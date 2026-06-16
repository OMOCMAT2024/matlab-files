% clear, clc, close all
import casadi.*

veh = veh_params();

% track_name = 'MyTrack';
% Track options: Austin, BrandsHatch, Budapest, Catalunya, Hockenheim, IMS,
% Melbourne, MexicoCity, Montreal, Monza, MoscowRaceway, Norisring,
% Nuerburgring, Oschersleben, Sakhir, SaoPaulo, Sepang, Shanghai,
% Silverstone, Sochi, Spa, Spielberg, Suzuka, YasMarina, Zandvoort

track = readtable(['tracks_smooth\' track_name '_smooth.csv']);

elapsedTime = 0; tic

%% Solver tolerance used only for soft complementarity-style constraints
OPT_e = 1e-3;

%% Vehicle model
% This file defines x,u,z,pv, dx,sf and output expressions.
% z is the lifted Chapter 9 tyre-load vector [Fz_fl Fz_fr Fz_rl Fz_rr].
run('tro_model_racecar.m');

f_dyn = Function('f_dyn',{x,u,z,pv},{dx,sf},{'x','u','z','pv'},{'dx','sf'});
f_sf  = Function('sf',{x,pv},{sf},{'x','pv'},{'sf'});

%% Optimal Control Problem - Path Constraints
% The corrected Chapter 9 load equations are enforced by h_load = 0.
% Friction is written without division by Fz^2 to improve IPOPT conditioning.

mux = veh.mux;
muy = veh.muy;

mufl_ineq = (fzfl^2 - (fxfl^2/mux^2 + fyfl^2/muy^2)) / veh.Fz_s^2;
mufr_ineq = (fzfr^2 - (fxfr^2/mux^2 + fyfr^2/muy^2)) / veh.Fz_s^2;
murl_ineq = (fzrl^2 - (fxrl^2/mux^2 + fyrl^2/muy^2)) / veh.Fz_s^2;
murr_ineq = (fzrr^2 - (fxrr^2/mux^2 + fyrr^2/muy^2)) / veh.Fz_s^2;

% Brake and throttle overlap, kept from original TOTPT.
BrTh = Tt_n * Tb_n;

% Engine power limits on driven rear wheels, kept from original TOTPT.
pwr_rl = Trl*omega_rl/(Tt_s*omega_s);
pwr_rr = Trr*omega_rr/(Tt_s*omega_s);
pwr_lb = -500e3/(Tt_s*omega_s);
pwr_ub =  150e3/(Tt_s*omega_s);

h_path = [mufl_ineq; mufr_ineq; murl_ineq; murr_ineq; h_load; BrTh; pwr_rl; pwr_rr];
h_lb = [0; 0; 0; 0; zeros(4,1); -OPT_e; pwr_lb; pwr_lb];
h_ub = [inf; inf; inf; inf; zeros(4,1);  OPT_e; pwr_ub; pwr_ub];
nh = length(h_lb);

h_eq = Function('h_eq', {x, u, z, pv}, {h_path}, {'x','u','z','pv'}, {'h'});

%% NLP - Collocation, polynomial coefficients

d = 3;
tau = collocation_points(d, 'legendre');
[C,D,B] = collocation_coeff(tau);
B = B(:);

%% NLP - Collocation, discretisation

N = length(track.s) - 1;
h = diff(track.s)';

step = 1:N+1;
step_col = kron(ones(1,N),tau) + kron(1:N,ones(1,d));
step_full = kron(ones(1,N),[0 tau]) + kron(1:N,ones(1,d+1));
step_full = [step_full N+1];

kappa_interp = interpolant('kappa_interp','linear',{step},track.kappa);
kappa_knot = kappa_interp(step);
kappa_col = kappa_interp(step_col);

wr_interp = interpolant('wr_interp','linear',{step},track.wr);
wr_knot = wr_interp(step); %#ok<NASGU>
wr_col = wr_interp(step_col); %#ok<NASGU>

wl_interp = interpolant('wl_interp','linear',{step},track.wl);
wl_knot = wl_interp(step); %#ok<NASGU>
wl_col = wl_interp(step_col); %#ok<NASGU>

s_interp = interpolant('s_interp','linear',{step},track.s); %#ok<NASGU>

%% NLP - formulation

Xk  = SX.sym('Xk',  nx,N+1);   % States at knots
Uk  = SX.sym('Uk',  nu,N+1);   % Inputs at knots
Zk  = SX.sym('Zk',  nz,N+1);   % Algebraic tyre-load variables at knots
Xkj = SX.sym('Xkj', nx,N*d);   % States at collocation points
Zkj = SX.sym('Zkj', nz,N*d);   % Algebraic tyre-load variables at collocation points

duk = diff(Uk')' ./ repmat(h,nu,1);
dzk = diff(Zk')' ./ repmat(h,nz,1);

%% NLP - boundary conditions

ux_start = 1;
Xi = [ux_start; nan; nan; nan; nan; nan; nan; nan; nan; nan; nan; nan; nan];
Xf = [nan;      nan; nan; nan; nan; nan; nan; nan; nan; nan; nan; nan; nan];

x0_min = x_min;
x0_max = x_max;
xf_min = x_min;
xf_max = x_max;
mask_i = ~isnan(Xi);
mask_f = ~isnan(Xf);
x0_min(mask_i) = max(x_min(mask_i), Xi(mask_i)./x_s(mask_i) - OPT_e);
x0_max(mask_i) = min(x_max(mask_i), Xi(mask_i)./x_s(mask_i) + OPT_e);
xf_min(mask_f) = max(x_min(mask_f), Xf(mask_f)./x_s(mask_f) - OPT_e);
xf_max(mask_f) = min(x_max(mask_f), Xf(mask_f)./x_s(mask_f) + OPT_e);

gb = [{Xk(:,1)}; {Xk(:,end)}];
lbg = [x0_min; xf_min];
ubg = [x0_max; xf_max];

%% NLP - initial guesses

ux_0 = ux_start*ones(1,N+1);
uy_0 = zeros(1,N+1);
r_0 = OPT_e*ones(1,N+1);
phi_0 = zeros(1,N+1);
phidot_0 = zeros(1,N+1);
theta_0 = zeros(1,N+1);
thetadot_0 = zeros(1,N+1);
n_0 = OPT_e*ones(1,N+1);
xi_0 = zeros(1,N+1);
omega_fl_0 = ux_0/veh.rw;
omega_fr_0 = ux_0/veh.rw;
omega_rl_0 = ux_0/veh.rw;
omega_rr_0 = ux_0/veh.rw;

Tt_0 = zeros(1,N+1);
Tb_0 = zeros(1,N+1);
delta_0 = OPT_e*ones(1,N+1);

Z10_0 = veh.m*veh.g*veh.lr/veh.l;
Z20_0 = veh.m*veh.g*veh.lf/veh.l;
fzfl_0 = 0.5*Z10_0*ones(1,N+1);
fzfr_0 = 0.5*Z10_0*ones(1,N+1);
fzrl_0 = 0.5*Z20_0*ones(1,N+1);
fzrr_0 = 0.5*Z20_0*ones(1,N+1);

x0 = [ux_0; uy_0; r_0; phi_0; phidot_0; theta_0; thetadot_0; n_0; xi_0; ...
      omega_fl_0; omega_fr_0; omega_rl_0; omega_rr_0] ./ x_s;
xc0_mat = kron(x0(:,1:end-1),ones(1,d));
xc0 = xc0_mat(:);

u0 = [Tt_0; Tb_0; delta_0] ./ u_s;
z0 = [fzfl_0; fzfr_0; fzrl_0; fzrr_0] ./ z_s;
zc0_mat = kron(z0(:,1:end-1),ones(1,d));
zc0 = zc0_mat(:);

%% NLP - collocation constraints, path constraints and objective function

dt = SX.zeros(1,N);
gck = {};
ghc = {};
J = 0;

for k = 0:N-1
    idx = d*k + (1:d);
    Xpoly = [Xk(:,k+1) Xkj(:,idx)];
    dPi = Xpoly*C;
    Xk_end = Xpoly*D;

    dtk = 0;
    for j = 1:d
        cidx = d*k + j;
        [dXkj_j, Qkj_j] = f_dyn(Xkj(:,cidx), Uk(:,k+1), Zkj(:,cidx), kappa_col(cidx));
        gck = [gck(:); {h(k+1)*dXkj_j - dPi(:,j)}];
        h_col_j = h_eq(Xkj(:,cidx), Uk(:,k+1), Zkj(:,cidx), kappa_col(cidx));
        ghc = [ghc(:); {h_col_j}];
        dtk = dtk + B(j)*Qkj_j*h(k+1);
    end

    gck = [gck(:); {Xk_end - Xk(:,k+2)}];
    J = J + dtk + sumsqr(rdu.*duk(:,k+1)) + sumsqr(rdy.*dzk(:,k+1));
    dt(k+1) = dtk;
end

%% Bounds at knot and collocation points
n_idx = 8;

Xk_min = repmat(x_min,N+1,1);
Xk_min(n_idx:nx:end) = full(-wr_interp(step) + (veh.wt/2 + veh.ws)) / n_s;
Xk_max = repmat(x_max,N+1,1);
Xk_max(n_idx:nx:end) = full( wl_interp(step) - (veh.wt/2 + veh.ws)) / n_s;

Xkj_min = repmat(x_min,N*d,1);
Xkj_min(n_idx:nx:end) = full(-wr_interp(step_col) + (veh.wt/2 + veh.ws)) / n_s;
Xkj_max = repmat(x_max,N*d,1);
Xkj_max(n_idx:nx:end) = full( wl_interp(step_col) - (veh.wt/2 + veh.ws)) / n_s;

Zk_min = repmat(z_min,N+1,1);
Zk_max = repmat(z_max,N+1,1);
Zkj_min = repmat(z_min,N*d,1);
Zkj_max = repmat(z_max,N*d,1);

%% Path constraints at knot points
ghk = {};
for k = 1:N+1
    ghk = [ghk(:); {h_eq(Xk(:,k), Uk(:,k), Zk(:,k), kappa_knot(k))}];
end

%% Rate of input constraints
gduk = {};
for k = 1:N
    sf_k = f_sf(Xk(:,k), kappa_knot(k));
    dsdt_k = 1 ./ sf_k;
    gduk = [gduk(:); {duk(:,k) * dsdt_k}];
end

%% NLP - define variables

w = [Xk(:); Uk(:); Zk(:); Xkj(:); Zkj(:)];
lbw = [Xk_min; repmat(u_min(:),N+1,1); Zk_min; Xkj_min; Zkj_min];
ubw = [Xk_max; repmat(u_max(:),N+1,1); Zk_max; Xkj_max; Zkj_max];
w0 = [x0(:); u0(:); z0(:); xc0(:); zc0(:)];

assert(size(lbw,1) == size(w,1));
assert(size(ubw,1) == size(w,1));
assert(size(w0,1) == size(w,1));

g = [gb(:); gck(:); ghk(:); ghc(:); gduk(:)];
g = vertcat(g{:});

lbg = [lbg; zeros((d+1)*N*nx,1); repmat(h_lb,N+1,1); repmat(h_lb,N*d,1); repmat(duk_lb,N,1)];
ubg = [ubg; zeros((d+1)*N*nx,1); repmat(h_ub,N+1,1); repmat(h_ub,N*d,1); repmat(duk_ub,N,1)];

assert(size(lbg,1) == size(g,1));
assert(size(ubg,1) == size(g,1));

%% NLP - output functions
f_dt_opt = Function('f_dt_opt',{w},{dt});
f_aero   = Function('f_aero',{x,u,z,pv},{y_aero});
f_tyre   = Function('f_tyre',{x,u,z,pv},{y_tyre});
f_slip   = Function('f_slip',{x,u,z,pv},{y_slip});
f_acc    = Function('f_acc',{x,u,z,pv},{y_acc});
f_torque = Function('f_torque',{x,u,z,pv},{y_torque});
f_mu     = Function('f_mu',{x,u,z,pv},{y_mu});
f_pwr    = Function('f_pwr',{x,u,z,pv},{y_pwr});
f_xdot   = Function('f_xdot',{x,u,z,pv},{y_xdot});
f_load   = Function('f_load',{x,u,z,pv},{y_load});

%% NLP - solve
opts = struct;
opts.ipopt.max_iter = 1000;
opts.ipopt.tol = 1e-3;

nlp = struct('f', J, 'x', w, 'g', g);
solver = nlpsol('solver', 'ipopt', nlp, opts);

elapsedTime(end+1) = toc;
tic;

sol_nlp = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

elapsedTime(end+1) = toc;
tic;

%% Postprocessing - Collect NLP variables

sol = struct;
sol.w_opt = full(sol_nlp.x);

idx_x0  = 0;
idx_u0  = idx_x0 + nx*(N+1);
idx_z0  = idx_u0 + nu*(N+1);
idx_xc0 = idx_z0 + nz*(N+1);
idx_zc0 = idx_xc0 + nx*N*d;

sol.x_opt  = reshape(sol.w_opt(idx_x0  + (1:nx*(N+1))),nx,N+1) .* x_s;
sol.u_opt  = reshape(sol.w_opt(idx_u0  + (1:nu*(N+1))),nu,N+1) .* u_s;
sol.z_opt  = reshape(sol.w_opt(idx_z0  + (1:nz*(N+1))),nz,N+1) .* z_s;
sol.xc_opt = reshape(sol.w_opt(idx_xc0 + (1:nx*N*d)),nx,N*d) .* x_s;
sol.zc_opt = reshape(sol.w_opt(idx_zc0 + (1:nz*N*d)),nz,N*d) .* z_s;

%% Postprocessing - Reconstruct solution at all interpolation points

sol.x_full = kron(sol.x_opt(:,1:N), [1 zeros(1,d)]) + reshape([zeros(nx,N); reshape(sol.xc_opt,nx*d,N)],nx,[]);
sol.x_full = [sol.x_full sol.x_opt(:,N+1)];
sol.z_full = kron(sol.z_opt(:,1:N), [1 zeros(1,d)]) + reshape([zeros(nz,N); reshape(sol.zc_opt,nz*d,N)],nz,[]);
sol.z_full = [sol.z_full sol.z_opt(:,N+1)];
sol.u_full = interp1(step',sol.u_opt',step_full)';

%% Postprocessing - Output variables at knot points

sol.aero   = zeros(N+1,4);
sol.tyre   = zeros(N+1,12);
sol.slip   = zeros(N+1,8);
sol.acc    = zeros(N+1,7);
sol.torque = zeros(N+1,4);
sol.mu     = zeros(N+1,4);
sol.pwr    = zeros(N+1,4);
sol.xdot   = zeros(N+1,nx);
sol.load   = zeros(N+1,11);

for k = 1:N+1
    xk_scaled = sol.x_opt(:,k)./x_s;
    uk_scaled = sol.u_opt(:,k)./u_s;
    zk_scaled = sol.z_opt(:,k)./z_s;
    kap = kappa_knot(k);

    sol.aero(k,:)   = full(f_aero(xk_scaled,uk_scaled,zk_scaled,kap))';
    sol.tyre(k,:)   = full(f_tyre(xk_scaled,uk_scaled,zk_scaled,kap))';
    sol.slip(k,:)   = full(f_slip(xk_scaled,uk_scaled,zk_scaled,kap))';
    sol.acc(k,:)    = full(f_acc(xk_scaled,uk_scaled,zk_scaled,kap))';
    sol.torque(k,:) = full(f_torque(xk_scaled,uk_scaled,zk_scaled,kap))';
    sol.mu(k,:)     = full(f_mu(xk_scaled,uk_scaled,zk_scaled,kap))';
    sol.pwr(k,:)    = full(f_pwr(xk_scaled,uk_scaled,zk_scaled,kap))';
    sol.xdot(k,:)   = (full(f_xdot(xk_scaled,uk_scaled,zk_scaled,kap)) .* x_s)';
    sol.load(k,:)   = full(f_load(xk_scaled,uk_scaled,zk_scaled,kap))';
end

%% Postprocessing - lap time

dt_opt_val = full(f_dt_opt(sol.w_opt));
sol.t_opt = [0 cumsum(dt_opt_val)];

%% Postprocessing - Reconstruct optimal racetrack

n_opt = sol.x_opt(n_idx,:)';
wr_opt = track.wr + n_opt;
wl_opt = track.wl - n_opt;

nor = [-sin(track.phi) cos(track.phi)];
track.x_opt = track.x + diag(nor(:,1))*n_opt;
track.y_opt = track.y + diag(nor(:,2))*n_opt;
[s_opt, ~] = calc_length_xy([track.x_opt,track.y_opt]);
[phi_opt, kappa_opt] = calc_head_curv_num([track.x_opt,track.y_opt]);

v_opt = sqrt(sol.x_opt(1,:)'.^2 + sol.x_opt(2,:)'.^2);
sol.track_opt = [track.x_opt track.y_opt wr_opt wl_opt s_opt phi_opt kappa_opt v_opt];

%% Save data
track_opt_table = array2table(sol.track_opt);
track_opt_table.Properties.VariableNames = {'x','y','wr','wl','s','phi','kappa','v'};

if ~exist(fullfile(pwd,'tracks_mintime'), 'dir')
   mkdir('tracks_mintime')
end
writetable(track_opt_table,['tracks_mintime\' track_name '_mintime' '.csv'])

if ~exist(fullfile(pwd,'tro_sol'), 'dir')
   mkdir('tro_sol')
end
save(['tro_sol\sol_' track_name '.mat'],'sol')
