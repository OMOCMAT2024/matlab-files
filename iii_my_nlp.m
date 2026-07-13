import casadi.*

veh = my_params();

track = track_interp_table;

elapsedTime = 0; tic %Start time counter

%%

OPT_e = 1e-5;

% Preserve the existing non-periodic formulation by default when my_nlp.m
% is run directly instead of through my_oc.m.
if ~exist('use_periodic_state_bc','var')
    use_periodic_state_bc = false;
end

%% Vehicle model

run('my_model.m');

f_dyn = Function('f_dyn',{x,u,z,pv},{dx,sf},{'x','u','z','pv'},{'dx','sf'});

f_sf = Function('sf',{x,pv},{sf},{'x','pv'},{'sf'});

%% Path Constraints

pw_drive_rear_normalized = Tt * (omega_rl+omega_rr)/2 / (Tt_s*omega_s);

% Available power is increased only in Rosberg's distance-based KERS-on
% intervals. The final interval uses its independently tunable parameter.
kers_power_add_available = veh.kers_power_add * kers_on + ...
    (veh.kers_power_add_last - veh.kers_power_add) * kers_last_on;
power_motor_available = veh.power_motor_max + kers_power_add_available;
pw_drive_rear_ub_normalized = ...
    power_motor_available / (Tt_s*omega_s); % [Watts] / [Watts]

BrTh = Tt_n*Tb_n;

Here_body = pw_drive_rear_normalized - pw_drive_rear_ub_normalized;

slip_ratio_value = veh.slip_ratio_value;

slip_body_fl = ((1-slip_ratio_value)*(ux-0.5*r*veh.t1) - omega_fl*veh.rw) / (omega_s*veh.rw);

slip_body_fr = ((1-slip_ratio_value)*(ux+0.5*r*veh.t1) - omega_fr*veh.rw) / (omega_s*veh.rw);

slip_body_rl = ((1-slip_ratio_value)*(ux-0.5*r*veh.t2) - omega_rl*veh.rw) / (omega_s*veh.rw);

slip_body_rr = ((1-slip_ratio_value)*(ux+0.5*r*veh.t2) - omega_rr*veh.rw) / (omega_s*veh.rw);

h_path = [BrTh; Here_body; slip_body_fl; slip_body_fr; slip_body_rl; slip_body_rr; h_load];

h_lb = [0-OPT_e; -inf; -inf; -inf; -inf; -inf; zeros(4,1)-OPT_e];

h_ub = [0+OPT_e; 0; 0; 0; 0; 0; zeros(4,1)+OPT_e];

nh = length(h_lb);

% Path constraints

h_eq = Function('h_eq', {x, u, z, pv}, {h_path}, {'x','u','z','pv'}, {'h'});

%% NLP - colloc, poly coefficients

d = 3;

tau = collocation_points(d, 'legendre');

[C,D,B] = collocation_coeff(tau);

B = B(:);

%% NLP - discretisation

N = length(track.s) - 1;

h = diff(track.s)';

step = 1:N+1;

step_col = kron(ones(1,N),tau) + kron(1:N,ones(1,d));

step_full = kron(ones(1,N),[0 tau]) + kron(1:N,ones(1,d+1));

step_full = [step_full N+1];

kappa_interp = interpolant('kappa_interp','linear',{step},track.curvature);

kappa_knot = kappa_interp(step);

kappa_col = kappa_interp(step_col);

% Cd and KERS flags are piecewise constant in distance because DRS/KERS
% are on/off states. Keep the left-knot value throughout each finite
% element to avoid non-physical intermediate values at transitions.
Cd_knot = reshape(track.Cd, 1, N+1);
Cd_col = reshape(repmat(Cd_knot(1:N), d, 1), 1, N*d);
kers_on_knot = reshape(track.kers_on, 1, N+1);
kers_on_col = reshape(repmat(kers_on_knot(1:N), d, 1), 1, N*d);
kers_last_on_knot = reshape(track.kers_last_on, 1, N+1);
kers_last_on_col = reshape( ...
    repmat(kers_last_on_knot(1:N), d, 1), 1, N*d);

wr_interp = interpolant('wr_interp','linear',{step},track.wr);

wr_knot = wr_interp(step);

wr_col = wr_interp(step_col);

wl_interp = interpolant('wl_interp','linear',{step},track.wl);

wl_knot = wl_interp(step);

wl_col = wl_interp(step_col);

s_interp = interpolant('s_interp','linear',{step},track.s);

%% NLP - formulation

Xk = SX.sym('Xk', nx,N+1); % States at knots

Uk = SX.sym('Uk', nu,N+1); % Inputs at knots

Zk = SX.sym('Zk', nz,N+1); % Algebraic tyre-load variables at knots

Xkj = SX.sym('Xkj', nx,N*d); % States at collocation points

Zkj = SX.sym('Zkj', nz,N*d); % Algebraic tyre-load variables at collocation points

% Knot-to-knot scaled input jump. The physical input is u = u_s .* Uk,
% so du_jump is dimensionless. It is used below for the Betts-style
% algebraic variable rate constraints.
du_jump = Uk(:,2:N+1) - Uk(:,1:N);

% Optional distance-rate quantities retained only for diagnostics or optional
% smoothing; the hard rate constraints below are imposed in time, not in s.
duk = du_jump ./ repmat(h,nu,1);

dzk = diff(Zk')' ./ repmat(h,nz,1);

%% NLP - boundary conditions

% if vehicle starts at center of the track, n0 = 0;

% ux_start = 5; % m/s
ux_start = 43.95; % m/s

% uy_start = 1e-9; % m/s
uy_start = -0.00108; % m/s

% Theta_start = deg2rad(1e-9);
Theta_start = -0.00385; % [rad]

Thetadot_start = 0.000287; % [rad/s]

% r_start = 1e-9;
r_start = 0.00289;

Phidot_start = 1e-5; % [rad/s]

% Phi_start = deg2rad(1e-9);
Phi_start = 6.487e-6; % [rad]

n_start = -4; % [m]

xi_start = deg2rad(0.754);

omegaFL_start = 126.879; % [rad/s]
omegaFR_start = 126.879; % [rad/s]
omegaRL_start = 130.008; % [rad/s]
omegaRR_start = 130.008; % [rad/s]

% Xi = [ux_start; uy_start; r_start; Phi_start; nan; Theta_start; nan; nan; nan; ux_start/veh.rw; ux_start/veh.rw; ux_start/veh.rw; ux_start/veh.rw];
Xi = [ux_start; nan; r_start; Phi_start; nan; Theta_start; nan; n_start; xi_start; omegaFL_start; omegaFR_start; omegaRL_start; omegaRR_start];

%-final

ux_end = nan;

% ux_end = ux_start;

Xf = [ux_end; nan; nan; nan; nan; nan; nan; nan; nan; nan; nan; nan; nan];

x0_min = max(x_min, Xi./x_s-OPT_e);

x0_max = min(x_max, Xi./x_s+OPT_e);

xf_min = max(x_min, Xf./x_s-OPT_e);

xf_max = min(x_max, Xf./x_s+OPT_e);

%-boundary constraints

% gb = [{Xk(:,1)}; {Xk(:,end)}]; lbg = [x0_min; xf_min]; ubg = [x0_max; xf_max];

gb = [{Xk(:,1)}; {Xk(:,end)}];

lbg = [x0_min; xf_min];

ubg = [x0_max; xf_max];

%% NLP - initial guesses (constant values)

%-States

% ux_0 = ux_start*ones(1,N+1);
ux_0 = 5*ones(1,N+1);

uy_0 = OPT_e*ones(1,N+1);

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

%-Inputs

Tt_0 = 100*ones(1,N+1);

Tb_0 = OPT_e*ones(1,N+1);

delta_0 = OPT_e*ones(1,N+1);

Z10_0 = veh.m*veh.g*veh.lr/veh.l;

Z20_0 = veh.m*veh.g*veh.lf/veh.l;

fzfl_0 = 0.5*Z10_0*ones(1,N+1);

fzfr_0 = 0.5*Z10_0*ones(1,N+1);

fzrl_0 = 0.5*Z20_0*ones(1,N+1);

fzrr_0 = 0.5*Z20_0*ones(1,N+1);

% Collect initial guesses

x0 = [ux_0; uy_0; r_0; phi_0; phidot_0; theta_0; thetadot_0; n_0; xi_0; ...

omega_fl_0; omega_fr_0; omega_rl_0; omega_rr_0] ./ x_s;

xc0_mat = kron(x0(:,1:end-1),ones(1,d));

xc0 = xc0_mat(:);

u0 = [Tt_0; Tb_0; delta_0] ./ u_s;

z0 = [fzfl_0; fzfr_0; fzrl_0; fzrr_0] ./ z_s;

zc0_mat = kron(z0(:,1:end-1),ones(1,d));

zc0 = zc0_mat(:);

%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % ==== Better initial guess for drift window ====

% % build a knot-window too (length N+1)

%

% index_1_corners_a = 82 - haha; % 2m interval

% index_2_corners_a = 1023 - haha; % 2m interval

%

% index_1_corners_b = 1024 - haha; % 2m interval

% index_2_corners_b = 1965 - haha; % 2m interval

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%

% k_trans = 1; % transition width in *intervals* (tune 5~20), 8

% beta_ref = deg2rad(45); % [rad] target drift magnitude in strong corners deg2rad(10, 15, 25, 40)

% kappa0 = 0.01; % [1/m] curvature scale for "how corner-ish" it is (tune)

%

% Wbeta = 50.0; % drift weight 200.0 75.0 90.0 45.0 22.5 50.0 0.0

%

% % ==== Drift window (smooth) built from knot indices ====

% k1_corners_a = index_1_corners_a - 1; % interval index start (your loop uses k=0..N-1)

% k2_corners_a = index_2_corners_a - 1; % interval index end

%

% k1_corners_b = index_1_corners_b - 1;

% k2_corners_b = index_2_corners_b - 1;

%

% sig = @(x) 0.5*(tanh(x) + 1); % smooth step in [0,1]

%% NLP - formulation, collocation constraints, path constraints and objective function

%

% w_drift_corners_a = zeros(1,N); % one weight per interval k

% w_drift_corners_b = zeros(1,N);

%

% for kk = 0:N-1

% w_drift_corners_a(kk+1) = sig((kk - k1_corners_a)/k_trans) * sig((k2_corners_a - kk)/k_trans);

% w_drift_corners_b(kk+1) = sig((kk - k1_corners_b)/k_trans) * sig((k2_corners_b - kk)/k_trans);

% end

%

% % Monitored variables

% dt = SX.zeros(1,N);

%

% % Loop over discretisation points to create collocation constraints

% gck = {}; % Collector for collocation constraints

% J = 0; % Initialise objective function

%

% for k = 0:N-1

% % Collocation constraints

% X = [Xk(:,k+1) Xkj(:,d*k+(1:d))]; % Concatenate states

%

% %-Dynamics

% %%-calculate derivatives of the approximating polynomial at the collocation points

% dPi = X*C;

%

% %%-calculate derivatives of the system at the collocation points

% [dXkj, Qk] = f_dyn(Xkj(:,d*k+(1:d)), Uk(:,k+1), kappa_col(d*k+(1:d)));

%

% %-State of approximating polynomial the end of the colocation interval

% Xk_end = X*D;

%

% gck = [gck(:); {h(k+1)*dXkj(:) - dPi(:)}; {Xk_end-Xk(:,k+2)}]; % Add collocation constraints

%

% % Integrate quadrature function

% dtk = Qk*B*h(k+1);

%

% wk_corners_a = w_drift_corners_a(k+1); % numeric scalar in [0,1]

% wk_corners_b = w_drift_corners_b(k+1);

%

% % % smoothness penalties (keep them; optionally relax them in-window)

% % rdu_k_corners_a = (1 - 0.8*wk_corners_a) * rdu; % optional: 80% less smoothing in drift window

% % rdy_k_corners_a = (1 - 0.8*wk_corners_a) * rdy;

%

% % ---- Drift encouragement (simple & smooth) ----

% beta_c = Xkj(2, d*k+1:(d*k+d)); % beta at collocation points (1xd)

% kappa_c = kappa_col_shifted(d*k+1:(d*k+d));

%

% % % beta_des = beta_ref * tanh(kappa_c / kappa0); % (1xd) smooth signed target

% beta_des = (-1) * beta_ref * tanh(kappa_c / kappa0); % (1xd) smooth signed target

% % % beta_des = (1) * beta_ref * tanh(kappa_c / kappa0); % this seems to have problems

%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% beta_des_mat(k+1, :) = full(beta_des); % beta_des is 1x3 casadi.DM

% % beta_des_mat is now a 1x3 double row vector

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%

% J_drift_k = ( (beta_c - beta_des).^2 ) * B * h(k+1);

%

% % J = J + wk * Wbeta * J_drift_k;

%

% time_weight_when_drifting = 1.0;

% % if k >= k1_corners_a && k <= k2_corners_a

% % % if k >= 198 && k <= 291

% %

% % % J = J + wk * Wbeta * J_drift_k;

% %

% % % J = J + wk * Wbeta * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % J = J + 1 * wk_corners_a * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % % J = J + 1 * wk_corners_a * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk;

% %

% % elseif k >= k1_corners_b && k <= k2_corners_b

% %

% % J = J + 1 * wk_corners_b * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % % J = J + 1 * wk_corners_b * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk;

% %

% % % % elseif k >= k1_corners_c && k <= k2_corners_c

% % % %

% % % % J = J + 1 * wk_corners_c * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % %

% % % elseif k >= k1_corners_d && k <= k2_corners_d

% % %

% % % J = J + 1 * wk_corners_d * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % %

% % % elseif k >= k1_corners_e && k <= k2_corners_e

% % %

% % % J = J + 1 * wk_corners_e * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % %

% % % elseif k >= k1_corners_f && k <= k2_corners_f

% % %

% % % J = J + 1 * wk_corners_f * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % %

% % % elseif k >= k1_corners_g && k <= k2_corners_g

% % %

% % % J = J + 1 * wk_corners_g * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % % %

% % % % % elseif k >= k1_corners_h && k <= k2_corners_h

% % % % %

% % % % % J = J + 1 * wk_corners_h * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% % % % %

% % % % % elseif k >= k1_corners_i && k <= k2_corners_i

% % % % %

% % % % % J = J + 1 * wk_corners_i * Wbeta*(1/200) * J_drift_k + time_weight_when_drifting*dtk + sumsqr(rdu.*duk(:,k+1));

% %

% % else

% %

% % J = J + dtk + sumsqr(rdu.*duk(:,k+1));

% % % J = J + dtk;

% %

% % end

%

% % J = J + dtk + sumsqr(rdu.*duk(:,k+1));

% % J = J + dtk + (1e-2)*sumsqr(rdu.*duk(:,k+1));

% J = J + dtk + (1e-6)*sumsqr(rdu.*duk(:,k+1));

%

% % Collect contribution to time of each interval

% dt(k+1) = dtk;

%

% end

%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure()

% plot(track.curvature)

% hold on

% plot(beta_des_mat(:,1))

% % hold on

% % plot(sol.x_opt(2, :))

% % hold on

% % plot(sol.u_opt(3,:))

% hold on

% plot(w_drift_corners_a*beta_ref)

% hold on

% plot(w_drift_corners_b*beta_ref)

% % legend('track.curvature at knots', '\beta desired at collocation points', '\beta solution', '\delta solution', 'wdrift window scaled', 'wdrift window scaled', 'wdrift window scaled', 'wdrift window scaled', 'wdrift window scaled', 'wdrift window scaled', 'wdrift window scaled', 'wdrift window scaled', 'wdrift window scaled')

% legend('curvature', 'beta desired', 'w drift corners a * beta ref', 'w drift corners b * beta ref')

% xlabel('Index')

% % ylabel('\beta desired, rad/s')

% title('\beta desired with wdrift window')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectorised collocation transcription. This is mathematically the same as
% the previous nested k/j loops, but avoids repeated cell-array growth and
% uses CasADi mapped functions for all collocation points at once.

% Repeat each knot input over the d collocation points of its interval:
% [U_1 ... U_N] -> [U_1 repeated d times, U_2 repeated d times, ...].
Uk_col = reshape(repmat(Uk(:,1:N), d, 1), nu, N*d);

% Per-collocation interval lengths and quadrature weights in the same column
% order as Xkj/Zkj: [k1j1 k1j2 ... k1jd k2j1 ...].
h_col = reshape(repmat(h, d, 1), 1, N*d);
B_col = repmat(B.', 1, N);
kappa_col_vec = reshape(kappa_col, 1, N*d);
pv_col = [kappa_col_vec; Cd_col; kers_on_col; kers_last_on_col];

% Collocation polynomial values stacked interval-by-interval:
% column k of Xpoly_stack is [Xk(:,k); Xkj(:,d*(k-1)+1); ...; Xkj(:,d*k)].
Xpoly_stack = [Xk(:,1:N); reshape(Xkj, nx*d, N)];

% Linear polynomial maps for all intervals at once.
% dPi_stack(:,k) = [dPi(:,1); ...; dPi(:,d)] for interval k.
C_blk = kron(full(C.'), eye(nx));
D_blk = kron(full(D.'), eye(nx));
dPi_stack = C_blk * Xpoly_stack;
dPi = reshape(dPi_stack, nx, N*d);
Xk_end = D_blk * Xpoly_stack;

% Dynamics and path constraints at all collocation points.
f_dyn_col = f_dyn.map(N*d);
h_eq_col = h_eq.map(N*d);
[dXkj_col, Qkj_col] = f_dyn_col(Xkj, Uk_col, Zkj, pv_col);
h_col_path = h_eq_col(Xkj, Uk_col, Zkj, pv_col);

% Collocation residuals arranged by finite element/stage.
% Each column k is exactly the same interval-k residual block as the old
% transcription:
% [colloc j=1..d residuals; interval-end residual].
% Keeping this matrix form lets the final NLP constraint vector be stacked
% stage-by-stage, which makes the Jacobian/KKT pattern more local.
coll_res = dXkj_col .* repmat(h_col, nx, 1) - dPi;
end_res = Xk_end - Xk(:,2:N+1);
gck_stage = [reshape(coll_res, nx*d, N); end_res];

% Collocation path constraints arranged by finite element/stage.
% Each column k is [h(k,j=1); ...; h(k,j=d)].
ghc_stage = reshape(h_col_path, nh*d, N);

% Interval time increments and objective.
dt_col = Qkj_col .* B_col .* h_col;
dt = ones(1,d) * reshape(dt_col, d, N);
% Minimum-time objective. The input-rate limits are enforced as hard
% Betts-style constraints below, so there is no active soft penalty on duk.
% Optional extra smoothing can be re-enabled deliberately if desired, e.g.
% J = dt * ones(N,1) + 1e-6*sumsqr(repmat(rdu,1,N).*duk);
% J = dt * ones(N,1) + sumsqr(repmat(rdu,1,N).*duk) + sumsqr(repmat(rdy,1,N).*dzk);
J = dt * ones(N,1);

%% Bounds at knot and collocation points

n_idx = 8; % this is because state variable 'n_n' is the eighth state in current model formulation

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

kappa_knot_vec = reshape(kappa_knot, 1, N+1);
pv_knot = [kappa_knot_vec; Cd_knot; ...
    kers_on_knot; kers_last_on_knot];
h_eq_knot = h_eq.map(N+1);
h_knot_path = h_eq_knot(Xk, Uk, Zk, pv_knot);
ghk = {h_knot_path(:)};

%% Rate of input constraints, Betts Eq. (4.261)-(4.262) style

% Betts writes the algebraic-variable rate limit over a time interval as
% 0 <= u_{k+1} - u_k - h_k*r_L,
% 0 >= u_{k+1} - u_k - h_k*r_U.
% This NLP is transcribed with distance s as the independent variable, so
% h_k must be the interval time, not the distance mesh spacing h(k). The
% interval time dt(k) has already been computed above by collocation
% quadrature from sf = dt/ds:
% dt(k) = integral_{s_k}^{s_{k+1}} sf ds.
%
% The controls Uk are scaled, i.e. physical u = u_s .* Uk, and my_model.m
% defines duk_lb/duk_ub as physical input-rate bounds divided by u_s.
% Therefore dt(k)*duk_lb/duk_ub has the same dimensionless units as
% Uk(:,k+1)-Uk(:,k).
dt_rate_mat = repmat(dt, nu, 1);
rate_lb_mat = repmat(duk_lb(:), 1, N);
rate_ub_mat = repmat(duk_ub(:), 1, N);

gdu_lower_mat = du_jump - rate_lb_mat .* dt_rate_mat; % >= 0
gdu_upper_mat = du_jump - rate_ub_mat .* dt_rate_mat; % <= 0
gdu_stage = [gdu_lower_mat; gdu_upper_mat];

%% NLP - define NLP variables

% Stage-ordered decision vector for a more banded/stage-structured KKT
% matrix. This is only a permutation of the same variables:
% stage k = [Xk(:,k); Uk(:,k); Zk(:,k); Xkj(:,k,1:d); Zkj(:,k,1:d)]
% followed by the final knot [Xk(:,N+1); Uk(:,N+1); Zk(:,N+1)].
% The mathematical NLP is unchanged; only the variable ordering seen by
% Ipopt's sparse linear algebra is changed.
stage_nw = nx + nu + nz + d*nx + d*nz;
final_nw = nx + nu + nz;

Xk_min_mat = reshape(Xk_min, nx, N+1);
Xk_max_mat = reshape(Xk_max, nx, N+1);
Zk_min_mat = reshape(Zk_min, nz, N+1);
Zk_max_mat = reshape(Zk_max, nz, N+1);
Xkj_min_stage = reshape(Xkj_min, nx*d, N);
Xkj_max_stage = reshape(Xkj_max, nx*d, N);
Zkj_min_stage = reshape(Zkj_min, nz*d, N);
Zkj_max_stage = reshape(Zkj_max, nz*d, N);

w_stage = [Xk(:,1:N); Uk(:,1:N); Zk(:,1:N); reshape(Xkj,nx*d,N); reshape(Zkj,nz*d,N)];
w = [w_stage(:); Xk(:,N+1); Uk(:,N+1); Zk(:,N+1)];

lbw_stage = [Xk_min_mat(:,1:N); repmat(u_min(:),1,N); Zk_min_mat(:,1:N); Xkj_min_stage; Zkj_min_stage];
ubw_stage = [Xk_max_mat(:,1:N); repmat(u_max(:),1,N); Zk_max_mat(:,1:N); Xkj_max_stage; Zkj_max_stage];
lbw = [lbw_stage(:); Xk_min_mat(:,N+1); u_min(:); Zk_min_mat(:,N+1)];
ubw = [ubw_stage(:); Xk_max_mat(:,N+1); u_max(:); Zk_max_mat(:,N+1)];

w0_stage = [x0(:,1:N); u0(:,1:N); z0(:,1:N); reshape(xc0_mat,nx*d,N); reshape(zc0_mat,nz*d,N)];
w0 = [w0_stage(:); x0(:,N+1); u0(:,N+1); z0(:,N+1)];

assert(size(lbw,1) == size(w,1));

assert(size(ubw,1) == size(w,1));

assert(size(w0,1) == size(w,1));

% Stage-ordered constraint vector. This is the same set of constraints as
% before, but permuted from
% [boundary; all dynamics; all knot path; all collocation path; all rates]
% to
% [initial boundary/path; interval-local blocks; final boundary].
% Each interval block contains only variables from the current interval and
% the next knot, improving locality of the Jacobian and therefore the KKT
% matrix used by Ipopt.
g_interval_stage = [gck_stage; ghc_stage; gdu_stage; h_knot_path(:,2:N+1)];

lbg_rate_stage = [zeros(nu,1); -inf(nu,1)];
ubg_rate_stage = [ inf(nu,1); zeros(nu,1)];
lbg_interval_stage = [zeros((d+1)*nx,1); repmat(h_lb,d,1); lbg_rate_stage; h_lb];
ubg_interval_stage = [zeros((d+1)*nx,1); repmat(h_ub,d,1); ubg_rate_stage; h_ub];

if use_periodic_state_bc
    % Flying-lap formulation: close all scaled vehicle states over one lap.
    % Equality in scaled coordinates is exactly equality in physical units.
    g = [h_knot_path(:,1); g_interval_stage(:); Xk(:,N+1)-Xk(:,1)];
    lbg = [h_lb; repmat(lbg_interval_stage,N,1); zeros(nx,1)];
    ubg = [h_ub; repmat(ubg_interval_stage,N,1); zeros(nx,1)];
else
    % Existing non-periodic formulation, unchanged.
    g = [Xk(:,1); h_knot_path(:,1); g_interval_stage(:); Xk(:,N+1)];
    lbg = [x0_min; h_lb; repmat(lbg_interval_stage,N,1); xf_min];
    ubg = [x0_max; h_ub; repmat(ubg_interval_stage,N,1); xf_max];
end

assert(size(lbg,1) == size(g,1));

assert(size(ubg,1) == size(g,1));

assert(size(w,1) == stage_nw*N + final_nw);

%% constrain the initial states of speed and wheel angular rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~use_periodic_state_bc
    delta_ux = 0.1; % [m/s]

    % Correspondingly enforce the initial longitudinal speed only for the
    % existing non-periodic formulation.
    lbw(1) = (ux_start-delta_ux) / x_s(1);
    ubw(1) = (ux_start+delta_ux) / x_s(1);
end

% % enforce initial rotation rate of wheels
%
% lbw(10) = (ux_start-delta_ux) / veh.rw / x_s(10);
%
% ubw(10) = (ux_start+delta_ux) / veh.rw / x_s(10);
%
% lbw(11) = (ux_start-delta_ux) / veh.rw / x_s(11);
%
% ubw(11) = (ux_start+delta_ux) / veh.rw / x_s(11);
%
% lbw(12) = (ux_start-delta_ux) / veh.rw / x_s(12);
%
% ubw(12) = (ux_start+delta_ux) / veh.rw / x_s(12);
%
% lbw(13) = (ux_start-delta_ux) / veh.rw / x_s(13);
%
% ubw(13) = (ux_start+delta_ux) / veh.rw / x_s(13);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NLP - output functions

f_dt_opt = Function('f_dt_opt',{w},{dt});

f_aero = Function('f_aero',{x,u,z,pv},{y_aero});

f_tyre = Function('f_tyre',{x,u,z,pv},{y_tyre});

f_slip = Function('f_slip',{x,u,z,pv},{y_slip});

f_acc = Function('f_acc',{x,u,z,pv},{y_acc});

f_torque = Function('f_torque',{x,u,z,pv},{y_torque});

f_mu = Function('f_mu',{x,u,z,pv},{y_mu});

f_pwr = Function('f_pwr',{x,u,z,pv},{y_pwr});

f_xdot = Function('f_xdot',{x,u,z,pv},{y_xdot});

f_load = Function('f_load',{x,u,z,pv},{y_load});

% f_control_rate = Function('f_control_rate',{x,u,z,pv},{gdu_stage}); %
% this will give error

%% NLP - solve

opts = struct;

opts.ipopt.max_iter = 1500; %max number of iterations. % 1500

opts.ipopt.mu_init = 1e0;

opts.ipopt.tol = 1e-6;

% opts.ipopt.dual_inf_tol = 1e-4;

% opts.ipopt.constr_viol_tol = 1e-6;

opts.ipopt.compl_inf_tol = 1e-8;

%-nlp problem [objective function 'f', decision variables 'x', constraints 'g']

nlp = struct('f', J, 'x', w, 'g', g);

% Create solver

solver = nlpsol('solver', 'ipopt', nlp, opts);

elapsedTime(end+1) = toc;

tic;

% Solve the NLP

sol_nlp = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

elapsedTime(end+1) = toc;

tic

%% Postprocessing - Collect NLP variables

sol = struct;

sol.w_opt = full(sol_nlp.x);

% Unpack the stage-ordered decision vector back to the original arrays used
% by the rest of the postprocessing/plotting code.
stage_opt = reshape(sol.w_opt(1:stage_nw*N), stage_nw, N);
final_opt = sol.w_opt(stage_nw*N + (1:final_nw));

row_x = 1:nx;
row_u = nx + (1:nu);
row_z = nx + nu + (1:nz);
row_xc = nx + nu + nz + (1:nx*d);
row_zc = nx + nu + nz + nx*d + (1:nz*d);

sol.x_opt = [stage_opt(row_x,:), final_opt(1:nx)] .* x_s;
sol.u_opt = [stage_opt(row_u,:), final_opt(nx + (1:nu))] .* u_s;
sol.z_opt = [stage_opt(row_z,:), final_opt(nx + nu + (1:nz))] .* z_s;
sol.xc_opt = reshape(stage_opt(row_xc,:), nx, N*d) .* x_s;
sol.zc_opt = reshape(stage_opt(row_zc,:), nz, N*d) .* z_s;

%% Postprocessing - Reconstruct solution at all interpolation points

sol.x_full = kron(sol.x_opt(:,1:N), [1 zeros(1,d)]) + reshape([zeros(nx,N); reshape(sol.xc_opt,nx*d,N)],nx,[]);

sol.x_full = [sol.x_full sol.x_opt(:,N+1)];

sol.z_full = kron(sol.z_opt(:,1:N), [1 zeros(1,d)]) + reshape([zeros(nz,N); reshape(sol.zc_opt,nz*d,N)],nz,[]);

sol.z_full = [sol.z_full sol.z_opt(:,N+1)];

sol.u_full = interp1(step',sol.u_opt',step_full)';

%% Postprocessing - Output variables at knot points

sol.aero = zeros(N+1,4);

sol.tyre = zeros(N+1,12);

sol.slip = zeros(N+1,8);

sol.acc = zeros(N+1,7);

sol.torque = zeros(N+1,4);

sol.mu = zeros(N+1,4);

sol.pwr = zeros(N+1,4);

sol.xdot = zeros(N+1,nx);

sol.load = zeros(N+1,11);

% Distance-based DRS/KERS quantities used at the knot points.
sol.Cd = full(Cd_knot);
sol.kers_on = full(kers_on_knot);
sol.kers_last_on = full(kers_last_on_knot);
sol.kers_power_add_W = veh.kers_power_add * sol.kers_on + ...
    (veh.kers_power_add_last - veh.kers_power_add) * sol.kers_last_on;
sol.power_motor_available_W = ...
    veh.power_motor_max + sol.kers_power_add_W;

for k = 1:N+1

xk_scaled = sol.x_opt(:,k)./x_s;

uk_scaled = sol.u_opt(:,k)./u_s;

zk_scaled = sol.z_opt(:,k)./z_s;

pv_here = [kappa_knot(k); Cd_knot(k); ...
    kers_on_knot(k); kers_last_on_knot(k)];

sol.aero(k,:) = full(f_aero(xk_scaled,uk_scaled,zk_scaled,pv_here))';

sol.tyre(k,:) = full(f_tyre(xk_scaled,uk_scaled,zk_scaled,pv_here))';

sol.slip(k,:) = full(f_slip(xk_scaled,uk_scaled,zk_scaled,pv_here))';

sol.acc(k,:) = full(f_acc(xk_scaled,uk_scaled,zk_scaled,pv_here))';

sol.torque(k,:) = full(f_torque(xk_scaled,uk_scaled,zk_scaled,pv_here))';

sol.mu(k,:) = full(f_mu(xk_scaled,uk_scaled,zk_scaled,pv_here))';

sol.pwr(k,:) = full(f_pwr(xk_scaled,uk_scaled,zk_scaled,pv_here))';

sol.xdot(k,:) = (full(f_xdot(xk_scaled,uk_scaled,zk_scaled,pv_here)) .* x_s)';

sol.load(k,:) = full(f_load(xk_scaled,uk_scaled,zk_scaled,pv_here))';

end

%% Postprocessing - Get lap time

dt_opt_val = full(f_dt_opt(sol.w_opt));

sol.t_opt = [0 cumsum(dt_opt_val)];

%% Postprocessing - Reconstruct optimal raceline

n_opt = sol.x_opt(n_idx,:)';

wr_opt = track.wr + n_opt;

wl_opt = track.wl - n_opt;

% nor = [-sin(track.phi) cos(track.phi)];

nor = [-sin(track.Theta) cos(track.Theta)];

track.x_opt = track.x + diag(nor(:,1))*n_opt;

track.y_opt = track.y + diag(nor(:,2))*n_opt;

[s_opt, ~] = calc_length_xy([track.x_opt,track.y_opt]);

[phi_opt, kappa_opt] = calc_head_curv_num([track.x_opt,track.y_opt]);

v_opt = sqrt(sol.x_opt(1,:)'.^2 + sol.x_opt(2,:)'.^2);

sol.track_opt = [track.x_opt track.y_opt wr_opt wl_opt s_opt phi_opt kappa_opt v_opt];

%% Save data for my visualization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = sol.t_opt.'; % [s]

carX = sol.track_opt(:, 1); % [m]

carY = sol.track_opt(:, 2);

carZ = sol.track_opt(:, 1) * 0;

% speed = sol.x_opt(1, :).'; % [m/s]

% sideslip_angle = sol.x_opt(2, :).'; % [rad]

speed = ((sol.x_opt(1, :).').^2 + (sol.x_opt(2, :).').^2).^0.5; % [m/s]

sideslip_angle = atan(sol.x_opt(2, :).' ./ sol.x_opt(1, :).'); % [rad]

% heading_unwraped = unwrap(track.Theta) + sol.x_opt(5,:).'; % [rad]

heading_unwraped = unwrap(track.Theta) + sol.x_opt(9,:).'; % [rad]

heading = wrapToPi(heading_unwraped); % wrap to [-pi, pi]

heading4Carsim = heading_unwraped;

% Extract states

pitch = sol.x_opt(6,:).'; % [rad]

roll = sol.x_opt(4,:).'; % [rad]

leftX = P_L_interp_here(:,1); % [m]

leftY = P_L_interp_here(:,2);

leftZ = P_L_interp_here(:,1) * 0;

rightX = P_R_interp_here(:,1);

rightY = P_R_interp_here(:,2);

rightZ = P_R_interp_here(:,1) * 0;

deltaFL = sol.u_opt(3,:).'; % [rad]

deltaFR = deltaFL;

omegaFL = sol.x_opt(10,:).'; % [rad/s]

omegaFR = sol.x_opt(11,:).';

omegaRL = sol.x_opt(12,:).';

omegaRR = sol.x_opt(13,:).';

Tdrive_here = sol.u_opt(1,:).';

Tbrake_here = (-1) * sol.u_opt(2,:).';

Tdrive_norm = Tdrive_here / veh.Tt_max;

Tbrake_norm = Tbrake_here / abs(veh.Tb_min);

cd('D:\Applications\TOTPT-main\oc\my_replay')

save('data4replay', 't', 'carX', 'carY', 'carZ', 'speed', 'sideslip_angle', 'heading', 'roll','pitch', 'leftX', 'leftY', 'leftZ', 'rightX', 'rightY', 'rightZ','deltaFL', 'deltaFR', 'omegaFL', 'omegaFR', 'omegaRL', 'omegaRR', 'Tdrive_norm', 'Tbrake_norm')

cd('D:\Applications\TOTPT-main\oc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug = 0;