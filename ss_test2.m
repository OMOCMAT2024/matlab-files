%% steady_state_circle_max_ay.m
% Steady-state constant-radius circular test:
% compute the maximum achievable lateral acceleration for a given radius
% using the same vehicle model structure as your current CasADi OCP.
%
% REQUIREMENT:
%   - my_params.m must be on the MATLAB path
%   - CasADi must be installed and on the MATLAB path
%
% NOTES:
%   - R_input > 0 : left turn
%   - R_input < 0 : right turn
%   - steady state is enforced by xdot = 0 (time-domain equilibrium)
%   - path centreline condition is enforced by n = 0
%   - your rear-axle power constraint is preserved
%   - explicit wheel-load positivity constraints are added

clear
close all
clc
import casadi.*

veh = my_params();

%% ========================= USER INPUTS ==================================
R_input = 40.0;         % [m], + left turn, - right turn
OPT_e   = 1e-6;

% Increase this if you do not want the auxiliary ax_bar / ay_bar states
% to artificially cap the solution.
AXAY_AUX_LIMIT_G = 5.0; % [g]

% Enforce positive wheel loads
FZ_MIN = 1.0;           % [N]

% Keep your current rear-axle power limit
USE_POWER_LIMIT = true;

DO_PLOTS = true;

if abs(R_input) < 1e-12
    error('R_input must be nonzero.');
end

kappa_val = 1 / R_input;
R_abs     = abs(R_input);
turn_sign = sign(R_input);

%% ========================= SCALING ======================================
V_s      = 100;
beta_s   = 1;
gamma_s  = 1;
n_s      = 5;
xi_s     = 1;
ax_bar_s = 2*9.8;
ay_bar_s = 2*9.8;
omega_s  = V_s / veh.rw;

Tt_s     = 2e3;
Tb_s     = 4e3;
delta_s  = pi/8;

%% ========================= STATES / INPUTS ==============================
nx = 11;
nu = 4;

% States (scaled)
V_n        = SX.sym('V_n');
beta_n     = SX.sym('beta_n');
gamma_n    = SX.sym('gamma_n');
n_n        = SX.sym('n_n');
xi_n       = SX.sym('xi_n');
ax_bar_n   = SX.sym('ax_bar_n');
ay_bar_n   = SX.sym('ay_bar_n');
omega_fl_n = SX.sym('omega_fl_n');
omega_fr_n = SX.sym('omega_fr_n');
omega_rl_n = SX.sym('omega_rl_n');
omega_rr_n = SX.sym('omega_rr_n');

x = [V_n; beta_n; gamma_n; n_n; xi_n; ax_bar_n; ay_bar_n; ...
     omega_fl_n; omega_fr_n; omega_rl_n; omega_rr_n];

x_s = [V_s; beta_s; gamma_s; n_s; xi_s; ax_bar_s; ay_bar_s; ...
       omega_s; omega_s; omega_s; omega_s];

% Inputs (scaled)
Tt_n              = SX.sym('Tt_n');
Tb_n              = SX.sym('Tb_n');
delta_n           = SX.sym('delta_n');
omega_rear_axle_n = SX.sym('omega_rear_axle_n');

u = [Tt_n; Tb_n; delta_n; omega_rear_axle_n];
u_s = [Tt_s; Tb_s; delta_s; omega_s];

% Unscaled physical variables
V               = V_s * V_n;
beta            = beta_s * beta_n;
gamma           = gamma_s * gamma_n;
n               = n_s * n_n;
xi              = xi_s * xi_n;
ax_bar          = ax_bar_s * ax_bar_n;
ay_bar          = ay_bar_s * ay_bar_n;
omega_fl        = omega_s * omega_fl_n;
omega_fr        = omega_s * omega_fr_n;
omega_rl        = omega_s * omega_rl_n;
omega_rr        = omega_s * omega_rr_n;

Tt              = Tt_s * Tt_n;
Tb              = Tb_s * Tb_n;
delta           = delta_s * delta_n;
omega_rear_axle = omega_s * omega_rear_axle_n;

%% ========================= BOUNDS =======================================
omega_max = veh.V_max / veh.rw;

x_min = [   -1e-3;
            -pi/4;
            -pi/2;
            -4;
            -pi/4;
            -2*9.8;
            -2*9.8;
             0;
             0;
             0;
             0] ./ x_s;

x_max = [veh.V_max;
          pi/4;
          pi/2;
          4;
          pi/4;
          2*9.8;
          2*9.8;
          omega_max;
          omega_max;
          omega_max;
          omega_max] ./ x_s;

% Force forward motion
x_min(1) = max(x_min(1), 0.5 / V_s);

% Relax auxiliary-state bounds so they do not clip the result
aux_lim = AXAY_AUX_LIMIT_G * veh.g;
x_min(6) = -aux_lim / ax_bar_s;
x_max(6) =  aux_lim / ax_bar_s;
x_min(7) = -aux_lim / ay_bar_s;
x_max(7) =  aux_lim / ay_bar_s;

u_min = [veh.Tt_min;
         veh.Tb_min;
         veh.delta_min;
         veh.omega_rear_axle_min] ./ u_s;

u_max = [veh.Tt_max;
         veh.Tb_max;
         veh.delta_max;
         veh.omega_rear_axle_max] ./ u_s;

%% ========================= MODEL EQUATIONS ==============================
cosd = cos(delta);
sind = sin(delta);
cosb = cos(beta);
sinb = sin(beta);

% Vehicle velocities
vx = V*cosb;
vy = V*sinb;

% Aero
f_drag = 0.5 * veh.drag_coeff * vx^2;
f_lift = 0.5 * veh.lift_coeff * vx^2;

% Tyre contact-patch velocities
vxfl = (vx - 0.5*veh.wt*gamma)*cosd + (vy + veh.lf*gamma)*sind + 1e-12;
vxfr = (vx + 0.5*veh.wt*gamma)*cosd + (vy + veh.lf*gamma)*sind + 1e-12;
vxrl = (vx - 0.5*veh.wt*gamma) + 1e-12;
vxrr = (vx + 0.5*veh.wt*gamma) + 1e-12;

vyfl = (vy + veh.lf*gamma)*cosd - (vx - 0.5*veh.wt*gamma)*sind;
vyfr = (vy + veh.lf*gamma)*cosd - (vx + 0.5*veh.wt*gamma)*sind;
vyrl = (vy - veh.lr*gamma);
vyrr = (vy - veh.lr*gamma);

% Wheel loads
Fx_tires = ax_bar * veh.m + f_drag;
Fy_tires = ay_bar * veh.m;

Fdown   = f_lift;
Fdrag   = f_drag;
acom    = veh.lf;
bcom    = veh.lr;
mass    = veh.m;
gravity = veh.g;
hcom    = veh.hc;
twF_h   = veh.twF_h;
twR_h   = veh.twR_h;
dr      = veh.dr;

acop = acom;
hcop = hcom;

[fzfl, fzfr, fzrl, fzrr] = wheel_load( ...
    Fx_tires, Fy_tires, Fdown, Fdrag, ...
    bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h);

fzf = fzfl + fzfr;
fzr = fzrl + fzrr;

% Slips
SXfl = (veh.rw*omega_fl - vxfl)/vxfl;
SXfr = (veh.rw*omega_fr - vxfr)/vxfr;
SXrl = (veh.rw*omega_rl - vxrl)/vxrl;
SXrr = (veh.rw*omega_rr - vxrr)/vxrr;

syfl = atan(vyfl/vxfl);
syfr = atan(vyfr/vxfr);
syrl = atan(vyrl/vxrl);
syrr = atan(vyrr/vxrr);

tan_syfl = vyfl/vxfl;
tan_syfr = vyfr/vxfr;
tan_syrl = vyrl/vxrl;
tan_syrr = vyrr/vxrr;

[fxfl, fyfl, mufl] = tire_force_tanEllipse(veh, SXfl, tan_syfl, fzfl);
[fxfr, fyfr, mufr] = tire_force_tanEllipse(veh, SXfr, tan_syfr, fzfr);
[fxrl, fyrl, murl] = tire_force_tanEllipse(veh, SXrl, tan_syrl, fzrl);
[fxrr, fyrr, murr] = tire_force_tanEllipse(veh, SXrr, tan_syrr, fzrr);

% Wheel torques
dfzf = 0.5*(fzfr - fzfl);
dfzr = 0.5*(fzrr - fzrl);

kb  = 0.65;
Tf  = kb*Tb;
Tfl = 0.5*Tf*(1 - dfzf/fzf);
Tfr = 0.5*Tf*(1 + dfzf/fzf);

Tr  = Tt + (1-kb)*Tb;
Trl = 0.5*Tr*(1 - dfzr/fzr);
Trr = 0.5*Tr*(1 + dfzr/fzr);

% Force resultants and yaw moment
X_force = (fxfl + fxfr)*cosd - (fyfl + fyfr)*sind + (fxrl + fxrr) - f_drag;
Y_force = (fxfl + fxfr)*sind + (fyfl + fyfr)*cosd + (fyrl + fyrr);

ax = X_force / veh.m;
ay = Y_force / veh.m;

Mz = (-1)*fxfl*cosd*(veh.wt/2) + fxfl*sind*veh.lf + fyfl*sind*(veh.wt/2) + fyfl*cosd*veh.lf + ...
          fxfr*cosd*(veh.wt/2) + fxfr*sind*veh.lf - fyfr*sind*(veh.wt/2) + fyfr*cosd*veh.lf + ...
     (-1)*fxrl*(veh.wt/2) - fyrl*veh.lr + ...
          fxrr*(veh.wt/2) - fyrr*veh.lr;

% Change of variable helper
chi = xi + beta;
sf  = (1 - n*kappa_val)/(V*cos(chi) + 1e-12) + 1e-12;

% Time-domain state derivatives (scaled)
xdot = [
    (X_force*cosb + Y_force*sinb) / veh.m;
    (Y_force*cosb - X_force*sinb) / (V*veh.m) - gamma;
    Mz / veh.Iz;
    V*sin(chi);
    gamma - kappa_val/sf;
    (ax - ax_bar) / veh.load_transfer_time_delay;
    (ay - ay_bar) / veh.load_transfer_time_delay;
    (Tfl - fxfl*veh.rw) / veh.Iw;
    (Tfr - fxfr*veh.rw) / veh.Iw;
    (Trl - fxrl*veh.rw) / veh.Iw;
    (Trr - fxrr*veh.rw) / veh.Iw
] ./ x_s;

%% ========================= PATH CONSTRAINTS =============================
% Same sign convention as your current model
BrTh = Tt_n * Tb_n;

Surething = omega_rear_axle_n - (omega_rl_n + omega_rr_n)/2;

power_E_max = power_E_max_func(omega_rear_axle);

pw_drive_rear_normalized    = Tt * omega_rear_axle / (Tt_s*omega_s);
pw_drive_rear_ub_normalized = power_E_max / (Tt_s*omega_s);

pw_brake_rear_normalized = Tb * (omega_rl + omega_rr)/2 / (Tt_s*omega_s);
pw_brake_rear_lb_normalized = (-500e3*2) / (Tt_s*omega_s);

Here_body = pw_drive_rear_normalized - pw_drive_rear_ub_normalized;

if USE_POWER_LIMIT
    h    = [BrTh; Surething; Here_body; pw_brake_rear_normalized];
    h_lb = [0-OPT_e; 0-OPT_e; -inf; pw_brake_rear_lb_normalized];
    h_ub = [0+OPT_e; 0+OPT_e;    0; 0];
else
    h    = [BrTh; Surething; pw_brake_rear_normalized];
    h_lb = [0-OPT_e; 0-OPT_e; pw_brake_rear_lb_normalized];
    h_ub = [0+OPT_e; 0+OPT_e; 0];
end

% Additional physically important constraints
wheel_load_con = [fzfl; fzfr; fzrl; fzrr];
wheel_load_lb  = FZ_MIN * ones(4,1);
wheel_load_ub  = inf(4,1);

%% ========================= OBJECTIVE ===================================
% For fixed radius, maximizing V is equivalent to maximizing V^2/|R|
% and is numerically simpler.
reg = 1e-7*(beta_n^2 + gamma_n^2 + delta_n^2) + ...
      1e-8*(Tt_n^2 + Tb_n^2) + ...
      1e-8*((omega_fl_n - V_n)^2 + (omega_fr_n - V_n)^2 + ...
            (omega_rl_n - V_n)^2 + (omega_rr_n - V_n)^2);

J = -V_n + reg;

%% ========================= NLP =========================================
% Enforce:
%   1) steady state      xdot = 0
%   2) on centreline     n = 0
%   3) path constraints  h
%   4) positive loads    Fz >= FZ_MIN
g = [xdot; n_n; h; wheel_load_con];

lbg = [zeros(nx,1); 0; h_lb; wheel_load_lb];
ubg = [zeros(nx,1); 0; h_ub; wheel_load_ub];

w   = [x; u];
lbw = [x_min; u_min];
ubw = [x_max; u_max];

nlp = struct('f', J, 'x', w, 'g', g);

opts = struct;
opts.ipopt.max_iter       = 3000;
opts.ipopt.tol            = 1e-8;
opts.ipopt.acceptable_tol = 1e-6;
opts.ipopt.compl_inf_tol  = 1e-8;
opts.ipopt.mu_init        = 1e-2;
opts.ipopt.print_level    = 0;
opts.print_time           = false;

solver = nlpsol('solver', 'ipopt', nlp, opts);

%% ========================= EVALUATION FUNCTIONS =========================
a_lat_geom_signed = V^2 / R_input;   % signed path-normal accel
a_lat_geom_mag    = V^2 / R_abs;     % magnitude
chi_expr          = xi + beta;

f_basic = Function('f_basic', {w}, {[
    V;
    beta;
    gamma;
    n;
    xi;
    chi_expr;
    ax;
    ay;
    ax_bar;
    ay_bar;
    delta;
    Tt;
    Tb;
    omega_fl;
    omega_fr;
    omega_rl;
    omega_rr;
    omega_rear_axle;
    a_lat_geom_signed;
    a_lat_geom_mag
]});

f_tyre = Function('f_tyre', {w}, { ...
    [fxfl; fxfr; fxrl; fxrr; fyfl; fyfr; fyrl; fyrr]});

f_fz = Function('f_fz', {w}, { ...
    [fzfl; fzfr; fzrl; fzrr]});

f_slip = Function('f_slip', {w}, { ...
    [SXfl; SXfr; SXrl; SXrr; syfl; syfr; syrl; syrr]});

f_mu = Function('f_mu', {w}, { ...
    [mufl; mufr; murl; murr]});

f_cons = Function('f_cons', {w}, {g});

%% ========================= MULTI-START =================================
V_guess_max = min(veh.V_max, max(10, sqrt(AXAY_AUX_LIMIT_G*veh.g*R_abs)));
V_guess_min = 5;

V_guess_list = unique([linspace(V_guess_min, V_guess_max, 8), ...
                       0.7*V_guess_max, 0.9*V_guess_max]);

wheelbase = veh.lf + veh.lr;

best.found   = false;
best.metric  = -inf;

for i = 1:length(V_guess_list)
    V0 = V_guess_list(i);

    beta0  = turn_sign * deg2rad(2.0);
    gamma0 = V0 / R_input;
    n0     = 0.0;
    xi0    = -beta0;
    axb0   = 0.0;
    ayb0   = V0^2 / R_input;

    omega0 = max(0.1, V0 / veh.rw);
    delta0 = turn_sign * atan(wheelbase / R_abs);

    Tt0 = 100.0;
    Tb0 = 0.0;
    omega_rear_axle0 = omega0;

    x0 = [V0; beta0; gamma0; n0; xi0; axb0; ayb0; ...
          omega0; omega0; omega0; omega0] ./ x_s;

    u0 = [Tt0; Tb0; delta0; omega_rear_axle0] ./ u_s;

    x0 = min(max(x0, lbw(1:nx) + 1e-6), ubw(1:nx) - 1e-6);
    u0 = min(max(u0, lbw(nx+1:end) + 1e-6), ubw(nx+1:end) - 1e-6);

    w0 = [x0; u0];

    try
        sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
        stats = solver.stats();

        ok = false;
        if isfield(stats, 'success') && stats.success
            ok = true;
        elseif isfield(stats, 'return_status')
            rs = stats.return_status;
            if contains(rs, 'Solve_Succeeded') || contains(rs, 'Solved_To_Acceptable_Level')
                ok = true;
            end
        end

        if ~ok
            continue
        end

        w_opt = full(sol.x);
        g_opt = full(f_cons(w_opt));

        total_viol = max([lbg - g_opt; g_opt - ubg; 0]);
        if total_viol > 1e-5
            continue
        end

        basic = full(f_basic(w_opt));
        basic = basic(:);
        V_opt_here = basic(1);

        if ~best.found || V_opt_here > best.metric
            best.found     = true;
            best.w_opt     = w_opt;
            best.stats     = stats;
            best.metric    = V_opt_here;
            best.total_viol = total_viol;
        end

    catch
        % try next initial guess
    end
end

if ~best.found
    error(['No feasible steady-state solution found. ', ...
           'Try a larger radius, looser AXAY_AUX_LIMIT_G, ', ...
           'or check whether your steering/power/load limits make this radius infeasible.']);
end

%% ========================= POSTPROCESS =================================
basic = full(f_basic(best.w_opt));
basic = basic(:);

tyreF = full(f_tyre(best.w_opt));
Fz    = full(f_fz(best.w_opt));
slip  = full(f_slip(best.w_opt));
muRat = full(f_mu(best.w_opt));
g_opt = full(f_cons(best.w_opt));

V_opt               = basic(1);
beta_opt            = basic(2);
gamma_opt           = basic(3);
n_opt               = basic(4);
xi_opt              = basic(5);
chi_opt             = basic(6);
ax_opt              = basic(7);
ay_opt              = basic(8);
ax_bar_opt          = basic(9);
ay_bar_opt          = basic(10);
delta_opt           = basic(11);
Tt_opt              = basic(12);
Tb_opt              = basic(13);
omega_fl_opt        = basic(14);
omega_fr_opt        = basic(15);
omega_rl_opt        = basic(16);
omega_rr_opt        = basic(17);
omega_rear_axle_opt = basic(18);
a_lat_signed_opt    = basic(19);
a_lat_mag_opt       = basic(20);

fprintf('\n============================================================\n');
fprintf('Steady-state circular test result\n');
fprintf('============================================================\n');
fprintf('Radius R                         = %+12.6f m\n', R_input);
fprintf('Curvature kappa                  = %+12.6f 1/m\n', kappa_val);
fprintf('Vehicle speed V                  = %+12.6f m/s\n', V_opt);
fprintf('Vehicle speed                    = %+12.6f km/h\n', 3.6*V_opt);
fprintf('Sideslip beta                    = %+12.6f rad  (%+10.6f deg)\n', beta_opt, rad2deg(beta_opt));
fprintf('Yaw rate gamma                   = %+12.6f rad/s\n', gamma_opt);
fprintf('Course-angle error xi            = %+12.6f rad  (%+10.6f deg)\n', xi_opt, rad2deg(xi_opt));
fprintf('chi = xi + beta                  = %+12.6e rad  (%+10.6e deg)\n', chi_opt, rad2deg(chi_opt));
fprintf('Steer angle delta                = %+12.6f rad  (%+10.6f deg)\n', delta_opt, rad2deg(delta_opt));
fprintf('n                                = %+12.6e m\n', n_opt);
fprintf('ax (body frame)                  = %+12.6f m/s^2  (%+10.6f g)\n', ax_opt, ax_opt/veh.g);
fprintf('ay (body frame)                  = %+12.6f m/s^2  (%+10.6f g)\n', ay_opt, ay_opt/veh.g);
fprintf('ax_bar                           = %+12.6f m/s^2  (%+10.6f g)\n', ax_bar_opt, ax_bar_opt/veh.g);
fprintf('ay_bar                           = %+12.6f m/s^2  (%+10.6f g)\n', ay_bar_opt, ay_bar_opt/veh.g);
fprintf('V^2/R (signed)                   = %+12.6f m/s^2  (%+10.6f g)\n', a_lat_signed_opt, a_lat_signed_opt/veh.g);
fprintf('|V^2/R|                          = %+12.6f m/s^2  (%+10.6f g)\n', a_lat_mag_opt, a_lat_mag_opt/veh.g);
fprintf('Drive torque Tt                  = %+12.6f Nm\n', Tt_opt);
fprintf('Brake torque Tb                  = %+12.6f Nm\n', Tb_opt);
fprintf('omega_rear_axle                  = %+12.6f rad/s\n', omega_rear_axle_opt);
fprintf('omega_rear_axle                  = %+12.6f rpm\n', omega_rear_axle_opt*30/pi);
fprintf('Minimum wheel load               = %+12.6f N\n', min(Fz));
fprintf('Maximum constraint violation     = %.3e\n', best.total_viol);
fprintf('============================================================\n\n');

result = struct();
result.R               = R_input;
result.kappa           = kappa_val;
result.V               = V_opt;
result.beta            = beta_opt;
result.gamma           = gamma_opt;
result.n               = n_opt;
result.xi              = xi_opt;
result.chi             = chi_opt;
result.ax              = ax_opt;
result.ay              = ay_opt;
result.ax_bar          = ax_bar_opt;
result.ay_bar          = ay_bar_opt;
result.delta           = delta_opt;
result.Tt              = Tt_opt;
result.Tb              = Tb_opt;
result.omega_fl        = omega_fl_opt;
result.omega_fr        = omega_fr_opt;
result.omega_rl        = omega_rl_opt;
result.omega_rr        = omega_rr_opt;
result.omega_rear_axle = omega_rear_axle_opt;
result.a_lat_signed    = a_lat_signed_opt;
result.a_lat_mag       = a_lat_mag_opt;
result.tyre_forces     = tyreF;
result.wheel_loads     = Fz;
result.slips           = slip;
result.mu_ratio        = muRat;
result.constraint_vals = g_opt;
result.w_opt           = best.w_opt;
result.solver_stats    = best.stats;

save('steady_state_circle_result.mat', 'result');

%% ========================= OPTIONAL PLOTS ===============================
if DO_PLOTS
    figure('Name','Steady-State Circle Test','Color','w');

    subplot(2,2,1)
    bar([tyreF(1:4), tyreF(5:8)])
    grid on
    xlabel('Wheel')
    ylabel('Force [N]')
    legend('Fx','Fy','Location','best')
    title('Tyre Forces')
    xticklabels({'FL','FR','RL','RR'})

    subplot(2,2,2)
    bar(Fz)
    grid on
    xlabel('Wheel')
    ylabel('Fz [N]')
    title('Wheel Loads')
    xticklabels({'FL','FR','RL','RR'})

    subplot(2,2,3)
    bar([slip(1:4), rad2deg(slip(5:8))])
    grid on
    xlabel('Wheel')
    ylabel('Slip')
    legend('\kappa','\alpha [deg]','Location','best')
    title('Slip Quantities')
    xticklabels({'FL','FR','RL','RR'})

    subplot(2,2,4)
    bar(muRat)
    grid on
    xlabel('Wheel')
    ylabel('\mu ratio')
    title('Embedded Tyre Saturation Ratio')
    xticklabels({'FL','FR','RL','RR'})
end

%% ========================= LOCAL FUNCTIONS ==============================
function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load( ...
    Fx_tires, Fy_tires, Fdown, Fdrag, ...
    bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h)

Fz_FL = ((bcom*(mass*gravity+Fdown) - (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) ...
         /(acom+bcom) + ((Fy_tires*hcom)*(-dr))/(twF_h*(-dr) - twR_h*(1-dr))) / 2;

Fz_FR = ((bcom*(mass*gravity+Fdown) - (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) ...
         /(acom+bcom) - ((Fy_tires*hcom)*(-dr))/(twF_h*(-dr) - twR_h*(1-dr))) / 2;

Fz_RL = ((acom*(mass*gravity+Fdown) + (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) ...
         /(acom+bcom) + ((Fy_tires*hcom)*(1-dr))/(twR_h*(1-dr) - (-dr)*twF_h)) / 2;

Fz_RR = ((acom*(mass*gravity+Fdown) + (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) ...
         /(acom+bcom) - ((Fy_tires*hcom)*(1-dr))/(twR_h*(1-dr) - (-dr)*twF_h)) / 2;
end

function [Fx, Fy, mu_ratio] = tire_force_tanEllipse(veh, kappa, tan_alpha, Fz)

Qx   = veh.Qx;
Qy   = veh.Qy;
eps1 = veh.eps1;
eps2 = veh.eps2;

k_max         = tan(pi/(2*veh.Cx))/veh.Bx;
alpha_peak    = tan(pi/(2*veh.Cy))/veh.By;
tan_alpha_max = tan(alpha_peak);

[mu_x_max, mu_y_max] = muxmax_muymax_func(veh, Fz);

k_n = kappa / k_max;
a_n = (-tan_alpha) / tan_alpha_max;

rou = sqrt(k_n.^2 + a_n.^2 + eps1) + eps2;

Sx = pi/(2*atan(Qx));
Sy = pi/(2*atan(Qy));

my_effective_mu_scaling = veh.my_effective_mu_scaling;

Fx = (my_effective_mu_scaling * mu_x_max) .* Fz .* sin(Qx * atan(Sx * rou)) .* (k_n ./ rou);
Fy = (my_effective_mu_scaling * mu_y_max) .* Fz .* sin(Qy * atan(Sy * rou)) .* (a_n ./ rou);

mu_ratio = hypot( ...
    Fx ./ ((my_effective_mu_scaling * mu_x_max).*Fz), ...
    Fy ./ ((my_effective_mu_scaling * mu_y_max).*Fz));
end

function [mu_x_max, mu_y_max] = muxmax_muymax_func(veh, Fz)

Fz1 = veh.Fz1;
Fz2 = veh.Fz2;
rou_smooth = veh.rou_smooth;

mu_x_max1 = (veh.d1x*Fz1 + veh.d2x)/Fz1;
mu_x_max2 = (veh.d1x*Fz2 + veh.d2x)/Fz2;
mu_y_max1 = (veh.d1y*Fz1 + veh.d2y)/Fz1;
mu_y_max2 = (veh.d1y*Fz2 + veh.d2y)/Fz2;

mu_x_max_raw = (Fz - Fz1) * (mu_x_max2 - mu_x_max1) / (Fz2 - Fz1) + mu_x_max1;
mu_y_max_raw = (Fz - Fz1) * (mu_y_max2 - mu_y_max1) / (Fz2 - Fz1) + mu_y_max1;

mu_x_max = mu_x_max_raw + smooth_ramp_func(mu_x_max_raw, rou_smooth);
mu_y_max = mu_y_max_raw + smooth_ramp_func(mu_y_max_raw, rou_smooth);
end

function smooth_ramp = smooth_ramp_func(x, rou_smooth)
smooth_ramp = (-x) .* (atan2(-x, rou_smooth)/pi + 1/2) + rou_smooth/pi;
end

function power_E_max = power_E_max_func(omega_rear_axle)
power_E_max = ...
    (0.5*tanh(2.5*omega_rear_axle - 581.776417331) - 0.5*tanh(2.5*omega_rear_axle - 707.7571987)) .* ...
        (1286.78621646*omega_rear_axle + 5.91995335186*omega_rear_axle.^2 - 0.0212738654539*omega_rear_axle.^3 - 33396.5368135) ...
    - 1.0*(0.5*tanh(2.5*omega_rear_axle - 581.776417331) - 0.5*tanh(2.5*omega_rear_axle - 442.079344477)) .* ...
        (1565.43335336*omega_rear_axle + 8.76142302002*omega_rear_axle.^2 - 0.038302838172*omega_rear_axle.^3 - 33396.5368135) ...
    + (0.5*tanh(2.5*omega_rear_axle - 267.852862491) - 0.5*tanh(2.5*omega_rear_axle - 338.2421031)) .* ...
        (2692.54536777*omega_rear_axle + 25.9197938624*omega_rear_axle.^2 - 0.194902000291*omega_rear_axle.^3 - 33396.5368135) ...
    - 1.0*(0.5*tanh(2.5*omega_rear_axle - 267.852862491) - 0.5*tanh(2.5*omega_rear_axle - 173.147743253)) .* ...
        (3400.12124349*omega_rear_axle + 41.3327490565*omega_rear_axle.^2 - 0.39247357235*omega_rear_axle.^3 - 33396.5368135) ...
    + (0.5*tanh(2.5*omega_rear_axle - 338.2421031) - 0.5*tanh(2.5*omega_rear_axle - 442.079344477)) .* ...
        (2060.11029302*omega_rear_axle + 15.1735230258*omega_rear_axle.^2 - 0.0872968601385*omega_rear_axle.^3 - 33396.5368135) ...
    - (0.5*tanh(2.5*omega_rear_axle - 173.147743253) - 0.5*tanh(2.5*omega_rear_axle - 110.814555682)) .* ...
        (5259.85606727*omega_rear_axle + 98.9129613268*omega_rear_axle.^2 - 1.45294370534*omega_rear_axle.^3 - 33396.5368135) ...
    - 1.0*(0.5*tanh(2.5*omega_rear_axle - 110.814555682) - 0.5) .* ...
        (8218.52510512*omega_rear_axle + 241.486721989*omega_rear_axle.^2 - 5.54254037986*omega_rear_axle.^3 - 33396.5368135) ...
    + (0.5*tanh(2.5*omega_rear_axle - 707.7571987) + 0.5) .* ...
        (1001.87734615*omega_rear_axle + 3.588678869*omega_rear_axle.^2 - 0.0100408592098*omega_rear_axle.^3 - 33396.5368135);
end

function veh = my_params()

veh = struct();

veh.mu = 1.0;
veh.g = 9.80665;

veh.m = 1496;
veh.Iz = 2241;
veh.lf = 1.217;
veh.lr = 1.231;
veh.l = veh.lf+veh.lr;
veh.hc = 0.45;
veh.Iw = 15/2;
veh.rw = 0.32;

rho = 1.2250;
A = 2.2;
Cd = 0.32;
Cl = 0.0;

veh.drag_coeff = rho*Cd*A;
veh.lift_coeff = rho*Cl*A;

veh.wt = 1.59;
veh.ws = 0.2;

veh.twF_h = veh.wt / 2; % front-wheel to car-centre-line distance
veh.twR_h = veh.wt / 2; %  rear-wheel to car-centre-line distance

veh.Ix = 550;
veh.Iy = 2000;
veh.Kphi = 75000; % Nm/rad
veh.Cphi = 4500; % Nm/rad
veh.Ktheta = 120000; % Nm/rad
veh.Ctheta = 10000; % Nm/rad

veh.roll_pitch_ode_initial_condition = [0 0 0 0].';


%% Brake and torque distribution
% Tf = kt*Tt + kb*Tb
% Tr = (1-kt)*Tt + (%1*kb)*Tb
% kt = 0

veh.kb = 0.65;

%% Roll stifness

veh.dr = 0.5;

veh.load_transfer_time_delay = 1e-1; % [s]

%% Tire parameters

veh.Bx = 18;
veh.Cx = 1.3;
veh.d1x = 0.95;
veh.d2x = 320;

veh.By = 13;
veh.Cy = 1.5;
veh.d1y = 0.95;
veh.d2y = 320;

% % veh.Bx=80; 
% % veh.Cx=2.8; 
% % veh.d1x=0.94582704; 
% % veh.d2x=4.4982565;
% % 
% % veh.By=14.014685;
% % veh.Cy=1.5002282; 
% % veh.d1y=0.88447465; 
% % veh.d2y=128.8949;

veh.mux = 1.0;
veh.muy = 1.0;

veh.Fz1 = 2000;   % N
veh.Fz2 = 6000;   % N
veh.Qx  = 1.35;
veh.Qy  = 1.5;
veh.eps1 = 1e-5;
veh.eps2 = 1e-5;

veh.my_effective_mu_scaling = 1.0; % 0.9749

% smoothing parameters for soft nonnegative clamp
veh.rou_smooth = 1e-2;
% veh.offset     = 0.00318298;
veh.offset     = veh.rou_smooth / pi;


%% Scaling factors
veh.V_s = 100;
veh.beta_s = 1;
veh.gamma_s = 1;
veh.ax_bar_s = veh.g;
veh.ay_bar_s = veh.g;
veh.s_s = 30;
veh.n_s = 0.1;
veh.xi_s = 1;
veh.Tt_s = 1000*2;
% veh.Tb_s = 1000*4;
veh.Tb_s = (-1)*4e3;
veh.delta_s = pi/8;
veh.omega_s = veh.V_s/veh.rw;


%% State an input limits
% states
veh.V_min = 10/3.6;

veh.V_max = 250/3.6;
veh.beta_min = -pi/4;
veh.beta_max =  pi/4;
veh.gamma_min = -pi/2;
veh.gamma_max =  pi/2;
veh.ax_bar_min = -3*veh.g;
veh.ax_bar_max =  3*veh.g;
veh.ay_bar_min = -3*veh.g;
veh.ay_bar_max =  3*veh.g;
veh.s_min =   0;
veh.s_max = 500;
veh.n_min = -4;
veh.n_max =  4;
veh.xi_min = -pi/4;
veh.xi_max =  pi/4;
veh.n_min = -4;
veh.n_max =  4;

% inputs
veh.Tt_min =  0;
veh.Tt_max =  3500; % 2 wheels
veh.Tb_min = -1200*4; % 4 wheels
veh.Tb_max = 0;
veh.delta_min = -deg2rad(43);
veh.delta_max =  deg2rad(43);

veh.omega_rear_axle_min = max(3.69825, 0);
veh.omega_rear_axle_max = min(366.5, veh.V_max/veh.rw);

veh.omega_min = veh.omega_rear_axle_min;
veh.omega_max = veh.omega_rear_axle_max;

% rate of inputs
veh.Tt_dot_min = -4000;
veh.Tt_dot_max =  2000;
veh.Tb_dot_min = -1500*4;
veh.Tb_dot_max =  1500*4;

veh.delta_dot_min = -deg2rad(90);
veh.delta_dot_max =  deg2rad(90);

veh.omega_rear_axle_dot_min = (-((50-12)/veh.rw) / 6) * 1.5;
veh.omega_rear_axle_dot_max = (-1) * veh.omega_rear_axle_dot_min;

end
