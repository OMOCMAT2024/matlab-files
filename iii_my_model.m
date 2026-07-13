import casadi.*

veh = my_params();

%% Vehicle Model - state variables
% Corrected Chapter 9 model with velocity at the vehicle invariant point M.
% States are scaled. Physical variables are recovered by multiplying by x_s.

nx = 13;

% VIP longitudinal velocity u [m/s]
ux_n = SX.sym('ux_n');
ux_s = veh.ux_s;
ux = ux_s * ux_n;

% VIP lateral velocity v [m/s]
uy_n = SX.sym('uy_n');
uy_s = veh.uy_s;
uy = uy_s * uy_n;

% yaw rate r [rad/s]
r_n = SX.sym('r_n');
r_s = veh.r_s;
r = r_s * r_n;

% roll angle phi [rad]
phi_n = SX.sym('phi_n');
phi_s = veh.phi_s;
phi = phi_s * phi_n;

% roll rate phidot [rad/s]
phidot_n = SX.sym('phidot_n');
phidot_s = veh.phidot_s;
phidot = phidot_s * phidot_n;

% pitch angle theta [rad]
theta_n = SX.sym('theta_n');
theta_s = veh.theta_s;
theta = theta_s * theta_n;

% pitch rate thetadot [rad/s]
thetadot_n = SX.sym('thetadot_n');
thetadot_s = veh.thetadot_s;
thetadot = thetadot_s * thetadot_n;

% lateral distance to centreline [m] - left of centreline => n > 0
n_n = SX.sym('n_n');
n_s = veh.n_s;
n = n_s * n_n;

% course angle to centreline tangent direction [rad]
xi_n = SX.sym('xi_n');
xi_s = veh.xi_s;
xi = xi_s * xi_n;

omega_s = veh.omega_s;

omega_fl_n = SX.sym('omega_fl_n');
omega_fr_n = SX.sym('omega_fr_n');
omega_rl_n = SX.sym('omega_rl_n');
omega_rr_n = SX.sym('omega_rr_n');
omega_fl = omega_s * omega_fl_n;
omega_fr = omega_s * omega_fr_n;
omega_rl = omega_s * omega_rl_n;
omega_rr = omega_s * omega_rr_n;

x = [ux_n; uy_n; r_n; phi_n; phidot_n; theta_n; thetadot_n; n_n; xi_n; ...
omega_fl_n; omega_fr_n; omega_rl_n; omega_rr_n];

x_s = [ux_s; uy_s; r_s; phi_s; phidot_s; theta_s; thetadot_s; n_s; xi_s; ...
omega_s; omega_s; omega_s; omega_s];

x_min = [veh.ux_min; veh.uy_min; veh.r_min; veh.phi_min; veh.phidot_min; ...
veh.theta_min; veh.thetadot_min; veh.n_min; veh.xi_min; ...
veh.omega_min; veh.omega_min; veh.omega_min; veh.omega_min] ./ x_s;
x_max = [veh.ux_max; veh.uy_max; veh.r_max; veh.phi_max; veh.phidot_max; ...
veh.theta_max; veh.thetadot_max; veh.n_max; veh.xi_max; ...
veh.omega_max; veh.omega_max; veh.omega_max; veh.omega_max] ./ x_s;

%% Vehicle model - control variables
nu = 3;

Tt_n = SX.sym('Tt_n');
Tt_s = veh.Tt_s;
Tt = Tt_s * Tt_n;

% Braking torque is non-positive as in the original convention.
Tb_n = SX.sym('Tb_n');
Tb_s = veh.Tb_s;
Tb = Tb_s * Tb_n;

delta_n = SX.sym('delta_n');
delta_s = veh.delta_s;
delta = delta_s * delta_n;

u = [Tt_n; Tb_n; delta_n];
u_s = [Tt_s; Tb_s; delta_s];

u_min = [veh.Tt_min; veh.Tb_min; veh.delta_min] ./ u_s;
u_max = [veh.Tt_max; veh.Tb_max; veh.delta_max] ./ u_s;

duk_lb = [veh.Tt_dot_min; veh.Tb_dot_min; veh.delta_dot_min] ./ u_s;
duk_ub = [veh.Tt_dot_max; veh.Tb_dot_max; veh.delta_dot_max] ./ u_s;

% Regularisation factors for knot-to-knot rates used in tro_code.m.
rdu = [1; 1; 10];
rdy = zeros(4,1);

%% Algebraic / auxiliary variables
% Lift only the four tyre vertical loads. This breaks the tyre-load algebraic
% loop without lifting ax, ay, r_dot, phi_ddot or theta_ddot.
nz = 4;

fzfl_n = SX.sym('fzfl_n');
fzfr_n = SX.sym('fzfr_n');
fzrl_n = SX.sym('fzrl_n');
fzrr_n = SX.sym('fzrr_n');

Fz_s = veh.Fz_s;
fzfl = Fz_s * fzfl_n;
fzfr = Fz_s * fzfr_n;
fzrl = Fz_s * fzrl_n;
fzrr = Fz_s * fzrr_n;

z = [fzfl_n; fzfr_n; fzrl_n; fzrr_n];
z_s = [Fz_s; Fz_s; Fz_s; Fz_s];
z_min = [veh.Fz_min; veh.Fz_min; veh.Fz_min; veh.Fz_min] ./ z_s;
z_max = [veh.Fz_max; veh.Fz_max; veh.Fz_max; veh.Fz_max] ./ z_s;

%% Vechicle model - parameter variables
%Symbolic variables that are not decision variables of the NLP - defined as such because their value changes (it is like a variable value parameter)

kappa = SX.sym('kappa'); % kappa > 0 for left turns
Cd = SX.sym('Cd');           % distance-varying aerodynamic drag coefficient
kers_on = SX.sym('kers_on'); % distance-varying KERS on/off state
kers_last_on = SX.sym('kers_last_on'); % final KERS interval indicator

pv = [kappa; Cd; kers_on; kers_last_on]; % collect variable parameters

%% Vehicle model - equations
%Calculate variables of the system other than the states and inputs

cosd = cos(delta);
sind = sin(delta);

% % longitudinal and lateral velocities
% vx = V*cosb;
% vy = V*sinb;

% Aerodynamics. With this sign convention, f_lift < 0 is downforce.
f_drag = 0.5 * veh.rho * Cd * veh.A * ux.^2;
f_lift = 0.5 * veh.lift_coeff * ux.^2;
Za1 = veh.aero_front_frac * f_lift;
Za2 = (1 - veh.aero_front_frac) * f_lift;

% Tire velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tyre contact-point velocities in wheel frames. Roll and pitch affect the
% Chapter 9 body EOM and vertical loads; they do not directly modify the tyre
% constitutive equations in this small-angle rigid-tyre model.
vx_eps = veh.vx_eps;

vxfl = (ux - 0.5*veh.t1*r)*cosd + (uy + veh.lf*r)*sind + vx_eps;
vxfr = (ux + 0.5*veh.t1*r)*cosd + (uy + veh.lf*r)*sind + vx_eps;
vxrl = ux - 0.5*veh.t2*r + vx_eps;
vxrr = ux + 0.5*veh.t2*r + vx_eps;

vyfl = (uy + veh.lf*r)*cosd - (ux - 0.5*veh.t1*r)*sind;
vyfr = (uy + veh.lf*r)*cosd - (ux + 0.5*veh.t1*r)*sind;
vyrl = uy - veh.lr*r;
vyrr = uy - veh.lr*r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tyre forces [N]

fzf = fzfl + fzfr;
fzr = fzrl + fzrr;

% Tyre longitudinal slips
kappa_fl = (veh.rw*omega_fl-vxfl)/vxfl;
kappa_fr = (veh.rw*omega_fr-vxfr)/vxfr;
kappa_rl = (veh.rw*omega_rl-vxrl)/vxrl;
kappa_rr = (veh.rw*omega_rr-vxrr)/vxrr;

% Tyre lateral slips
tan_alpha_fl = vyfl/vxfl;
tan_alpha_fr = vyfr/vxfr;
tan_alpha_rl = vyrl/vxrl;
tan_alpha_rr = vyrr/vxrr;

% alpha_fl = atan(tan_alpha_fl);
% alpha_fr = atan(tan_alpha_fr);
% alpha_rl = atan(tan_alpha_rl);
% alpha_rr = atan(tan_alpha_rr);
alpha_fl = atan2(vyfl, vxfl);
alpha_fr = atan2(vyfr, vxfr);
alpha_rl = atan2(vyrl, vxrl);
alpha_rr = atan2(vyrr, vxrr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fxfl, fyfl, mufl_ratio] = tire_force_tanEllipse(veh, kappa_fl, alpha_fl, fzfl);
[fxfr, fyfr, mufr_ratio] = tire_force_tanEllipse(veh, kappa_fr, alpha_fr, fzfr);
[fxrl, fyrl, murl_ratio] = tire_force_tanEllipse(veh, kappa_rl, alpha_rl, fzrl);
[fxrr, fyrr, murr_ratio] = tire_force_tanEllipse(veh, kappa_rr, alpha_rr, fzrr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% monitoring
mufl = mufl_ratio;
mufr = mufr_ratio;
murl = murl_ratio;
murr = murr_ratio;

%% Corrected Chapter 9 force and moment resultants
fxf = fxfl + fxfr;
fxr = fxrl + fxrr;
fyf = fyfl + fyfr;
fyr = fyrl + fyrr;

X1 = fxf*cosd - fyf*sind;
X2 = fxr;
Y1 = fyf*cosd + fxf*sind;
Y2 = fyr;

X = X1 + X2 - f_drag;
Y = Y1 + Y2;

N = veh.lf*Y1 - veh.lr*Y2 ...
+ 0.5*veh.t1*((fxfr - fxfl)*cosd + (fyfl - fyfr)*sind) ...
+ 0.5*veh.t2*(fxrr - fxrl);

H = veh.H;
LM = -veh.kphi*phi - veh.cphi*phidot + veh.m*veh.g*H*phi;
MM = -veh.ktheta*theta - veh.ctheta*thetadot;

% Corrected Chapter 9 compact EOM, solved explicitly after X,Y,N are known.
ax = X ./ veh.m;
ay = Y ./ veh.m;

theta_ddot = (MM - H*X) ./ veh.Iy;

Djr = veh.Ix*veh.Iz - veh.Jzx^2;
Aphi = LM + H*Y;
Ar = N - H*X*phi;
phi_ddot = (veh.Iz*Aphi + veh.Jzx*Ar) ./ Djr;
r_dot = (veh.Jzx*Aphi + veh.Ix*Ar) ./ Djr;

ux_dot = ax + uy*r - H*theta_ddot;
uy_dot = ay - ux*r + H*phi_ddot;

%% Corrected Chapter 9 tyre-load equations, enforced as path constraints
Z10 = veh.m*veh.g*veh.lr/veh.l;
Z20 = veh.m*veh.g*veh.lf/veh.l;

DeltaZ = (-X*veh.h + veh.Iy*theta_ddot) ./ veh.l;
Z1 = Z10 + DeltaZ - Za1;
Z2 = Z20 - DeltaZ - Za2;

DeltaZ1 = (Y1*veh.q1 + veh.kphi1*phi + veh.cphi1*phidot) ./ veh.t1;
DeltaZ2 = (Y2*veh.q2 + veh.kphi2*phi + veh.cphi2*phidot) ./ veh.t2;

fzfl_ch9 = 0.5*Z1 - DeltaZ1;
fzfr_ch9 = 0.5*Z1 + DeltaZ1;
fzrl_ch9 = 0.5*Z2 - DeltaZ2;
fzrr_ch9 = 0.5*Z2 + DeltaZ2;

load_res = [fzfl - fzfl_ch9;
fzfr - fzfr_ch9;
fzrl - fzrl_ch9;
fzrr - fzrr_ch9];

h_load = load_res ./ Fz_s;

%% Wheel torque

% dfzf = 0.5*(fzfr-fzfl);
% dfzr = 0.5*(fzrr-fzrl);

kb = veh.kb;
Tf = kb*Tb;
Tbr = (1-kb) * Tb;

% Tfl = 0.5*Tf*(1-dfzf/fzf);
% Tfr = 0.5*Tf*(1+dfzf/fzf);
Tfl = 0.5*Tf;
Tfr = 0.5*Tf;

% Trl = 0.5*Tt + 0.5*Tbr*(1-dfzr/fzr);
% Trr = 0.5*Tt + 0.5*Tbr*(1+dfzr/fzr);

kd = veh.kd;
Trl = (Tt - kd*(omega_rl-omega_rr))/2 + 0.5*Tbr;
Trr = (Tt + kd*(omega_rl-omega_rr))/2 + 0.5*Tbr;

% Tbf = kb * Tb;
% Tbr = (1-kb) * Tb;
%
% % user-chosen constants, tuned small enough
% tuning_const_f = 0.5;
% tuning_const_r = 0.5;
%
% ay_bar_max_sol_predicted = veh.ay_bar_max;
% Kf = tuning_const_f * 0.5/ay_bar_max_sol_predicted;
% Kr = tuning_const_r * 0.5/ay_bar_max_sol_predicted;
%
% zf = Kf * ay_bar;
% zr = Kr * ay_bar;
%
% Tfl = (0.5 - zf) * Tbf;
% Tfr = (0.5 + zf) * Tbf;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trl = 0.5*Tt + (0.5 - zr) * Tbr;
% Trr = 0.5*Tt + (0.5 + zr) * Tbr;
% % Trl = 0.5*Td_state + (0.5 - zr) * Tbr;
% % Trr = 0.5*Td_state + (0.5 + zr) * Tbr;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Vehicle model - state derivatives
% Define derivatives of the state-space model

% % Change of independent variable
% chi = xi+beta;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % sf = (1-n*kappa)/(V*cos(chi));
% sf = (1-n*kappa)/(V*cos(chi) + 1e-12) + 1e-12;
% % sf = 1/(V*cos(chi) + 1e-12) + 1e-12; % for straight test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Change of independent variable and state derivatives
% sf = dt/ds.
sf_den = ux*cos(xi) - uy*sin(xi) + vx_eps;
sf = (1 - n*kappa) ./ sf_den;

n_dot = ux*sin(xi) + uy*cos(xi);
xi_dot = r - kappa ./ sf;

%%
% xdot = [
%
% (X*cosb + Y*sinb) / veh.m;
% (Y*cosb - X*sinb)/(V*veh.m) - gamma;
%
% Mz/veh.Iz;
% V*sin(chi);
%
% gamma - kappa/sf;
%
% (ax - ax_bar) / veh.acceleration_delay; % in this way ax will 'lead' ax_bar
% (ay - ay_bar) / veh.acceleration_delay; % in this way ay will 'lead' ay_bar
%
% (Tfl-fxfl*veh.rw)/veh.Iw;
% (Tfr-fxfr*veh.rw)/veh.Iw;
% (Trl-fxrl*veh.rw)/veh.Iw;
% (Trr-fxrr*veh.rw)/veh.Iw;
%
% theta_pitch_dot;
% phi_roll_dot;
%
% % ( -veh.m*hcom*ax_bar + Mtheta_ext ...
% % - veh.Ctheta*theta_pitch_dot - veh.Ktheta*theta_pitch ) / veh.Iy;
% %
% % ( veh.m*hcom*ay_bar ...
% % - veh.Cphi*phi_roll_dot - veh.Kphi*phi_roll ) / veh.Ix;
%
% ( -m_s*h_s*ax_bar + Mtheta_ext ...
% - veh.Ctheta*theta_pitch_dot - veh.Ktheta*theta_pitch ) / Iy_sprung;
%
% ( m_s*h_roll_s*ay_bar ...
% - veh.Cphi*phi_roll_dot - veh.Kphi*phi_roll ) / Ix_sprung;
%
% ]./x_s;
%
% dx = sf*xdot;

xdot = [
ux_dot;
uy_dot;
r_dot;
phidot;
phi_ddot;
thetadot;
theta_ddot;
n_dot;
xi_dot;
(Tfl - fxfl*veh.rw) ./ veh.Iw;
(Tfr - fxfr*veh.rw) ./ veh.Iw;
(Trl - fxrl*veh.rw) ./ veh.Iw;
(Trr - fxrr*veh.rw) ./ veh.Iw;
] ./ x_s;

dx = sf .* xdot;

%-Collect output variables
y_aero = [f_drag; f_lift; Za1; Za2];
y_tyre = [fxfl; fxfr; fxrl; fxrr; fyfl; fyfr; fyrl; fyrr; fzfl; fzfr; fzrl; fzrr];
y_slip = [kappa_fl; kappa_fr; kappa_rl; kappa_rr; alpha_fl; alpha_fr; alpha_rl; alpha_rr];
y_acc = [ax; ay; ux_dot; uy_dot; r_dot; phi_ddot; theta_ddot];
y_torque = [Tfl; Tfr; Trl; Trr];
y_mu = [mufl; mufr; murl; murr];
y_pwr = [Tfl*omega_fl; Tfr*omega_fr; Trl*omega_rl; Trr*omega_rr];
y_xdot = xdot;
y_load = [h_load; fzfl_ch9; fzfr_ch9; fzrl_ch9; fzrr_ch9; DeltaZ; DeltaZ1; DeltaZ2];

%% Tyre model

% function [Fx, Fy, mu_ratio] = tire_force_tanEllipse(veh, kappa, tan_alpha, Fz)
function [Fx, Fy, mu_ratio] = tire_force_tanEllipse(veh, kappa, alpha, Fz)

Qx = veh.Qx;
Qy = veh.Qy;
eps1 = veh.eps1;
eps2 = veh.eps2;

% k_max = tan(pi/(2*veh.Cx))/veh.Bx;
% alpha_peak = tan(pi/(2*veh.Cy))/veh.By;
% tan_alpha_max = tan(alpha_peak);
k_max = 0.105;
alpha_peak = 0.148352986419518; % deg2rad(8.5)
tan_alpha_max = tan(alpha_peak);

[mu_x_max, mu_y_max] = muxmax_muymax_func(veh, Fz);

k_n = kappa / k_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a_n = (-tan_alpha) / tan_alpha_max;
a_n = (-tan(alpha)) / tan_alpha_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rou = (k_n.^2 + a_n.^2 + eps1).^0.5 + eps2;

% Sx = pi/(2*atan(Qx));
% Sy = pi/(2*atan(Qy));
Sx = pi/(2*atan2(Qx, 1));
Sy = pi/(2*atan2(Qy, 1));

my_effective_mu_scaling = veh.my_effective_mu_scaling;

% % Fx = (my_effective_mu_scaling * mu_x_max) .* Fz .* sin(Qx * atan(Sx * rou)) .* (k_n ./ rou);
% % Fy = (my_effective_mu_scaling * mu_y_max) .* Fz .* sin(Qy * atan(Sy * rou)) .* (a_n ./ rou);
% Fx = (1.04*my_effective_mu_scaling * mu_x_max) .* Fz .* sin(Qx * atan2(Sx * rou, 1)) .* (k_n ./ rou);
% Fy = (my_effective_mu_scaling * mu_y_max) .* Fz .* sin(Qy * atan2(Sy * rou, 1)) .* (a_n ./ rou);
Fx = (1.08*my_effective_mu_scaling * mu_x_max) .* Fz .* sin(Qx * atan2(Sx * rou, 1)) .* (k_n ./ rou);
Fy = (my_effective_mu_scaling * mu_y_max) .* Fz .* sin(Qy * atan2(Sy * rou, 1)) .* (a_n ./ rou);

% for debugging: should be <= 1 by construction
mu_ratio = hypot(Fx./((my_effective_mu_scaling * mu_x_max).*Fz), Fy./((my_effective_mu_scaling * mu_y_max).*Fz));

end

function [mu_x_max, mu_y_max] = muxmax_muymax_func(veh, Fz)

Fz1 = veh.Fz1;
Fz2 = veh.Fz2;
rou_smooth = veh.rou_smooth;

% mu_x_max1 = (veh.d1x*Fz1 + veh.d2x)/Fz1; mu_x_max2 = (veh.d1x*Fz2 + veh.d2x)/Fz2;
% mu_y_max1 = (veh.d1y*Fz1 + veh.d2y)/Fz1; mu_y_max2 = (veh.d1y*Fz2 + veh.d2y)/Fz2;
mu_x_max1 = 1.75; 
mu_x_max2 = 1.40;
mu_y_max1 = 1.80; 
mu_y_max2 = 1.45;

mu_x_max_raw = (Fz - Fz1) * (mu_x_max2 - mu_x_max1) / (Fz2 - Fz1) + mu_x_max1;
mu_y_max_raw = (Fz - Fz1) * (mu_y_max2 - mu_y_max1) / (Fz2 - Fz1) + mu_y_max1;

mu_x_max = mu_x_max_raw + smooth_ramp_func(mu_x_max_raw, rou_smooth);
mu_y_max = mu_y_max_raw + smooth_ramp_func(mu_y_max_raw, rou_smooth);

end

function smooth_ramp = smooth_ramp_func(x, rou_smooth)

smooth_ramp = (-x) .* ( atan2(-x, rou_smooth)/pi + 1/2 ) + rou_smooth / pi;

end