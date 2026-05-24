import casadi.*

veh = my_params();

%% Vehicle Model - state variables

nx = 11 + 4; % number of state variables

% longitudinal velocity [m/s]
V_n = SX.sym('v_n');
V_s = veh.V_s;
V = V_s * V_n;

% sideslip angle [rad]
beta_n = SX.sym('beta_n');
beta_s = veh.beta_s;
beta = beta_s * beta_n;

% yaw rate [rad/s]
gamma_n = SX.sym('gamma_n');
gamma_s = veh.gamma_s;
gamma = gamma_s * gamma_n;

% lateral distance to centreline [m] - left of centreline => n > 0; right => n < 0
n_n = SX.sym('n_n');
% n_s = 5;
n_s = veh.n_s;
n = n_s * n_n;

% course angle to centreline tangent direction [rad]
xi_n = SX.sym('xi_n');
xi_s = veh.xi_s;
xi = xi_s * xi_n;

ax_bar_n = SX.sym('ax_bar_n');
% ax_bar_s = 2*9.8;
ax_bar_s = veh.ax_bar_s;
ax_bar = ax_bar_s * ax_bar_n;

ay_bar_n = SX.sym('ay_bar_n');
% ay_bar_s = 2*9.8;
ay_bar_s = veh.ay_bar_s;
ay_bar = ay_bar_s * ay_bar_n;

omega_s = V_s/veh.rw;

% angular velocity front left tyre [rad/s]
omega_fl_n = SX.sym('omega_fl_n');
omega_fl = omega_s * omega_fl_n;

% angular velocity front right tyre [rad/s]
omega_fr_n = SX.sym('omega_fr_n');
omega_fr = omega_s * omega_fr_n;

% angular velocity rear left tyre [rad/s]
omega_rl_n = SX.sym('omega_rl_n');
omega_rl = omega_s * omega_rl_n;

% angular velocity rear right tyre [rad/s]
omega_rr_n = SX.sym('omega_rr_n');
omega_rr = omega_s * omega_rr_n;

% sprung mass pitch angle [rad]
theta_pitch_s = veh.theta_pitch_s;
theta_pitch_n = SX.sym('theta_pitch_n');
theta_pitch = theta_pitch_s * theta_pitch_n;

% sprung mass roll angle [rad]
phi_roll_s = veh.phi_roll_s;
phi_roll_n = SX.sym('phi_roll_n');
phi_roll = phi_roll_s * phi_roll_n;

% sprung mass pitch angle rate [rad/s]
theta_pitch_dot_s = veh.theta_pitch_dot_s;
theta_pitch_dot_n = SX.sym('theta_pitch_dot_n');
theta_pitch_dot = theta_pitch_dot_s * theta_pitch_dot_n;

% sprung mass roll angle rate [rad/s]
phi_roll_dot_s = veh.phi_roll_dot_s;
phi_roll_dot_n = SX.sym('phi_roll_dot_n');
phi_roll_dot = phi_roll_dot_s * phi_roll_dot_n;

% states vector (scaled)
x = [V_n beta_n gamma_n n_n xi_n ax_bar_n ay_bar_n omega_fl_n omega_fr_n omega_rl_n omega_rr_n theta_pitch_n phi_roll_n theta_pitch_dot_n phi_roll_dot_n]';

% scaling factors
x_s = [V_s beta_s gamma_s n_s xi_s ax_bar_s ay_bar_s omega_s omega_s omega_s omega_s theta_pitch_s phi_roll_s theta_pitch_dot_s phi_roll_dot_s]';

% state limits
x_min = [ -1e-3; veh.beta_min; veh.gamma_min; veh.n_min; veh.xi_min; veh.ax_bar_min; veh.ay_bar_min; veh.omega_min; veh.omega_min; veh.omega_min; veh.omega_min; veh.theta_pitch_min; veh.phi_roll_min; veh.theta_pitch_dot_min; veh.phi_roll_dot_min]./x_s;
x_max = [veh.V_max; veh.beta_max; veh.gamma_max; veh.n_max; veh.xi_max; veh.ax_bar_max; veh.ay_bar_max; veh.omega_max; veh.omega_max; veh.omega_max; veh.omega_max; veh.theta_pitch_max; veh.phi_roll_max; veh.theta_pitch_dot_max; veh.phi_roll_dot_max]./x_s;

%% Vehicle model - control variables (inputs)

nu = 3; % number of control variables

% driving torque (Nm)
Tt_n = SX.sym('Tt_n');
Tt_s = veh.Tt_s;
Tt = Tt_s * Tt_n;

% braking torque [Nm]
Tb_n = SX.sym('Tb_n');
Tb_s = veh.Tb_s;
Tb = Tb_s * Tb_n;

% steer angle [rad]
delta_n = SX.sym('delta_n');
delta_s = veh.delta_s;
delta = delta_s * delta_n;

% inputs vector (scaled)
u = [Tt_n Tb_n delta_n]';

% scaling factors for inputs
u_s = [Tt_s Tb_s delta_s]';

% input limits
u_min = [veh.Tt_min veh.Tb_min veh.delta_min]'./u_s;
u_max = [veh.Tt_max veh.Tb_max veh.delta_max]'./u_s;

%Constraints on the rate of inputs (units of each input per second)
duk_lb = [veh.Tt_dot_min veh.Tb_dot_min veh.delta_dot_min]'./u_s; %lower bound
duk_ub = [veh.Tt_dot_max veh.Tb_dot_max veh.delta_dot_max]'./u_s; %upper bound

% % Regularisation factors
rdu = [1; 1; 10]; %

%% Vehicle model - additional variables

nz = 0; % Number of aux variables

%% Vechicle model - parameter variables
%Symbolic variables that are not decision variables of the NLP - defined as such because their value changes (it is like a variable value parameter)

kappa = SX.sym('kappa'); % kappa > 0 for left turns

pv = kappa; % collect variable parameters

%% Vehicle model - equations
%Calculate variables of the system other than the states and inputs

cosd = cos(delta);
sind = sin(delta);
cosb = cos(beta);
sinb = sin(beta);

% longitudinal and lateral velocities
vx = V*cosb;
vy = V*sinb;

% Aero forces
Fdrag = 0.5*veh.drag_coeff*vx^2;
Fdown = 0.5*veh.downforce_coeff*vx^2;

% Tire velocities
vxfl = (vx-0.5*veh.wt*gamma)*cosd + (vy+veh.acom*gamma)*sind + 1e-12;

vxfr = (vx+0.5*veh.wt*gamma)*cosd + (vy+veh.acom*gamma)*sind + 1e-12;

vxrl = (vx-0.5*veh.wt*gamma) + 1e-12;
vxrr = (vx+0.5*veh.wt*gamma) + 1e-12;

vyfl = (vy+veh.acom*gamma)*cosd - (vx-0.5*veh.wt*gamma)*sind;
vyfr = (vy+veh.acom*gamma)*cosd - (vx+0.5*veh.wt*gamma)*sind;

vyrl = (vy-veh.bcom*gamma);
vyrr = (vy-veh.bcom*gamma);

% vertical tyre forces [N]
acom = veh.acom;
bcom = veh.bcom;
hcom = veh.hcom;

hcop = veh.hcop;
acop = veh.acop;

mass = veh.m;
gravity = veh.g;

dr = veh.dr;
twF_h = veh.wt/2;
twR_h = veh.wt/2;

% External pitch moment (front-unloading)
Mtheta_ext = -(Fdrag * hcop + Fdown * (acop - acom));

% Dynamic load
Mphi_sus2unspr_dyn = veh.Cphi * phi_roll_dot + veh.Kphi * phi_roll;
Mtheta_sus2unspr_dyn = veh.Ctheta * theta_pitch_dot + veh.Ktheta * theta_pitch;

dis_Kphi_ratio = veh.dis_Kphi_ratio;
dis_Cphi_dot_ratio = veh.dis_Cphi_dot_ratio;

Kphi_f = veh.Kphi * dis_Kphi_ratio;
Kphi_r = veh.Kphi * (1-dis_Kphi_ratio);

Cphi_f = veh.Cphi * dis_Cphi_dot_ratio;
Cphi_r = veh.Cphi * (1 - dis_Cphi_dot_ratio);

wt = veh.wt;

% Sprung/unsprung + roll-centre parameters for upgraded wheel-load model.
% IMPORTANT: existing parameters such as acom, bcom, hcom, hcop, acop, wt
% keep their original definitions. In particular, the existing aero pitch
% moment Mtheta_ext = -(Fdrag*hcop + Fdown*(acop - acom)) is not redefined.
%
% If the new fields are not yet added in my_params(), the fallback values
% reduce this model close to the previous lumped model:
%   m_uf = m_ur = 0, m_s = m, h_s = hcom, h_RC_f = h_RC_r = 0,
%   h_uF = h_uR = 0, wt_f = wt_r = wt, Ix_sprung = Ix, Iy_sprung = Iy.
m_uf = getvehfield_any(veh, {'m_uf','m_unsprung_f','muf','m_uns_f'}, 0);       % total front unsprung mass, both front wheels [kg]
m_ur = getvehfield_any(veh, {'m_ur','m_unsprung_r','mur','m_uns_r'}, 0);       % total rear unsprung mass, both rear wheels [kg]
m_s  = getvehfield_any(veh, {'m_s','m_sprung','ms'}, mass - m_uf - m_ur);     % sprung mass [kg]

h_s    = getvehfield_any(veh, {'h_s','h_sprung','hs'}, hcom);                % sprung-mass CG height above ground/reference plane [m]
h_RC_f = getvehfield_any(veh, {'h_RC_f','h_rc_f','hRCf'}, 0);                % front roll-centre height [m]
h_RC_r = getvehfield_any(veh, {'h_RC_r','h_rc_r','hRCr'}, 0);                % rear roll-centre height [m]
h_uF   = getvehfield_any(veh, {'h_uF','h_unsprung_f','huF'}, 0);             % front unsprung-mass CG height [m]
h_uR   = getvehfield_any(veh, {'h_uR','h_unsprung_r','huR'}, 0);             % rear unsprung-mass CG height [m]

tF = getvehfield_any(veh, {'wt_f','tF','track_f'}, wt);                      % front full track width [m]
tR = getvehfield_any(veh, {'wt_r','tR','track_r'}, wt);                      % rear full track width [m]

% Front/rear distribution of sprung-mass lateral inertia force.
% A simple default is static axle distribution: front = b/L, rear = a/L.
%
% IMPORTANT CONSISTENCY NOTE:
% The same lambdaF_s_lat/lambdaR_s_lat is used both for:
%   1) sprung geometric load transfer through front/rear roll centres, and
%   2) the effective roll-centre height seen by the sprung-mass lateral
%      inertia force in the roll equation.
% This keeps the steady-state identity
%   M_geom_F + M_geom_R + M_el_F + M_el_R = m_s*ay_bar*h_s
% when phi_dot = phi_ddot = 0, ignoring unsprung mass.
L = acom + bcom;
lambdaF_s_lat = getvehfield_any(veh, {'lambdaF_s_lat','lambda_f_s_lat'}, bcom/L);
lambdaR_s_lat = 1 - lambdaF_s_lat;

% Effective roll-centre height for the sprung-mass lateral inertia force.
% If lambdaF_s_lat = b/L, this is the usual geometric roll-axis height
% at the sprung-mass CG station. If lambdaF_s_lat is tuned to represent
% lateral force distribution, this force-weighted effective height keeps
% the wheel-load decomposition and roll equation mutually consistent.
h_RC_eff_s = lambdaF_s_lat*h_RC_f + lambdaR_s_lat*h_RC_r;
h_roll_s   = h_s - h_RC_eff_s;
Ix_sprung  = getvehfield_any(veh, {'Ix_sprung','Ix_s','Ixs'}, veh.Ix);
Iy_sprung  = getvehfield_any(veh, {'Iy_sprung','Iy_s','Iys'}, veh.Iy);

[Fz_FL, Fz_FR, Fz_RL, Fz_RR] = ...
wheel_load_new( ...
Mtheta_sus2unspr_dyn, Fdown, ax_bar, ...
bcom, mass, gravity, acom, ...
ay_bar, phi_roll, phi_roll_dot, ...
Kphi_f, Cphi_f, Kphi_r, Cphi_r, ...
m_uf, m_ur, m_s, ...
h_RC_f, h_RC_r, h_uF, h_uR, ...
lambdaF_s_lat, tF, tR);

fzfl = Fz_FL;
fzfr = Fz_FR;
fzrl = Fz_RL;
fzrr = Fz_RR;

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

alpha_fl = atan2(vyfl, vxfl);
alpha_fr = atan2(vyfr, vxfr);
alpha_rl = atan2(vyrl, vxrl);
alpha_rr = atan2(vyrr, vxrr);

[fxfl, fyfl, mufl_ratio] = tire_force_tanEllipse(veh, kappa_fl, alpha_fl, fzfl);
[fxfr, fyfr, mufr_ratio] = tire_force_tanEllipse(veh, kappa_fr, alpha_fr, fzfr);
[fxrl, fyrl, murl_ratio] = tire_force_tanEllipse(veh, kappa_rl, alpha_rl, fzrl);
[fxrr, fyrr, murr_ratio] = tire_force_tanEllipse(veh, kappa_rr, alpha_rr, fzrr);

% optional monitoring
% for my model
mufl = mufl_ratio;
mufr = mufr_ratio;
murl = murl_ratio;
murr = murr_ratio;

kb = veh.kb;

Tbf = kb * Tb;
Tbr = (1-kb) * Tb;

% user-chosen constants, tuned small enough
tuning_const_f = 0.5;
tuning_const_r = 0.5;

ay_bar_max_sol_predicted = veh.ay_bar_max;
Kf = tuning_const_f * 0.5/ay_bar_max_sol_predicted;
Kr = tuning_const_r * 0.5/ay_bar_max_sol_predicted;

zf = Kf * ay_bar;
zr = Kr * ay_bar;

Tfl = (0.5 - zf) * Tbf;
Tfr = (0.5 + zf) * Tbf;

Trl = 0.5*Tt + (0.5 - zr) * Tbr;
Trr = 0.5*Tt + (0.5 + zr) * Tbr;

%% Vehicle model - state derivatives
% Define derivatives of the state-space model
X = (fxfl+fxfr)*cosd - (fyfl+fyfr)*sind + (fxrl+fxrr) - Fdrag;

Y = (fxfl+fxfr)*sind + (fyfl+fyfr)*cosd + (fyrl+fyrr);

ax = X / veh.m;
ay = Y / veh.m;

Mz = (-1)*fxfl*cosd*(veh.wt/2) + fxfl*sind*veh.acom + fyfl*sind*(veh.wt/2) + fyfl*cosd*veh.acom +...
fxfr*cosd*(veh.wt/2) + fxfr*sind*veh.acom - fyfr*sind*(veh.wt/2) + fyfr*cosd*veh.acom +...
(-1)*fxrl*(veh.wt/2) - fyrl*veh.bcom +...
fxrr*(veh.wt/2) - fyrr*veh.bcom;

% Change of independent variable
chi = xi+beta;
sf = (1-n*kappa)/(V*cos(chi) + 1e-12) + 1e-12;

xdot = [

(X*cosb + Y*sinb) / veh.m;
(Y*cosb - X*sinb)/(V*veh.m) - gamma;

Mz/veh.Iz;
V*sin(chi);

gamma - kappa/sf;

(ax - ax_bar) / veh.acceleration_delay; % in this way ax will 'lead' ax_bar
(ay - ay_bar) / veh.acceleration_delay; % in this way ay will 'lead' ay_bar

(Tfl-fxfl*veh.rw)/veh.Iw;
(Tfr-fxfr*veh.rw)/veh.Iw;
(Trl-fxrl*veh.rw)/veh.Iw;
(Trr-fxrr*veh.rw)/veh.Iw;

theta_pitch_dot;
phi_roll_dot;

( -m_s*h_s*ax_bar + Mtheta_ext ...
- veh.Ctheta*theta_pitch_dot - veh.Ktheta*theta_pitch ) / Iy_sprung;

( m_s*h_roll_s*ay_bar ...
- veh.Cphi*phi_roll_dot - veh.Kphi*phi_roll ) / Ix_sprung;

]./x_s;

dx = sf*xdot;

%-Collect output variables
y_aero = [Fdrag; Fdown];
y_tyre = [fxfl; fxfr; fxrl; fxrr; fyfl; fyfr; fyrl; fyrr; fzfl; fzfr; fzrl; fzrr];
y_slip = [kappa_fl; kappa_fr; kappa_rl; kappa_rr; alpha_fl; alpha_fr; alpha_rl; alpha_rr];
y_acc = [ax; ay];
y_torque = [Tfl; Tfr; Trl; Trr];
y_mu = [mufl; mufr; murl; murr];
y_pwr = [Tfl*omega_fl; Tfr*omega_fr; Trl*omega_rl; Trr*omega_rr];
y_xdot = xdot;

%% Tyre model

% MY MODEL----------------------------------------------------------------
function [Fx, Fy, mu_ratio] = tire_force_tanEllipse(veh, kappa, alpha, Fz)

Qx = veh.Qx;
Qy = veh.Qy;
eps1 = veh.eps1;
eps2 = veh.eps2;

k_max = tan(pi/(2*veh.Cx))/veh.Bx;
alpha_peak = tan(pi/(2*veh.Cy))/veh.By;
tan_alpha_max = tan(alpha_peak);

[mu_x_max, mu_y_max] = muxmax_muymax_func(veh, Fz);

k_n = kappa / k_max;

a_n = (-tan(alpha)) / tan_alpha_max;

rou = (k_n.^2 + a_n.^2 + eps1).^0.5 + eps2;

Sx = pi/(2*atan2(Qx, 1));
Sy = pi/(2*atan2(Qy, 1));

my_effective_mu_scaling = veh.my_effective_mu_scaling;

Fx = (my_effective_mu_scaling * mu_x_max) .* Fz .* sin(Qx * atan2(Sx * rou, 1)) .* (k_n ./ rou);
Fy = (my_effective_mu_scaling * mu_y_max) .* Fz .* sin(Qy * atan2(Sy * rou, 1)) .* (a_n ./ rou);

% for debugging: should be <= 1 by construction
mu_ratio = hypot(Fx./((my_effective_mu_scaling * mu_x_max).*Fz), Fy./((my_effective_mu_scaling * mu_y_max).*Fz));

end

function [mu_x_max, mu_y_max] = muxmax_muymax_func(veh, Fz)

Fz1 = veh.Fz1;
Fz2 = veh.Fz2;
rou_smooth = veh.rou_smooth;

mu_x_max1 = (veh.d1x*Fz1 + veh.d2x)/Fz1; mu_x_max2 = (veh.d1x*Fz2 + veh.d2x)/Fz2;
mu_y_max1 = (veh.d1y*Fz1 + veh.d2y)/Fz1; mu_y_max2 = (veh.d1y*Fz2 + veh.d2y)/Fz2;

mu_x_max_raw = (Fz - Fz1) * (mu_x_max2 - mu_x_max1) / (Fz2 - Fz1) + mu_x_max1;
mu_y_max_raw = (Fz - Fz1) * (mu_y_max2 - mu_y_max1) / (Fz2 - Fz1) + mu_y_max1;

mu_x_max = mu_x_max_raw + smooth_ramp_func(mu_x_max_raw, rou_smooth);
mu_y_max = mu_y_max_raw + smooth_ramp_func(mu_y_max_raw, rou_smooth);

end

function smooth_ramp = smooth_ramp_func(x, rou_smooth)

smooth_ramp = (-x) .* ( atan2(-x, rou_smooth)/pi + 1/2 ) + rou_smooth / pi;

end

function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load_new( ...
Mtheta_sus2unspr, Fdown, ax_bar, ...
bcom, mass, gravity, acom, ...
ay_bar, phi_roll, phi_roll_dot, ...
Kphi_f, Cphi_f, Kphi_r, Cphi_r, ...
m_uf, m_ur, m_s, ...
h_RC_f, h_RC_r, h_uF, h_uR, ...
lambdaF_s_lat, tF, tR)

% Upgraded wheel-load model:
%   1) front/rear total normal load from an upgraded longitudinal
%      sprung/unsprung pitch-load model,
%   2) lateral load transfer split into:
%        - sprung geometric contribution through roll centres,
%        - unsprung direct contribution,
%        - elastic roll spring/damper contribution,
%   3) wheel loads using the physically consistent moment relation:
%        M_LT_axle = (Fz_right - Fz_left)*track/2.
%
% Sign convention kept consistent with the original code:
% positive longitudinal pitch-transfer moment increases front axle load;
% positive lateral transfer increases right-side load and decreases
% left-side load.

% 1. Front/rear total load distribution with upgraded longitudinal transfer.
% Keep existing parameter definitions. In particular, acop/hcop are not
% reinterpreted here; their pitch moment still enters via Mtheta_ext in the
% sprung pitch equation outside this function.
L = acom + bcom;

% Static/base vertical load.
%
% The existing acom/bcom are kept as the ORIGINAL TOTAL-VEHICLE CoM
% geometry. After splitting off front/rear unsprung masses, the sprung-mass
% static front fraction is therefore not generally bcom/L.  It must be
% chosen so that, when ax_bar=0, theta=0, theta_dot=0 and Fdown=0, the
% total vehicle still has the original static axle loads:
%   Fz_front = bcom*mass*gravity/L,
%   Fz_rear  = acom*mass*gravity/L.
%
% With front unsprung weight m_uf*gravity sitting directly on the front
% axle, this gives:
%   lambdaF_s_long*m_s*gravity + m_uf*gravity = bcom*mass*gravity/L.
% Hence:
lambdaF_s_long = (bcom*mass - L*m_uf) / (L*m_s);
lambdaR_s_long = 1 - lambdaF_s_long;
%
% Aero vertical force keeps the original front/rear reference split bcom/L
% here, while the existing aero pitch moment convention
%   Mtheta_ext = -(Fdrag*hcop + Fdown*(acop - acom))
% continues to enter through the pitch dynamics outside this function.
% Therefore acop/hcop are not reinterpreted in this wheel-load function.
Fz_front_base = lambdaF_s_long * m_s * gravity + m_uf * gravity + bcom * Fdown / L;
Fz_rear_base  = lambdaR_s_long * m_s * gravity + m_ur * gravity + acom * Fdown / L;

% Longitudinal load-transfer moment.
% Mtheta_sus2unspr is the elastic pitch spring/damper moment
%   Ktheta*theta + Ctheta*theta_dot.
% The unsprung longitudinal inertia contribution bypasses sprung pitch
% dynamics and directly transfers load between rear/front axles.
% With positive ax_bar (acceleration), this term is negative: front unloads,
% rear loads.
M_pitch_el  = Mtheta_sus2unspr;
M_pitch_uns = -(m_uf*h_uF + m_ur*h_uR) * ax_bar;
M_pitch_LT  = M_pitch_el + M_pitch_uns;

Fz_front_total = Fz_front_base + M_pitch_LT / L;
Fz_rear_total  = Fz_rear_base  - M_pitch_LT / L;

% 2. Lateral load-transfer moments at each axle.
lambdaR_s_lat = 1 - lambdaF_s_lat;

% Sprung-mass geometric load transfer through roll centres.
% Fy_sF and Fy_sR are approximated from the sprung-mass lateral inertia
% force distribution. This avoids an implicit tyre-force/Fz loop.
M_geom_F = lambdaF_s_lat * m_s * ay_bar * h_RC_f;
M_geom_R = lambdaR_s_lat * m_s * ay_bar * h_RC_r;

% Unsprung-mass direct lateral load transfer.
M_uns_F = m_uf * ay_bar * h_uF;
M_uns_R = m_ur * ay_bar * h_uR;

% Elastic roll contribution from suspension roll stiffness/damping.
M_el_F = Kphi_f * phi_roll + Cphi_f * phi_roll_dot;
M_el_R = Kphi_r * phi_roll + Cphi_r * phi_roll_dot;

% Total lateral load-transfer moments.
M_LT_F = M_geom_F + M_uns_F + M_el_F;
M_LT_R = M_geom_R + M_uns_R + M_el_R;

% 3. Final wheel loads.
% Because M_LT = (Fz_right - Fz_left)*track/2,
% the per-wheel increment is M_LT/track.
Fz_FL = 0.5 * Fz_front_total - M_LT_F / tF;
Fz_FR = 0.5 * Fz_front_total + M_LT_F / tF;
Fz_RL = 0.5 * Fz_rear_total  - M_LT_R / tR;
Fz_RR = 0.5 * Fz_rear_total  + M_LT_R / tR;

end

function val = getvehfield_any(veh, names, default_val)
%GETVEHFIELD_ANY Return the first existing field from a list of candidate names.
% This helper lets the upgraded model run even before all new fields are
% added to my_params(). Existing original parameters are not redefined.

val = default_val;
for ii = 1:numel(names)
    if isfield(veh, names{ii})
        val = veh.(names{ii});
        return;
    end
end

end
