function veh = my_params()

veh = struct();

veh.g = 9.80665;

veh.m = 2010;

veh.Ix = 624.4;

veh.Iy = 2521.3;

veh.Iz = 2742.2;

veh.Jzx = 0; % product of inertia Jzx [kg m^2]

veh.lf = 1.469;

veh.lr = 1.430;

veh.l = veh.lf + veh.lr;

veh.wt = 1.634; % distance between the center of FL and FR wheels (assuming front and the rear the same)

veh.t1 = veh.wt; % front track [m]

veh.t2 = veh.wt; % rear track [m]

veh.hc = 0.453;

veh.Iw = 15/2;

veh.rw = 0.346; % wheel radius

veh.ws = 0.275; % width of one tyre (assuming all four tyres the same width)

veh.h = veh.hc; % Chapter 9 CG height h [m]

veh.d = veh.h * 0.7; % VIP / equivalent no-roll-center height d [m]

veh.H = veh.h - veh.d; % H = h - d [m]

% q1 and q2 are the heights/moment arms used by the Chapter 9 lateral

% load-transfer equations. If only one no-roll-center height is available,

% q1=q2=d is the consistent two-axle first estimate.

veh.q1 = veh.d;

veh.q2 = veh.d;

rho = 1.206;

A = 2.8;

Cd = 0.3003;

Cl = 0.0; % assuming no downforce but enters VD EOM (due to uncertainty)

% negative means downforce in this sign convention

% veh.hcop = veh.hc;

% veh.acop = veh.lf;

veh.drag_coeff = rho*Cd*A;

veh.lift_coeff = rho*Cl*A;

veh.aero_front_frac = veh.lr/veh.l; % fraction of vertical aero force on front axle

veh.delay_Tdstate = 1e-1 * 0.5; % [s] 1e-1

veh.acceleration_delay = 1e-1; % [s]

Kw_FL = 30*1000;

Kw_FR = 30*1000;

Kw_RL = 43*1000;

Kw_RR = 43*1000;

Cw_FL = 1.01*1000;

Cw_FR = 1.01*1000;

Cw_RL = 1.83*1000;

Cw_RR = 1.83*1000;

% Make sure the anti-roll-bar stiffness is 0 for now, otherwise ODEs may

% have to change too?

veh.Karb_f = 0; % [Nm/(rad)]

veh.Karb_r = 0; % [Nm/(rad)]

% Pitch stiffness/damping. In Chapter 9 these can be obtained from the

% vertical/suspension stiffness and damping distribution; here they are

% explicit parameters for the OCP implementation.

veh.ktheta = (Kw_FL+Kw_FR)*(veh.lf)^2 + (Kw_RL+Kw_RR)*(veh.lr)^2; % Nm/rad

veh.ctheta = (Cw_FL+Cw_FR)*(veh.lf)^2 + (Cw_RL+Cw_RR)*(veh.lr)^2; % Nms/rad

% Roll stiffness/damping split by axle and total roll stiffness/damping.

veh.Kphi_f_spr = (Kw_FL+Kw_FR)*(veh.t1/2)^2;

veh.Kphi_r_spr = (Kw_RL+Kw_RR)*(veh.t2/2)^2;

veh.cphi1 = (Cw_FL+Cw_FR)*(veh.t1/2)^2;

veh.cphi2 = (Cw_RL+Cw_RR)*(veh.t2/2)^2;

veh.kphi1 = veh.Kphi_f_spr + veh.Karb_f;

veh.kphi2 = veh.Kphi_r_spr + veh.Karb_r;

veh.kphi = veh.kphi1 + veh.kphi2;

veh.cphi = veh.cphi1 + veh.cphi2;

% these parameters are for post-calculation of roll and pitch angles (outdated)

% initial condition of [roll, dot_roll, pitch, dot_pitch]

veh.roll_pitch_ode_initial_condition = [0 0 0 0]; % [rad, rad/s, rad, rad/s]

%% Brake and torque distribution

% veh.dr = 0.5;

veh.kb = 0.667;

veh.slip_ratio_value = 0.12; % this is the absolute value of the minimum negative slip ratio of each tyre targets

%% Tyre parameters

veh.Bx = 28.488818;

veh.Cx = 1.194995;

veh.d1x = 1.163672;

veh.d2x = 279.438470;

veh.Qx = 1.569014;

veh.By = 12.274751;

veh.Cy = 1.366372;

veh.d1y = 0.942175;

veh.d2y = 76.064551;

veh.Qy = 1.453161;

%% Don't change these for now

veh.Fz1 = 2000;

veh.Fz2 = 6000;

veh.eps1 = 1e-5;

veh.eps2 = 1e-5;

veh.my_effective_mu_scaling = 1.0; % 1.0 0.925 0.9 0.89 0.88 0.87 0.86 0.85 0.7 0.6 0.4

veh.rou_smooth = 1e-2;

veh.offset = veh.rou_smooth / pi;

veh.mux = 1.0; % this is for the old tyre model

veh.muy = 1.0; % this is for the old tyre model

%% Scaling factors

veh.V_s = 100;

veh.ux_s = veh.V_s;

veh.uy_s = veh.V_s;

veh.r_s = 1;

veh.s_s = 30;

veh.n_s = 5;

veh.xi_s = 1;

veh.Tt_s = 2000*2;

veh.Tb_s = 3000*2;

veh.delta_s = pi/8;

veh.Td_state_s = veh.Tt_s;

veh.theta_s = deg2rad(30);

veh.phi_s = deg2rad(30);

veh.thetadot_s = deg2rad(30);

veh.phidot_s = deg2rad(30);

veh.Fz_s = veh.m*veh.g/4;

% Small smooth numerical guards used only to avoid zero divisions in tyre-slip

% and curvilinear-coordinate denominators. Normal race solutions should stay

% far from these values.

veh.vx_eps = 1e-3;

veh.slip_eps = 1e-8;

%% State and control bounds

% states

veh.V_min = 10/3.6;

veh.V_max = 210/3.6; % [m/s]

veh.ux_min = veh.V_min; 

veh.ux_max = veh.V_max;

veh.uy_min = -veh.V_max;

veh.uy_max = veh.V_max;

% % veh.beta_min = -deg2rad(30);

% veh.beta_min = -deg2rad(50);

% veh.beta_max = (-1) * veh.beta_min;

veh.r_min = -pi/2;

veh.r_max = pi/2;

veh.n_min = -4;

veh.n_max = 4;

veh.xi_min = -pi/3; % xi is the angle between track reference line heading and vehicle's front

veh.xi_max = pi/3;

% veh.ax_bar_min = -1.5*veh.g;

% veh.ax_bar_max = 1.0*veh.g;

% veh.ay_bar_min = -1.1*veh.g;

% veh.ay_bar_max = (-1)*veh.ay_bar_min;

veh.theta_min = deg2rad(-30);

veh.theta_max = (-1) * veh.theta_min;

veh.phi_min = deg2rad(-30);

veh.phi_max = (-1) * veh.phi_min;

veh.thetadot_min = deg2rad(-30);

veh.thetadot_max = (-1) * veh.thetadot_min;

veh.phidot_min = deg2rad(-30);

veh.phidot_max = (-1) * veh.phidot_min;

% estimate gear ratio (engine_rad/s / wheel_rad/s)

veh.omega_motor_max = 18500 * (2*pi/60); % [rad/s]

% veh.gr = veh.omega_motor_max / (veh.V_max/veh.rw);

veh.gr = 10.49;

veh.power_motor_max = 250 * 10^3; % [Watts]

veh.tc_motor = 0.05; % [second]

% inputs

veh.Tt_min = 0;

veh.Tt_max = 373 * veh.gr;

% veh.Td_state_min = 0;

% veh.Td_state_max = veh.Tt_max;

veh.Tb_min = -2500*2 / (1-veh.kb);

veh.Tb_max = 0;

veh.delta_min = -deg2rad(43);

veh.delta_max = deg2rad(43);

% % omega_rear_axle_min and omega_rear_axle_max are determined from

% % power_engine_max VS omega_rear_axle curve

veh.omega_rear_axle_min = max(3.69827, veh.V_min/veh.rw); % [rad/s]

veh.omega_rear_axle_max = min(366.518, veh.V_max/veh.rw); % [rad/s]

veh.omega_min = veh.omega_rear_axle_min;

veh.omega_max = veh.omega_rear_axle_max;

veh.omega_s = veh.omega_max;

veh.omega_rear_axle_s = veh.omega_s;

% rate of inputs

veh.Tt_dot_min = -3000 * veh.gr;

veh.Tt_dot_max = 3000 * veh.gr; % [Nm/s]

veh.Tb_dot_min = -50000*2;

veh.Tb_dot_max = 50000*2;

veh.delta_dot_min = -deg2rad(90);

veh.delta_dot_max = deg2rad(90);

% omega_rear_axle_dot_min is estimated by the maximum deceleration of the

% car

estimated_seconds_from_max2min_speed = 2; % [s]

veh.omega_rear_axle_dot_min = (-1) * (veh.V_max - veh.V_min) / veh.rw / estimated_seconds_from_max2min_speed; % [rad/s^2]

veh.omega_rear_axle_dot_max = (-1) * veh.omega_rear_axle_dot_min;

veh.Fz_min = 50; % positive lower bound for friction/load equations [N]

veh.Fz_max = 4.0*veh.m*veh.g; % generous per-tyre upper bound [N]

%% Some small numerical values for protection

% Smoothly protect the tyre-slip divisions from zero wheel-frame forward speed.

% In normal race operation vx*_raw remains positive and this equals vx*_raw to

% numerical precision. The small guard is only for robustness.

veh.vx_eps = 1e-12;

end