function veh = veh_params()
%VEH_PARAMS Vehicle and solver scaling parameters for the Chapter-9 TOTPT model.
%
% This file keeps the original TOTPT race-car parameters where possible and
% adds the inertial, roll/pitch, no-roll-center, and tyre-load parameters
% required by the corrected Chapter 9 equations of motion and tyre-load
% equations.  The numerical roll/pitch parameters below are reasonable
% placeholders for testing the implementation; replace them with measured or
% identified data before using the model for conclusions.

veh = struct();

veh.mu = 1.0;
veh.g = 9.81;

%% Basic mass and geometry
veh.m = 1250;
veh.Iz = 1050;
veh.lf = 1.4;                  % a1: CG/VIP to front axle [m]
veh.lr = 1.4;                  % a2: CG/VIP to rear axle [m]
veh.l = veh.lf + veh.lr;
veh.hc = 0.35;                 % CG height h [m]
veh.Iw = 1.2;
veh.rw = 0.3;

% Front and rear tracks used by Chapter 9 tyre-load equations.
veh.wt = 1.5;                  % kept for compatibility with original plotting/code
veh.t1 = veh.wt;               % front track [m]
veh.t2 = veh.wt;               % rear track [m]
veh.ws = 0.2;                  % safety margin to track bounds [m]

%% Chapter 9 roll/pitch/yaw parameters
veh.Ix = 450;                  % roll inertia Jx [kg m^2]
veh.Iy = 1600;                 % pitch inertia Jy [kg m^2]
veh.Jzx = 0;                   % product of inertia Jzx [kg m^2]

veh.h = veh.hc;                % Chapter 9 CG height h [m]
veh.d = 0.10;                  % VIP / equivalent no-roll-center height d [m]
veh.H = veh.h - veh.d;         % H = h - d [m]

% q1 and q2 are the heights/moment arms used by the Chapter 9 lateral
% load-transfer equations.  If only one no-roll-center height is available,
% q1=q2=d is the consistent two-axle first estimate.
veh.q1 = veh.d;
veh.q2 = veh.d;

% Roll stiffness/damping split by axle and total roll stiffness/damping.
veh.kphi1 = 35000;             % front roll stiffness contribution [N m/rad]
veh.kphi2 = 45000;             % rear roll stiffness contribution [N m/rad]
veh.cphi1 = 2500;              % front roll damping contribution [N m s/rad]
veh.cphi2 = 3000;              % rear roll damping contribution [N m s/rad]
veh.kphi = veh.kphi1 + veh.kphi2;
veh.cphi = veh.cphi1 + veh.cphi2;

% Pitch stiffness/damping.  In Chapter 9 these can be obtained from the
% vertical/suspension stiffness and damping distribution; here they are
% explicit parameters for the OCP implementation.
veh.ktheta = 90000;            % pitch stiffness [N m/rad]
veh.ctheta = 6000;             % pitch damping [N m s/rad]

%% Aerodynamics
rho = 1.2;
A = 1.6;
Cd = 0.3;
Cl = -0.6;                     % negative means downforce in this sign convention
veh.drag_coeff = rho*Cd*A;
veh.lift_coeff = rho*Cl*A;
veh.aero_front_frac = veh.lr/veh.l;  % fraction of vertical aero force on front axle

%% Brake and torque distribution
veh.kb = 0.65;                 % front brake distribution fraction

%% Tire parameters
veh.Bx = 18;
veh.Cx = 1.3;
veh.d1x = 0.95;
veh.d2x = 320;

veh.By = 13;
veh.Cy = 1.5;
veh.d1y = 0.95;
veh.d2y = 320;

veh.mux = 1.0;
veh.muy = 1.0;

%% Scaling factors
veh.V_s = 100;
veh.ux_s = veh.V_s;
veh.uy_s = veh.V_s;
veh.r_s = 1;
veh.phi_s = 0.1;
veh.phidot_s = 1;
veh.theta_s = 0.05;
veh.thetadot_s = 1;
veh.n_s = 5;
veh.xi_s = 1;
veh.Tt_s = 1000*2;
veh.Tb_s = 1000*4;
veh.delta_s = pi/8;
veh.omega_s = veh.V_s/veh.rw;
veh.Fz_s = veh.m*veh.g/4;

% Small smooth numerical guards used only to avoid zero divisions in tyre-slip
% and curvilinear-coordinate denominators.  Normal race solutions should stay
% far from these values.
veh.vx_eps = 1e-3;
veh.slip_eps = 1e-8;

% Kept for backward compatibility with scripts using the original names.
veh.beta_s = 1;
veh.gamma_s = veh.r_s;
veh.ax_bar_s = veh.g;
veh.ay_bar_s = veh.g;
veh.s_s = 30;

%% State and input limits
veh.V_min = 0;
veh.V_max = 250/3.6;
veh.ux_min = 0.5;              % avoid tyre-slip singularity at vx ~= 0
veh.ux_max = veh.V_max;
veh.uy_min = -veh.V_max;
veh.uy_max =  veh.V_max;
veh.r_min = -pi/2;
veh.r_max =  pi/2;
veh.phi_min = -15*pi/180;
veh.phi_max =  15*pi/180;
veh.phidot_min = -2.5;
veh.phidot_max =  2.5;
veh.theta_min = -8*pi/180;
veh.theta_max =  8*pi/180;
veh.thetadot_min = -2.5;
veh.thetadot_max =  2.5;
veh.n_min = -4;
veh.n_max =  4;
veh.xi_min = -pi/4;
veh.xi_max =  pi/4;
veh.omega_min = 0;
veh.omega_max = veh.V_max/veh.rw;

veh.Fz_min = 50;               % positive lower bound for friction/load equations [N]
veh.Fz_max = 4.0*veh.m*veh.g;  % generous per-tyre upper bound [N]

% inputs
veh.Tt_min =  0;
veh.Tt_max =  1200*2;          % 2 driven rear wheels
veh.Tb_min = -1200*4;          % 4 braking wheels; negative as in original TOTPT
veh.Tb_max = 0;
veh.delta_min = -pi/4;
veh.delta_max =  pi/4;

% rate of inputs
veh.Tt_dot_min = -1500*2;
veh.Tt_dot_max =  1500*2;
veh.Tb_dot_min = -1500*4;
veh.Tb_dot_max =  1500*4;
veh.delta_dot_min = -22.5*pi/180;
veh.delta_dot_max =  22.5*pi/180;

end
