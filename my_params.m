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

%% Brake and torque distribution
% Tf = kt*Tt + kb*Tb
% Tr = (1-kt)*Tt + (%1*kb)*Tb
% kt = 0

veh.kb = 0.65;

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

veh.omega_min = 0;
veh.omega_max = veh.V_max/veh.rw;

% rate of inputs
veh.Tt_dot_min = -4000;
veh.Tt_dot_max =  2000;
veh.Tb_dot_min = -1500*4;
veh.Tb_dot_max =  1500*4;

veh.delta_dot_min = -deg2rad(90);
veh.delta_dot_max =  deg2rad(90);




end
