function veh = my_params()
%MY_PARAMS Default parameter set for the reduced vehicle comparison scripts.
%
% This file is provided as a runnable fallback when your own my_params.m is
% not available on the MATLAB path.  The values are conservative demo values
% and are intended for code verification, not vehicle identification.
%
% The original MATLAB model uses these fields for scaling, limits, tyre
% forces, wheel-load transfer, pitch/roll dynamics, and input scaling.

veh = struct();

%% Basic vehicle constants
veh.m  = 1450.0;      % total mass [kg]
veh.g  = 9.81;        % gravity [m/s^2]
veh.rw = 0.31;        % tyre rolling radius [m]
veh.Iw = 1.25;        % wheel spin inertia [kg m^2]
veh.Ix = 650.0;       % roll inertia used by original reduced model [kg m^2]
veh.Iy = 1500.0;      % pitch inertia used by original reduced model [kg m^2]
veh.Iz = 2200.0;      % total yaw inertia [kg m^2]
veh.dr = 3.15;        % final-drive placeholder, not used by this comparison

%% Geometry
veh.acom = 1.35;      % CoM to front axle [m]
veh.bcom = 1.35;      % CoM to rear axle [m]
veh.wt   = 1.60;      % nominal full track width [m]
veh.tF   = veh.wt;    % front track [m]
veh.tR   = veh.wt;    % rear track [m]
veh.hcom = 0.55;      % CoM height [m]
veh.hcop = 0.30;      % aero drag CoP height [m]
veh.acop = 0.05;      % aero downforce CoP longitudinal location [m]

%% Sprung/unsprung and roll-centre load-transfer parameters
veh.m_uf = 84.0;                         % total front unsprung mass [kg]
veh.m_ur = 84.0;                         % total rear unsprung mass [kg]
veh.m_s  = veh.m - veh.m_uf - veh.m_ur;  % sprung mass [kg]
veh.h_s = veh.hcom;                      % sprung mass CG height [m]
veh.h_RC_f = 0.055;                      % front roll centre height [m]
veh.h_RC_r = 0.075;                      % rear roll centre height [m]
veh.h_uF = veh.rw;                       % front unsprung mass CG height [m]
veh.h_uR = veh.rw;                       % rear unsprung mass CG height [m]
veh.lambdaF_s_lat = veh.bcom/(veh.acom + veh.bcom);
veh.Ix_sprung = veh.Ix;
veh.Iy_sprung = veh.Iy;
veh.Izz_sprung = veh.Iz - veh.m_uf*(veh.acom^2 + (veh.tF/2)^2) - veh.m_ur*(veh.bcom^2 + (veh.tR/2)^2);

%% Pitch/roll reduced dynamics from CarSim-equivalent suspension/ARB data
% IMPORTANT:
%   - The four k_w/c_w values below must be WHEEL-CENTRE vertical rates.
%   - If CarSim gives spring/damper element rates, first convert them by
%       k_w = k_s * MR^2,  c_w = c_s * MR^2,
%     where MR = spring displacement / wheel vertical displacement.
%   - CarSim anti-roll-bar values are assumed here to be axle roll
%     stiffnesses directly in [N*m/rad], as requested.
%
% Replace these demo values by your CarSim linearized/equivalent values.
veh.k_w_FL = 11660.0;      % front-left wheel vertical stiffness [N/m]
veh.k_w_FR = 11660.0;      % front-right wheel vertical stiffness [N/m]
veh.k_w_RL = 11660.0;      % rear-left wheel vertical stiffness [N/m]
veh.k_w_RR = 11660.0;      % rear-right wheel vertical stiffness [N/m]

veh.c_w_FL = 1235.0;       % front-left wheel vertical damping [N*s/m]
veh.c_w_FR = 1235.0;       % front-right wheel vertical damping [N*s/m]
veh.c_w_RL = 1235.0;       % rear-left wheel vertical damping [N*s/m]
veh.c_w_RR = 1235.0;       % rear-right wheel vertical damping [N*s/m]

veh.Karb_F = 17576.0;      % front ARB roll stiffness from CarSim [N*m/rad]
veh.Karb_R = 17576.0;      % rear ARB roll stiffness from CarSim [N*m/rad]

% Build the reduced pitch/roll constants used by my_model_sprung_unsprung.m.
veh = set_pitch_roll_from_corner_rates_and_arb(veh);

veh.acceleration_delay = 0.10;

%% Aero
veh.rho = 1.225;
veh.CdA = 0.70;
veh.drag_coeff = veh.rho * veh.CdA;
veh.downforce_coeff = 1.80;

%% Tyre tan-ellipse parameters
veh.Qx = 1.35;
veh.Qy = 1.35;
veh.eps1 = 1.0e-8;
veh.eps2 = 1.0e-6;
veh.Bx = 12.0;
veh.Cx = 1.65;
veh.By = 11.0;
veh.Cy = 1.35;
veh.Fz1 = 1000.0;
veh.Fz2 = 6000.0;
veh.rou_smooth = 0.02;
veh.d1x = 1.10;
veh.d2x = 0.0;
veh.d1y = 1.15;
veh.d2y = 0.0;
veh.my_effective_mu_scaling = 1.0;

%% Input scaling, limits, and demo command magnitudes
veh.Tt_s = 1800.0;
veh.Tb_s = 4200.0;
veh.delta_s = 1.0;
veh.T_drive_max = veh.Tt_s;
veh.T_brake_max = veh.Tb_s;
veh.kb = 0.65;
veh.brake_bias_front = veh.kb;
veh.steer_max = deg2rad(8.0);
veh.Tt_min = 0.0;
veh.Tt_max = veh.Tt_s;
veh.Tb_min = -veh.Tb_s;
veh.Tb_max = 0.0;
veh.delta_min = -veh.steer_max;
veh.delta_max = veh.steer_max;
veh.Tt_dot_min = -5000.0;
veh.Tt_dot_max = 5000.0;
veh.Tb_dot_min = -10000.0;
veh.Tb_dot_max = 10000.0;
veh.delta_dot_min = -2.0;
veh.delta_dot_max = 2.0;

%% State scaling and bounds used by the original CasADi file
veh.V_s = 30.0;
veh.beta_s = 0.10;
veh.gamma_s = 0.50;
veh.n_s = 5.0;
veh.xi_s = 0.10;
veh.ax_bar_s = 2.0*veh.g;
veh.ay_bar_s = 2.0*veh.g;
veh.theta_pitch_s = 0.05;
veh.phi_roll_s = 0.05;
veh.theta_pitch_dot_s = 0.50;
veh.phi_roll_dot_s = 0.50;

veh.beta_min = -0.5;
veh.beta_max = 0.5;
veh.gamma_min = -2.0;
veh.gamma_max = 2.0;
veh.n_min = -10.0;
veh.n_max = 10.0;
veh.xi_min = -1.0;
veh.xi_max = 1.0;
veh.ax_bar_min = -3.0*veh.g;
veh.ax_bar_max = 3.0*veh.g;
veh.ay_bar_min = -3.0*veh.g;
veh.ay_bar_max = 3.0*veh.g;
veh.V_max = 80.0;
veh.omega_min = -10.0;
veh.omega_max = veh.V_max/veh.rw*3.0;
veh.theta_pitch_min = -0.3;
veh.theta_pitch_max = 0.3;
veh.phi_roll_min = -0.3;
veh.phi_roll_max = 0.3;
veh.theta_pitch_dot_min = -5.0;
veh.theta_pitch_dot_max = 5.0;
veh.phi_roll_dot_min = -5.0;
veh.phi_roll_dot_max = 5.0;
end

function veh = set_pitch_roll_from_corner_rates_and_arb(veh)
%SET_PITCH_ROLL_FROM_CORNER_RATES_AND_ARB
% Convert CarSim-equivalent corner wheel rates and front/rear ARB roll
% stiffnesses into the reduced model's lumped pitch/roll constants.
%
% Sign convention:
%   The scalar stiffnesses are positive restoring coefficients. The model
%   applies them as -Ktheta*theta and -Kphi*phi in the pitch/roll ODEs.
%
% Required fields:
%   k_w_FL,k_w_FR,k_w_RL,k_w_RR     [N/m]
%   c_w_FL,c_w_FR,c_w_RL,c_w_RR     [N*s/m]
%   Karb_F,Karb_R                   [N*m/rad]
%   acom,bcom,tF,tR                 [m]

kwFL = veh.k_w_FL;  kwFR = veh.k_w_FR;
kwRL = veh.k_w_RL;  kwRR = veh.k_w_RR;
cwFL = veh.c_w_FL;  cwFR = veh.c_w_FR;
cwRL = veh.c_w_RL;  cwRR = veh.c_w_RR;

a  = veh.acom;
b  = veh.bcom;
tF = veh.tF;
tR = veh.tR;

% Pitch: front and rear axle vertical motions from pitch.
veh.Ktheta = (kwFL + kwFR)*a^2 + (kwRL + kwRR)*b^2;
veh.Ctheta = (cwFL + cwFR)*a^2 + (cwRL + cwRR)*b^2;

% Roll from suspension vertical wheel rates.
veh.Kphi_F_spr = (kwFL + kwFR)*(tF/2)^2;
veh.Kphi_R_spr = (kwRL + kwRR)*(tR/2)^2;
veh.Cphi_F     = (cwFL + cwFR)*(tF/2)^2;
veh.Cphi_R     = (cwRL + cwRR)*(tR/2)^2;

% Roll from anti-roll bars. CarSim values are assumed [N*m/rad].
veh.Kphi_F_arb = veh.Karb_F;
veh.Kphi_R_arb = veh.Karb_R;

% Total front/rear roll stiffness used for elastic lateral load transfer.
veh.Kphi_F_total = veh.Kphi_F_spr + veh.Kphi_F_arb;
veh.Kphi_R_total = veh.Kphi_R_spr + veh.Kphi_R_arb;

% Lumped total coefficients used in the single roll DOF equation.
veh.Kphi = veh.Kphi_F_total + veh.Kphi_R_total;
veh.Cphi = veh.Cphi_F + veh.Cphi_R;

% Distribution ratios kept for backward compatibility with the existing model.
veh.dis_Kphi_ratio = veh.Kphi_F_total / max(veh.Kphi, eps);
veh.dis_Cphi_dot_ratio = veh.Cphi_F / max(veh.Cphi, eps);
end

