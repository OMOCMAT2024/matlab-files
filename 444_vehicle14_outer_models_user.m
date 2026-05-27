function [IN, aux] = vehicle14_outer_models_user(mode, t, x, S, par, args, road, meta)
%VEHICLE14_OUTER_MODELS_USER User-editable outer models for the 14-DOF enak/CasADi vehicle.
%
% This is the main file you should edit when you want to change:
%   1) control signals: total driving torque, total braking torque, steering;
%   2) powertrain model and torque split;
%   3) braking model and brake torque split;
%   4) tyre/contact model;
%   5) aero model.
%
% The generated file vehicle14_enak_matrices_casadi.m should remain untouched.
% This file converts your outer-skeleton models into the 34-input vector IN
% expected by the generated enak/CasADi equations of motion.
%
% Calling convention used by run_vehicle14_matlab_casadi_demo.m:
%   [IN0, ctrl]  = vehicle14_outer_models_user('control', t, x, [], par, args, road, meta)
%   [IN,  aux]   = vehicle14_outer_models_user('full',    t, x, S,  par, args, road, meta)
%
% The 'control' mode returns only delta and delta_dot in IN.  This is needed
% before the sensor evaluation because the front wheel-plate geometry depends
% on steering.  The 'full' mode returns all tyre, powertrain/brake, and aero
% inputs after the sensor vector S has been computed.

if nargin < 1 || isempty(mode)
    mode = 'full';
end
mode = lower(char(mode));

cfg = user_outer_model_options(par, args);
ctrl = user_control_signals(t, x, par, args, cfg);

IN = zeros(meta.nin,1);
IN(meta.input_index.delta) = ctrl.delta;
IN(meta.input_index.delta_dot) = ctrl.delta_dot;

aux = ctrl;
aux.cfg = cfg;
aux.contact = empty_contact_diag();

if strcmp(mode, 'control')
    return;
elseif ~strcmp(mode, 'full')
    error('Unknown vehicle14_outer_models_user mode: %s. Use control or full.', mode);
end

if isempty(S)
    error('Full outer-model mode requires the sensor vector S from the generated enak/CasADi model.');
end

qv = x(1:14);
uv = x(15:28);
in = meta.input_index;

% -------------------------------------------------------------------------
% 1. Powertrain and braking model: total torques -> four wheel torques.
% -------------------------------------------------------------------------
[wheel_torque, torque_diag] = user_powertrain_brake_model(ctrl, x, S, par, args, cfg, meta);
IN(in.TFL) = wheel_torque.FL;
IN(in.TFR) = wheel_torque.FR;
IN(in.TRL) = wheel_torque.RL;
IN(in.TRR) = wheel_torque.RR;

% -------------------------------------------------------------------------
% 2. Tyre/contact model: road contact -> body-frame wheel-centre wrenches.
% -------------------------------------------------------------------------
names = {'FL','FR','RL','RR'};
idx = meta.sensor_index;
contact = empty_contact_diag();
contact.names = names;

for ii = 1:4
    nm = names{ii};
    p_wc_N = [S(idx.([nm '_Wc_x'])); S(idx.([nm '_Wc_y'])); S(idx.([nm '_Wc_z']))];
    v_wc_N = [S(idx.([nm '_Wc_vx'])); S(idx.([nm '_Wc_vy'])); S(idx.([nm '_Wc_vz']))];
    head_N = [S(idx.([nm '_head_Nx'])); S(idx.([nm '_head_Ny'])); S(idx.([nm '_head_Nz']))];
    spin_N = [S(idx.([nm '_spin_Nx'])); S(idx.([nm '_spin_Ny'])); S(idx.([nm '_spin_Nz']))];
    carrier_omega_N = [S(idx.([nm '_carrier_wx'])); S(idx.([nm '_carrier_wy'])); S(idx.([nm '_carrier_wz']))];
    omega = uv(10+ii);

    ck = contact_kinematics_from_pose(par, road, p_wc_N, v_wc_N, head_N, spin_N, carrier_omega_N, omega);
    [F_B, M_B, d] = user_tyre_contact_wrench(qv, par, road, p_wc_N, v_wc_N, head_N, spin_N, omega, carrier_omega_N, ck, cfg);

    IN(in.(['FBx' nm])) = F_B(1);
    IN(in.(['FBy' nm])) = F_B(2);
    IN(in.(['FBz' nm])) = F_B(3);
    IN(in.(['MBx' nm])) = M_B(1);
    IN(in.(['MBy' nm])) = M_B(2);
    IN(in.(['MBz' nm])) = M_B(3);

    contact.Fx(ii,1) = d.Fx;
    contact.Fy(ii,1) = d.Fy;
    contact.Fz(ii,1) = d.Fz;
    contact.kappa(ii,1) = d.kappa;
    contact.alpha(ii,1) = d.alpha;
    contact.compression(ii,1) = d.compression;
    contact.road_z_at_wheel(ii,1) = d.road_patch_z;
    contact.contact_x(ii,1) = d.contact_x;
    contact.contact_y(ii,1) = d.contact_y;
    contact.contact_z(ii,1) = d.contact_z;
    contact.mu_ratio(ii,1) = d.mu_ratio;
end

% -------------------------------------------------------------------------
% 3. Aero model.
% -------------------------------------------------------------------------
[FdragF, FdownF, FdragR, FdownR, aero_diag] = user_aero_model(t, x, S, par, args, road, meta, cfg);
IN(in.FdragF) = FdragF;
IN(in.FdownF) = FdownF;
IN(in.FdragR) = FdragR;
IN(in.FdownR) = FdownR;

aux.contact = contact;
aux.torque = torque_diag;
aux.aero = aero_diag;
aux.wheel_torque = wheel_torque;
end

% =========================================================================
% USER-EDITABLE SECTION
% =========================================================================

function cfg = user_outer_model_options(par, args)
%USER_OUTER_MODEL_OPTIONS All editable model choices and constants in one place.
%
% This section is inspired by the user's older MATLAB model:
%   - define total Tt, total Tb, and delta as user signals;
%   - split brake torque with front bias kb and optional lateral redistribution;
%   - use a combined-slip tan-ellipse tyre model when tyre parameters exist;
%   - use smooth fallbacks so the demo still runs with the default vehicle14 parameters.

cfg.control_profile = getfield_default(args, 'maneuver_profile', 'aggressive-3d');
cfg.drop_settle_time = getfield_default(args, 'settle_time', 2.0);
cfg.maneuver_time = getfield_default(args, 'maneuver_time', 16.0);

cfg.T_drive_max = getfield_default(par, 'T_drive_max', 1900.0);     % total positive rear driving torque [Nm]
cfg.T_brake_max = getfield_default(par, 'T_brake_max', 4500.0);     % total positive braking magnitude [Nm]
cfg.steer_max = getfield_default(par, 'steer_max', deg2rad(7.0));

cfg.throttle_level = getfield_default(args, 'throttle', 0.65);
cfg.throttle_start = getfield_default(args, 'throttle_start', 0.20);
cfg.throttle_ramp = getfield_default(args, 'throttle_ramp', 0.45);
cfg.steer_deg = getfield_default(args, 'steer_deg', 7.0);
cfg.steer_start = getfield_default(args, 'steer_start', 1.0);
cfg.steer_ramp = getfield_default(args, 'steer_ramp', 0.5);
cfg.steer_freq = getfield_default(args, 'steer_freq', 0.55);
cfg.brake_start = getfield_default(args, 'brake_start', 11.5);
cfg.brake_level = getfield_default(args, 'brake_level', 0.28);
cfg.brake_ramp = getfield_default(args, 'brake_ramp', 0.35);
cfg.brake_hold = getfield_default(args, 'brake_hold', 1.15);

cfg.brake_bias_front = getfield_default(par, 'brake_bias_front', 0.65);
cfg.enable_lateral_brake_redistribution = true;
cfg.brake_lateral_tuning_f = 0.50;       % inspired by tuning_const_f in the old file
cfg.brake_lateral_tuning_r = 0.50;       % inspired by tuning_const_r in the old file
cfg.ay_bar_max_for_split = 2.0*getfield_default(par, 'g', 9.81);
cfg.max_left_right_brake_bias = 0.35;    % clamp zf/zr so each side remains positive

cfg.enable_rear_diff = true;
cfg.diff_gain = getfield_default(par, 'diff_gain', 35.0);
cfg.diff_width = getfield_default(par, 'diff_width', 15.0);

cfg.tyre_model = 'tan-ellipse';          % 'tan-ellipse' or 'linear-saturated'
cfg.my_effective_mu_scaling = getfield_default(par, 'my_effective_mu_scaling', 1.0);
cfg.v_eps = getfield_default(par, 'v_eps', 1.0);
cfg.kt = getfield_default(par, 'kt', 220000.0);
cfg.ct = getfield_default(par, 'ct', 1400.0);
cfg.normal_smooth_eps = getfield_default(par, 'normal_smooth_eps', 30.0);
cfg.mu = getfield_default(par, 'mu', 1.15);
cfg.Ck = getfield_default(par, 'Ck', 90000.0);
cfg.Ca = getfield_default(par, 'Ca', 85000.0);

% Tan-ellipse parameters.  If your my_params() has the old fields, these
% values can be replaced by par/veh fields.  The fallback values are moderate
% and solver-friendly for a runnable demonstration.
cfg.Bx = getfield_default(par, 'Bx', 10.0);
cfg.Cx = getfield_default(par, 'Cx', 1.65);
cfg.By = getfield_default(par, 'By', 7.5);
cfg.Cy = getfield_default(par, 'Cy', 1.30);
cfg.Qx = getfield_default(par, 'Qx', 1.60);
cfg.Qy = getfield_default(par, 'Qy', 1.50);
cfg.eps1 = getfield_default(par, 'eps1', 1.0e-8);
cfg.eps2 = getfield_default(par, 'eps2', 1.0e-6);
cfg.Fz1 = getfield_default(par, 'Fz1', 1000.0);
cfg.Fz2 = getfield_default(par, 'Fz2', 6000.0);
cfg.d1x = getfield_default(par, 'd1x', cfg.mu);
cfg.d2x = getfield_default(par, 'd2x', 0.0);
cfg.d1y = getfield_default(par, 'd1y', cfg.mu);
cfg.d2y = getfield_default(par, 'd2y', 0.0);
cfg.rou_smooth = getfield_default(par, 'rou_smooth', 1.0e-3);

cfg.disable_aero = getfield_default(args, 'disable_aero', false);
cfg.aero_h_min_multiplier = 0.20;
end

function ctrl = user_control_signals(t_abs, x, par, args, cfg)
%USER_CONTROL_SIGNALS Define total drive torque, total braking torque, and steering.
%
% This is the direct analogue of the old MATLAB variables:
%   Tt    = total driving torque, positive [Nm]
%   Tb    = total braking torque magnitude, positive here [Nm]
%   delta = steering angle [rad]
%
% Edit this function when you want to impose another maneuver, controller,
% MPC reference, or measured actuator trace.

if t_abs < cfg.drop_settle_time
    throttle = 0.0;
    brake = 0.0;
    delta = 0.0;
else
    tau = t_abs - cfg.drop_settle_time;
    profile = lower(cfg.control_profile);

    if strcmp(profile, 'straight-launch')
        throttle = cfg.throttle_level * smooth_step(tau, cfg.throttle_start, cfg.throttle_start + cfg.throttle_ramp);
        brake = 0.0;
        delta = 0.0;
    elseif strcmp(profile, 'creep-steer')
        throttle = cfg.throttle_level * smooth_step(tau, cfg.throttle_start, cfg.throttle_start + cfg.throttle_ramp);
        brake = 0.0;
        if cfg.brake_start >= 0.0 && tau >= cfg.brake_start
            bw = smooth_step(tau, cfg.brake_start, cfg.brake_start + cfg.brake_ramp);
            brake = cfg.brake_level*bw;
            throttle = throttle*(1.0 - bw);
        end
        steer_amp = deg2rad(cfg.steer_deg);
        env = smooth_step(tau, cfg.steer_start, cfg.steer_start + cfg.steer_ramp);
        delta = steer_amp*env*sin(2*pi*cfg.steer_freq*max(tau - cfg.steer_start, 0.0));
    elseif strcmp(profile, 'source-demo')
        throttle = 0.55*smooth_step(tau, 0.20, 0.80);
        brake = 0.0;
        delta = cfg.steer_max*smooth_step(tau, 0.80, 1.20)*sin(2*pi*0.35*max(tau - 0.80, 0.0));
    elseif strcmp(profile, 'aggressive-3d')
        launch = smooth_step(tau, 0.10, 0.55);
        throttle = cfg.throttle_level * launch;
        brake_window = 0.0;
        if cfg.brake_start >= 0.0
            brake_window = smooth_step(tau, cfg.brake_start, cfg.brake_start + cfg.brake_ramp) * ...
                (1.0 - smooth_step(tau, cfg.brake_start + cfg.brake_hold, cfg.brake_start + cfg.brake_hold + cfg.brake_ramp));
        end
        brake = cfg.brake_level * brake_window;
        throttle = throttle * (1.0 - 0.85*brake_window);
        steer_amp = deg2rad(cfg.steer_deg);
        env = smooth_step(tau, cfg.steer_start, cfg.steer_start + cfg.steer_ramp) * ...
            (1.0 - smooth_step(tau, cfg.maneuver_time - 1.0, cfg.maneuver_time - 0.2));
        tt = max(tau - cfg.steer_start, 0.0);
        raw = 0.78*sin(2*pi*cfg.steer_freq*tt) + 0.28*sin(2*pi*1.75*cfg.steer_freq*tt + 0.45);
        delta = steer_amp * env * min(max(raw, -1.0), 1.0);
    else
        error('Unknown control_profile/maneuver_profile: %s', profile);
    end
end

% Finite-difference steering rate.  This keeps the user signal definition
% simple and avoids adding delta as an independent state in this demo.
h = 1.0e-5;
delta_p = steering_only(t_abs + h, cfg);
delta_m = steering_only(t_abs - h, cfg);
delta_dot = (delta_p - delta_m)/(2*h);

ctrl.throttle = throttle;
ctrl.brake = brake;
ctrl.T_drive_total = cfg.T_drive_max * throttle;
ctrl.T_brake_total = cfg.T_brake_max * brake;
ctrl.delta = delta;
ctrl.delta_dot = delta_dot;
end

function delta = steering_only(t_abs, cfg)
%STEERING_ONLY Helper used only for delta_dot finite difference.
if t_abs < cfg.drop_settle_time
    delta = 0.0;
    return;
end
tau = t_abs - cfg.drop_settle_time;
profile = lower(cfg.control_profile);
if strcmp(profile, 'straight-launch')
    delta = 0.0;
elseif strcmp(profile, 'creep-steer')
    steer_amp = deg2rad(cfg.steer_deg);
    env = smooth_step(tau, cfg.steer_start, cfg.steer_start + cfg.steer_ramp);
    delta = steer_amp*env*sin(2*pi*cfg.steer_freq*max(tau - cfg.steer_start, 0.0));
elseif strcmp(profile, 'source-demo')
    delta = cfg.steer_max*smooth_step(tau, 0.80, 1.20)*sin(2*pi*0.35*max(tau - 0.80, 0.0));
elseif strcmp(profile, 'aggressive-3d')
    steer_amp = deg2rad(cfg.steer_deg);
    env = smooth_step(tau, cfg.steer_start, cfg.steer_start + cfg.steer_ramp) * ...
        (1.0 - smooth_step(tau, cfg.maneuver_time - 1.0, cfg.maneuver_time - 0.2));
    tt = max(tau - cfg.steer_start, 0.0);
    raw = 0.78*sin(2*pi*cfg.steer_freq*tt) + 0.28*sin(2*pi*1.75*cfg.steer_freq*tt + 0.45);
    delta = steer_amp * env * min(max(raw, -1.0), 1.0);
else
    error('Unknown control_profile/maneuver_profile: %s', profile);
end
end

function [wheel_torque, diag] = user_powertrain_brake_model(ctrl, x, S, par, args, cfg, meta)
%USER_POWERTRAIN_BRAKE_MODEL Split total Tt/Tb into wheel torques.
%
% Sign convention for wheel torques sent to the enak model:
%   positive torque drives wheel spin forward;
%   negative torque brakes the wheel.
%
% This follows the spirit of the old MATLAB code:
%   kb = front brake bias;
%   Tbf = kb*Tb; Tbr = (1-kb)*Tb;
%   optional left/right brake redistribution based on lateral acceleration.

uv = x(15:28);
Tt = ctrl.T_drive_total;   % positive drive torque magnitude
Tb = ctrl.T_brake_total;   % positive brake torque magnitude
kb = cfg.brake_bias_front;

% Rear-drive open-differential-like split with smooth speed-sensitive bias.
omRL = uv(13);
omRR = uv(14);
if cfg.enable_rear_diff
    diff_bias = cfg.diff_gain * cfg.diff_width * tanh((omRL - omRR)/cfg.diff_width);
else
    diff_bias = 0.0;
end
Tdrive_RL = 0.5*Tt - diff_bias;
Tdrive_RR = 0.5*Tt + diff_bias;

% Approximate lateral acceleration proxy for brake redistribution.  This is
% deliberately algebraic in the current state and therefore compatible with
% ode15i.  You can replace it with your own controller signal if desired.
vx = uv(1);
r = uv(6);
ay_proxy = vx*r;
if cfg.enable_lateral_brake_redistribution
    zf = cfg.brake_lateral_tuning_f * 0.5/cfg.ay_bar_max_for_split * ay_proxy;
    zr = cfg.brake_lateral_tuning_r * 0.5/cfg.ay_bar_max_for_split * ay_proxy;
    zf = min(max(zf, -cfg.max_left_right_brake_bias), cfg.max_left_right_brake_bias);
    zr = min(max(zr, -cfg.max_left_right_brake_bias), cfg.max_left_right_brake_bias);
else
    zf = 0.0;
    zr = 0.0;
end

Tbf = kb*Tb;
Tbr = (1.0-kb)*Tb;

wheel_torque.FL = -(0.5 - zf)*Tbf;
wheel_torque.FR = -(0.5 + zf)*Tbf;
wheel_torque.RL = Tdrive_RL - (0.5 - zr)*Tbr;
wheel_torque.RR = Tdrive_RR - (0.5 + zr)*Tbr;

diag.T_drive_total = Tt;
diag.T_brake_total = Tb;
diag.brake_bias_front = kb;
diag.zf = zf;
diag.zr = zr;
diag.ay_proxy = ay_proxy;
diag.diff_bias = diff_bias;
end

function [F_B, M_B, diag] = user_tyre_contact_wrench(qv, par, road, p_wc_N, v_wc_N, wheel_heading_N, wheel_spin_axis_N, omega, carrier_omega_N, ck, cfg)
%USER_TYRE_CONTACT_WRENCH 3D road/contact/tyre model external to enak.
R_NB = body_rotation_from_q(qv);
R_BN = R_NB.';
geom = oriented_tyre_contact_geometry(par, road, p_wc_N, wheel_heading_N, wheel_spin_axis_N);

n_road_N = geom.n_road_N;
t_long_N = geom.t_long_N;
t_lat_N = geom.t_lat_N;
r_wc_to_contact_N = geom.r_wc_to_contact_N;
compression = geom.compression;
compression_rate = ck.compression_rate_normal;
compression_rate_radial = ck.compression_rate_radial;
v_contact_geom_N = ck.contact_point_velocity_geom_N;
v_contact_carrier_N = ck.contact_patch_velocity_carrier_N;
v_spin_contact_N = ck.contact_patch_velocity_spin_N;
v_contact_material_N = ck.contact_patch_velocity_material_N;

Fz = smoothplus(cfg.kt*compression + cfg.ct*compression_rate, cfg.normal_smooth_eps);
Vx = dot(v_contact_carrier_N, t_long_N);
Vy = dot(v_contact_carrier_N, t_lat_N);
Vx_material = dot(v_contact_material_N, t_long_N);
Vy_material = dot(v_contact_material_N, t_lat_N);
rolling_speed = -dot(v_spin_contact_N, t_long_N);
Vden = sqrt(Vx*Vx + cfg.v_eps*cfg.v_eps);
kappa = (rolling_speed - Vx)/Vden;
alpha = atan2(Vy, Vden);

if strcmpi(cfg.tyre_model, 'tan-ellipse')
    [Fx, Fy, mu_ratio] = tyre_force_tan_ellipse(cfg, kappa, alpha, Fz);
else
    [Fx, Fy, mu_ratio] = tyre_force_linear_saturated(cfg, kappa, alpha, Fz);
end

F_contact_N = Fx*t_long_N + Fy*t_lat_N + Fz*n_road_N;
F_B = R_BN*F_contact_N;
r_contact_B = R_BN*r_wc_to_contact_N;
M_B = cross(r_contact_B, F_B);

diag.s = geom.s;
diag.n = geom.n;
diag.compression = compression;
diag.compression_rate = compression_rate;
diag.compression_rate_radial = compression_rate_radial;
diag.Fz = Fz;
diag.Fx = Fx;
diag.Fy = Fy;
diag.mu_ratio = mu_ratio;
diag.kappa = kappa;
diag.alpha = alpha;
diag.Vx_contact_carrier = Vx;
diag.Vy_contact_carrier = Vy;
diag.Vx_contact_material = Vx_material;
diag.Vy_contact_material = Vy_material;
diag.rolling_speed = rolling_speed;
diag.contact_x = geom.p_contact_N(1);
diag.contact_y = geom.p_contact_N(2);
diag.contact_z = geom.p_contact_N(3);
diag.road_patch_x = geom.p_road_N(1);
diag.road_patch_y = geom.p_road_N(2);
diag.road_patch_z = geom.p_road_N(3);
diag.road_normal_velocity_residual = ck.road_normal_velocity_residual;
diag.road_plane_residual_rate = ck.road_plane_residual_rate;
diag.contact_vx = v_contact_geom_N(1);
diag.contact_vy = v_contact_geom_N(2);
diag.contact_vz = v_contact_geom_N(3);
end

function [Fx, Fy, mu_ratio] = tyre_force_tan_ellipse(cfg, kappa, alpha, Fz)
%TYRE_FORCE_TAN_ELLIPSE Combined-slip tyre inspired by the uploaded old MATLAB model.
Qx = cfg.Qx;
Qy = cfg.Qy;
eps1 = cfg.eps1;
eps2 = cfg.eps2;

k_max = tan(pi/(2*cfg.Cx))/cfg.Bx;
alpha_peak = tan(pi/(2*cfg.Cy))/cfg.By;
tan_alpha_max = tan(alpha_peak);

[mu_x_max, mu_y_max] = muxmax_muymax_func(cfg, Fz);

k_n = kappa / k_max;
a_n = (-tan(alpha)) / tan_alpha_max;
rou = sqrt(k_n.^2 + a_n.^2 + eps1) + eps2;

Sx = pi/(2*atan2(Qx, 1));
Sy = pi/(2*atan2(Qy, 1));
mu_scale = cfg.my_effective_mu_scaling;

Fx = (mu_scale * mu_x_max) .* Fz .* sin(Qx * atan2(Sx * rou, 1)) .* (k_n ./ rou);
Fy = (mu_scale * mu_y_max) .* Fz .* sin(Qy * atan2(Sy * rou, 1)) .* (a_n ./ rou);

denx = max(abs((mu_scale * mu_x_max).*Fz), 1.0e-9);
deny = max(abs((mu_scale * mu_y_max).*Fz), 1.0e-9);
mu_ratio = hypot(Fx./denx, Fy./deny);
end

function [Fx, Fy, mu_ratio] = tyre_force_linear_saturated(cfg, kappa, alpha, Fz)
%TYRE_FORCE_LINEAR_SATURATED Fallback from the previous executable demo.
Fx0 = cfg.Ck*kappa;
Fy0 = -cfg.Ca*alpha;
limit = cfg.mu*Fz + 1.0e-9;
norm_tan = sqrt(Fx0*Fx0 + Fy0*Fy0 + 1.0e-12);
scale = (1.0 + (norm_tan/limit)^4)^(-0.25);
Fx = Fx0*scale;
Fy = Fy0*scale;
mu_ratio = norm_tan/max(limit, 1.0e-9);
end

function [mu_x_max, mu_y_max] = muxmax_muymax_func(cfg, Fz)
%MUXMAX_MUYMAX_FUNC Load-sensitive friction estimate from the old MATLAB model.
Fz1 = cfg.Fz1;
Fz2 = cfg.Fz2;
rou_smooth = cfg.rou_smooth;

mu_x_max1 = (cfg.d1x*Fz1 + cfg.d2x)/Fz1;
mu_x_max2 = (cfg.d1x*Fz2 + cfg.d2x)/Fz2;
mu_y_max1 = (cfg.d1y*Fz1 + cfg.d2y)/Fz1;
mu_y_max2 = (cfg.d1y*Fz2 + cfg.d2y)/Fz2;

mu_x_max_raw = (Fz - Fz1) * (mu_x_max2 - mu_x_max1) / (Fz2 - Fz1) + mu_x_max1;
mu_y_max_raw = (Fz - Fz1) * (mu_y_max2 - mu_y_max1) / (Fz2 - Fz1) + mu_y_max1;

mu_x_max = mu_x_max_raw + smooth_ramp_func(mu_x_max_raw, rou_smooth);
mu_y_max = mu_y_max_raw + smooth_ramp_func(mu_y_max_raw, rou_smooth);
end

function smooth_ramp = smooth_ramp_func(x, rou_smooth)
smooth_ramp = (-x) .* ( atan2(-x, rou_smooth)/pi + 1/2 ) + rou_smooth / pi;
end

function [FdragF, FdownF, FdragR, FdownR, diag] = user_aero_model(t_abs, x, S, par, args, road, meta, cfg)
%USER_AERO_MODEL Drag/downforce model.  Edit here for aero maps.
qv = x(1:14);
uv = x(15:28);

if t_abs < cfg.drop_settle_time || cfg.disable_aero
    FdragF = 0.0;
    FdownF = 0.0;
    FdragR = 0.0;
    FdownR = 0.0;
    diag.enabled = false;
    return;
end

vx = uv(1);
qdyn = 0.5*getfield_default(par, 'rho', 1.225)*vx*vx;
[p_af_N,~] = body_point_state_num(qv, uv, [getfield_default(par, 'aero_front_x', 1.25);0;getfield_default(par, 'aero_z', 0.05)]);
[p_ar_N,~] = body_point_state_num(qv, uv, [getfield_default(par, 'aero_rear_x', -1.20);0;getfield_default(par, 'aero_z', 0.05)]);
[hF,~,~,~] = road_height_normal_at_point(road, p_af_N);
[hR,~,~,~] = road_height_normal_at_point(road, p_ar_N);

multF_raw = 1.0 + getfield_default(par, 'aero_h_sens_front', 1.0)*(getfield_default(par, 'aero_h_ref', 0.080) - hF);
multR_raw = 1.0 + getfield_default(par, 'aero_h_sens_rear', 1.2)*(getfield_default(par, 'aero_h_ref', 0.080) - hR);
multF = smooth_max(cfg.aero_h_min_multiplier, multF_raw, 1.0e-3);
multR = smooth_max(cfg.aero_h_min_multiplier, multR_raw, 1.0e-3);

D = 0.5*getfield_default(par, 'rho', 1.225)*getfield_default(par, 'CdA', 0.65)*vx*smooth_abs(vx, 0.25);
FdragF = -0.5*D;
FdragR = -0.5*D;
FdownF = -qdyn*getfield_default(par, 'ClA_front', 0.48)*multF;
FdownR = -qdyn*getfield_default(par, 'ClA_rear', 0.72)*multR;

diag.enabled = true;
diag.hF = hF;
diag.hR = hR;
diag.multF = multF;
diag.multR = multR;
diag.drag_total = D;
end

% =========================================================================
% GEOMETRY, CONTACT, ROAD, AND MATH HELPERS
% =========================================================================

function ck = contact_kinematics_from_pose(par, road, p_wc_N, v_wc_N, head_N, spin_N, carrier_omega_N, omega)
head_N = normalize_num(head_N, [1;0;0]);
spin_N = normalize_num(spin_N, [0;1;0]);
geom0 = oriented_tyre_contact_geometry(par, road, p_wc_N, head_N, spin_N);
head_dot_N = cross(carrier_omega_N, head_N);
spin_dot_N = cross(carrier_omega_N, spin_N);
h = 1.0e-5;
geomp = oriented_tyre_contact_geometry(par, road, p_wc_N + h*v_wc_N, head_N + h*head_dot_N, spin_N + h*spin_dot_N);
geomm = oriented_tyre_contact_geometry(par, road, p_wc_N - h*v_wc_N, head_N - h*head_dot_N, spin_N - h*spin_dot_N);
ck.geom = geom0;
ck.contact_point_velocity_geom_N = (geomp.p_contact_N - geomm.p_contact_N)/(2*h);
ck.contact_unloaded_velocity_geom_N = (geomp.p_contact_unloaded_N - geomm.p_contact_unloaded_N)/(2*h);
ck.radial_down_rate_N = (geomp.radial_down_N - geomm.radial_down_N)/(2*h);
ck.road_normal_rate_N = (geomp.n_road_N - geomm.n_road_N)/(2*h);
ck.t_long_rate_N = (geomp.t_long_N - geomm.t_long_N)/(2*h);
ck.t_lat_rate_N = (geomp.t_lat_N - geomm.t_lat_N)/(2*h);
ck.compression_rate_normal = (geomp.compression_normal - geomm.compression_normal)/(2*h);
ck.compression_rate_radial = (geomp.compression_radial - geomm.compression_radial)/(2*h);
ck.contact_gap_rate_normal = (geomp.contact_gap_normal - geomm.contact_gap_normal)/(2*h);
ck.contact_gap_rate_radial = (geomp.contact_gap_radial - geomm.contact_gap_radial)/(2*h);
ck.road_plane_residual_rate = (geomp.road_plane_residual - geomm.road_plane_residual)/(2*h);
r_wc_to_contact_N = geom0.r_wc_to_contact_N;
spin_axis_N = geom0.spin_axis_N;
ck.contact_patch_velocity_carrier_N = v_wc_N + cross(carrier_omega_N, r_wc_to_contact_N);
ck.contact_patch_velocity_spin_N = cross(omega*spin_axis_N, r_wc_to_contact_N);
ck.contact_patch_velocity_material_N = ck.contact_patch_velocity_carrier_N + ck.contact_patch_velocity_spin_N;
ck.carrier_omega_N = carrier_omega_N;
ck.road_normal_velocity_residual = dot(ck.contact_point_velocity_geom_N, geom0.n_road_N);
end

function geom = oriented_tyre_contact_geometry(par, road, p_wc_N, wheel_heading_N, wheel_spin_axis_N)
[s,n,p_road_N,e_long_N,e_lat_N,n_road_N] = road_project_point(road, p_wc_N);
wheel_heading_N = normalize_num(wheel_heading_N, e_long_N);
spin_axis_N = normalize_num(wheel_spin_axis_N, e_lat_N);
t_long_raw = cross(spin_axis_N, n_road_N);
if norm(t_long_raw) < 1.0e-10
    t_long_raw = wheel_heading_N - dot(wheel_heading_N,n_road_N)*n_road_N;
end
t_long_N = normalize_num(t_long_raw, e_long_N);
if dot(t_long_N, wheel_heading_N) < 0.0
    t_long_N = -t_long_N;
end
t_lat_N = normalize_num(cross(n_road_N, t_long_N), e_lat_N);
normal_in_radial_plane = n_road_N - dot(n_road_N, spin_axis_N)*spin_axis_N;
gamma = norm(normal_in_radial_plane);
if gamma < 1.0e-10
    gamma = 1.0;
    radial_down_N = -n_road_N;
    degenerate = true;
else
    radial_down_N = -normal_in_radial_plane/gamma;
    degenerate = false;
end
contact_to_wc_unit_N = -radial_down_N;
r_wc_to_unloaded_N = par.rw*radial_down_N;
p_contact_unloaded_N = p_wc_N + r_wc_to_unloaded_N;
signed_center_distance = dot(p_wc_N - p_road_N, n_road_N);
radial_depth_normal = par.rw*gamma;
compression_normal = radial_depth_normal - signed_center_distance;
distance_center_to_road_along_radial = signed_center_distance/gamma;
compression_radial = par.rw - distance_center_to_road_along_radial;
r_wc_to_contact_N = distance_center_to_road_along_radial*radial_down_N;
p_contact_N = p_wc_N + r_wc_to_contact_N;
road_plane_residual = dot(p_contact_N - p_road_N, n_road_N);
geom.s = s;
geom.n = n;
geom.p_road_N = p_road_N;
geom.e_long_N = e_long_N;
geom.e_lat_N = e_lat_N;
geom.n_road_N = n_road_N;
geom.t_long_N = t_long_N;
geom.t_lat_N = t_lat_N;
geom.spin_axis_N = spin_axis_N;
geom.radial_down_N = radial_down_N;
geom.contact_to_wc_unit_N = contact_to_wc_unit_N;
geom.r_wc_to_contact_N = r_wc_to_contact_N;
geom.p_contact_N = p_contact_N;
geom.p_contact_road_N = p_contact_N;
geom.r_wc_to_unloaded_N = r_wc_to_unloaded_N;
geom.p_contact_unloaded_N = p_contact_unloaded_N;
geom.signed_center_distance = signed_center_distance;
geom.radial_projection_norm = gamma;
geom.radial_depth = radial_depth_normal;
geom.radial_depth_normal = radial_depth_normal;
geom.distance_center_to_road_along_radial = distance_center_to_road_along_radial;
geom.compression = compression_normal;
geom.compression_normal = compression_normal;
geom.compression_radial = compression_radial;
geom.contact_gap_normal = -compression_normal;
geom.contact_gap_radial = -compression_radial;
geom.road_plane_residual = road_plane_residual;
geom.degenerate_contact = degenerate;
end

function [h,e_z,s,n] = road_height_normal_at_point(road,p_N)
[s,n,p_road,~,~,e_z] = road_project_point(road,p_N);
h = dot(p_N-p_road,e_z);
end

function [pN,vN] = body_point_state_num(qv, uv, r_B)
R = body_rotation_from_q(qv);
pP = qv(1:3);
vP = R*uv(1:3);
omega = R*uv(4:6);
rN = R*r_B;
pN = pP + rN;
vN = vP + cross(omega,rN);
end

function R = body_rotation_from_q(qv)
R = rot_z(qv(4))*rot_y(qv(5))*rot_x(qv(6));
end

function R = rot_x(a)
ca = cos(a);
sa = sin(a);
R = [1 0 0; 0 ca -sa; 0 sa ca];
end

function R = rot_y(a)
ca = cos(a);
sa = sin(a);
R = [ca 0 sa; 0 1 0; -sa 0 ca];
end

function R = rot_z(a)
ca = cos(a);
sa = sin(a);
R = [ca -sa 0; sa ca 0; 0 0 1];
end

function [s,n,p_road,e_long,e_lat,e_z] = road_project_point(road,p)
s = p(1);
for kk = 1:max(0,road.projection_iterations)
    [f,~,~,~,~,~,~] = projection_residual_given_s(road,p,s);
    if abs(f) < road.projection_tol
        break;
    end
    ds = max(1.0e-4, 1.0e-6*(1.0+abs(s)));
    fp = projection_residual_given_s(road,p,s+ds);
    fm = projection_residual_given_s(road,p,s-ds);
    dfds = (fp-fm)/(2*ds);
    if ~isfinite(dfds) || abs(dfds) < 1.0e-12
        break;
    end
    step = min(max(f/dfds, -5.0), 5.0);
    s = s - step;
    if abs(step) < road.projection_tol
        break;
    end
end
[~,n,p_road,e_long,e_lat,e_z,~] = projection_residual_given_s(road,p,s);
end

function [f,n,p_road,e_long,e_lat,e_z,p_s] = projection_residual_given_s(road,p,s)
[C,~,~,e_y,~,~] = centerline_quantities(road,s);
n = dot(p-C,e_y);
[p_road,e_long,e_lat,e_z,p_s,~] = road_surface_at(road,s,n);
r = p-p_road;
f = dot(r,p_s);
end

function [p_road,e_long,e_lat,e_z,p_s,p_n] = road_surface_at(road,s,n)
[C,Cs,ex,ey,dey,ezc] = centerline_quantities(road,s);
p_road = C + n*ey;
p_s = Cs + n*dey;
p_n = ey;
e_long = normalize_num(p_s,ex);
e_lat = normalize_num(p_n,ey);
e_z = normalize_num(cross(e_long,e_lat),ezc);
if e_z(3) < 0.0
    e_z = -e_z;
    e_lat = -e_lat;
end
end

function [C,Cs,ex,ey,dey,ezc] = centerline_quantities(road,s)
if isfield(road,'type') && strcmpi(road.type,'flat')
    C = [s;0;0];
    Cs = [1;0;0];
    ex = [1;0;0];
    ey = [0;1;0];
    dey = [0;0;0];
    ezc = [0;0;1];
    return;
end
ke = 2*pi/road.elev_wavelength;
kb = 2*pi/road.bank_wavelength;
zc = road.elev_amp*(1-cos(ke*s));
dz = road.elev_amp*ke*sin(ke*s);
d2z = road.elev_amp*ke*ke*cos(ke*s);
bank = road.bank_amp*sin(kb*s);
db = road.bank_amp*kb*cos(kb*s);
C = [s;0;zc];
Cs = [1;0;dz];
L = sqrt(1+dz*dz);
ex = [1/L;0;dz/L];
ey0 = [0;1;0];
ez0 = [-dz/L;0;1/L];
dez = [-d2z/(L^3);0;-dz*d2z/(L^3)];
cb = cos(bank);
sb = sin(bank);
ey = cb*ey0 + sb*ez0;
dey = -sb*db*ey0 + cb*db*ez0 + sb*dez;
ezc = normalize_num(cross(ex,ey),[0;0;1]);
end

function y = smooth_step(t, t0, t1)
if t <= t0
    y = 0.0;
elseif t >= t1
    y = 1.0;
else
    a = (t-t0)/(t1-t0);
    y = a*a*(3.0 - 2.0*a);
end
end

function y = smoothplus(x,epsv)
y = 0.5*(x + sqrt(x*x + epsv*epsv));
end

function y = smooth_max(a,b,epsv)
d = b - a;
y = a + smoothplus(d,epsv);
end

function y = smooth_abs(x,epsv)
y = sqrt(x*x + epsv*epsv);
end

function y = normalize_num(v,fallback)
n = norm(v);
if n < 1.0e-10
    y = fallback/norm(fallback);
else
    y = v/n;
end
end

function val = getfield_default(s, name, default_val)
if isstruct(s) && isfield(s, name)
    val = s.(name);
else
    val = default_val;
end
end

function d = empty_contact_diag()
d.names = {'FL','FR','RL','RR'};
d.Fx = zeros(4,1);
d.Fy = zeros(4,1);
d.Fz = zeros(4,1);
d.kappa = zeros(4,1);
d.alpha = zeros(4,1);
d.compression = zeros(4,1);
d.road_z_at_wheel = zeros(4,1);
d.contact_x = zeros(4,1);
d.contact_y = zeros(4,1);
d.contact_z = zeros(4,1);
d.mu_ratio = zeros(4,1);
end
