function run_vehicle14_matlab_casadi_demo(varargin)
%RUN_VEHICLE14_MATLAB_CASADI_DEMO Integrate, plot, and render the 14-DOF MATLAB/CasADi demo.
%
% This version lets you choose which outer/component model supplies the Enak input vector IN:
%
%   run_vehicle14_matlab_casadi_demo
%   run_vehicle14_matlab_casadi_demo('outer')
%   run_vehicle14_matlab_casadi_demo('original')
%   run_vehicle14_matlab_casadi_demo('simulate','outer')
%   run_vehicle14_matlab_casadi_demo('simulate','original')
%   run_vehicle14_matlab_casadi_demo('plots-only','outer')
%   run_vehicle14_matlab_casadi_demo('render-only','outer')
%   run_vehicle14_matlab_casadi_demo('simulate','outer','flat')
%   run_vehicle14_matlab_casadi_demo('simulate','outer','ribbon')
%   run_vehicle14_matlab_casadi_demo('original','flat')
%
% Solver policy:
%   1) Try MATLAB ode15i on the implicit residual M(x)*xdot-F(x)=0 first.
%   2) If ode15i fails, fall back to ode15s on the explicit RHS xdot=M\F.
%
% The generated Enak/CasADi file only supplies M, F, and sensors.  The model_source
% decides how tyre/contact/powertrain/brake/aero/control inputs are produced:
%   'outer'    -> vehicle14_outer_models_user.m
%   'original' -> built-in original demo component model in this runner

% addpath('D:/Applications/casadi-windows-matlabR2016a-v3.6.5');  % edit if CasADi is not already on the MATLAB path

[mode, model_source, road_type] = parse_run_mode_source_road(varargin{:});

import casadi.*
[ver_ok, msg] = check_casadi_available();
if ~ver_ok, error('%s', msg); end

[mf_fun, meta] = vehicle14_build_casadi_function();
par = vehicle14_params_default();
args = vehicle14_demo_options();
args = apply_model_source_to_args(args, model_source, road_type);
road = vehicle14_make_demo_road(args);

fprintf('vehicle14 MATLAB/CasADi: mode = %s, model_source = %s, road_type = %s\n', mode, args.model_source, args.road_type);

if strcmpi(mode, 'render-only')
    render_vehicle14_matlab_animation(args.output_mat, args.animation_mp4, args.animation_gif);
    return;
elseif strcmpi(mode, 'plots-only')
    if ~exist(args.output_mat, 'file')
        error('Cannot run plots-only because %s does not exist yet. Run run_vehicle14_matlab_casadi_demo(''simulate'',''%s'') first.', args.output_mat, args.model_source);
    end
    Sload = load(args.output_mat, 't', 'x', 'par', 'args', 'road', 'meta');
    vehicle14_save_plots_matlab(Sload.t, Sload.x, mf_fun, Sload.par, Sload.args, Sload.road, Sload.meta, Sload.args.plot_dir);
    return;
elseif ~strcmpi(mode, 'simulate')
    error('Unknown mode: %s. Use simulate, plots-only, or render-only.', mode);
end

x0 = vehicle14_initial_state(mf_fun, par, args, road, meta);
xdot0 = vehicle14_rhs_full_numeric(0, x0, mf_fun, par, args, road, meta);
tspan = [0 args.t_final];
opts = odeset('RelTol', args.rtol, 'AbsTol', args.atol, 'MaxStep', args.max_step);
solver_used = 'ode15i';
try
    fprintf('Trying implicit MATLAB solver ode15i on M(x)*xdot-F(x)=0...\n');
    [t,x] = ode15i(@(tt,xx,xxp) vehicle14_residual_full_numeric(tt,xx,xxp,mf_fun,par,args,road,meta), tspan, x0, xdot0, opts);
catch ME
    warning('ode15i failed: %s\nFalling back to ode15s with explicit RHS M\\F.', ME.message);
    solver_used = 'ode15s-explicit';
    [t,x] = ode15s(@(tt,xx) vehicle14_rhs_full_numeric(tt,xx,mf_fun,par,args,road,meta), tspan, x0, opts);
end

model_source_saved = args.model_source; %#ok<NASGU>
save(args.output_mat, 't', 'x', 'par', 'args', 'road', 'meta', 'solver_used', 'model_source_saved', '-v7.3');
fprintf('Saved solution to %s using %s, model_source=%s, road_type=%s.\n', args.output_mat, solver_used, args.model_source, args.road_type);
vehicle14_save_plots_matlab(t, x, mf_fun, par, args, road, meta, args.plot_dir);
try
    render_vehicle14_matlab_animation(args.output_mat, args.animation_mp4, args.animation_gif);
catch ME
    warning('Animation rendering failed: %s', ME.message);
end
end

function [mode, model_source, road_type] = parse_run_mode_source_road(varargin)
%PARSE_RUN_MODE_SOURCE_ROAD Flexible parser for mode/source/road options.
%
% Valid modes:
%   simulate, plots-only, render-only
% Valid model sources:
%   outer, original
% Valid road types:
%   ribbon, 3d, 3d-ribbon, road-3d, flat
%
% Examples:
%   run_vehicle14_matlab_casadi_demo
%   run_vehicle14_matlab_casadi_demo('outer')
%   run_vehicle14_matlab_casadi_demo('original','flat')
%   run_vehicle14_matlab_casadi_demo('simulate','outer','flat')
%   run_vehicle14_matlab_casadi_demo('plots-only','outer','ribbon')

valid_modes = {'simulate','plots-only','render-only'};
valid_sources = {'outer','original'};
valid_ribbon_aliases = {'ribbon','3d','3d-ribbon','road-3d','demo-ribbon'};
valid_flat_aliases = {'flat','flat-road','plane'};

mode = 'simulate';
model_source = 'outer';
road_type = 'ribbon';

for ii = 1:nargin
    if isempty(varargin{ii})
        continue;
    end
    token = lower(char(varargin{ii}));
    if any(strcmp(token, valid_modes))
        mode = token;
    elseif any(strcmp(token, valid_sources))
        model_source = token;
    elseif any(strcmp(token, valid_ribbon_aliases))
        road_type = 'ribbon';
    elseif any(strcmp(token, valid_flat_aliases))
        road_type = 'flat';
    else
        error('Unknown argument: %s. Use mode={simulate,plots-only,render-only}, source={outer,original}, road={ribbon,flat}.', token);
    end
end
end

function args = apply_model_source_to_args(args, model_source, road_type)
args.model_source = lower(char(model_source));
args.road_type = lower(char(road_type));
args.output_mat = sprintf('vehicle14_matlab_casadi_solution_%s_%s.mat', args.model_source, args.road_type);
args.plot_dir = sprintf('vehicle14_matlab_casadi_plots_%s_%s', args.model_source, args.road_type);
args.animation_mp4 = sprintf('vehicle14_matlab_casadi_animation_%s_%s.mp4', args.model_source, args.road_type);
args.animation_gif = sprintf('vehicle14_matlab_casadi_animation_%s_%s.gif', args.model_source, args.road_type);
end

function [ok,msg] = check_casadi_available()
ok = true; msg = '';
try
    import casadi.*
    SX.sym('probe',1,1);
catch ME
    ok = false;
    msg = sprintf('CasADi is not on the MATLAB path or failed to initialize: %s', ME.message);
end
end


% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_build_casadi_function.m
% -------------------------------------------------------------------------
function [mf_fun, meta] = vehicle14_build_casadi_function()
%VEHICLE14_BUILD_CASADI_FUNCTION Build the numeric CasADi Function used by MATLAB solvers.
import casadi.*
meta = vehicle14_Enak_meta();
q  = SX.sym('q',  meta.nq, 1);
u  = SX.sym('u',  meta.nu, 1);
P  = SX.sym('P',  meta.np, 1);
IN = SX.sym('IN', meta.nin, 1);
[M,F,S] = vehicle14_enak_matrices_casadi(q,u,P,IN);
mf_fun = Function('vehicle14_mf', {q,u,P,IN}, {M,F,S}, {'q','u','P','IN'}, {'M','F','S'});
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_Enak_meta.m
% -------------------------------------------------------------------------
function meta = vehicle14_Enak_meta()
% Metadata matching vehicle14_Enak_matrices_casadi.m.
meta.q_names = {'X', 'Y', 'Z', 'psi', 'theta', 'phi', 'zFL', 'zFR', 'zRL', 'zRR', 'chiFL', 'chiFR', 'chiRL', 'chiRR'};
meta.u_names = {'vx', 'vy', 'vz', 'p', 'q_pitch', 'r', 'dzFL', 'dzFR', 'dzRL', 'dzRR', 'omFL', 'omFR', 'omRL', 'omRR'};
meta.p_names = {'mB', 'mW', 'Ixx', 'Iyy', 'Izz', 'Ixz', 'Iw', 'lf', 'lr', 'tf', 'tr', 'hc', 'rw', 'g', 'ksf', 'ksr', 'csf', 'csr', 'karbf', 'karbr', 'carbf', 'carbr', 'Fpre_FL', 'Fpre_FR', 'Fpre_RL', 'Fpre_RR', 'rack_m_per_rad', 'kcx1_f', 'kcx2_f', 'kcx1_r', 'kcx2_r', 'kcy1_f', 'kcy2_f', 'kcy1_r', 'kcy2_r', 'kcz1_f', 'kcz1_r', 'kcz2_f', 'kcz2_r', 'kcx_y_f', 'kcx_zy_f', 'kcx_y2_f', 'kcy_y_f', 'kcy_zy_f', 'kcy_y2_f', 'kcz_y_f', 'kcz_zy_f', 'kcz_y2_f', 'toe0_f', 'toe1_f', 'toe2_f', 'toe0_r', 'toe1_r', 'toe2_r', 'toe_ysr_f', 'toe_zysr_f', 'toe_ysr2_f', 'pitch0_f', 'pitch1_f', 'pitch2_f', 'pitch0_r', 'pitch1_r', 'pitch2_r', 'pitch_ysr_f', 'pitch_zysr_f', 'pitch_ysr2_f', 'cam0_f', 'cam1_f', 'cam2_f', 'cam0_r', 'cam1_r', 'cam2_r', 'cam_ysr_f', 'cam_zysr_f', 'cam_ysr2_f', 'aero_front_x', 'aero_rear_x', 'aero_z'};
meta.input_names = {'delta', 'delta_dot', 'FBxFL', 'FByFL', 'FBzFL', 'MBxFL', 'MByFL', 'MBzFL', 'FBxFR', 'FByFR', 'FBzFR', 'MBxFR', 'MByFR', 'MBzFR', 'FBxRL', 'FByRL', 'FBzRL', 'MBxRL', 'MByRL', 'MBzRL', 'FBxRR', 'FByRR', 'FBzRR', 'MBxRR', 'MByRR', 'MBzRR', 'TFL', 'TFR', 'TRL', 'TRR', 'FdragF', 'FdownF', 'FdragR', 'FdownR'};
meta.sensor_names = {'FL_Wc_x', 'FL_Wc_y', 'FL_Wc_z', 'FL_Wc_vx', 'FL_Wc_vy', 'FL_Wc_vz', 'FL_Vx_wc', 'FL_Vy_wc', 'FL_x_body', 'FL_y_body', 'FL_z_body', 'FL_toe', 'FL_camber', 'FL_pitch', 'FL_rack', 'FL_head_Nx', 'FL_head_Ny', 'FL_head_Nz', 'FL_spin_Nx', 'FL_spin_Ny', 'FL_spin_Nz', 'FL_carrier_wx', 'FL_carrier_wy', 'FL_carrier_wz', 'FR_Wc_x', 'FR_Wc_y', 'FR_Wc_z', 'FR_Wc_vx', 'FR_Wc_vy', 'FR_Wc_vz', 'FR_Vx_wc', 'FR_Vy_wc', 'FR_x_body', 'FR_y_body', 'FR_z_body', 'FR_toe', 'FR_camber', 'FR_pitch', 'FR_rack', 'FR_head_Nx', 'FR_head_Ny', 'FR_head_Nz', 'FR_spin_Nx', 'FR_spin_Ny', 'FR_spin_Nz', 'FR_carrier_wx', 'FR_carrier_wy', 'FR_carrier_wz', 'RL_Wc_x', 'RL_Wc_y', 'RL_Wc_z', 'RL_Wc_vx', 'RL_Wc_vy', 'RL_Wc_vz', 'RL_Vx_wc', 'RL_Vy_wc', 'RL_x_body', 'RL_y_body', 'RL_z_body', 'RL_toe', 'RL_camber', 'RL_pitch', 'RL_rack', 'RL_head_Nx', 'RL_head_Ny', 'RL_head_Nz', 'RL_spin_Nx', 'RL_spin_Ny', 'RL_spin_Nz', 'RL_carrier_wx', 'RL_carrier_wy', 'RL_carrier_wz', 'RR_Wc_x', 'RR_Wc_y', 'RR_Wc_z', 'RR_Wc_vx', 'RR_Wc_vy', 'RR_Wc_vz', 'RR_Vx_wc', 'RR_Vy_wc', 'RR_x_body', 'RR_y_body', 'RR_z_body', 'RR_toe', 'RR_camber', 'RR_pitch', 'RR_rack', 'RR_head_Nx', 'RR_head_Ny', 'RR_head_Nz', 'RR_spin_Nx', 'RR_spin_Ny', 'RR_spin_Nz', 'RR_carrier_wx', 'RR_carrier_wy', 'RR_carrier_wz', 'aero_front_z', 'aero_rear_z', 'vx_body'};
meta.nq = 14; meta.nu = 14; meta.nx = 28;
meta.np = 78; meta.nin = 34; meta.ns = 99;
meta.input_index = struct();
meta.input_index.delta = 1;
meta.input_index.delta_dot = 2;
meta.input_index.FBxFL = 3;
meta.input_index.FByFL = 4;
meta.input_index.FBzFL = 5;
meta.input_index.MBxFL = 6;
meta.input_index.MByFL = 7;
meta.input_index.MBzFL = 8;
meta.input_index.FBxFR = 9;
meta.input_index.FByFR = 10;
meta.input_index.FBzFR = 11;
meta.input_index.MBxFR = 12;
meta.input_index.MByFR = 13;
meta.input_index.MBzFR = 14;
meta.input_index.FBxRL = 15;
meta.input_index.FByRL = 16;
meta.input_index.FBzRL = 17;
meta.input_index.MBxRL = 18;
meta.input_index.MByRL = 19;
meta.input_index.MBzRL = 20;
meta.input_index.FBxRR = 21;
meta.input_index.FByRR = 22;
meta.input_index.FBzRR = 23;
meta.input_index.MBxRR = 24;
meta.input_index.MByRR = 25;
meta.input_index.MBzRR = 26;
meta.input_index.TFL = 27;
meta.input_index.TFR = 28;
meta.input_index.TRL = 29;
meta.input_index.TRR = 30;
meta.input_index.FdragF = 31;
meta.input_index.FdownF = 32;
meta.input_index.FdragR = 33;
meta.input_index.FdownR = 34;
meta.sensor_index = struct();
meta.sensor_index.FL_Wc_x = 1;
meta.sensor_index.FL_Wc_y = 2;
meta.sensor_index.FL_Wc_z = 3;
meta.sensor_index.FL_Wc_vx = 4;
meta.sensor_index.FL_Wc_vy = 5;
meta.sensor_index.FL_Wc_vz = 6;
meta.sensor_index.FL_Vx_wc = 7;
meta.sensor_index.FL_Vy_wc = 8;
meta.sensor_index.FL_x_body = 9;
meta.sensor_index.FL_y_body = 10;
meta.sensor_index.FL_z_body = 11;
meta.sensor_index.FL_toe = 12;
meta.sensor_index.FL_camber = 13;
meta.sensor_index.FL_pitch = 14;
meta.sensor_index.FL_rack = 15;
meta.sensor_index.FL_head_Nx = 16;
meta.sensor_index.FL_head_Ny = 17;
meta.sensor_index.FL_head_Nz = 18;
meta.sensor_index.FL_spin_Nx = 19;
meta.sensor_index.FL_spin_Ny = 20;
meta.sensor_index.FL_spin_Nz = 21;
meta.sensor_index.FL_carrier_wx = 22;
meta.sensor_index.FL_carrier_wy = 23;
meta.sensor_index.FL_carrier_wz = 24;
meta.sensor_index.FR_Wc_x = 25;
meta.sensor_index.FR_Wc_y = 26;
meta.sensor_index.FR_Wc_z = 27;
meta.sensor_index.FR_Wc_vx = 28;
meta.sensor_index.FR_Wc_vy = 29;
meta.sensor_index.FR_Wc_vz = 30;
meta.sensor_index.FR_Vx_wc = 31;
meta.sensor_index.FR_Vy_wc = 32;
meta.sensor_index.FR_x_body = 33;
meta.sensor_index.FR_y_body = 34;
meta.sensor_index.FR_z_body = 35;
meta.sensor_index.FR_toe = 36;
meta.sensor_index.FR_camber = 37;
meta.sensor_index.FR_pitch = 38;
meta.sensor_index.FR_rack = 39;
meta.sensor_index.FR_head_Nx = 40;
meta.sensor_index.FR_head_Ny = 41;
meta.sensor_index.FR_head_Nz = 42;
meta.sensor_index.FR_spin_Nx = 43;
meta.sensor_index.FR_spin_Ny = 44;
meta.sensor_index.FR_spin_Nz = 45;
meta.sensor_index.FR_carrier_wx = 46;
meta.sensor_index.FR_carrier_wy = 47;
meta.sensor_index.FR_carrier_wz = 48;
meta.sensor_index.RL_Wc_x = 49;
meta.sensor_index.RL_Wc_y = 50;
meta.sensor_index.RL_Wc_z = 51;
meta.sensor_index.RL_Wc_vx = 52;
meta.sensor_index.RL_Wc_vy = 53;
meta.sensor_index.RL_Wc_vz = 54;
meta.sensor_index.RL_Vx_wc = 55;
meta.sensor_index.RL_Vy_wc = 56;
meta.sensor_index.RL_x_body = 57;
meta.sensor_index.RL_y_body = 58;
meta.sensor_index.RL_z_body = 59;
meta.sensor_index.RL_toe = 60;
meta.sensor_index.RL_camber = 61;
meta.sensor_index.RL_pitch = 62;
meta.sensor_index.RL_rack = 63;
meta.sensor_index.RL_head_Nx = 64;
meta.sensor_index.RL_head_Ny = 65;
meta.sensor_index.RL_head_Nz = 66;
meta.sensor_index.RL_spin_Nx = 67;
meta.sensor_index.RL_spin_Ny = 68;
meta.sensor_index.RL_spin_Nz = 69;
meta.sensor_index.RL_carrier_wx = 70;
meta.sensor_index.RL_carrier_wy = 71;
meta.sensor_index.RL_carrier_wz = 72;
meta.sensor_index.RR_Wc_x = 73;
meta.sensor_index.RR_Wc_y = 74;
meta.sensor_index.RR_Wc_z = 75;
meta.sensor_index.RR_Wc_vx = 76;
meta.sensor_index.RR_Wc_vy = 77;
meta.sensor_index.RR_Wc_vz = 78;
meta.sensor_index.RR_Vx_wc = 79;
meta.sensor_index.RR_Vy_wc = 80;
meta.sensor_index.RR_x_body = 81;
meta.sensor_index.RR_y_body = 82;
meta.sensor_index.RR_z_body = 83;
meta.sensor_index.RR_toe = 84;
meta.sensor_index.RR_camber = 85;
meta.sensor_index.RR_pitch = 86;
meta.sensor_index.RR_rack = 87;
meta.sensor_index.RR_head_Nx = 88;
meta.sensor_index.RR_head_Ny = 89;
meta.sensor_index.RR_head_Nz = 90;
meta.sensor_index.RR_spin_Nx = 91;
meta.sensor_index.RR_spin_Ny = 92;
meta.sensor_index.RR_spin_Nz = 93;
meta.sensor_index.RR_carrier_wx = 94;
meta.sensor_index.RR_carrier_wy = 95;
meta.sensor_index.RR_carrier_wz = 96;
meta.sensor_index.aero_front_z = 97;
meta.sensor_index.aero_rear_z = 98;
meta.sensor_index.vx_body = 99;
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_params_default.m
% -------------------------------------------------------------------------
function par = vehicle14_params_default()
% Default parameters copied from the uploaded Python VehicleParams dataclass.
par.mB = 1250.0;
par.mW = 42.0;
par.Ixx = 650.0;
par.Iyy = 1500.0;
par.Izz = 2200.0;
par.Ixz = 0.0;
par.Iw = 1.25;
par.lf = 1.35;
par.lr = 1.35;
par.tf = 1.6;
par.tr = 1.6;
par.hc = 0.55;
par.rw = 0.31;
par.g = 9.81;
par.ksf = 36000.0;
par.ksr = 39000.0;
par.csf = 3600.0;
par.csr = 3900.0;
par.karbf = 12000.0;
par.karbr = 10000.0;
par.carbf = 700.0;
par.carbr = 600.0;
par.Fpre_FL = NaN;
par.Fpre_FR = NaN;
par.Fpre_RL = NaN;
par.Fpre_RR = NaN;
par.rack_m_per_rad = 0.05;
par.kcx1_f = 0.0;
par.kcx2_f = 0.0;
par.kcx1_r = 0.0;
par.kcx2_r = 0.0;
par.kcy1_f = 0.0;
par.kcy2_f = 0.0;
par.kcy1_r = 0.0;
par.kcy2_r = 0.0;
par.kcz1_f = 1.0;
par.kcz1_r = 1.0;
par.kcz2_f = 0.0;
par.kcz2_r = 0.0;
par.kcx_y_f = 0.0;
par.kcx_zy_f = 0.0;
par.kcx_y2_f = 0.0;
par.kcy_y_f = 0.0;
par.kcy_zy_f = 0.0;
par.kcy_y2_f = 0.0;
par.kcz_y_f = 0.0;
par.kcz_zy_f = 0.0;
par.kcz_y2_f = 0.0;
par.toe0_f = 0.0;
par.toe1_f = 0.0;
par.toe2_f = 0.0;
par.toe0_r = 0.0;
par.toe1_r = 0.0;
par.toe2_r = 0.0;
par.toe_ysr_f = 20.0;
par.toe_zysr_f = 0.0;
par.toe_ysr2_f = 0.0;
par.pitch0_f = 0.0;
par.pitch1_f = 0.0;
par.pitch2_f = 0.0;
par.pitch0_r = 0.0;
par.pitch1_r = 0.0;
par.pitch2_r = 0.0;
par.pitch_ysr_f = 0.0;
par.pitch_zysr_f = 0.0;
par.pitch_ysr2_f = 0.0;
par.cam0_f = 0.0;
par.cam1_f = -0.06;
par.cam2_f = 0.0;
par.cam0_r = 0.0;
par.cam1_r = -0.04;
par.cam2_r = 0.0;
par.cam_ysr_f = 0.0;
par.cam_zysr_f = 0.0;
par.cam_ysr2_f = 0.0;
par.kt = 220000.0;
par.ct = 1400.0;
par.normal_smooth_eps = 30.0;
par.mu = 1.15;
par.Ck = 90000.0;
par.Ca = 85000.0;
par.v_eps = 1.0;
par.T_drive_max = 1900.0;
par.T_brake_max = 4500.0;
par.brake_bias_front = 0.65;
par.steer_max = 0.12217304763960307;
par.diff_gain = 35.0;
par.diff_width = 15.0;
par.rho = 1.225;
par.CdA = 0.65;
par.ClA_front = 0.48;
par.ClA_rear = 0.72;
par.aero_front_x = 1.25;
par.aero_rear_x = -1.2;
par.aero_z = 0.05;
par.aero_h_ref = 0.08;
par.aero_h_sens_front = 1.0;
par.aero_h_sens_rear = 1.2;
par.corner_names = {'FL','FR','RL','RR'};
par.p_names = {'mB', 'mW', 'Ixx', 'Iyy', 'Izz', 'Ixz', 'Iw', 'lf', 'lr', 'tf', 'tr', 'hc', 'rw', 'g', 'ksf', 'ksr', 'csf', 'csr', 'karbf', 'karbr', 'carbf', 'carbr', 'Fpre_FL', 'Fpre_FR', 'Fpre_RL', 'Fpre_RR', 'rack_m_per_rad', 'kcx1_f', 'kcx2_f', 'kcx1_r', 'kcx2_r', 'kcy1_f', 'kcy2_f', 'kcy1_r', 'kcy2_r', 'kcz1_f', 'kcz1_r', 'kcz2_f', 'kcz2_r', 'kcx_y_f', 'kcx_zy_f', 'kcx_y2_f', 'kcy_y_f', 'kcy_zy_f', 'kcy_y2_f', 'kcz_y_f', 'kcz_zy_f', 'kcz_y2_f', 'toe0_f', 'toe1_f', 'toe2_f', 'toe0_r', 'toe1_r', 'toe2_r', 'toe_ysr_f', 'toe_zysr_f', 'toe_ysr2_f', 'pitch0_f', 'pitch1_f', 'pitch2_f', 'pitch0_r', 'pitch1_r', 'pitch2_r', 'pitch_ysr_f', 'pitch_zysr_f', 'pitch_ysr2_f', 'cam0_f', 'cam1_f', 'cam2_f', 'cam0_r', 'cam1_r', 'cam2_r', 'cam_ysr_f', 'cam_zysr_f', 'cam_ysr2_f', 'aero_front_x', 'aero_rear_x', 'aero_z'};
par.P = [1250.0; 42.0; 650.0; 1500.0; 2200.0; 0.0; 1.25; 1.35; 1.35; 1.6; 1.6; 0.55; 0.31; 9.81; 36000.0; 39000.0; 3600.0; 3900.0; 12000.0; 10000.0; 700.0; 600.0; 3065.625; 3065.625; 3065.625; 3065.625; 0.05; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 20.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; -0.06; 0.0; 0.0; -0.04; 0.0; 0.0; 0.0; 0.0; 1.25; -1.2; 0.05];
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_demo_options.m
% -------------------------------------------------------------------------
function args = vehicle14_demo_options()
%VEHICLE14_DEMO_OPTIONS Options matching the uploaded Python drop-then-maneuver demo.
args.drop_height = 0.25;
args.settle_time = 2.0;
args.maneuver_time = 16.0;
args.t_final = args.settle_time + args.maneuver_time;
args.max_step = 0.006;
args.rtol = 1e-7;
args.atol = 1e-9;
args.suspension_init = 'zero-force';
args.initial_x = 0.0;
args.initial_y = 0.0;
args.initial_vx = 0.0;
args.initial_vz = 0.0;
args.initial_wheel_omega = 0.0;
args.maneuver_profile = 'aggressive-3d';
args.throttle = 0.65;
args.throttle_start = 0.2;
args.throttle_ramp = 0.45;
args.steer_deg = 7.0;
args.steer_start = 1.0;
args.steer_ramp = 0.5;
args.steer_freq = 0.55;
args.brake_start = 11.5;
args.brake_level = 0.28;
args.brake_ramp = 0.35;
args.brake_hold = 1.15;
args.disable_aero = false;
args.road_type = 'ribbon';   % 'ribbon' for the 3D demo road, or 'flat' for a pure flat road
args.road_width = 12.0;
args.road_elev_amp = 1.0;
args.road_elev_wavelength = 80.0;
args.road_bank_deg = 2.0;
args.road_bank_wavelength = 60.0;
args.road_projection_iterations = 0;
args.output_mat = 'vehicle14_matlab_casadi_solution.mat';
args.plot_dir = 'vehicle14_matlab_casadi_plots';
args.animation_mp4 = 'vehicle14_matlab_casadi_animation.mp4';
args.animation_gif = 'vehicle14_matlab_casadi_animation.gif';
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_make_demo_road.m
% -------------------------------------------------------------------------
function road = vehicle14_make_demo_road(args)
%VEHICLE14_MAKE_DEMO_ROAD Build either the 3D ribbon road or a pure flat road.
road.type = lower(getfield_default_runner(args, 'road_type', 'ribbon'));
road.width = args.road_width;
road.projection_iterations = args.road_projection_iterations;
road.projection_tol = 1e-11;

if strcmp(road.type, 'flat')
    % Pure horizontal plane: z=0, no banking, no elevation variation.
    road.elev_amp = 0.0;
    road.elev_wavelength = 1.0;
    road.bank_amp = 0.0;
    road.bank_wavelength = 1.0;
    road.projection_iterations = 0;
elseif strcmp(road.type, 'ribbon')
    % Original 3D straight ribbon demo road.
    road.elev_amp = args.road_elev_amp;
    road.elev_wavelength = args.road_elev_wavelength;
    road.bank_amp = deg2rad(args.road_bank_deg);
    road.bank_wavelength = args.road_bank_wavelength;
else
    error('Unknown road.type: %s. Use ribbon or flat.', road.type);
end
end

function [throttle, brake, delta, delta_dot] = vehicle14_command_inputs(t_abs, par, args)
%VEHICLE14_COMMAND_INPUTS Driver command profiles ported from the uploaded Python demo.
[throttle, brake, delta] = command_no_derivative(t_abs, par, args);
h = 1.0e-5;
[~,~,delta_p] = command_no_derivative(t_abs + h, par, args);
[~,~,delta_m] = command_no_derivative(t_abs - h, par, args);
delta_dot = (delta_p - delta_m)/(2*h);
end

function [throttle, brake, delta] = command_no_derivative(t_abs, par, args)
if t_abs < args.settle_time
    throttle = 0; brake = 0; delta = 0; return;
end
tau = t_abs - args.settle_time;
profile = args.maneuver_profile;
if strcmp(profile, 'straight-launch')
    throttle = args.throttle * smooth_step(tau, args.throttle_start, args.throttle_start + args.throttle_ramp);
    brake = 0; delta = 0; return;
elseif strcmp(profile, 'creep-steer')
    throttle = args.throttle * smooth_step(tau, args.throttle_start, args.throttle_start + args.throttle_ramp);
    brake = 0;
    if args.brake_start >= 0 && tau >= args.brake_start
        bw = smooth_step(tau, args.brake_start, args.brake_start + args.brake_ramp);
        brake = args.brake_level*bw;
        throttle = throttle*(1-bw);
    end
    steer_amp = deg2rad(args.steer_deg);
    env = smooth_step(tau, args.steer_start, args.steer_start + args.steer_ramp);
    delta = steer_amp*env*sin(2*pi*args.steer_freq*max(tau-args.steer_start,0));
    return;
elseif strcmp(profile, 'source-demo')
    throttle = 0.17 * (1.0 - smooth_step(tau, 2.2, 2.7));
    brake = 0.0;
    if 3.2 <= tau && tau <= 4.8
        brake = 0.20 * smooth_step(tau, 3.2, 3.5) * (1.0 - smooth_step(tau, 4.4, 4.8));
    end
    envelope = smooth_step(tau, 0.8, 1.2) * (1.0 - smooth_step(tau, 5.4, 6.0));
    delta = par.steer_max * envelope * sin(2*pi*0.32*max(tau-1.0,0));
    return;
elseif strcmp(profile, 'aggressive-3d')
    launch = smooth_step(tau, 0.10, 0.55);
    throttle = args.throttle * launch;
    brake_window = 0;
    if args.brake_start >= 0
        brake_window = smooth_step(tau, args.brake_start, args.brake_start + args.brake_ramp) * ...
            (1 - smooth_step(tau, args.brake_start + args.brake_hold, args.brake_start + args.brake_hold + args.brake_ramp));
    end
    brake = args.brake_level * brake_window;
    throttle = throttle * (1 - 0.85*brake_window);
    steer_amp = deg2rad(args.steer_deg);
    env = smooth_step(tau, args.steer_start, args.steer_start + args.steer_ramp) * ...
        (1 - smooth_step(tau, args.maneuver_time - 1.0, args.maneuver_time - 0.2));
    tt = max(tau - args.steer_start, 0);
    raw = 0.78*sin(2*pi*args.steer_freq*tt) + 0.28*sin(2*pi*1.75*args.steer_freq*tt + 0.45);
    delta = steer_amp * env * min(max(raw, -1), 1);
    return;
else
    error('Unknown maneuver_profile: %s', profile);
end
end

function y = smooth_step(t, t0, t1)
if t <= t0
    y = 0;
elseif t >= t1
    y = 1;
else
    a = (t-t0)/(t1-t0);
    y = a*a*(3 - 2*a);
end
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_initial_state.m
% -------------------------------------------------------------------------
function x0 = vehicle14_initial_state(mf_fun, par, args, road, meta)
%VEHICLE14_INITIAL_STATE Sensor-based drop initial condition with requested tyre gap.
q0 = zeros(14,1);
u0 = zeros(14,1);
q0(1) = args.initial_x;
q0(2) = args.initial_y;
q0(7:10) = initial_suspension_coordinates(par, args.suspension_init);
u0(1) = args.initial_vx;
u0(3) = args.initial_vz;
u0(11:14) = args.initial_wheel_omega;
q0(3) = par.hc + args.drop_height;
for k = 1:8
    gap = min_contact_gap(q0, u0, mf_fun, par, road, meta);
    err = gap - args.drop_height;
    if abs(err) < 1e-10, break; end
    dZ = 1e-5;
    qp = q0; qm = q0; qp(3) = qp(3) + dZ; qm(3) = qm(3) - dZ;
    gp = min_contact_gap(qp, u0, mf_fun, par, road, meta);
    gm = min_contact_gap(qm, u0, mf_fun, par, road, meta);
    dgap = (gp-gm)/(2*dZ);
    if ~isfinite(dgap) || abs(dgap) < 1e-9
        q0(3) = q0(3) - err;
    else
        q0(3) = q0(3) - min(max(err/dgap, -0.5), 0.5);
    end
end
x0 = [q0; u0];
end

function z = initial_suspension_coordinates(par, mode)
Fpre = resolved_preloads(par);
if strcmp(mode, 'zero-force')
    z = [-Fpre(1)/par.ksf; -Fpre(2)/par.ksf; -Fpre(3)/par.ksr; -Fpre(4)/par.ksr];
elseif strcmp(mode, 'nominal') || strcmp(mode, 'static-source')
    z = zeros(4,1);
else
    error('Unknown suspension_init: %s', mode);
end
end

function Fpre = resolved_preloads(par)
total_w = par.mB*par.g;
front_axle = total_w*par.lr/(par.lf+par.lr);
rear_axle = total_w*par.lf/(par.lf+par.lr);
defaults = [0.5*front_axle; 0.5*front_axle; 0.5*rear_axle; 0.5*rear_axle];
user = [par.Fpre_FL; par.Fpre_FR; par.Fpre_RL; par.Fpre_RR];
Fpre = defaults;
for i=1:4
    if ~isnan(user(i)), Fpre(i)=user(i); end
end
end

function gap = min_contact_gap(q0, u0, mf_fun, par, road, meta)
IN0 = zeros(meta.nin,1);
[~,~,Sdm] = mf_fun(q0, u0, par.P, IN0);
S = full(Sdm);
idx = meta.sensor_index;
names = {'FL','FR','RL','RR'};
gaps = zeros(4,1);
for i=1:4
    nm = names{i};
    p_wc = [S(idx.([nm '_Wc_x'])); S(idx.([nm '_Wc_y'])); S(idx.([nm '_Wc_z']))];
    head = [S(idx.([nm '_head_Nx'])); S(idx.([nm '_head_Ny'])); S(idx.([nm '_head_Nz']))];
    spin = [S(idx.([nm '_spin_Nx'])); S(idx.([nm '_spin_Ny'])); S(idx.([nm '_spin_Nz']))];
    geom = oriented_tyre_contact_geometry_local(par, road, p_wc, head, spin);
    gaps(i) = -geom.compression;
end
gap = min(gaps);
end

% Local minimal geometry copies needed before the main contact file is loaded.
function geom = oriented_tyre_contact_geometry_local(par, road, p_wc_N, wheel_heading_N, wheel_spin_axis_N)
[s,n,p_road_N,e_long_N,e_lat_N,n_road_N] = road_project_point_local(road, p_wc_N);
wheel_heading_N = normalize_num_local(wheel_heading_N, e_long_N);
spin_axis_N = normalize_num_local(wheel_spin_axis_N, e_lat_N);
t_long_raw = cross(spin_axis_N, n_road_N);
if norm(t_long_raw) < 1e-10
    t_long_raw = wheel_heading_N - dot(wheel_heading_N,n_road_N)*n_road_N;
end
t_long_N = normalize_num_local(t_long_raw, e_long_N);
if dot(t_long_N, wheel_heading_N) < 0, t_long_N = -t_long_N; end
t_lat_N = normalize_num_local(cross(n_road_N,t_long_N), e_lat_N);
normal_in_radial_plane = n_road_N - dot(n_road_N, spin_axis_N)*spin_axis_N;
gamma = norm(normal_in_radial_plane);
if gamma < 1e-10
    gamma = 1.0; radial_down_N = -n_road_N;
else
    radial_down_N = -normal_in_radial_plane/gamma;
end
signed_center_distance = dot(p_wc_N - p_road_N, n_road_N);
radial_depth_normal = par.rw*gamma;
compression_normal = radial_depth_normal - signed_center_distance;
distance_center_to_road_along_radial = signed_center_distance/gamma;
compression_radial = par.rw - distance_center_to_road_along_radial;
r_wc_to_contact_N = distance_center_to_road_along_radial*radial_down_N;
p_contact_N = p_wc_N + r_wc_to_contact_N;
geom.s=s; geom.n=n; geom.p_road_N=p_road_N; geom.e_long_N=e_long_N; geom.e_lat_N=e_lat_N; geom.n_road_N=n_road_N;
geom.t_long_N=t_long_N; geom.t_lat_N=t_lat_N; geom.spin_axis_N=spin_axis_N; geom.radial_down_N=radial_down_N;
geom.r_wc_to_contact_N=r_wc_to_contact_N; geom.p_contact_N=p_contact_N;
geom.compression=compression_normal; geom.compression_normal=compression_normal; geom.compression_radial=compression_radial;
geom.contact_gap_normal=-compression_normal; geom.contact_gap_radial=-compression_radial;
end

function [s,n,p_road,e_long,e_lat,e_z] = road_project_point_local(road,p)
s = p(1);
[C,~,~,e_y,~,~] = centerline_quantities_local(road,s);
n = dot(p-C,e_y);
[p_road,e_long,e_lat,e_z] = road_surface_at_local(road,s,n);
end

function [p_road,e_long,e_lat,e_z] = road_surface_at_local(road,s,n)
[C,Cs,ex,ey,dey,ezc] = centerline_quantities_local(road,s);
p_road = C + n*ey; ps = Cs + n*dey; pn = ey;
e_long = normalize_num_local(ps,ex); e_lat = normalize_num_local(pn,ey); e_z = normalize_num_local(cross(e_long,e_lat),ezc);
if e_z(3)<0, e_z=-e_z; e_lat=-e_lat; end
end

function [C,Cs,ex,ey,dey,ezc] = centerline_quantities_local(road,s)
if isfield(road,'type') && strcmpi(road.type,'flat')
    C = [s;0;0];
    Cs = [1;0;0];
    ex = [1;0;0];
    ey = [0;1;0];
    dey = [0;0;0];
    ezc = [0;0;1];
    return;
end
ke=2*pi/road.elev_wavelength; kb=2*pi/road.bank_wavelength;
zc=road.elev_amp*(1-cos(ke*s)); dz=road.elev_amp*ke*sin(ke*s); d2z=road.elev_amp*ke*ke*cos(ke*s);
bank=road.bank_amp*sin(kb*s); db=road.bank_amp*kb*cos(kb*s);
C=[s;0;zc]; Cs=[1;0;dz]; L=sqrt(1+dz*dz); ex=[1/L;0;dz/L]; ey0=[0;1;0]; ez0=[-dz/L;0;1/L];
dez=[-d2z/(L^3);0;-dz*d2z/(L^3)]; cb=cos(bank); sb=sin(bank); ey=cb*ey0+sb*ez0; dey=-sb*db*ey0+cb*db*ez0+sb*dez; ezc=normalize_num_local(cross(ex,ey),[0;0;1]);
end

function y = normalize_num_local(v,fallback)
n=norm(v); if n<1e-10, y=fallback/norm(fallback); else, y=v/n; end
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_core_numeric.m
% -------------------------------------------------------------------------
function [M,F,aux] = vehicle14_core_numeric(t, x, mf_fun, par, args, road, meta)
%VEHICLE14_CORE_NUMERIC Evaluate M(x)*xdot = F(x) using selected input/component model.
q = x(1:14);
u = x(15:28);
source = lower(getfield_default_runner(args, 'model_source', 'outer'));

if strcmp(source, 'outer')
    % New user-editable outer-model path.  First ask only for steering so
    % the generated sensor geometry uses steering-consistent front K/C maps.
    [IN0, ctrl] = vehicle14_outer_models_user('control', t, x, [], par, args, road, meta);
    [~,~,Sdm] = mf_fun(q, u, par.P, IN0);
    S = full(Sdm);

    % Now compute the complete IN vector from the user-editable outer model.
    [IN, outer] = vehicle14_outer_models_user('full', t, x, S, par, args, road, meta);
    [Mdm,Fdm,~] = mf_fun(q, u, par.P, IN);
    M = full(Mdm);
    F = full(Fdm);

    if nargout > 2
        aux.S = S;
        aux.IN = IN;
        aux.control = ctrl;
        aux.outer = outer;
        aux.contact = outer.contact;
        aux.throttle = outer.throttle;
        aux.brake = outer.brake;
        aux.delta = outer.delta;
        aux.delta_dot = outer.delta_dot;
        aux.T_drive_total = outer.T_drive_total;
        aux.T_brake_total = outer.T_brake_total;
        aux.model_source = source;
    end

elseif strcmp(source, 'original')
    % Original built-in demo component path from the previous compact runner.
    [throttle, brake, delta, delta_dot] = vehicle14_command_inputs(t, par, args);
    IN0 = zeros(meta.nin,1);
    IN0(meta.input_index.delta) = delta;
    IN0(meta.input_index.delta_dot) = delta_dot;
    [~,~,Sdm] = mf_fun(q, u, par.P, IN0);
    S = full(Sdm);
    [IN, cdiag] = vehicle14_component_inputs_numeric(t, x, S, par, args, road, meta, throttle, brake, delta, delta_dot);
    [Mdm,Fdm,~] = mf_fun(q, u, par.P, IN);
    M = full(Mdm);
    F = full(Fdm);

    if nargout > 2
        aux.S = S;
        aux.IN = IN;
        aux.contact = cdiag;
        aux.throttle = throttle;
        aux.brake = brake;
        aux.delta = delta;
        aux.delta_dot = delta_dot;
        aux.T_drive_total = throttle * par.T_drive_max;
        aux.T_brake_total = brake * par.T_brake_max;
        aux.model_source = source;
    end
else
    error('Unknown args.model_source: %s. Use outer or original.', source);
end
end

function val = getfield_default_runner(s, name, default_val)
if isstruct(s) && isfield(s, name)
    val = s.(name);
else
    val = default_val;
end
end
function [IN, diag] = vehicle14_component_inputs_numeric(t_abs, x, S, par, args, road, meta, throttle, brake, delta, delta_dot)
%VEHICLE14_COMPONENT_INPUTS_NUMERIC Numeric 3D-ribbon tyre/contact layer outside Enak.
qv = x(1:14);
uv = x(15:28);
idx = meta.sensor_index;
in = meta.input_index;
IN = zeros(meta.nin,1);
IN(in.delta) = delta;
IN(in.delta_dot) = delta_dot;

T_drive = throttle * par.T_drive_max;
T_brake = brake * par.T_brake_max;
kb = par.brake_bias_front;
omRL = uv(13); omRR = uv(14);
diff_bias = par.diff_gain * par.diff_width * tanh((omRL - omRR)/par.diff_width);
wheel_torque.FL = -0.5 * kb * T_brake;
wheel_torque.FR = -0.5 * kb * T_brake;
wheel_torque.RL = +0.5 * T_drive - diff_bias - 0.5*(1-kb)*T_brake;
wheel_torque.RR = +0.5 * T_drive + diff_bias - 0.5*(1-kb)*T_brake;

names = {'FL','FR','RL','RR'};
diag.names = names;
for i = 1:4
    nm = names{i};
    p_wc_N = [S(idx.([nm '_Wc_x'])); S(idx.([nm '_Wc_y'])); S(idx.([nm '_Wc_z']))];
    v_wc_N = [S(idx.([nm '_Wc_vx'])); S(idx.([nm '_Wc_vy'])); S(idx.([nm '_Wc_vz']))];
    head_N = [S(idx.([nm '_head_Nx'])); S(idx.([nm '_head_Ny'])); S(idx.([nm '_head_Nz']))];
    spin_N = [S(idx.([nm '_spin_Nx'])); S(idx.([nm '_spin_Ny'])); S(idx.([nm '_spin_Nz']))];
    carrier_omega_N = [S(idx.([nm '_carrier_wx'])); S(idx.([nm '_carrier_wy'])); S(idx.([nm '_carrier_wz']))];
    omega = uv(10+i);
    ck = contact_kinematics_from_pose(par, road, p_wc_N, v_wc_N, head_N, spin_N, carrier_omega_N, omega);
    [F_B, M_B, d] = ribbon_road_tyre_contact_to_body_wrench(qv, par, road, p_wc_N, v_wc_N, head_N, spin_N, omega, carrier_omega_N, ck);
    IN(in.(['FBx' nm])) = F_B(1); IN(in.(['FBy' nm])) = F_B(2); IN(in.(['FBz' nm])) = F_B(3);
    IN(in.(['MBx' nm])) = M_B(1); IN(in.(['MBy' nm])) = M_B(2); IN(in.(['MBz' nm])) = M_B(3);
    IN(in.(['T' nm])) = wheel_torque.(nm);
    diag.Fx(i,1) = d.Fx; diag.Fy(i,1) = d.Fy; diag.Fz(i,1) = d.Fz;
    diag.kappa(i,1) = d.kappa; diag.alpha(i,1) = d.alpha; diag.compression(i,1) = d.compression;
    diag.road_z_at_wheel(i,1) = d.road_patch_z; diag.contact_x(i,1)=d.contact_x; diag.contact_y(i,1)=d.contact_y; diag.contact_z(i,1)=d.contact_z;
end

% Aero: zero during drop/settle; smooth placeholder aero during maneuver.
if t_abs < args.settle_time || args.disable_aero
    IN(in.FdragF) = 0; IN(in.FdragR) = 0; IN(in.FdownF) = 0; IN(in.FdownR) = 0;
else
    vx = uv(1); qdyn = 0.5*par.rho*vx*vx;
    [p_af_N,~] = body_point_state_num(qv, uv, [par.aero_front_x;0;par.aero_z]);
    [p_ar_N,~] = body_point_state_num(qv, uv, [par.aero_rear_x;0;par.aero_z]);
    [hF,~,~,~] = road_height_normal_at_point(road, p_af_N);
    [hR,~,~,~] = road_height_normal_at_point(road, p_ar_N);
    multF_raw = 1 + par.aero_h_sens_front*(par.aero_h_ref - hF);
    multR_raw = 1 + par.aero_h_sens_rear*(par.aero_h_ref - hR);
    multF = smooth_max(0.20, multF_raw, 1e-3);
    multR = smooth_max(0.20, multR_raw, 1e-3);
    D = 0.5*par.rho*par.CdA*vx*smooth_abs(vx,0.25);
    IN(in.FdragF) = -0.5*D; IN(in.FdragR) = -0.5*D;
    IN(in.FdownF) = -qdyn*par.ClA_front*multF;
    IN(in.FdownR) = -qdyn*par.ClA_rear*multR;
end
diag.throttle = throttle; diag.brake = brake; diag.delta = delta; diag.delta_dot = delta_dot;
end

function ck = contact_kinematics_from_pose(par, road, p_wc_N, v_wc_N, head_N, spin_N, carrier_omega_N, omega)
head_N = normalize_num(head_N, [1;0;0]); spin_N = normalize_num(spin_N, [0;1;0]);
geom0 = oriented_tyre_contact_geometry(par, road, p_wc_N, head_N, spin_N);
head_dot_N = cross(carrier_omega_N, head_N);
spin_dot_N = cross(carrier_omega_N, spin_N);
h = 1e-5;
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
if norm(t_long_raw) < 1e-10
    t_long_raw = wheel_heading_N - dot(wheel_heading_N,n_road_N)*n_road_N;
end
t_long_N = normalize_num(t_long_raw, e_long_N);
if dot(t_long_N, wheel_heading_N) < 0, t_long_N = -t_long_N; end
t_lat_N = normalize_num(cross(n_road_N, t_long_N), e_lat_N);
normal_in_radial_plane = n_road_N - dot(n_road_N, spin_axis_N)*spin_axis_N;
gamma = norm(normal_in_radial_plane);
if gamma < 1e-10
    gamma = 1.0; radial_down_N = -n_road_N; degenerate = true;
else
    radial_down_N = -normal_in_radial_plane/gamma; degenerate = false;
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
geom.s=s; geom.n=n; geom.p_road_N=p_road_N; geom.e_long_N=e_long_N; geom.e_lat_N=e_lat_N; geom.n_road_N=n_road_N;
geom.t_long_N=t_long_N; geom.t_lat_N=t_lat_N; geom.spin_axis_N=spin_axis_N; geom.radial_down_N=radial_down_N; geom.contact_to_wc_unit_N=contact_to_wc_unit_N;
geom.r_wc_to_contact_N=r_wc_to_contact_N; geom.p_contact_N=p_contact_N; geom.p_contact_road_N=p_contact_N; geom.r_wc_to_unloaded_N=r_wc_to_unloaded_N; geom.p_contact_unloaded_N=p_contact_unloaded_N;
geom.signed_center_distance=signed_center_distance; geom.radial_projection_norm=gamma; geom.radial_depth=radial_depth_normal; geom.radial_depth_normal=radial_depth_normal;
geom.distance_center_to_road_along_radial=distance_center_to_road_along_radial; geom.compression=compression_normal; geom.compression_normal=compression_normal; geom.compression_radial=compression_radial;
geom.contact_gap_normal=-compression_normal; geom.contact_gap_radial=-compression_radial; geom.road_plane_residual=road_plane_residual; geom.degenerate_contact=degenerate;
end

function [F_B, M_B, diag] = ribbon_road_tyre_contact_to_body_wrench(qv, par, road, p_wc_N, v_wc_N, wheel_heading_N, wheel_spin_axis_N, omega, carrier_omega_N, ck)
R_NB = body_rotation_from_q(qv); R_BN = R_NB.';
geom = oriented_tyre_contact_geometry(par, road, p_wc_N, wheel_heading_N, wheel_spin_axis_N);
n_road_N = geom.n_road_N; t_long_N = geom.t_long_N; t_lat_N = geom.t_lat_N; spin_axis_N = geom.spin_axis_N; r_wc_to_contact_N = geom.r_wc_to_contact_N;
compression = geom.compression;
compression_rate = ck.compression_rate_normal;
compression_rate_radial = ck.compression_rate_radial;
v_contact_geom_N = ck.contact_point_velocity_geom_N;
v_contact_carrier_N = ck.contact_patch_velocity_carrier_N;
v_spin_contact_N = ck.contact_patch_velocity_spin_N;
v_contact_material_N = ck.contact_patch_velocity_material_N;
Fz = smoothplus(par.kt*compression + par.ct*compression_rate, par.normal_smooth_eps);
Vx = dot(v_contact_carrier_N, t_long_N); Vy = dot(v_contact_carrier_N, t_lat_N);
Vx_material = dot(v_contact_material_N, t_long_N); Vy_material = dot(v_contact_material_N, t_lat_N);
rolling_speed = -dot(v_spin_contact_N, t_long_N);
Vden = sqrt(Vx*Vx + par.v_eps*par.v_eps);
kappa = (rolling_speed - Vx)/Vden;
alpha = atan2(Vy, Vden);
Fx0 = par.Ck*kappa; Fy0 = -par.Ca*alpha;
limit = par.mu*Fz + 1e-9;
norm_tan = sqrt(Fx0*Fx0 + Fy0*Fy0 + 1e-12);
scale = (1 + (norm_tan/limit)^4)^(-0.25);
Fx = Fx0*scale; Fy = Fy0*scale;
F_contact_N = Fx*t_long_N + Fy*t_lat_N + Fz*n_road_N;
F_B = R_BN*F_contact_N;
r_contact_B = R_BN*r_wc_to_contact_N;
M_B = cross(r_contact_B, F_B);
diag.s=geom.s; diag.n=geom.n; diag.compression=compression; diag.compression_rate=compression_rate; diag.compression_rate_radial=compression_rate_radial;
diag.Fz=Fz; diag.Fx=Fx; diag.Fy=Fy; diag.kappa=kappa; diag.alpha=alpha; diag.Vx_contact_carrier=Vx; diag.Vy_contact_carrier=Vy; diag.Vx_contact_material=Vx_material; diag.Vy_contact_material=Vy_material; diag.rolling_speed=rolling_speed;
diag.contact_x=geom.p_contact_N(1); diag.contact_y=geom.p_contact_N(2); diag.contact_z=geom.p_contact_N(3); diag.road_patch_x=geom.p_road_N(1); diag.road_patch_y=geom.p_road_N(2); diag.road_patch_z=geom.p_road_N(3);
diag.road_normal_velocity_residual=ck.road_normal_velocity_residual; diag.road_plane_residual_rate=ck.road_plane_residual_rate; diag.contact_vx=v_contact_geom_N(1); diag.contact_vy=v_contact_geom_N(2); diag.contact_vz=v_contact_geom_N(3);
end

function [h,e_z,s,n] = road_height_normal_at_point(road,p_N)
[s,n,p_road,~,~,e_z] = road_project_point(road,p_N);
h = dot(p_N-p_road,e_z);
end

function [pN,vN] = body_point_state_num(qv, uv, r_B)
R = body_rotation_from_q(qv); pP = qv(1:3); vP = R*uv(1:3); omega = R*uv(4:6); rN = R*r_B;
pN = pP + rN; vN = vP + cross(omega,rN);
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
if n < 1e-10
    y = fallback/norm(fallback);
else
    y = v/n;
end
end

function [s,n,p_road,e_long,e_lat,e_z] = road_project_point(road,p)
s = p(1);
for k=1:max(0,road.projection_iterations)
    [f,~,~,~,~,~,~] = projection_residual_given_s(road,p,s);
    if abs(f) < road.projection_tol, break; end
    ds=max(1e-4,1e-6*(1+abs(s)));
    fp=projection_residual_given_s(road,p,s+ds); fm=projection_residual_given_s(road,p,s-ds);
    dfds=(fp-fm)/(2*ds);
    if ~isfinite(dfds) || abs(dfds)<1e-12, break; end
    step=min(max(f/dfds,-5),5); s=s-step; if abs(step)<road.projection_tol, break; end
end
[~,n,p_road,e_long,e_lat,e_z,~] = projection_residual_given_s(road,p,s);
end

function [f,n,p_road,e_long,e_lat,e_z,p_s] = projection_residual_given_s(road,p,s)
[C,~,~,e_y,~,~] = centerline_quantities(road,s);
n = dot(p-C,e_y);
[p_road,e_long,e_lat,e_z,p_s,~] = road_surface_at(road,s,n);
r = p-p_road; f=dot(r,p_s);
end

function [p_road,e_long,e_lat,e_z,p_s,p_n] = road_surface_at(road,s,n)
[C,Cs,ex,ey,dey,ezc] = centerline_quantities(road,s);
p_road = C + n*ey; p_s = Cs + n*dey; p_n = ey;
e_long=normalize_num(p_s,ex); e_lat=normalize_num(p_n,ey); e_z=normalize_num(cross(e_long,e_lat),ezc);
if e_z(3)<0, e_z=-e_z; e_lat=-e_lat; end
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
ke=2*pi/road.elev_wavelength; kb=2*pi/road.bank_wavelength;
zc=road.elev_amp*(1-cos(ke*s)); dz=road.elev_amp*ke*sin(ke*s); d2z=road.elev_amp*ke*ke*cos(ke*s);
bank=road.bank_amp*sin(kb*s); db=road.bank_amp*kb*cos(kb*s);
C=[s;0;zc]; Cs=[1;0;dz]; L=sqrt(1+dz*dz); ex=[1/L;0;dz/L]; ey0=[0;1;0]; ez0=[-dz/L;0;1/L];
dez=[-d2z/(L^3);0;-dz*d2z/(L^3)]; cb=cos(bank); sb=sin(bank); ey=cb*ey0+sb*ez0; dey=-sb*db*ey0+cb*db*ez0+sb*dez; ezc=normalize_num(cross(ex,ey),[0;0;1]);
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_rhs_full_numeric.m
% -------------------------------------------------------------------------
function dx = vehicle14_rhs_full_numeric(t, x, mf_fun, par, args, road, meta)
%VEHICLE14_RHS_FULL_NUMERIC Explicit RHS xdot = M\F for MATLAB explicit/stiff ODE solvers.
[M,F] = vehicle14_core_numeric(t, x, mf_fun, par, args, road, meta);
dx = M \ F;
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_residual_full_numeric.m
% -------------------------------------------------------------------------
function res = vehicle14_residual_full_numeric(t, x, xdot, mf_fun, par, args, road, meta)
%VEHICLE14_RESIDUAL_FULL_NUMERIC Implicit residual M(x)*xdot - F(x) for ode15i.
[M,F] = vehicle14_core_numeric(t, x, mf_fun, par, args, road, meta);
res = M*xdot(:) - F;
end

% -------------------------------------------------------------------------
% Local helper from original package file: vehicle14_save_plots_matlab.m
% -------------------------------------------------------------------------
function vehicle14_save_plots_matlab(t, x, mf_fun, par, args, road, meta, plot_dir)
%VEHICLE14_SAVE_PLOTS_MATLAB Generate diagnostic PNG plots from the MATLAB/CasADi solution.
if nargin < 8 || isempty(plot_dir), plot_dir = args.plot_dir; end
if ~exist(plot_dir, 'dir'), mkdir(plot_dir); end
N = numel(t);
throttle=zeros(N,1); brake=zeros(N,1); delta=zeros(N,1); Fz=zeros(4,N); Fx=zeros(4,N); Fy=zeros(4,N); comp=zeros(4,N); roadz=zeros(4,N);
for k=1:N
    [~,~,aux] = vehicle14_core_numeric(t(k), x(k,:).', mf_fun, par, args, road, meta);
    throttle(k)=aux.throttle; brake(k)=aux.brake; delta(k)=aux.delta;
    Fz(:,k)=aux.contact.Fz; Fx(:,k)=aux.contact.Fx; Fy(:,k)=aux.contact.Fy; comp(:,k)=aux.contact.compression; roadz(:,k)=aux.contact.road_z_at_wheel;
end
names={'FL','FR','RL','RR'};
fig=figure('Visible','off'); plot3(x(:,1),x(:,2),x(:,3),'LineWidth',1.2); grid on; xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); title('Chassis CoM trajectory'); saveas(fig, fullfile(plot_dir,'trajectory_3d.png')); close(fig);
fig=figure('Visible','off'); plot(t, rad2deg(x(:,6)),'DisplayName','roll [deg]'); hold on; plot(t, rad2deg(x(:,5)),'DisplayName','pitch [deg]'); plot(t, x(:,15),'DisplayName','vx [m/s]'); plot(t, x(:,20),'DisplayName','yaw rate r [rad/s]'); grid on; xlabel('time [s]'); title(sprintf('Vehicle response channels (%s)', getfield_default_runner(args,'model_source','outer'))); legend show; saveas(fig, fullfile(plot_dir,'response_channels.png')); close(fig);
fig=figure('Visible','off'); plot(t, throttle,'DisplayName','throttle'); hold on; plot(t, brake,'DisplayName','brake'); plot(t, rad2deg(delta),'DisplayName','steering [deg]'); grid on; xlabel('time [s]'); title(sprintf('Driver commands (%s)', getfield_default_runner(args,'model_source','outer'))); legend show; saveas(fig, fullfile(plot_dir,'driver_commands.png')); close(fig);
fig=figure('Visible','off'); hold on; for i=1:4, plot(t,Fz(i,:),'DisplayName',names{i}); end; grid on; xlabel('time [s]'); ylabel('Fz [N]'); title('Wheel normal loads'); legend show; saveas(fig, fullfile(plot_dir,'wheel_normal_loads.png')); close(fig);
fig=figure('Visible','off'); hold on; for i=1:4, plot(t,comp(i,:),'DisplayName',names{i}); end; grid on; xlabel('time [s]'); ylabel('compression [m]'); title('Tyre compression'); legend show; saveas(fig, fullfile(plot_dir,'tyre_compression.png')); close(fig);
fig=figure('Visible','off'); plot(x(:,1),x(:,3),'DisplayName','CoM Z'); hold on; plot(x(:,1),mean(roadz,1),'DisplayName','mean road Z at wheels'); grid on; xlabel('X [m]'); ylabel('Z [m]'); title('3D ribbon road excitation'); legend show; saveas(fig, fullfile(plot_dir,'road_profile_response.png')); close(fig);
fig=figure('Visible','off'); plot(t, x(:,7:10)); grid on; xlabel('time [s]'); ylabel('z_i [m]'); title('Suspension/wheel-centre heave coordinates'); legend(names{:}); saveas(fig, fullfile(plot_dir,'suspension_travels.png')); close(fig);
fprintf('Saved MATLAB diagnostic plots to %s\n', plot_dir);
end

% -------------------------------------------------------------------------
% Local animation renderer, compact package version
% -------------------------------------------------------------------------
function render_vehicle14_matlab_animation(mat_file, out_mp4, out_gif)
%RENDER_VEHICLE14_MATLAB_ANIMATION Simple MATLAB MP4/GIF renderer for the saved solution.
if nargin < 1 || isempty(mat_file), mat_file = 'vehicle14_matlab_casadi_solution.mat'; end
if nargin < 2 || isempty(out_mp4), out_mp4 = 'vehicle14_matlab_casadi_animation.mp4'; end
if nargin < 3 || isempty(out_gif), out_gif = 'vehicle14_matlab_casadi_animation.gif'; end
load(mat_file, 't', 'x', 'par', 'args', 'road', 'meta');
[mf_fun, meta2] = vehicle14_build_casadi_function(); %#ok<ASGLU>
fps = 20; playback_speed = 1.0;
nFrames = min(450, max(2, ceil((t(end)-t(1))*fps/playback_speed)+1));
tFrames = linspace(t(1), t(end), nFrames);
xFrames = interp1(t, x, tFrames, 'linear');
fig = figure('Color','w','Position',[100 100 950 720]);
try
    v = VideoWriter(out_mp4, 'MPEG-4');
catch
    [p0,n0,~] = fileparts(out_mp4);
    out_mp4 = fullfile(p0, [n0 '.avi']);
    v = VideoWriter(out_mp4, 'Motion JPEG AVI');
end
v.FrameRate = fps; open(v);
for k = 1:nFrames
    clf(fig); ax = axes('Parent',fig); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); view(ax, 3);
    q = xFrames(k,1:14).'; u = xFrames(k,15:28).';
    xk = [q; u];
    source = lower(getfield_default_runner(args, 'model_source', 'outer'));
    if strcmp(source, 'outer')
        [IN0,~] = vehicle14_outer_models_user('control', tFrames(k), xk, [], par, args, road, meta);
    else
        [~,~,delta,delta_dot] = vehicle14_command_inputs(tFrames(k), par, args);
        IN0=zeros(meta.nin,1); IN0(meta.input_index.delta)=delta; IN0(meta.input_index.delta_dot)=delta_dot;
    end
    [~,~,Sdm] = mf_fun(q,u,par.P,IN0); S=full(Sdm); idx=meta.sensor_index;
    draw_local_road(ax, road, q(1), 10.0, road.width/2, 0.8, 1.2);
    draw_chassis(ax,q,par);
    names={'FL','FR','RL','RR'};
    for i=1:4
        nm=names{i};
        wc=[S(idx.([nm '_Wc_x'])); S(idx.([nm '_Wc_y'])); S(idx.([nm '_Wc_z']))];
        head=[S(idx.([nm '_head_Nx'])); S(idx.([nm '_head_Ny'])); S(idx.([nm '_head_Nz']))];
        spin=[S(idx.([nm '_spin_Nx'])); S(idx.([nm '_spin_Ny'])); S(idx.([nm '_spin_Nz']))];
        draw_wheel(ax,wc,head,spin,par.rw);
    end
    xlabel(ax,'X [m]'); ylabel(ax,'Y [m]'); zlabel(ax,'Z [m]');
    title(ax,sprintf('vehicle14 MATLAB/CasADi | t = %.2f s | roll %.2f deg | pitch %.2f deg', tFrames(k), rad2deg(q(6)), rad2deg(q(5))));
    xlim(ax,[q(1)-4 q(1)+8]); ylim(ax,[q(2)-4 q(2)+4]); zlim(ax,[q(3)-1.0 q(3)+1.8]);
    frame = getframe(fig); writeVideo(v,frame);
    [A,map] = rgb2ind(frame2im(frame),256);
    if k == 1
        imwrite(A,map,out_gif,'gif','LoopCount',Inf,'DelayTime',1/fps);
    else
        imwrite(A,map,out_gif,'gif','WriteMode','append','DelayTime',1/fps);
    end
end
close(v); close(fig);
fprintf('Saved animation: %s and %s\n', out_mp4, out_gif);
end

function draw_chassis(ax,q,par)
R=body_rotation_from_q(q); P=q(1:3); Lf=par.lf; Lr=par.lr; W=max(par.tf,par.tr)/2; H=0.18;
ptsB=[ Lf W H; Lf -W H; -Lr -W H; -Lr W H; Lf W -H; Lf -W -H; -Lr -W -H; -Lr W -H].';
pts=(P + R*ptsB).';
edges=[1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;3 7;4 8];
for e=1:size(edges,1), plot3(ax,pts(edges(e,:),1),pts(edges(e,:),2),pts(edges(e,:),3),'k-','LineWidth',1.6); end
end

function draw_wheel(ax,wc,head,spin,rw)
spin=normalize_num(spin,[0;1;0]); e1=head-dot(head,spin)*spin; e1=normalize_num(e1,[1;0;0]); e3=normalize_num(cross(spin,e1),[0;0;1]);
th=linspace(0,2*pi,40); C=wc + rw*(e1*cos(th)+e3*sin(th)); plot3(ax,C(1,:),C(2,:),C(3,:),'LineWidth',1.4);
plot3(ax,[wc(1)-0.12*spin(1), wc(1)+0.12*spin(1)], [wc(2)-0.12*spin(2), wc(2)+0.12*spin(2)], [wc(3)-0.12*spin(3), wc(3)+0.12*spin(3)], 'LineWidth',2.0);
end

function draw_local_road(ax,road,x0,xhalf,nhalf,ds,dn)
sVals=(x0-xhalf):ds:(x0+xhalf); nVals=(-nhalf):dn:nhalf; [SS,NN]=meshgrid(sVals,nVals); X=zeros(size(SS));Y=X;Z=X;
for i=1:numel(SS), [p,~,~,~]=road_surface_at_local(road,SS(i),NN(i)); X(i)=p(1);Y(i)=p(2);Z(i)=p(3); end
surf(ax,X,Y,Z,'FaceAlpha',0.25,'EdgeAlpha',0.25,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.4 0.4 0.4]);
end
% Uses local road_surface_at_local/centerline_quantities_local helpers defined in this runner.

