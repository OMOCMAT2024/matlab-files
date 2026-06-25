function [sim_ode15s, sim_fixed, sim_ode15i_dae, model, ctrl] = run_oc_simple_time_sim_arbitrary(profile, user_opts)
%RUN_OC_SIMPLE_TIME_SIM_ARBITRARY  Time simulation using your my_model.m only.
%
% This file deliberately does NOT run my_nlp.m, does NOT use result.mat, and
% does NOT use any optimal-control solution. It only reuses:
%   - my_params.m
%   - my_model.m
%   - the same scaled state/control/load definitions
%   - the same y_xdot time-domain derivative
%   - the same h_load tyre vertical-load algebraic equations
%
% Example:
%   addpath('D:\casadi_3_7_2')
%   addpath('D:\Applications\TOTPT-main\oc')
%   [sim_ode15s, sim_fixed] = run_oc_simple_time_sim_arbitrary('mild_sine');
%
% Optional short run:
%   opts = struct('t_end',3,'fixed_dt',1e-3,'do_plot',true);
%   [sim_ode15s, sim_fixed] = run_oc_simple_time_sim_arbitrary('mild_sine',opts);

if nargin < 1 || isempty(profile)
    profile = 'mild_sine';
end
if nargin < 2 || isempty(user_opts)
    user_opts = struct();
end

fprintf('\n=== Simple time simulation from my_model.m only ===\n');
fprintf('Profile: %s\n', profile);
fprintf('Simulator version: 2026-06-25-v8-command-actual-power-limit-no-clamps.\n');

model = build_model_functions_from_my_model();
veh = model.veh;
opts = default_options(user_opts, veh);
ctrl = make_arbitrary_controls(profile, opts, veh);

x0_scaled = make_initial_state_scaled(veh, model);
z0_static_scaled = static_load_guess_scaled(veh, model);

% For the reduced-ODE integrators, improve z0 by solving h_load = 0 once.
% For the true DAE integrator, use z0_static_scaled and let decic/ode15i
% determine a consistent algebraic initial condition without calling our
% tyre-load Newton solve.
u0_scaled = make_applied_control_scaled(ctrl.t(1), x0_scaled, ctrl, veh, model, opts);
kappa0 = eval_kappa(opts, 0, x0_scaled);
z0_reduced_scaled = solve_loads_newton(model, x0_scaled, u0_scaled, kappa0, z0_static_scaled, opts, true);

% Variable-step stiff integration using ode15s. This is an ODE in x only,
% with the tyre-load algebraic equations solved inside the RHS.
fprintf('Running ode15s variable-step stiff simulation...\n');
t_ode_start = tic;
z_cache = z0_reduced_scaled;
odefun = @(t,xs) rhs_scaled(t, xs, ctrl, veh, model, opts, z_cache);
% Use nested wrapper to update the z cache after each accepted/internal call.
    function dxs = rhs_cached(t, xs)
        [dxs, z_cache] = rhs_scaled(t, xs, ctrl, veh, model, opts, z_cache);
    end

tspan = 0:opts.output_dt:opts.t_end;
ode_opts = odeset('RelTol', opts.ode_reltol, 'AbsTol', opts.ode_abstol, ...
                  'MaxStep', opts.ode_max_step);
[t_ode, x_ode] = ode15s(@rhs_cached, tspan, x0_scaled, ode_opts);
time_ode15s = toc(t_ode_start);
fprintf('  ode15s wall time: %.3f s, output points: %d.\n', time_ode15s, numel(t_ode));
sim_ode15s = postprocess_solution('ode15s', t_ode(:), x_ode.', ctrl, veh, model, opts, z0_reduced_scaled);

% Fixed-step RK4 integration.
fprintf('Running fixed-step RK4 simulation... fixed_dt = %.4g s, stages/step = 4.\n', opts.fixed_dt);
t_fixed_start = tic;
sim_fixed = simulate_fixed_rk4(x0_scaled, z0_reduced_scaled, ctrl, veh, model, opts);

time_fixed = toc(t_fixed_start);
fprintf('  fixed RK4 wall time: %.3f s, steps: %d, nominal RHS calls: %d.\n', time_fixed, numel(sim_fixed.t)-1, 4*(numel(sim_fixed.t)-1));

% True DAE integration using ode15i.  Here the four tyre vertical loads are
% solver variables, not quantities solved by our RHS Newton reduction.  MATLAB
% will solve the coupled nonlinear DAE residual internally.
if opts.run_dae
    fprintf('Running ode15i true DAE simulation for [x; z]...\n');
    t_dae_start = tic;
    sim_ode15i_dae = simulate_ode15i_dae(x0_scaled, z0_static_scaled, ctrl, veh, model, opts);
    time_dae = toc(t_dae_start);
    fprintf('  ode15i DAE wall time: %.3f s, output points: %d.\n', time_dae, numel(sim_ode15i_dae.t));
else
    sim_ode15i_dae = [];
end

if isempty(sim_ode15i_dae)
    fprintf('Done. Max |h_load| ode15s = %.3e, fixed = %.3e\n', ...
        max(abs(sim_ode15s.h_load(:))), max(abs(sim_fixed.h_load(:))));
else
    fprintf('Done. Max |h_load| ode15s = %.3e, fixed = %.3e, ode15i DAE = %.3e\n', ...
        max(abs(sim_ode15s.h_load(:))), max(abs(sim_fixed.h_load(:))), max(abs(sim_ode15i_dae.h_load(:))));
end

if opts.do_plot
    plot_simple_results(sim_ode15s, sim_fixed, sim_ode15i_dae, ctrl, veh);
end
end

%% ------------------------------------------------------------------------
function model = build_model_functions_from_my_model()
% Build only single-output CasADi functions. No NLP, no collocation, no IPOPT.
import casadi.*
veh = my_params(); %#ok<NASGU>
run('my_model.m');

required = {'x','u','z','pv','y_xdot','h_load','y_acc','y_torque','y_mu','y_load', ...
            'x_s','u_s','z_s','x_min','x_max','z_min','z_max','nx','nz'};
for i = 1:numel(required)
    assert(exist(required{i}, 'var') == 1, 'my_model.m did not define %s.', required{i});
end

model = struct();
model.veh = veh;
model.nx = nx;
model.nz = nz;
model.x_s = full(x_s(:));
model.u_s = full(u_s(:));
model.z_s = full(z_s(:));
model.x_min = full(x_min(:));
model.x_max = full(x_max(:));
model.z_min = full(z_min(:));
model.z_max = full(z_max(:));

% Single-output functions only. This avoids the previous multi-output and
% concatenated-output issues completely.
model.f_xdot   = casadi.Function('sim_f_xdot',   {x,u,z,pv}, {y_xdot});
model.f_hload  = casadi.Function('sim_f_hload',  {x,u,z,pv}, {h_load});
model.f_hload_jac_z = casadi.Function('sim_f_hload_jac_z', {x,u,z,pv}, {jacobian(h_load, z)});
model.f_acc    = casadi.Function('sim_f_acc',    {x,u,z,pv}, {y_acc});
model.f_torque = casadi.Function('sim_f_torque', {x,u,z,pv}, {y_torque});
model.f_mu     = casadi.Function('sim_f_mu',     {x,u,z,pv}, {y_mu});
model.f_load   = casadi.Function('sim_f_load',   {x,u,z,pv}, {y_load});

fprintf('Built CasADi functions from my_model.m. nx = %d, nz = %d.\n', model.nx, model.nz);
end

%% ------------------------------------------------------------------------
function opts = default_options(user_opts, veh)
opts = struct();
opts.t_end = 15.0;
opts.output_dt = 0.02;
opts.fixed_dt = 1e-3;
opts.ode_reltol = 1e-6;
opts.ode_abstol = 1e-8;
opts.ode_max_step = 0.02;
opts.dae_reltol = 1e-6;
opts.dae_abstol = 1e-8;
opts.dae_max_step = 0.02;
opts.do_plot = true;
opts.run_dae = true;
% Fixed RK4 should not silently clamp states by default, because ode15s and
% ode15i do not impose the NLP state bounds during simulation.  Set this true
% only as a numerical guard for aggressive arbitrary inputs.
opts.enforce_state_bounds_fixed = false;
opts.kappa_constant = 0.0;       % straight road by default
opts.load_tol = 1e-8;
opts.load_max_iter = 12;
opts.load_fd_rel = 1e-6;
opts.warn_load_residual = 1e-5;
opts.min_power_omega = max(1.0, veh.omega_rear_axle_min);

names = fieldnames(user_opts);
for i = 1:numel(names)
    opts.(names{i}) = user_opts.(names{i});
end
end

%% ------------------------------------------------------------------------
% function ctrl = make_arbitrary_controls(profile, opts, veh)
% t = 0:opts.output_dt:opts.t_end;
% Tt_max = veh.Tt_max;
% 
% switch lower(profile)
%     case 'straight_accel'
%         Tt = 0.45*Tt_max*(1 - exp(-t/2.0));
%         delta = zeros(size(t));
%     case 'mild_sine'
%         Tt = 0.35*Tt_max*(1 - exp(-t/2.0)).*(1 + 0.08*sin(2*pi*0.20*t));
%         delta = deg2rad(3.0)*sin(2*pi*0.18*t).*(1 - exp(-t/1.0));
%     case 'slalom'
%         Tt = 0.25*Tt_max*(1 - exp(-t/1.5));
%         delta = deg2rad(6.0)*sin(2*pi*0.35*t).*(1 - exp(-t/1.0));
%     otherwise
%         error('Unknown profile "%s". Use straight_accel, mild_sine, or slalom.', profile);
% end
% 
% ctrl = struct();
% ctrl.t = t(:).';
% ctrl.Tt_cmd = min(max(Tt(:).', veh.Tt_min), veh.Tt_max);
% ctrl.delta_cmd = min(max(delta(:).', veh.delta_min), veh.delta_max);
% ctrl.Tb_cmd = zeros(size(ctrl.t)); % brake torque fixed zero as requested
% end









% function ctrl = make_arbitrary_controls(profile, opts, veh)
% % Generates open-loop control signals for vehicle simulations.
% % Integrates original legacy profiles with new multi-phase steering tests.
% 
% t = 0:opts.output_dt:opts.t_end;
% Tt_max = veh.Tt_max;
% 
% % Pre-allocate outputs with zeros to ensure clean initialization
% Tt = zeros(size(t));
% delta = zeros(size(t));
% 
% % --- Phased Profile Configuration Defaults ---
% % These variables read customizable fields from "opts". If they do not exist,
% % they gracefully fall back to safe default values.
% t_spinup = iffield(opts, 't_spinup', 3.0); % Time to spin up the vehicle speed [s]
% t_hold = iffield(opts, 't_hold', 2.0); % Time to maintain constant torque [s]
% Tt_target = iffield(opts, 'Tt_target', 0.25 * Tt_max); % Dynamic equilibrium target torque
% 
% % Compute absolute transition time markers
% t1 = t_spinup;
% t2 = t1 + t_hold;
% 
% switch lower(profile)
% % =========================================================================
% % LEGACY PROFILES (Maintained exactly as originally provided)
% % =========================================================================
% case 'straight_accel'
% Tt = 0.45*Tt_max*(1 - exp(-t/2.0));
% delta = zeros(size(t));
% 
% case 'mild_sine'
% Tt = 0.35*Tt_max*(1 - exp(-t/2.0)).*(1 + 0.08*sin(2*pi*0.20*t));
% delta = deg2rad(3.0)*sin(2*pi*0.18*t).*(1 - exp(-t/1.0));
% 
% case 'slalom'
% Tt = 0.25*Tt_max*(1 - exp(-t/1.5));
% delta = deg2rad(6.0)*sin(2*pi*0.35*t).*(1 - exp(-t/1.0));
% 
% % =========================================================================
% % NEW PHASED TEST PROFILES (Torque drop -> Steering activation)
% % =========================================================================
% case 'step_steer'
% % Setup common multi-phase torque profile
% Tt = compute_phased_torque(t, t1, t2, Tt_target);
% 
% % Extract custom steering step parameters
% mag_step = iffield(opts, 'step_mag', deg2rad(5.0)); % Step input magnitude [rad]
% dur_step = iffield(opts, 'step_dur', 2.0); % Duration to hold the step [s]
% 
% % Step is active from t2 up to (t2 + dur_step)
% idx_step = (t >= t2 & t < (t2 + dur_step));
% delta(idx_step) = mag_step;
% 
% case 'sine_steer'
% % Setup common multi-phase torque profile
% Tt = compute_phased_torque(t, t1, t2, Tt_target);
% 
% % Extract custom sinusoidal parameters
% mag_sine = iffield(opts, 'sine_mag', deg2rad(8.0)); % Peak amplitude [rad]
% freq_sine = iffield(opts, 'sine_freq', 0.5); % Frequency [Hz]
% tau_sine = iffield(opts, 'sine_tau', 0.4); % Steering ramp-in/steepness factor [s]
% 
% % Shift time index to start cleanly at t2 (the moment torque cuts to 0)
% idx_steer = (t >= t2);
% t_rel = t(idx_steer) - t2;
% 
% % Sinusoidal command with smooth exponential steepness curve
% delta(idx_steer) = mag_sine * sin(2*pi*freq_sine * t_rel) .* (1 - exp(-t_rel / tau_sine));
% 
% case 'linear_steer'
% % Setup common multi-phase torque profile
% Tt = compute_phased_torque(t, t1, t2, Tt_target);
% 
% % Extract custom linear ramp parameters
% slope_linear = iffield(opts, 'linear_slope', deg2rad(3.0)); % Ramp rate [rad/s]
% 
% % Linear slope starts incrementing exactly at t2
% idx_steer = (t >= t2);
% t_rel = t(idx_steer) - t2;
% delta(idx_steer) = slope_linear * t_rel;
% 
% % Strictly enforce vehicle max physical limits to avoid numerical issues
% delta = max(min(delta, veh.delta_max), veh.delta_min);
% 
% otherwise
% error('Unknown profile "%s". Use straight_accel, mild_sine, slalom, step_steer, sine_steer, or linear_steer.', profile);
% end
% 
% % Check total duration window for steering inputs (Capping total user execution if requested)
% if isfield(opts, 't_steer')
% idx_past_steer_window = (t >= (t2 + opts.t_steer));
% delta(idx_past_steer_window) = 0;
% end
% 
% % Assign back to the expected output control struct format
% % ctrl.Tt = Tt;
% % ctrl.delta = delta;
% 
% ctrl = struct();
% ctrl.t = t(:).';
% ctrl.Tt_cmd = min(max(Tt(:).', veh.Tt_min), veh.Tt_max);
% ctrl.delta_cmd = min(max(delta(:).', veh.delta_min), veh.delta_max);
% ctrl.Tb_cmd = zeros(size(ctrl.t)); % brake torque fixed zero as requested
% 
% end
% 
% % =========================================================================
% % LOCAL HELPER FUNCTIONS
% % =========================================================================
% 
% function Tt = compute_phased_torque(t, t1, t2, Tt_target)
% % Computes the vectorized phased torque sequence across 3 specific timelines:
% % Phase 1: Exponential spin up to target torque
% % Phase 2: Constant torque hold for dynamic equilibrium
% % Phase 3: Immediate torque cut to 0
% 
% Tt = zeros(size(t));
% 
% % Phase 1: Smooth exponential buildup up to t1
% % idx_p1 = (t >= 0 & t < t1);
% % tau_torque = t1 / 3; % Reaches roughly 95% of target by t1 limit
% % Tt(idx_p1) = Tt_target * (1 - exp(-t(idx_p1) / tau_torque));
% idx_p1 = (t >= 0 & t < t1);
% Tt(idx_p1) = Tt_target * (t(idx_p1) / t1);
% 
% % Phase 2: Constant torque boundary from t1 to t2
% idx_p2 = (t >= t1 & t < t2);
% Tt(idx_p2) = Tt_target;
% 
% % Phase 3: Immediate torque drop to 0 starting at t2
% idx_p3 = (t >= t2);
% Tt(idx_p3) = 0;
% end
% 
% function val = iffield(strct, field, default_val)
% % Safe check utility function to keep defaults if options structure fields are omitted
% if isfield(strct, field)
% val = strct.(field);
% else
% val = default_val;
% end
% end










function ctrl = make_arbitrary_controls(profile, opts, veh)
% Generates open-loop control signals for vehicle simulations.
% Integrates original legacy profiles with new multi-phase steering tests.
%
% New phased-test timing:
%
%   Phase 1: 0  -> t1
%           torque ramps up, steering = 0
%
%   Phase 2: t1 -> t2
%           torque held constant, steering = 0
%
%   Phase 3: t2 -> t_steer_start
%           torque = 0, steering = 0
%
%   Phase 4: t >= t_steer_start
%           torque = 0, steering profile starts
%
% where
%
%   t_steer_start = t2 + t_steer_delay
%
% You can set:
%
%   opts.t_steer_delay = 0.7;
%
% For backward compatibility, opts.delta_t is also accepted if
% opts.t_steer_delay is not provided.

t = 0:opts.output_dt:opts.t_end;
Tt_max = veh.Tt_max;

% Pre-allocate outputs with zeros to ensure clean initialization
Tt = zeros(size(t));
delta = zeros(size(t));

% -------------------------------------------------------------------------
% Phased Profile Configuration Defaults
% -------------------------------------------------------------------------
t_spinup = iffield(opts, 't_spinup', 3.0);              % [s]
t_hold = iffield(opts, 't_hold', 2.0);                  % [s]
Tt_target = iffield(opts, 'Tt_target', 0.25 * Tt_max);  % [Nm], or whatever unit your model uses

% New delay between torque cut and steering start
if isfield(opts, 't_steer_delay')
    t_steer_delay = opts.t_steer_delay;
elseif isfield(opts, 'delta_t')
    % Backward-compatible name, because you mentioned delta_t = 0.7 s.
    % However, t_steer_delay is clearer because delta_t often means time step.
    t_steer_delay = opts.delta_t;
else
    t_steer_delay = 0.7;
end

if ~isscalar(t_steer_delay) || ~isfinite(t_steer_delay) || t_steer_delay < 0
    error('opts.t_steer_delay must be a finite nonnegative scalar.');
end

% Compute absolute transition time markers
t1 = t_spinup;
t2 = t1 + t_hold;                 % torque becomes zero from here
t_steer_start = t2 + t_steer_delay; % steering starts from here

% Reference time used by optional steering-window cut-off
t_steer_window_ref = t2;

switch lower(profile)

% =========================================================================
% LEGACY PROFILES
% =========================================================================
case 'straight_accel'
    Tt = 0.45*Tt_max*(1 - exp(-t/2.0));
    delta = zeros(size(t));

case 'mild_sine'
    Tt = 0.35*Tt_max*(1 - exp(-t/2.0)).*(1 + 0.08*sin(2*pi*0.20*t));
    delta = deg2rad(3.0)*sin(2*pi*0.18*t).*(1 - exp(-t/1.0));

case 'slalom'
    Tt = 0.25*Tt_max*(1 - exp(-t/1.5));
    delta = deg2rad(6.0)*sin(2*pi*0.35*t).*(1 - exp(-t/1.0));

% =========================================================================
% PHASED TEST PROFILES
% Torque drop -> zero-input waiting period -> steering activation
% =========================================================================
case 'step_steer'
    % Torque:
    %   0 -> t1             ramp up
    %   t1 -> t2            hold
    %   t >= t2             zero
    Tt = compute_phased_torque(t, t1, t2, Tt_target);

    % Steering starts only after the zero-torque / zero-steering delay
    t_steer_window_ref = t_steer_start;

    mag_step = iffield(opts, 'step_mag', deg2rad(5.0)); % [rad]
    dur_step = iffield(opts, 'step_dur', 2.0);          % [s]

    idx_step = (t >= t_steer_start & t < (t_steer_start + dur_step));
    delta(idx_step) = mag_step;

case 'sine_steer'
    Tt = compute_phased_torque(t, t1, t2, Tt_target);

    % Steering starts only after the zero-torque / zero-steering delay
    t_steer_window_ref = t_steer_start;

    mag_sine = iffield(opts, 'sine_mag', deg2rad(8.0)); % [rad]
    freq_sine = iffield(opts, 'sine_freq', 0.5);        % [Hz]
    tau_sine = iffield(opts, 'sine_tau', 0.4);          % [s]

    idx_steer = (t >= t_steer_start);
    t_rel = t(idx_steer) - t_steer_start;

    delta(idx_steer) = mag_sine ...
        .* sin(2*pi*freq_sine * t_rel) ...
        .* (1 - exp(-t_rel / tau_sine));

case 'linear_steer'
    Tt = compute_phased_torque(t, t1, t2, Tt_target);

    % Steering starts only after the zero-torque / zero-steering delay
    t_steer_window_ref = t_steer_start;

    slope_linear = iffield(opts, 'linear_slope', deg2rad(3.0)); % [rad/s]

    idx_steer = (t >= t_steer_start);
    t_rel = t(idx_steer) - t_steer_start;

    delta(idx_steer) = slope_linear * t_rel;

    % Strictly enforce vehicle max physical limits to avoid numerical issues
    % delta = max(min(delta, veh.delta_max), veh.delta_min);
    % No steering command clamp here: ctrl.delta_cmd is the raw user command.

otherwise
    error(['Unknown profile "%s". Use straight_accel, mild_sine, slalom, ', ...
           'step_steer, sine_steer, or linear_steer.'], profile);
end

% -------------------------------------------------------------------------
% Optional steering duration cap
%
% For phased steering tests, opts.t_steer is now measured from the delayed
% steering start time, not from the torque-cut time.
%
% Example:
%   t2 = 5.0
%   opts.t_steer_delay = 0.7
%   opts.t_steer = 2.0
%
% Then:
%   torque becomes zero at 5.0 s
%   steering starts at 5.7 s
%   steering is forced back to zero at 7.7 s
% -------------------------------------------------------------------------
if isfield(opts, 't_steer')
    idx_past_steer_window = (t >= (t_steer_window_ref + opts.t_steer));
    delta(idx_past_steer_window) = 0;
end

% Assign back to the expected output control struct format
ctrl = struct();
ctrl.t = t(:).';
% ctrl.Tt_cmd = min(max(Tt(:).', veh.Tt_min), veh.Tt_max);
% ctrl.delta_cmd = min(max(delta(:).', veh.delta_min), veh.delta_max);
ctrl.Tt_cmd = Tt(:).';      % raw command, not clipped by NLP bounds
ctrl.delta_cmd = delta(:).'; % raw command, not clipped by NLP bounds
ctrl.Tb_cmd = zeros(size(ctrl.t)); % raw brake command; currently zero for these arbitrary profiles

end

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function Tt = compute_phased_torque(t, t1, t2, Tt_target)
% Computes the vectorized phased torque sequence:
%
% Phase 1: 0  <= t < t1   linear torque buildup
% Phase 2: t1 <= t < t2   constant torque hold
% Phase 3: t  >= t2       immediate torque cut to zero

Tt = zeros(size(t));

% Phase 1: Linear buildup up to t1
idx_p1 = (t >= 0 & t < t1);
Tt(idx_p1) = Tt_target * (t(idx_p1) / t1);

% Phase 2: Constant torque from t1 to t2
idx_p2 = (t >= t1 & t < t2);
Tt(idx_p2) = Tt_target;

% Phase 3: Immediate torque drop to zero starting at t2
idx_p3 = (t >= t2);
Tt(idx_p3) = 0;

end

function val = iffield(strct, field, default_val)
% Safe check utility function to keep defaults if options structure fields are omitted
if isfield(strct, field)
    val = strct.(field);
else
    val = default_val;
end
end











%% ------------------------------------------------------------------------
function x0_scaled = make_initial_state_scaled(veh, model)
ux0 = 5.0;
x0_phys = [ ...
    ux0;                    % ux
    0;                      % uy
    0;                      % r
    0;                      % phi
    0;                      % phidot
    0;                      % theta
    0;                      % thetadot
    0;                      % n
    0;                      % xi
    ux0/veh.rw;             % omega_fl
    ux0/veh.rw;             % omega_fr
    ux0/veh.rw;             % omega_rl
    ux0/veh.rw];            % omega_rr
x0_scaled = x0_phys(:) ./ model.x_s(:);
% x0_scaled = min(max(x0_scaled, model.x_min), model.x_max);
end

%% ------------------------------------------------------------------------
function z0_scaled = static_load_guess_scaled(veh, model)
Z_front = veh.m*veh.g*veh.lr/veh.l;
Z_rear  = veh.m*veh.g*veh.lf/veh.l;
z_phys = [0.5*Z_front; 0.5*Z_front; 0.5*Z_rear; 0.5*Z_rear];
z0_scaled = z_phys ./ model.z_s(:);
% z0_scaled = min(max(z0_scaled, model.z_min), model.z_max);
end

%% ------------------------------------------------------------------------
function [dxs, z_new] = rhs_scaled(t, xs, ctrl, veh, model, opts, z_guess)
xs = xs(:);
kappa = eval_kappa(opts, t, xs);
u_scaled = make_applied_control_scaled(t, xs, ctrl, veh, model, opts);
z_new = solve_loads_newton(model, xs, u_scaled, kappa, z_guess, opts, false);
dxs = full(model.f_xdot(xs, u_scaled, z_new, kappa));
dxs = dxs(:);
end

%% ------------------------------------------------------------------------
function u_cmd_phys = make_command_control_phys(t, ctrl)
% Interpolate the raw arbitrary command histories.  These are intentionally
% not clipped by the NLP control bounds; they are the commands the user asked
% the open-loop simulator to replay.
Tt_cmd = interp1(ctrl.t, ctrl.Tt_cmd, t, 'linear', 'extrap');
Tb_cmd = interp1(ctrl.t, ctrl.Tb_cmd, t, 'linear', 'extrap');
delta_cmd = interp1(ctrl.t, ctrl.delta_cmd, t, 'linear', 'extrap');
u_cmd_phys = [Tt_cmd; Tb_cmd; delta_cmd];
end

%% ------------------------------------------------------------------------
function [u_phys, u_diag] = make_applied_control_phys(t, xs, ctrl, veh, model, opts)
% Convert raw commands to the actual controls seen by my_model.m.
%
% Deliberate policy in this version:
%   - no state clamp here;
%   - no steering clamp here;
%   - no brake clamp here;
%   - no lower/upper Tt clamp except the requested drive-power/torque cap.
%
% The NLP power constraint is
%   Tt*((omega_rl + omega_rr)/2) <= veh.power_motor_max
% for the usual forward-driving case where rear wheel speeds are positive.
% Therefore the simulation applies
%   Tt_actual = min(Tt_cmd, Td_max, Pmax/omega_rear_avg)
% with only a small denominator guard to avoid division by zero at very low
% wheel speed.  In this model Td_max is veh.Tt_max because Td_state is not an
% active control/state in the uploaded my_model.m.
u_cmd_phys = make_command_control_phys(t, ctrl);
Tt_cmd = u_cmd_phys(1);
Tb_cmd = u_cmd_phys(2);
delta_cmd = u_cmd_phys(3);

omega_rl = xs(12) * model.x_s(12);
omega_rr = xs(13) * model.x_s(13);
omega_rear = 0.5*(omega_rl + omega_rr);

% Denominator guard only.  Do not use abs(omega_rear): the NLP formula uses
% the signed rear-axle average, and normal forward driving keeps it positive.
omega_for_power = max(omega_rear, opts.min_power_omega);
Tt_power_max = veh.power_motor_max / omega_for_power;

Td_max = veh.Tt_max;
if isfield(veh, 'Td_state_max')
    Td_max = veh.Td_state_max;
end

% Old bounded-control policy, intentionally disabled:
% Tt = min(max(Tt_cmd, veh.Tt_min), min(veh.Tt_max, Tt_power_max));
% Tb = 0.0;
% delta = min(max(delta_cmd, veh.delta_min), veh.delta_max);

Tt = min([Tt_cmd; Td_max; Tt_power_max]);
Tb = Tb_cmd;
delta = delta_cmd;

u_phys = [Tt; Tb; delta];

if nargout > 1
    u_diag = struct();
    u_diag.u_cmd_phys = u_cmd_phys;
    u_diag.omega_rl = omega_rl;
    u_diag.omega_rr = omega_rr;
    u_diag.omega_rear = omega_rear;
    u_diag.omega_for_power = omega_for_power;
    u_diag.Tt_power_max = Tt_power_max;
    u_diag.Td_max = Td_max;
    u_diag.Tt_upper_bound = min(Td_max, Tt_power_max);
    u_diag.power_rear = Tt * omega_rear;
    u_diag.power_motor_max = veh.power_motor_max;
end
end

%% ------------------------------------------------------------------------
function u_scaled = make_applied_control_scaled(t, xs, ctrl, veh, model, opts)
u_phys = make_applied_control_phys(t, xs, ctrl, veh, model, opts);
u_scaled = u_phys ./ model.u_s(:);
end

%% ------------------------------------------------------------------------
function kappa = eval_kappa(opts, t, xs) %#ok<INUSD>
% The simple simulator uses a straight road by default. You can replace this
% with a function handle through opts.kappa_fun if needed.
if isfield(opts, 'kappa_fun') && ~isempty(opts.kappa_fun)
    kappa = opts.kappa_fun(t, xs);
else
    kappa = opts.kappa_constant;
end
end

%% ------------------------------------------------------------------------
function z = solve_loads_newton(model, xs, us, kappa, z_guess, opts, do_warn)
% Solve the four lifted tyre-load equations h_load(x,u,z,kappa)=0.
%
% This is the core DAE reduction used by both integrators.  z is not
% integrated as a state; at every RHS/stage evaluation it is recomputed so
% that the exact same Chapter-9 load residuals from my_model.m are satisfied.
%
% Version v2 uses the analytic CasADi Jacobian dh_load/dz instead of finite
% differencing the four algebraic variables.  This is much faster for RK4,
% because RK4 evaluates the RHS four times per fixed step.
z = z_guess(:);
% z = min(max(z, model.z_min), model.z_max);

r = full(model.f_hload(xs, us, z, kappa));
r = r(:);
if norm(r, inf) < opts.load_tol
    return;
end

for iter = 1:opts.load_max_iter
    r = full(model.f_hload(xs, us, z, kappa));
    r = r(:);
    rnorm = norm(r, inf);
    if rnorm < opts.load_tol
        return;
    end

    J = full(model.f_hload_jac_z(xs, us, z, kappa));
    if any(~isfinite(J(:))) || any(~isfinite(r(:)))
        break;
    end

    % Four-by-four Newton system.  Use pinv only as a last-resort guard.
    if rcond(J) < 1e-12
        step = -pinv(J) * r;
    else
        step = -J \ r;
    end
    if any(~isfinite(step))
        break;
    end

    % Damped Newton with simple backtracking and z-bounds preservation.
    alpha = 1.0;
    accepted = false;
    for bt = 1:8
        z_trial = z + alpha*step;
        % z_trial = min(max(z_trial, model.z_min), model.z_max);
        r_trial = full(model.f_hload(xs, us, z_trial, kappa));
        r_trial = r_trial(:);
        if norm(r_trial, inf) <= 0.9*rnorm || norm(r_trial, inf) < opts.load_tol
            z = z_trial;
            accepted = true;
            break;
        end
        alpha = 0.5*alpha;
    end

    if ~accepted
        % Do not take a full bad Newton step; take a guarded small step.
        % z = min(max(z + 0.10*step, model.z_min), model.z_max);
        z = z + 0.10*step;
    end
end

r = full(model.f_hload(xs, us, z, kappa));
rnorm = norm(r(:), inf);
if do_warn || rnorm > opts.warn_load_residual
    warning('Tyre-load Newton residual max |h_load| = %.3e.', rnorm);
end
end

%% ------------------------------------------------------------------------
function sim = simulate_fixed_rk4(x0_scaled, z0_scaled, ctrl, veh, model, opts)
dt = opts.fixed_dt;
t = 0:dt:opts.t_end;
X = zeros(model.nx, numel(t));
Z = zeros(model.nz, numel(t));
X(:,1) = x0_scaled(:);
Z(:,1) = z0_scaled(:);
z_guess = z0_scaled(:);

for k = 1:numel(t)-1
    tk = t(k);
    xk = X(:,k);

    [k1, z1] = rhs_scaled(tk,          xk,             ctrl, veh, model, opts, z_guess);
    [k2, z2] = rhs_scaled(tk + 0.5*dt, xk + 0.5*dt*k1, ctrl, veh, model, opts, z1);
    [k3, z3] = rhs_scaled(tk + 0.5*dt, xk + 0.5*dt*k2, ctrl, veh, model, opts, z2);
    [k4, z4] = rhs_scaled(tk + dt,     xk + dt*k3,     ctrl, veh, model, opts, z3);

    X(:,k+1) = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    % if opts.enforce_state_bounds_fixed
    %     X(:,k+1) = min(max(X(:,k+1), model.x_min), model.x_max);
    % end
    Z(:,k+1) = z4;
    z_guess = z4;
end

sim = postprocess_solution('fixed_rk4', t(:), X, ctrl, veh, model, opts, z0_scaled);
end

%% ------------------------------------------------------------------------
function sim = simulate_ode15i_dae(x0_scaled, z0_scaled, ctrl, veh, model, opts)
%SIMULATE_ODE15I_DAE  True index-1 DAE integration of [x; z].
%
% Residual form for ode15i:
%   0 = xdot - f_xdot(x,u,z,kappa)
%   0 = h_load(x,u,z,kappa)
%
% Unlike ode15s/RK4 above, this function does not call solve_loads_newton in
% the RHS.  The four lifted tyre loads z are part of the solver state vector.
% MATLAB's DAE solver handles the nonlinear algebraic coupling internally.
% A consistent initial condition is found by decic, with x(0) fixed and z(0)
% allowed to move from the static-load guess.

t0 = 0.0;
y0_guess  = [x0_scaled(:); z0_scaled(:)];
yp0_guess = zeros(model.nx + model.nz, 1);

% Give decic a good derivative guess for the differential states.  This does
% not solve h_load manually; it only improves decic's starting point.
u0 = make_applied_control_scaled(t0, x0_scaled, ctrl, veh, model, opts);
k0 = eval_kappa(opts, t0, x0_scaled);
try
    yp0_guess(1:model.nx) = full(model.f_xdot(x0_scaled(:), u0, z0_scaled(:), k0));
catch
    yp0_guess(1:model.nx) = 0;
end

fixed_y0 = [true(model.nx,1); false(model.nz,1)];   % keep physical initial x fixed; let loads become consistent
fixed_yp0 = [false(model.nx,1); true(model.nz,1)];  % z_dot is absent from residual, so fix it to zero

resfun = @(t,y,yp) dae_residual_scaled(t, y, yp, ctrl, veh, model, opts);
decic_opts = odeset('RelTol', opts.dae_reltol, 'AbsTol', opts.dae_abstol);
[y0_cons, yp0_cons] = decic(resfun, t0, y0_guess, fixed_y0, yp0_guess, fixed_yp0, decic_opts);
res0_cons = resfun(t0, y0_cons, yp0_cons);
res0_cons_norm = norm(res0_cons, inf);
if res0_cons_norm > 100*opts.dae_abstol
    warning('decic returned a relatively large initial DAE residual: %.3e.', res0_cons_norm);
end

% Integrate the DAE.  ode15i uses variable steps internally; tspan below only
% fixes the requested output times for easier comparison with the other runs.
tspan = 0:opts.output_dt:opts.t_end;
dae_opts = odeset('RelTol', opts.dae_reltol, 'AbsTol', opts.dae_abstol, ...
                  'MaxStep', opts.dae_max_step);
[t_dae, y_dae] = ode15i(resfun, tspan, y0_cons, yp0_cons, dae_opts);

Y = y_dae.';
X = Y(1:model.nx,:);
Z = Y(model.nx+1:model.nx+model.nz,:);

% Do not re-solve tyre loads here.  Use the algebraic variables returned by
% ode15i so that the diagnostic residual reflects the true DAE solution.
sim = postprocess_dae_solution('ode15i_DAE', t_dae(:), X, Z, ctrl, veh, model, opts);

% Useful diagnostics for users comparing the reduced-ODE and true-DAE paths.
sim.y0_guess = y0_guess;
sim.y0_consistent = y0_cons;
sim.yp0_consistent = yp0_cons;
sim.decic_residual_inf = res0_cons_norm;
sim.initial_z_shift_scaled = y0_cons(model.nx+1:end) - z0_scaled(:);
sim.initial_z_shift_phys = sim.initial_z_shift_scaled .* model.z_s(:);
end

%% ------------------------------------------------------------------------
function res = dae_residual_scaled(t, y, yp, ctrl, veh, model, opts)
y = y(:);
yp = yp(:);
xs = y(1:model.nx);
zs = y(model.nx+1:model.nx+model.nz);
xsdot = yp(1:model.nx);

kap = eval_kappa(opts, t, xs);
us = make_applied_control_scaled(t, xs, ctrl, veh, model, opts);

f = full(model.f_xdot(xs, us, zs, kap));
h = full(model.f_hload(xs, us, zs, kap));

res = [xsdot(:) - f(:); h(:)];
end

%% ------------------------------------------------------------------------
function sim = postprocess_solution(name, t, X, ctrl, veh, model, opts, z0_scaled)
N = numel(t);
Z = zeros(model.nz, N);
U = zeros(3, N);
Ucmd = zeros(3, N);
H = zeros(model.nz, N);
ACC = zeros(7, N);
TORQUE = zeros(4, N);
MU = zeros(4, N);
LOAD = zeros(11, N);
omega_rear = zeros(1, N);
omega_for_power = zeros(1, N);
Tt_power_max = zeros(1, N);
Td_max_hist = zeros(1, N);
Tt_upper_bound = zeros(1, N);
power_rear = zeros(1, N);
z_guess = z0_scaled(:);

for k = 1:N
    xs = X(:,k);
    kap = eval_kappa(opts, t(k), xs);
    [u_phys, u_diag] = make_applied_control_phys(t(k), xs, ctrl, veh, model, opts);
    us = u_phys ./ model.u_s(:);
    zs = solve_loads_newton(model, xs, us, kap, z_guess, opts, false);
    z_guess = zs;

    Z(:,k) = zs;
    U(:,k) = us;
    Ucmd(:,k) = u_diag.u_cmd_phys;
    omega_rear(k) = u_diag.omega_rear;
    omega_for_power(k) = u_diag.omega_for_power;
    Tt_power_max(k) = u_diag.Tt_power_max;
    Td_max_hist(k) = u_diag.Td_max;
    Tt_upper_bound(k) = u_diag.Tt_upper_bound;
    power_rear(k) = u_diag.power_rear;
    H(:,k) = full(model.f_hload(xs, us, zs, kap));
    ACC(:,k) = full(model.f_acc(xs, us, zs, kap));
    TORQUE(:,k) = full(model.f_torque(xs, us, zs, kap));
    MU(:,k) = full(model.f_mu(xs, us, zs, kap));
    LOAD(:,k) = full(model.f_load(xs, us, zs, kap));
end

sim = struct();
sim.name = name;
sim.t = t(:).';
sim.x_scaled = X;
sim.x_phys = X .* model.x_s(:);
sim.z_scaled = Z;
sim.z_phys = Z .* model.z_s(:);
sim.u_scaled = U;
sim.u_phys = U .* model.u_s(:);
sim.u_cmd_phys = Ucmd;
sim.omega_rear = omega_rear;
sim.omega_for_power = omega_for_power;
sim.Tt_power_max = Tt_power_max;
sim.Td_max = Td_max_hist;
sim.Tt_upper_bound = Tt_upper_bound;
sim.power_rear = power_rear;
sim.h_load = H;
sim.acc = ACC;
sim.torque = TORQUE;
sim.mu = MU;
sim.load = LOAD;
end

%% ------------------------------------------------------------------------
function sim = postprocess_dae_solution(name, t, X, Z, ctrl, veh, model, opts)
N = numel(t);
U = zeros(3, N);
Ucmd = zeros(3, N);
H = zeros(model.nz, N);
ACC = zeros(7, N);
TORQUE = zeros(4, N);
MU = zeros(4, N);
LOAD = zeros(11, N);
omega_rear = zeros(1, N);
omega_for_power = zeros(1, N);
Tt_power_max = zeros(1, N);
Td_max_hist = zeros(1, N);
Tt_upper_bound = zeros(1, N);
power_rear = zeros(1, N);

for k = 1:N
    xs = X(:,k);
    zs = Z(:,k);
    kap = eval_kappa(opts, t(k), xs);
    [u_phys, u_diag] = make_applied_control_phys(t(k), xs, ctrl, veh, model, opts);
    us = u_phys ./ model.u_s(:);

    U(:,k) = us;
    Ucmd(:,k) = u_diag.u_cmd_phys;
    omega_rear(k) = u_diag.omega_rear;
    omega_for_power(k) = u_diag.omega_for_power;
    Tt_power_max(k) = u_diag.Tt_power_max;
    Td_max_hist(k) = u_diag.Td_max;
    Tt_upper_bound(k) = u_diag.Tt_upper_bound;
    power_rear(k) = u_diag.power_rear;
    H(:,k) = full(model.f_hload(xs, us, zs, kap));
    ACC(:,k) = full(model.f_acc(xs, us, zs, kap));
    TORQUE(:,k) = full(model.f_torque(xs, us, zs, kap));
    MU(:,k) = full(model.f_mu(xs, us, zs, kap));
    LOAD(:,k) = full(model.f_load(xs, us, zs, kap));
end

sim = struct();
sim.name = name;
sim.t = t(:).';
sim.x_scaled = X;
sim.x_phys = X .* model.x_s(:);
sim.z_scaled = Z;
sim.z_phys = Z .* model.z_s(:);
sim.u_scaled = U;
sim.u_phys = U .* model.u_s(:);
sim.u_cmd_phys = Ucmd;
sim.omega_rear = omega_rear;
sim.omega_for_power = omega_for_power;
sim.Tt_power_max = Tt_power_max;
sim.Td_max = Td_max_hist;
sim.Tt_upper_bound = Tt_upper_bound;
sim.power_rear = power_rear;
sim.h_load = H;
sim.acc = ACC;
sim.torque = TORQUE;
sim.mu = MU;
sim.load = LOAD;
end

%% ------------------------------------------------------------------------

function plot_command_vs_actual_controls(sim_ode15s, sim_fixed, sim_ode15i_dae, ctrl, veh)
% Plot raw command histories together with the actual controls used in the
% vehicle dynamics.  The only intentional command-vs-actual difference should
% be drive torque when the motor power/torque cap is active.
has_dae = ~isempty(sim_ode15i_dae);

figure('Name','Simple time simulation: command vs actual controls');
tiledlayout(4,1);

nexttile;
plot(ctrl.t, ctrl.Tt_cmd, 'k', 'LineWidth', 1.4); hold on; grid on;
plot(sim_ode15s.t, sim_ode15s.u_phys(1,:), 'LineWidth', 1.1);
plot(sim_fixed.t, sim_fixed.u_phys(1,:), '--', 'LineWidth', 1.0);
if has_dae, plot(sim_ode15i_dae.t, sim_ode15i_dae.u_phys(1,:), ':', 'LineWidth', 1.4); end
plot(sim_ode15s.t, sim_ode15s.Tt_upper_bound, '-.', 'LineWidth', 1.0);
ylabel('T_t [Nm]');
title('Drive torque: command vs actual applied to dynamics');
if has_dae
    legend('command','ode15s actual','RK4 actual','ode15i actual','min(Td max, Pmax/omega)', 'Location','best');
else
    legend('command','ode15s actual','RK4 actual','min(Td max, Pmax/omega)', 'Location','best');
end

nexttile;
plot(ctrl.t, ctrl.Tb_cmd, 'k', 'LineWidth', 1.4); hold on; grid on;
plot(sim_ode15s.t, sim_ode15s.u_phys(2,:), 'LineWidth', 1.1);
plot(sim_fixed.t, sim_fixed.u_phys(2,:), '--', 'LineWidth', 1.0);
if has_dae, plot(sim_ode15i_dae.t, sim_ode15i_dae.u_phys(2,:), ':', 'LineWidth', 1.4); end
ylabel('T_b [Nm]');
title('Brake torque: command vs actual applied to dynamics');

nexttile;
plot(ctrl.t, rad2deg(ctrl.delta_cmd), 'k', 'LineWidth', 1.4); hold on; grid on;
plot(sim_ode15s.t, rad2deg(sim_ode15s.u_phys(3,:)), 'LineWidth', 1.1);
plot(sim_fixed.t, rad2deg(sim_fixed.u_phys(3,:)), '--', 'LineWidth', 1.0);
if has_dae, plot(sim_ode15i_dae.t, rad2deg(sim_ode15i_dae.u_phys(3,:)), ':', 'LineWidth', 1.4); end
ylabel('\delta [deg]');
title('Steering: command vs actual applied to dynamics');

nexttile;
plot(sim_ode15s.t, sim_ode15s.power_rear/veh.power_motor_max, 'LineWidth', 1.1); hold on; grid on;
plot(sim_fixed.t, sim_fixed.power_rear/veh.power_motor_max, '--', 'LineWidth', 1.0);
if has_dae, plot(sim_ode15i_dae.t, sim_ode15i_dae.power_rear/veh.power_motor_max, ':', 'LineWidth', 1.4); end
yline(1.0, 'k-.', 'P_{max}');
xlabel('t [s]'); ylabel('T_t \omega_{rear}/P_{max} [-]');
title('Rear drive power usage after actual torque limiting');
end

%% ------------------------------------------------------------------------
function plot_simple_results(sim_ode15s, sim_fixed, sim_ode15i_dae, ctrl, veh)
has_dae = ~isempty(sim_ode15i_dae);

figure('Name','Simple time simulation: speed comparison');
plot(sim_ode15s.t, 3.6*sim_ode15s.x_phys(1,:), 'LineWidth', 1.2); hold on; grid on;
plot(sim_fixed.t, 3.6*sim_fixed.x_phys(1,:), '--', 'LineWidth', 1.0);
if has_dae
    plot(sim_ode15i_dae.t, 3.6*sim_ode15i_dae.x_phys(1,:), ':', 'LineWidth', 1.5);
    legend('ode15s reduced ODE','fixed RK4 reduced ODE','ode15i true DAE');
else
    legend('ode15s reduced ODE','fixed RK4 reduced ODE');
end
xlabel('t [s]'); ylabel('u_x [km/h]'); title('Longitudinal speed');

figure('Name','Simple time simulation: selected states comparison');
tiledlayout(3,1);
nexttile;
plot(sim_ode15s.t, 3.6*sim_ode15s.x_phys(1,:), 'LineWidth', 1.1); hold on; grid on;
plot(sim_fixed.t, 3.6*sim_fixed.x_phys(1,:), '--');
if has_dae, plot(sim_ode15i_dae.t, 3.6*sim_ode15i_dae.x_phys(1,:), ':', 'LineWidth', 1.4); end
xlabel('t [s]'); ylabel('u_x [km/h]'); title('Speed');
nexttile;
plot(sim_ode15s.t, rad2deg(sim_ode15s.x_phys(3,:)), 'LineWidth', 1.1); hold on; grid on;
plot(sim_fixed.t, rad2deg(sim_fixed.x_phys(3,:)), '--');
if has_dae, plot(sim_ode15i_dae.t, rad2deg(sim_ode15i_dae.x_phys(3,:)), ':', 'LineWidth', 1.4); end
xlabel('t [s]'); ylabel('r [deg/s]'); title('Yaw rate');
nexttile;
plot(sim_ode15s.t, rad2deg(sim_ode15s.x_phys(9,:)), 'LineWidth', 1.1); hold on; grid on;
plot(sim_fixed.t, rad2deg(sim_fixed.x_phys(9,:)), '--');
if has_dae, plot(sim_ode15i_dae.t, rad2deg(sim_ode15i_dae.x_phys(9,:)), ':', 'LineWidth', 1.4); end
xlabel('t [s]'); ylabel('\xi [deg]'); title('Course angle');
if has_dae
    legend('ode15s reduced ODE','fixed RK4 reduced ODE','ode15i true DAE', 'Location','best');
else
    legend('ode15s reduced ODE','fixed RK4 reduced ODE', 'Location','best');
end

% Old command-only plot, intentionally replaced by command-vs-actual plot.
% figure('Name','Simple time simulation: controls');
% plot(ctrl.t, ctrl.Tt_cmd, 'LineWidth', 1.2); hold on; grid on;
% yyaxis right; plot(ctrl.t, rad2deg(ctrl.delta_cmd), 'LineWidth', 1.2);
% xlabel('t [s]'); yyaxis left; ylabel('T_t command [Nm]'); yyaxis right; ylabel('\delta command [deg]');
% title('Arbitrary command histories');
plot_command_vs_actual_controls(sim_ode15s, sim_fixed, sim_ode15i_dae, ctrl, veh);

figure('Name','Simple time simulation: tyre loads comparison');
tiledlayout(2,2);
labels = {'FL','FR','RL','RR'};
for i = 1:4
    nexttile;
    plot(sim_ode15s.t, sim_ode15s.z_phys(i,:), 'LineWidth', 1.1); hold on; grid on;
    plot(sim_fixed.t, sim_fixed.z_phys(i,:), '--');
    if has_dae, plot(sim_ode15i_dae.t, sim_ode15i_dae.z_phys(i,:), ':', 'LineWidth', 1.4); end
    xlabel('t [s]'); ylabel('F_z [N]'); title(labels{i});
end
if has_dae
    legend('ode15s reduced ODE','fixed RK4 reduced ODE','ode15i true DAE', 'Location','best');
else
    legend('ode15s reduced ODE','fixed RK4 reduced ODE', 'Location','best');
end

figure('Name','Simple time simulation: load residual comparison');
semilogy(sim_ode15s.t, max(abs(sim_ode15s.h_load),[],1) + eps, 'LineWidth', 1.2); hold on; grid on;
semilogy(sim_fixed.t, max(abs(sim_fixed.h_load),[],1) + eps, '--', 'LineWidth', 1.0);
if has_dae
    semilogy(sim_ode15i_dae.t, max(abs(sim_ode15i_dae.h_load),[],1) + eps, ':', 'LineWidth', 1.5);
    legend('ode15s reduced ODE','fixed RK4 reduced ODE','ode15i true DAE');
else
    legend('ode15s reduced ODE','fixed RK4 reduced ODE');
end
xlabel('t [s]'); ylabel('max |h_{load}| [-]'); title('Tyre-load algebraic residual');

if has_dae
    % Compare the two reduced-ODE results to the DAE result at common DAE times.
    tx = sim_ode15i_dae.t;
    ux_ode = interp1(sim_ode15s.t, sim_ode15s.x_phys(1,:), tx, 'linear', 'extrap');
    ux_rk4 = interp1(sim_fixed.t, sim_fixed.x_phys(1,:), tx, 'linear', 'extrap');
    ux_dae = sim_ode15i_dae.x_phys(1,:);
    figure('Name','Simple time simulation: speed difference to DAE');
    plot(tx, 3.6*(ux_ode - ux_dae), 'LineWidth', 1.2); hold on; grid on;
    plot(tx, 3.6*(ux_rk4 - ux_dae), '--', 'LineWidth', 1.0);
    xlabel('t [s]'); ylabel('\Delta u_x [km/h]');
    legend('ode15s reduced ODE - ode15i DAE','fixed RK4 reduced ODE - ode15i DAE');
    title('Speed difference relative to true DAE solution');
end
end
