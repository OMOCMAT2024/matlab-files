function [sim_ode15s, sim_fixed, model, ctrl] = run_oc_simple_time_sim_arbitrary(profile, user_opts)
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
fprintf('Simulator version: 2026-06-23-v2-analytic-load-jacobian.\n');

model = build_model_functions_from_my_model();
veh = model.veh;
opts = default_options(user_opts, veh);
ctrl = make_arbitrary_controls(profile, opts, veh);

x0_scaled = make_initial_state_scaled(veh, model);
z0_scaled = static_load_guess_scaled(veh, model);

% Improve z0 by solving the same h_load = 0 equations at t = 0.
u0_scaled = make_applied_control_scaled(ctrl.t(1), x0_scaled, ctrl, veh, model, opts);
kappa0 = eval_kappa(opts, 0, x0_scaled);
z0_scaled = solve_loads_newton(model, x0_scaled, u0_scaled, kappa0, z0_scaled, opts, true);

% Variable-step stiff integration using ode15s. This is an ODE in x only,
% with the tyre-load algebraic equations solved inside the RHS.
fprintf('Running ode15s variable-step stiff simulation...\n');
t_ode_start = tic;
z_cache = z0_scaled;
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
sim_ode15s = postprocess_solution('ode15s', t_ode(:), x_ode.', ctrl, veh, model, opts, z0_scaled);

% Fixed-step RK4 integration.
fprintf('Running fixed-step RK4 simulation... fixed_dt = %.4g s, stages/step = 4.\n', opts.fixed_dt);
t_fixed_start = tic;
sim_fixed = simulate_fixed_rk4(x0_scaled, z0_scaled, ctrl, veh, model, opts);

time_fixed = toc(t_fixed_start);
fprintf('  fixed RK4 wall time: %.3f s, steps: %d, nominal RHS calls: %d.\n', time_fixed, numel(sim_fixed.t)-1, 4*(numel(sim_fixed.t)-1));

fprintf('Done. Max |h_load| ode15s = %.3e, fixed = %.3e\n', ...
    max(abs(sim_ode15s.h_load(:))), max(abs(sim_fixed.h_load(:))));

if opts.do_plot
    plot_simple_results(sim_ode15s, sim_fixed, ctrl, veh);
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
opts.do_plot = true;
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
function ctrl = make_arbitrary_controls(profile, opts, veh)
t = 0:opts.output_dt:opts.t_end;
Tt_max = veh.Tt_max;

switch lower(profile)
    case 'straight_accel'
        Tt = 0.45*Tt_max*(1 - exp(-t/2.0));
        delta = zeros(size(t));
    case 'mild_sine'
        Tt = 0.35*Tt_max*(1 - exp(-t/2.0)).*(1 + 0.08*sin(2*pi*0.20*t));
        delta = deg2rad(3.0)*sin(2*pi*0.18*t).*(1 - exp(-t/1.0));
    case 'slalom'
        Tt = 0.25*Tt_max*(1 - exp(-t/1.5));
        delta = deg2rad(6.0)*sin(2*pi*0.35*t).*(1 - exp(-t/1.0));
    otherwise
        error('Unknown profile "%s". Use straight_accel, mild_sine, or slalom.', profile);
end

ctrl = struct();
ctrl.t = t(:).';
ctrl.Tt_cmd = min(max(Tt(:).', veh.Tt_min), veh.Tt_max);
ctrl.delta_cmd = min(max(delta(:).', veh.delta_min), veh.delta_max);
ctrl.Tb_cmd = zeros(size(ctrl.t)); % brake torque fixed zero as requested
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
x0_scaled = min(max(x0_scaled, model.x_min), model.x_max);
end

%% ------------------------------------------------------------------------
function z0_scaled = static_load_guess_scaled(veh, model)
Z_front = veh.m*veh.g*veh.lr/veh.l;
Z_rear  = veh.m*veh.g*veh.lf/veh.l;
z_phys = [0.5*Z_front; 0.5*Z_front; 0.5*Z_rear; 0.5*Z_rear];
z0_scaled = z_phys ./ model.z_s(:);
z0_scaled = min(max(z0_scaled, model.z_min), model.z_max);
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
function u_scaled = make_applied_control_scaled(t, xs, ctrl, veh, model, opts)
Tt_cmd = interp1(ctrl.t, ctrl.Tt_cmd, t, 'linear', 'extrap');
delta_cmd = interp1(ctrl.t, ctrl.delta_cmd, t, 'linear', 'extrap');

omega_rl = xs(12) * model.x_s(12);
omega_rr = xs(13) * model.x_s(13);
omega_rear = 0.5*(omega_rl + omega_rr);
omega_for_power = max(abs(omega_rear), opts.min_power_omega);
Tt_power_max = veh.power_motor_max / omega_for_power;

Tt = min(max(Tt_cmd, veh.Tt_min), min(veh.Tt_max, Tt_power_max));
Tb = 0.0;
delta = min(max(delta_cmd, veh.delta_min), veh.delta_max);

u_phys = [Tt; Tb; delta];
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
z = min(max(z, model.z_min), model.z_max);

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
        z_trial = min(max(z_trial, model.z_min), model.z_max);
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
        z = min(max(z + 0.10*step, model.z_min), model.z_max);
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
    X(:,k+1) = min(max(X(:,k+1), model.x_min), model.x_max);
    Z(:,k+1) = z4;
    z_guess = z4;
end

sim = postprocess_solution('fixed_rk4', t(:), X, ctrl, veh, model, opts, z0_scaled);
end

%% ------------------------------------------------------------------------
function sim = postprocess_solution(name, t, X, ctrl, veh, model, opts, z0_scaled)
N = numel(t);
Z = zeros(model.nz, N);
U = zeros(3, N);
H = zeros(model.nz, N);
ACC = zeros(7, N);
TORQUE = zeros(4, N);
MU = zeros(4, N);
LOAD = zeros(11, N);
z_guess = z0_scaled(:);

for k = 1:N
    xs = X(:,k);
    kap = eval_kappa(opts, t(k), xs);
    us = make_applied_control_scaled(t(k), xs, ctrl, veh, model, opts);
    zs = solve_loads_newton(model, xs, us, kap, z_guess, opts, false);
    z_guess = zs;

    Z(:,k) = zs;
    U(:,k) = us;
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
sim.h_load = H;
sim.acc = ACC;
sim.torque = TORQUE;
sim.mu = MU;
sim.load = LOAD;
end

%% ------------------------------------------------------------------------
function plot_simple_results(sim_ode15s, sim_fixed, ctrl, veh)
figure('Name','Simple time simulation: speed');
plot(sim_ode15s.t, 3.6*sim_ode15s.x_phys(1,:), 'LineWidth', 1.2); hold on; grid on;
plot(sim_fixed.t, 3.6*sim_fixed.x_phys(1,:), '--');
xlabel('t [s]'); ylabel('u_x [km/h]'); legend('ode15s','fixed RK4'); title('Longitudinal speed');

figure('Name','Simple time simulation: controls');
plot(ctrl.t, ctrl.Tt_cmd, 'LineWidth', 1.2); hold on; grid on;
yyaxis right; plot(ctrl.t, rad2deg(ctrl.delta_cmd), 'LineWidth', 1.2);
xlabel('t [s]'); yyaxis left; ylabel('T_t command [Nm]'); yyaxis right; ylabel('\delta command [deg]');
title('Arbitrary command histories');

figure('Name','Simple time simulation: tyre loads');
plot(sim_ode15s.t, sim_ode15s.z_phys.', 'LineWidth', 1.0); grid on;
xlabel('t [s]'); ylabel('F_z [N]'); legend('FL','FR','RL','RR'); title('ode15s tyre vertical loads');

figure('Name','Simple time simulation: load residual');
semilogy(sim_ode15s.t, max(abs(sim_ode15s.h_load),[],1) + eps, 'LineWidth', 1.2); hold on; grid on;
semilogy(sim_fixed.t, max(abs(sim_fixed.h_load),[],1) + eps, '--');
xlabel('t [s]'); ylabel('max |h_{load}| [-]'); legend('ode15s','fixed RK4'); title('Tyre-load algebraic residual');
end
