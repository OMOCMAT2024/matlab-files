clear
close all
clc
format long

%% =========================================================================
% USER INPUTS
% =========================================================================
% Vehicle / gearbox
gear_ratios    = [5.250, 3.360, 2.172, 1.720, 1.316, 1.000, 0.822, 0.640];
final_drive    = 3.150;
engine_max_rpm = 7000;
axle_rpm_max   = 3500;

% Raw engine power data
rpm_samples = [1000, 2000, 3000, 4500, 5000, 6000, 7000];
kw_samples  = [26.2, 100.5, 169.6, 269.1, 293.2, 323.6, 322.5];

% Fit / smoothing settings
poly_degree              = 3;     % polynomial degree in engine-speed domain
k_branchblend            = 5.0;   % larger = sharper shifts (I used 5.0)
num_interpolation_points = 100000;

% Optional output controls
save_symbolic_text_files = true; % derivative expressions can be very large (false)
zoom_half_width          = 1.0;   % rad/s, for zoomed plots around each shift

%% =========================================================================
% BASIC CONSTANTS
% =========================================================================
rpm_to_rads  = pi / 30;
total_ratios = gear_ratios .* final_drive;

%% =========================================================================
% INPUT CHECKS
% =========================================================================
validateattributes(gear_ratios, {'double'}, {'real','finite','vector','positive'});
validateattributes(final_drive, {'double'}, {'real','finite','scalar','positive'});
validateattributes(engine_max_rpm, {'double'}, {'real','finite','scalar','positive'});
validateattributes(axle_rpm_max, {'double'}, {'real','finite','scalar','positive'});
validateattributes(rpm_samples, {'double'}, {'real','finite','vector','positive'});
validateattributes(kw_samples, {'double'}, {'real','finite','vector','numel',numel(rpm_samples)});
validateattributes(poly_degree, {'double'}, {'real','finite','scalar','integer','>=',1,'<',numel(rpm_samples)});
validateattributes(k_branchblend, {'double'}, {'real','finite','scalar','positive'});
validateattributes(num_interpolation_points, {'double'}, {'real','finite','scalar','integer','>=',10});
validateattributes(save_symbolic_text_files, {'logical','numeric'}, {'scalar'});
validateattributes(zoom_half_width, {'double'}, {'real','finite','scalar','positive'});

assert(all(diff(total_ratios) < 0), ...
    'total_ratios must be strictly decreasing from 1st gear to top gear.');
assert(all(diff(rpm_samples) > 0), ...
    'rpm_samples must be strictly increasing.');

%% =========================================================================
% NUMERICAL GRID
% =========================================================================
axle_rads_max  = axle_rpm_max * rpm_to_rads;
axle_rads_fine = linspace(0.1, axle_rads_max, num_interpolation_points);

%% =========================================================================
% BUILD AUTOMATIC SYMBOLIC + NUMERIC POWER MAP
% =========================================================================
model = build_power_map_from_raw_engine_data( ...
    rpm_samples, kw_samples, total_ratios, engine_max_rpm, poly_degree, k_branchblend);

% Final symbolic variable with the exact desired name
syms omega_rear_axle real
power_E_max  = subs(model.power_E_max_sym,  model.omega_a_sym, omega_rear_axle);
dpower_E_max = subs(model.dpower_E_max_sym, model.omega_a_sym, omega_rear_axle);
d2power_E_max = subs(model.d2power_E_max_sym, model.omega_a_sym, omega_rear_axle);

%% =========================================================================
% EXACT PIECEWISE 8-SPEED REFERENCE
% =========================================================================
[omega_e_piecewise, P_piecewise] = evaluate_exact_piecewise_power( ...
    axle_rads_fine, model.total_ratios, model.shifts_rads, model.engine_poly_coeffs);

%% =========================================================================
% FAST NUMERICAL EVALUATION OF SMOOTH MODEL AND DERIVATIVES
% =========================================================================
[P_smooth, omega_e_smooth, W_eval, ~] = evaluate_smooth_power_map(axle_rads_fine, model);

% Symbolic exact derivative evaluations
dP_smooth  = arrayfun(model.dpower_E_max_fun,  axle_rads_fine);
d2P_smooth = arrayfun(model.d2power_E_max_fun, axle_rads_fine);

%% =========================================================================
% INDEPENDENT FINITE-DIFFERENCE CHECKS
% =========================================================================
% Independent checks from the smooth power curve itself
dP_fd_from_P   = gradient(P_smooth, axle_rads_fine);
d2P_fd_from_P  = gradient(dP_fd_from_P, axle_rads_fine);

% Additional second-derivative check from dP_smooth itself
d2P_fd_from_dP = gradient(dP_smooth, axle_rads_fine);

% Ignore boundary points for FD comparisons
fd_idx = 3:(numel(axle_rads_fine)-2);

%% =========================================================================
% SUMMARY PRINTING
% =========================================================================
fprintf('\n============================================================\n');
fprintf('AUTOMATICALLY GENERATED POWER MAP SUMMARY\n');
fprintf('============================================================\n');

fprintf('\nEngine power polynomial fit in SI units:\n');
print_polynomial(model.engine_poly_coeffs, 'omega_e', 'P_engine');

fprintf('\nGear-specific rear-axle branch polynomials:\n');
for i = 1:numel(model.total_ratios)
    print_polynomial(model.branch_poly_coeffs(i,:), 'omega_a', sprintf('P_branch_%d', i));
end

fprintf('\nShift centers T_i [rad/s]:\n');
disp(model.shifts_rads(:).')

fprintf('k_branchblend = %.6g\n', model.k_branchblend);
fprintf('Max |sum(weights)-1| = %.3e\n', max(abs(sum(W_eval,1) - 1)));
fprintf('Min weight          = %.3e\n', min(W_eval(:)));
fprintf('Max weight          = %.3e\n', max(W_eval(:)));

fprintf('\nIndependent finite-difference consistency checks:\n');
fprintf('Max |dP_exact - FD(P)|         on interior = %.3e\n', ...
    max(abs(dP_smooth(fd_idx)  - dP_fd_from_P(fd_idx))));
fprintf('Max |d2P_exact - FD(FD(P))|    on interior = %.3e\n', ...
    max(abs(d2P_smooth(fd_idx) - d2P_fd_from_P(fd_idx))));
fprintf('Max |d2P_exact - FD(dP_exact)| on interior = %.3e\n', ...
    max(abs(d2P_smooth(fd_idx) - d2P_fd_from_dP(fd_idx))));

% Optional: report low-speed zero crossing of the raw polynomial fit
positive_real_roots = roots(model.engine_poly_coeffs);
positive_real_roots = positive_real_roots(abs(imag(positive_real_roots)) < 1e-10 & real(positive_real_roots) > 0);
positive_real_roots = sort(real(positive_real_roots));

if ~isempty(positive_real_roots)
    omega_e_zero = positive_real_roots(1);
    omega_a_zero_first = omega_e_zero / total_ratios(1);
    fprintf('\nLowest positive zero of fitted P_engine(omega_e):\n');
    fprintf('omega_e ~= %.6f rad/s  (%.3f rpm)\n', omega_e_zero, omega_e_zero * 30/pi);
    fprintf('Equivalent rear axle speed in 1st gear: omega_a ~= %.6f rad/s\n', omega_a_zero_first);
end

fprintf('\nThe symbolic expressions have been generated automatically.\n');
fprintf('Inspect them using:\n');
fprintf('  disp(power_E_max)\n');
fprintf('  disp(dpower_E_max)\n');
fprintf('  disp(d2power_E_max)\n');

%% =========================================================================
% OPTIONAL: SAVE SYMBOLIC EXPRESSIONS TO TEXT FILES
% =========================================================================
if save_symbolic_text_files
    expr_char = char(vpa(power_E_max, 12));
    fid = fopen('power_E_max_symbolic.txt', 'w');
    fprintf(fid, '%s\n', expr_char);
    fclose(fid);

    expr_char_d1 = char(vpa(dpower_E_max, 12));
    fid = fopen('dpower_E_max_symbolic.txt', 'w');
    fprintf(fid, '%s\n', expr_char_d1);
    fclose(fid);

    expr_char_d2 = char(vpa(d2power_E_max, 12));
    fid = fopen('d2power_E_max_symbolic.txt', 'w');
    fprintf(fid, '%s\n', expr_char_d2);
    fclose(fid);

    fprintf('\nSymbolic expressions saved to text files.\n');
end

%% =========================================================================
% ENGINE-SPEED-DOMAIN FIT PLOTS
% =========================================================================
omega_e_query = linspace(min(model.omega_engine_samples), max(model.omega_engine_samples), 1200);
P_engine_fit  = polyval(model.engine_poly_coeffs, omega_e_query);

figure('Color', 'w', 'Position', [80, 80, 1300, 520]);
tiledlayout(1,2, 'TileSpacing', 'Loose', 'Padding', 'Compact');

nexttile; hold on; grid on; box on;
plot(omega_e_query, P_engine_fit/1000, 'b', 'LineWidth', 2.5);
plot(model.omega_engine_samples, model.power_samples_W/1000, 'ko', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 7);
xlabel('Engine speed \omega_e (rad/s)');
ylabel('Power (kW)');
title('Raw engine data and fitted engine-power polynomial');
legend('Polynomial fit', 'Raw data', 'Location', 'best');

nexttile; hold on; grid on; box on;
torque_fit = P_engine_fit ./ omega_e_query;
torque_raw = model.power_samples_W ./ model.omega_engine_samples;
plot(omega_e_query, torque_fit, 'r', 'LineWidth', 2.5);
plot(model.omega_engine_samples, torque_raw, 'ks', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 7);
xlabel('Engine speed \omega_e (rad/s)');
ylabel('Torque (N·m)');
title('Torque implied by fitted power curve');
legend('Implied torque from fit', 'Raw-point torque', 'Location', 'best');

%% =========================================================================
% REAR-AXLE DOMAIN: EXACT PIECEWISE VS SMOOTH power_E_max
% =========================================================================
figure('Color', 'w', 'Position', [80, 80, 1400, 560]);
tiledlayout(1,2, 'TileSpacing', 'Loose', 'Padding', 'Compact');

nexttile; hold on; grid on; box on;
plot(axle_rads_fine, omega_e_smooth,    'b',   'LineWidth', 2.5);
plot(axle_rads_fine, omega_e_piecewise, 'k--', 'LineWidth', 1.2);
for i = 1:numel(model.shifts_rads)
    xline(model.shifts_rads(i), 'k:');
end
xlabel('Rear axle speed \omega_a (rad/s)');
ylabel('Engine speed \omega_e (rad/s)');
title('Engine speed map');
legend('Smooth branch-blend map', 'Exact piecewise 8-speed', 'Location', 'northwest');

nexttile; hold on; grid on; box on;
plot(axle_rads_fine, P_smooth/1000,     'g',   'LineWidth', 2.8);
plot(axle_rads_fine, P_piecewise/1000,  'k--', 'LineWidth', 1.2);
for i = 1:numel(model.shifts_rads)
    xline(model.shifts_rads(i), 'k:');
end
xlabel('Rear axle speed \omega_a (rad/s)');
ylabel('Engine power (kW)');
title('Automatically generated power\_E\_max(\omega_a)');
legend('Smooth automatically generated power\_E\_max', ...
       'Exact piecewise 8-speed reference', ...
       'Location', 'southeast');

%% =========================================================================
% FULL-RANGE PLOTS OF POWER AND ITS 1ST / 2ND DERIVATIVES
% =========================================================================
figure('Color', 'w', 'Position', [60, 60, 1450, 900]);
t_der = tiledlayout(3,1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Power
nexttile; hold on; grid on; box on;
plot(axle_rads_fine, P_smooth/1000, 'g', 'LineWidth', 2.5);
plot(axle_rads_fine, P_piecewise/1000, 'k--', 'LineWidth', 1.1);
for i = 1:numel(model.shifts_rads)
    xline(model.shifts_rads(i), 'k:');
end
xlabel('Rear axle speed \omega_a (rad/s)');
ylabel('Power (kW)');
title('Smooth power\_E\_max(\omega_a)');
legend('Smooth power\_E\_max', 'Exact piecewise 8-speed', 'Location', 'best');

% First derivative
nexttile; hold on; grid on; box on;
plot(axle_rads_fine, dP_smooth, 'b', 'LineWidth', 2.2);
plot(axle_rads_fine, dP_fd_from_P, 'r--', 'LineWidth', 1.0);
for i = 1:numel(model.shifts_rads)
    xline(model.shifts_rads(i), 'k:');
end
xlabel('Rear axle speed \omega_a (rad/s)');
ylabel('dP/d\omega_a  (W / (rad/s))');
title('First derivative of smooth power\_E\_max');
legend('Symbolic exact derivative', 'Finite-difference check from P', 'Location', 'best');

% Second derivative
nexttile; hold on; grid on; box on;
plot(axle_rads_fine, d2P_smooth, 'm', 'LineWidth', 2.2);
plot(axle_rads_fine, d2P_fd_from_P, 'r--', 'LineWidth', 1.0);
for i = 1:numel(model.shifts_rads)
    xline(model.shifts_rads(i), 'k:');
end
xlabel('Rear axle speed \omega_a (rad/s)');
ylabel('d^2P/d\omega_a^2  (W / (rad/s)^2)');
title('Second derivative of smooth power\_E\_max');
legend('Symbolic exact derivative', 'Finite-difference check from P', 'Location', 'best');

title(t_der, 'Smooth Power Map and Automatically Generated Derivatives', 'FontSize', 16);

%% =========================================================================
% ZOOM AROUND EACH SHIFT: POWER
% =========================================================================
nShifts = numel(model.shifts_rads);

figure('Color', 'w', 'Position', [90, 50, 1350, 250*nShifts]);
tiledlayout(nShifts, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for i = 1:nShifts
    Tshift = model.shifts_rads(i);

    nexttile; hold on; grid on; box on;
    plot(axle_rads_fine, P_smooth/1000,    'g',   'LineWidth', 2.8);
    plot(axle_rads_fine, P_piecewise/1000, 'k--', 'LineWidth', 1.2);
    xline(Tshift, 'k:', sprintf('shift %d -> %d', i, i+1), ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelOrientation', 'horizontal');

    xlim([max(axle_rads_fine(1), Tshift - zoom_half_width), ...
          min(axle_rads_fine(end), Tshift + zoom_half_width)]);

    xlabel('Rear axle speed \omega_a (rad/s)');
    ylabel('Power (kW)');
    title(sprintf('Zoom around shift %d -> %d', i, i+1));
    legend('Smooth power\_E\_max', 'Exact piecewise', 'Location', 'best');
end

%% =========================================================================
% ZOOMED DERIVATIVE PLOTS AROUND EACH SHIFT
% =========================================================================
figure('Color', 'w', 'Position', [70, 40, 1500, 320*nShifts]);
t_zoom_der = tiledlayout(nShifts, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for i = 1:nShifts
    Tshift = model.shifts_rads(i);
    xL = max(axle_rads_fine(1), Tshift - zoom_half_width);
    xU = min(axle_rads_fine(end), Tshift + zoom_half_width);

    % First derivative zoom
    nexttile; hold on; grid on; box on;
    plot(axle_rads_fine, dP_smooth, 'b', 'LineWidth', 2.2);
    plot(axle_rads_fine, dP_fd_from_P, 'r--', 'LineWidth', 1.0);
    xline(Tshift, 'k:', sprintf('shift %d -> %d', i, i+1), ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelOrientation', 'horizontal');
    xlim([xL, xU]);
    xlabel('Rear axle speed \omega_a (rad/s)');
    ylabel('dP/d\omega_a');
    title(sprintf('Zoomed 1st derivative near shift %d -> %d', i, i+1));
    legend('Symbolic exact derivative', 'FD check from P', 'Location', 'best');

    % Second derivative zoom
    nexttile; hold on; grid on; box on;
    plot(axle_rads_fine, d2P_smooth, 'm', 'LineWidth', 2.2);
    plot(axle_rads_fine, d2P_fd_from_P, 'r--', 'LineWidth', 1.0);
    xline(Tshift, 'k:', sprintf('shift %d -> %d', i, i+1), ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelOrientation', 'horizontal');
    xlim([xL, xU]);
    xlabel('Rear axle speed \omega_a (rad/s)');
    ylabel('d^2P/d\omega_a^2');
    title(sprintf('Zoomed 2nd derivative near shift %d -> %d', i, i+1));
    legend('Symbolic exact derivative', 'FD check from P', 'Location', 'best');
end

title(t_zoom_der, 'Zoomed Derivative Behaviour Around Shift Points', 'FontSize', 16);

%% =========================================================================
% OPTIONAL: SHOW GEAR WEIGHTS
% =========================================================================
figure('Color', 'w', 'Position', [100, 100, 1250, 550]);
hold on; grid on; box on;
for i = 1:size(W_eval,1)
    plot(axle_rads_fine, W_eval(i,:), 'LineWidth', 2.0);
end
for i = 1:numel(model.shifts_rads)
    xline(model.shifts_rads(i), 'k:');
end
xlabel('Rear axle speed \omega_a (rad/s)');
ylabel('Weight');
title('Smooth gear weights (partition of unity)');
legend(compose('w_%d', 1:size(W_eval,1)), 'Location', 'eastoutside');

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function model = build_power_map_from_raw_engine_data( ...
    rpm_samples, kw_samples, total_ratios, engine_max_rpm, poly_degree, k_branchblend)

    rpm_to_rads = pi / 30;

    % Raw engine data in SI units
    omega_engine_samples = rpm_samples * rpm_to_rads;
    power_samples_W      = kw_samples * 1000;

    % Polynomial fit in engine-speed domain
    engine_poly_coeffs = polyfit(omega_engine_samples, power_samples_W, poly_degree);

    % Shift locations in rear-axle-speed domain
    shifts_rads = (engine_max_rpm ./ total_ratios(1:end-1)) * rpm_to_rads;
    assert(all(diff(shifts_rads) > 0), 'Computed shift locations must be strictly increasing.');

    % Symbolics
    syms omega_a omega_e real
    P_engine_sym = poly2sym(engine_poly_coeffs, omega_e);

    nGears = numel(total_ratios);

    % Smooth logistic/tanh gate
    sig = @(z) sym(1)/2 * (1 + tanh(z/2));

    % Smooth cumulative shift gates
    s_shift_sym = sym(zeros(nGears-1,1));
    for j = 1:nGears-1
        s_shift_sym(j) = sig(k_branchblend * (omega_a - shifts_rads(j)));
    end

    % Partition-of-unity gear weights
    w_gear_sym = sym(zeros(nGears,1));
    w_gear_sym(1) = 1 - s_shift_sym(1);
    for i = 2:nGears-1
        w_gear_sym(i) = s_shift_sym(i-1) - s_shift_sym(i);
    end
    w_gear_sym(nGears) = s_shift_sym(end);

    % Gear-wise engine-speed and power branches
    omega_e_branch_sym = sym(zeros(nGears,1));
    P_branch_sym       = sym(zeros(nGears,1));
    branch_poly_coeffs = zeros(nGears, poly_degree+1);

    for i = 1:nGears
        omega_e_branch_sym(i) = total_ratios(i) * omega_a;
        P_branch_sym(i)       = subs(P_engine_sym, omega_e, omega_e_branch_sym(i));
        branch_poly_coeffs(i,:) = transform_poly_to_axle_domain(engine_poly_coeffs, total_ratios(i));
    end

    % Final smooth symbolic maps
    omega_e_smooth_sym = sum(w_gear_sym .* omega_e_branch_sym);
    power_E_max_sym    = sum(w_gear_sym .* P_branch_sym);

    % First and second derivatives
    dpower_E_max_sym  = diff(power_E_max_sym, omega_a);
    d2power_E_max_sym = diff(dpower_E_max_sym, omega_a);

    % Function handles
    omega_e_smooth_fun = matlabFunction(omega_e_smooth_sym, 'Vars', omega_a);
    power_E_max_fun    = matlabFunction(power_E_max_sym,    'Vars', omega_a);
    dpower_E_max_fun   = matlabFunction(dpower_E_max_sym,   'Vars', omega_a);
    d2power_E_max_fun  = matlabFunction(d2power_E_max_sym,  'Vars', omega_a);

    % Package
    model = struct();
    model.total_ratios         = total_ratios;
    model.shifts_rads          = shifts_rads;
    model.k_branchblend        = k_branchblend;

    model.omega_engine_samples = omega_engine_samples;
    model.power_samples_W      = power_samples_W;

    model.engine_poly_coeffs   = engine_poly_coeffs;
    model.branch_poly_coeffs   = branch_poly_coeffs;

    model.omega_a_sym          = omega_a;
    model.omega_e_sym          = omega_e;
    model.P_engine_sym         = P_engine_sym;
    model.s_shift_sym          = s_shift_sym;
    model.w_gear_sym           = w_gear_sym;
    model.omega_e_branch_sym   = omega_e_branch_sym;
    model.P_branch_sym         = P_branch_sym;

    model.omega_e_smooth_sym   = omega_e_smooth_sym;
    model.power_E_max_sym      = power_E_max_sym;
    model.dpower_E_max_sym     = dpower_E_max_sym;
    model.d2power_E_max_sym    = d2power_E_max_sym;

    model.omega_e_smooth_fun   = omega_e_smooth_fun;
    model.power_E_max_fun      = power_E_max_fun;
    model.dpower_E_max_fun     = dpower_E_max_fun;
    model.d2power_E_max_fun    = d2power_E_max_fun;
end

function coeffs_axle = transform_poly_to_axle_domain(coeffs_engine, ratio)
    % If
    %   P_engine(omega_e) = a_n*omega_e^n + ... + a_1*omega_e + a_0
    % and
    %   omega_e = ratio * omega_a
    % then
    %   P_branch(omega_a) = a_n*(ratio^n)*omega_a^n + ... + a_1*ratio*omega_a + a_0

    n = numel(coeffs_engine) - 1;
    coeffs_axle = zeros(size(coeffs_engine));

    for k = 1:numel(coeffs_engine)
        power_of_term = n - (k - 1);
        coeffs_axle(k) = coeffs_engine(k) * ratio^power_of_term;
    end
end

function [omega_e_piecewise, P_piecewise] = evaluate_exact_piecewise_power( ...
    omega_a_vec, total_ratios, shifts_rads, engine_poly_coeffs)

    omega_a_vec = omega_a_vec(:).';
    nGears = numel(total_ratios);

    omega_e_piecewise = zeros(size(omega_a_vec));
    P_piecewise       = zeros(size(omega_a_vec));

    bounds = [0, shifts_rads, inf];

    for i = 1:nGears
        idx = omega_a_vec >= bounds(i) & omega_a_vec < bounds(i+1);
        omega_e_piecewise(idx) = total_ratios(i) * omega_a_vec(idx);
        P_piecewise(idx)       = polyval(engine_poly_coeffs, omega_e_piecewise(idx));
    end
end

function [P_smooth, omega_e_smooth, W_eval, P_branch_eval] = evaluate_smooth_power_map(omega_a_vec, model)

    omega_a_vec = omega_a_vec(:).';
    nPts   = numel(omega_a_vec);
    nGears = numel(model.total_ratios);

    % Smooth cumulative gates
    S_eval = zeros(nGears-1, nPts);
    for j = 1:nGears-1
        S_eval(j,:) = 0.5 * (1 + tanh(model.k_branchblend * (omega_a_vec - model.shifts_rads(j)) / 2));
    end

    % Partition-of-unity weights
    W_eval = zeros(nGears, nPts);
    W_eval(1,:) = 1 - S_eval(1,:);
    for i = 2:nGears-1
        W_eval(i,:) = S_eval(i-1,:) - S_eval(i,:);
    end
    W_eval(nGears,:) = S_eval(end,:);

    % Engine-speed branches
    omega_e_branch = model.total_ratios(:) * omega_a_vec;
    omega_e_smooth = sum(W_eval .* omega_e_branch, 1);

    % Power branches evaluated directly in axle domain
    P_branch_eval = zeros(nGears, nPts);
    for i = 1:nGears
        P_branch_eval(i,:) = polyval(model.branch_poly_coeffs(i,:), omega_a_vec);
    end

    % Final smooth power map
    P_smooth = sum(W_eval .* P_branch_eval, 1);
end

function print_polynomial(coeffs, varname, funcname)
    n = numel(coeffs) - 1;

    fprintf('%s(%s) = ', funcname, varname);

    for k = 1:numel(coeffs)
        c = coeffs(k);
        p = n - (k - 1);

        if k == 1
            if c < 0
                fprintf('-');
                c = abs(c);
            end
        else
            if c >= 0
                fprintf(' + ');
            else
                fprintf(' - ');
                c = abs(c);
            end
        end

        if p > 1
            fprintf('%.12g*%s^%d', c, varname, p);
        elseif p == 1
            fprintf('%.12g*%s', c, varname);
        else
            fprintf('%.12g', c);
        end
    end

    fprintf('\n');
end