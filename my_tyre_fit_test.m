function fit_my_tire_to_weber_brush_ONEFILE()
% Fit your tire model parameters (Bx,Cx,d1x,d2x,By,Cy,d1y,d2y)
% so that |F| matches Weber/Supra combined-slip brush model |F|
% over:
%   kappa in [-1.5, 1.5]
%   alpha in [-45, 45] deg
%
% References (from your uploaded PDFs):
% - Weber Thesis Table 3.2 (Takumi Supra parameters) :contentReference[oaicite:2]{index=2}
% - Combined-slip brush eqs (16)-(19) in drifting paper :contentReference[oaicite:3]{index=3}
%
% Notes:
% - Brush uses kappa/(kappa+1) and tan(alpha)/(kappa+1), singular at kappa=-1.
%   We handle kappa near -1 robustly by treating it as fully sliding (|F|=mu*Fz).
% - Your tire model here is copied from the ACTIVE part of tire_force_tanEllipse
%   subfunction in your my_model.m (same constants and formulas).

clc; close all;

%% Load your baseline parameters (if available)
if exist('my_params','file') == 2
    veh0 = my_params();
else
    warning('my_params.m not found. Using a minimal default veh struct.');
    veh0 = default_veh_struct();
end

%% Supra (Takumi) rear brush parameters (Weber Table 3.2)
% Rear axle stiffnesses (axle-lumped):
CxR_axle = 213e3;   % N/rad
CyR_axle = 185e3;   % N/rad
% Rear affine friction parameters:
muR0 = 0.814;
muR1 = 8.2e-6;      % 1/N

% Convert axle-lumped stiffness to per-wheel (consistent with your per-wheel tire model)
Cx_brush = CxR_axle/2;
Cy_brush = CyR_axle/2;

% Choose DeltaFz,long for mu_r = muR0 + muR1*DeltaFz,long
DeltaFz_long = 0;            % default
mu_brush = muR0 + muR1*DeltaFz_long;

fprintf('Brush target (per-wheel): Cx=%.3g N/rad, Cy=%.3g N/rad, mu=%.6g\n', Cx_brush, Cy_brush, mu_brush);

%% Slip grids requested by you
kappa_vec = linspace(-1.5, 1.5, 61);
alpha_deg_vec = linspace(-45, 45, 61);
[KAP, ALPdeg] = ndgrid(kappa_vec, alpha_deg_vec);
TANALP = tand(ALPdeg);

%% Choose representative Fz levels (per wheel) for fitting
FzR_wheel_static = rear_static_Fz_per_wheel(veh0);
Fz_levels = FzR_wheel_static * [0.7, 1.0, 1.3];  % you can edit
fprintf('Using Fz levels per wheel [N]: %s\n', mat2str(Fz_levels,4));

%% Precompute brush target magnitude surfaces
Fmag_brush = cell(numel(Fz_levels),1);
for i = 1:numel(Fz_levels)
    Fz = Fz_levels(i);
    Fmag_brush{i} = brush_combined_mag(KAP, TANALP, Fz, Cx_brush, Cy_brush, mu_brush);
end

%% Parameter vector to tune
% p = [Bx Cx d1x d2x By Cy d1y d2y]
p0 = [veh0.Bx, veh0.Cx, veh0.d1x, veh0.d2x, veh0.By, veh0.Cy, veh0.d1y, veh0.d2y];

% Bounds (IMPORTANT: Cx,Cy must be > 1 to avoid tan(pi/(2*C)) singularity)
lb = [  0.5, 1.05, 0.00, -5000,  0.5, 1.05, 0.00, -5000 ];
ub = [ 80.0, 2.80, 3.00,  5000, 80.0, 2.80, 3.00,  5000 ];
p0 = max(lb, min(ub, p0));

%% Solve (lsqnonlin preferred, otherwise fminsearch fallback)
objfun = @(p) residual_vector_mag_only( ...
    p, veh0, KAP, TANALP, Fz_levels, Fmag_brush, mu_brush, Cx_brush, Cy_brush);

use_lsq = (exist('lsqnonlin','file') == 2);
if use_lsq
    opts = optimoptions('lsqnonlin',...
        'Display','iter',...
        'MaxFunctionEvaluations', 8000,...
        'MaxIterations', 300,...
        'FunctionTolerance', 1e-10,...
        'StepTolerance', 1e-12);
    p_fit = lsqnonlin(objfun, p0, lb, ub, opts);
else
    warning('lsqnonlin not found (Optimization Toolbox). Using fminsearch fallback (slower).');
    costfun = @(p) sum(objfun(project_to_bounds(p,lb,ub)).^2);
    p_fit = fminsearch(costfun, p0);
    p_fit = project_to_bounds(p_fit, lb, ub);
end

veh_fit = veh0;
veh_fit.Bx  = p_fit(1); veh_fit.Cx  = p_fit(2); veh_fit.d1x = p_fit(3); veh_fit.d2x = p_fit(4);
veh_fit.By  = p_fit(5); veh_fit.Cy  = p_fit(6); veh_fit.d1y = p_fit(7); veh_fit.d2y = p_fit(8);

fprintf('\n=== FIT DONE ===\n');
fprintf('Fitted tire params:\n');
fprintf('  Bx=%.8g, Cx=%.8g, d1x=%.8g, d2x=%.8g\n', veh_fit.Bx, veh_fit.Cx, veh_fit.d1x, veh_fit.d2x);
fprintf('  By=%.8g, Cy=%.8g, d1y=%.8g, d2y=%.8g\n', veh_fit.By, veh_fit.Cy, veh_fit.d1y, veh_fit.d2y);

%% Plot results
plot_fit(veh0, veh_fit, KAP, ALPdeg, TANALP, Fz_levels, Fmag_brush, mu_brush, Cx_brush, Cy_brush);

end

%% ========================= Residual (MAGNITUDE ONLY) =========================
function r = residual_vector_mag_only(p, veh0, KAP, TANALP, Fz_levels, Fmag_brush, mu_brush, Cx_brush, Cy_brush)
veh = veh0;
veh.Bx  = p(1); veh.Cx  = p(2); veh.d1x = p(3); veh.d2x = p(4);
veh.By  = p(5); veh.Cy  = p(6); veh.d1y = p(7); veh.d2y = p(8);

r_all = [];

for i = 1:numel(Fz_levels)
    Fz = Fz_levels(i);

    % Brush target magnitude
    Fmag_t = Fmag_brush{i};

    % Your tire model magnitude (copied from your my_model.m tire_force_tanEllipse)
    [Fx_u, Fy_u] = my_tire_force_tanEllipse_like_my_model(veh, KAP, TANALP, Fz);
    Fmag_u = hypot(Fx_u, Fy_u);

    denom = max(mu_brush*Fz, 1.0);
    res = (Fmag_u - Fmag_t) / denom;

    % Weight knee region heavier (most important for drifting transitions)
    s = Fmag_t / denom; % in [0,1]
    w = ones(size(s));
    w(s > 0.2 & s < 0.98) = 4.0;
    w(s >= 0.98) = 1.5;

    % Also downweight near kappa = -1 where brush mapping is singular (we saturate there anyway)
    den = 1 + KAP;
    w(abs(den) < 0.02) = 0.5;

    r_all = [r_all; sqrt(w(:)).*res(:)];
end

% Soft penalties to prevent negative mu_max over the interpolation range
% (since your model defines mu_max = (d1*Fz + d2)/Fz)
Fz1 = 2000; Fz2 = 6000;
Fz_test = linspace(Fz1, Fz2, 7);
mu_x = (veh.d1x*Fz_test + veh.d2x)./Fz_test;
mu_y = (veh.d1y*Fz_test + veh.d2y)./Fz_test;
pen = [min(mu_x)-0.05, min(mu_y)-0.05]; % want >=0.05
pen = max(0, -pen);                     % positive if violated
r_all = [r_all; 5.0*pen(:)];

% Optional: prevent extreme k_max / tan_alpha_max (numerical stability)
[~, ~, info] = internal_slip_scales(veh);
r_all = [r_all; 0.05*max(0, info.k_max_min - info.k_max);
               0.05*max(0, info.tan_alpha_max_min - info.tan_alpha_max)];

r = r_all;
end

%% ========================= Brush model magnitude =========================
function Fmag = brush_combined_mag(kappa, tan_alpha, Fz, Cx, Cy, mu)
% From combined-slip brush model equations (16)-(19) :contentReference[oaicite:4]{index=4}
% Using:
%   sx = Cx * (kappa/(1+kappa))
%   sy = Cy * (tan(alpha)/(1+kappa))
%   f  = sqrt(sx^2 + sy^2)
%   F  = polynomial in f up to saturation at mu*Fz
% Magnitude is F.

den = 1 + kappa;

% Handle singularity at kappa=-1 robustly:
% if den <= 0 (or very small), treat as fully sliding.
den_safe = den;
den_safe(den_safe < 1e-3) = 1e-3;

sx = Cx .* (kappa ./ den_safe);
sy = Cy .* (tan_alpha ./ den_safe);
f  = sqrt(sx.^2 + sy.^2 + 1e-12);

Fmax = mu .* Fz;
r = f ./ (3*Fmax);

g = 3*r - 3*r.^2 + r.^3;
g(r > 1) = 1;
Fmag = Fmax .* g;

% For kappa <= -1, fully sliding anyway:
Fmag(den <= 0) = Fmax;
end

%% ========================= Your tire model (copied from my_model.m) =========================
function [Fx, Fy] = my_tire_force_tanEllipse_like_my_model(veh, kappa, tan_alpha, Fz)
% This matches the ACTIVE part of tire_force_tanEllipse in your my_model.m

Fz1 = 2000;   % N
Fz2 = 6000;   % N
Qx  = 1.35;
Qy  = 1.5;
eps1 = 1e-5;
eps2 = 1e-5;

% mu_x_max1/2, mu_y_max1/2 as in your file
mu_x_max1 = (veh.d1x*Fz1 + veh.d2x)/Fz1;  mu_x_max2 = (veh.d1x*Fz2 + veh.d2x)/Fz2;
mu_y_max1 = (veh.d1y*Fz1 + veh.d2y)/Fz1;  mu_y_max2 = (veh.d1y*Fz2 + veh.d2y)/Fz2;

% Guard Cx,Cy > 1 to avoid tan(pi/(2*C)) singular
Cx = max(veh.Cx, 1.01);
Cy = max(veh.Cy, 1.01);

k_max = tan(pi/(2*Cx))/veh.Bx;
alpha_peak = tan(pi/(2*Cy))/veh.By;
tan_alpha_max = tan(alpha_peak);

% Clamp Fz within [Fz1,Fz2] for stability (your file extrapolates; clamp is safer for fitting)
Fz_clamped = min(max(Fz, Fz1), Fz2);

mu_x_max = (Fz_clamped - Fz1) * (mu_x_max2 - mu_x_max1) / (Fz2 - Fz1) + mu_x_max1;
mu_y_max = (Fz_clamped - Fz1) * (mu_y_max2 - mu_y_max1) / (Fz2 - Fz1) + mu_y_max1;

k_n = kappa / k_max;
a_n = (-tan_alpha) / tan_alpha_max;

rou = (k_n.^2 + a_n.^2 + eps1).^0.5 + eps2;

Sx = pi/(2*atan(Qx));
Sy = pi/(2*atan(Qy));

my_scaling = 1.0;

Fx = (my_scaling * mu_x_max) .* Fz .* sin(Qx * atan(Sx * rou)) .* (k_n ./ rou);
Fy = (my_scaling * mu_y_max) .* Fz .* sin(Qy * atan(Sy * rou)) .* (a_n ./ rou);
end

%% ========================= Plots =========================
function plot_fit(veh0, veh_fit, KAP, ALPdeg, TANALP, Fz_levels, Fmag_brush, mu_brush, Cx_brush, Cy_brush)

levels = [0.2 0.4 0.6 0.8 0.9 0.95 1.0];

for i = 1:numel(Fz_levels)
    Fz = Fz_levels(i);
    denom = max(mu_brush*Fz, 1.0);

    Fmag_t = Fmag_brush{i};

    [Fx0, Fy0] = my_tire_force_tanEllipse_like_my_model(veh0,   KAP, TANALP, Fz);
    [Fx1, Fy1] = my_tire_force_tanEllipse_like_my_model(veh_fit, KAP, TANALP, Fz);

    Fmag0 = hypot(Fx0, Fy0);
    Fmag1 = hypot(Fx1, Fy1);

    figure('Name', sprintf('Contours Fz=%.0fN',Fz));
    hold on; grid on;
    contour(KAP, ALPdeg, Fmag_t/denom, levels, 'LineWidth', 1.6);
    contour(KAP, ALPdeg, Fmag0/denom, levels, '--', 'LineWidth', 1.2);
    contour(KAP, ALPdeg, Fmag1/denom, levels, '-.', 'LineWidth', 1.2);
    xlabel('\kappa'); ylabel('\alpha [deg]');
    title(sprintf('Normalized |F| contours (Brush solid, Your0 dashed, YourFit dashdot), Fz=%.0fN',Fz));
    legend('Brush','Your0','YourFit','Location','best');
    xlim([min(KAP(:)) max(KAP(:))]); ylim([min(ALPdeg(:)) max(ALPdeg(:))]);

    % 1D cuts: alpha=0 and kappa=0
    figure('Name', sprintf('Cuts Fz=%.0fN',Fz));

    subplot(1,2,1);
    [~, idx_a0] = min(abs(ALPdeg(1,:) - 0));
    plot(KAP(:,idx_a0), Fmag_t(:,idx_a0)/denom, 'LineWidth',1.6); hold on; grid on;
    plot(KAP(:,idx_a0), Fmag0(:,idx_a0)/denom, '--', 'LineWidth',1.2);
    plot(KAP(:,idx_a0), Fmag1(:,idx_a0)/denom, '-.', 'LineWidth',1.2);
    xlabel('\kappa'); ylabel('|F|/(\mu F_z)'); title('\alpha = 0 deg');
    legend('Brush','Your0','YourFit','Location','best');

    subplot(1,2,2);
    [~, idx_k0] = min(abs(KAP(:,1) - 0));
    plot(ALPdeg(idx_k0,:), Fmag_t(idx_k0,:)/denom, 'LineWidth',1.6); hold on; grid on;
    plot(ALPdeg(idx_k0,:), Fmag0(idx_k0,:)/denom, '--', 'LineWidth',1.2);
    plot(ALPdeg(idx_k0,:), Fmag1(idx_k0,:)/denom, '-.', 'LineWidth',1.2);
    xlabel('\alpha [deg]'); ylabel('|F|/(\mu F_z)'); title('\kappa = 0');
    legend('Brush','Your0','YourFit','Location','best');
end

end

%% ========================= Small helpers =========================
function p = project_to_bounds(p, lb, ub)
p = max(lb, min(ub, p));
end

function veh = default_veh_struct()
veh = struct();
veh.g  = 9.80665;
veh.m  = 1496;
veh.lf = 1.217;
veh.lr = 1.231;
veh.l  = veh.lf + veh.lr;

% Your initial tire params from my_params.m (defaults)
veh.Bx  = 18;   veh.Cx  = 1.3;   veh.d1x = 0.95; veh.d2x = 320;
veh.By  = 13;   veh.Cy  = 1.5;   veh.d1y = 0.95; veh.d2y = 320;
end

function FzR_wheel = rear_static_Fz_per_wheel(veh)
% In your my_model.m static rear axle load is m*g*lf/l (then /2 per wheel)
m = veh.m; g = veh.g;

if isfield(veh,'lf') && isfield(veh,'l')
    FzR_axle = m*g * veh.lf/veh.l;
elseif isfield(veh,'a') && isfield(veh,'L')
    FzR_axle = m*g * veh.a/veh.L;
else
    warning('Cannot infer lf/l. Using 0.5*m*g for rear axle load.');
    FzR_axle = 0.5*m*g;
end

FzR_wheel = 0.5 * FzR_axle;
end

function [k_max, tan_alpha_max, info] = internal_slip_scales(veh)
Cx = max(veh.Cx, 1.01);
Cy = max(veh.Cy, 1.01);
k_max = tan(pi/(2*Cx))/veh.Bx;
alpha_peak = tan(pi/(2*Cy))/veh.By;
tan_alpha_max = tan(alpha_peak);
info.k_max = k_max;
info.tan_alpha_max = tan_alpha_max;
info.k_max_min = 0.02;
info.tan_alpha_max_min = 0.02;
end