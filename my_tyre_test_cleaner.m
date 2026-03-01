%% Tyre model comparison (clean version)
% clear; close all; clc;

%% ======================
%% User settings
%% ======================
Fz = 3.0625e+03;   % N

kappa_vec     = linspace(-1.0, 1.0, 161*6);
tan_alpha_vec = linspace(-tan(deg2rad(70)), tan(deg2rad(70)), 161*6);
[KAPPA, TANALPHA] = meshgrid(kappa_vec, tan_alpha_vec);

% Plot switches
DO_PLOT_OVERLAY        = true;  % overlay 3 surfaces: model1 + monotone + old
DO_PLOT_OLD_VIOLATION  = true;  % old model surface with mufl_ineq<=0 indication

% Monotone cap smoothness
kcap = 30;

%% ======================
%% Tyre/vehicle parameters
%% ======================
veh = struct();

veh.Bx  = 18;
veh.Cx  = 1.3;
veh.d1x = 0.95;
veh.d2x = 320;

veh.By  = 13;
veh.Cy  = 1.5;
veh.d1y = 0.95;
veh.d2y = 320;

veh.mux = 1.0;
veh.muy = 1.0;

%% ======================
%% Model parameters (put all "magic numbers" in one place)
%% ======================
par = struct();
par.Fz1  = 2000;
par.Fz2  = 6000;
par.Qx   = 1.35;
par.Qy   = 1.5;
par.eps1 = 1e-5;
par.eps2 = 1e-5;

% Old model eps to avoid divide by zero
par.eps_s = 1e-12;

%% ======================
%% Evaluate all models (ONLY ONCE each)
%% ======================
[Fx1, Fy1, meta1] = eval_model1(veh, KAPPA, TANALPHA, Fz, par);
Ftot1 = hypot(Fx1, Fy1);

[Fx1m, Fy1m, meta1m] = eval_model1_monotone(veh, KAPPA, TANALPHA, Fz, par, kcap);
Ftot1m = hypot(Fx1m, Fy1m);

old = eval_old_model(veh, KAPPA, TANALPHA, Fz, par);
Ftot_old = old.Ftot;           % hypot(old.fx, old.fy)
mufl_ineq = old.mufl_ineq;     % 1 - (fx^2/(mux^2 Fz^2) + fy^2/(muy^2 Fz^2))

%% ======================
%% Plot: overlay 3 surfaces
%% ======================
if DO_PLOT_OVERLAY
    figure;
    h1 = surf(KAPPA, TANALPHA, Ftot1,    "EdgeColor","none"); hold on;
    h2 = surf(KAPPA, TANALPHA, Ftot1m,   "EdgeColor","none");
    h3 = surf(KAPPA, TANALPHA, Ftot_old, "EdgeColor","none");

    h1.FaceColor = [0.10 0.45 0.80]; h1.FaceAlpha = 0.35; h1.DisplayName = sprintf("Model 1 (scale=%.3f)", meta1.my_scaling);
    h2.FaceColor = [0.20 0.70 0.25]; h2.FaceAlpha = 0.35; h2.DisplayName = sprintf("Monotone (scale=%.3f)", meta1m.my_scaling);
    h3.FaceColor = [0.85 0.33 0.10]; h3.FaceAlpha = 0.35; h3.DisplayName = "Old model";

    xlabel('\kappa (slip ratio)');
    ylabel('tan(\alpha)');
    zlabel('F_{tot} [N]');
    title('Overlay: Model 1 vs Monotone Model 1 vs Old model');
    grid on; view(45, 30);
    legend("show", "Location", "best");
    colorbar;

    % Optional: show old-model boundary mufl_ineq = 0 on same axes
    % contour3(KAPPA, TANALPHA, mufl_ineq, [0 0], 'k', 'LineWidth', 2, 'DisplayName','old: mufl\_ineq=0');
end

%% ======================
%% Plot: old model + mufl_ineq<=0 region
%% ======================
mask_feas = (mufl_ineq >= 0);

figure; hold on;

% Base surface (whole old model), muted
hBase = surf(KAPPA, TANALPHA, Ftot_old, "EdgeColor","none");
hBase.FaceColor = [0.75 0.75 0.75];
hBase.FaceAlpha = 0.20;

% Feasible region overlay (GREEN)
Zfeas = Ftot_old;
Zfeas(~mask_feas) = NaN;                 % show only feasible region
hFeas = surf(KAPPA, TANALPHA, Zfeas, "EdgeColor","none");
hFeas.FaceColor = [0.10 0.75 0.10];
hFeas.FaceAlpha = 0.90;

% (Optional) boundary mufl_ineq = 0 drawn ON the Ftot surface
C = contourc(kappa_vec, tan_alpha_vec, mufl_ineq, [0 0]);
idx = 1;
while idx < size(C,2)
    n  = C(2,idx);
    xC = C(1, idx+1:idx+n);
    yC = C(2, idx+1:idx+n);
    zC = interp2(KAPPA, TANALPHA, Ftot_old, xC, yC, "linear");
    plot3(xC, yC, zC, 'k-', 'LineWidth', 2.5);
    idx = idx + n + 1;
end

xlabel('\kappa (slip ratio)');
ylabel('tan(\alpha)');
zlabel('F_{tot} [N]');
title('Old model: feasible region (mufl\_ineq \ge 0) highlighted');
grid on; view(45,30);
legend([hBase hFeas], {'F_{tot} (old model)', 'Feasible: mufl\_ineq \ge 0'}, 'Location','best');
camlight headlight; lighting gouraud;


%% =======================================================================
%% Local functions
%% =======================================================================

function [Fx, Fy, meta] = eval_model1(veh, kappa, tan_alpha, Fz, par)
    [mu_x_max, mu_y_max, k_max, tan_alpha_max] = common_load_dependent_mu(veh, Fz, par);

    k_n = kappa ./ k_max;
    a_n = (-tan_alpha) ./ tan_alpha_max;
    rou = sqrt(k_n.^2 + a_n.^2 + par.eps1) + par.eps2;

    Sx = pi/(2*atan(par.Qx));
    Sy = pi/(2*atan(par.Qy));

    my_scaling = auto_scaling_match_old_model(veh, Fz, mu_x_max, mu_y_max, par.Qx, par.Qy);

    gx = sin(par.Qx * atan(Sx * rou));
    gy = sin(par.Qy * atan(Sy * rou));

    Fx = (my_scaling * mu_x_max) .* Fz .* gx .* (k_n ./ rou);
    Fy = (my_scaling * mu_y_max) .* Fz .* gy .* (a_n ./ rou);

    meta = struct("my_scaling", my_scaling, "mu_x_max", mu_x_max, "mu_y_max", mu_y_max);
end

function [Fx, Fy, meta] = eval_model1_monotone(veh, kappa, tan_alpha, Fz, par, kcap)
    [mu_x_max, mu_y_max, k_max, tan_alpha_max] = common_load_dependent_mu(veh, Fz, par);

    k_n = kappa ./ k_max;
    a_n = (-tan_alpha) ./ tan_alpha_max;
    rou = sqrt(k_n.^2 + a_n.^2 + par.eps1) + par.eps2;

    Sx = pi/(2*atan(par.Qx));
    Sy = pi/(2*atan(par.Qy));

    % Same scaling as Model 1 (recommended for consistent magnitude)
    my_scaling = auto_scaling_match_old_model(veh, Fz, mu_x_max, mu_y_max, par.Qx, par.Qy);

    theta_x = par.Qx * atan(Sx * rou);
    theta_y = par.Qy * atan(Sy * rou);

    theta_x_cap = smoothmin(theta_x, pi/2, kcap);
    theta_y_cap = smoothmin(theta_y, pi/2, kcap);

    gx = sin(theta_x_cap);
    gy = sin(theta_y_cap);

    Fx = (my_scaling * mu_x_max) .* Fz .* gx .* (k_n ./ rou);
    Fy = (my_scaling * mu_y_max) .* Fz .* gy .* (a_n ./ rou);

    meta = struct("my_scaling", my_scaling, "mu_x_max", mu_x_max, "mu_y_max", mu_y_max);
end

function out = eval_old_model(veh, kappa, tan_alpha, Fz, par)
    s = sqrt(kappa.^2 + tan_alpha.^2 + par.eps_s);

    Fx_mag = (veh.d1x*Fz + veh.d2x) .* sin(veh.Cx * atan(veh.Bx * s));
    Fy_mag = -(veh.d1y*Fz + veh.d2y) .* sin(veh.Cy * atan(veh.By * s));

    fx = (kappa    ./ s) .* Fx_mag;
    fy = (tan_alpha ./ s) .* Fy_mag;

    Ftot = hypot(fx, fy);

    mufl_ineq = 1 - ( (fx.^2)./(veh.mux^2 * Fz^2) + (fy.^2)./(veh.muy^2 * Fz^2) );

    out = struct("fx", fx, "fy", fy, "Ftot", Ftot, "mufl_ineq", mufl_ineq);
end

function [mu_x_max, mu_y_max, k_max, tan_alpha_max] = common_load_dependent_mu(veh, Fz, par)
    Fz1 = par.Fz1; Fz2 = par.Fz2;

    mu_x_max1 = (veh.d1x*Fz1 + veh.d2x)/Fz1;  mu_x_max2 = (veh.d1x*Fz2 + veh.d2x)/Fz2;
    mu_y_max1 = (veh.d1y*Fz1 + veh.d2y)/Fz1;  mu_y_max2 = (veh.d1y*Fz2 + veh.d2y)/Fz2;

    k_max = tan(pi/(2*veh.Cx))/veh.Bx;
    alpha_peak = tan(pi/(2*veh.Cy))/veh.By;
    tan_alpha_max = tan(alpha_peak);

    mu_x_max = (Fz - Fz1) .* (mu_x_max2 - mu_x_max1) ./ (Fz2 - Fz1) + mu_x_max1;
    mu_y_max = (Fz - Fz1) .* (mu_y_max2 - mu_y_max1) ./ (Fz2 - Fz1) + mu_y_max1;
end

function m = smoothmin(a, b, k)
    m0 = min(a, b);
    m  = m0 - (1./k).*log(exp(-k.*(a - m0)) + exp(-k.*(b - m0)));
end

function my_scaling = auto_scaling_match_old_model(veh, Fz, mu_x_max, mu_y_max, Qx, Qy)
    old_peak_fac_x = peak_sin_factor(veh.Cx);
    old_peak_fac_y = peak_sin_factor(veh.Cy);

    new_peak_fac_x = peak_sin_factor(Qx);
    new_peak_fac_y = peak_sin_factor(Qy);

    old_peak_Fx = abs((veh.d1x * Fz + veh.d2x) * old_peak_fac_x);
    old_peak_Fy = abs((veh.d1y * Fz + veh.d2y) * old_peak_fac_y);

    new_peak_Fx = abs((mu_x_max * Fz) * new_peak_fac_x);
    new_peak_Fy = abs((mu_y_max * Fz) * new_peak_fac_y);

    sx = old_peak_Fx / max(new_peak_Fx, eps);
    sy = old_peak_Fy / max(new_peak_Fy, eps);

    my_scaling = sqrt(sx * sy);
    my_scaling = min(max(my_scaling, 0.05), 20);
end

function pf = peak_sin_factor(C)
    if C >= 1
        pf = 1;
    else
        pf = sin(C*pi/2);
    end
end
