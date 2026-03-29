% import casadi.*
% 
% veh = my_params();
% 
% %% Vehicle Model - state variables
% 
% nx = 15; % number of state variables
% 
% % longitudinal velocity [m/s]
% V_n = SX.sym('v_n');
% V_s = 100;
% V = V_s * V_n;
% 
% % sideslip angle [rad]
% beta_n = SX.sym('beta_n');
% beta_s = 1;
% beta = beta_s * beta_n;
% 
% % yaw rate [rad/s]
% gamma_n = SX.sym('gamma_n');
% gamma_s = 1;
% gamma = gamma_s * gamma_n;
% 
% % lateral distance to centreline [m] - left of centreline => n > 0; right => n < 0
% n_n = SX.sym('n_n');
% n_s = 5;
% n = n_s * n_n;
% 
% % course angle to centreline tangent direction [rad]
% xi_n = SX.sym('xi_n');
% xi_s = 1;
% xi = xi_s * xi_n;
% 
% ax_bar_n = SX.sym('ax_bar_n');
% ax_bar_s = 2*9.8;
% ax_bar = ax_bar_s * ax_bar_n;
% 
% ay_bar_n = SX.sym('ay_bar_n');
% ay_bar_s = 2*9.8;
% ay_bar = ay_bar_s * ay_bar_n;
% 
% omega_s = V_s/veh.rw;
% 
% % angular velocity front left tyre [rad/s]
% omega_fl_n = SX.sym('omega_fl_n');
% omega_fl = omega_s * omega_fl_n;
% 
% % angular velocity front right tyre [rad/s]
% omega_fr_n = SX.sym('omega_fr_n');
% omega_fr = omega_s * omega_fr_n;
% 
% % angular velocity rear left tyre [rad/s]
% omega_rl_n = SX.sym('omega_rl_n');
% omega_rl = omega_s * omega_rl_n;
% 
% % angular velocity rear right tyre [rad/s]
% omega_rr_n = SX.sym('omega_rr_n');
% omega_rr = omega_s * omega_rr_n;
% 
% % NEW: chassis roll / pitch states
% phi_n = SX.sym('phi_n');
% phi_s = 0.2;              % [rad] scaling
% 
% phi_dot_n = SX.sym('phi_dot_n');
% phi_dot_s = 1.0;          % [rad/s] scaling
% 
% theta_n = SX.sym('theta_n');
% theta_s = 0.2;            % [rad] scaling
% 
% theta_dot_n = SX.sym('theta_dot_n');
% theta_dot_s = 1.0;        % [rad/s] scaling
% 
% phi = phi_s * phi_n;
% phi_dot = phi_dot_s * phi_dot_n;
% theta = theta_s * theta_n;
% theta_dot = theta_dot_s * theta_dot_n;
% 
% % states vector (scaled)
% x = [V_n beta_n gamma_n n_n xi_n ax_bar_n ay_bar_n ...
%      omega_fl_n omega_fr_n omega_rl_n omega_rr_n ...
%      phi_n phi_dot_n theta_n theta_dot_n]';
% 
% % scaling factors
% x_s = [V_s beta_s gamma_s n_s xi_s ax_bar_s ay_bar_s ...
%        omega_s omega_s omega_s omega_s ...
%        phi_s phi_dot_s theta_s theta_dot_s]';
% 
% % state limits
% omega_max = veh.V_max/veh.rw;
% x_min = [   -1e-3; -pi/4; -pi/2; -4; -pi/4; -2*9.8; -2*9.8; ...
%                  0;     0;     0;     0; ...
%             -pi/4;   -10; -pi/4;   -10] ./ x_s;
% 
% x_max = [veh.V_max;  pi/4;  pi/2;  4;  pi/4;  2*9.8;  2*9.8; ...
%           omega_max; omega_max; omega_max; omega_max; ...
%              pi/4;    10;    pi/4;    10] ./ x_s;
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Keep these for monitoring / exact steady-state matching
% Fx_tires = ax_bar * veh.m + f_drag;
% Fy_tires = ay_bar * veh.m;
% 
% Fdown = f_lift;
% Fdrag = f_drag;
% 
% acom = veh.lf;
% bcom = veh.lr;
% mass = veh.m;
% gravity = veh.g;
% hcom = veh.hc;
% twF_h = veh.twF_h; % front-wheel to car-centre-line distance
% twR_h = veh.twR_h;
% dr = veh.dr;
% 
% acop = acom;
% hcop = hcom;
% 
% % NEW: suspension-transmitted chassis moments
% % Sign convention used here:
% %   Mphi_sus   > 0  => increases right-side loads, decreases left-side loads
% %   Mtheta_sus > 0  => increases front axle load, decreases rear axle load
% Mphi_sus   = veh.Cphi   * phi_dot   + veh.Kphi   * phi;
% Mtheta_sus = veh.Ctheta * theta_dot + veh.Ktheta * theta;
% 
% [fzfl, fzfr, fzrl, fzrr] = wheel_load( ...
%     Mtheta_sus, Mphi_sus, Fdown, Fdrag, ...
%     bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = (fxfl+fxfr)*cosd - (fyfl+fyfr)*sind + (fxrl+fxrr) - f_drag;
% Y = (fxfl+fxfr)*sind + (fyfl+fyfr)*cosd + (fyrl+fyrr);
% 
% ax = X / veh.m;
% ay = Y / veh.m;
% 
% Mz = (-1)*fxfl*cosd* (veh.wt/2) + fxfl*sind*veh.lf + fyfl*sind*(veh.wt/2) + fyfl*cosd*veh.lf + ...
%           fxfr*cosd* (veh.wt/2) + fxfr*sind*veh.lf - fyfr*sind*(veh.wt/2) + fyfr*cosd*veh.lf + ...
%      (-1)*fxrl*(veh.wt/2) - fyrl*veh.lr + ...
%           fxrr*(veh.wt/2) - fyrr*veh.lr;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % NEW: chassis roll / pitch dynamics
% % Use ax_bar, ay_bar as excitations so that, at steady state, the new
% % wheel-load equations collapse exactly to your current algebraic
% % wheel_load(...) equations.
% phi_ddot   = (-mass*hcom*ay_bar - Mphi_sus)   / veh.Ix;
% theta_ddot = (-mass*hcom*ax_bar - Mtheta_sus) / veh.Iy;
% 
% % Change of independent variable
% chi = xi + beta;
% sf = (1 - n*kappa) / (V*cos(chi) + 1e-12) + 1e-12;
% 
% xdot = [
%     (X*cosb + Y*sinb) / veh.m;
%     (Y*cosb - X*sinb) / (V*veh.m) - gamma;
%     Mz / veh.Iz;
%     V*sin(chi);
%     gamma - kappa/sf;
% 
%     (ax - ax_bar) / veh.load_transfer_time_delay;
%     (ay - ay_bar) / veh.load_transfer_time_delay;
% 
%     (Tfl - fxfl*veh.rw) / veh.Iw;
%     (Tfr - fxfr*veh.rw) / veh.Iw;
%     (Trl - fxrl*veh.rw) / veh.Iw;
%     (Trr - fxrr*veh.rw) / veh.Iw;
% 
%     phi_dot;
%     phi_ddot;
%     theta_dot;
%     theta_ddot;
%     ] ./ x_s;
% 
% dx = sf * xdot;
% 
% 
% 
% 
% 
% 
% 
% 
% function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load( ...
%     Mtheta_sus, Mphi_sus, Fdown, Fdrag, ...
%     bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h)
% 
% L = acom + bcom;
% Nz_total = mass*gravity + Fdown;
% 
% % Positive Mtheta_sus loads the FRONT axle, unloads the REAR axle
% % This external pitch term is required for exact steady-state recovery of
% % your current algebraic wheel_load(...) equations.
% Mtheta_ext = -(Fdrag*hcop + Fdown*(acop-acom));
% 
% front_axle_total = ( bcom*Nz_total + Mtheta_sus + Mtheta_ext ) / L;
% rear_axle_total  = ( acom*Nz_total - Mtheta_sus - Mtheta_ext ) / L;
% 
% % Positive Mphi_sus loads RIGHT side, unloads LEFT side
% Dphi = dr*twF_h + (1-dr)*twR_h;
% 
% dF_front = dr      * Mphi_sus / Dphi;
% dF_rear  = (1-dr)  * Mphi_sus / Dphi;
% 
% Fz_FL = 0.5 * (front_axle_total - dF_front);
% Fz_FR = 0.5 * (front_axle_total + dF_front);
% 
% Fz_RL = 0.5 * (rear_axle_total  - dF_rear);
% Fz_RR = 0.5 * (rear_axle_total  + dF_rear);
% 
% end















% 
% %% compare_wheel_load_models.m
% % Compare:
% %   1) original wheel-load equations
% %   2) new dynamic wheel-load equations
% %   3) new quasi-static wheel-load equations
% %
% % Assumptions:
% %   - Uses your current my_params()
% %   - Uses the sign conventions agreed earlier:
% %       Mtheta_sus > 0  => increases FRONT axle load, decreases REAR axle load
% %       Mphi_sus   > 0  => increases RIGHT-side loads, decreases LEFT-side loads
% %
% % Notes:
% %   - The "new dynamic" wheel loads are obtained by integrating roll/pitch ODEs
% %     driven by ax_bar and ay_bar.
% %   - The "new quasi-static" wheel loads use:
% %         Mphi_sus_qs   = -m*h*ay_bar
% %         Mtheta_sus_qs = -m*h*ax_bar
% %     which should collapse exactly to your original wheel-load equations.
% %
% % Replace the demo input signals (V, ax_bar, ay_bar) below with your own data
% % if desired.
% 
% clear
% close all
% clc
% 
% %% ------------------------------------------------------------------------
% % Parameters
% % -------------------------------------------------------------------------
% veh = my_params();
% 
% % Consistency check
% assert(isfield(veh,'roll_pitch_ode_initial_condition'), ...
%     'veh.roll_pitch_ode_initial_condition is missing.');
% assert(numel(veh.roll_pitch_ode_initial_condition) == 4, ...
%     'veh.roll_pitch_ode_initial_condition must be a 4x1 vector [phi; phi_dot; theta; theta_dot].');
% 
% %% ------------------------------------------------------------------------
% % USER INPUTS / DEMO INPUTS
% % Replace these with your own time histories if you already have them
% % -------------------------------------------------------------------------
% dt = 0.01;
% t_end = 12.0;
% t = (0:dt:t_end).';
% 
% % Demo vehicle speed used only for aero drag / lift
% % If you have your own solved V array, replace this.
% V = 28 + 5*sin(2*pi*0.12*t) + 2*sin(2*pi*0.35*t + 0.4);  % [m/s]
% V = max(V, veh.V_min);
% 
% % Demo filtered accelerations used by your original wheel-load equations
% % Replace with your own ax_bar, ay_bar if available
% ax_bar = 0.8*sin(2*pi*0.18*t) - 0.35*sin(2*pi*0.55*t + 0.3);   % [m/s^2]
% ay_bar = 4.5*sin(2*pi*0.16*t + 0.5) + 1.0*sin(2*pi*0.48*t);    % [m/s^2]
% 
% % Optional: if you want to use your own arrays from workspace, uncomment and adapt:
% % load('your_data.mat','t','V','ax_bar','ay_bar');
% 
% assert(isvector(t) && isvector(V) && isvector(ax_bar) && isvector(ay_bar), ...
%     't, V, ax_bar, ay_bar must be vectors.');
% assert(numel(t)==numel(V) && numel(t)==numel(ax_bar) && numel(t)==numel(ay_bar), ...
%     't, V, ax_bar, ay_bar must have the same length.');
% 
% t      = t(:);
% V      = V(:);
% ax_bar = ax_bar(:);
% ay_bar = ay_bar(:);
% 
% %% ------------------------------------------------------------------------
% % Aero quantities
% % -------------------------------------------------------------------------
% % In your full model, f_drag uses vx = V*cos(beta).
% % Here, for this standalone comparison, we take beta = 0 so vx = V.
% vx_for_aero = V;
% 
% f_drag = 0.5 * veh.drag_coeff * vx_for_aero.^2;   % [N]
% f_lift = 0.5 * veh.lift_coeff * vx_for_aero.^2;   % [N]
% 
% %% ------------------------------------------------------------------------
% % ORIGINAL wheel-load equations
% % -------------------------------------------------------------------------
% Fx_tires = veh.m * ax_bar + f_drag;   % exactly as in your current model
% Fy_tires = veh.m * ay_bar;            % exactly as in your current model
% 
% Fz_FL_orig = zeros(size(t));
% Fz_FR_orig = zeros(size(t));
% Fz_RL_orig = zeros(size(t));
% Fz_RR_orig = zeros(size(t));
% 
% for k = 1:numel(t)
%     [Fz_FL_orig(k), Fz_FR_orig(k), Fz_RL_orig(k), Fz_RR_orig(k)] = ...
%         wheel_load_original( ...
%             Fx_tires(k), Fy_tires(k), f_lift(k), f_drag(k), ...
%             veh.lr, veh.m, veh.g, veh.hc, veh.hc, veh.lf, veh.lf, ...
%             veh.dr, veh.twF_h, veh.twR_h);
% end
% 
% %% ------------------------------------------------------------------------
% % NEW DYNAMIC wheel-load equations
% % Integrate roll/pitch ODEs driven by ax_bar and ay_bar
% % -------------------------------------------------------------------------
% z0 = veh.roll_pitch_ode_initial_condition(:);  % [phi; phi_dot; theta; theta_dot]
% 
% ode_fun = @(tt,zz) roll_pitch_ode(tt, zz, t, ax_bar, ay_bar, veh);
% 
% % Use ode45 for the comparison script
% opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
% [t_ode, z_ode] = ode45(ode_fun, [t(1), t(end)], z0, opts);
% 
% % Interpolate ODE result back to the original time grid
% phi       = interp1(t_ode, z_ode(:,1), t, 'pchip');
% phi_dot   = interp1(t_ode, z_ode(:,2), t, 'pchip');
% theta     = interp1(t_ode, z_ode(:,3), t, 'pchip');
% theta_dot = interp1(t_ode, z_ode(:,4), t, 'pchip');
% 
% % Suspension moments from dynamic roll/pitch states
% Mphi_sus_dyn   = veh.Cphi   * phi_dot   + veh.Kphi   * phi;
% Mtheta_sus_dyn = veh.Ctheta * theta_dot + veh.Ktheta * theta;
% 
% Fz_FL_new_dyn = zeros(size(t));
% Fz_FR_new_dyn = zeros(size(t));
% Fz_RL_new_dyn = zeros(size(t));
% Fz_RR_new_dyn = zeros(size(t));
% 
% for k = 1:numel(t)
%     [Fz_FL_new_dyn(k), Fz_FR_new_dyn(k), Fz_RL_new_dyn(k), Fz_RR_new_dyn(k)] = ...
%         wheel_load_new( ...
%             Mtheta_sus_dyn(k), Mphi_sus_dyn(k), f_lift(k), f_drag(k), ...
%             veh.lr, veh.m, veh.g, veh.hc, veh.hc, veh.lf, veh.lf, ...
%             veh.dr, veh.twF_h, veh.twR_h);
% end
% 
% %% ------------------------------------------------------------------------
% % NEW QUASI-STATIC wheel-load equations
% % These should match the original equations up to numerical roundoff
% % -------------------------------------------------------------------------
% Mphi_sus_qs   = -veh.m * veh.hc * ay_bar;
% Mtheta_sus_qs = -veh.m * veh.hc * ax_bar;
% 
% Fz_FL_new_qs = zeros(size(t));
% Fz_FR_new_qs = zeros(size(t));
% Fz_RL_new_qs = zeros(size(t));
% Fz_RR_new_qs = zeros(size(t));
% 
% for k = 1:numel(t)
%     [Fz_FL_new_qs(k), Fz_FR_new_qs(k), Fz_RL_new_qs(k), Fz_RR_new_qs(k)] = ...
%         wheel_load_new( ...
%             Mtheta_sus_qs(k), Mphi_sus_qs(k), f_lift(k), f_drag(k), ...
%             veh.lr, veh.m, veh.g, veh.hc, veh.hc, veh.lf, veh.lf, ...
%             veh.dr, veh.twF_h, veh.twR_h);
% end
% 
% %% ------------------------------------------------------------------------
% % Errors: original vs new quasi-static
% % -------------------------------------------------------------------------
% err_FL = Fz_FL_new_qs - Fz_FL_orig;
% err_FR = Fz_FR_new_qs - Fz_FR_orig;
% err_RL = Fz_RL_new_qs - Fz_RL_orig;
% err_RR = Fz_RR_new_qs - Fz_RR_orig;
% 
% fprintf('\n=============================================================\n');
% fprintf('Max abs difference: NEW quasi-static vs ORIGINAL\n');
% fprintf('Fz_FL: %.12e N\n', max(abs(err_FL)));
% fprintf('Fz_FR: %.12e N\n', max(abs(err_FR)));
% fprintf('Fz_RL: %.12e N\n', max(abs(err_RL)));
% fprintf('Fz_RR: %.12e N\n', max(abs(err_RR)));
% fprintf('=============================================================\n\n');
% 
% %% ------------------------------------------------------------------------
% % Sanity checks
% % -------------------------------------------------------------------------
% sum_orig = Fz_FL_orig + Fz_FR_orig + Fz_RL_orig + Fz_RR_orig;
% sum_dyn  = Fz_FL_new_dyn + Fz_FR_new_dyn + Fz_RL_new_dyn + Fz_RR_new_dyn;
% sum_qs   = Fz_FL_new_qs  + Fz_FR_new_qs  + Fz_RL_new_qs  + Fz_RR_new_qs;
% 
% figure('Name','Total vertical load check');
% plot(t, sum_orig, 'LineWidth', 1.5); hold on; grid on;
% plot(t, sum_dyn,  '--', 'LineWidth', 1.5);
% plot(t, sum_qs,   ':', 'LineWidth', 1.8);
% xlabel('t [s]');
% ylabel('Fz_{FL}+Fz_{FR}+Fz_{RL}+Fz_{RR} [N]');
% legend('Original','New dynamic','New quasi-static','Location','best');
% title('Total vertical load');
% 
% %% ------------------------------------------------------------------------
% % Plot roll / pitch states and moments
% % -------------------------------------------------------------------------
% figure('Name','Roll and Pitch States');
% tiledlayout(2,2);
% 
% nexttile;
% plot(t, phi, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('\phi [rad]');
% title('Roll angle');
% 
% nexttile;
% plot(t, phi_dot, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('\phi_{dot} [rad/s]');
% title('Roll rate');
% 
% nexttile;
% plot(t, theta, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('\theta [rad]');
% title('Pitch angle');
% 
% nexttile;
% plot(t, theta_dot, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('\theta_{dot} [rad/s]');
% title('Pitch rate');
% 
% figure('Name','Suspension Moments');
% tiledlayout(2,1);
% 
% nexttile;
% plot(t, Mphi_sus_dyn, 'LineWidth', 1.5); hold on; grid on;
% plot(t, Mphi_sus_qs, '--', 'LineWidth', 1.5);
% xlabel('t [s]');
% ylabel('M_{\phi,sus} [Nm]');
% legend('Dynamic','Quasi-static','Location','best');
% title('Roll suspension moment');
% 
% nexttile;
% plot(t, Mtheta_sus_dyn, 'LineWidth', 1.5); hold on; grid on;
% plot(t, Mtheta_sus_qs, '--', 'LineWidth', 1.5);
% xlabel('t [s]');
% ylabel('M_{\theta,sus} [Nm]');
% legend('Dynamic','Quasi-static','Location','best');
% title('Pitch suspension moment');
% 
% %% ------------------------------------------------------------------------
% % Plot wheel loads
% % -------------------------------------------------------------------------
% plot_corner_compare(t, Fz_FL_orig, Fz_FL_new_dyn, Fz_FL_new_qs, 'FL');
% plot_corner_compare(t, Fz_FR_orig, Fz_FR_new_dyn, Fz_FR_new_qs, 'FR');
% plot_corner_compare(t, Fz_RL_orig, Fz_RL_new_dyn, Fz_RL_new_qs, 'RL');
% plot_corner_compare(t, Fz_RR_orig, Fz_RR_new_dyn, Fz_RR_new_qs, 'RR');
% 
% %% ------------------------------------------------------------------------
% % Plot quasi-static matching errors
% % -------------------------------------------------------------------------
% figure('Name','Quasi-static matching error');
% tiledlayout(2,2);
% 
% nexttile;
% plot(t, err_FL, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('Error [N]');
% title('Fz_{FL,new,qs} - Fz_{FL,orig}');
% 
% nexttile;
% plot(t, err_FR, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('Error [N]');
% title('Fz_{FR,new,qs} - Fz_{FR,orig}');
% 
% nexttile;
% plot(t, err_RL, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('Error [N]');
% title('Fz_{RL,new,qs} - Fz_{RL,orig}');
% 
% nexttile;
% plot(t, err_RR, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('Error [N]');
% title('Fz_{RR,new,qs} - Fz_{RR,orig}');
% 
% %% ========================================================================
% % Local functions
% % ========================================================================
% 
% function dz = roll_pitch_ode(tt, z, t, ax_bar, ay_bar, veh)
%     % State order:
%     % z = [phi; phi_dot; theta; theta_dot]
% 
%     phi       = z(1);
%     phi_dot   = z(2);
%     theta     = z(3);
%     theta_dot = z(4);
% 
%     % Interpolate excitation signals
%     ax_bar_t = interp1(t, ax_bar, tt, 'pchip');
%     ay_bar_t = interp1(t, ay_bar, tt, 'pchip');
% 
%     % Sign convention:
%     %   Mphi_sus   > 0 loads RIGHT side
%     %   Mtheta_sus > 0 loads FRONT axle
%     %
%     % To recover the original wheel-load equations at steady state:
%     %   Mphi_sus_ss   = -m*h*ay_bar
%     %   Mtheta_sus_ss = -m*h*ax_bar
% 
%     phi_ddot = ( -veh.m*veh.hc*ay_bar_t - veh.Cphi*phi_dot - veh.Kphi*phi ) / veh.Ix;
%     theta_ddot = ( -veh.m*veh.hc*ax_bar_t - veh.Ctheta*theta_dot - veh.Ktheta*theta ) / veh.Iy;
% 
%     dz = [phi_dot; phi_ddot; theta_dot; theta_ddot];
% end
% 
% function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load_original( ...
%     Fx_tires, Fy_tires, Fdown, Fdrag, ...
%     bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h)
% 
%     Fz_FL = ((bcom*(mass*gravity+Fdown) ...
%             - (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
%             + ((Fy_tires*hcom)*(-dr)) / (twF_h*(-dr) - twR_h*(1-dr))) / 2;
% 
%     Fz_FR = ((bcom*(mass*gravity+Fdown) ...
%             - (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
%             - ((Fy_tires*hcom)*(-dr)) / (twF_h*(-dr) - twR_h*(1-dr))) / 2;
% 
%     Fz_RL = ((acom*(mass*gravity+Fdown) ...
%             + (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
%             + ((Fy_tires*hcom)*(1-dr)) / (twR_h*(1-dr) - (-dr)*twF_h)) / 2;
% 
%     Fz_RR = ((acom*(mass*gravity+Fdown) ...
%             + (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
%             - ((Fy_tires*hcom)*(1-dr)) / (twR_h*(1-dr) - (-dr)*twF_h)) / 2;
% end
% 
% function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load_new( ...
%     Mtheta_sus, Mphi_sus, Fdown, Fdrag, ...
%     bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h)
% 
%     L = acom + bcom;
%     Nz_total = mass*gravity + Fdown;
% 
%     % This external pitch term is the one needed so that the new quasi-static
%     % equations collapse exactly to the original equations.
%     Mtheta_ext = -(Fdrag*hcop + Fdown*(acop-acom));
% 
%     % Positive Mtheta_sus loads FRONT axle, unloads REAR axle
%     Fz_front_total = ( bcom*Nz_total + Mtheta_sus + Mtheta_ext ) / L;
%     Fz_rear_total  = ( acom*Nz_total - Mtheta_sus - Mtheta_ext ) / L;
% 
%     % Positive Mphi_sus loads RIGHT side, unloads LEFT side
%     Dphi = dr*twF_h + (1-dr)*twR_h;
% 
%     dF_front = dr     * Mphi_sus / Dphi;   % FR - FL
%     dF_rear  = (1-dr) * Mphi_sus / Dphi;   % RR - RL
% 
%     Fz_FL = 0.5 * (Fz_front_total - dF_front);
%     Fz_FR = 0.5 * (Fz_front_total + dF_front);
% 
%     Fz_RL = 0.5 * (Fz_rear_total  - dF_rear);
%     Fz_RR = 0.5 * (Fz_rear_total  + dF_rear);
% end
% 
% function plot_corner_compare(t, Fz_orig, Fz_new_dyn, Fz_new_qs, corner_name)
%     figure('Name',['Wheel Load Comparison - ' corner_name]);
%     plot(t, Fz_orig, 'LineWidth', 1.6); hold on; grid on;
%     plot(t, Fz_new_dyn, '--', 'LineWidth', 1.6);
%     plot(t, Fz_new_qs, ':', 'LineWidth', 2.0);
%     xlabel('t [s]');
%     ylabel(['Fz_' corner_name ' [N]']);
%     legend('Original','New dynamic','New quasi-static','Location','best');
%     title(['Wheel load comparison - ' corner_name]);
% end










% 
% %% compare_wheel_load_models_checked.m
% % Compare:
% %   1) original wheel-load equations
% %   2) new dynamic wheel-load equations
% %   3) new quasi-static wheel-load equations
% %
% % This version is checked for:
% %   - exact quasi-static collapse back to your original wheel_load(...)
% %   - dynamic pitch ODE including the direct aerodynamic pitch moment
% %   - aero consistency via vx = V*cos(beta)
% %   - optional quasi-static-consistent initial conditions for roll/pitch
% 
% clear
% close all
% clc
% 
% %% ------------------------------------------------------------------------
% % Parameters
% % -------------------------------------------------------------------------
% veh = my_params();
% 
% assert(isfield(veh,'roll_pitch_ode_initial_condition'), ...
%     'veh.roll_pitch_ode_initial_condition is missing.');
% assert(numel(veh.roll_pitch_ode_initial_condition) == 4, ...
%     'veh.roll_pitch_ode_initial_condition must be a 4x1 vector [phi; phi_dot; theta; theta_dot].');
% 
% %% ------------------------------------------------------------------------
% % USER INPUTS / DEMO INPUTS
% % Replace these with your own time histories if desired
% % -------------------------------------------------------------------------
% dt = 0.01;
% t_end = 12.0;
% t = (0:dt:t_end).';
% 
% % Demo vehicle speed
% V = 28 + 5*sin(2*pi*0.12*t) + 2*sin(2*pi*0.35*t + 0.4);  % [m/s]
% V = max(V, veh.V_min);
% 
% % Demo sideslip angle for aero consistency
% % Replace with your own beta if available; otherwise keep zeros
% beta = zeros(size(t));   % [rad]
% 
% % Demo filtered accelerations used by your original wheel-load equations
% ax_bar = 0.8*sin(2*pi*0.18*t) - 0.35*sin(2*pi*0.55*t + 0.3);   % [m/s^2]
% ay_bar = 4.5*sin(2*pi*0.16*t + 0.5) + 1.0*sin(2*pi*0.48*t);    % [m/s^2]
% 
% % Optional: load your own solved data
% % load('your_data.mat','t','V','beta','ax_bar','ay_bar');
% 
% assert(isvector(t) && isvector(V) && isvector(beta) && isvector(ax_bar) && isvector(ay_bar), ...
%     't, V, beta, ax_bar, ay_bar must be vectors.');
% assert(numel(t)==numel(V) && numel(t)==numel(beta) && ...
%        numel(t)==numel(ax_bar) && numel(t)==numel(ay_bar), ...
%     't, V, beta, ax_bar, ay_bar must have the same length.');
% 
% t      = t(:);
% V      = V(:);
% beta   = beta(:);
% ax_bar = ax_bar(:);
% ay_bar = ay_bar(:);
% 
% %% ------------------------------------------------------------------------
% % Geometry / aero points used exactly like your current my_model.m
% % -------------------------------------------------------------------------
% acom = veh.lf;
% bcom = veh.lr;
% hcom = veh.hc;
% 
% % In your current my_model.m these are hard-coded as:
% acop = acom;
% hcop = hcom;
% 
% %% ------------------------------------------------------------------------
% % Aero quantities
% % -------------------------------------------------------------------------
% % Consistent with your full model: vx = V*cos(beta)
% vx_for_aero = V .* cos(beta);
% 
% f_drag = 0.5 * veh.drag_coeff * vx_for_aero.^2;   % [N]
% f_lift = 0.5 * veh.lift_coeff * vx_for_aero.^2;   % [N]
% 
% %% ------------------------------------------------------------------------
% % ORIGINAL wheel-load equations
% % -------------------------------------------------------------------------
% Fx_tires = veh.m * ax_bar + f_drag;   % exactly as in your current model
% Fy_tires = veh.m * ay_bar;            % exactly as in your current model
% 
% Fz_FL_orig = zeros(size(t));
% Fz_FR_orig = zeros(size(t));
% Fz_RL_orig = zeros(size(t));
% Fz_RR_orig = zeros(size(t));
% 
% for k = 1:numel(t)
%     [Fz_FL_orig(k), Fz_FR_orig(k), Fz_RL_orig(k), Fz_RR_orig(k)] = ...
%         wheel_load_original( ...
%             Fx_tires(k), Fy_tires(k), f_lift(k), f_drag(k), ...
%             bcom, veh.m, veh.g, hcom, hcop, acop, acom, ...
%             veh.dr, veh.twF_h, veh.twR_h);
% end
% 
% %% ------------------------------------------------------------------------
% % NEW QUASI-STATIC wheel-load equations
% % These should match the original equations to numerical roundoff
% % -------------------------------------------------------------------------
% Mphi_sus_qs   = -veh.m * hcom * ay_bar;
% Mtheta_sus_qs = -veh.m * hcom * ax_bar;
% 
% Fz_FL_new_qs = zeros(size(t));
% Fz_FR_new_qs = zeros(size(t));
% Fz_RL_new_qs = zeros(size(t));
% Fz_RR_new_qs = zeros(size(t));
% 
% for k = 1:numel(t)
%     [Fz_FL_new_qs(k), Fz_FR_new_qs(k), Fz_RL_new_qs(k), Fz_RR_new_qs(k)] = ...
%         wheel_load_new( ...
%             Mtheta_sus_qs(k), Mphi_sus_qs(k), f_lift(k), f_drag(k), ...
%             bcom, veh.m, veh.g, hcom, hcop, acop, acom, ...
%             veh.dr, veh.twF_h, veh.twR_h);
% end
% 
% %% ------------------------------------------------------------------------
% % NEW DYNAMIC wheel-load equations
% % Integrate roll/pitch ODEs driven by ax_bar, ay_bar, and aero pitch moment
% % -------------------------------------------------------------------------
% use_quasi_static_initial_condition = true;
% 
% if use_quasi_static_initial_condition
%     phi0       = Mphi_sus_qs(1)   / veh.Kphi;
%     phi_dot0   = 0;
%     theta0     = Mtheta_sus_qs(1) / veh.Ktheta;
%     theta_dot0 = 0;
%     z0 = [phi0; phi_dot0; theta0; theta_dot0];
% else
%     z0 = veh.roll_pitch_ode_initial_condition(:);  % [phi; phi_dot; theta; theta_dot]
% end
% 
% ode_fun = @(tt,zz) roll_pitch_ode_checked(tt, zz, ...
%     t, V, beta, ax_bar, ay_bar, veh, hcom, hcop, acop, acom);
% 
% opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
% [t_ode, z_ode] = ode45(ode_fun, [t(1), t(end)], z0, opts);
% 
% phi       = interp1(t_ode, z_ode(:,1), t, 'linear', 'extrap');
% phi_dot   = interp1(t_ode, z_ode(:,2), t, 'linear', 'extrap');
% theta     = interp1(t_ode, z_ode(:,3), t, 'linear', 'extrap');
% theta_dot = interp1(t_ode, z_ode(:,4), t, 'linear', 'extrap');
% 
% Mphi_sus_dyn   = veh.Cphi   * phi_dot   + veh.Kphi   * phi;
% Mtheta_sus_dyn = veh.Ctheta * theta_dot + veh.Ktheta * theta;
% 
% Fz_FL_new_dyn = zeros(size(t));
% Fz_FR_new_dyn = zeros(size(t));
% Fz_RL_new_dyn = zeros(size(t));
% Fz_RR_new_dyn = zeros(size(t));
% 
% for k = 1:numel(t)
%     [Fz_FL_new_dyn(k), Fz_FR_new_dyn(k), Fz_RL_new_dyn(k), Fz_RR_new_dyn(k)] = ...
%         wheel_load_new( ...
%             Mtheta_sus_dyn(k), Mphi_sus_dyn(k), f_lift(k), f_drag(k), ...
%             bcom, veh.m, veh.g, hcom, hcop, acop, acom, ...
%             veh.dr, veh.twF_h, veh.twR_h);
% end
% 
% %% ------------------------------------------------------------------------
% % Errors: original vs new quasi-static
% % -------------------------------------------------------------------------
% err_FL = Fz_FL_new_qs - Fz_FL_orig;
% err_FR = Fz_FR_new_qs - Fz_FR_orig;
% err_RL = Fz_RL_new_qs - Fz_RL_orig;
% err_RR = Fz_RR_new_qs - Fz_RR_orig;
% 
% fprintf('\n=============================================================\n');
% fprintf('Max abs difference: NEW quasi-static vs ORIGINAL\n');
% fprintf('Fz_FL: %.12e N\n', max(abs(err_FL)));
% fprintf('Fz_FR: %.12e N\n', max(abs(err_FR)));
% fprintf('Fz_RL: %.12e N\n', max(abs(err_RL)));
% fprintf('Fz_RR: %.12e N\n', max(abs(err_RR)));
% fprintf('=============================================================\n\n');
% 
% %% ------------------------------------------------------------------------
% % Sanity checks
% % -------------------------------------------------------------------------
% sum_orig = Fz_FL_orig + Fz_FR_orig + Fz_RL_orig + Fz_RR_orig;
% sum_dyn  = Fz_FL_new_dyn + Fz_FR_new_dyn + Fz_RL_new_dyn + Fz_RR_new_dyn;
% sum_qs   = Fz_FL_new_qs  + Fz_FR_new_qs  + Fz_RL_new_qs  + Fz_RR_new_qs;
% 
% figure('Name','Total vertical load check');
% plot(t, sum_orig, 'LineWidth', 1.5); hold on; grid on;
% plot(t, sum_dyn,  '--', 'LineWidth', 1.5);
% plot(t, sum_qs,   ':', 'LineWidth', 1.8);
% xlabel('t [s]');
% ylabel('Fz_{FL}+Fz_{FR}+Fz_{RL}+Fz_{RR} [N]');
% legend('Original','New dynamic','New quasi-static','Location','best');
% title('Total vertical load');
% 
% %% ------------------------------------------------------------------------
% % Plot roll / pitch states and moments
% % -------------------------------------------------------------------------
% figure('Name','Roll and Pitch States');
% tiledlayout(2,2);
% 
% nexttile;
% plot(t, phi, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('\phi [rad]');
% title('Roll angle');
% 
% nexttile;
% plot(t, phi_dot, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('\phi_{dot} [rad/s]');
% title('Roll rate');
% 
% nexttile;
% plot(t, theta, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('\theta [rad]');
% title('Pitch angle');
% 
% nexttile;
% plot(t, theta_dot, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('\theta_{dot} [rad/s]');
% title('Pitch rate');
% 
% figure('Name','Suspension Moments');
% tiledlayout(2,1);
% 
% nexttile;
% plot(t, Mphi_sus_dyn, 'LineWidth', 1.5); hold on; grid on;
% plot(t, Mphi_sus_qs, '--', 'LineWidth', 1.5);
% xlabel('t [s]');
% ylabel('M_{\phi,sus} [Nm]');
% legend('Dynamic','Quasi-static','Location','best');
% title('Roll suspension moment');
% 
% nexttile;
% plot(t, Mtheta_sus_dyn, 'LineWidth', 1.5); hold on; grid on;
% plot(t, Mtheta_sus_qs, '--', 'LineWidth', 1.5);
% xlabel('t [s]');
% ylabel('M_{\theta,sus} [Nm]');
% legend('Dynamic','Quasi-static','Location','best');
% title('Pitch suspension moment');
% 
% %% ------------------------------------------------------------------------
% % Plot wheel loads
% % -------------------------------------------------------------------------
% plot_corner_compare(t, Fz_FL_orig, Fz_FL_new_dyn, Fz_FL_new_qs, 'FL');
% plot_corner_compare(t, Fz_FR_orig, Fz_FR_new_dyn, Fz_FR_new_qs, 'FR');
% plot_corner_compare(t, Fz_RL_orig, Fz_RL_new_dyn, Fz_RL_new_qs, 'RL');
% plot_corner_compare(t, Fz_RR_orig, Fz_RR_new_dyn, Fz_RR_new_qs, 'RR');
% 
% %% ------------------------------------------------------------------------
% % Plot quasi-static matching errors
% % -------------------------------------------------------------------------
% figure('Name','Quasi-static matching error');
% tiledlayout(2,2);
% 
% nexttile;
% plot(t, err_FL, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('Error [N]');
% title('Fz_{FL,new,qs} - Fz_{FL,orig}');
% 
% nexttile;
% plot(t, err_FR, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('Error [N]');
% title('Fz_{FR,new,qs} - Fz_{FR,orig}');
% 
% nexttile;
% plot(t, err_RL, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('Error [N]');
% title('Fz_{RL,new,qs} - Fz_{RL,orig}');
% 
% nexttile;
% plot(t, err_RR, 'LineWidth', 1.5); grid on;
% xlabel('t [s]'); ylabel('Error [N]');
% title('Fz_{RR,new,qs} - Fz_{RR,orig}');
% 
% %% ========================================================================
% % Local functions
% % ========================================================================
% 
% function dz = roll_pitch_ode_checked(tt, z, t, V, beta, ax_bar, ay_bar, veh, hcom, hcop, acop, acom)
%     % State order:
%     % z = [phi; phi_dot; theta; theta_dot]
% 
%     phi       = z(1);
%     phi_dot   = z(2);
%     theta     = z(3);
%     theta_dot = z(4);
% 
%     % Interpolate excitation signals
%     V_t      = interp1(t, V,      tt, 'linear', 'extrap');
%     beta_t   = interp1(t, beta,   tt, 'linear', 'extrap');
%     ax_bar_t = interp1(t, ax_bar, tt, 'linear', 'extrap');
%     ay_bar_t = interp1(t, ay_bar, tt, 'linear', 'extrap');
% 
%     % Aero at current instant
%     vx_t    = V_t * cos(beta_t);
%     Fdrag_t = 0.5 * veh.drag_coeff * vx_t^2;
%     Fdown_t = 0.5 * veh.lift_coeff * vx_t^2;
% 
%     % Positive Mtheta_sus loads FRONT axle.
%     % This same external pitch moment must appear in the pitch ODE if the
%     % dynamic pitch state is to be physically consistent with wheel_load_new.
%     Mtheta_ext_t = -(Fdrag_t*hcop + Fdown_t*(acop-acom));
% 
%     % Positive Mphi_sus loads RIGHT side.
%     % To recover the original wheel-load equations at steady state:
%     %   Mphi_sus_ss   = -m*h*ay_bar
%     %   Mtheta_sus_ss = -m*h*ax_bar
%     %
%     % Therefore:
%     %   Ix*phi_ddot   = -m*h*ay_bar - Cphi*phi_dot - Kphi*phi
%     %   Iy*theta_ddot = -m*h*ax_bar + Mtheta_ext - Ctheta*theta_dot - Ktheta*theta
% 
%     phi_ddot   = ( -veh.m*hcom*ay_bar_t ...
%                    - veh.Cphi*phi_dot - veh.Kphi*phi ) / veh.Ix;
% 
%     theta_ddot = ( -veh.m*hcom*ax_bar_t + Mtheta_ext_t ...
%                    - veh.Ctheta*theta_dot - veh.Ktheta*theta ) / veh.Iy;
% 
%     dz = [phi_dot; phi_ddot; theta_dot; theta_ddot];
% end
% 
% function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load_original( ...
%     Fx_tires, Fy_tires, Fdown, Fdrag, ...
%     bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h)
% 
%     Fz_FL = ((bcom*(mass*gravity+Fdown) ...
%             - (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
%             + ((Fy_tires*hcom)*(-dr)) / (twF_h*(-dr) - twR_h*(1-dr))) / 2;
% 
%     Fz_FR = ((bcom*(mass*gravity+Fdown) ...
%             - (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
%             - ((Fy_tires*hcom)*(-dr)) / (twF_h*(-dr) - twR_h*(1-dr))) / 2;
% 
%     Fz_RL = ((acom*(mass*gravity+Fdown) ...
%             + (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
%             + ((Fy_tires*hcom)*(1-dr)) / (twR_h*(1-dr) - (-dr)*twF_h)) / 2;
% 
%     Fz_RR = ((acom*(mass*gravity+Fdown) ...
%             + (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
%             - ((Fy_tires*hcom)*(1-dr)) / (twR_h*(1-dr) - (-dr)*twF_h)) / 2;
% end
% 
% function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load_new( ...
%     Mtheta_sus, Mphi_sus, Fdown, Fdrag, ...
%     bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h)
% 
%     L = acom + bcom;
%     Nz_total = mass*gravity + Fdown;
% 
%     % Positive Mtheta_sus loads FRONT axle, unloads REAR axle
%     % This external moment is required for exact quasi-static recovery of
%     % the original wheel_load(...) equations.
%     Mtheta_ext = -(Fdrag*hcop + Fdown*(acop-acom));
% 
%     Fz_front_total = ( bcom*Nz_total + Mtheta_sus + Mtheta_ext ) / L;
%     Fz_rear_total  = ( acom*Nz_total - Mtheta_sus - Mtheta_ext ) / L;
% 
%     % Positive Mphi_sus loads RIGHT side, unloads LEFT side
%     Dphi = dr*twF_h + (1-dr)*twR_h;
% 
%     dF_front = dr     * Mphi_sus / Dphi;   % FR - FL
%     dF_rear  = (1-dr) * Mphi_sus / Dphi;   % RR - RL
% 
%     Fz_FL = 0.5 * (Fz_front_total - dF_front);
%     Fz_FR = 0.5 * (Fz_front_total + dF_front);
% 
%     Fz_RL = 0.5 * (Fz_rear_total  - dF_rear);
%     Fz_RR = 0.5 * (Fz_rear_total  + dF_rear);
% end
% 
% function plot_corner_compare(t, Fz_orig, Fz_new_dyn, Fz_new_qs, corner_name)
%     figure('Name',['Wheel Load Comparison - ' corner_name]);
%     plot(t, Fz_orig, 'LineWidth', 1.6); hold on; grid on;
%     plot(t, Fz_new_dyn, '--', 'LineWidth', 1.6);
%     plot(t, Fz_new_qs, ':', 'LineWidth', 2.0);
%     xlabel('t [s]');
%     ylabel(['Fz_' corner_name ' [N]']);
%     legend('Original','New dynamic','New quasi-static','Location','best');
%     title(['Wheel load comparison - ' corner_name]);
% end










%% compare_wheel_load_models_fully_checked.m
% Compare:
%   1) original wheel-load equations
%   2) new dynamic wheel-load equations
%   3) new quasi-static wheel-load equations
%
% Fully checked points:
%   - original equations are exactly your current equations
%   - new quasi-static equations collapse exactly to the original equations
%   - new dynamic equations are now consistent with the new wheel-load definition
%   - no double-counting of external pitch moment
%   - quasi-static initial conditions are defined consistently

clear
close all
clc

%% ------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
veh = my_params();

assert(isfield(veh,'roll_pitch_ode_initial_condition'), ...
    'veh.roll_pitch_ode_initial_condition is missing.');
assert(numel(veh.roll_pitch_ode_initial_condition) == 4, ...
    'veh.roll_pitch_ode_initial_condition must be a 4x1 vector [phi; phi_dot; theta; theta_dot].');

%% ------------------------------------------------------------------------
% USER INPUTS / DEMO INPUTS
% Replace these with your own time histories if desired
% -------------------------------------------------------------------------
dt = 0.01;
t_end = 12.0;
t = (0:dt:t_end).';

% Demo vehicle speed
V = 28 + 5*sin(2*pi*0.12*t) + 2*sin(2*pi*0.35*t + 0.4);  % [m/s]
V = max(V, veh.V_min);

% Demo sideslip angle for aero consistency
% Replace with your own beta if available
beta = zeros(size(t));   % [rad]

% Demo filtered accelerations used by your original wheel-load equations
ax_bar = 0.8*sin(2*pi*0.18*t) - 0.35*sin(2*pi*0.55*t + 0.3);   % [m/s^2]
ay_bar = 4.5*sin(2*pi*0.16*t + 0.5) + 1.0*sin(2*pi*0.48*t);    % [m/s^2]

% Optional: load your own solved data
% load('your_data.mat','t','V','beta','ax_bar','ay_bar');

assert(isvector(t) && isvector(V) && isvector(beta) && isvector(ax_bar) && isvector(ay_bar), ...
    't, V, beta, ax_bar, ay_bar must be vectors.');
assert(numel(t)==numel(V) && numel(t)==numel(beta) && ...
       numel(t)==numel(ax_bar) && numel(t)==numel(ay_bar), ...
    't, V, beta, ax_bar, ay_bar must have the same length.');

t      = t(:);
V      = V(:);
beta   = beta(:);
ax_bar = ax_bar(:);
ay_bar = ay_bar(:);

%% ------------------------------------------------------------------------
% Geometry / aero points used exactly like your current my_model.m
% -------------------------------------------------------------------------
acom = veh.lf;
bcom = veh.lr;
hcom = veh.hc;

% In your current my_model.m these are hard-coded as:
acop = acom;
hcop = hcom;

%% ------------------------------------------------------------------------
% Aero quantities
% -------------------------------------------------------------------------
% Consistent with your full model: vx = V*cos(beta)
vx_for_aero = V .* cos(beta);

f_drag = 0.5 * veh.drag_coeff * vx_for_aero.^2;   % [N]
f_lift = 0.5 * veh.lift_coeff * vx_for_aero.^2;   % [N]

% External pitch moment that is not generated by suspension deflection itself
% Positive Mtheta_sus loads FRONT axle, so this term is front-unloading.
Mtheta_ext = -(f_drag .* hcop + f_lift .* (acop - acom));

%% ------------------------------------------------------------------------
% ORIGINAL wheel-load equations
% -------------------------------------------------------------------------
Fx_tires = veh.m * ax_bar + f_drag;   % exactly as in your current model
Fy_tires = veh.m * ay_bar;            % exactly as in your current model

Fz_FL_orig = zeros(size(t));
Fz_FR_orig = zeros(size(t));
Fz_RL_orig = zeros(size(t));
Fz_RR_orig = zeros(size(t));

for k = 1:numel(t)
    [Fz_FL_orig(k), Fz_FR_orig(k), Fz_RL_orig(k), Fz_RR_orig(k)] = ...
        wheel_load_original( ...
            Fx_tires(k), Fy_tires(k), f_lift(k), f_drag(k), ...
            bcom, veh.m, veh.g, hcom, hcop, acop, acom, ...
            veh.dr, veh.twF_h, veh.twR_h);
end

%% ------------------------------------------------------------------------
% NEW QUASI-STATIC wheel-load equations
% These should match the original equations to numerical roundoff
% -------------------------------------------------------------------------
Mphi_sus_qs   = -veh.m * hcom * ay_bar;
Mtheta_sus_qs = -veh.m * hcom * ax_bar + Mtheta_ext;

Fz_FL_new_qs = zeros(size(t));
Fz_FR_new_qs = zeros(size(t));
Fz_RL_new_qs = zeros(size(t));
Fz_RR_new_qs = zeros(size(t));

for k = 1:numel(t)
    [Fz_FL_new_qs(k), Fz_FR_new_qs(k), Fz_RL_new_qs(k), Fz_RR_new_qs(k)] = ...
        wheel_load_new( ...
            Mtheta_sus_qs(k), Mphi_sus_qs(k), f_lift(k), ...
            bcom, veh.m, veh.g, acom, veh.dr, veh.twF_h, veh.twR_h);
end

%% ------------------------------------------------------------------------
% NEW DYNAMIC wheel-load equations
% Integrate roll/pitch ODEs driven by ax_bar, ay_bar, and external pitch moment
% -------------------------------------------------------------------------
use_quasi_static_initial_condition = true;

if use_quasi_static_initial_condition
    phi0       = Mphi_sus_qs(1)   / veh.Kphi;
    phi_dot0   = 0;
    theta0     = Mtheta_sus_qs(1) / veh.Ktheta;
    theta_dot0 = 0;
    z0 = [phi0; phi_dot0; theta0; theta_dot0];
else
    z0 = veh.roll_pitch_ode_initial_condition(:);  % [phi; phi_dot; theta; theta_dot]
end

ode_fun = @(tt,zz) roll_pitch_ode_fully_checked(tt, zz, ...
    t, V, beta, ax_bar, ay_bar, veh, hcom, hcop, acop, acom);

opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
[t_ode, z_ode] = ode45(ode_fun, [t(1), t(end)], z0, opts);

phi       = interp1(t_ode, z_ode(:,1), t, 'linear', 'extrap');
phi_dot   = interp1(t_ode, z_ode(:,2), t, 'linear', 'extrap');
theta     = interp1(t_ode, z_ode(:,3), t, 'linear', 'extrap');
theta_dot = interp1(t_ode, z_ode(:,4), t, 'linear', 'extrap');

Mphi_sus_dyn   = veh.Cphi   * phi_dot   + veh.Kphi   * phi;
Mtheta_sus_dyn = veh.Ctheta * theta_dot + veh.Ktheta * theta;

Fz_FL_new_dyn = zeros(size(t));
Fz_FR_new_dyn = zeros(size(t));
Fz_RL_new_dyn = zeros(size(t));
Fz_RR_new_dyn = zeros(size(t));

for k = 1:numel(t)
    [Fz_FL_new_dyn(k), Fz_FR_new_dyn(k), Fz_RL_new_dyn(k), Fz_RR_new_dyn(k)] = ...
        wheel_load_new( ...
            Mtheta_sus_dyn(k), Mphi_sus_dyn(k), f_lift(k), ...
            bcom, veh.m, veh.g, acom, veh.dr, veh.twF_h, veh.twR_h);
end

%% ------------------------------------------------------------------------
% Errors: original vs new quasi-static
% -------------------------------------------------------------------------
err_FL = Fz_FL_new_qs - Fz_FL_orig;
err_FR = Fz_FR_new_qs - Fz_FR_orig;
err_RL = Fz_RL_new_qs - Fz_RL_orig;
err_RR = Fz_RR_new_qs - Fz_RR_orig;

fprintf('\n=============================================================\n');
fprintf('Max abs difference: NEW quasi-static vs ORIGINAL\n');
fprintf('Fz_FL: %.12e N\n', max(abs(err_FL)));
fprintf('Fz_FR: %.12e N\n', max(abs(err_FR)));
fprintf('Fz_RL: %.12e N\n', max(abs(err_RL)));
fprintf('Fz_RR: %.12e N\n', max(abs(err_RR)));
fprintf('=============================================================\n\n');

%% ------------------------------------------------------------------------
% Sanity checks
% -------------------------------------------------------------------------
sum_orig = Fz_FL_orig + Fz_FR_orig + Fz_RL_orig + Fz_RR_orig;
sum_dyn  = Fz_FL_new_dyn + Fz_FR_new_dyn + Fz_RL_new_dyn + Fz_RR_new_dyn;
sum_qs   = Fz_FL_new_qs  + Fz_FR_new_qs  + Fz_RL_new_qs  + Fz_RR_new_qs;

figure('Name','Total vertical load check');
plot(t, sum_orig, 'LineWidth', 1.5); hold on; grid on;
plot(t, sum_dyn,  '--', 'LineWidth', 1.5);
plot(t, sum_qs,   ':', 'LineWidth', 1.8);
xlabel('t [s]');
ylabel('Fz_{FL}+Fz_{FR}+Fz_{RL}+Fz_{RR} [N]');
legend('Original','New dynamic','New quasi-static','Location','best');
title('Total vertical load');

%% ------------------------------------------------------------------------
% Plot roll / pitch states and moments
% -------------------------------------------------------------------------
figure('Name','Roll and Pitch States');
tiledlayout(2,2);

nexttile;
plot(t, phi, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('\phi [rad]');
title('Roll angle');

nexttile;
plot(t, phi_dot, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('\phi_{dot} [rad/s]');
title('Roll rate');

nexttile;
plot(t, theta, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('\theta [rad]');
title('Pitch angle');

nexttile;
plot(t, theta_dot, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('\theta_{dot} [rad/s]');
title('Pitch rate');

figure('Name','Suspension Moments');
tiledlayout(2,1);

nexttile;
plot(t, Mphi_sus_dyn, 'LineWidth', 1.5); hold on; grid on;
plot(t, Mphi_sus_qs, '--', 'LineWidth', 1.5);
xlabel('t [s]');
ylabel('M_{\phi,sus} [Nm]');
legend('Dynamic','Quasi-static','Location','best');
title('Roll suspension moment');

nexttile;
plot(t, Mtheta_sus_dyn, 'LineWidth', 1.5); hold on; grid on;
plot(t, Mtheta_sus_qs, '--', 'LineWidth', 1.5);
xlabel('t [s]');
ylabel('M_{\theta,sus} [Nm]');
legend('Dynamic','Quasi-static','Location','best');
title('Pitch suspension moment');

%% ------------------------------------------------------------------------
% Plot wheel loads
% -------------------------------------------------------------------------
plot_corner_compare(t, Fz_FL_orig, Fz_FL_new_dyn, Fz_FL_new_qs, 'FL');
plot_corner_compare(t, Fz_FR_orig, Fz_FR_new_dyn, Fz_FR_new_qs, 'FR');
plot_corner_compare(t, Fz_RL_orig, Fz_RL_new_dyn, Fz_RL_new_qs, 'RL');
plot_corner_compare(t, Fz_RR_orig, Fz_RR_new_dyn, Fz_RR_new_qs, 'RR');

%% ------------------------------------------------------------------------
% Plot quasi-static matching errors
% -------------------------------------------------------------------------
figure('Name','Quasi-static matching error');
tiledlayout(2,2);

nexttile;
plot(t, err_FL, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('Error [N]');
title('Fz_{FL,new,qs} - Fz_{FL,orig}');

nexttile;
plot(t, err_FR, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('Error [N]');
title('Fz_{FR,new,qs} - Fz_{FR,orig}');

nexttile;
plot(t, err_RL, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('Error [N]');
title('Fz_{RL,new,qs} - Fz_{RL,orig}');

nexttile;
plot(t, err_RR, 'LineWidth', 1.5); grid on;
xlabel('t [s]'); ylabel('Error [N]');
title('Fz_{RR,new,qs} - Fz_{RR,orig}');

%% ========================================================================
% Local functions
% ========================================================================

function dz = roll_pitch_ode_fully_checked(tt, z, t, V, beta, ax_bar, ay_bar, veh, hcom, hcop, acop, acom)
    % State order:
    % z = [phi; phi_dot; theta; theta_dot]

    phi       = z(1);
    phi_dot   = z(2);
    theta     = z(3);
    theta_dot = z(4);

    % Interpolate signals
    V_t      = interp1(t, V,      tt, 'linear', 'extrap');
    beta_t   = interp1(t, beta,   tt, 'linear', 'extrap');
    ax_bar_t = interp1(t, ax_bar, tt, 'linear', 'extrap');
    ay_bar_t = interp1(t, ay_bar, tt, 'linear', 'extrap');

    % Aero at current instant
    vx_t    = V_t * cos(beta_t);
    Fdrag_t = 0.5 * veh.drag_coeff * vx_t^2;
    Fdown_t = 0.5 * veh.lift_coeff * vx_t^2;

    % External pitch moment (front-unloading)
    Mtheta_ext_t = -(Fdrag_t * hcop + Fdown_t * (acop - acom));

    % Positive Mphi_sus loads RIGHT side.
    % Positive Mtheta_sus loads FRONT axle.
    %
    % Consistent dynamic definitions:
    %   Ix*phi_ddot   = -m*h*ay_bar - Cphi*phi_dot - Kphi*phi
    %   Iy*theta_ddot = -m*h*ax_bar + Mtheta_ext - Ctheta*theta_dot - Ktheta*theta

    phi_ddot = ( -veh.m*hcom*ay_bar_t ...
                 - veh.Cphi*phi_dot - veh.Kphi*phi ) / veh.Ix;

    theta_ddot = ( -veh.m*hcom*ax_bar_t + Mtheta_ext_t ...
                   - veh.Ctheta*theta_dot - veh.Ktheta*theta ) / veh.Iy;

    dz = [phi_dot; phi_ddot; theta_dot; theta_ddot];
end

function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load_original( ...
    Fx_tires, Fy_tires, Fdown, Fdrag, ...
    bcom, mass, gravity, hcom, hcop, acop, acom, dr, twF_h, twR_h)

    Fz_FL = ((bcom*(mass*gravity+Fdown) ...
            - (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
            + ((Fy_tires*hcom)*(-dr)) / (twF_h*(-dr) - twR_h*(1-dr))) / 2;

    Fz_FR = ((bcom*(mass*gravity+Fdown) ...
            - (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
            - ((Fy_tires*hcom)*(-dr)) / (twF_h*(-dr) - twR_h*(1-dr))) / 2;

    Fz_RL = ((acom*(mass*gravity+Fdown) ...
            + (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
            + ((Fy_tires*hcom)*(1-dr)) / (twR_h*(1-dr) - (-dr)*twF_h)) / 2;

    Fz_RR = ((acom*(mass*gravity+Fdown) ...
            + (Fx_tires*hcom + Fdrag*(hcop-hcom) + Fdown*(acop-acom))) / (acom+bcom) ...
            - ((Fy_tires*hcom)*(1-dr)) / (twR_h*(1-dr) - (-dr)*twF_h)) / 2;
end

function [Fz_FL, Fz_FR, Fz_RL, Fz_RR] = wheel_load_new( ...
    Mtheta_sus, Mphi_sus, Fdown, ...
    bcom, mass, gravity, acom, dr, twF_h, twR_h)

    L = acom + bcom;
    Nz_total = mass*gravity + Fdown;

    % Positive Mtheta_sus loads FRONT axle, unloads REAR axle
    Fz_front_total = ( bcom * Nz_total + Mtheta_sus ) / L;
    Fz_rear_total  = ( acom * Nz_total - Mtheta_sus ) / L;

    % Positive Mphi_sus loads RIGHT side, unloads LEFT side
    Dphi = dr*twF_h + (1-dr)*twR_h;

    dF_front = dr     * Mphi_sus / Dphi;   % FR - FL
    dF_rear  = (1-dr) * Mphi_sus / Dphi;   % RR - RL

    Fz_FL = 0.5 * (Fz_front_total - dF_front);
    Fz_FR = 0.5 * (Fz_front_total + dF_front);

    Fz_RL = 0.5 * (Fz_rear_total  - dF_rear);
    Fz_RR = 0.5 * (Fz_rear_total  + dF_rear);
end

function plot_corner_compare(t, Fz_orig, Fz_new_dyn, Fz_new_qs, corner_name)
    figure('Name',['Wheel Load Comparison - ' corner_name]);
    plot(t, Fz_orig, 'LineWidth', 1.6); hold on; grid on;
    plot(t, Fz_new_dyn, '--', 'LineWidth', 1.6);
    plot(t, Fz_new_qs, ':', 'LineWidth', 2.0);
    xlabel('t [s]');
    ylabel(['Fz_' corner_name ' [N]']);
    legend('Original','New dynamic','New quasi-static','Location','best');
    title(['Wheel load comparison - ' corner_name]);
end