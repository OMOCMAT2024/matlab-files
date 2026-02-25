close all
clear
clc

%% Parameter setup
%
sa = 3; % do not change this!

nx = 14; % total number of states of my ribbon model

shift_X = 0; % shift the X coordinates so that (0,0) is somewhere within the track

p_x = 1e-1; % smooth parameter for raw x coordinates
p_y = 1e-1; % smooth parameter for raw y coordinates
p_z = (1e-4)/2;  % smooth parameter for raw z coordinates

p = 1e-1; % smooth parameter for smoothed raw coordinates

N_integration_grid = 10^8 * 1;

N_nodes = 5452*2; % total number of nodes including both the start and end points (no intermediate two collocation points)

%% Radua IIA
%
[A_matrix_RadauIIA, ~, c_vec_RadauIIA] = RadauIIA_tabular();

c1_RadauIIA = c_vec_RadauIIA(1);
c2_RadauIIA = c_vec_RadauIIA(2);
c3_RadauIIA = c_vec_RadauIIA(3);

%% Load Google Earth hand picked data
%
% load C_Cartesian_coords_GE_raw
% load CL_Cartesian_coords_GE_raw
% load CR_Cartesian_coords_GE_raw

load C_Cartesian_coords_Blender_raw
load CL_Cartesian_coords_Blender_raw
load CR_Cartesian_coords_Blender_raw
xEast_C_GE_raw = xEast_C_Blender_raw;
yNorth_C_GE_raw = yNorth_C_Blender_raw;
zUp_C_GE_raw = zUp_C_Blender_raw;
xEast_CL_GE_raw = xEast_CL_Blender_raw;
yNorth_CL_GE_raw = yNorth_CL_Blender_raw;
zUp_CL_GE_raw = zUp_CL_Blender_raw;
xEast_CR_GE_raw = xEast_CR_Blender_raw;
yNorth_CR_GE_raw = yNorth_CR_Blender_raw;
zUp_CR_GE_raw = zUp_CR_Blender_raw;

% Assemble raw points
rawC = [ xEast_C_GE_raw(:)+shift_X,  yNorth_C_GE_raw(:),  zUp_C_GE_raw(:) ];
rawL = [ xEast_CL_GE_raw(:)+shift_X, yNorth_CL_GE_raw(:), zUp_CL_GE_raw(:) ];
rawR = [ xEast_CR_GE_raw(:)+shift_X, yNorth_CR_GE_raw(:), zUp_CR_GE_raw(:) ];

figure()
plot(rawC(:,1), rawC(:,2), 'k-o')
hold on
plot(rawL(:,1), rawL(:,2), 'r-o')
hold on
plot(rawR(:,1), rawR(:,2), 'b-o')
hold on
plot(rawC(1,1), rawC(1,2), 'k-*')
hold on
plot(rawC(end,1), rawC(end,2), 'k-s')
xlabel('X, m')
ylabel('Y, m')
title('Google Earth raw data with X shift')
axis equal; grid on

%% Piecewise-length parameterization of raw Google Earth data
%
rawdC = diff(rawC,1,1);
rawUC = [0; cumsum( sqrt(sum(rawdC.^2,2)) )];

rawdL = diff(rawL,1,1);
rawUL = [0; cumsum( sqrt(sum(rawdL.^2,2)) )];

rawdR = diff(rawR,1,1);
rawUR = [0; cumsum( sqrt(sum(rawdR.^2,2)) )];

%% Smooth raw z data
%
my_pp_Cz = myspcsp(rawUC.', zUp_C_GE_raw.', p_z);
Cz_raw_smoothed = ppval(my_pp_Cz, rawUC);

my_pp_CLz = myspcsp(rawUL.', zUp_CL_GE_raw.', p_z);
CLz_raw_smoothed = ppval(my_pp_CLz, rawUL);

my_pp_CRz = myspcsp(rawUR.', zUp_CR_GE_raw.', p_z);
CRz_raw_smoothed = ppval(my_pp_CRz, rawUR);

figure()
subplot(3,1,1)
plot(rawUC, zUp_C_GE_raw(:))
hold on
plot(rawUC, Cz_raw_smoothed(:))
xlabel('rawUC, m')
ylabel('zUP C GE raw, m')
legend('raw Cz', 'Cz smoothed')
title('zUP C GE raw')
subplot(3,1,2)
plot(rawUL, zUp_CL_GE_raw(:))
hold on
plot(rawUL, CLz_raw_smoothed(:))
xlabel('rawUL, m')
ylabel('zUP L GE raw, m')
legend('raw CLz', 'CLz smoothed')
title('zUP L GE raw')
subplot(3,1,3)
plot(rawUR, zUp_CR_GE_raw(:))
hold on
plot(rawUR, CRz_raw_smoothed(:))
xlabel('rawUR, m')
ylabel('zUP R GE raw, m')
legend('raw CRz', 'CRz smoothed')
title('zUP R GE raw')

%% Smooth raw x data
%
my_pp_Cx = myspcsp(rawUC.', xEast_C_GE_raw.', p_x);
Cx_raw_smoothed = ppval(my_pp_Cx, rawUC);

my_pp_CLx = myspcsp(rawUL.', xEast_CL_GE_raw.', p_x);
CLx_raw_smoothed = ppval(my_pp_CLx, rawUL);

my_pp_CRx = myspcsp(rawUR.', xEast_CR_GE_raw.', p_x);
CRx_raw_smoothed = ppval(my_pp_CRx, rawUR);

figure()
subplot(3,1,1)
plot(rawUC, xEast_C_GE_raw(:))
hold on
plot(rawUC, Cx_raw_smoothed(:))
xlabel('rawUC, m')
ylabel('xEast C GE raw, m')
legend('raw Cx', 'Cx smoothed')
title('xEast C GE raw')
subplot(3,1,2)
plot(rawUL, xEast_CL_GE_raw(:))
hold on
plot(rawUL, CLx_raw_smoothed(:))
xlabel('rawUL, m')
ylabel('xEast L GE raw, m')
legend('raw CLx', 'CLx smoothed')
title('xEast L GE raw')
subplot(3,1,3)
plot(rawUR, xEast_CR_GE_raw(:))
hold on
plot(rawUR, CRx_raw_smoothed(:))
xlabel('rawUR, m')
ylabel('xEast R GE raw, m')
legend('raw CRx', 'CRz smoothed')
title('xEast R GE raw')

%% Smooth raw y data
%
my_pp_Cy = myspcsp(rawUC.', yNorth_C_GE_raw.', p_y);
Cy_raw_smoothed = ppval(my_pp_Cy, rawUC);

my_pp_CLy = myspcsp(rawUL.', yNorth_CL_GE_raw.', p_y);
CLy_raw_smoothed = ppval(my_pp_CLy, rawUL);

my_pp_CRy = myspcsp(rawUR.', yNorth_CR_GE_raw.', p_y);
CRy_raw_smoothed = ppval(my_pp_CRy, rawUR);

figure()
subplot(3,1,1)
plot(rawUC, yNorth_C_GE_raw(:))
hold on
plot(rawUC, Cy_raw_smoothed(:))
xlabel('rawUC, m')
ylabel('yNorth C GE raw, m')
legend('raw Cy', 'Cy smoothed')
title('yNorth C GE raw')
subplot(3,1,2)
plot(rawUL, yNorth_CL_GE_raw(:))
hold on
plot(rawUL, CLy_raw_smoothed(:))
xlabel('rawUL, m')
ylabel('yNorth L GE raw, m')
legend('raw CLy', 'CLy smoothed')
title('yNorth L GE raw')
subplot(3,1,3)
plot(rawUR, yNorth_CR_GE_raw(:))
hold on
plot(rawUR, CRy_raw_smoothed(:))
xlabel('rawUR, m')
ylabel('yNorth R GE raw, m')
legend('raw CRy', 'CRy smoothed')
title('yNorth R GE raw')

%% Replace the raw Google Earth data with smoothed data
%
rawC(:, 3) = Cz_raw_smoothed;
rawL(:, 3) = CLz_raw_smoothed;
rawR(:, 3) = CRz_raw_smoothed;

rawC(:, 1) = Cx_raw_smoothed  + shift_X;
rawL(:, 1) = CLx_raw_smoothed + shift_X;
rawR(:, 1) = CRx_raw_smoothed + shift_X;

rawC(:, 2) = Cy_raw_smoothed;
rawL(:, 2) = CLy_raw_smoothed;
rawR(:, 2) = CRy_raw_smoothed;

% Piecewise-length parameterization
rawdC = diff(rawC,1,1);
rawUC = [0; cumsum( sqrt(sum(rawdC.^2,2)) )];

rawdL = diff(rawL,1,1);
rawUL = [0; cumsum( sqrt(sum(rawdL.^2,2)) )];

rawdR = diff(rawR,1,1);
rawUR = [0; cumsum( sqrt(sum(rawdR.^2,2)) )];

%% Smooth the smoothed raw coordinates
%
my_pp_C = myspcsp(rawUC.', rawC.', p);
my_pp_L = myspcsp(rawUL.', rawL.', p);
my_pp_R = myspcsp(rawUR.', rawR.', p);

%% Totoal length computation
%
[C_length_piecewise, C_length_cumtrapz, uC_of_S, h_vec, S_all_points] = Compute_length_of_spline(my_pp_C, N_integration_grid, N_nodes, c1_RadauIIA, c2_RadauIIA, sa);
[L_length_piecewise, L_length_cumtrapz, uL_of_S, ~, ~]                = Compute_length_of_spline(my_pp_L, N_integration_grid, N_nodes, c1_RadauIIA, c2_RadauIIA, sa);
[R_length_piecewise, R_length_cumtrapz, uR_of_S, ~, ~]                = Compute_length_of_spline(my_pp_R, N_integration_grid, N_nodes, c1_RadauIIA, c2_RadauIIA, sa);

%% Compute C_true and unit tangent vector of C_true
%
[C_true, T_C] = Compute_C_true_and_T_C(my_pp_C, uC_of_S);

disp('Compute_C_true_and_T_C finsished.')

C_true_diff_SF = C_true(end,:) - C_true(1,:)
T_C_diff_SF = T_C(end,:) - T_C(1,:)

% debug = 1;

%% Compute intersection points: between a plane perpendicular to T_C and L spline (and R spline)
%
[P_L, P_R, U_L, U_R, FL_at_UL, FR_at_UR] = Compute_intersection_points(C_true, T_C, my_pp_L, my_pp_R, uL_of_S, uR_of_S);

P_L_SF_diff = P_L(1,:)-P_L(end,:)
P_R_SF_diff = P_R(1,:)-P_R(end,:)

figure()
subplot(3,1,1)
plot(P_L(:, 1), '-o')
title('PL X')
subplot(3,1,2)
plot(P_L(:, 2), '-o')
title('PL Y')
subplot(3,1,3)
plot(P_L(:, 3), '-o')
title('PL Z')

figure()
subplot(3,1,1)
plot(P_R(:, 1), '-o')
title('PR X')
subplot(3,1,2)
plot(P_R(:, 2), '-o')
title('PR Y')
subplot(3,1,3)
plot(P_R(:, 3), '-o')
title('PR Z')

figure()
subplot(2,1,1)
plot(U_L, FL_at_UL)
xlabel('UL')
ylabel('FL at UL')
title('FL at UL')
subplot(2,1,2)
plot(U_R, FR_at_UR)
xlabel('UR')
ylabel('FR at UR')
title('FR at UR')

disp('Compute_intersection_points finished.')

% debug = 1;

%% Compute C_ture_project
%
% C_true_project = P_L + (dot((C_true-P_L), (P_R-P_L))/dot((P_R-P_L), (P_R-P_L)))*(P_R-P_L);

C_true_project = nan(length(P_L), 3);
for index = 1:length(P_L)
    
C_true_project(index, :) = P_L(index, :) + (dot((C_true(index, :)-P_L(index, :)), (P_R(index, :)-P_L(index, :))) / norm(P_R(index, :)-P_L(index, :))) * ((P_R(index, :)-P_L(index, :)) / norm((P_R(index, :)-P_L(index, :))));

end

C_true_project_diff_SF = C_true_project(1,:) - C_true_project(end,:)

figure;
% Plot Left intersection points
plot(P_L(:,1), P_L(:,2), 'r-o', 'LineWidth', 1.5);
hold on;
% Plot Right intersection points
plot(P_R(:,1), P_R(:,2), 'b-o', 'LineWidth', 1.5);
% Draw only sampled C_true_project
plot( ...
    C_true_project(:,1), ...
    C_true_project(:,2), ...
    'k-o', ...              % black circles connected by lines
    'LineWidth', 0.5, ...
    'MarkerSize', 3, ...
    'MarkerFaceColor','k' ...
);
% plot(C_true(:,1),C_true(:,2), 'g-o', 'LineWidth', 0.5, 'MarkerSize', 3, 'MarkerFaceColor','g')
% legend('P\_L', 'P\_R', 'C\_true\_project', 'C\_true')
legend('P\_L', 'P\_R', 'C\_true\_project')
xlabel('X'); ylabel('Y');
axis equal;
grid on;
hold off;

disp('C_ture_project finished.')

% debug = 1;

%% Compute unit N_C vector
%
N_C = (P_L - P_R) ./ vecnorm((P_L - P_R),2,2);

%% Compute unit M_C vector
%
M_C = cross(T_C, N_C, 2);
M_C = M_C ./ vecnorm(M_C,2,2);

%% Compute Euler angles
%-- direct “inverse‐Euler” extraction from [t,n,m]
Theta = unwrap( atan2( T_C(:,2), T_C(:,1) ) );    % θ = atan2(t_y, t_x)
Mu    = -asin(    T_C(:,3)       );               % μ = -asin(t_z)
Phi   = atan2( N_C(:,3), M_C(:,3)    );           % φ = atan2(n_z, m_z)

Theta_diff_SF = Theta(1,:) - Theta(end,:) - (2*pi)  % here since Theta has been unwrapped, the initial and end values may have difference of 2pi*k
Mu_diff_SF    = Mu(1,:)    - Mu(end,:)
Phi_diff_SF   = Phi(1,:)   - Phi(end,:)

cos_Theta = cos(Theta);
sin_Theta = sin(Theta);
cos_Mu = cos(Mu);
sin_Mu = sin(Mu);
cos_Phi = cos(Phi);
sin_Phi = sin(Phi);

figure()
subplot(3,1,1)
plot(S_all_points, Theta, '-o')
title('Theta')
xlabel('SC')
subplot(3,1,2)
plot(S_all_points, Mu, '-o')
title('Mu')
xlabel('SC')
subplot(3,1,3)
plot(S_all_points, Phi, '-o')
title('Phi')
xlabel('SC')

% fname1 = fullfile('Theta.txt');
%     fid1 = fopen(fname1, 'w');
%     for i = 1:length(Theta)
%         fprintf(fid1, '%0.17g\n', Theta(i));
%     end
% fclose(fid1);
%     
% fname2 = fullfile('Mu.txt');
%     fid2 = fopen(fname2, 'w');
%     for i = 1:length(Mu)
%         fprintf(fid2, '%0.17g\n', Mu(i));
%     end
% fclose(fid2);
% 
% fname3 = fullfile('Phi.txt');
%     fid3 = fopen(fname3, 'w');
%     for i = 1:length(Phi)
%         fprintf(fid3, '%0.17g\n', Phi(i));
%     end
% fclose(fid3);
% 
% fname4 = fullfile('S_all_points.txt');
%     fid4 = fopen(fname4, 'w');
%     for i = 1:length(S_all_points)
%         fprintf(fid4, '%0.17g\n', S_all_points(i));
%     end
% fclose(fid4);

%% Estimate dEuler_angles/ds and d^2Euler_angles/ds^2
%
[dTheta, dMu, dPhi, ddTheta, ddMu, ddPhi] = Estimate_dEulerds_ddEulerds2(S_all_points, Theta, Mu, Phi);

dTheta_diff_SF = dTheta(1) - dTheta(end)
dMu_diff_SF    = dMu(1)    - dMu(end)
dPhi_diff_SF   = dPhi(1)   - dPhi(end)

figure()
subplot(3,1,1)
plot(S_all_points, dTheta, '-o')
xlabel('SC')
title('dTheta')
subplot(3,1,2)
plot(S_all_points, dMu, '-o')
xlabel('SC')
title('dMu')
subplot(3,1,3)
plot(S_all_points, dPhi, '-o')
xlabel('SC')
title('dPhi')

figure()
subplot(3,1,1)
plot(S_all_points, ddTheta, '-o')
xlabel('SC')
title('ddTheta')
subplot(3,1,2)
plot(S_all_points, ddMu, '-o')
xlabel('SC')
title('ddMu')
subplot(3,1,3)
plot(S_all_points, ddPhi, '-o')
xlabel('SC')
title('ddPhi')

disp('Estimate_dEulerds_ddEulerds2 finished.')

% debug = 1;

%% Compute nl and nr
%
nl = vecnorm((P_L - C_true_project),2,2);
nr = (-1) * vecnorm((P_R - C_true_project),2,2);

nl_diff_SF = nl(1) - nl(end)
nr_diff_SF = nr(1) - nr(end)

figure()
plot(S_all_points, nl, '-r')
hold on
plot(S_all_points, nr, '-b')
xlabel('SC')
ylabel('nl, nr')
legend('nl', 'nr')

%% Estimate dnl/ds and dnr/ds
%
[dnl, dnr] = Estimate_dnlds_dnrds(S_all_points, nl, nr);

figure()
subplot(2,1,1)
plot(S_all_points, dnl, '-r')
xlabel('SC')
title('dnl')
subplot(2,1,2)
plot(S_all_points, dnr, '-b')
xlabel('SC')
title('dnr')

disp('Estimate_dnlds_dnrds finished.')

% debug = 1;

%%
%
% All_X_table = [C_true_project(:,1),...
%                C_true_project(:,2),...
%                C_true_project(:,3),...
%                cos_Theta,...
%                sin_Theta,...
%                cos_Mu,...
%                sin_Mu,...
%                cos_Phi,...
%                sin_Phi,...
%                dTheta,...
%                dMu,...
%                dPhi,...
%                nl,...
%                nr];
All_X_table = [C_true(:,1),...
               C_true(:,2),...
               C_true(:,3),...
               cos_Theta,...
               sin_Theta,...
               cos_Mu,...
               sin_Mu,...
               cos_Phi,...
               sin_Phi,...
               dTheta,...
               dMu,...
               dPhi,...
               nl,...
               nr];
           
[D_bar_kron_I, D_bar] = Compute_Dbar_kron_I(N_nodes, sa, nx, A_matrix_RadauIIA, h_vec);

ddTheta_from_DbartimesdTheta = D_bar * dTheta;
ddMu_from_DbartimesdMu       = D_bar * dMu;
ddPhi_from_DbartimesdPhi     = D_bar * dPhi;

dnl_from_Dbartimesnl = D_bar * nl;
dnr_from_Dbartimesnr = D_bar * nr;

ddTheta_from_DbartimesdTheta_with_t1 = [ddTheta(1); ddTheta_from_DbartimesdTheta];
ddMu_from_DbartimesdMu_with_t1       = [ddMu(1);    ddMu_from_DbartimesdMu];
ddPhi_from_DbartimesdPhi_with_t1     = [ddPhi(1); ddPhi_from_DbartimesdPhi];

dnl_from_Dbartimesnl_with_t1 = [dnl(1); dnl_from_Dbartimesnl];
dnr_from_Dbartimesnr_with_t1 = [dnr(1); dnr_from_Dbartimesnr];

All_Fs_table_new = ODEs_rhs(cos_Theta, sin_Theta, cos_Mu, sin_Mu, cos_Phi, sin_Phi,...
                                            dTheta, dMu, dPhi,...
                                            ddTheta_from_DbartimesdTheta_with_t1, ddMu_from_DbartimesdMu_with_t1, ddPhi_from_DbartimesdPhi_with_t1,...
                                            dnl_from_Dbartimesnl_with_t1, dnr_from_Dbartimesnr_with_t1);

discrete_ODE_residuals_new = discrete_ODE_residuals(All_Fs_table_new, All_X_table, D_bar_kron_I);

discrete_ODE_residuals_reshaped_new = reshape(discrete_ODE_residuals_new, nx, []);

for index_state = 1:nx
    
    figure()
    plot(discrete_ODE_residuals_reshaped_new(index_state, :))
    xlabel('Index of collocation points')
    ylabel(['Discrete ODE equation ', num2str(index_state), ' residuals'])
    grid on
    
end

disp('Verify new ODE residuals finished.')

% debug = 1;

%% A quick test of system dynamics
%
All_Fs_table = ODEs_rhs(cos_Theta, sin_Theta, cos_Mu, sin_Mu, cos_Phi, sin_Phi,...
                                            dTheta, dMu, dPhi, ddTheta, ddMu, ddPhi, dnl, dnr);
F1  = All_Fs_table(:, 1);  
F2  = All_Fs_table(:, 2);
F3  = All_Fs_table(:, 3);
F4  = All_Fs_table(:, 4);
F5  = All_Fs_table(:, 5);
F6  = All_Fs_table(:, 6);
F7  = All_Fs_table(:, 7);
F8  = All_Fs_table(:, 8);
F9  = All_Fs_table(:, 9);
F10 = All_Fs_table(:, 10);
F11 = All_Fs_table(:, 11);
F12 = All_Fs_table(:, 12);
F13 = All_Fs_table(:, 13);
F14 = All_Fs_table(:, 14);

X_estimate = cumtrapz(S_all_points, F1) + C_true(1, 1);
Y_estimate = cumtrapz(S_all_points, F2) + C_true(1, 2);
Z_estimate = cumtrapz(S_all_points, F3) + C_true(1, 3);

cos_Theta_estimate = cumtrapz(S_all_points, F4) + cos_Theta(1);
sin_Theta_estimate = cumtrapz(S_all_points, F5) + sin_Theta(1);

cos_Mu_estimate = cumtrapz(S_all_points, F6) + cos_Mu(1);
sin_Mu_estimate = cumtrapz(S_all_points, F7) + sin_Mu(1);

cos_Phi_estimate = cumtrapz(S_all_points, F8) + cos_Phi(1);
sin_Phi_estimate = cumtrapz(S_all_points, F9) + sin_Phi(1);

d_Theta_estimate = cumtrapz(S_all_points, F10) + dTheta(1);
d_Mu_estimate    = cumtrapz(S_all_points, F11) + dMu(1);
d_Phi_estimate   = cumtrapz(S_all_points, F12) + dPhi(1);

nl_estimate = cumtrapz(S_all_points, F13) + nl(1);
nr_estimate = cumtrapz(S_all_points, F14) + nr(1);

figure()
subplot(3,1,1)
plot(S_all_points, C_true(:, 1)-X_estimate)
xlabel('SC')
title('Difference between X true and X estimate')
subplot(3,1,2)
plot(S_all_points, C_true(:, 2)-Y_estimate)
xlabel('SC')
title('Difference between Y true and Y estimate')
subplot(3,1,3)
plot(S_all_points, C_true(:, 3)-Z_estimate)
xlabel('SC')
title('Difference between Z true and Z estimate')

figure()
subplot(3,1,1)
plot(S_all_points, cos_Theta-cos_Theta_estimate)
hold on
plot(S_all_points, sin_Theta-sin_Theta_estimate)
legend('cos', 'sin')
title('Difference between cos Theta and its estimate and sin Theta and its estimate')
xlabel('SC')
subplot(3,1,2)
plot(S_all_points, cos_Mu-cos_Mu_estimate)
hold on
plot(S_all_points, sin_Mu-sin_Mu_estimate)
legend('cos', 'sin')
title('Difference between cos Mu and its estimate and sin Mu and its estimate')
xlabel('SC')
subplot(3,1,3)
plot(S_all_points, cos_Phi-cos_Phi_estimate)
hold on
plot(S_all_points, sin_Phi-sin_Phi_estimate)
legend('cos', 'sin')
title('Difference between cos Phi and its estimate and sin Phi and its estimate')
xlabel('SC')

figure()
subplot(3,1,1)
plot(S_all_points, dTheta-d_Theta_estimate)
xlabel('SC')
title('Difference between dTheta and its estimate')
subplot(3,1,2)
plot(S_all_points, dMu-d_Mu_estimate)
xlabel('SC')
title('Difference between dMu and its estimate')
subplot(3,1,3)
plot(S_all_points, dPhi-d_Phi_estimate)
xlabel('SC')
title('Difference between dPhi and its estimate')

figure()
subplot(2,1,1)
plot(S_all_points, nl-nl_estimate)
xlabel('SC')
title('Difference between nl and its estimate')
subplot(2,1,2)
plot(S_all_points, nr-nr_estimate)
xlabel('SC')
title('Difference between nr and its estimate')

disp('A quick test of system dynamics finished.')

%% Verify tracking errors
%
All_X_table_for_tracking =...
              [C_true_project(:,1),...
               C_true_project(:,2),...
               C_true_project(:,3),...
               cos_Theta,...
               sin_Theta,...
               cos_Mu,...
               sin_Mu,...
               cos_Phi,...
               sin_Phi,...
               dTheta,...
               dMu,...
               dPhi,...
               nl,...
               nr];

x_C_m = All_X_table_for_tracking(:, 1);
y_C_m = All_X_table_for_tracking(:, 2);
z_C_m = All_X_table_for_tracking(:, 3);

x_CL_m = P_L(:, 1);
y_CL_m = P_L(:, 2);
z_CL_m = P_L(:, 3);

x_CR_m = P_R(:, 1);
y_CR_m = P_R(:, 2);
z_CR_m = P_R(:, 3);

[tracking_error_value_1, tracking_error_value_2, tracking_error_value_3, x_CL, y_CL, z_CL, x_n_hat, y_n_hat, z_n_hat] = tracking_error(All_X_table_for_tracking, x_C_m, y_C_m, z_C_m, x_CL_m, y_CL_m, z_CL_m, x_CR_m, y_CR_m, z_CR_m);

figure()
subplot(3,1,1)
plot(tracking_error_value_1)
title('C tracking error')
subplot(3,1,2)
plot(tracking_error_value_2)
title('L tracking error')
subplot(3,1,3)
plot(tracking_error_value_3)
title('R tracking error')

figure()
subplot(3,1,1)
plot(P_L(:,1) - x_CL)
title('Difference between PLx and xCL')
subplot(3,1,2)
plot(P_L(:,2) - y_CL)
title('Difference between PLy and yCL')
subplot(3,1,3)
plot(P_L(:,3) - z_CL)
title('Difference between PLz and zCL')

figure()
subplot(3,1,1)
plot(N_C(:,1) - x_n_hat)
title('Difference between NCx and xnhat')
subplot(3,1,2)
plot(N_C(:,2) - y_n_hat)
title('Difference between NCy and ynhat')
subplot(3,1,3)
plot(N_C(:,3) - z_n_hat)
title('Difference between NCz and znhat')

disp('Verify tracking errors finished.')

%% Angular velocity (independent variable is distance) of the ribbon model
%
OmegaB = NaN(length(S_all_points), 3);

% Jacobian K(θ,μ,φ) for z–y–x Euler sequence (eq 4.25):
for index_points = 1:length(S_all_points)
    
    Mu_i  = Mu(index_points);
    Phi_i = Phi(index_points);
    
    K_i = [
                -sin(Mu_i),         0,   1;
        cos(Mu_i)*sin(Phi_i),  cos(Phi_i), 0;
        cos(Mu_i)*cos(Phi_i), -sin(Phi_i), 0
        ];
    
    dTheta_i = dTheta(index_points);
    dMu_i    = dMu(index_points);
    dPhi_i   = dPhi(index_points);
    
    OmegaB_i_temp = K_i * [dTheta_i, dMu_i, dPhi_i].'; % (eq 4.22)
    
    OmegaB_i = OmegaB_i_temp.';
    
    OmegaB(index_points, :) = OmegaB_i;
    
end

figure()
subplot(3,1,1)
plot(S_all_points, OmegaB(:, 1))
grid on
xlabel('SC, m')
ylabel('Omega B x')
subplot(3,1,2)
plot(S_all_points, OmegaB(:, 2))
grid on
xlabel('SC, m')
ylabel('Omega B y')
subplot(3,1,3)
plot(S_all_points, OmegaB(:, 3))
grid on
xlabel('SC, m')
ylabel('Omega B z')

%% Plot 3D track
%
figure()
plot3(P_L(:,1), P_L(:,2), P_L(:,3)*0, '-k')
hold on
plot3(P_R(:,1), P_R(:,2), P_R(:,3)*0, '-k')
hold on
plot3(P_L(:,1), P_L(:,2), P_L(:,3)-min(P_L(:,3)), '-r')
hold on
plot3(P_R(:,1), P_R(:,2), P_R(:,3)-min(P_R(:,3)), '-b')
title('Processed 3D track with X shifted with lowest elevation 0m')
axis equal; grid on;

%% Save
%
save('S_all_points', 'S_all_points')
save('h_vec', 'h_vec')
save('P_L', 'P_L')
save('P_R', 'P_R')
save('C_true', 'C_true')
save('C_true_project', 'C_true_project')
save('nl', 'nl')
save('nr', 'nr')
save('Theta', 'Theta')
save('Mu', 'Mu')
save('Phi', 'Phi')
save('T_C', 'T_C')
save('N_C', 'N_C')
save('M_C', 'M_C')
save('OmegaB', 'OmegaB')

%%
%
function pp = myspcsp(x,Y,p,varargin)
% SMOOTHING CUBIC SPLINE WITH PERIODIC CONDITIONS
%
% SYNTAX:
%   pp = spcsp(x,Y,p);
%   pp = spcsp(x,Y,p,w);
%
% spcsp(x,Y,p,w) returns the ppform of a cubic smoothing spline for the
% given data x,Y, with periodic conditions. The smoothing spline 
% approximates, at the data site x(j), the given data value Y(:,j), where 
% j=1:length(x). The data values may be scalars, vectors, matrices, 
% or even ND-arrays. If weigth vector w is not supplied, default weigths to
% 1 are used.
%
% INPUTS:
% x = row vector of breaks (parametric variable of the d-dimensional curve)
% Y = d-by-n matrix of the values to approximate
%       d: dimension of the space
%       n: number of points
% p = smoothing parameter
% w = row vector of weights
%
% The smoothing spline f minimizes
%
%    p * sum_j w(j)^(-2) |Y(:,j) - f(x(j))|^2  +  (1-p) * integral |D^2 f|^2
% 
% where the sum is over j=1:length(X).
%
% OUTPUTS:
% pp = ppform of the approximating spline.

% ################# ALGORITHM #############################################
    
x=x(:)';
n = numel(x); % number of breaks
[dim,n2] = size(Y); % space dimension and number of points
Y = Y(:,1:end-1);

%check inputs
if ~isempty(varargin)
    w=(varargin{1}(:))';
else
    w=ones(size(x));
end
if n~=n2
    error('size of x and Y do not match');
end
if n~=numel(w)
    error('numel(x)~=numel(w), check dimensions');
end

% h array of sub-intervals between breaks
h = x(2:end)-x(1:end-1); % n-1 array

% S matrix
aux = [ [h(1:end-1), 0] ; 2*(h+[h(end), h(1:end-1)]) ; [0, h(1:end-1)] ]';
S = spdiags( aux , [-1 0 1] , n-1 , n-1 );
S(1,n-1) = h(end);
S(n-1,1) = h(end);

% V matrix
aux = [ [1./h(1:end-1), 0] ; -1./h-1./[h(end), h(1:end-1)] ; [0, 1./h(1:end-1)] ]';
V = spdiags( aux , [-1 0 1] , n-1 , n-1 );
V(1,n-1) = 1/h(end);
V(n-1,1) = 1/h(end);

% generate weight matrix from w, matrix W:=W^-2
W = spdiags( w(1:end-1)'.^(-2) , 0 , n-1 , n-1 );

% constructing the quadratic form: F(a) = 1/2 a'Ua - v'a
AUX = S\V; %AUX = S^-1 V
U = 2*p*W + 12*(1-p)*V'*AUX;
v = 2*p*W*Y';
h = h';

% minimizing and computing coefficients
a = U\v;
c = 3*AUX*a;
d = ([c(2:end,:);c(1,:)]-c)./(3*h);
b = ([a(2:end,:);a(1,:)]-a)./h - c.*h - d.*(h.^2);
coefs = cat(3,d',c',b',a');

% build spline in pp form: breaks=x
pp = mkpp(x,coefs,dim);

end

%%
%
function [length_integral, length_cumtrapz, u_of_S, h_vec, S_all_points] = Compute_length_of_spline(my_pp, N_uFine, N_Sgrid, c1_RadauIIA, c2_RadauIIA, sa)

% Compute true length of the curve
% 1) Differentiate the pp-splines
dpp = fnder(my_pp, 1);

ppd = fnder(my_pp,1);
integrand = @(s) sqrt( sum( ppval(ppd,s).^2 , 1 ) );
breaks = my_pp.breaks;       % your x-knots (rawUC)
length_integral = 0;
for k = 1 : numel(breaks)-1
  length_integral = length_integral + integral( integrand, breaks(k), breaks(k+1), ...
                    'RelTol',1e-15, 'AbsTol',1e-15 );
end
fprintf('Spline length ≈ %.12g m\n', length_integral);


% 3) NUMERICAL length via cumtrapz on a fine uniform grid
uMin  = breaks(1);
uMax  = breaks(end);
uFine = linspace(uMin, uMax, N_uFine);

dxyz    = ppval(dpp, uFine);
dxyz = dxyz.';
speed = sqrt(dxyz(:,1).^2 + dxyz(:,2).^2 + dxyz(:,3).^2);

length_cumtrapz = cumtrapz(uFine.', speed);

% Map a uniform grid of the true curve length back to the u domain
% True length of the curve
Length_true  = length_cumtrapz(end);

% Define a uniform arc-length grid
S_true_uniform_grid  = linspace(0, Length_true, N_Sgrid);

h_vec = diff(S_true_uniform_grid);

S_all_points = augment_S_nodes(S_true_uniform_grid, c1_RadauIIA, c2_RadauIIA, sa);

u_of_S = interp1(length_cumtrapz, uFine, S_all_points.', 'pchip');

% debug = 1;

end

%%
%
function [A_matrix, b_vec, c_vec] = RadauIIA_tabular()

% the following values come from paper "Stiff differential equations solved
% by Radau methods, Ernst Hairer , Gerhard Wanner

a11 = (88-7*(6^0.5)) / 360;

a12 = (296-169*(6^0.5)) / 1800;

a13 = (-2+3*(6^0.5)) / 225;

a21 = (296+169*(6^0.5)) / 1800;

a22 = (88+7*(6^0.5)) / 360;

a23 = (-2-3*(6^0.5)) / 225;

a31 = (16-(6^0.5)) / 36;

a32 = (16+(6^0.5)) / 36;

a33 = 1/9;

% A matrix
A_matrix = [a11, a12, a13;
            a21, a22, a23;
            a31, a32, a33];

% b_vec
b_vec = [a31, a32, a33];

% c_vec
c1 = (4-6^0.5) / 10;
c2 = (4+6^0.5) / 10;
c3 = 1;
c_vec = [c1; c2; c3];

end

%%
%
function S_all_points = augment_S_nodes(S_nodes, c1, c2, sa)
    % S_nodes:    (N×1) sorted vector
    % c1, c2:     c1 and c2 in Radau IIA method
    % sa:         number of collocation points per interval
    
    N = numel(S_nodes);
    if N<2
        S_all_points = S_nodes;
        return;
    end

    h = diff(S_nodes);                 % (N−1)×1 intervals
    L = sa*(N-1) + 1;                   % total length = 3*(N−1)+1
    S_all_points = zeros(L,1);         % pre‐allocate

    idx = 1;
    for i = 1:N-1
        S0 = S_nodes(i);
        S_all_points(idx)   = S0;                  % original node
        S_all_points(idx+1) = S0 + c1 * h(i);      % first extra
        S_all_points(idx+2) = S0 + c2 * h(i);      % second extra
        idx = idx + 3;
    end

    S_all_points(end) = S_nodes(end);   % append the very last original node
end

%%
%
function [C_true, T_C] = Compute_C_true_and_T_C(my_pp, uC_of_S)

C_true = ppval(my_pp, uC_of_S)';    % N×3

% 1) form the derivative‐splines (order = 1)
dspC = fnder(my_pp, 1);

% 2) evaluate them on uC_of_S
dC_du = ppval(dspC, uC_of_S);

% 3) dC/du
dC_du = dC_du';    % size N×3

% ds/du
norm_dC_du = vecnorm(dC_du,2,2);

% unit tangent vector of spline C, dC/ds = (dC/du) / (ds/du)
T_C = dC_du ./ norm_dC_du;

% debug = 1;

end

%%
% %
% function [P_L, P_R, U_L, U_R, FL_at_UL, FR_at_UR] = Compute_intersection_points(C_true, T_C, my_pp_L, my_pp_R, uL_of_S, uR_of_S)
% 
% N_all_points = length(C_true(:, 1));
% 
% % pre‐allocate
% P_L = zeros(N_all_points,3);
% P_R = zeros(N_all_points,3);
% U_L = zeros(N_all_points,1);
% U_R = zeros(N_all_points,1);
% 
% FL_at_UL = zeros(N_all_points,1);
% FR_at_UR = zeros(N_all_points,1);
% 
% U_L_max = uL_of_S(end);
% U_R_max = uR_of_S(end);
% 
% for i = 1:N_all_points
%   Ci = C_true(i,:).';    Ti = T_C(i,:).';
% 
%   % define the function
% %   fL = @(u) ([fnval(spLx,u); fnval(spLy,u); fnval(spLz,u)] - Ci).' * Ti;
% %   fR = @(u) ([fnval(spRx,u); fnval(spRy,u); fnval(spRz,u)] - Ci).' * Ti;
%   fL = @(u) (ppval(my_pp_L,u) - Ci).' * Ti;
%   fR = @(u) (ppval(my_pp_R,u) - Ci).' * Ti;
% 
%   % initial guess
%   u0_L = uL_of_S(i);
%   u0_R = uR_of_S(i);
% 
%   % left: try single‐guess fzero first
%   try
%     uLroot = fzero(fL, u0_L);
%   catch
%     % if that fails, find a valid bracket by sampling
%     wL   = 0.1*U_L_max;
%     samps_L = linspace(max(0,u0_L-wL), min(U_L_max,u0_L+wL), 21);
%     vals_L  = arrayfun(fL, samps_L);
%     idx_L   = find(vals_L(1:end-1).*vals_L(2:end) < 0, 1);
%     if isempty(idx_L)
%       % last‐ditch: sample the whole domain
%       samps_L = linspace(0, U_L_max, 200);
%       vals_L  = arrayfun(fL, samps_L);
%       idx_L   = find(vals_L(1:end-1).*vals_L(2:end) < 0, 1);
%       if isempty(idx_L)
%         error('No sign change for left‐rail at station %d.', i);
%       end
%     end
%     a_L = samps_L(idx_L);    b_L = samps_L(idx_L+1);
%     uLroot = fzero(fL, [a_L, b_L]);
%   end
%   
%   % right: try single‐guess fzero first
%   try
%     uRroot = fzero(fR, u0_R);
%   catch
%     % if that fails, find a valid bracket by sampling
%     wR   = 0.1*U_R_max;
%     samps_R = linspace(max(0,u0_R-wR), min(U_R_max,u0_R+wR), 21);
%     vals_R  = arrayfun(fR, samps_R);
%     idx_R   = find(vals_R(1:end-1).*vals_R(2:end) < 0, 1);
%     if isempty(idx_R)
%       % last‐ditch: sample the whole domain
%       samps_R = linspace(0, U_R_max, 200);
%       vals_R  = arrayfun(fR, samps_R);
%       idx_R   = find(vals_R(1:end-1).*vals_R(2:end) < 0, 1);
%       if isempty(idx_R)
%         error('No sign change for right‐rail at station %d.', i);
%       end
%     end
%     a_R = samps_R(idx_R);    b_R = samps_R(idx_R+1);
%     uRroot = fzero(fR, [a_R, b_R]);
%   end
% 
%   % evaluate and store
% %   P_L(i,:) = [ fnval(spLx,uLroot), fnval(spLy,uLroot), fnval(spLz,uLroot) ];
%   P_L(i,:) = ppval(my_pp_L,uLroot);
%   U_L(i)  = uLroot;
%   
%   FL_at_UL(i) = fL(uLroot); 
%   
% %   P_R(i,:) = [ fnval(spRx,uRroot), fnval(spRy,uRroot), fnval(spRz,uRroot) ];
%   P_R(i,:) = ppval(my_pp_R,uRroot);
%   U_R(i)  = uRroot;
%   
%   FR_at_UR(i) = fR(uRroot); 
%   
% %   debug = 1;
%   
% end
% 
% end


% function [P_L, P_R, U_L, U_R, FL_at_UL, FR_at_UR] = ...
%     Compute_intersection_points(C_true, T_C, my_pp_L, my_pp_R, uL_of_S, uR_of_S)
% 
% N = size(C_true,1);
% P_L = zeros(N,3);  P_R = zeros(N,3);
% U_L = zeros(N,1);  U_R = zeros(N,1);
% FL_at_UL = zeros(N,1); FR_at_UR = zeros(N,1);
% 
% U_L_max = my_pp_L.breaks(end);
% U_R_max = my_pp_R.breaks(end);
% 
% % --- your track specifics ---
% ds_nom      = 1.0;     % m (Δs between stations, for reference)
% W_half      = 6.0;     % m (half track width near hairpin)
% search_s    = 2.0*W_half;  % m (± arclength window in hairpins)
% alpha       = 0.5;     % weight for continuity (m/m)
% dist_cap    = 1.5*W_half;  % m (reject far intersections)
% angle_thr   = 25*pi/180;   % rad, hairpin detector
% grid_n      = 41;      % samples to bracket zeros
% % --------------------------------
% 
% % derivatives for speed (m per unit-u)
% dppL = fnder(my_pp_L,1);
% dppR = fnder(my_pp_R,1);
% 
% for i = 1:N
%     Ci = C_true(i,:).';  Ti = T_C(i,:).';
% 
%     fL = @(u) (ppval(my_pp_L,u) - Ci).' * Ti;
%     fR = @(u) (ppval(my_pp_R,u) - Ci).' * Ti;
% 
%     % clamp seeds
%     u0L = min(max(uL_of_S(i),0), U_L_max);
%     u0R = min(max(uR_of_S(i),0), U_R_max);
% 
%     % hairpin detector (tangent change)
%     if i>1 && i<N
%         ang = acos(max(-1,min(1, dot(T_C(i-1,:), T_C(i+1,:)) )));
%         isHairpin = ang > angle_thr;
%     else
%         isHairpin = false;
%     end
% 
%     % local u-windows based on arclength window
%     vL = ppval(dppL,u0L); speedL = max(norm(vL), eps);      % m per u
%     vR = ppval(dppR,u0R); speedR = max(norm(vR), eps);
%     duWinL = search_s / speedL; duWinR = search_s / speedR;
%     uLoL = max(0, u0L - duWinL); uHiL = min(U_L_max, u0L + duWinL);
%     uLoR = max(0, u0R - duWinR); uHiR = min(U_R_max, u0R + duWinR);
% 
%     % ===== LEFT =====
%     uLroot = NaN; pL = [NaN;NaN;NaN];
%     if ~isHairpin
%         try
%             uLroot = fzero(fL, u0L); pL = ppval(my_pp_L,uLroot);
%         catch, end
%     end
%     if isnan(uLroot) || norm(pL-Ci)>dist_cap
%         [uLroot,pL] = bestRootWindow(my_pp_L,dppL,fL,Ci,uLoL,uHiL,u0L,alpha,grid_n,dist_cap);
%     end
%     if isnan(uLroot)   % last ditch: closest point in same window
%         [uLroot,pL] = closestPointOnCurve(my_pp_L,Ci,uLoL,uHiL);
%     end
%     P_L(i,:) = pL.'; U_L(i)=uLroot; FL_at_UL(i)=fL(uLroot);
% 
%     % ===== RIGHT =====
%     uRroot = NaN; pR = [NaN;NaN;NaN];
%     if ~isHairpin
%         try
%             uRroot = fzero(fR, u0R); pR = ppval(my_pp_R,uRroot);
%         catch, end
%     end
%     if isnan(uRroot) || norm(pR-Ci)>dist_cap
%         [uRroot,pR] = bestRootWindow(my_pp_R,dppR,fR,Ci,uLoR,uHiR,u0R,alpha,grid_n,dist_cap);
%     end
%     if isnan(uRroot)
%         [uRroot,pR] = closestPointOnCurve(my_pp_R,Ci,uLoR,uHiR);
%     end
%     P_R(i,:) = pR.'; U_R(i)=uRroot; FR_at_UR(i)=fR(uRroot);
% end
% end
% 
% % ---- helpers ----
% function [uBest,pBest] = bestRootWindow(pp,dpp,f,Ci,uLo,uHi,uSeed,alpha,grid_n,dist_cap)
%     grid = linspace(uLo,uHi,max(3,grid_n));
%     vals = arrayfun(f, grid);
%     idx  = find(vals(1:end-1).*vals(2:end) <= 0);
%     uC = [];
%     for k = idx(:).'
%         try, uC(end+1,1)=fzero(f,[grid(k),grid(k+1)]); catch, end %#ok<AGROW>
%     end
%     if isempty(uC), uBest=NaN; pBest=[NaN;NaN;NaN]; return; end
%     uC = unique(uC);
%     P   = ppval(pp,uC).';                     % M×3
%     d   = vecnorm(P - Ci.',2,2);              % radial distance (m)
%     speedSeed = max(norm(ppval(dpp,uSeed)), eps);
%     ds  = abs(uC - uSeed) * speedSeed;        % arclength deviation (m)
%     d(d>dist_cap) = inf;                      % width sanity
%     score = d + alpha*ds;
%     [minval,j] = min(score);
%     if ~isfinite(minval), uBest=NaN; pBest=[NaN;NaN;NaN]; return; end
%     uBest = uC(j); pBest = P(j,:).';
% end
% 
% function [uBest,pBest] = closestPointOnCurve(pp,Ci,uLo,uHi)
%     obj = @(u) sum((ppval(pp,u)-Ci).^2,1);
%     uBest = fminbnd(obj,uLo,uHi);
%     pBest = ppval(pp,uBest);
% end



% function [P_L, P_R, U_L, U_R, FL_at_UL, FR_at_UR] = ...
%     Compute_intersection_points(C_true, T_C, my_pp_L, my_pp_R, uL_of_S, uR_of_S)
% % Robust plane–rail intersections with hairpin fallback.
% % Signature is unchanged from your original function.
% 
% N = size(C_true,1);
% P_L = zeros(N,3);  P_R = zeros(N,3);
% U_L = zeros(N,1);  U_R = zeros(N,1);
% FL_at_UL = zeros(N,1); FR_at_UR = zeros(N,1);
% 
% U_L_max = my_pp_L.breaks(end);
% U_R_max = my_pp_R.breaks(end);
% 
% % ======= Tunables (set for Δs≈1 m, half-width≈6 m) =======
% W_half      = 6.0;             % [m] half track width near hairpin
% search_s    = 18; % 2.0*W_half      % [m] ± arclength window in hairpins
% dist_cap    = 1.5*W_half;      % [m] reject intersections farther than this
% angle_thr   = 25*pi/180;       % [rad] hairpin detector (tangent bend)
% grid_n      = 41;              % samples to bracket sign changes
% alpha_seed  = 0.30;            % weight on distance from uSeed (m per m)
% alpha_prev  = 2.50;  % 1.60          % weight on distance from previous accepted
% jump_cap_m  = 6.0;             % [m] hard cap from previous accepted
% side_tol_m  = 0.10;  % 0.20          % [m] side tolerance in XY test
% % =========================================================
% 
% % Derivatives for local speed (m per unit-u)
% dppL = fnder(my_pp_L,1);
% dppR = fnder(my_pp_R,1);
% 
% % Rolling continuity seeds = previously accepted u
% uPrevL = min(max(uL_of_S(1),0), U_L_max);
% uPrevR = min(max(uR_of_S(1),0), U_R_max);
% 
% for i = 1:N
%     Ci = C_true(i,:).';  Ti = T_C(i,:).';
% 
%     % fast evaluators (3D plane condition)
%     fL = @(u) (ppval(my_pp_L,u) - Ci).' * Ti;
%     fR = @(u) (ppval(my_pp_R,u) - Ci).' * Ti;
% 
%     % clamp initial seeds
%     u0L = min(max(uL_of_S(i),0), U_L_max);
%     u0R = min(max(uR_of_S(i),0), U_R_max);
% 
%     % simple hairpin detector from tangent change
%     if i>1 && i<N
%         c = max(-1,min(1, dot(T_C(i-1,:), T_C(i+1,:)) ));
%         ang = acos(c);
%         isHairpin = ang > angle_thr;
%     else
%         isHairpin = false;
%     end
% 
%     % local 2D normal (rotate tangent +90° in XY)
%     txy = T_C(i,1:2);  nxy = [-txy(2), txy(1)];
%     nxy = nxy / max(norm(nxy), eps);
% 
%     % arclength-to-parameter windows centered at u0*
%     speedL = max(norm(ppval(dppL,u0L)), eps);
%     speedR = max(norm(ppval(dppR,u0R)), eps);
%     duWinL = search_s / speedL;  duWinR = search_s / speedR;
%     uLoL = max(0, u0L - duWinL);  uHiL = min(U_L_max, u0L + duWinL);
%     uLoR = max(0, u0R - duWinR);  uHiR = min(U_R_max, u0R + duWinR);
% 
%     % ===== LEFT =====
%     uLroot = NaN;  pL = [NaN;NaN;NaN];
%     if ~isHairpin
%         try, uLroot = fzero(fL, u0L);  pL = ppval(my_pp_L,uLroot); catch, end
%     end
%     if isnan(uLroot) || norm(pL - Ci) > dist_cap
%         [uLroot,pL] = bestRootWindow(my_pp_L,dppL,fL,Ci, ...
%                         uLoL,uHiL,u0L,uPrevL,nxy,+1, ...
%                         alpha_seed,alpha_prev,grid_n,dist_cap,jump_cap_m,side_tol_m);
%     end
%     if isnan(uLroot)
%         [uLroot,pL] = closestPointOnCurve(my_pp_L,Ci,uLoL,uHiL);
%     end
%     P_L(i,:) = pL.';  U_L(i) = uLroot;  FL_at_UL(i) = fL(uLroot);
%     uPrevL = U_L(i);   % continuity update
% 
%     % ===== RIGHT =====
%     uRroot = NaN;  pR = [NaN;NaN;NaN];
%     if ~isHairpin
%         try, uRroot = fzero(fR, u0R);  pR = ppval(my_pp_R,uRroot); catch, end
%     end
%     if isnan(uRroot) || norm(pR - Ci) > dist_cap
%         [uRroot,pR] = bestRootWindow(my_pp_R,dppR,fR,Ci, ...
%                         uLoR,uHiR,u0R,uPrevR,nxy,-1, ...
%                         alpha_seed,alpha_prev,grid_n,dist_cap,jump_cap_m,side_tol_m);
%     end
%     if isnan(uRroot)
%         [uRroot,pR] = closestPointOnCurve(my_pp_R,Ci,uLoR,uHiR);
%     end
%     P_R(i,:) = pR.';  U_R(i) = uRroot;  FR_at_UR(i) = fR(uRroot);
%     uPrevR = U_R(i);   % continuity update
% end
% end
% 
% % ===== helpers (subfunctions) ============================================
% 
% function [uBest,pBest] = bestRootWindow(pp,dpp,f,Ci, ...
%         uLo,uHi,uSeed,uPrev,nxy,sideSign, ...
%         alpha_seed,alpha_prev,grid_n,dist_cap,jump_cap_m,side_tol_m)
% 
%     % 1) bracket all sign changes in [uLo,uHi]
%     grid = linspace(uLo,uHi,max(3,grid_n));
%     vals = arrayfun(f, grid);
%     idx  = find(vals(1:end-1).*vals(2:end) <= 0);
%     uC = [];
%     for k = idx(:).'
%         try, uC(end+1,1) = fzero(f,[grid(k),grid(k+1)]); catch, end %#ok<AGROW>
%     end
%     % widen once if nothing found
%     if isempty(uC)
%         grid = linspace(pp.breaks(1), pp.breaks(end), max(201,5*grid_n));
%         vals = arrayfun(f, grid);
%         idx  = find(vals(1:end-1).*vals(2:end) <= 0);
%         for k = idx(:).'
%             try, uC(end+1,1) = fzero(f,[grid(k),grid(k+1)]); catch, end %#ok<AGROW>
%         end
%     end
%     if isempty(uC), uBest=NaN; pBest=[NaN;NaN;NaN]; return; end
% 
%     % de-dup roots (tolerant)
%     if exist('uniquetol','file')
%         uC = uniquetol(uC, 1e-10);
%     else
%         epsu = 1e-10*max(1,abs(uC));
%         [uC,~] = unique(round(uC./epsu).*epsu);
%     end
% 
%     % 2) evaluate candidates and score
%     P   = ppval(pp,uC).';                      % M×3
%     d   = vecnorm(P - Ci.',2,2);               % distance to center (m)
% 
%     % side filter in XY (left: +1 keeps sXY>=-tol; right: -1 keeps sXY<=+tol)
%     sXY = (P(:,1)-Ci(1))*nxy(1) + (P(:,2)-Ci(2))*nxy(2);
%     if sideSign > 0
%         keep = (sXY >= -side_tol_m);
%     else
%         keep = (sXY <= +side_tol_m);
%     end
%     if any(keep), uC=uC(keep); P=P(keep,:); d=d(keep); else
%         % don’t drop all—just penalize
%         d = d + 50;
%     end
% 
%     % continuity penalties (convert Δu to meters using local speeds)
%     speedSeed = max(norm(ppval(dpp,uSeed)), eps);
%     speedPrev = max(norm(ppval(dpp,uPrev)), eps);
%     ds_seed   = abs(uC - uSeed) * speedSeed;
%     ds_prev   = abs(uC - uPrev) * speedPrev;
% 
%     % hard cap on jumps from previous accepted
%     ds_prev(ds_prev > jump_cap_m) = inf;
% 
%     % distance sanity
%     d(d > dist_cap) = inf;
% 
%     % composite score
%     score = d + alpha_seed*ds_seed + alpha_prev*ds_prev;
% 
%     [minval,j] = min(score);
%     if ~isfinite(minval), uBest=NaN; pBest=[NaN;NaN;NaN]; return; end
%     uBest = uC(j);  pBest = P(j,:).';
% end
% 
% function [uBest,pBest] = closestPointOnCurve(pp,Ci,uLo,uHi)
%     obj = @(u) sum((ppval(pp,u)-Ci).^2,1);
%     uBest = fminbnd(obj, uLo, uHi);
%     pBest = ppval(pp, uBest);
% end




% function [P_L, P_R, U_L, U_R, FL_at_UL, FR_at_UR] = ...
%     Compute_intersection_points(C_true, T_C, my_pp_L, my_pp_R, uL_of_S, uR_of_S)
% 
% N = size(C_true,1);
% P_L = zeros(N,3);  P_R = zeros(N,3);
% U_L = zeros(N,1);  U_R = zeros(N,1);
% FL_at_UL = zeros(N,1); FR_at_UR = zeros(N,1);
% 
% U_L_max = my_pp_L.breaks(end);
% U_R_max = my_pp_R.breaks(end);
% 
% % ---- tunables ----
% W_half      = 6.0;           % m
% search_s    = 15.0;          % ± arclength window (m)
% dist_cap    = 9.0;           % m
% angle_thr   = 25*pi/180;     % rad (hairpin detector)
% grid_n      = 41;
% 
% alpha_seed  = 0.30;          % weight on |Δu| from uSeed (m/m)
% alpha_prev  = 2.20;          % weight on |Δu| from previous accepted (m/m)
% jump_cap_m  = 3.0;           % hard cap on arclength jump (m)
% 
% beta_lat_prev = 1.00;        % weight on |sXY - sXY_prev| (m/m)
% beta_lat_nom  = 0.50;        % weight on ||sXY|-W_half|  (m/m)
% side_tol_m    = 0.10;        % side tolerance (m)
% % ------------------
% 
% dppL = fnder(my_pp_L,1);
% dppR = fnder(my_pp_R,1);
% 
% % rolling continuity seeds
% uPrevL = min(max(uL_of_S(1),0), U_L_max); sPrevL = +W_half;
% uPrevR = min(max(uR_of_S(1),0), U_R_max); sPrevR = -W_half;
% 
% for i = 1:N
%     Ci = C_true(i,:).';                     % 3x1
%     Ti = T_C(i,:).';                        % 3x1
%     fL = @(u) (ppval(my_pp_L,u) - Ci).' * Ti;
%     fR = @(u) (ppval(my_pp_R,u) - Ci).' * Ti;
% 
%     u0L = min(max(uL_of_S(i),0), U_L_max);
%     u0R = min(max(uR_of_S(i),0), U_R_max);
% 
%     % hairpin detector
%     if i>1 && i<N
%         ang = acos(max(-1,min(1, dot(T_C(i-1,:), T_C(i+1,:)) )));
%         isHairpin = ang > angle_thr;
%     else
%         isHairpin = false;
%     end
% 
%     % 2D normal (+90° rotation of tangent in XY)
%     txy = T_C(i,1:2);
%     nxy = [-txy(2), txy(1)];
%     nxy = nxy / max(norm(nxy), eps);        % 1x2
% 
%     % convert ±search_s (m) to ±Δu around seed
%     speedL = max(norm(ppval(dppL,u0L)), eps);
%     speedR = max(norm(ppval(dppR,u0R)), eps);
%     duWinL = search_s / speedL;  duWinR = search_s / speedR;
%     uLoL = max(0, u0L - duWinL);  uHiL = min(U_L_max, u0L + duWinL);
%     uLoR = max(0, u0R - duWinR);  uHiR = min(U_R_max, u0R + duWinR);
% 
%     % ===== LEFT =====
%     uLroot = NaN; pL = [NaN;NaN;NaN];
%     if ~isHairpin
%         try
%             uLroot = fzero(fL, u0L);
%             pL = ppval(my_pp_L,uLroot); pL = pL(:);      % force 3x1
%         catch
%         end
%     end
%     if isnan(uLroot) || norm(pL - Ci) > dist_cap
%         [uLroot,pL,sXY_L] = bestRootWindow( ...
%             my_pp_L,dppL,fL,Ci,uLoL,uHiL,u0L,uPrevL,nxy,+1, ...
%             alpha_seed,alpha_prev,beta_lat_prev,beta_lat_nom, ...
%             grid_n,dist_cap,jump_cap_m,W_half,side_tol_m, sPrevL);
%     else
%         sXY_L = dot(pL(1:2)-Ci(1:2), nxy(:));            % scalar
%     end
%     if isnan(uLroot)
%         [uLroot,pL] = closestPointOnCurve(my_pp_L,Ci,uLoL,uHiL);
%         pL = pL(:);  sXY_L = dot(pL(1:2)-Ci(1:2), nxy(:));
%     end
%     P_L(i,:) = pL.';  U_L(i)=uLroot;  FL_at_UL(i)=fL(uLroot);
%     uPrevL = U_L(i);  sPrevL = sXY_L;
% 
%     % ===== RIGHT =====
%     uRroot = NaN; pR = [NaN;NaN;NaN];
%     if ~isHairpin
%         try
%             uRroot = fzero(fR, u0R);
%             pR = ppval(my_pp_R,uRroot); pR = pR(:);
%         catch
%         end
%     end
%     if isnan(uRroot) || norm(pR - Ci) > dist_cap
%         [uRroot,pR,sXY_R] = bestRootWindow( ...
%             my_pp_R,dppR,fR,Ci,uLoR,uHiR,u0R,uPrevR,nxy,-1, ...
%             alpha_seed,alpha_prev,beta_lat_prev,beta_lat_nom, ...
%             grid_n,dist_cap,jump_cap_m,W_half,side_tol_m, sPrevR);
%     else
%         sXY_R = dot(pR(1:2)-Ci(1:2), nxy(:));
%     end
%     if isnan(uRroot)
%         [uRroot,pR] = closestPointOnCurve(my_pp_R,Ci,uLoR,uHiR);
%         pR = pR(:);  sXY_R = dot(pR(1:2)-Ci(1:2), nxy(:));
%     end
%     P_R(i,:) = pR.';  U_R(i)=uRroot;  FR_at_UR(i)=fR(uRroot);
%     uPrevR = U_R(i);  sPrevR = sXY_R;
% end
% end
% 
% % ==== helpers =============================================================
% 
% function [uBest,pBest,sXY_best] = bestRootWindow(pp,dpp,f,Ci, ...
%         uLo,uHi,uSeed,uPrev,nxy,sideSign, ...
%         alpha_seed,alpha_prev,beta_lat_prev,beta_lat_nom, ...
%         grid_n,dist_cap,jump_cap_m,W_half,side_tol_m, sPrev)
% 
%     % 1) bracket sign-changes
%     grid = linspace(uLo,uHi,max(3,grid_n));
%     vals = arrayfun(f, grid);
%     idx  = find(vals(1:end-1).*vals(2:end) <= 0);
%     uC = [];
%     for k = idx(:).'
%         try, uC(end+1,1) = fzero(f,[grid(k),grid(k+1)]); catch, end %#ok<AGROW>
%     end
%     if isempty(uC)
%         big = linspace(pp.breaks(1), pp.breaks(end), max(201,5*grid_n));
%         vals = arrayfun(f, big);
%         idx  = find(vals(1:end-1).*vals(2:end) <= 0);
%         for k = idx(:).'
%             try, uC(end+1,1) = fzero(f,[big(k),big(k+1)]); catch, end %#ok<AGROW>
%         end
%     end
%     if isempty(uC), uBest=NaN; pBest=[NaN;NaN;NaN]; sXY_best=NaN; return; end
% 
%     % de-dupe
%     if exist('uniquetol','file')
%         uC = uniquetol(uC, 1e-10);
%     else
%         epsu = 1e-10*max(1,abs(uC));
%         [uC,~] = unique(round(uC./epsu).*epsu);
%     end
% 
%     % 2) evaluate and score
%     P   = ppval(pp,uC).';                                % Mx3
%     d   = vecnorm(P - Ci.',2,2);                         % Mx1
% 
%     sXY = (P(:,1)-Ci(1))*nxy(1) + (P(:,2)-Ci(2))*nxy(2); % signed lateral (m)
%     if sideSign > 0, keep = (sXY >= -side_tol_m);
%     else,            keep = (sXY <= +side_tol_m);
%     end
%     if any(keep), uC=uC(keep); P=P(keep,:); d=d(keep); sXY=sXY(keep);
%     else, d = d + 50;  % penalize if none strictly satisfy side
%     end
% 
%     speedSeed = max(norm(ppval(dpp,uSeed)), eps);
%     speedPrev = max(norm(ppval(dpp,uPrev)), eps);
%     ds_seed   = abs(uC - uSeed) * speedSeed;
%     ds_prev   = abs(uC - uPrev) * speedPrev;
% 
%     ds_prev(ds_prev > jump_cap_m) = inf;
%     d(d > dist_cap)               = inf;
% 
%     lat_prev = abs(sXY - sPrev);
%     lat_nom  = abs(abs(sXY) - W_half);
% 
%     score = d + alpha_seed*ds_seed + alpha_prev*ds_prev ...
%               + beta_lat_prev*lat_prev + beta_lat_nom*lat_nom;
% 
%     [minval,j] = min(score);
%     if ~isfinite(minval), uBest=NaN; pBest=[NaN;NaN;NaN]; sXY_best=NaN; return; end
%     uBest = uC(j);  pBest = P(j,:).';  sXY_best = sXY(j);
% end
% 
% function [uBest,pBest] = closestPointOnCurve(pp,Ci,uLo,uHi)
%     obj = @(u) sum((ppval(pp,u)-Ci).^2,1);   % scalar objective
%     uBest = fminbnd(obj, uLo, uHi);
%     pBest = ppval(pp, uBest);
% end






function [P_L, P_R, U_L, U_R, FL_at_UL, FR_at_UR] = ...
    Compute_intersection_points(C_true, T_C, my_pp_L, my_pp_R, uL_of_S, uR_of_S)

N = size(C_true,1);
P_L = zeros(N,3);  P_R = zeros(N,3);
U_L = zeros(N,1);  U_R = zeros(N,1);
FL_at_UL = zeros(N,1); FR_at_UR = zeros(N,1);

U_L_max = my_pp_L.breaks(end);
U_R_max = my_pp_R.breaks(end);

% ---- tunables (tightened) ----
W_half      = 6.0;           % m
search_s    = 15.0;          % ± arclength window (m)
dist_cap    = 9.0;           % m
angle_thr   = 25*pi/180;     % rad
grid_n      = 41;

alpha_seed  = 0.30;          % m/m
alpha_prev  = 3.50;          % ↑ stronger continuity in u (m/m)
jump_cap_m  = 1.2;           % ↓ smaller allowed jump along rail (m)

beta_lat_prev = 2.50;        % ↑ keep signed lateral close to previous (m/m)
beta_lat_nom  = 0.50;        % keep |lateral| near W_half (m/m)
side_tol_m    = 0.05;        % stricter side test (m)

gamma_xyz_prev = 1.50;       % NEW: weight on ||P - P_prev|| (m/m)
step_cap_m     = 1.50;       % NEW: hard cap on that 3D step (m)
% -----------------------------

dppL = fnder(my_pp_L,1);
dppR = fnder(my_pp_R,1);

% rolling continuity seeds
uPrevL = min(max(uL_of_S(1),0), U_L_max); sPrevL = +W_half; pPrevL = ppval(my_pp_L,uPrevL).';
uPrevR = min(max(uR_of_S(1),0), U_R_max); sPrevR = -W_half; pPrevR = ppval(my_pp_R,uPrevR).';

for i = 1:N
    Ci = C_true(i,:).';  Ti = T_C(i,:).';
    fL = @(u) (ppval(my_pp_L,u) - Ci).' * Ti;
    fR = @(u) (ppval(my_pp_R,u) - Ci).' * Ti;

    u0L = min(max(uL_of_S(i),0), U_L_max);
    u0R = min(max(uR_of_S(i),0), U_R_max);

    % hairpin detector
    if i>1 && i<N
        ang = acos(max(-1,min(1, dot(T_C(i-1,:), T_C(i+1,:)) )));
        isHairpin = ang > angle_thr;
    else
        isHairpin = false;
    end

    % 2D normal (+90° rotation)
    txy = T_C(i,1:2);  nxy = [-txy(2), txy(1)];
    nxy = nxy / max(norm(nxy), eps);

    % ±search_s (m) -> ±Δu around seed
    speedL = max(norm(ppval(dppL,u0L)), eps);
    speedR = max(norm(ppval(dppR,u0R)), eps);
    duWinL = search_s / speedL;  duWinR = search_s / speedR;
    uLoL = max(0, u0L - duWinL);  uHiL = min(U_L_max, u0L + duWinL);
    uLoR = max(0, u0R - duWinR);  uHiR = min(U_R_max, u0R + duWinR);

    % ===== LEFT =====
    uLroot = NaN; pL = [NaN;NaN;NaN];
    if ~isHairpin
        try, uLroot = fzero(fL, u0L);  pL = ppval(my_pp_L,uLroot); pL=pL(:); catch, end
    end
    if isnan(uLroot) || norm(pL - Ci) > dist_cap
        [uLroot,pL,sXY_L] = bestRootWindow( ...
            my_pp_L,dppL,fL,Ci,uLoL,uHiL,u0L,uPrevL,nxy,+1, ...
            alpha_seed,alpha_prev,beta_lat_prev,beta_lat_nom, ...
            gamma_xyz_prev,grid_n,dist_cap,jump_cap_m,step_cap_m, ...
            W_half,side_tol_m, sPrevL, pPrevL);
    else
        sXY_L = dot(pL(1:2)-Ci(1:2), nxy(:));
    end
    if isnan(uLroot)
        [uLroot,pL] = closestPointOnCurve(my_pp_L,Ci,uLoL,uHiL); pL=pL(:);
        sXY_L = dot(pL(1:2)-Ci(1:2), nxy(:));
    end
    P_L(i,:) = pL.'; U_L(i)=uLroot; FL_at_UL(i)=fL(uLroot);
    uPrevL = U_L(i); sPrevL = sXY_L; pPrevL = pL;

    % ===== RIGHT =====
    uRroot = NaN; pR = [NaN;NaN;NaN];
    if ~isHairpin
        try, uRroot = fzero(fR, u0R);  pR = ppval(my_pp_R,uRroot); pR=pR(:); catch, end
    end
    if isnan(uRroot) || norm(pR - Ci) > dist_cap
        [uRroot,pR,sXY_R] = bestRootWindow( ...
            my_pp_R,dppR,fR,Ci,uLoR,uHiR,u0R,uPrevR,nxy,-1, ...
            alpha_seed,alpha_prev,beta_lat_prev,beta_lat_nom, ...
            gamma_xyz_prev,grid_n,dist_cap,jump_cap_m,step_cap_m, ...
            W_half,side_tol_m, sPrevR, pPrevR);
    else
        sXY_R = dot(pR(1:2)-Ci(1:2), nxy(:));
    end
    if isnan(uRroot)
        [uRroot,pR] = closestPointOnCurve(my_pp_R,Ci,uLoR,uHiR); pR=pR(:);
        sXY_R = dot(pR(1:2)-Ci(1:2), nxy(:));
    end
    P_R(i,:) = pR.'; U_R(i)=uRroot; FR_at_UR(i)=fR(uRroot);
    uPrevR = U_R(i); sPrevR = sXY_R; pPrevR = pR;
end
end

% ==== helpers =============================================================
function [uBest,pBest,sXY_best] = bestRootWindow(pp,dpp,f,Ci, ...
        uLo,uHi,uSeed,uPrev,nxy,sideSign, ...
        alpha_seed,alpha_prev,beta_lat_prev,beta_lat_nom, ...
        gamma_xyz_prev,grid_n,dist_cap,jump_cap_m,step_cap_m, ...
        W_half,side_tol_m, sPrev, pPrev)

    grid = linspace(uLo,uHi,max(3,grid_n));
    vals = arrayfun(f, grid);
    idx  = find(vals(1:end-1).*vals(2:end) <= 0);
    uC = [];
    for k = idx(:).'   % refine each bracket
        try, uC(end+1,1) = fzero(f,[grid(k),grid(k+1)]); catch, end %#ok<AGROW>
    end
    if isempty(uC)
        big = linspace(pp.breaks(1), pp.breaks(end), max(201,5*grid_n));
        vals = arrayfun(f, big);
        idx  = find(vals(1:end-1).*vals(2:end) <= 0);
        for k = idx(:).'
            try, uC(end+1,1) = fzero(f,[big(k),big(k+1)]); catch, end %#ok<AGROW>
        end
    end
    if isempty(uC), uBest=NaN; pBest=[NaN;NaN;NaN]; sXY_best=NaN; return; end

    % de-dup
    if exist('uniquetol','file'), uC = uniquetol(uC,1e-10);
    else
        epsu = 1e-10*max(1,abs(uC)); [uC,~] = unique(round(uC./epsu).*epsu);
    end

    P   = ppval(pp,uC).';                              % Mx3
    d   = vecnorm(P - Ci.',2,2);                       % center distance (m)

    % side in XY
    sXY = (P(:,1)-Ci(1))*nxy(1) + (P(:,2)-Ci(2))*nxy(2);
    keep = sideSign>0 & (sXY >= -side_tol_m) | sideSign<0 & (sXY <= +side_tol_m);
    if any(keep), uC=uC(keep); P=P(keep,:); d=d(keep); sXY=sXY(keep);
    else, d = d + 50; end

    % continuity in u -> meters
    speedSeed = max(norm(ppval(dpp,uSeed)), eps);
    speedPrev = max(norm(ppval(dpp,uPrev)), eps);
    ds_seed   = abs(uC - uSeed) * speedSeed;
    ds_prev   = abs(uC - uPrev) * speedPrev;
    ds_prev(ds_prev > jump_cap_m) = inf;

    % NEW: continuity in 3D space to previous accepted rail point
    dxyz_prev = vecnorm(P - pPrev.',2,2);
    dxyz_prev(dxyz_prev > step_cap_m) = inf;

    % sanity caps
    d(d > dist_cap) = inf;

    % lateral consistency
    lat_prev = abs(sXY - sPrev);
    lat_nom  = abs(abs(sXY) - W_half);

    score = d ...
          + alpha_seed*ds_seed + alpha_prev*ds_prev ...
          + beta_lat_prev*lat_prev + beta_lat_nom*lat_nom ...
          + gamma_xyz_prev*dxyz_prev;     % NEW term

    [minval,j] = min(score);
    if ~isfinite(minval), uBest=NaN; pBest=[NaN;NaN;NaN]; sXY_best=NaN; return; end
    uBest = uC(j);  pBest = P(j,:).';  sXY_best = sXY(j);
end

function [uBest,pBest] = closestPointOnCurve(pp,Ci,uLo,uHi)
    obj = @(u) sum((ppval(pp,u)-Ci).^2,1);
    uBest = fminbnd(obj, uLo, uHi);
    pBest = ppval(pp, uBest);
end




















%%
%
function [dTheta, dMu, dPhi, ddTheta, ddMu, ddPhi] = Estimate_dEulerds_ddEulerds2(S_all_points, Theta, Mu, Phi)

p = 1;
% spTh = myspcsp(S_all_points.', Theta.', p);
spTh = csaps(S_all_points.', Theta.', p);
spMu = myspcsp(S_all_points.', Mu.',    p);
spPh = myspcsp(S_all_points.', Phi.',   p);

% 2) form the derivative splines
dspTh  = fnder(spTh,1);
dspMu  = fnder(spMu,1);
dspPh  = fnder(spPh,1);

% 3) and even the second‐derivative splines
ddspTh = fnder(spTh,2);
ddspMu = fnder(spMu,2);
ddspPh = fnder(spPh,2);

% 4) now sample
dTheta    = ppval(dspTh,  S_all_points);
dMu       = ppval(dspMu,  S_all_points);
dPhi      = ppval(dspPh,  S_all_points);

ddTheta   = ppval(ddspTh, S_all_points);
ddMu      = ppval(ddspMu, S_all_points);
ddPhi     = ppval(ddspPh, S_all_points);

% debug = 1;

end

%%
%
function [dnl, dnr] = Estimate_dnlds_dnrds(S_all_points, nl, nr)

% 1) fit (smoothing) splines
p = 1;
spnl = myspcsp(S_all_points.', nl.', p);
spnr = myspcsp(S_all_points.', nr.', p);

% 2) form the derivative splines
dspnl  = fnder(spnl,1);
dspnr  = fnder(spnr,1);

% 4) now sample
dnl  = ppval(dspnl,  S_all_points);
dnr  = ppval(dspnr,  S_all_points);

end

%%
%
function All_Fs_table = ODEs_rhs(cos_Theta, sin_Theta, cos_Mu, sin_Mu, cos_Phi, sin_Phi,...
                                            dTheta, dMu, dPhi, ddTheta, ddMu, ddPhi, dnl, dnr)


F1 = cos_Theta .* cos_Mu;
F2 = sin_Theta .* cos_Mu;
F3 = -sin_Mu;

F4 = -sin_Theta .* dTheta;
F5 =  cos_Theta .* dTheta;

F6 = -sin_Mu .* dMu;
F7 =  cos_Mu .* dMu;

F8 = -sin_Phi .* dPhi;
F9 =  cos_Phi .* dPhi;

F10 = ddTheta;
F11 = ddMu;
F12 = ddPhi;

F13 = dnl;
F14 = dnr;

All_Fs_table = [F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14];

end

%%
%
function discrete_ODE_residuals = discrete_ODE_residuals(All_Fs_table, All_X_table, D_bar_kron_I)


All_Fs_table_except_t1 = All_Fs_table(2:end, :);

All_Fs_table_except_t1_reshaped = reshape(All_Fs_table_except_t1.', [], 1);

All_X_table_reshaped = reshape(All_X_table.', [], 1);

discrete_ODE_residuals = D_bar_kron_I * All_X_table_reshaped - All_Fs_table_except_t1_reshaped;



end

%%
%
function [D_bar_kron_I, D_bar] = Compute_Dbar_kron_I(N_nodes, sa, nx, A_matrix_RadauIIA, h_vec)

%
nmesh = N_nodes;

ninterval = nmesh-1;
n_all_collocation_points = sa * ninterval; % number of all collocation points

% D matrix
D = A_matrix_RadauIIA \ [-1*ones(sa, 1), speye(sa)];

% D_bar
D_bar   = sparse(n_all_collocation_points, n_all_collocation_points + 1); 
for i = 1:ninterval
    % Which rows to fill for the i-th block:
    row_idx = (i-1)*sa + (1:sa);
    
    % Which columns to fill:
    %   Block 1 -> columns 1..(s+1)
    %   Block 2 -> columns s+1..(2s+1)
    %   etc.
    col_start = (i-1)*sa + 1;
    col_idx   = col_start : (col_start + sa);  % s+1 columns
    
    % Scale D by 1/h_vec(i) and place it in D_bar
    D_bar(row_idx, col_idx) = D / h_vec(i);
end

D_bar_kron_I = kron(D_bar, speye(nx));

end

%%
%
function [tracking_penalty_1, tracking_penalty_2, tracking_penalty_3, x_CL, y_CL, z_CL, x_n_hat, y_n_hat, z_n_hat] = tracking_error(All_X_table, x_C_m, y_C_m, z_C_m, x_CL_m, y_CL_m, z_CL_m, x_CR_m, y_CR_m, z_CR_m)

x_C = All_X_table(:,1);
y_C = All_X_table(:,2);
z_C = All_X_table(:,3);

c_theta = All_X_table(:,4);
s_theta = All_X_table(:,5);

c_mu    = All_X_table(:,6);
s_mu    = All_X_table(:,7);

c_phi   = All_X_table(:,8);
s_phi   = All_X_table(:,9);

dtheta = All_X_table(:,10);
dmu    = All_X_table(:,11);
dphi   = All_X_table(:,12);

nl = All_X_table(:,13);
nr = All_X_table(:,14);

% coordinates of vector n_hat, Giacomo thesis Eq(4.19) matrix second column
x_n_hat = c_theta.*s_mu.*s_phi - s_theta.*c_phi;
y_n_hat = s_theta.*s_mu.*s_phi + c_theta.*c_phi;
z_n_hat = c_mu.*s_phi;

% coordinates of the left boundary of ribbon
x_CL = x_C + x_n_hat.*nl;
y_CL = y_C + y_n_hat.*nl;
z_CL = z_C + z_n_hat.*nl;

% coordinates of the right boundary of ribbon
x_CR = x_C + x_n_hat.*nr;
y_CR = y_C + y_n_hat.*nr;
z_CR = z_C + z_n_hat.*nr;

r_x_C = (x_C - x_C_m);
r_y_C = (y_C - y_C_m);
r_z_C = (z_C - z_C_m);

r_x_CL = (x_CL - x_CL_m);
r_y_CL = (y_CL - y_CL_m);
r_z_CL = (z_CL - z_CL_m);

r_x_CR = (x_CR - x_CR_m);
r_y_CR = (y_CR - y_CR_m);
r_z_CR = (z_CR - z_CR_m);

%
tracking_penalty_1 = (r_x_C.^2  + r_y_C.^2  + r_z_C.^2);
tracking_penalty_2 = (r_x_CL.^2 + r_y_CL.^2 + r_z_CL.^2);
tracking_penalty_3 = (r_x_CR.^2 + r_y_CR.^2 + r_z_CR.^2);

% tracking_penalty = tracking_penalty_1 + tracking_penalty_2 + tracking_penalty_3;

end
