%% Load left/right boundaries and plot with equal axis scaling (3D)
% clc; clear; close all;

% --- set paths (change if needed) ---
leftPath  = "edge_loop_00.csv";

% --- read CSVs (expected columns: X, Y(up), Z, PoT) ---
L = readmatrix(leftPath);

% drop any NaN/Inf rows
L = L(all(isfinite(L),2),:);

% take first 3 columns as coordinates
XL = L(:,1); YL = L(:,2); ZL = L(:,3);

% --- plot ---
figure('Color','w'); hold on;
plot3(XL, YL, ZL, '-', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
view(3);