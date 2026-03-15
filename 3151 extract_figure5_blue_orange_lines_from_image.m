
%% extract_figure5_blue_orange_lines_from_image.m
% Extract the blue "Grip Driving" line and the orange
% "Handbrake Destabilization" line from Figure 5 in one run.
%
% IMPORTANT
% - This is an approximate digitization from the published figure image.
% - Figure 5 only contains East/North plan-view information.
%   True Z is not available from the figure, so z_m is written as NaN.
% - This script uses deterministic HSV color thresholding + row-wise median
%   centerline extraction on the right-hand ROI, where the blue/orange
%   segments appear in Figure 5.
%
% OUTPUTS
%   figure5_grip_driving_blue_line_extracted.csv
%   figure5_handbrake_orange_line_extracted.csv
%   figure5_blue_orange_extracted_combined.csv
%   figure5_grip_driving_blue_line_pixels.csv
%   figure5_handbrake_orange_line_pixels.csv
%   figure5_blue_orange_overlay_full_preview.png
%
% REQUIRES
%   MATLAB base functions only (no special toolbox assumed).
%
% AUTHOR NOTE
%   This script is designed to reproduce the provided CSV files closely.
%   If you need exact byte-for-byte reproduction, use the companion script:
%       extract_figure5_blue_orange_lines_exact.m

clear; close all; clc;

%% ========================= USER SETTINGS ===============================
imgFile = 'figure5_plot_crop.png';

csvBlue       = 'figure5_grip_driving_blue_line_extracted.csv';
csvOrange     = 'figure5_handbrake_orange_line_extracted.csv';
csvBluePix    = 'figure5_grip_driving_blue_line_pixels.csv';
csvOrangePix  = 'figure5_handbrake_orange_line_pixels.csv';
csvCombined   = 'figure5_blue_orange_extracted_combined.csv';
pngOverlay    = 'figure5_blue_orange_overlay_full_preview.png';

% Figure 5 plot calibration
east_min  = -780.0;
east_max  = -620.0;
north_min = -380.0;
north_max = -275.0;

% Right-side ROI where the blue/orange lines live in Figure 5
x_roi_min = 1101;  % MATLAB 1-based column index

% HSV thresholds (MATLAB hue is in [0,1])
% Blue: roughly OpenCV hue 100..120 -> MATLAB hue 100/180 .. 120/180
blue_h_min    = 100/180;
blue_h_max    = 120/180;
blue_sat_min  = 100/255;
blue_val_min  =  70/255;

% Orange: roughly OpenCV hue 8..25 -> MATLAB hue 8/180 .. 25/180
orange_h_min   =  8/180;
orange_h_max   = 25/180;
orange_sat_min = 100/255;
orange_val_min =  70/255;

blue_smooth_window   = 9;  % odd integer
orange_smooth_window = 7;  % odd integer
ds = 0.5;                  % output arc-length spacing [m]

%% ========================= LOAD IMAGE ==================================
I = imread(imgFile);
if size(I,3) ~= 3
    error('Input image must be RGB.');
end

[Himg, Wimg, ~] = size(I);
Ihsv = rgb2hsv(I);
H = Ihsv(:,:,1);
S = Ihsv(:,:,2);
V = Ihsv(:,:,3);

roiMask = false(Himg, Wimg);
roiMask(:, x_roi_min:end) = true;

%% ========================= COLOR MASKS =================================
blueMask = ...
    (H >= blue_h_min) & (H <= blue_h_max) & ...
    (S >  blue_sat_min) & (V > blue_val_min) & roiMask;

orangeMask = ...
    (H >= orange_h_min) & (H <= orange_h_max) & ...
    (S >  orange_sat_min) & (V > orange_val_min) & roiMask;

%% ========================= EXTRACT CENTERLINES ==========================
[xb_pix, yb_pix] = rowwiseCenterline(blueMask, blue_smooth_window);
[xo_pix, yo_pix] = rowwiseCenterline(orangeMask, orange_smooth_window);

%% ========================= SAVE PIXEL CSVs =============================
TbPix = table((1:numel(xb_pix))', xb_pix(:), yb_pix(:), ...
    'VariableNames', {'point_id','x_pix','y_pix'});
ToPix = table((1:numel(xo_pix))', xo_pix(:), yo_pix(:), ...
    'VariableNames', {'point_id','x_pix','y_pix'});

writetable(TbPix, csvBluePix);
writetable(ToPix, csvOrangePix);

%% ========================= PIXELS -> WORLD =============================
[east_b, north_b] = pix2world(xb_pix, yb_pix, Wimg, Himg, ...
    east_min, east_max, north_min, north_max);

[east_o, north_o] = pix2world(xo_pix, yo_pix, Wimg, Himg, ...
    east_min, east_max, north_min, north_max);

%% ========================= RESAMPLE BY ARC LENGTH ======================
[east_b_r, north_b_r, s_b] = resampleByArcLength(east_b, north_b, ds);
[east_o_r, north_o_r, s_o] = resampleByArcLength(east_o, north_o, ds);

%% ========================= WRITE WORLD CSVs ============================
z_b = nan(size(s_b));
z_o = nan(size(s_o));

Tb = table((1:numel(s_b))', east_b_r(:), north_b_r(:), z_b(:), s_b(:), ...
    'VariableNames', {'point_id','east_m','north_m','z_m','s_m'});

To = table((1:numel(s_o))', east_o_r(:), north_o_r(:), z_o(:), s_o(:), ...
    'VariableNames', {'point_id','east_m','north_m','z_m','s_m'});

writetable(Tb, csvBlue);
writetable(To, csvOrange);

curve_id_b = ones(height(Tb),1);
curve_id_o = 2*ones(height(To),1);

curve_name_b = repmat("grip_driving_blue", height(Tb), 1);
curve_name_o = repmat("handbrake_destabilization_orange", height(To), 1);

Tbc = addvars(Tb, curve_name_b, curve_id_b, 'Before', 1, ...
    'NewVariableNames', {'curve_name','curve_id'});
Toc = addvars(To, curve_name_o, curve_id_o, 'Before', 1, ...
    'NewVariableNames', {'curve_name','curve_id'});

Tcombined = [Tbc; Toc];
writetable(Tcombined, csvCombined);

%% ========================= OVERLAY PREVIEW =============================
fig = figure('Color','w','Position',[100 100 1410 899]);
imshow(I); hold on;
plot(xb_pix, yb_pix, 'c-', 'LineWidth', 2);
plot(xo_pix, yo_pix, 'y-', 'LineWidth', 2);
hold off;
exportgraphics(gca, pngOverlay, 'Resolution', 150);

disp('Done.');
disp(['Wrote: ', csvBlue]);
disp(['Wrote: ', csvOrange]);
disp(['Wrote: ', csvCombined]);
disp(['Wrote: ', csvBluePix]);
disp(['Wrote: ', csvOrangePix]);
disp(['Wrote: ', pngOverlay]);

%% ========================= LOCAL FUNCTIONS =============================
function [x_full_smooth, y_full] = rowwiseCenterline(mask, smoothWindow)
    [H, ~] = size(mask);

    y_list = [];
    x_list = [];

    for y = 1:H
        cols = find(mask(y,:));
        if numel(cols) >= 2
            y_list(end+1,1) = y; %#ok<AGROW>
            x_list(end+1,1) = median(cols);
        end
    end

    if isempty(y_list)
        error('No pixels found for one of the color masks.');
    end

    y_full = (min(y_list):max(y_list))';
    x_full = interp1(y_list, x_list, y_full, 'linear');

    % Simple moving-average smoothing
    smoothWindow = max(3, round(smoothWindow));
    if mod(smoothWindow,2) == 0
        smoothWindow = smoothWindow + 1;
    end
    kernel = ones(smoothWindow,1) / smoothWindow;
    pad = floor(smoothWindow/2);

    x_pad = [repmat(x_full(1), pad, 1); x_full; repmat(x_full(end), pad, 1)];
    x_full_smooth = conv(x_pad, kernel, 'valid');
end

function [east, north] = pix2world(xpix, ypix, Wimg, Himg, ...
        east_min, east_max, north_min, north_max)

    east = east_min + (xpix - 1) / (Wimg - 1) * (east_max - east_min);
    north = north_max - (ypix - 1) / (Himg - 1) * (north_max - north_min);
end

function [xr, yr, s_new] = resampleByArcLength(x, y, ds)
    pts = [x(:), y(:)];
    d = sqrt(sum(diff(pts,1,1).^2, 2));
    s = [0; cumsum(d)];

    if s(end) <= 0
        xr = x(:);
        yr = y(:);
        s_new = s(:);
        return;
    end

    s_new = (0:ds:s(end))';
    if abs(s_new(end) - s(end)) > 1e-9
        s_new = [s_new; s(end)];
    end

    xr = interp1(s, x(:), s_new, 'linear');
    yr = interp1(s, y(:), s_new, 'linear');
end
