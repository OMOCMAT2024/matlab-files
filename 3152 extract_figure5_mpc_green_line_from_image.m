%% extract_figure5_mpc_green_line_from_image.m
% MATLAB script to extract the MPC green curve from Figure 5 in one run.
%
% Paper:
%   Beyond the stable handling limits: nonlinear model predictive control
%   for highly transient autonomous drifting
%
% WHAT THIS SCRIPT DOES
% 1) Reads the cropped Figure 5 image.
% 2) Thresholds the green MPC curve.
% 3) Converts extracted pixel coordinates to East/North.
% 4) Resamples the curve to uniform arc-length spacing.
% 5) Writes CSV and an overlay preview.
%
% IMPORTANT
% - This is an approximate digitization from the published figure.
% - Figure 5 is 2D only, so true Z cannot be recovered. z_m is NaN.
% - The extraction excludes the lower-left legend box region.
% - Because this is an image-based extraction, tiny numerical differences
%   can occur if you change the crop, resolution, or thresholds.
%
% REQUIRED FILE
%   figure5_plot_crop.png
%
% OUTPUT FILES
%   figure5_mpc_green_line_extracted_from_image.csv
%   figure5_mpc_green_line_pixels_from_image.csv
%   figure5_mpc_green_line_overlay_preview_from_image_matlab.png

clear; close all; clc;

%% User settings
imgFile   = 'figure5_plot_crop.png';
outCsv    = 'figure5_mpc_green_line_extracted_from_image.csv';
outPixCsv = 'figure5_mpc_green_line_pixels_from_image.csv';
outPng    = 'figure5_mpc_green_line_overlay_preview_from_image_matlab.png';
ds_target = 0.5;   % [m]

east_min  = -780.0;
east_max  = -620.0;
north_min = -380.0;
north_max = -275.0;

% Thresholds tuned for the provided Figure 5 crop
g_min = 90;
green_ratio = 1.35;

% Legend exclusion region (MATLAB 1-based indexing)
legend_row_start = 701;
legend_col_end   = 370;

%% Read image
assert(isfile(imgFile), 'Image file not found: %s', imgFile);
I = imread(imgFile);
[H, W, ~] = size(I);

R = double(I(:,:,1));
G = double(I(:,:,2));
B = double(I(:,:,3));

%% Green mask
greenMask = (G > g_min) & (G > green_ratio*R) & (G > green_ratio*B);

% Exclude lower-left legend region
greenMask(legend_row_start:end, 1:legend_col_end) = false;

%% Extract one y-value per x-column using column-wise median
x_present = find(any(greenMask, 1));
assert(~isempty(x_present), 'No green pixels found. Check image and thresholds.');

y_med = nan(size(x_present));
for k = 1:numel(x_present)
    xk = x_present(k);
    rows = find(greenMask(:, xk));
    y_med(k) = median(rows);
end

%% Fill internal gaps and smooth
x_full = (x_present(1):x_present(end)).';
y_full = interp1(x_present(:), y_med(:), x_full, 'linear', 'extrap');

% Smooth deterministically
if exist('sgolayfilt','file') == 2
    y_smooth = sgolayfilt(y_full, 3, 31);
else
    % Fallback if Signal Processing Toolbox is unavailable
    y_smooth = smoothdata(y_full, 'sgolay', 31);
end

%% Save extracted pixel curve
TP = table(x_full, y_smooth, 'VariableNames', {'x_pix','y_pix'});
writetable(TP, outPixCsv);

%% Pixel -> East/North
east_m  = east_min  + ((x_full - 1) ./ (W-1)) .* (east_max  - east_min);
north_m = north_max - ((y_smooth - 1) ./ (H-1)) .* (north_max - north_min);

%% Arc length and resample
ds = sqrt(diff(east_m).^2 + diff(north_m).^2);
s  = [0; cumsum(ds)];

s_uniform = (0:ds_target:s(end)).';
east_u  = interp1(s, east_m,  s_uniform, 'linear');
north_u = interp1(s, north_m, s_uniform, 'linear');

T = table((1:numel(s_uniform)).', east_u, north_u, ...
    nan(numel(s_uniform),1), s_uniform, ...
    'VariableNames', {'point_id','east_m','north_m','z_m','s_m'});
writetable(T, outCsv);

%% Overlay preview
figure('Color','w');
imshow(I); hold on;
plot(x_full, y_smooth, '-', 'Color', [0 1 1], 'LineWidth', 2.0);
title('Figure 5 MPC green line: image-based extraction overlay');
legend('Extracted MPC green line', 'Location', 'best');
exportgraphics(gca, outPng, 'Resolution', 180);

fprintf('Wrote: %s\n', outCsv);
fprintf('Wrote: %s\n', outPixCsv);
fprintf('Wrote: %s\n', outPng);
