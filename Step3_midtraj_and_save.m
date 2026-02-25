clear
close all
clc






%%
% test_mid_traj.m
% Robust centerline from two edge loops; outputs edge_loop_mid.csv and plots.

%% ------------ Settings ------------
file0 = 'edgeloop00_resampled.csv';
file1 = 'edgeloop01_resampled.csv';
outname = 'midtraj_generated.csv';

% Set to 0 to use max(#points) of the two files; otherwise force N samples.
N_common = 0;

% Continuity window (how far we search along the other edge, in segments)
search_window = 40;   % increase to 60–80 if your edges are very uneven

%% ------------ Load CSVs (2D or 3D) ------------
P0 = read_numeric_xyz(file0);
P1 = read_numeric_xyz(file1);
D  = min(size(P0,2), size(P1,2));
P0 = P0(:,1:D);
P1 = P1(:,1:D);

if size(P0,1) < 2 || size(P1,1) < 2
    error('Each CSV must have at least two numeric rows.');
end

%% ------------ Resample both by arc length to a common N ------------
if N_common <= 0
    N_common = max(size(P0,1), size(P1,1));
end
P0r = resample_by_arclength(P0, N_common);
P1r = resample_by_arclength(P1, N_common);

%% ------------ Coarse alignment (orientation + circular shift) ------------
[P1a, stats] = best_align_loop(P0r, P1r);
fprintf('Alignment: orientation=%s, shift=%d, mean pairwise distance=%.6g\n', ...
    stats.orientation, stats.shift, stats.mean_dist);

%% ------------ Map each point on P0r to nearest on P1a (continuous) ------------
% Close P1a into segments
P1c = [P1a; P1a(1,:)];
A   = P1c(1:end-1,:);     % segment starts
B   = P1c(2:end,:);       % segment ends
AB  = B - A;
AB2 = sum(AB.^2, 2);
M   = size(A,1);

Q1  = zeros(size(P0r));   % nearest points on P1a
idx = zeros(N_common,1);  % segment index used
tt  = zeros(N_common,1);  % local parameter in [0,1]

% First point: global search
[pd2, idx(1), tt(1), Q1(1,:)] = project_to_segments(P0r(1,:), A, AB, AB2, 1:M); %#ok<ASGLU>

% Subsequent points: local window around previous segment index
for i = 2:N_common
    cand = circular_window(idx(i-1), search_window, M);
    [pd2, idx(i), tt(i), Q1(i,:)] = project_to_segments(P0r(i,:), A, AB, AB2, cand); %#ok<ASGLU>
end

%% ------------ Midpoints ------------
Pm = 0.5 * (P0r + Q1);

%% ------------ Write CSV ------------
if D == 2
    T = array2table(Pm, 'VariableNames', {'x','y'});
else
    T = array2table(Pm, 'VariableNames', {'x','y','z'});
end
writetable(T, outname);
fprintf('Wrote %s (%d pts, %d dims)\n', outname, size(Pm,1), size(Pm,2));

%% ------------ Plot ------------
figure('Name','Edge loops and midpoint','Color','w');
hold on; grid on; axis equal;

if D == 2
    plot(P0r(:,1), P0r(:,2), '-', 'LineWidth', 1.5);
    plot(P1a(:,1), P1a(:,2), '-', 'LineWidth', 1.5);
    plot(Pm(:,1),  Pm(:,2),  '-', 'LineWidth', 2.0);
    xlabel('X'); ylabel('Y');
else
    plot3(P0r(:,1), P0r(:,2), P0r(:,3), '-', 'LineWidth', 1.5);
    plot3(P1a(:,1), P1a(:,2), P1a(:,3), '-', 'LineWidth', 1.5);
    plot3(Pm(:,1),  Pm(:,2),  Pm(:,3),  '-', 'LineWidth', 2.0);
    xlabel('X'); ylabel('Y'); zlabel('Z'); view(3);
end
legend({'edge\_loop\_00','edge\_loop\_01 (aligned)','midpoint'}, 'Location','best');
title('Edge loops and midpoint trajectory');

figure()
plot3(P0(:,1), P0(:,2), P0(:,3), '-o')
hold on
plot3(P1(:,1), P1(:,2), P1(:,3), '-o')
hold on
plot3(Pm(:,1),  Pm(:,2),  Pm(:,3),  '-o');
hold on
plot3(P0(1,1), P0(1,2), P0(1,3), 'or', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
hold on
plot3(P1(1,1), P1(1,2), P1(1,3), 'or', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
hold on
plot3(Pm(1,1), Pm(1,2), Pm(1,3), 'or', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
axis equal
grid minor


%% Save

cd('C:\Shanghai\Step2')

% center line
xEast_C_Blender_raw = Pm(:, 1);
yNorth_C_Blender_raw = Pm(:, 2);
zUp_C_Blender_raw = Pm(:, 3);
save('C_Cartesian_coords_Blender_raw', 'xEast_C_Blender_raw', 'yNorth_C_Blender_raw', 'zUp_C_Blender_raw')

% left boundary
xEast_CL_Blender_raw = P0(:, 1);
yNorth_CL_Blender_raw = P0(:, 2);
zUp_CL_Blender_raw = P0(:, 3);
save('CL_Cartesian_coords_Blender_raw', 'xEast_CL_Blender_raw', 'yNorth_CL_Blender_raw', 'zUp_CL_Blender_raw')

% right boundary
xEast_CR_Blender_raw = P1(:, 1);
yNorth_CR_Blender_raw = P1(:, 2);
zUp_CR_Blender_raw = P1(:, 3);
save('CR_Cartesian_coords_Blender_raw', 'xEast_CR_Blender_raw', 'yNorth_CR_Blender_raw', 'zUp_CR_Blender_raw')


%% ======================= Helpers =======================

function A = read_numeric_xyz(fname)
    % Reads CSV, keeps first up to 3 numeric columns, drops rows with NaNs.
    A = readmatrix(fname);
    if isempty(A)
        error('File %s is empty or unreadable.', fname);
    end
    A = A(all(~isnan(A),2), :);
    if size(A,2) < 2
        error('File %s must contain at least two numeric columns.', fname);
    end
    A = A(:, 1:min(3, size(A,2)));
end

function P = resample_by_arclength(P, N)
    % Remove exact/near duplicates to keep s strictly increasing
    dif = diff(P,1,1);
    seg = sqrt(sum(dif.^2, 2));
    keep = [true; seg > 1e-12];
    P = P(keep,:);
    if size(P,1) < 2
        P = repmat(P(1,:), N, 1);
        return;
    end

    % Cumulative arc length
    dif = diff(P,1,1);
    seg = sqrt(sum(dif.^2, 2));
    s   = [0; cumsum(seg)];
    L   = s(end);
    if L <= eps
        P = repmat(P(1,:), N, 1);
        return;
    end

    s_new = linspace(0, L, N).';
    P = interp1(s, P, s_new, 'linear');
end

function [P2_best, stats] = best_align_loop(P1, P2)
    % Try forward vs reversed and a small circular-shift search.
    N = size(P1,1);
    cands = {P2, flipud(P2)};
    oris  = {'forward','reversed'};
    best.mean = inf; best.P = []; best.shift = 0; best.ori = oris{1};

    for o = 1:2
        Q = cands{o};
        % initial shift by nearest to P1(1,:)
        d = sqrt(sum((Q - P1(1,:)).^2, 2));
        [~,k0] = min(d);
        win = unique(mod((k0-5):(k0+5) - 1, N) + 1);
        for k = win
            Qs = circshift(Q, -(k-1), 1);
            m  = mean(sqrt(sum((P1 - Qs).^2, 2)));
            if m < best.mean
                best.mean = m; best.P = Qs; best.shift = k-1; best.ori = oris{o};
            end
        end
    end

    P2_best = best.P;
    stats = struct('mean_dist', best.mean, 'shift', best.shift, 'orientation', best.ori);
end

function [dmin, idx_best, t_best, q_best] = project_to_segments(p, A, AB, AB2, seg_idx)
    % Project point p onto chosen segments of a closed polyline.
    Ai   = A(seg_idx,:);
    ABi  = AB(seg_idx,:);
    AB2i = AB2(seg_idx,:);
    AB2i(AB2i < eps) = eps;          % protect zero-length segments

    AP   = bsxfun(@minus, p, Ai);    % p - Ai
    ti   = sum(AP .* ABi, 2) ./ AB2i;
    ti   = max(0, min(1, ti));       % clamp to [0,1]
    qi   = Ai + bsxfun(@times, ti, ABi);
    di2  = sum((qi - p).^2, 2);

    [dmin, k] = min(di2);
    idx_best  = seg_idx(k);
    t_best    = ti(k);
    q_best    = qi(k,:);
end

function idxs = circular_window(center, W, M)
    % Return a wrapped index window of segments around 'center'.
    if center <= 0, center = 1; end
    rng  = (center-W):(center+W);
    idxs = mod(rng-1, M) + 1;
    idxs = unique(idxs);
end





