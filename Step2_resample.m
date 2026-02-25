clear 
close all
clc

% % function resample_shanghai_boundaries()
    % Target spacing (meters)
    ds = 1.0;

    % ---- Resample both sides (returns resampled + raw for plotting) ----
    [L_new, L_raw] = resample_boundary_csv('edgeloop00_toprocess.csv', ...
    'edgeloop00_resampled.csv', ds);

    [R_new, R_raw] = resample_boundary_csv('edgeloop01_toprocess.csv', ...
    'edgeloop01_resampled.csv', ds);

    % ---- 3D plot (Z kept; spacing computed in XY unless you switch to 3D) ----
    figure('Name','Shanghai F1 Track Boundaries (Resampled vs Raw)','Color','w');
    hold on; grid on;

    % raw points (light, small)
    plot3(L_raw(:,1), L_raw(:,2), L_raw(:,3), '.', 'Color',[0.7 0.7 1.0], 'MarkerSize',6);
    plot3(R_raw(:,1), R_raw(:,2), R_raw(:,3), '.', 'Color',[1.0 0.7 0.7], 'MarkerSize',6);

    % resampled curves (bold)
    plot3(L_new(:,1), L_new(:,2), L_new(:,3), '-', 'LineWidth',2.2, 'Color',[0.0 0.2 1.0]);
    plot3(R_new(:,1), R_new(:,2), R_new(:,3), '-', 'LineWidth',2.2, 'Color',[1.0 0.1 0.1]);

    % endpoints (pinned)
    plot3(L_new(1,1), L_new(1,2), L_new(1,3), 'o', 'MarkerSize',8, 'MarkerFaceColor',[0 0.2 1], 'MarkerEdgeColor','k');
    plot3(L_new(end,1), L_new(end,2), L_new(end,3), 's', 'MarkerSize',8, 'MarkerFaceColor',[0 0.2 1], 'MarkerEdgeColor','k');
    plot3(R_new(1,1), R_new(1,2), R_new(1,3), 'o', 'MarkerSize',8, 'MarkerFaceColor',[1 0.1 0.1], 'MarkerEdgeColor','k');
    plot3(R_new(end,1), R_new(end,2), R_new(end,3), 's', 'MarkerSize',8, 'MarkerFaceColor',[1 0.1 0.1], 'MarkerEdgeColor','k');

    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title(sprintf('Resampled at ds = %.3f m (endpoints preserved)', ds));
    axis equal; view(3); box on; legend( ...
        'Left raw','Right raw','Left resampled','Right resampled', ...
        'Left start','Left finish','Right start','Right finish', ...
        'Location','bestoutside');

    % Optional: also show top view
    figure('Name','Top View (XY)','Color','w');
    hold on; grid on; axis equal; box on;
    plot(L_new(:,1), L_new(:,2), '-o', 'LineWidth',2.0, 'Color',[0.0 0.2 1.0]);
    plot(R_new(:,1), R_new(:,2), '-o', 'LineWidth',2.0, 'Color',[1.0 0.1 0.1]);
    xlabel('X (m)'); ylabel('Y (m)');
    title('Top View of Resampled Boundaries');
% end

function [XYZ_new, XYZ_raw] = resample_boundary_csv(infile, outfile, ds)
    % ---- Read CSV (x,y,z). Works with/without header row. ----
    T = readmatrix(infile);
    if size(T,2) < 3
        error('Input "%s" must have ≥3 columns (x,y,z).', infile);
    end
    x = T(:,1); y = T(:,2); z = T(:,3);
    XYZ_raw = [x y z];

    % ---- Clean: drop NaNs and collapse consecutive duplicates ----
    good = ~(isnan(x) | isnan(y) | isnan(z));
    x = x(good); y = y(good); z = z(good);
    if numel(x) < 2, error('"%s" has <2 valid points after cleaning.', infile); end
    keep = [true; (diff(x)~=0 | diff(y)~=0 | diff(z)~=0)];
    x = x(keep); y = y(keep); z = z(keep);

    % ---- Arc length (change to 3D if desired) ----
    dx = diff(x); dy = diff(y); dz = diff(z);
    seglen = hypot(dx, dy);                        % XY ground distance
    % seglen = sqrt(dx.^2 + dy.^2 + dz.^2);        % <- use for true 3D spacing

    s = [0; cumsum(seglen)];
    L = s(end);  if ~isfinite(L) || L<=0, error('Non-positive length in "%s".', infile); end

    % Ensure strictly increasing s for interp
    [s_unique, iu] = unique(s, 'stable');
    x = x(iu); y = y(iu); z = z(iu);

    % ---- Targets: include 0 and L; force exact L ----
    if ds <= 0, error('ds must be > 0'); end
    s_target = 0:ds:L;
    if isempty(s_target) || abs(s_target(end)-L) > 1e-12
        s_target(end+1) = L;
    end

    % ---- Interpolate (shape-preserving) ----
    x_new = interp1(s_unique, x, s_target, 'pchip');
    y_new = interp1(s_unique, y, s_target, 'pchip');
    z_new = interp1(s_unique, z, s_target, 'pchip');

    % ---- HARD-PIN endpoints to original values ----
    x_new(1)   = x(1);   y_new(1)   = y(1);   z_new(1)   = z(1);
    x_new(end) = x(end); y_new(end) = y(end); z_new(end) = z(end);

    XYZ_new = [x_new(:) y_new(:) z_new(:)];

    % ---- Write CSV (full precision) ----
    fid = fopen(outfile, 'w');  assert(fid~=-1, 'Cannot open "%s".', outfile);
    fprintf(fid, 'x,y,z\n');
    fprintf(fid, '%.15g,%.15g,%.15g\n', XYZ_new.');
    fclose(fid);

    fprintf('Wrote %-34s  N=%d  L=%.3f m  ds~%.3f m\n', ...
        outfile, size(XYZ_new,1), L, mean(diff(s_target)));
end
