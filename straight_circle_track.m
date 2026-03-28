clear
close all
clc

%% ========================= USER INPUTS ==================================
L_straight = 80;      % [m] straight length
R_circle   = 50;      % [m] circle radius
nLaps      = 3.0;     % [-] number of circular laps (can be non-integer)
ds_target  = 0.5;     % [m] target centerline arc-length spacing
trackWidth = 8.0;     % [m] total track width (equal left/right widths)

%% ========================= BUILD TRACK ==================================
[s, kappa, x, y, psi, xLeft, yLeft, xRight, yRight, totalLength, info] = ...
    build_straight_circle_track(L_straight, R_circle, nLaps, ds_target, trackWidth);

%% ========================= PRINT INFO ===================================
fprintf('Total travelled distance      = %.12f m\n', totalLength);
fprintf('Analytical total distance     = %.12f m\n', L_straight + 2*pi*R_circle*nLaps);
fprintf('Number of sample points       = %d\n', numel(s));
fprintf('Target ds                     = %.12f m\n', ds_target);
fprintf('Actual straight ds            = %.12f m\n', info.ds_straight);
fprintf('Actual full-circle ds         = %.12f m\n', info.ds_circle);
fprintf('Min/Max spacing in s          = %.12f / %.12f m\n', min(diff(s)), max(diff(s)));

%% ========================= OPTIONAL PLOTS ===============================
figure('Color','w');
plot(x, y, 'k-o', 'LineWidth', 1.5); hold on
plot(xLeft,  yLeft,  'r-o', 'LineWidth', 1.0);
plot(xRight, yRight, 'b-o', 'LineWidth', 1.0);
plot(0,0,'ko','MarkerFaceColor','y','MarkerSize',6)
axis equal
grid on
xlabel('x [m] (East)')
ylabel('y [m] (North)')
legend('Centerline','Left boundary','Right boundary','Circle center','Location','best')
title('Straight + repeated circle track')

figure('Color','w');
subplot(3,1,1)
plot(s, kappa, 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('\kappa [1/m]')
title('Curvature')

subplot(3,1,2)
plot(s, psi, 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('\psi [rad]')
title('Heading in [-\pi,\pi]')

subplot(3,1,3)
plot(s(2:end), diff(s), 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('\Deltas [m]')
title('Arc-length spacing along centerline')

%% ========================= LOCAL FUNCTION ===============================
function [s, kappa, x, y, psi, xLeft, yLeft, xRight, yRight, totalLength, info] = ...
    build_straight_circle_track(L_straight, R_circle, nLaps, ds_target, trackWidth)

    % ------------------------- checks ------------------------------------
    validateattributes(L_straight, {'numeric'}, {'real','scalar','nonnegative'});
    validateattributes(R_circle,   {'numeric'}, {'real','scalar','positive'});
    validateattributes(nLaps,      {'numeric'}, {'real','scalar','nonnegative'});
    validateattributes(ds_target,  {'numeric'}, {'real','scalar','positive'});
    validateattributes(trackWidth, {'numeric'}, {'real','scalar','nonnegative'});

    halfWidth   = trackWidth / 2;
    circLength  = 2*pi*R_circle;
    totalLength = L_straight + nLaps*circLength;

    % ------------------------- straight ----------------------------------
    % Choose number of straight segments so spacing is close to ds_target.
    if L_straight > 0
        nSeg_straight = max(1, round(L_straight / ds_target));
        ds_straight   = L_straight / nSeg_straight;
        s_straight    = (0:nSeg_straight).' * ds_straight;

        x_straight    = R_circle * ones(size(s_straight));
        y_straight    = -L_straight + s_straight;
        psi_straight  = (pi/2) * ones(size(s_straight));
        kap_straight  = zeros(size(s_straight));
    else
        nSeg_straight = 0;
        ds_straight   = NaN;
        s_straight    = 0;

        x_straight    = R_circle;
        y_straight    = 0;
        psi_straight  = pi/2;
        kap_straight  = 0;
    end

    % ------------------------- circle ------------------------------------
    % Choose one fixed angular grid for a full lap. Reuse it for every full
    % lap so all full laps overlap exactly.
    if nLaps > 0
        nSeg_perLap = max(3, round(circLength / ds_target));
        dtheta      = 2*pi / nSeg_perLap;
        ds_circle   = R_circle * dtheta;

        nFull = floor(nLaps);
        frac  = nLaps - nFull;

        s_circle   = [];
        theta_all  = [];

        % Full laps: identical theta grid every lap
        for j = 0:nFull-1
            theta_j = (0:nSeg_perLap).' * dtheta + 2*pi*j;
            s_j     = L_straight + j*circLength + (0:nSeg_perLap).' * ds_circle;

            % remove duplicated lap start point
            theta_j = theta_j(2:end);
            s_j     = s_j(2:end);

            theta_all = [theta_all; theta_j]; %#ok<AGROW>
            s_circle  = [s_circle;  s_j];     %#ok<AGROW>
        end

        % Fractional last lap: exact end point, spacing still close to ds_target
        if frac > 0
            theta_end = 2*pi*frac;
            nSeg_frac = max(1, round((R_circle*theta_end) / ds_target));

            theta_f = linspace(0, theta_end, nSeg_frac+1).' + 2*pi*nFull;
            s_f     = linspace(L_straight + nFull*circLength, totalLength, nSeg_frac+1).';

            % remove duplicated start point
            theta_f = theta_f(2:end);
            s_f     = s_f(2:end);

            theta_all = [theta_all; theta_f]; %#ok<AGROW>
            s_circle  = [s_circle;  s_f];     %#ok<AGROW>
        end

        x_circle   = R_circle * cos(theta_all);
        y_circle   = R_circle * sin(theta_all);
        psi_circle = theta_all + pi/2;
        kap_circle = (1/R_circle) * ones(size(theta_all));
    else
        nSeg_perLap = 0;
        ds_circle   = NaN;

        s_circle    = [];
        x_circle    = [];
        y_circle    = [];
        psi_circle  = [];
        kap_circle  = [];
    end

    % ------------------------- combine -----------------------------------
    s        = [s_straight;   s_circle];
    x        = [x_straight;   x_circle];
    y        = [y_straight;   y_circle];
    psi_raw  = [psi_straight; psi_circle];
    kappa    = [kap_straight; kap_circle];

    % Wrap heading to [-pi, pi]
    psi = atan2(sin(psi_raw), cos(psi_raw));

    % ------------------------- boundaries --------------------------------
    % Left normal for increasing s: [-sin(psi), cos(psi)]
    nx = -sin(psi);
    ny =  cos(psi);

    xLeft  = x + halfWidth * nx;
    yLeft  = y + halfWidth * ny;
    xRight = x - halfWidth * nx;
    yRight = y - halfWidth * ny;

    % ------------------------- info --------------------------------------
    info = struct();
    info.ds_target      = ds_target;
    info.ds_straight    = ds_straight;
    info.ds_circle      = ds_circle;
    info.nSeg_straight  = nSeg_straight;
    info.nSeg_perLap    = nSeg_perLap;
    info.circleCenter   = [0, 0];
end



% clear
% close all
% clc

%% ========================= USER INPUTS ==================================
vehicleWidth     = 1.85;   % [m] overall vehicle width without mirrors
ds_target        = 0.20;   % [m] target spacing along smooth reference path
coneBaseDiameter = 0.185;  % [m] cone base-circle diameter
nConesPerRow     = 5;      % figure-consistent cone rows
nCirclePts       = 120;    % drawing resolution for cone base circles
tolInside        = 1e-9;   % numerical tolerance for containment check

%% ========================= BUILD EXACT CORRIDOR =========================
track = build_iso3888_2_corridor(vehicleWidth, coneBaseDiameter, ...
                                 nConesPerRow, nCirclePts);

%% ========================= BUILD SMOOTH REFERENCE PATH ==================
ref = build_smooth_reference_path(track, ds_target);

%% ========================= VERIFY CONTAINMENT ===========================
report = verify_reference_against_corridor(track, ref, tolInside);

%% ========================= PRINT RESULTS ================================
fprintf('ISO 3888-2 corridor + smooth reference path\n');
fprintf('-------------------------------------------\n');
fprintf('Longitudinal corridor length     = %.12f m\n', track.X5);
fprintf('Smooth centerline length         = %.12f m\n', ref.totalDistance);
fprintf('Target ds                        = %.12f m\n', ds_target);
fprintf('Actual ds                        = %.12f m\n', ref.ds_actual);
fprintf('b1, b3, b5                       = [%.12f  %.12f  %.12f] m\n', ...
        track.b1, track.b3, track.b5);
fprintf('\n');
fprintf('Centerline inside corridor?      %s\n', pass_fail(report.centerline.pass));
fprintf('Left ref boundary inside?        %s\n', pass_fail(report.left.pass));
fprintf('Right ref boundary inside?       %s\n', pass_fail(report.right.pass));
fprintf('\n');
fprintf('Worst centerline upper margin    = %.12e m\n', report.centerline.minUpperMargin);
fprintf('Worst centerline lower margin    = %.12e m\n', report.centerline.minLowerMargin);
fprintf('Worst left upper margin          = %.12e m\n', report.left.minUpperMargin);
fprintf('Worst left lower margin          = %.12e m\n', report.left.minLowerMargin);
fprintf('Worst right upper margin         = %.12e m\n', report.right.minUpperMargin);
fprintf('Worst right lower margin         = %.12e m\n', report.right.minLowerMargin);

if report.allPass
    fprintf('\nVERIFICATION RESULT: PASS\n');
else
    fprintf('\nVERIFICATION RESULT: FAIL\n');

    if ~report.centerline.pass
        print_fail_info('CENTERLINE', ref.x, ref.y, report.centerline);
    end
    if ~report.left.pass
        print_fail_info('LEFT REFERENCE BOUNDARY', ref.xLeft, ref.yLeft, report.left);
    end
    if ~report.right.pass
        print_fail_info('RIGHT REFERENCE BOUNDARY', ref.xRight, ref.yRight, report.right);
    end
end

%% ========================= OPTIONAL PLOT ================================
figure('Color','w'); hold on

patch(track.polyX, track.polyY, [0.92 0.92 0.92], ...
    'EdgeColor', 'k', 'LineWidth', 1.2, 'DisplayName', 'Exact corridor');

plot(track.upperBoundary(:,1), track.upperBoundary(:,2), 'k-', ...
    'LineWidth', 1.8, 'DisplayName', 'Exact corridor boundary');
plot(track.lowerBoundary(:,1), track.lowerBoundary(:,2), 'k-', ...
    'LineWidth', 1.8, 'HandleVisibility', 'off');

for i = 1:size(track.coneCenters,1)
    xc = track.coneCenters(i,1);
    yc = track.coneCenters(i,2);
    plot(xc + track.coneBaseRadius*cos(track.thetaCircle), ...
         yc + track.coneBaseRadius*sin(track.thetaCircle), ...
         'b-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
plot(track.coneCenters(:,1), track.coneCenters(:,2), 'b.', ...
    'MarkerSize', 10, 'DisplayName', 'Cone centers');

plot(ref.x, ref.y, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Smooth centerline');
plot(ref.xLeft, ref.yLeft, '--', 'Color', [0.85 0.2 0.2], ...
    'LineWidth', 1.0, 'DisplayName', 'Smooth left reference boundary');
plot(ref.xRight, ref.yRight, '--', 'Color', [0.2 0.6 0.2], ...
    'LineWidth', 1.0, 'DisplayName', 'Smooth right reference boundary');

axis equal
grid on
xlabel('x [m] (longitudinal / East)')
ylabel('y [m] (lateral / North)')
title('Exact corridor paired with smooth reference path')
legend('Location', 'bestoutside')

figure('Color','w');
subplot(3,1,1)
plot(ref.s, ref.kappa, 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('\kappa [1/m]')
title('Smooth centerline curvature')

subplot(3,1,2)
plot(ref.s, ref.psi, 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('\psi [rad]')
title('Smooth centerline heading')

subplot(3,1,3)
plot(ref.s, ref.width, 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('width [m]')
title('Smooth reference width profile')

%% ========================= LOCAL FUNCTIONS ==============================
function track = build_iso3888_2_corridor(vehicleWidth, coneBaseDiameter, ...
                                          nConesPerRow, nCirclePts)

    validateattributes(vehicleWidth,     {'numeric'}, {'real','scalar','positive'});
    validateattributes(coneBaseDiameter, {'numeric'}, {'real','scalar','positive'});
    validateattributes(nConesPerRow,     {'numeric'}, {'real','scalar','integer','>=',2});
    validateattributes(nCirclePts,       {'numeric'}, {'real','scalar','integer','>=',20});

    % ISO 3888-2:2011 dimensions
    track.L1 = 12.0;
    track.L2 = 13.5;
    track.L3 = 11.0;
    track.L4 = 12.5;
    track.L5 = 12.0;
    track.laneOffset = 1.0;

    track.b1 = 1.1*vehicleWidth + 0.25;
    track.b3 = vehicleWidth + 1.0;
    track.b5 = max(1.3*vehicleWidth + 0.25, 3.0);

    track.coneBaseDiameter = coneBaseDiameter;
    track.coneBaseRadius   = coneBaseDiameter/2.0;
    track.nConesPerRow     = nConesPerRow;
    track.nCirclePts       = nCirclePts;

    track.X0 = 0.0;
    track.X1 = track.L1;
    track.X2 = track.L1 + track.L2;
    track.X3 = track.L1 + track.L2 + track.L3;
    track.X4 = track.L1 + track.L2 + track.L3 + track.L4;
    track.X5 = track.L1 + track.L2 + track.L3 + track.L4 + track.L5;

    track.y1U =  0.5*track.b1;
    track.y1L = -0.5*track.b1;

    track.y3U = track.laneOffset + 0.5*track.b3;
    track.y3L = track.laneOffset - 0.5*track.b3;

    track.y5U =  0.5*track.b5;
    track.y5L = -0.5*track.b5;

    poly = [ ...
        track.X0, track.y1L;
        track.X2, track.y1L;
        track.X2, track.y3L;
        track.X3, track.y3L;
        track.X3, track.y5L;
        track.X5, track.y5L;
        track.X5, track.y5U;
        track.X4, track.y5U;
        track.X4, track.y3U;
        track.X1, track.y3U;
        track.X1, track.y1U;
        track.X0, track.y1U;
        track.X0, track.y1L];

    track.polyX = poly(:,1);
    track.polyY = poly(:,2);

    track.upperBoundary = [ ...
        track.X0, track.y1U;
        track.X1, track.y1U;
        track.X1, track.y3U;
        track.X4, track.y3U;
        track.X4, track.y5U;
        track.X5, track.y5U ];

    track.lowerBoundary = [ ...
        track.X0, track.y1L;
        track.X2, track.y1L;
        track.X2, track.y3L;
        track.X3, track.y3L;
        track.X3, track.y5L;
        track.X5, track.y5L ];

    % Figure-consistent cone rows
    xRow1 = linspace(track.X0 + track.coneBaseRadius, ...
                     track.X1 - track.coneBaseRadius, nConesPerRow);
    xRow3 = linspace(track.X2 + track.coneBaseRadius, ...
                     track.X3 - track.coneBaseRadius, nConesPerRow);
    xRow5 = linspace(track.X4 + track.coneBaseRadius, ...
                     track.X5 - track.coneBaseRadius, nConesPerRow);

    yRow1Top = (track.y1U + track.coneBaseRadius) * ones(size(xRow1));
    yRow1Bot = (track.y1L - track.coneBaseRadius) * ones(size(xRow1));

    yRow3Top = (track.y3U + track.coneBaseRadius) * ones(size(xRow3));
    yRow3Bot = (track.y3L - track.coneBaseRadius) * ones(size(xRow3));

    yRow5Top = (track.y5U + track.coneBaseRadius) * ones(size(xRow5));
    yRow5Bot = (track.y5L - track.coneBaseRadius) * ones(size(xRow5));

    track.coneCenters = [ ...
        xRow1(:), yRow1Top(:);
        xRow1(:), yRow1Bot(:);
        xRow3(:), yRow3Top(:);
        xRow3(:), yRow3Bot(:);
        xRow5(:), yRow5Top(:);
        xRow5(:), yRow5Bot(:)];

    track.thetaCircle = linspace(0, 2*pi, nCirclePts);
end

function ref = build_smooth_reference_path(track, ds_target)

    validateattributes(ds_target, {'numeric'}, {'real','scalar','positive'});

    dx_dense = min(0.005, max(ds_target/50, 1e-4));
    X_dense = (track.X0:dx_dense:track.X5).';
    if abs(X_dense(end) - track.X5) > 1e-14
        X_dense = [X_dense; track.X5];
    end

    [~, dYdX_dense, ~, ~] = eval_smooth_reference_profile(X_dense, track);

    dS_dX = sqrt(1.0 + dYdX_dense.^2);
    S_dense = cumtrapz(X_dense, dS_dX);
    totalDistance = S_dense(end);

    nSeg = max(1, round(totalDistance / ds_target));
    s = linspace(0.0, totalDistance, nSeg + 1).';
    ds_actual = totalDistance / nSeg;

    x = interp1(S_dense, X_dense, s, 'pchip');

    [y, dYdX, d2YdX2, width] = eval_smooth_reference_profile(x, track);

    psi = atan2(dYdX, ones(size(dYdX)));
    psi = atan2(sin(psi), cos(psi));

    kappa = d2YdX2 ./ (1.0 + dYdX.^2).^(3/2);

    nx = -sin(psi);
    ny =  cos(psi);

    halfW = 0.5 * width;
    xLeft  = x + halfW .* nx;
    yLeft  = y + halfW .* ny;
    xRight = x - halfW .* nx;
    yRight = y - halfW .* ny;

    ref = struct();
    ref.s             = s;
    ref.x             = x;
    ref.y             = y;
    ref.psi           = psi;
    ref.kappa         = kappa;
    ref.width         = width;
    ref.xLeft         = xLeft;
    ref.yLeft         = yLeft;
    ref.xRight        = xRight;
    ref.yRight        = yRight;
    ref.totalDistance = totalDistance;
    ref.ds_actual     = ds_actual;
end

function [Y, dYdX, d2YdX2, W] = eval_smooth_reference_profile(X, track)

    Y       = zeros(size(X));
    dYdX    = zeros(size(X));
    d2YdX2  = zeros(size(X));
    W       = zeros(size(X));

    idx1 = (X <= track.X1);
    Y(idx1)      = 0.0;
    dYdX(idx1)   = 0.0;
    d2YdX2(idx1) = 0.0;
    W(idx1)      = track.b1;

    idx2 = (X > track.X1) & (X <= track.X2);
    if any(idx2)
        u = (X(idx2) - track.X1) / track.L2;
        [h, dh, d2h] = quintic_smoothstep(u);

        Y(idx2)      = track.laneOffset * h;
        dYdX(idx2)   = track.laneOffset * dh / track.L2;
        d2YdX2(idx2) = track.laneOffset * d2h / (track.L2^2);
        W(idx2)      = track.b1 + (track.b3 - track.b1) * h;
    end

    idx3 = (X > track.X2) & (X <= track.X3);
    Y(idx3)      = track.laneOffset;
    dYdX(idx3)   = 0.0;
    d2YdX2(idx3) = 0.0;
    W(idx3)      = track.b3;

    idx4 = (X > track.X3) & (X <= track.X4);
    if any(idx4)
        u = (X(idx4) - track.X3) / track.L4;
        [h, dh, d2h] = quintic_smoothstep(u);

        Y(idx4)      = track.laneOffset * (1.0 - h);
        dYdX(idx4)   = -track.laneOffset * dh / track.L4;
        d2YdX2(idx4) = -track.laneOffset * d2h / (track.L4^2);
        W(idx4)      = track.b3 + (track.b5 - track.b3) * h;
    end

    idx5 = (X > track.X4);
    Y(idx5)      = 0.0;
    dYdX(idx5)   = 0.0;
    d2YdX2(idx5) = 0.0;
    W(idx5)      = track.b5;
end

function report = verify_reference_against_corridor(track, ref, tol)

    report.centerline = containment_report(track, ref.x,      ref.y,      tol);
    report.left       = containment_report(track, ref.xLeft,  ref.yLeft,  tol);
    report.right      = containment_report(track, ref.xRight, ref.yRight, tol);

    report.allPass = report.centerline.pass && report.left.pass && report.right.pass;
end

function rep = containment_report(track, x, y, tol)

    x = x(:);
    y = y(:);

    upper = zeros(size(x));
    lower = zeros(size(x));

    idxU1 = (x <= track.X1);
    idxU2 = (x >  track.X1) & (x <= track.X4);
    idxU3 = (x >  track.X4);

    upper(idxU1) = track.y1U;
    upper(idxU2) = track.y3U;
    upper(idxU3) = track.y5U;

    idxL1 = (x <= track.X2);
    idxL2 = (x >  track.X2) & (x <= track.X3);
    idxL3 = (x >  track.X3);

    lower(idxL1) = track.y1L;
    lower(idxL2) = track.y3L;
    lower(idxL3) = track.y5L;

    marginUpper = upper - y;
    marginLower = y - lower;
    marginXMin  = x - track.X0;
    marginXMax  = track.X5 - x;

    insidePiecewise = (marginUpper >= -tol) & ...
                      (marginLower >= -tol) & ...
                      (marginXMin  >= -tol) & ...
                      (marginXMax  >= -tol);

    [inPoly, onPoly] = inpolygon(x, y, track.polyX, track.polyY);
    insidePolygon = inPoly | onPoly;

    rep.pass = all(insidePiecewise) && all(insidePolygon);

    rep.minUpperMargin = min(marginUpper);
    rep.minLowerMargin = min(marginLower);
    rep.minXMinMargin  = min(marginXMin);
    rep.minXMaxMargin  = min(marginXMax);

    rep.numPiecewiseFail = sum(~insidePiecewise);
    rep.numPolygonFail   = sum(~insidePolygon);

    rep.firstPiecewiseFailIndex = first_fail_index(insidePiecewise);
    rep.firstPolygonFailIndex   = first_fail_index(insidePolygon);
end

function idx = first_fail_index(mask)
    k = find(~mask, 1, 'first');
    if isempty(k)
        idx = NaN;
    else
        idx = k;
    end
end

function [h, dh, d2h] = quintic_smoothstep(u)
    u = min(max(u, 0.0), 1.0);
    h   = 10*u.^3 - 15*u.^4 + 6*u.^5;
    dh  = 30*u.^2 - 60*u.^3 + 30*u.^4;
    d2h = 60*u - 180*u.^2 + 120*u.^3;
end

function txt = pass_fail(tf)
    if tf
        txt = 'PASS';
    else
        txt = 'FAIL';
    end
end

function print_fail_info(nameStr, x, y, rep)
    fprintf('\n%s failure details:\n', nameStr);

    if ~isnan(rep.firstPiecewiseFailIndex)
        k = rep.firstPiecewiseFailIndex;
        fprintf('  First piecewise fail index = %d\n', k);
        fprintf('  Coordinates                = [%.12f, %.12f]\n', x(k), y(k));
    end

    if ~isnan(rep.firstPolygonFailIndex)
        k = rep.firstPolygonFailIndex;
        fprintf('  First polygon fail index   = %d\n', k);
        fprintf('  Coordinates                = [%.12f, %.12f]\n', x(k), y(k));
    end
end