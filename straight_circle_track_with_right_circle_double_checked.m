clear
close all
clc

%% ========================= USER INPUTS ==================================
L_straight      = 80;    % [m] initial straight length
R_circle        = 50;    % [m] radius of both circles
nLapsPerCircle  = 3;     % [-] repeated full laps on each circle; must be integer
                         %     for a continuous left-circle -> right-circle switch

ds_target       = 0.5;   % [m] target centerline arc-length spacing
trackWidth      = 8.0;   % [m] total track width, equal left/right widths

%% ========================= BUILD TRACK ==================================
[s, kappa, x, y, psi, xLeft, yLeft, xRight, yRight, totalLength, info] = ...
    build_straight_two_tangent_circles_track(L_straight, R_circle, ...
                                             nLapsPerCircle, ds_target, ...
                                             trackWidth);

%% ========================= PRINT INFO ===================================
fprintf('Straight + left repeated circle + right repeated circle track\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Total travelled distance      = %.12f m\n', totalLength);
fprintf('Analytical total distance     = %.12f m\n', ...
        L_straight + 2*nLapsPerCircle*(2*pi*R_circle));
fprintf('Number of sample points       = %d\n', numel(s));
fprintf('Target ds                     = %.12f m\n', ds_target);
fprintf('Actual straight ds            = %.12f m\n', info.ds_straight);
fprintf('Actual circle ds              = %.12f m\n', info.ds_circle);
if numel(s) >= 2
    fprintf('Min/Max spacing in s          = %.12f / %.12f m\n', ...
            min(diff(s)), max(diff(s)));
else
    fprintf('Min/Max spacing in s          = not available: only one sample point\n');
end
fprintf('Left circle center            = [%.12f, %.12f] m\n', ...
        info.leftCircleCenter(1), info.leftCircleCenter(2));
fprintf('Right circle center           = [%.12f, %.12f] m\n', ...
        info.rightCircleCenter(1), info.rightCircleCenter(2));
fprintf('Common tangent point          = [%.12f, %.12f] m\n', ...
        info.commonTangentPoint(1), info.commonTangentPoint(2));
fprintf('Left full-lap overlap error   = %.3e m\n', info.leftLapOverlapError);
fprintf('Right full-lap overlap error  = %.3e m\n', info.rightLapOverlapError);

%% ========================= OPTIONAL PLOTS ===============================
figure('Color','w');
plot(x, y, 'k-', 'LineWidth', 1.8); hold on
plot(xLeft,  yLeft,  'r-', 'LineWidth', 1.0);
plot(xRight, yRight, 'b-', 'LineWidth', 1.0);
plot(info.leftCircleCenter(1),  info.leftCircleCenter(2),  'ko', ...
    'MarkerFaceColor','y', 'MarkerSize', 6);
plot(info.rightCircleCenter(1), info.rightCircleCenter(2), 'ko', ...
    'MarkerFaceColor','c', 'MarkerSize', 6);
plot(info.commonTangentPoint(1), info.commonTangentPoint(2), 'ks', ...
    'MarkerFaceColor','m', 'MarkerSize', 6);
axis equal
grid on
xlabel('x [m] (East)')
ylabel('y [m] (North)')
legend('Centerline', 'Left boundary', 'Right boundary', ...
       'Left circle center', 'Right circle center', ...
       'Common tangent point', 'Location', 'bestoutside')
title('Straight + tangent left and right repeated circles')

figure('Color','w');
subplot(3,1,1)
plot(s, kappa, 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('\kappa [1/m]')
title('Curvature: left circle positive, right circle negative')

subplot(3,1,2)
plot(s, psi, 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('\psi [rad]')
title('Heading in [-\pi, \pi]')

subplot(3,1,3)
plot(s(2:end), diff(s), 'LineWidth', 1.4)
grid on
xlabel('s [m]')
ylabel('\Deltas [m]')
title('Arc-length spacing along centerline')

figure('Color','w');
plot(s, x, 'LineWidth', 1.2); hold on
plot(s, y, 'LineWidth', 1.2);
grid on
xlabel('s [m]')
ylabel('position [m]')
legend('x(s)', 'y(s)', 'Location', 'best')
title('Centerline coordinates versus arc length')

%% ========================= LOCAL FUNCTIONS ===============================
function [s, kappa, x, y, psi, xLeft, yLeft, xRight, yRight, totalLength, info] = ...
    build_straight_two_tangent_circles_track(L_straight, R_circle, ...
                                             nLapsPerCircle, ds_target, ...
                                             trackWidth)

    % Geometry convention:
    %   - The initial straight is the vertical line x = R_circle.
    %   - The vehicle travels upward/north along the straight.
    %   - The original/left circle has center [0, 0].
    %   - The new/right circle has center [2*R_circle, 0].
    %   - Both circles have radius R_circle.
    %   - Both circles and the straight are tangent at [R_circle, 0].
    %   - The left circle is traversed counter-clockwise, kappa = +1/R.
    %   - The right circle is traversed clockwise,       kappa = -1/R.
    %
    % The clockwise direction on the right circle is not optional if the
    % path must be tangent to the incoming straight with heading +pi/2.

    % ------------------------- checks ------------------------------------
    validateattributes(L_straight,     {'numeric'}, {'real','scalar','nonnegative'});
    validateattributes(R_circle,       {'numeric'}, {'real','scalar','positive'});
    validateattributes(nLapsPerCircle, {'numeric'}, {'real','scalar','nonnegative'});
    validateattributes(ds_target,      {'numeric'}, {'real','scalar','positive'});
    validateattributes(trackWidth,     {'numeric'}, {'real','scalar','nonnegative'});

    nLapsRounded = round(nLapsPerCircle);
    if abs(nLapsPerCircle - nLapsRounded) > 1e-12
        error(['nLapsPerCircle must be an integer for this continuous two-circle track. ', ...
               'Reason: the path can switch from the left circle to the right circle ', ...
               'only at the common tangent point [R_circle, 0].']);
    end
    nLapsPerCircle = nLapsRounded;

    halfWidth   = trackWidth / 2;
    if halfWidth >= R_circle
        warning(['trackWidth/2 is greater than or equal to R_circle. ', ...
                 'The boundary curves are still computed, but the circular ', ...
                 'inner boundary is geometrically degenerate or inverted.']);
    end

    circLength  = 2*pi*R_circle;
    totalLength = L_straight + 2*nLapsPerCircle*circLength;

    leftCenter  = [0, 0];
    rightCenter = [2*R_circle, 0];
    tangentPt   = [R_circle, 0];

    % ------------------------- straight ----------------------------------
    % Choose the number of straight segments so the spacing is close to
    % ds_target, while still ending exactly at the common tangent point.
    if L_straight > 0
        nSeg_straight = max(1, round(L_straight / ds_target));
        ds_straight   = L_straight / nSeg_straight;
        s_straight    = (0:nSeg_straight).' * ds_straight;

        x_straight    = R_circle * ones(size(s_straight));
        y_straight    = -L_straight + s_straight;
        psi_straight  = (pi/2) * ones(size(s_straight));
        kap_straight  = zeros(size(s_straight));
        seg_straight  = zeros(size(s_straight));
    else
        nSeg_straight = 0;
        ds_straight   = NaN;
        s_straight    = 0;

        x_straight    = R_circle;
        y_straight    = 0;
        psi_straight  = pi/2;
        kap_straight  = 0;
        seg_straight  = 0;
    end

    % ------------------------- circle base grid ---------------------------
    % Use one fixed angular grid for one full lap. Reuse exactly the same
    % local grid for every full lap, so repeated laps overlap point-for-point.
    if nLapsPerCircle > 0
        nSeg_perLap = max(8, round(circLength / ds_target));
        dtheta      = 2*pi / nSeg_perLap;
        ds_circle   = R_circle * dtheta;

        theta_base = (0:nSeg_perLap).' * dtheta;
        theta_keep = theta_base(2:end);  % remove duplicated lap-start point

        % Left circle, counter-clockwise, center [0, 0].
        x_left_base   = R_circle * cos(theta_keep);
        y_left_base   = R_circle * sin(theta_keep);
        psi_left_base = theta_keep + pi/2;
        kap_left_base = (1/R_circle) * ones(size(theta_keep));
        seg_left_base = ones(size(theta_keep));

        % Right circle, clockwise, center [2R, 0]. Let lambda = theta_keep.
        % At lambda = 0, the point is [R, 0] and the heading is +pi/2.
        lambda_keep    = theta_keep;
        x_right_base   = 2*R_circle - R_circle * cos(lambda_keep);
        y_right_base   = R_circle * sin(lambda_keep);
        psi_right_base = pi/2 - lambda_keep;
        kap_right_base = -(1/R_circle) * ones(size(lambda_keep));
        seg_right_base = -ones(size(lambda_keep));

        % Build repeated left-circle laps.
        s_left  = zeros(nLapsPerCircle*numel(theta_keep), 1);
        x_left  = zeros(size(s_left));
        y_left  = zeros(size(s_left));
        psi_l   = zeros(size(s_left));
        kap_l   = zeros(size(s_left));
        seg_l   = zeros(size(s_left));

        for j = 0:nLapsPerCircle-1
            idx = j*numel(theta_keep) + (1:numel(theta_keep));
            s_left(idx) = L_straight + j*circLength + theta_keep*R_circle;
            x_left(idx) = x_left_base;
            y_left(idx) = y_left_base;
            psi_l(idx)  = psi_left_base;
            kap_l(idx)  = kap_left_base;
            seg_l(idx)  = seg_left_base;
        end

        % Build repeated right-circle laps. They begin after all left laps.
        s_right  = zeros(nLapsPerCircle*numel(theta_keep), 1);
        x_right  = zeros(size(s_right));
        y_right  = zeros(size(s_right));
        psi_r    = zeros(size(s_right));
        kap_r    = zeros(size(s_right));
        seg_r    = zeros(size(s_right));

        s_right_offset = L_straight + nLapsPerCircle*circLength;
        for j = 0:nLapsPerCircle-1
            idx = j*numel(theta_keep) + (1:numel(theta_keep));
            s_right(idx) = s_right_offset + j*circLength + theta_keep*R_circle;
            x_right(idx) = x_right_base;
            y_right(idx) = y_right_base;
            psi_r(idx)   = psi_right_base;
            kap_r(idx)   = kap_right_base;
            seg_r(idx)   = seg_right_base;
        end
    else
        nSeg_perLap = 0;
        ds_circle   = NaN;

        s_left  = [];
        x_left  = [];
        y_left  = [];
        psi_l   = [];
        kap_l   = [];
        seg_l   = [];

        s_right = [];
        x_right = [];
        y_right = [];
        psi_r   = [];
        kap_r   = [];
        seg_r   = [];
    end

    % ------------------------- combine -----------------------------------
    s       = [s_straight;   s_left;  s_right];
    x       = [x_straight;   x_left;  x_right];
    y       = [y_straight;   y_left;  y_right];
    psi_raw = [psi_straight; psi_l;   psi_r];
    kappa   = [kap_straight; kap_l;   kap_r];
    segment = [seg_straight; seg_l;   seg_r];

    % Wrap heading to [-pi, pi].
    psi = atan2(sin(psi_raw), cos(psi_raw));

    % ------------------------- boundaries --------------------------------
    % Left normal for increasing s: [-sin(psi), cos(psi)].
    nx = -sin(psi);
    ny =  cos(psi);

    xLeft  = x + halfWidth * nx;
    yLeft  = y + halfWidth * ny;
    xRight = x - halfWidth * nx;
    yRight = y - halfWidth * ny;

    % ------------------------- self-checks --------------------------------
    if any(diff(s) <= 0)
        error('Internal error: s is not strictly increasing.');
    end

    if numel(s) >= 2
        if abs(s(end) - totalLength) > 100*eps(max(1, totalLength))
            error('Internal error: final s does not match totalLength.');
        end
    end

    % Check tangency and exact geometry at the common point.
    if abs(norm(tangentPt - leftCenter)  - R_circle) > 1e-12 || ...
       abs(norm(tangentPt - rightCenter) - R_circle) > 1e-12 || ...
       abs(norm(rightCenter - leftCenter) - 2*R_circle) > 1e-12
        error('Internal error: circle centers/radii are not mutually tangent as intended.');
    end

    if abs(x_straight(end) - tangentPt(1)) > 1e-12 || ...
       abs(y_straight(end) - tangentPt(2)) > 1e-12
        error('Internal error: straight does not end at the common tangent point.');
    end

    if nLapsPerCircle > 0
        leftEndIdx = numel(s_straight) + nLapsPerCircle*nSeg_perLap;
        if hypot(x(leftEndIdx) - tangentPt(1), y(leftEndIdx) - tangentPt(2)) > 1e-10
            error('Internal error: left-circle section does not end at the common tangent point.');
        end

        if hypot(x(end) - tangentPt(1), y(end) - tangentPt(2)) > 1e-10
            error('Internal error: right-circle section does not end at the common tangent point.');
        end
    end

    if any(~isfinite(s)) || any(~isfinite(x)) || any(~isfinite(y)) || ...
       any(~isfinite(psi)) || any(~isfinite(kappa)) || ...
       any(~isfinite(xLeft)) || any(~isfinite(yLeft)) || ...
       any(~isfinite(xRight)) || any(~isfinite(yRight))
        error('Internal error: generated track contains NaN or Inf values.');
    end

    % ------------------------- info --------------------------------------
    info = struct();
    info.ds_target          = ds_target;
    info.ds_straight        = ds_straight;
    info.ds_circle          = ds_circle;
    info.nSeg_straight      = nSeg_straight;
    info.nSeg_perLap        = nSeg_perLap;
    info.nLapsPerCircle     = nLapsPerCircle;
    info.leftCircleCenter   = leftCenter;
    info.rightCircleCenter  = rightCenter;
    info.commonTangentPoint = tangentPt;
    info.segment            = segment;
    info.psi_unwrapped      = psi_raw;
    info.headingNote        = 'Output psi is wrapped to [-pi, pi]. Use info.psi_unwrapped if interpolating heading directly versus s.';
    info.segmentMeaning     = '0 = straight, +1 = left CCW circle, -1 = right CW circle';

    if nLapsPerCircle >= 2
        nPer = nSeg_perLap;
        i1L = numel(s_straight) + (1:nPer);
        i2L = numel(s_straight) + nPer + (1:nPer);
        info.leftLapOverlapError = max(hypot(x(i1L) - x(i2L), y(i1L) - y(i2L)));

        i1R = numel(s_straight) + nLapsPerCircle*nPer + (1:nPer);
        i2R = numel(s_straight) + nLapsPerCircle*nPer + nPer + (1:nPer);
        info.rightLapOverlapError = max(hypot(x(i1R) - x(i2R), y(i1R) - y(i2R)));
    else
        info.leftLapOverlapError  = 0;
        info.rightLapOverlapError = 0;
    end
end
