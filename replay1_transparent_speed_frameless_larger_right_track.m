% Replay planar-track / planar-vehicle data with four wheels.

% Output: ONE combined video

% - main 3D chase-camera view

% - small 2D top-down inset at upper-left

% - upper-middle HUD:

% * virtual steering wheel driven by deltaFL

% * drive torque bar driven by Tdrive_norm

% * brake torque bar driven by Tbrake_norm

% - top-right speed mini-plot:

% * animated speed curve

% * moving current-point marker

% * current speed value shown beside the marker

%

% REQUIRED VARIABLES IN WORKSPACE:

% t

% carX, carY, carZ

% speed

% heading, roll, pitch

% leftX, leftY, leftZ

% rightX, rightY, rightZ

% deltaFL, deltaFR

% omegaFL, omegaFR, omegaRL, omegaRR

% Tdrive_norm, Tbrake_norm

%

% USER CONVENTIONS:

% Global frame:

% X = East

% Y = North

% Z = Up

%

% Heading:

% positive CCW from +X (East), in [-pi, pi]

%

% Roll:

% positive when car rolls to the right

%

% Pitch:

% positive when car front pitches down

%

% BODY FRAME USED HERE:

% +x forward

% +y left

% +z up

%

% ROTATION USED:

% R = Rz(heading) * Ry(pitch) * Rx(roll)

%

% Steering convention assumed here:

% positive delta = left turn = CCW about body +z

%

% Wheel angular rate convention:

% omega in rad/s

% If the visual spin direction looks opposite to what you want,

% flip wheelSpinSign below.

% close all; clear; clc;

% cd('D:\App\my_oc_test\oc_main\oc')
cd('D:\Applications\TOTPT-main\oc')

veh = my_params();

% cd('D:\App\my_oc_test\oc_main\oc\my_replay')
cd('D:\Applications\TOTPT-main\oc\my_replay')

load data4replay.mat

% load dataofpowertrain.mat

%% ========================= USER SETTINGS =========================

% ----- angle units -----

ANGLES_IN_DEG = false; % heading/roll/pitch in deg?

STEER_ANGLES_IN_DEG = false; % deltaFL/deltaFR in deg?

% ----- video settings -----

videoFPS = 60;

videoQuality = 95;

videoProfile = 'MPEG-4'; % or 'Motion JPEG AVI'

if strcmpi(videoProfile, 'Motion JPEG AVI')

combinedFile = 'track_replay_combined_wheels.avi';

else

combinedFile = 'track_replay_combined_wheels.mp4';

end

% ----- vehicle body dimensions [m] -----

% carLength = (veh.lf + veh.lr) * 1.3;

% carWidth = veh.wt * 1.1;

% carHeight = 1.35;

carLength = 4.38;

carWidth = 1.865;

carHeight = 1.285;

% ----- wheel / axle geometry [m] -----

lf = veh.lf; % COM -> front axle

lr = veh.lr; % COM -> rear axle

frontTrack = veh.t1;

rearTrack = veh.t2;

wheelRadius = veh.rw;

wheelWidth = veh.ws;

% Top-view wheel rectangle dimensions

wheelLength2D = wheelRadius * 2;

wheelWidth2D = veh.ws;

% ----- planar 3D visualization settings -----

hcg = veh.hc; % nominal COM height above flat ground

bodyStaticClearance = 0.10; % nominal body clearance at zero roll/pitch

USE_AUTO_NONPENETRATION_LIFT = true;

minBodyClearance = 0.01;

% Optional XY offset of drawn body relative to COM in BODY frame

carBodyOffsetBodyXY = [0 0];

% ----- transparency / appearance -----

carBodyAlpha = 0.2; % transparent so wheels are visible 0.42

wheelColor3D = [0.15 0.15 0.15];

wheelEdgeColor3D = [0.05 0.05 0.05];

wheelMarkerColor = [0.95 0.80 0.10];

carFaceColor = [0.10 0.45 0.85];

carEdgeColor = [0.05 0.05 0.05];

noseColor = [0.95 0.15 0.15];

roadColor = [0.82 0.82 0.82];

wheelColor2D = [0.20 0.20 0.20];

trailColor = [0.85 0.15 0.15];

% ----- wheel spin sign -----

wheelSpinSign = 1; % use -1 if you want the opposite visual spin direction

% ----- top-down options -----

topdownTrailLength_s = 4.0;

showSpeedText = false; % upper-left time/speed/sideslip text disabled

showHeadingArrow2D = true;

% Larger fixed-size marker for the car position in the circuit inset
trackMarkerSize = 6;

% Display-only rotation of the 2D circuit inset
trackDisplayRotation_deg = -27.0; % positive = counterclockwise

% ----- 3D chase camera options -----

% chaseDist = 10.0;

chaseDist = 12.0;

chaseHeight = 4.0;

lookAhead = 15.0;

camSide = 0.0;

cameraViewAngle_deg = 35;

CAMERA_FOLLOWS_ROLL = false;

% ----- interpolation method -----

interpMethod = 'pchip'; % 'linear', 'pchip', 'makima'

% ----- optional road surface strip -----

DRAW_ROAD_SURFACE_IF_SAME_LENGTH = true;

% ----- figure sizes -----

figPos = [120 60 1450 920];

% ----- combined layout -----

mainAxPos = [0.04 0.06 0.93 0.90]; % main 3D axes

insetAxPos = [0.72 0.49 0.275 0.46]; % moved farther right and enlarged vertically; right edge = 0.995, top edge = 0.95

% ----- HUD layout / animation -----

hudAxPos = [0.34 0.72 0.32 0.23]; % upper-middle overlay, aligned with other displays

steeringWheelRatio = 16; % steering wheel angle = ratio * deltaFL

steeringWheelSign = 1; % use -1 if visual direction is opposite

driveBarColor = [0.18 0.72 0.18];

brakeBarColor = [0.85 0.20 0.20];

hudFrameColor = [0.10 0.10 0.10];

hudTextColor = [0.05 0.05 0.05];

hudBgColor = [1 1 1];

% ----- speed mini-plot layout / appearance -----

speedAxPos = [0.04 0.72 0.26 0.23]; % upper-left overlay, aligned with other displays

speedCurveColor = [0.10 0.40 0.90];

speedRefColor = [0.82 0.82 0.82];

speedMarkerColor = [0.95 0.45 0.10];

speedAxBgColor = [1 1 1];

speedTextColor = [0.05 0.05 0.05];

%% ========================= INPUT CHECKS =========================

requiredVars = {'leftX','leftY','rightX','rightY', 't','carX','carY','speed','heading','roll','pitch', 'deltaFL','deltaFR', 'omegaFL','omegaFR','omegaRL','omegaRR', 'Tdrive_norm','Tbrake_norm'};

for i = 1:numel(requiredVars)

if ~exist(requiredVars{i}, 'var')

error('Missing required variable: %s', requiredVars{i});

end

end

if ~exist('leftZ','var') || isempty(leftZ), leftZ = zeros(size(leftX)); end

if ~exist('rightZ','var') || isempty(rightZ), rightZ = zeros(size(rightX)); end

if ~exist('carZ','var') || isempty(carZ), carZ = zeros(size(carX)); end

% Force column vectors

leftX = leftX(:); leftY = leftY(:); leftZ = leftZ(:);

rightX = rightX(:); rightY = rightY(:); rightZ = rightZ(:);

t = t(:);

carX = carX(:);

carY = carY(:);

carZ = carZ(:);

speed = speed(:);

sideslip_angle = sideslip_angle(:);

heading = heading(:);

roll = roll(:);

pitch = pitch(:);

deltaFL = deltaFL(:);

deltaFR = deltaFR(:);

omegaFL = omegaFL(:);

omegaFR = omegaFR(:);

omegaRL = omegaRL(:);

omegaRR = omegaRR(:);

Tdrive_norm = Tdrive_norm(:);

Tbrake_norm = Tbrake_norm(:);

% omega_engine_query_rpm = omega_engine_query_rpm(:);

% gear_query = gear_query(:);

N = numel(t);

if N < 2

error('Time history must contain at least 2 samples.');

end

sameLengthVars = {'carX','carY','carZ','speed','heading','roll','pitch','deltaFL','deltaFR','omegaFL','omegaFR','omegaRL','omegaRR','Tdrive_norm','Tbrake_norm'};

for i = 1:numel(sameLengthVars)

tmp = eval(sameLengthVars{i});

if numel(tmp) ~= N

error('Array size mismatch: %s must have the same length as t.', sameLengthVars{i});

end

end

allFinite = all(isfinite(t)) && all(isfinite(carX)) && all(isfinite(carY)) && all(isfinite(carZ)) && all(isfinite(speed)) && all(isfinite(heading)) && all(isfinite(roll)) && all(isfinite(pitch)) && all(isfinite(deltaFL))&& all(isfinite(deltaFR)) && all(isfinite(omegaFL))&& all(isfinite(omegaFR)) && all(isfinite(omegaRL))&& all(isfinite(omegaRR)) && all(isfinite(Tdrive_norm)) && all(isfinite(Tbrake_norm)) && all(isfinite(leftX)) && all(isfinite(leftY)) && all(isfinite(leftZ)) && all(isfinite(rightX)) && all(isfinite(rightY)) && all(isfinite(rightZ));

if ~allFinite

error('All input arrays must be finite. NaN or Inf found.');

end

if any(diff(t) <= 0)

error('Time vector t must be strictly increasing.');

end

if ANGLES_IN_DEG

heading = deg2rad(heading);

roll = deg2rad(roll);

pitch = deg2rad(pitch);

end

if STEER_ANGLES_IN_DEG

deltaFL = deg2rad(deltaFL);

deltaFR = deg2rad(deltaFR);

end

%% ========================= FLAT GROUND CHECK =========================

allZraw = [leftZ; rightZ; carZ];

groundZ0 = mean(allZraw);

zSpread = max(allZraw) - min(allZraw);

if zSpread > 1e-8

warning('Z data are not exactly constant. A flat ground at groundZ0 = %.6f will be used for visualization.', groundZ0);

end

%% ========================= BUILD VIDEO TIMELINE =========================

dtVideo = 1 / videoFPS;

nFrames = floor((t(end) - t(1)) * videoFPS) + 1;

tv = t(1) + (0:nFrames-1).' * dtVideo;

heading_u = unwrap(heading);

xv = interp1(t, carX, tv, interpMethod);

yv = interp1(t, carY, tv, interpMethod);

vv = interp1(t, speed, tv, interpMethod);

ssav = interp1(t, sideslip_angle, tv, interpMethod);

yawv = interp1(t, heading_u, tv, interpMethod);

rollv = interp1(t, roll, tv, interpMethod);

pitchv = interp1(t, pitch, tv, interpMethod);

deltaFLv = interp1(t, deltaFL, tv, interpMethod);

deltaFRv = interp1(t, deltaFR, tv, interpMethod);

omegaFLv = interp1(t, omegaFL, tv, interpMethod);

omegaFRv = interp1(t, omegaFR, tv, interpMethod);

omegaRLv = interp1(t, omegaRL, tv, interpMethod);

omegaRRv = interp1(t, omegaRR, tv, interpMethod);

Tdrive_normv = interp1(t, Tdrive_norm, tv, interpMethod);

Tbrake_normv = interp1(t, Tbrake_norm, tv, interpMethod);

% Engine_RPM = interp1(t, omega_engine_query_rpm, tv, interpMethod);

% Gear_Selected = interp1(t, gear_query, tv, 'linear');

Tdrive_normv = max(0, min(1, Tdrive_normv));

Tbrake_normv = max(0, min(1, Tbrake_normv));

thetaFLv = cumtrapz(tv, omegaFLv);

thetaFRv = cumtrapz(tv, omegaFRv);

thetaRLv = cumtrapz(tv, omegaRLv);

thetaRRv = cumtrapz(tv, omegaRRv);

Nv = numel(tv);

trailSamplesTop = max(2, round(topdownTrailLength_s * videoFPS));

% Speed mini-plot limits

speedMin = min(vv);

speedMax = max(vv);

if abs(speedMax - speedMin) < 1e-9

speedPad = max(1.0, 0.15*max(abs(speedMax),1));

else

speedPad = 0.08 * (speedMax - speedMin);

end

speedYL = [speedMin - speedPad, speedMax + speedPad];

tSpan = tv(end) - tv(1);

if tSpan <= 0

tSpan = 1;

end

%% ========================= LOCAL GEOMETRY =========================

% Body geometry offset in BODY frame

carBoxCenterOffsetBody = [ ...

carBodyOffsetBodyXY(1), ...

carBodyOffsetBodyXY(2), ...

bodyStaticClearance + carHeight/2 - hcg];

[vertsBodyLocal, facesBodyLocal] = makeBoxAtOffset(carLength, carWidth, carHeight, carBoxCenterOffsetBody);

% 2D body footprint and nose already shifted by the chosen body offset

[carXY_local_raw, noseXY_local_raw] = make2DCarFootprint(carLength, carWidth);

carXY_local = carXY_local_raw + carBodyOffsetBodyXY;

noseXY_local_2D = noseXY_local_raw + carBodyOffsetBodyXY;

% Wheel centers in BODY frame

wheelCenterZBody = wheelRadius - hcg;

cFL = [ lf + carBodyOffsetBodyXY(1), frontTrack/2 + carBodyOffsetBodyXY(2), wheelCenterZBody];

cFR = [ lf + carBodyOffsetBodyXY(1), -frontTrack/2 + carBodyOffsetBodyXY(2), wheelCenterZBody];

cRL = [-lr + carBodyOffsetBodyXY(1), rearTrack/2 + carBodyOffsetBodyXY(2), wheelCenterZBody];

cRR = [-lr + carBodyOffsetBodyXY(1), -rearTrack/2 + carBodyOffsetBodyXY(2), wheelCenterZBody];

wheelCentersBody = [cFL; cFR; cRL; cRR];

% 2D wheel local rectangle in wheel frame

wheelRectLocal2D = make2DWheelRect(wheelLength2D, wheelWidth2D);

% 3D wheel local geometry (axis along local +y)

[XwLoc, YwLoc, ZwLoc] = makeWheelCylinderY(wheelRadius, wheelWidth, 28);

% Small marker line on the wheel so spin is visible

wheelMarkerLine = [0 0 0; wheelRadius 0 0];

% Wheel sample points used for anti-penetration check

wheelPenPtsLocal = makeWheelPenetrationPoints(wheelRadius, wheelWidth, 36);

%% ========================= PRECOMPUTE VISUAL HEIGHT HISTORY =========================

% Body uses yaw + pitch + roll

% Wheels use yaw only + steering (no body pitch/roll)

zVisHist = zeros(Nv,1);

extraLiftHist = zeros(Nv,1);

for k = 1:Nv

psi = yawv(k);

the = pitchv(k);

phi = rollv(k);

pNom = [xv(k); yv(k); groundZ0 + hcg];

Rbody = rotmZYX(psi, the, phi); % body yaw + pitch + roll

Ryaw = rotZ3(psi); % wheels yaw only

minZworld = inf;

% Body: full yaw + pitch + roll

VbodyW = (Rbody * vertsBodyLocal.').';

VbodyW(:,1) = VbodyW(:,1) + pNom(1);

VbodyW(:,2) = VbodyW(:,2) + pNom(2);

VbodyW(:,3) = VbodyW(:,3) + pNom(3);

minZworld = min(minZworld, min(VbodyW(:,3)));

% Wheels: yaw only + steering, no body pitch/roll

deltaVec = [deltaFLv(k); deltaFRv(k); 0; 0];

for iw = 1:4

cw = wheelCentersBody(iw,:).';

Rsteer = rotZ3(deltaVec(iw));

Pbody = (Rsteer * wheelPenPtsLocal.').';

Pbody(:,1) = Pbody(:,1) + cw(1);

Pbody(:,2) = Pbody(:,2) + cw(2);

Pbody(:,3) = Pbody(:,3) + cw(3);

Pworld = (Ryaw * Pbody.').';

Pworld(:,1) = Pworld(:,1) + pNom(1);

Pworld(:,2) = Pworld(:,2) + pNom(2);

Pworld(:,3) = Pworld(:,3) + pNom(3);

minZworld = min(minZworld, min(Pworld(:,3)));

end

extraLift = 0;

if USE_AUTO_NONPENETRATION_LIFT

desiredMinZ = groundZ0 + minBodyClearance;

if minZworld < desiredMinZ

extraLift = desiredMinZ - minZworld;

end

end

extraLiftHist(k) = extraLift;

zVisHist(k) = pNom(3) + extraLift;

end

%% ========================= AXIS LIMITS =========================

allX = [leftX; rightX; xv];

allY = [leftY; rightY; yv];

allZ = [leftZ; rightZ; zVisHist];

xMin = min(allX); xMax = max(allX);

yMin = min(allY); yMax = max(allY);

zMin = min(allZ); zMax = max(allZ);

xySpan = max([xMax - xMin, yMax - yMin, 1]);

marginXY = 0.08 * xySpan + max(carLength, carWidth);

xL = [xMin - marginXY, xMax + marginXY];

yL = [yMin - marginXY, yMax + marginXY];

zL = [groundZ0 - 0.10, max(groundZ0 + 2.5, max(zVisHist) + carHeight + 0.8)];

%% ========================= 2D CIRCUIT DISPLAY ROTATION =========================
% This rotates only the upper circuit inset. It does not change the vehicle data,
% the main 3D animation, or any simulation quantity.
trackDisplayCenter = [ ...
    0.5*(min([leftX; rightX]) + max([leftX; rightX])), ...
    0.5*(min([leftY; rightY]) + max([leftY; rightY]))];

trackDisplayRotation = deg2rad(trackDisplayRotation_deg);

[leftX_disp, leftY_disp] = rotate2DPoints( ...
    leftX, leftY, trackDisplayRotation, trackDisplayCenter);
[rightX_disp, rightY_disp] = rotate2DPoints( ...
    rightX, rightY, trackDisplayRotation, trackDisplayCenter);
[xv_disp, yv_disp] = rotate2DPoints( ...
    xv, yv, trackDisplayRotation, trackDisplayCenter);

allInsetX = [leftX_disp; rightX_disp; xv_disp];
allInsetY = [leftY_disp; rightY_disp; yv_disp];
insetSpan = max([max(allInsetX)-min(allInsetX), ...
                 max(allInsetY)-min(allInsetY), 1]);
insetMargin = 0.010 * insetSpan; % minimal display-only margin so the rotated circuit fills the larger inset
xL_inset = [min(allInsetX)-insetMargin, max(allInsetX)+insetMargin];
yL_inset = [min(allInsetY)-insetMargin, max(allInsetY)+insetMargin];

%% ========================= COMBINED VIDEO =========================

fprintf('Creating combined video: %s\n', combinedFile);

fig = figure('Name','Combined Replay','Color','w','Position',figPos);

set(fig, 'Renderer', 'opengl');

% ------------------------- Main 3D axes -------------------------

ax3 = axes('Parent', fig, 'Position', mainAxPos);

hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');

axis(ax3,'equal');

axis(ax3,'vis3d');

xlim(ax3, xL);

ylim(ax3, yL);

zlim(ax3, zL);

xlabel(ax3,'X East [m]');

ylabel(ax3,'Y North [m]');

zlabel(ax3,'Z Up [m]');

title(ax3,'3D Chase View');

view(ax3, 3);

camproj(ax3, 'perspective');

camva(ax3, cameraViewAngle_deg);

canDrawRoadSurface = DRAW_ROAD_SURFACE_IF_SAME_LENGTH && (numel(leftX) == numel(rightX));

if canDrawRoadSurface

Xroad = [leftX, rightX];

Yroad = [leftY, rightY];

Zroad = groundZ0 * ones(size(Xroad));

surf(ax3, Xroad, Yroad, Zroad, 'FaceColor', roadColor, 'EdgeColor', 'none', 'FaceAlpha', 1.0);

end

plot3(ax3, leftX, leftY, groundZ0*ones(size(leftX)), 'k-', 'LineWidth', 1.2);

plot3(ax3, rightX, rightY, groundZ0*ones(size(rightX)), 'k-', 'LineWidth', 1.2);

% Full racing line intentionally hidden so the future path is not shown.
% plot3(ax3, xv, yv, zVisHist, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);

hTrail3 = plot3(ax3, NaN, NaN, NaN, '-', 'Color', trailColor, 'LineWidth', 1.8);

% Car transform hierarchy:

% hCarYawTF : translation + yaw only

% hBodyTF : body-only pitch/roll relative to yaw frame

hCarYawTF = hgtransform('Parent', ax3);

hBodyTF = hgtransform('Parent', hCarYawTF);

% Body patch

patch('Parent', hBodyTF, 'Vertices', vertsBodyLocal, 'Faces', facesBodyLocal, 'FaceColor', carFaceColor, 'EdgeColor', carEdgeColor, 'LineWidth', 1.0, 'FaceAlpha', carBodyAlpha);

% Nose line

noseLocal3D = [ ...

carBoxCenterOffsetBody(1) + 0.10*carLength, carBoxCenterOffsetBody(2), carBoxCenterOffsetBody(3);

carBoxCenterOffsetBody(1) + 0.55*carLength, carBoxCenterOffsetBody(2), carBoxCenterOffsetBody(3)];

line('Parent', hBodyTF, 'XData', noseLocal3D(:,1), 'YData', noseLocal3D(:,2), 'ZData', noseLocal3D(:,3), 'Color', noseColor, 'LineWidth', 2.5);

% Wheel transform hierarchy

hWheelMount = gobjects(4,1);

hWheelSteer = gobjects(4,1);

hWheelSpin = gobjects(4,1);

for iw = 1:4

% Wheels are attached to yaw-only frame, not body pitch/roll frame

hWheelMount(iw) = hgtransform('Parent', hCarYawTF);

hWheelSteer(iw) = hgtransform('Parent', hWheelMount(iw));

hWheelSpin(iw) = hgtransform('Parent', hWheelSteer(iw));

surf('Parent', hWheelSpin(iw), 'XData', XwLoc, 'YData', YwLoc, 'ZData', ZwLoc, 'FaceColor', wheelColor3D, 'EdgeColor', wheelEdgeColor3D, 'FaceAlpha', 1.0);

line('Parent', hWheelSpin(iw), 'XData', wheelMarkerLine(:,1), 'YData', wheelMarkerLine(:,2), 'ZData', wheelMarkerLine(:,3), 'Color', wheelMarkerColor, 'LineWidth', 2.2);

end

for iw = 1:4

set(hWheelMount(iw), 'Matrix', makeTransform(eye(3), wheelCentersBody(iw,:).'));

end

light(ax3, 'Position', [1 1 2], 'Style', 'infinite');

lighting(ax3, 'gouraud');

material(ax3, 'dull');

% ------------------------- Inset 2D axes -------------------------

ax2 = axes('Parent', fig, 'Position', insetAxPos);

% Preserve the requested inset box size when axis equal is applied.
set(ax2, 'PositionConstraint', 'innerposition');

hold(ax2,'on'); grid(ax2,'off'); box(ax2,'off'); % transparent frameless track display

axis(ax2,'equal');

axis(ax2,'manual');

xlim(ax2, xL_inset);

ylim(ax2, yL_inset);

set(ax2, 'Color', 'none', 'LineWidth', 1.0, 'Layer', 'top', ...
    'XColor', 'none', 'YColor', 'none'); % transparent track display with no rectangular frame

set(ax2, 'XTick', [], 'YTick', []);

title(ax2, 'Top view', 'FontSize', 10);

plot(ax2, leftX_disp, leftY_disp, 'k-', 'LineWidth', 0.9);

plot(ax2, rightX_disp, rightY_disp, 'k-', 'LineWidth', 0.9);

% Full racing line intentionally hidden in the circuit inset.
% plot(ax2, xv_disp, yv_disp, '-', 'Color', [0.80 0.80 0.80], 'LineWidth', 0.7);

hTrail2 = plot(ax2, NaN, NaN, '-', 'Color', trailColor, 'LineWidth', 1.2);

hCar2D = patch(ax2, 'XData', NaN, 'YData', NaN, 'FaceColor', carFaceColor, 'EdgeColor', 'k', 'LineWidth', 0.9, 'FaceAlpha', carBodyAlpha);

if showHeadingArrow2D

hNose2D = plot(ax2, NaN, NaN, '-', 'Color', noseColor, 'LineWidth', 1.8);

end

hWheel2D = gobjects(4,1);

for iw = 1:4

hWheel2D(iw) = patch(ax2, 'XData', NaN, 'YData', NaN, 'FaceColor', wheelColor2D, 'EdgeColor', 'k', 'LineWidth', 0.7);

end

% Fixed-size current-position dot, kept visible at full-circuit scale
hTrackMarker = plot(ax2, NaN, NaN, 'o', ...
    'MarkerSize', trackMarkerSize, ...
    'MarkerFaceColor', trailColor, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.0);

% ------------------------- HUD axes -------------------------

axHUD = axes('Parent', fig, 'Position', hudAxPos);

hold(axHUD,'on');

axis(axHUD, [0 1 0 1]);

axis(axHUD, 'off');

set(axHUD, 'Color', 'none');

% HUD background panel

patch(axHUD, [0 1 1 0], [0 0 1 1], hudBgColor, ...
    'FaceAlpha', 0.70, 'EdgeColor', hudFrameColor, 'LineWidth', 1.2); % full axes for exact display alignment

% Steering wheel geometry

swCenter = [0.50, 0.56];

swRimR = 0.17;

swHubR = 0.040;

swSpokeR1 = 0.055;

swSpokeR2 = 0.145;

swBaseAngles = deg2rad([90, 210, 330]);

tt = linspace(0, 2*pi, 220);

% Rim and hub

plot(axHUD, swCenter(1) + swRimR*cos(tt), swCenter(2) + swRimR*sin(tt), '-', 'Color', hudFrameColor, 'LineWidth', 4.0);

patch(axHUD, swCenter(1) + swHubR*cos(tt), swCenter(2) + swHubR*sin(tt), hudFrameColor, 'EdgeColor', hudFrameColor);

% Three spokes

hSWSpokes = gobjects(3,1);

for ispk = 1:3

hSWSpokes(ispk) = line(axHUD, NaN, NaN, 'Color', hudFrameColor, 'LineWidth', 3.0);

end

% Rotating rim marker so the steering wheel animation is more obvious

swMarkerR1 = swRimR - 0.015;

swMarkerR2 = swRimR + 0.010;

hSWMarker = line(axHUD, NaN, NaN, 'Color', wheelMarkerColor, 'LineWidth', 3.0);

% Drive / brake bar geometry

barY0 = 0.18;

barH = 0.62;

barW = 0.065;

brakeBarXC = 0.18; % brake pedal on the left

driveBarXC = 0.82; % drive pedal on the right

driveBarX = driveBarXC + 0.5*barW*[-1 1 1 -1];

brakeBarX = brakeBarXC + 0.5*barW*[-1 1 1 -1];

barYFrame = [barY0 barY0 barY0+barH barY0+barH];

% Bar backgrounds

patch(axHUD, driveBarX, barYFrame, [1 1 1], 'FaceAlpha', 0.75, 'EdgeColor', hudFrameColor, 'LineWidth', 1.4);

patch(axHUD, brakeBarX, barYFrame, [1 1 1], 'FaceAlpha', 0.75, 'EdgeColor', hudFrameColor, 'LineWidth', 1.4);

% Dynamic bar fills

hDriveBarFill = patch(axHUD, driveBarX, [barY0 barY0 barY0 barY0], driveBarColor, 'EdgeColor', 'none', 'FaceAlpha', 0.95);

hBrakeBarFill = patch(axHUD, brakeBarX, [barY0 barY0 barY0 barY0], brakeBarColor, 'EdgeColor', 'none', 'FaceAlpha', 0.95);

% HUD labels

text(axHUD, driveBarXC, 0.10, 'TDrive', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'Color', hudTextColor, 'FontSize', 10);

text(axHUD, brakeBarXC, 0.10, 'TBrake', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'Color', hudTextColor, 'FontSize', 10);

title(axHUD, 'Steering wheel (ratio 16)', ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'Color', hudTextColor);

hSWText = text(axHUD, swCenter(1), 0.16, '', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'Color', hudTextColor, 'FontSize', 10);

% ------------------------- Speed mini-plot axes -------------------------

axSpd = axes('Parent', fig, 'Position', speedAxPos);

hold(axSpd,'on'); grid(axSpd,'off'); box(axSpd,'off'); % transparent frameless speed display

xlim(axSpd, [tv(1) tv(end)]);

ylim(axSpd, speedYL);

set(axSpd, 'Color', 'none', 'LineWidth', 1.0, 'Layer', 'top', ...
    'XColor', 'none', 'YColor', 'none'); % transparent speed-display background with no rectangular frame

title(axSpd, 'Speed', 'FontSize', 10);

set(axSpd, 'XTick', [], 'YTick', []);

% Full grey speed reference intentionally hidden so only completed motion is shown.

% plot(axSpd, tv, vv, '-', 'Color', speedRefColor, 'LineWidth', 1.0);

% Dynamic speed curve, cursor, marker, value text

hSpeedPast = plot(axSpd, NaN, NaN, '-', 'Color', speedCurveColor, 'LineWidth', 1.8);

hSpeedCursor = line(axSpd, [NaN NaN], speedYL, 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.0);

hSpeedMarker = plot(axSpd, NaN, NaN, 'o', 'MarkerSize', 6, 'MarkerFaceColor', speedMarkerColor, 'MarkerEdgeColor', [0.20 0.20 0.20], 'LineWidth', 1.0);

hSpeedVal = text(axSpd, NaN, NaN, '', 'FontWeight', 'bold', 'Color', speedTextColor, 'BackgroundColor', 'none', 'Margin', 2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Clipping', 'on');

% ------------------------- Figure annotation -------------------------

% Upper-left time/speed/sideslip annotation intentionally disabled.

% ------------------------- Video writer -------------------------

vw = VideoWriter(combinedFile, videoProfile);

vw.FrameRate = videoFPS;

vw.Quality = videoQuality;

open(vw);

worldUp = [0; 0; 1];

speedRange = speedYL(2) - speedYL(1);

for k = 1:Nv

% ===================== 3D update =====================

psi = yawv(k);

the = pitchv(k);

phi = rollv(k);

Rwb = rotmZYX(psi, the, phi); % full body orientation for camera only

Ryaw = rotZ3(psi); % yaw-only frame for wheels / chassis base

Rpr = rotmZYX(0, the, phi); % pitch/roll only, applied to body relative to yaw frame

pVis = [xv(k); yv(k); zVisHist(k)];

set(hCarYawTF, 'Matrix', makeTransform(Ryaw, pVis));

set(hBodyTF, 'Matrix', makeTransform(Rpr, [0;0;0]));

% Wheel steer + spin

deltaVec = [deltaFLv(k); deltaFRv(k); 0; 0];

thetaVec = [thetaFLv(k); thetaFRv(k); thetaRLv(k); thetaRRv(k)];

for iw = 1:4

set(hWheelSteer(iw), 'Matrix', makeTransform(rotZ3(deltaVec(iw)), [0;0;0]));

set(hWheelSpin(iw), 'Matrix', makeTransform(rotY3(wheelSpinSign * thetaVec(iw)), [0;0;0]));

end

% 3D trail

k0 = max(1, k - trailSamplesTop + 1);

% set(hTrail3, 'XData', xv(k0:k), 'YData', yv(k0:k), 'ZData', zVisHist(k0:k));

% Camera follows body orientation

fwd = Rwb(:,1);

left = Rwb(:,2);

up = Rwb(:,3);

camPos = pVis - chaseDist * fwd + chaseHeight * worldUp + camSide * left;

camTgt = pVis + lookAhead * fwd;

campos(ax3, camPos.');

camtarget(ax3, camTgt.');

if CAMERA_FOLLOWS_ROLL

camup(ax3, up.');

else

camup(ax3, worldUp.');

end

% ===================== 2D inset update =====================

xc = xv(k);

yc = yv(k);

R2 = rot2(psi);

% Body

carXY_world = (R2 * carXY_local.').';

carXY_world(:,1) = carXY_world(:,1) + xc;

carXY_world(:,2) = carXY_world(:,2) + yc;

[carXY_disp_x, carXY_disp_y] = rotate2DPoints( ...
    carXY_world(:,1), carXY_world(:,2), ...
    trackDisplayRotation, trackDisplayCenter);
set(hCar2D, 'XData', carXY_disp_x, 'YData', carXY_disp_y);

if showHeadingArrow2D

noseXY_world = (R2 * noseXY_local_2D.').';

noseXY_world(:,1) = noseXY_world(:,1) + xc;

noseXY_world(:,2) = noseXY_world(:,2) + yc;

[noseXY_disp_x, noseXY_disp_y] = rotate2DPoints( ...
    noseXY_world(:,1), noseXY_world(:,2), ...
    trackDisplayRotation, trackDisplayCenter);
set(hNose2D, 'XData', noseXY_disp_x, 'YData', noseXY_disp_y);

end

% Wheels

for iw = 1:4

cw = wheelCentersBody(iw,1:2).';

Rsteer2 = rot2(deltaVec(iw));

ptsBody = (Rsteer2 * wheelRectLocal2D.').';

ptsBody(:,1) = ptsBody(:,1) + cw(1);

ptsBody(:,2) = ptsBody(:,2) + cw(2);

ptsWorld = (R2 * ptsBody.').';

ptsWorld(:,1) = ptsWorld(:,1) + xc;

ptsWorld(:,2) = ptsWorld(:,2) + yc;

[wheelXY_disp_x, wheelXY_disp_y] = rotate2DPoints( ...
    ptsWorld(:,1), ptsWorld(:,2), ...
    trackDisplayRotation, trackDisplayCenter);
set(hWheel2D(iw), 'XData', wheelXY_disp_x, 'YData', wheelXY_disp_y);

end

% Larger current-position dot in the circuit inset
[xc_disp, yc_disp] = rotate2DPoints( ...
    xc, yc, trackDisplayRotation, trackDisplayCenter);
set(hTrackMarker, 'XData', xc_disp, 'YData', yc_disp);

% 2D trail

% Racing-line/trail display intentionally hidden.
% [trailX_disp, trailY_disp] = rotate2DPoints( ...
%     xv(k0:k), yv(k0:k), trackDisplayRotation, trackDisplayCenter);
% set(hTrail2, 'XData', trailX_disp, 'YData', trailY_disp);

% ===================== HUD update =====================

swAngle = steeringWheelSign * steeringWheelRatio * deltaFLv(k);

% Steering wheel spokes

for ispk = 1:3

ang = swBaseAngles(ispk) + swAngle;

x1 = swCenter(1) + swSpokeR1*cos(ang);

y1 = swCenter(2) + swSpokeR1*sin(ang);

x2 = swCenter(1) + swSpokeR2*cos(ang);

y2 = swCenter(2) + swSpokeR2*sin(ang);

set(hSWSpokes(ispk), 'XData', [x1 x2], 'YData', [y1 y2]);

end

% Steering wheel rim marker

angm = pi/2 + swAngle;

xm1 = swCenter(1) + swMarkerR1*cos(angm);

ym1 = swCenter(2) + swMarkerR1*sin(angm);

xm2 = swCenter(1) + swMarkerR2*cos(angm);

ym2 = swCenter(2) + swMarkerR2*sin(angm);

set(hSWMarker, 'XData', [xm1 xm2], 'YData', [ym1 ym2]);

% Drive / brake bars

driveVal = Tdrive_normv(k);

brakeVal = Tbrake_normv(k);

driveYTop = barY0 + barH * driveVal;

brakeYTop = barY0 + barH * brakeVal;

set(hDriveBarFill, 'YData', [barY0 barY0 driveYTop driveYTop]);

set(hBrakeBarFill, 'YData', [barY0 barY0 brakeYTop brakeYTop]);

set(hSWText, 'String', sprintf('\\delta_{Front} = %.1f^\\circ', rad2deg(deltaFLv(k))));

% ===================== Speed mini-plot update =====================

set(hSpeedPast, 'XData', tv(1:k), 'YData', vv(1:k));

set(hSpeedCursor, 'XData', [tv(k) tv(k)], 'YData', speedYL);

set(hSpeedMarker, 'XData', tv(k), 'YData', vv(k));

% Put current speed value beside the marker, while keeping it inside axes

xText = tv(k) + 0.025 * tSpan;

yText = vv(k) + 0.04 * speedRange;

if xText > tv(end) - 0.12*tSpan

xText = tv(k) - 0.17*tSpan;

end

if yText > speedYL(2) - 0.06*speedRange

yText = vv(k) - 0.08*speedRange;

end

if yText < speedYL(1) + 0.04*speedRange

yText = vv(k) + 0.04*speedRange;

end

set(hSpeedVal, 'Position', [xText yText 0], 'String', sprintf('%.2f m/s', vv(k)));

% ===================== Text =====================

% % % if showSpeedText

% % % ann.String = sprintf(['t = %.2f s\nspeed = %.2f m/s\n' ...

% % % 'visual COM height = %.3f m\nextra lift = %.3f m'], ...

% % % tv(k), vv(k), zVisHist(k) - groundZ0, extraLiftHist(k));

% % % end

% % if showSpeedText

% % ann.String = sprintf(['t = %.2f s\nspeed = %.2f m/s\nsideslip = %.2f deg\n'], ...

% % tv(k), vv(k), rad2deg(ssav(k)));

% % end

% if showSpeedText

% ann.String = sprintf(['t = %.2f s\nspeed = %.2f m/s\nsideslip = %.2f deg\nEngineRPM = %.2f RPM\n'], ...

% tv(k), vv(k), rad2deg(ssav(k)), Engine_RPM(k));

% end

% if showSpeedText

% ann.String = sprintf(['t = %.2f s\nspeed = %.2f m/s\nsideslip = %.2f deg\nEngineRPM = %.2f RPM\nGearPosition = %.0f \n'], ...

% tv(k), vv(k), rad2deg(ssav(k)), Engine_RPM(k), Gear_Selected(k));

% end

% Upper-left time/speed/sideslip annotation intentionally disabled.

drawnow limitrate;

frame = getframe(fig);

writeVideo(vw, frame);

end

close(vw);

fprintf('Combined video saved.\n');

disp('Done.');

disp(['Saved: ', combinedFile]);

%% ========================= LOCAL FUNCTIONS =========================

function [carXY_local, noseXY_local] = make2DCarFootprint(L, W)

carXY_local = [ ...

L/2, W/2;

L/2, -W/2;

-L/2, -W/2;

-L/2, W/2 ];

noseXY_local = [ ...

0.10*L, 0;

0.55*L, 0 ];

end

function wheelRect = make2DWheelRect(L, W)

wheelRect = [ ...

L/2, W/2;

L/2, -W/2;

-L/2, -W/2;

-L/2, W/2 ];

end

function [V, F] = makeBoxAtOffset(L, W, H, centerBody)

cx = centerBody(1);

cy = centerBody(2);

cz = centerBody(3);

x1 = cx - L/2; x2 = cx + L/2;

y1 = cy - W/2; y2 = cy + W/2;

z1 = cz - H/2; z2 = cz + H/2;

V = [ ...

x1 y1 z1;

x2 y1 z1;

x2 y2 z1;

x1 y2 z1;

x1 y1 z2;

x2 y1 z2;

x2 y2 z2;

x1 y2 z2 ];

F = [ ...

1 2 3 4;

5 6 7 8;

1 2 6 5;

2 3 7 6;

3 4 8 7;

4 1 5 8 ];

end

function [X, Y, Z] = makeWheelCylinderY(radius, width, nCirc)

% Cylinder with axis along local +y, centered at y=0

th = linspace(0, 2*pi, nCirc);

yEnds = [-width/2, width/2];

[TH, YY] = meshgrid(th, yEnds);

X = radius * cos(TH);

Y = YY;

Z = radius * sin(TH);

end

function P = makeWheelPenetrationPoints(radius, width, nCirc)

% Sample points on wheel cylinder for min-z checking

th = linspace(0, 2*pi, nCirc).';

x = radius * cos(th);

z = radius * sin(th);

P1 = [x, -width/2*ones(size(x)), z];

P2 = [x, width/2*ones(size(x)), z];

P = [P1; P2];

end

function [xr, yr] = rotate2DPoints(x, y, angle, centerXY)

originalSizeX = size(x);
originalSizeY = size(y);

xy = [x(:)-centerXY(1), y(:)-centerXY(2)];
xyRot = (rot2(angle) * xy.').';

xr = reshape(xyRot(:,1)+centerXY(1), originalSizeX);
yr = reshape(xyRot(:,2)+centerXY(2), originalSizeY);

end

function R = rotmZYX(yaw, pitch, roll)

cy = cos(yaw); sy = sin(yaw);

cp = cos(pitch); sp = sin(pitch);

cr = cos(roll); sr = sin(roll);

Rz = [ cy -sy 0;

sy cy 0;

0 0 1 ];

Ry = [ cp 0 sp;

0 1 0;

-sp 0 cp ];

Rx = [ 1 0 0;

0 cr -sr;

0 sr cr ];

R = Rz * Ry * Rx;

end

function R = rot2(a)

R = [cos(a), -sin(a);

sin(a), cos(a)];

end

function R = rotY3(a)

ca = cos(a); sa = sin(a);

R = [ ca 0 sa;

0 1 0;

-sa 0 ca ];

end

function R = rotZ3(a)

ca = cos(a); sa = sin(a);

R = [ ca -sa 0;

sa ca 0;

0 0 1 ];

end

function T = makeTransform(R, p)

T = eye(4);

T(1:3,1:3) = R;

T(1:3,4) = p(:);

end