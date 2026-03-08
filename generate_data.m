%% generate_demo_planar_track_data_with_wheels.m
% Artificial planar-track + vehicle + wheel data generator
%
% Variables created in workspace:
%   t
%   carX, carY, carZ
%   speed
%   heading, roll, pitch
%   leftX,leftY,leftZ
%   rightX,rightY,rightZ
%   deltaFL, deltaFR              % front steering angles [rad]
%   omegaFL, omegaFR, omegaRL, omegaRR   % wheel angular rates [rad/s]
%
% Global frame:
%   X = East
%   Y = North
%   Z = Up
%
% Heading convention:
%   positive CCW from +X (East), in [-pi, pi]
%
% Artificial steering convention used here:
%   positive steering = left turn = CCW about body +z

clear; clc; close all;

%% ========================= USER SETTINGS =========================

% Track geometry
Ntrack = 1400;
trackWidth = 10.0;   % [m]

% Vehicle history
Nstate = 2200;
nLaps  = 1.25;
meanDt = 0.04;       % mean time step [s], but nonuniform

% Lateral offset of vehicle relative to centerline
maxLateralOffset = 1.0;   % [m]

% Synthetic attitude
rollAmp_deg  = 4.0;
pitchAmp_deg = 2.5;

% Vehicle / wheel geometry used for artificial wheel data
lf = 1.40;           % COM -> front axle [m]
lr = 1.40;           % COM -> rear axle  [m]
wheelbase  = lf + lr;
frontTrack = 1.60;
rearTrack  = 1.58;
wheelRadius = 0.34;  % [m]

% Save optional MAT file
saveMatFile = true;
matFileName = 'demo_planar_track_vehicle_wheels.mat';

%% ========================= TRACK CENTERLINE =========================

thetaTrack = linspace(0, 2*pi, Ntrack+1).';
thetaTrack(end) = [];

[xc, yc, dxc, dyc] = centerlineShape(thetaTrack);

ds_dtheta = hypot(dxc, dyc);
tx = dxc ./ ds_dtheta;
ty = dyc ./ ds_dtheta;

nx = -ty;
ny =  tx;

leftX  = xc + 0.5*trackWidth*nx;
leftY  = yc + 0.5*trackWidth*ny;
leftZ  = zeros(size(leftX));

rightX = xc - 0.5*trackWidth*nx;
rightY = yc - 0.5*trackWidth*ny;
rightZ = zeros(size(rightX));

%% ========================= NONUNIFORM TIME VECTOR =========================

k = (0:Nstate-2).';
dt = meanDt * (1 ...
    + 0.18*sin(2*pi*k/max(1,numel(k)-1)) ...
    + 0.08*sin(2*pi*3*k/max(1,numel(k)-1) + 0.7) ...
    + 0.04*randn(size(k)));

dt = max(dt, 0.35*meanDt);
t = [0; cumsum(dt)];
Tend = t(end);

%% ========================= VEHICLE TRAJECTORY =========================

u = t / Tend;

thetaVeh_unwrapped = 2*pi*nLaps*u ...
                   + 0.08*sin(2*pi*2*u) ...
                   + 0.03*sin(2*pi*5*u + 0.4);

thetaVeh = mod(thetaVeh_unwrapped, 2*pi);

[xcVeh, ycVeh, dxcVeh, dycVeh] = centerlineShape(thetaVeh);

dsVeh = hypot(dxcVeh, dycVeh);
txVeh = dxcVeh ./ dsVeh;
tyVeh = dycVeh ./ dsVeh;

nxVeh = -tyVeh;
nyVeh =  txVeh;

latOffset = maxLateralOffset * sin(2*pi*1.4*u + 0.2) .* (0.7 + 0.3*cos(2*pi*0.7*u));

carX = xcVeh + latOffset .* nxVeh;
carY = ycVeh + latOffset .* nyVeh;
carZ = zeros(size(carX));   % planar

%% ========================= SPEED AND HEADING =========================

vx = gradient(carX, t);
vy = gradient(carY, t);

speed = hypot(vx, vy);

heading = atan2(vy, vx);
heading = atan2(sin(heading), cos(heading));   % wrap to [-pi, pi]

%% ========================= ARTIFICIAL ROLL / PITCH =========================

rollAmp  = deg2rad(rollAmp_deg);
pitchAmp = deg2rad(pitchAmp_deg);

roll  = rollAmp  * (0.60*sin(2*pi*2.2*u + 0.4) + 0.40*sin(2*pi*5.3*u));
pitch = pitchAmp * (0.75*sin(2*pi*1.7*u - 0.5) + 0.25*sin(2*pi*4.1*u + 1.0));

%% ========================= ARTIFICIAL FRONT STEERING =========================
% Build a plausible steering history from yaw rate and speed, then use a
% simple Ackermann split to get deltaFL / deltaFR.

yaw_u   = unwrap(heading);
yawRate = gradient(yaw_u, t);

speedSafe = max(speed, 1.0);
deltaEq = atan(wheelbase * yawRate ./ speedSafe);

% Clip to a reasonable range
deltaMax = deg2rad(26);
deltaEq = max(min(deltaEq, deltaMax), -deltaMax);

deltaFL = zeros(size(deltaEq));
deltaFR = zeros(size(deltaEq));

for i = 1:numel(deltaEq)
    deq = deltaEq(i);

    if abs(deq) < 1e-6
        deltaFL(i) = 0;
        deltaFR(i) = 0;
    else
        R = wheelbase / tan(abs(deq));  % positive turning radius magnitude

        din  = atan(wheelbase / max(R - frontTrack/2, 0.20));
        dout = atan(wheelbase / (R + frontTrack/2));

        if deq > 0
            % left turn positive
            deltaFL(i) =  din;   % inner
            deltaFR(i) =  dout;  % outer
        else
            deltaFL(i) = -dout;  % outer
            deltaFR(i) = -din;   % inner
        end
    end
end

%% ========================= ARTIFICIAL WHEEL ANGULAR RATES =========================
% Use a kinematic approximation at each wheel center.

% Wheel center locations in BODY frame: +x forward, +y left
pFL = [ lf;  frontTrack/2];
pFR = [ lf; -frontTrack/2];
pRL = [-lr;  rearTrack/2];
pRR = [-lr; -rearTrack/2];

uLong = speed;   % approximate longitudinal speed in body frame
rBody = yawRate;

omegaFL = zeros(size(t));
omegaFR = zeros(size(t));
omegaRL = zeros(size(t));
omegaRR = zeros(size(t));

for i = 1:numel(t)
    u_i = uLong(i);
    r_i = rBody(i);

    % Velocity of each wheel center in BODY frame
    vFL = [u_i - r_i*pFL(2); r_i*pFL(1)];
    vFR = [u_i - r_i*pFR(2); r_i*pFR(1)];
    vRL = [u_i - r_i*pRL(2); r_i*pRL(1)];
    vRR = [u_i - r_i*pRR(2); r_i*pRR(1)];

    eFL = [cos(deltaFL(i)); sin(deltaFL(i))];
    eFR = [cos(deltaFR(i)); sin(deltaFR(i))];
    eR  = [1; 0];

    vrollFL = dot(vFL, eFL);
    vrollFR = dot(vFR, eFR);
    vrollRL = dot(vRL, eR);
    vrollRR = dot(vRR, eR);

    omegaFL(i) = max(vrollFL / wheelRadius, 0.1);
    omegaFR(i) = max(vrollFR / wheelRadius, 0.1);
    omegaRL(i) = max(vrollRL / wheelRadius, 0.1);
    omegaRR(i) = max(vrollRR / wheelRadius, 0.1);
end

% Small smooth modulation so wheel spin is not too uniform
omegaFL = omegaFL .* (1 + 0.015*sin(2*pi*2.7*u + 0.1));
omegaFR = omegaFR .* (1 + 0.015*sin(2*pi*2.5*u + 0.4));
omegaRL = omegaRL .* (1 + 0.012*sin(2*pi*2.0*u + 0.2));
omegaRR = omegaRR .* (1 + 0.012*sin(2*pi*2.2*u + 0.6));

%% ========================= SANITY PLOTS =========================

figure('Color','w','Position',[100 100 1350 900]);

subplot(2,3,1); hold on; grid on; box on; axis equal;
plot(leftX, leftY, 'k-', 'LineWidth', 1.2);
plot(rightX, rightY, 'k-', 'LineWidth', 1.2);
plot(carX, carY, 'r-', 'LineWidth', 1.0);
xlabel('X East [m]');
ylabel('Y North [m]');
title('Track and COM trajectory');

subplot(2,3,2); hold on; grid on; box on;
plot(t, speed, 'b-');
xlabel('t [s]');
ylabel('speed [m/s]');
title('Speed');

subplot(2,3,3); hold on; grid on; box on;
plot(t, heading, 'k-');
xlabel('t [s]');
ylabel('heading [rad]');
title('Heading');

subplot(2,3,4); hold on; grid on; box on;
plot(t, roll, 'r-', 'DisplayName','roll');
plot(t, pitch,'b-', 'DisplayName','pitch');
xlabel('t [s]');
ylabel('angle [rad]');
title('Roll / pitch');
legend('Location','best');

subplot(2,3,5); hold on; grid on; box on;
plot(t, deltaFL, 'm-', 'DisplayName','deltaFL');
plot(t, deltaFR, 'c-', 'DisplayName','deltaFR');
xlabel('t [s]');
ylabel('steer [rad]');
title('Front steering');
legend('Location','best');

subplot(2,3,6); hold on; grid on; box on;
plot(t, omegaFL, 'DisplayName','omegaFL');
plot(t, omegaFR, 'DisplayName','omegaFR');
plot(t, omegaRL, 'DisplayName','omegaRL');
plot(t, omegaRR, 'DisplayName','omegaRR');
xlabel('t [s]');
ylabel('\omega [rad/s]');
title('Wheel angular rates');
legend('Location','best');

%% ========================= OPTIONAL SAVE =========================

if saveMatFile
    save(matFileName, ...
        't', ...
        'carX','carY','carZ', ...
        'speed', ...
        'heading','roll','pitch', ...
        'leftX','leftY','leftZ', ...
        'rightX','rightY','rightZ', ...
        'deltaFL','deltaFR', ...
        'omegaFL','omegaFR','omegaRL','omegaRR');
    fprintf('Saved demo data to: %s\n', matFileName);
end

fprintf('Demo data with wheels generated successfully.\n');

%% ========================= LOCAL FUNCTION =========================

function [x, y, dx, dy] = centerlineShape(theta)
% Closed wavy oval centerline and derivatives wrt theta.
    x = 110*cos(theta) + 18*cos(3*theta + 0.25);
    y =  72*sin(theta) + 10*sin(2*theta - 0.35);

    dx = -110*sin(theta) - 54*sin(3*theta + 0.25);
    dy =   72*cos(theta) + 20*cos(2*theta - 0.35);
end