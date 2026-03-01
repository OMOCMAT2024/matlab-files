% clear
% close all
% clc

x = -1:0.001:1;

rou_smooth = 1e-1; offset = 0.00318298 + 0.0285425 + 3.74306e-8;
% rou_smooth = 1e-2; offset = 0.00318298;
% rou_smooth = 1e-3; offset = 0.00318298-0.00286467;
% rou_smooth = 1e-4; offset = 0.00318298-0.00286467-0.000286479;

my_smooth_ramp_func1 = (0 - x) .* ( atan2(0 - x, rou_smooth)/pi + 1/2 ) + offset;

figure()
plot(x, my_smooth_ramp_func1)
xlabel('x')
ylabel('my smooth ramp func1')
grid minor
axis equal

my_smooth_ramp_func2 = (0 + x) .* ( atan2(0 + x, rou_smooth)/pi + 1/2 ) + offset;

figure()
plot(x, my_smooth_ramp_func2)
xlabel('x')
ylabel('my smooth ramp func2')
grid minor
axis equal