close all

veh = my_params();

track_opt = sol.track_opt;


%% States
s_end = track.s(end);

figure()
plot(sol.t_opt,sol.x_opt(1,:)*3.6), grid on, hold on
ylabel('ux [km/h]')
xlabel('Time, [s]')
% xlim([0 s_end])
title('ux')

figure()
plot(sol.t_opt,sol.x_opt(2,:)*3.6), grid on, hold on
ylabel('uy [km/h]')
xlabel('Time, [s]')
% xlim([0 s_end])
title('uy')

figure();plot(sol.t_opt, rad2deg(sol.x_opt(3,:))), grid on
ylabel('r [deg/s]')
xlabel('Time, [s]')
% xlim([0 s_end])
title('Yaw rate')

figure();plot(sol.t_opt, sol.x_opt(8,:)), grid on
ylabel('n [deg]')
xlabel('Time, [s]')
% xlim([0 s_end])
title('Lateral displacement')

figure();plot(sol.t_opt, rad2deg(sol.x_opt(9,:))), grid on
ylabel('\xi [deg]')
xlabel('Time, [s]')
% xlim([0 s_end])
title('Angle between track reference and vehicle front or total speed vector?')


%%
figure()
plot(sol.t_opt,sol.acc(:,1))
hold on
plot(sol.t_opt,sol.acc(:,2))
grid on
ylabel('[m/s2]')
xlabel('Time, [s]')
legend('ax', 'ay')
% xlim([0 s_end])
xlabel('Time, [s]')
title('Longitudinal and lateral acceleration')

figure()
scatter(sol.acc(:,2)/veh.g, sol.acc(:,1)/veh.g, 'o')
grid minor
axis equal
xlabel('ay, [g]')
ylabel('ax, [g]')
ylim([-1.5 0.7])

%% overlay GG with vairables of interest
figure()
scatter(sol.acc(:,2)/veh.g, sol.acc(:,1)/veh.g, 5, sol.x_opt(1,:)*3.6), hold on, grid minor
cb=colorbar;
ylabel(cb,'V [km/h]','FontSize',12)
axis equal
xlabel('ay, [g]')
ylabel('ax, [g]')
ylim([-1.5 0.7])
title('GG with speed overlay')

figure()
scatter(sol.acc(:,2)/veh.g, sol.acc(:,1)/veh.g, 5, sol.mu(:,1)), hold on, grid minor
cb=colorbar;
ylabel(cb,'FL tyre usage [-]','FontSize',12)
axis equal
xlabel('ay, [g]')
ylabel('ax, [g]')
ylim([-1.5 0.7])
title('GG with FL tyre usage')

figure()
scatter(sol.acc(:,2)/veh.g, sol.acc(:,1)/veh.g, 5, sol.mu(:,2)), hold on, grid minor
cb=colorbar;
ylabel(cb,'FR tyre usage [-]','FontSize',12)
axis equal
xlabel('ay, [g]')
ylabel('ax, [g]')
ylim([-1.5 0.7])
title('GG with FR tyre usage')

figure()
scatter(sol.acc(:,2)/veh.g, sol.acc(:,1)/veh.g, 5, sol.mu(:,3)), hold on, grid minor
cb=colorbar;
ylabel(cb,'RL tyre usage [-]','FontSize',12)
axis equal
xlabel('ay, [g]')
ylabel('ax, [g]')
ylim([-1.5 0.7])
title('GG with RL tyre usage')

figure()
scatter(sol.acc(:,2)/veh.g, sol.acc(:,1)/veh.g, 5, sol.mu(:,4)), hold on, grid minor
cb=colorbar;
ylabel(cb,'RR tyre usage [-]','FontSize',12)
axis equal
xlabel('ay, [g]')
ylabel('ax, [g]')
ylim([-1.5 0.7])
title('GG with RR tyre usage')


%%
figure()
plot(sol.t_opt,sol.x_opt(10,:))
hold on
plot(sol.t_opt,sol.x_opt(11,:))
hold on
plot(sol.t_opt,sol.x_opt(12,:))
hold on
plot(sol.t_opt,sol.x_opt(13,:))
grid on
ylabel('[rad/s]')
xlabel('Time, [s]')
legend('\omega_{fl}', '\omega_{fr}', '\omega_{rl}', '\omega_{ff}')
% xlim([0 s_end])
xlabel('Time, [s]')
title('Wheel rotation rate')


%%
figure()
plot(sol.t_opt,sol.x_opt(4,:))
hold on
plot(sol.t_opt,sol.x_opt(5,:))
hold on
plot(sol.t_opt,sol.x_opt(6,:))
hold on
plot(sol.t_opt,sol.x_opt(7,:))
grid on
ylabel('[rad or rad/s]')
xlabel('Time, [s]')
legend('\phi_{roll}', '\phi_{roll}_{dot}','\theta_{pitch}', '\theta_{pitch}_{dot}')
% xlim([0 s_end])
xlabel('Time, [s]')
title('Pitch and roll')



%% Controls

% re-compute motor rotation speed
omega_motor = (sol.x_opt(12,:)+sol.x_opt(13,:))/2 * veh.gr; % [rad/s]
% re-compute maximum available torque from motor
torque_motor_max = veh.power_motor_max ./ omega_motor; % [Nm]
% re-compute maximum avaialbe torque at the rear axle
torque_rear_axle_max = torque_motor_max * veh.gr; % [Nm]

% figure()
% plot(track.s,sol.u_opt(1,:)), grid on, hold on
% plot(track.s,sol.u_opt(2,:))
% ylabel('T_t/T_b [N-m]')
% xlabel('Distance, [m]')
% xlim([0 s_end])
% title('Total driving/braking torque')
figure()
plot(sol.t_opt,sol.u_opt(1,:)), grid on, hold on
plot(sol.t_opt,sol.u_opt(2,:))
plot(sol.t_opt,torque_rear_axle_max)
ylabel('T_t/T_b [N-m]')
xlabel('Time, [s]')
legend('Total drive torque', 'Total brake torque', 'Maximum available torque due to power limit')
% xlim([0 s_end])
title('Total driving/braking torque')

figure()
plot(sol.t_opt,rad2deg(sol.u_opt(3,:))), grid on,
ylabel('\delta [deg]')
xlabel('Time, [s]')
% xlim([0 s_end])
title('Front wheel(s) steering angle')


%% Output Variables
figure()
plot(sol.t_opt,sol.torque(:,1),'-'), grid on, hold on
plot(sol.t_opt,sol.torque(:,2),'--'),
plot(sol.t_opt,sol.torque(:,3),'-.'),
plot(sol.t_opt,sol.torque(:,4),':')
ylabel('T, [N-m]')
xlabel('Time, [s]')
% xlim([0 s_end])
legend('T_{fl}','T_{fr}','T_{rl}','T_{rr}')
xlabel('Time, [s]')
title('Torque on each wheel')



%% Trajectory

n_track = [-sin(track.Theta) cos(track.Theta)];
trackbnd_l = [track.x track.y] + track.wl.*n_track;
trackbnd_r = [track.x track.y] - track.wr.*n_track;

f4 = figure;
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,sol.x_opt(1,:)*3.6), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'V [km/h]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with speed overlay')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,rad2deg(atan(sol.x_opt(2,:)./sol.x_opt(1,:)))), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'\beta [deg]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with sideslip angle overlay')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,rad2deg(sol.u_opt(3,:))), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'\delta [deg]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with front wheel(s) steering angle overlay')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,sol.u_opt(1,:)), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'Tdrive [-]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with total drive torque overlay')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,sol.u_opt(2,:)), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'Tdbrake [-]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with total brake torque overlay')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% trajectory overlay with each tyre's friction usage
figure()
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,sol.mu(:,1)), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'FL tyre usage [-]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with FL tyre usage overlay')

figure()
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,sol.mu(:,2)), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'FR tyre usage [-]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with FR tyre usage overlay')

figure()
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,sol.mu(:,3)), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'RL tyre usage [-]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with RL tyre usage overlay')

figure()
plot(trackbnd_l(:,1),trackbnd_l(:,2),'-k'), grid on, hold on
plot(trackbnd_r(:,1),trackbnd_r(:,2),'-k')

scatter(sol.track_opt(:,1),sol.track_opt(:,2),5,sol.mu(:,4)), hold on

plot(track_opt(:,1),track_opt(:,2),'Color',[0.3 0.3 0.3])
cb=colorbar;
ylabel(cb,'RR tyre usage [-]','FontSize',12)

axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Trajectory with RR tyre usage overlay')






%%

figure()
plot(sol.t_opt,sol.tyre(:,9),'-'), grid on, hold on
plot(sol.t_opt,sol.tyre(:,10),'--'),
plot(sol.t_opt,sol.tyre(:,11),'-.'),
plot(sol.t_opt,sol.tyre(:,12),':')
ylabel('[N]')
% xlim([0 s_end])
legend('Fz_{fl}','Fz_{fr}','Fz_{rl}','Fz_{rr}')

figure()
plot(sol.t_opt,sol.tyre(:,1),'-'), grid on, hold on
plot(sol.t_opt,sol.tyre(:,2),'--'),
plot(sol.t_opt,sol.tyre(:,3),'-.'),
plot(sol.t_opt,sol.tyre(:,4),':')
ylabel('[N]')
% xlim([0 s_end])
legend('Fx_{fl}','Fx_{fr}','Fx_{rl}','Fx_{rr}')

figure()
plot(sol.t_opt,sol.tyre(:,5),'-'), grid on, hold on
plot(sol.t_opt,sol.tyre(:,6),'--'),
plot(sol.t_opt,sol.tyre(:,7),'-.'),
plot(sol.t_opt,sol.tyre(:,8),':')
ylabel('[N]')
% xlim([0 s_end])
legend('Fy_{fl}','Fy_{fr}','Fy_{rl}','Fy_{rr}')




figure()
plot(sol.t_opt,sol.slip(:,1),'-'), grid on, hold on
plot(sol.t_opt,sol.slip(:,2),'--'),
plot(sol.t_opt,sol.slip(:,3),'-.'),
plot(sol.t_opt,sol.slip(:,4),':')
ylabel('[-]')
% xlim([0 s_end])
legend('\lambda_{fl}','\lambda_{fr}','\lambda_{rl}','\lambda_{rr}')

figure()
plot((sol.t_opt).',sol.slip(:,5),'-'), grid on, hold on
plot((sol.t_opt).',sol.slip(:,6),'--'),
plot((sol.t_opt).',sol.slip(:,7),'-.'),
plot((sol.t_opt).',sol.slip(:,8),':')
ylabel('[rad]')
% xlim([0 s_end])
legend('\alpha_{fl}','\alpha_{fr}','\alpha_{rl}','\alpha_{rr}')


%% Control rate at knots
control_rate_T = (sol.u_opt(1, 2:end) - sol.u_opt(1, 1:end-1)) ./ (sol.t_opt(2:end) - sol.t_opt(1:end-1));
control_rate_Tb = (sol.u_opt(2, 2:end) - sol.u_opt(2, 1:end-1)) ./ (sol.t_opt(2:end) - sol.t_opt(1:end-1));
control_rate_delta = (sol.u_opt(3, 2:end) - sol.u_opt(3, 1:end-1)) ./ (sol.t_opt(2:end) - sol.t_opt(1:end-1));

figure();plot(sol.t_opt(2:end), control_rate_T)
ylabel('[Nm/s]')
title('Total drive torque rate')

figure();plot(sol.t_opt(2:end), control_rate_Tb)
ylabel('[Nm/s]')
title('Total brake torque rate')

figure();plot(sol.t_opt(2:end), control_rate_delta)
ylabel('[rad/s]')
title('Front wheel steering angle rate')


%%
debug = 0;