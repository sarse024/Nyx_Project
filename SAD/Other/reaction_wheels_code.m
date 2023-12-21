clear
clc
close all

% Reaction wheels - pyramid configuration
a = 1/sqrt(3);
A = [-a a a -a; -a -a a a; a a a a];
A = eye(3);
b = sqrt(3)/4;
A_ = [-b -b b; b -b b; b b b; -b b b];

h_r_0 = [0 0 0 0]';
M_c = [1 0 0]';

%% 3 reaction wheel

t0 = 0;
tf = 200;
t_step = 0.1;
A = eye(3);
A_inv = inv(A);
h_r_0 = [0 0 0]';
M_c = [1 0 0]';

%%
% Satellite data
Ix = 4.350; %[kg*m^2]
Iy = 4.3370 ; %[kg*m^2]
Iz = 3.6640; %[kg*m^2]
I_RW = 2*10^(-3); %[kg*m^2]
I = [Ix 0 0; 0 Iy 0; 0 0 Iz];
I_vector_RW  = [I_RW I_RW I_RW]';
I_inv = inv(I);
omega_0 = [0 0 0]'; %[rad/s]
w0_RW = [0 0 0];
RPM_max_RW = 5035;
wMax_RW = 2*pi * RPM_max_RW / 60;
M = 0;

% maximum torque m = 0.4 Nm
%The moment of inertia for the wheel is 4 · 10−5 kgm2
% saturation of Reaction Wheel 

tic
out = sim('reaction_wheel.slx');
toc

% plot figure of w
fig5 = figure();
wx_plot = plot(out.t, out.w(:,1), '-b', 'LineWidth', 1);
hold on
grid on
wy_plot = plot(out.t, out.w(:,2), '-c', 'LineWidth', 1);
wz_plot = plot(out.t, out.w(:,3), '-g', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' w [rad/s] ')
title('Angulary Velocity in time')

legend([wx_plot, wy_plot, wz_plot], 'wx', 'wy', 'wz');

% plot figure of RW angular velocity
fig5 = figure();
wx_plot = plot(out.t, out.wr(:,1), '-b', 'LineWidth', 2);
hold on
grid on
wy_plot = plot(out.t, out.wr(:,2), '-c', 'LineWidth', 2);
wz_plot = plot(out.t, out.wr(:,3), '-g', 'LineWidth', 2);
w_lim = plot(out.t, wMax_RW*ones(length(out.t)),'-r', 'LineWidth', 1);
w_lim = plot(out.t, -wMax_RW*ones(length(out.t)),'-r', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' w [rad/s] ')
title('Reaction Wheel Angulary Velocity in time')

legend([wx_plot, wy_plot, wz_plot, w_lim], 'wx', 'wy', 'wz', 'Saturation limit');
%%
% plot figure Torque due to RW
fig5 = figure();
wx_plot = plot(out.t, out.M_actuator(:,1), '-b', 'LineWidth', 1);
hold on
grid on
wy_plot = plot(out.t, out.M_actuator(:,2), '-c', 'LineWidth', 1);
wz_plot = plot(out.t, out.M_actuator(:,3), '-g', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' M [Nm] ')
title('Torque due to RW')

legend([wx_plot, wy_plot, wz_plot], 'wx', 'wy', 'wz');

