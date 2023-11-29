clear all;
close all;
clc;

%% ----- Enviroment data -----

% Universal constant
R_E = astroConstants(23); 
mu_E = astroConstants(13);
w_E = deg2rad(15.04)/60/60;

% Orbit satellite data
a = 20000 + R_E; %km
e = 0.001;
i = deg2rad(60); 
OM = deg2rad(40); 
om = deg2rad(20); 
th0 = deg2rad(0); 
kep = [a,e,i,OM,om,th0];

% Other orbit satellite data
p = a*(1-e^2);
n_e = sqrt(mu_E/a^3);
T = 2*pi/n_e;
DCM_eci_per = rotationMatrix(rad2deg(i),rad2deg(OM), rad2deg(om));
[r0, v0] = kep2carRAD(kep,mu_E);
car0 = [r0;v0];

% Sun data for ephemerides
date0 = [2023 12 25 0 0 0]; %date [year month day h m s]
MJD2000_initial = date2mjd2000(date0);
epsilon = deg2rad(23.45);

%% SIMULATION data

% Time definition
t0 = 0;
tf = 200;
t_step = 0.01;
ts_sensor = t_step;

% MAG-3 sensor data and error
FS = 800e-9; %Full scale value
inc_lin = 0.15/100; % Accuracy error
inc_acc = 0.75/100; % Linear error
inc_tot = sqrt(inc_lin^2 + inc_acc^2); %Total error

% Definition of Magnetic Field (Dipole model)
j_b = [0.01 0.05 0.01]';
H0 = 1e-9.*[-29404.8 -1450.9 4652.5]';

% ----------- Initialization of S/C ------------------
% Principal moments of inertia
Ix = 0.040; %kg m^2
Iy = 0.060; %kg m^2
Iz = 0.10; %kg m^2
I =  [Ix, 0, 0; 0, Iy, 0; 0, 0, Iz];
I_inv = inv(I);

% Initial Angular Velocity Conditions
wx0 = 1e-4; %rad/s
wy0 = 1e-4; %rad/s
wz0 = 0.055; %rad/s
w_r = 0; %rad/s
w0 = [wx0 wy0 wz0]'; %vector of Initial Angular Velocity

% Initial quaternion position
q0 = [0 0 0 1];

% Body configuration
BODY = [1, 0, 0;
        0, 1, 0;
        -1,0, 0;
        0, -1, 0;
        0, 0, +1;
        0, 0, -1;
        1, 0, 0;
        -1, 0, 0;
        +1, 0, 0;
        -1, 0, 0];

r = 1e-2*[10, 0, 0;
    0, 10, 0;
    -10,0, 0;
    0, -10, 0;
    0, 0, 15;
    0, 0, -15;
    0, 0, 45;
    0, 0, 45;
    0, 0, -45;
    0, 0, -45];

CMC_position = [0,0,0]; % applied rigid shift

% Solar pressure data
r1 = r(:,1) + CMC_position(1).*ones(10,1);
r2 = r(:,2) + CMC_position(2).*ones(10,1);
r3 = r(:,3)+ CMC_position(3).*ones(10,1);
A1 = ones(length(BODY),1);
A = 1e-2.*[6, 6, 6, 6, 4, 4, 12, 12, 12, 12]';
ps = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1]';
pd = 0.1.*A1;
Fe = 1358;
c = 3e8;

%% RUN SIMULATION
out = sim(['SAD_project.slx']);

%% PLOTTING RESULT

%%% DISTURBANCE PLOT %%%

% GRAVITY GRADIENT disturbance
fig1 = figure('WindowState', 'maximized');
tiledlayout(2,3);
nexttile
hold on;
grid on;
plot(out.t, out.M_GG(:,1), '-b', 'LineWidth', 1);
plot(out.t, out.M_GG(:,2), '-g', 'LineWidth', 1);
plot(out.t, out.M_GG(:,3), '-r', 'LineWidth', 1);


xlabel(' time [s] ')
ylabel(' M ')
title('GRAVITY GRADIENT')

legend('M_{GG}x', 'M_{GG}y', 'M_{GG}z');

% SOLAR PRESSURE disturbance
nexttile
hold on;
grid on;
plot(out.t, out.M_SRP(:,1), '-b', 'LineWidth', 1);
plot(out.t, out.M_SRP(:,2), '-g', 'LineWidth', 1);
plot(out.t, out.M_SRP(:,3), '-r', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' M [N]')
title('SOLAR PRESSURE')

legend('M_{SRP}x', 'M_{SRP}y', 'M_{SRP}z');

% MAGNETIC FIELD
% History of disturbance
nexttile
hold on;
grid on;
plot(out.t, out.M_b(:,1), '-b', 'LineWidth', 1);
plot(out.t, out.M_b(:,2), '-g', 'LineWidth', 1);
plot(out.t, out.M_b(:,3), '-r', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' M [N]')
title('MAGNETIC FIELD')

legend('M_{b}x', 'M_{b}y', 'M_{b}z');

% Comparison of all disturbance
nexttile([1 3])
hold on;
grid on;
semilogy(out.t, vecnorm(out.M_GG,2,2), '-b', 'LineWidth', 1);
semilogy(out.t, vecnorm(out.M_SRP,2,2), '-g', 'LineWidth', 1);
semilogy(out.t, vecnorm(out.M_b,2,2), '-r', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' M [N]')
title('History of disturbance')

legend('M_{GG}', 'M_{SRP}', 'M_{b}');

%%% OTHER PLOT %%%

% History of Magnetic Field B(r)
fig2 = figure;
hold on;
grid on;
plot(out.t, 1e9.*out.bn(:,1), '-b', 'LineWidth', 1);
plot(out.t, 1e9.*out.bn(:,2), '-g', 'LineWidth', 1);
plot(out.t, 1e9.*out.bn(:,3), '-r', 'LineWidth', 1);
plot(out.t, 1e9.*vecnorm(out.bn,2,2), '-r', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' B [nT]')
title('History of Magnetic Field B(r)')

legend('Bx', 'By', 'Bz', 'Norm(B)');

% ATTITUDE ERROR
% plot figure of w
fig3 = figure();
wx_plot = plot(out.t, out.wb(:,1), '-b', 'LineWidth', 1);
hold on
grid on
wy_plot = plot(out.t, out.wb(:,2), '-c', 'LineWidth', 1);
wz_plot = plot(out.t, out.wb(:,3), '-g', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' w [rad/s] ')
title('Attitude Error')

legend([wx_plot, wy_plot, wz_plot], 'wx', 'wy', 'wz');

% MAG-3 SENSOR
fig4 = figure;
hold on;
grid on;
plot(out.t, 1e9.*out.b_sensor(:,1), '-b', 'LineWidth', 1);
plot(out.t, 1e9.*out.b_sensor(:,2), '-g', 'LineWidth', 1);
plot(out.t, 1e9.*out.b_sensor(:,3), '-r', 'LineWidth', 1);
plot(out.t, 1e9.*vecnorm(out.b_sensor,2,2), '-c', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' B [nT]')
title('Magnetic Field B(r) by sensor')

legend('Bx', 'By', 'Bz', 'Norm(B)');

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
%% ANIMATION PLOT

% ancora da sistemare 
%{
    ground Track
    plot Orbit
    Animation of satellite
%}