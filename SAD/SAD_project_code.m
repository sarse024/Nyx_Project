clear all;
close all;
clc;

%% ----- Enviroment data -----

% Universal constant
R_E = astroConstants(23); 
mu_E = astroConstants(13);
w_E = deg2rad(15.04)/60/60;

% Orbit satellite data from space-track
%TLE 1 43792U 18099AL  23355.11707279  .00003582  00000-0  30335-3 0  9999
%TLE 2 43792  97.5473  50.7464 0011241 206.3209 153.7444 14.98567371275372
a = 575 + R_E; %km
e = 0.0011241;
i = deg2rad(97.55); 
OM = deg2rad(50.7464); 
om = deg2rad(206.3209); 
th0 = deg2rad(0); 
kep = [a,e,i,OM,om,th0];

% Other orbit satellite data
p = a*(1-e^2);
n_e = sqrt(mu_E/a^3);
period = 2*pi/n_e;
DCM_eci_per = rotationMatrix(rad2deg(i),rad2deg(OM), rad2deg(om));
[r0, v0] = kep2carRAD(kep,mu_E);
car0 = [r0;v0];

% Sun data for ephemerides DA CAMBIARE IN BASE AL BLOCCO !!!!!!
date0 = [2023 12 25 0 0 0]; %date [year month day h m s]
MJD2000_initial = date2mjd2000(date0);
epsilon = deg2rad(23.45);

%% SIMULATION data

% Time definition
t0 = 0;
tf = period;
t_step = 0.01; %maybe change
ts_sensor = t_step*10;

% MAG-3 sensor data and error
FS = 800e-9; %Full scale value
inc_lin = 0.15/100; % Accuracy error
inc_acc = 0.75/100; % Linear error
inc_tot = sqrt(inc_lin^2 + inc_acc^2); %Total error

% weigh for statistical method
alpha = [0.2 0.8];

% Definition of Magnetic Field (Dipole model)
j_b = [0.01 0.05 0.01]';
H0 = 1e-9.*[-29404.8 -1450.9 4652.5]';

% gyro sensor (Da sistmare)
bias = 0.05; %[deg/h]
bias = bias*(pi/180)*(1/3600); %[rad/s]
bias = bias*[1 1 1]';
ts = 0.1;
ARW = 0.007; %[deg/sqrt(h)]
ARW = ARW*(pi/180)*(1/60)*(1/sqrt(ts)); %[rad/sqrt(h)]
ARW = ARW*[1 1 1]';

% Horizon sensor 
horizon_cone_limit = deg2rad(67);

% Reaction Wheel DATA (Da sistemare)
I_RW = 1*10^(-5); %[kg*m^2]
I_vector_RW  = [I_RW I_RW I_RW]';
omega_0 = [0 0 0]'; %[rad/s]
matrix_RW_config = eye(3);
matrix_RW_config_inv = inv(matrix_RW_config);
w0_RW = [0 0 0];
h_r_0 = [0,0,0]';
RPM_max_RW = 5035;
wMax_RW = 2*pi * RPM_max_RW / 60;

% 

% ----------- Initialization of S/C ------------------
% Principal moments of inertia
h_dim = 0.66; %[m] heigth dimension
l_dim = 0.33; %[m] base dimension of square
mass = 50; %[Kg]
Ix = 4.350; %[kg*m^2]
Iy = 4.3370 ; %[kg*m^2]
Iz = 3.6640; %[kg*m^2]
I =  5/12.*diag([Ix Iy Iz]);
I_inv = inv(I);

% Initial Angular Velocity Conditions
wx0 = 1e-4; %rad/s
wy0 = 1e-4; %rad/s
wz0 = 0.055; %rad/s
w0 = [wx0 wy0 wz0]'; %vector of Initial Angular Velocity

% Initial quaternion position
q0 = [0 0 0 1];

% Body configuration 66x33x33 cm mass = 50 kg
BODY = [1, 0, 0;
        0, 1, 0;
        -1,0, 0;
        0, -1, 0;
        0, 0, +1;
        0, 0, -1];

% point of application of F SRP on each body face (DA SISTEMARE)
r = 1e-2*[10, 0, 0;
    0, 10, 0;
    -10,0, 0;
    0, -10, 0;
    0, 0, 15;
    0, 0, -15];

CMC_position = [0,0,0]; % applied rigid shift

% Solar pressure data
r1 = r(:,1) + CMC_position(1).*ones(6,1);
r2 = r(:,2) + CMC_position(2).*ones(6,1);
r3 = r(:,3)+ CMC_position(3).*ones(6,1);
one_vector_body = ones(length(BODY),1);
area_face = [h_dim*l_dim, h_dim*l_dim, h_dim*l_dim, h_dim*l_dim, l_dim^2, l_dim^2]';
ps = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]';
pd = 0.1.*one_vector_body;
Fe = 1358;
c = 3e8; % [m/s]

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