clear all;
close all;
clc;

%% ----- Enviroment data -----

% Universal constant
R_E = astroConstants(23); 
mu_E = astroConstants(13);
w_E = deg2rad(15.04)/60/60;

% Set up for better plot
setGraphicPlot;

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
timeoflaunch = date2mjd2000(date0);
epsilon = deg2rad(23.45);
InitialTrueAnomalyOfEarth = 0;

%% SIMULATION data

% Time definition
t0 = 0;
tf = 200;
t_step = 0.1; %maybe change
ts_sensor = t_step*10;

% MAG-3 sensor data and error
FS = 50000e-9; %Full scale value
inc_lin = 0.15/100; % Accuracy error
inc_acc = 0.75/100; % Linear error
inc_tot = sqrt(inc_lin^2 + inc_acc^2); %Total error

% weigh for statistical method
alpha = [0.2 0.8];

% Definition of Magnetic Field (Dipole model)
j_b = [0.01 0.05 0.01]';


% Definition of Magnetic Field (Order model)
A1 = importdata("igrfSg.txt");
A2 = importdata("igrfSh.txt");

% export data for dipole model
g_10 = A1(1,3) + A1(1,4)*(timeoflaunch - 7305)/365;
g_11 = A1(2,3) + A1(2,4)*(timeoflaunch - 7305)/365;
h_11 = A2(2,3) + A2(2,4)*(timeoflaunch - 7305)/365;

H0 = [g_10 g_11 h_11]';

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
N_1 = BODY(1,:);
N_2 = BODY(2,:);
N_3 = BODY(3,:);
N_4 = BODY(4,:);
N_5 = BODY(5,:);
N_6 = BODY(6,:);

% point of application of F SRP on each body face (DA SISTEMARE)
r = 0.5e-2*[33, 0, 0;
    0, 33, 0;
    -33,0, 0;
    0, -33, 0;
    0, 0, 15;
    0, 0, -15];
r_1 = r(1,:);
r_2 = r(2,:);
r_3 = r(3,:);
r_4 = r(4,:);
r_5 = r(5,:);
r_6 = r(6,:);


CMC_position = [0,0,15]*1e-2; % applied rigid shift

% Solar pressure data
r1 = r(:,1) + CMC_position(1).*ones(6,1);
r2 = r(:,2) + CMC_position(2).*ones(6,1);
r3 = r(:,3)+ CMC_position(3).*ones(6,1);
one_vector_body = ones(length(BODY),1);
area_face = [h_dim*l_dim, h_dim*l_dim, h_dim*l_dim, h_dim*l_dim, l_dim^2, l_dim^2]';
A_1 = area_face(1);
A_2 = area_face(2);
A_3 = area_face(3);
A_4 = area_face(4);
A_5 = area_face(5);
A_6 = area_face(6);
rho_s = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]';

rho_s1 = rho_s(1);
rho_s2 = rho_s(2);
rho_s3 = rho_s(3);
rho_s4 = rho_s(4);
rho_s5 = rho_s(5);
rho_s6 = rho_s(6);

rho_d1 = 0.1;
rho_d2 = 0.1;
rho_d3 = 0.1;
rho_d4 = 0.1;
rho_d5 = 0.1;
rho_d6 = 0.1;

Fe = 1358;
c = 3e8; % [m/s]

% Valuation of max disturbance
Mgg_MAX = 3/2 * mu_E / a^3 *abs(Ix - Iz)
Msrp = rho_s(1)*area_face(1)*(1+0.1)*(h_dim/4)
%j_b = 4*10^(-3).*[area_face(1) area_face(2) area_face(end)]./mass'; %estimate from Nasa table 




%% RUN SIMULATION
tic
out = sim(['SAD_project.slx']);
toc
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
axis tight

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
plot(out.t, 1e9.*out.bn(:,1), '-b', 'LineWidth', 2);
plot(out.t, 1e9.*out.bn(:,2), '-g', 'LineWidth', 2);
plot(out.t, 1e9.*out.bn(:,3), '-r', 'LineWidth', 2);
plot(out.t, 1e9.*vecnorm(out.bn,2,2), '-c', 'LineWidth', 2);

xlabel(' time [s] ')
ylabel(' B [nT]')
title('History of Magnetic Field B(r)')

legend('Bx', 'By', 'Bz', 'Norm(B)');
%%
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
title('TRUE Angulary Velocity in time')

legend([wx_plot, wy_plot, wz_plot], 'wx', 'wy', 'wz');
%%
% reaction wheel
% plot figure of RW angular velocity
fig6 = figure();
tiledlayout(1,2);
nexttile    
wx_plot = plot(out.t, out.wr(:,1), '-b', 'LineWidth', 2);
hold on
grid on
wy_plot = plot(out.t, out.wr(:,2), '-c', 'LineWidth', 2);
wz_plot = plot(out.t, out.wr(:,3), '-g', 'LineWidth', 2);
plot(out.t, wMax_RW*ones(length(out.t)),'-r', 'LineWidth', 1);
w_lim = plot(out.t, -wMax_RW*ones(length(out.t)),'-r', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' w [rad/s] ')
title('Reaction Wheel Angulary Velocity in time')

legend([wx_plot, wy_plot, wz_plot, w_lim(1)], 'wx', 'wy', 'wz', 'Saturation limit');
%
% plot figure Torque due to RW
nexttile    
Mx_plot = plot(out.t, out.M_actuator(:,1), '-b', 'LineWidth', 1);
hold on
grid on
My_plot = plot(out.t, out.M_actuator(:,2), '-c', 'LineWidth', 1);
Mz_plot = plot(out.t, out.M_actuator(:,3), '-g', 'LineWidth', 1);

xlabel(' time [s] ')
ylabel(' M [Nm] ')
title('Torque due to RW')

legend([Mx_plot, My_plot, Mz_plot], 'Mx', 'My', 'Mz');
%% ANIMATION PLOT

% ancora da sistemare 
%{
    ground Track
    plot Orbit
    Animation of satellite
%}