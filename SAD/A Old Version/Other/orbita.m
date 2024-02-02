clear
clc

mu_E = astroConstants(13); % [km^3/s^2]
R_E = astroConstants(23); % [km]
mu_S = astroConstants(4); % [km^3/s^2]

% Satellite orbit parameters
h = 1500; % [km]
a = R_E + h; % [km]
e = 0; 
T = 2*pi*sqrt((a^3)/mu_E); % [s]
n = (2*pi)/T; % [1/s]
i = deg2rad(30); % [rad]
OMEGA = deg2rad(45); % [rad]
omega = deg2rad(0); % [rad]
theta_0 = deg2rad(230); % [rad]

[r0, v0] = kep2car(a, e, i, OMEGA, omega, theta_0, mu_E);
%r0 = r0.*10^3;
%v0 = v0.*10^3;

%%
% Earth's orbit around the Sun
a_E = 149598023; % [km]
e_E = 0.0167086; 
T_E = 2*pi*sqrt((a_E^3)/mu_S); % [s]
n_E = (2*pi)/T_E; % [1/s]
i_E = deg2rad(7.155); % [rad]
OMEGA_E = deg2rad(-11.26064); % [rad]
omega_E = deg2rad(114.20783); % [rad]
theta_0_E = deg2rad(230); % [rad]

[r0_E, v0_E] = kep2car(a_E, e_E, i_E, OMEGA_E, omega_E, theta_0_E, mu_S);

%%
% Satellite data
Ix = 0.0700; %[kg*m^2]
Iy = 0.0550; %[kg*m^2]
Iz = 0.025; %[kg*m^2]
Ir = 0; %[kg*m^2]
I = [Ix 0 0; 0 Iy 0; 0 0 Iz];
I_inv = inv(I);
omega_0 = [0.45 0.52 0.55]'; %[rad/s]
A_0 = [1 0 0; 0 1 0; 0 0 1];
omega_r0 = 0;
rot = [0 0 Ir*omega_r0]';
M = 0;

%% 
% Orbit plot
out = sim('Prova_orbita.slx');
plot3(out.r_N(:,1), out.r_N(:,2), out.r_N(:,3))
xlabel('X [km]'); 
ylabel('Y [km]'); 
zlabel('Z [km]');
title('Satellite orbit around Earth');

%% 
% Earth's orbit plot
figure()
plot3(out.S_N(:,1), out.S_N(:,2), out.S_N(:,3))
xlabel('X [km]'); 
ylabel('Y [km]'); 
zlabel('Z [km]');
title('Earth orbit around Sun');
