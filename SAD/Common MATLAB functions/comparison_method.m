clear;
close all;
clc;

%% CAR METHOD

% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
mu_E = astroConstants(13); %[km^3/s^2]
R_E = astroConstants(23);
w_E = deg2rad(15.04)/60/60;
tetha_G0 = 0;
N = 10000;

% Set your initial time and state (Cartesian or Keplerian): t0, s0
e = 0.4;

a_deorbit = @(a) a*(1-e) - (250 + R_E);
a = fzero(a_deorbit, 10000);

i = deg2rad(87.9); % deg
OM = deg2rad(180); %deg
w = deg2rad(180); %deg
f = deg2rad(0); %deg
T_period = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]
kep0 = [a, e, i, OM, w, f];

[r,v] = kep2carRAD(kep0, mu_E);

s0 = [r;v];

% Set your integration time span: tspan
tF = 10000*T_period; % 
tspan = linspace(0, tF, N);

% parameters
parameters.rE = R_E;              %[Km]
parameters.wE = w_E;              % rad/s;
parameters.mu = mu_E;
parameters.drag.CD = 2.1;         %[]
parameters.drag.AM = 0.0043;      %[m^2/kg]

%parameters.drag.CD = 0;         %[]
%parameters.drag.AM = 0;      %[m^2/kg]
parameters.drag.rE = R_E;         %[Km]
parameters.j2 = 0;
parameters.j2 = astroConstants(9);
parameters.kep = kep0;

% Set ODE solver options (e.g., RelTol, AbsTol): options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

[T, S] = ode113( @(t,s) eq_motion_CAR( t, s, @(t,s) acc_pert_fun_CAR(t,s,parameters), parameters ), tspan, s0, options);
% Analyse and plot the results

figure()
Terra3d
plot3( S(:,1), S(:,2), S(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

kep_matrix = zeros(length(S), 6);

for i = 1:length(S) 
    kep = car2kepRAD(S(i,1:3)', S(i,4:6)', mu_E);
    kep_matrix(i,:) = kep;
end
altitude_CAR =  kep_matrix(:,1).*((ones(length(S),1) - kep_matrix(:,2)))- R_E*ones(length(S),1);

%% GAUSS METHOD

%state vector
s0 = kep0;

% Numerical integration of the equations of motion
[T, S] = ode113( @(t,s) eq_motion_GAUSS( t, s, @(t,s) acc_pert_fun_RWS(t,s,parameters), parameters ), tspan, s0, options );

%% figure and comparison

%adjust tspan
tspan = tspan./T_period;
figure();
tiledlayout(2,3);
title('Keplerian change by perturbation')

% a plot
nexttile
hold on;
grid on;
plot( tspan, kep_matrix(:,1), '-b' )
plot( tspan, S(:,1), '-r' )
xlabel('time [T]');
ylabel('a [km]');
title('a plot');


% e plot
nexttile
hold on;
grid on;
plot( tspan, kep_matrix(:,2), '-b' )
plot( tspan, S(:,2), '-r' )
xlabel('time [T]');
ylabel('e [-]');
title('e plot');
grid on;

% i plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,3)), '-' )
plot( tspan, rad2deg(S(:,3)), '-r' )
xlabel('time [T]');
ylabel('i [deg]');
title('i plot');
grid on;

% OM plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,4)), '-b' )
plot( tspan, rad2deg(S(:,4)), '-r' )
xlabel('time [T]');
ylabel('OM [deg]');
title('OM plot');
grid on;


% om plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,5)), '-b' )
plot( tspan, rad2deg(S(:,5)), '-r' )
xlabel('time [T]');
ylabel('om [deg]');
title('om plot');
grid on;


% th plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,6)), '-' )
plot( tspan, mod(rad2deg(S(:,6)), 360), '-r' )
xlabel('time [T]');
ylabel('f [deg]');
title('f plot');
legend

%% comparison error
figure()
grid on;
hold on;
err = (kep_matrix(:,1:5)-S(:,1:5));
semilogy( tspan, err, '-' )
legend('a','e','i','OM','om')
xlabel('time [T]');
ylabel('err [-]');
title('Error from GAUSS and CAR method')

figure()
grid on
hold on
plot(tspan, altitude_CAR)
xlabel('time [T]');
ylabel('r_p [Km]');
title('Altitude of pericenter')


%% figure of lab 5
index = fix(length(S)/10);

% a plot
figure();
tiledlayout(1,3);
nexttile
err = abs((kep_matrix(:,1)-S(:,1)))./kep0(1);
semilogy( tspan, err, '-b' )
grid on;
xlabel('time [T]');
ylabel('err [-]');
%ylabel('err = $\frac{\abs{a_{car} - a_{gauss}}}{a_0} $', 'interpreter','latex');
title('Relative error of a');

nexttile
hold on;
grid on;
plot( tspan, kep_matrix(:,1), '-b' )
plot( tspan, S(:,1), '-r' )
xlabel('time [T]');
ylabel('a [km]');
title('a plot');
legend('Cartesian', 'Gauss')

nexttile
hold on;
grid on;
plot( tspan(1:index), kep_matrix(1:index,1), '-b' )
plot( tspan(1:index), S(1:index,1), '-r' )
xlabel('time [T]');
ylabel('a [km]');
title('10 orbit study');
legend('Cartesian', 'Gauss')

% e plot
figure();
tiledlayout(1,3);
nexttile
err = abs((kep_matrix(:,2)-S(:,2)));
semilogy( tspan, err, '-b' )
grid on;
xlabel('time [T]');
ylabel('err [-]');
title('Relative error of e');

nexttile
hold on;
grid on;
plot( tspan, kep_matrix(:,2), '-b' )
plot( tspan, S(:,2), '-r' )
xlabel('time [T]');
ylabel('e [-]');
title('e plot');
legend('Cartesian', 'Gauss')

nexttile
hold on;
grid on;
plot( tspan(1:index), kep_matrix(1:index,2), '-b' )
plot( tspan(1:index), S(1:index,2), '-r' )
xlabel('time [T]');
ylabel('e [-]');
title('10 orbit study');
legend('Cartesian', 'Gauss')

% i plot
figure();
tiledlayout(1,3);
nexttile
err = abs((kep_matrix(:,3)-S(:,3)))./(2*pi);
semilogy( tspan, err, '-b' )
grid on;
xlabel('time [T]');
ylabel('err [-]');
title('Relative error of i');

nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,3)), '-b' )
plot( tspan, rad2deg(S(:,3)), '-r' )
xlabel('time [T]');
ylabel('i [deg]');
title('i plot');
legend('Cartesian', 'Gauss')

nexttile
hold on;
grid on;
plot( tspan(1:index), rad2deg(kep_matrix(1:index,3)), '-b' )
plot( tspan(1:index), rad2deg(S(1:index,3)), '-r' )
xlabel('time [T]');
ylabel('i [deg]');
title('10 orbit study');

% OM plot
figure();
tiledlayout(1,3);
nexttile
err = abs((kep_matrix(:,4)-S(:,4)))./(2*pi);
semilogy( tspan, err, '-b' )
grid on;
xlabel('time [T]');
ylabel('err [-]');
title('Relative error of \Omega');

nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,4)), '-b' )
plot( tspan, rad2deg(S(:,4)), '-r' )
xlabel('time [T]');
ylabel('\Omega [deg]');
title('\Omega plot');

nexttile
hold on;
grid on;
plot( tspan(1:index), rad2deg(kep_matrix(1:index,4)), '-b' )
plot( tspan(1:index), rad2deg(S(1:index,4)), '-r' )
xlabel('time [T]');
ylabel('\Omega [deg]');
title('10 orbit study');

% om plot
figure();
tiledlayout(1,3);
nexttile
err = abs((kep_matrix(:,5)-S(:,5)))./(2*pi);
semilogy( tspan, err, '-b' )
grid on;
xlabel('time [T]');
ylabel('err [-]');
title('Relative error of \omega');

nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,5)), '-b' )
plot( tspan, rad2deg(S(:,5)), '-r' )
xlabel('time [T]');
ylabel('\omega [deg]');
title('\omega plot');

nexttile
hold on;
grid on;
plot( tspan(1:index), rad2deg(kep_matrix(1:index,5)), '-b' )
plot( tspan(1:index), rad2deg(S(1:index,5)), '-r' )
xlabel('time [T]');
ylabel('om [deg]');
title('10 orbit study');

% th plot
figure();
tiledlayout(1,3);
nexttile
err = abs( mod(kep_matrix(:,6), 2*pi) - S(:,6)) ./ (2*pi);
semilogy( tspan, err, '-b' )
grid on;
xlabel('time [T]');
ylabel('err [-]');
title('Relative error of \theta');

nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,6)), '-b' )
plot( tspan, rad2deg(S(:,6)), '-r' )
xlabel('time [T]');
ylabel('\theta [deg]');
title('om plot');

nexttile
hold on;
grid on;
plot( tspan(1:index), rad2deg(kep_matrix(1:index,6)), '-b' )
plot( tspan(1:index), rad2deg(S(1:index,6)), '-r' )
xlabel('time [T]');
ylabel('th [deg]');
title('10 orbit study');
