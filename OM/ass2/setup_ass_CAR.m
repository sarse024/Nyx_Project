clear;
close all;
clc;
 

% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
mu_E = astroConstants(13); %[km^3/s^2]
R_E = astroConstants(23);
w_E = deg2rad(15.04)/60/60;
tetha_G0 = 0;
N = 1000;

% Set your initial time and state (Cartesian or Keplerian): t0, s0
a = 8350; % Km
e = 0.8;
i = deg2rad(60); % deg
OM = deg2rad(270); %deg
w = deg2rad(45); %deg
f = deg2rad(230); %deg
%T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]
kep0 = [a, e, i, OM, w, f];

[r,v] = kep2carRAD(kep0, mu_E);

s0 = [r;v];

% Set your integration time span: tspan
tF = 20*60*60; % 
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
parameters.j2 = 0.00108263;
parameters.kep = kep0;

% Set ODE solver options (e.g., RelTol, AbsTol): options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Numerical integration of the equations of motion
%[ T, S ] = ode113( @(t,s) eq_motion( t, s, @(t,s) acc_pert_fun(t,s,parameters), parameters), tspan, s0, options);

[T, S ] = ode113( @(t,s) eq_motion_CAR( t, s, @(t,s) acc_pert_fun_CAR(t,s,parameters), parameters ), tspan, s0, options );
% Analyse and plot the results

figure()
Terra3d
plot3( S(:,1), S(:,2), S(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;
%%
kep_matrix = zeros(length(S), 6);
for i = 1:length(S) 
    kep = car2kepRAD(S(i,1:3), S(i,4:6), mu_E);
    kep_matrix(i,:) = kep;
end

[dOM_sec, dw_sec, dM_sec] = j2_secular();

figure()
tiledlayout(2,3);
title('Keplerian change by perturbation')

% a plot
nexttile
plot( tspan./3600, kep_matrix(:,1), '-' )
xlabel('time [H]');
ylabel('a [km]');
title('a plot');
grid on;

% e plot
nexttile
plot( tspan./3600, kep_matrix(:,2), '-' )
xlabel('time [H]');
ylabel('e [-]');
title('e plot');
grid on;

% i plot
nexttile
plot( tspan./3600, rad2deg(kep_matrix(:,3)), '-' )
xlabel('time [H]');
ylabel('i [deg]');
title('i plot');
grid on;

% OM plot
nexttile
plot( tspan./3600, rad2deg(kep_matrix(:,4)), '-' )
xlabel('time [H]');
ylabel('OM [deg]');
title('OM plot');
grid on;


% om plot
nexttile
plot( tspan./3600, rad2deg(kep_matrix(:,5)), '-' )
xlabel('time [H]');
ylabel('om [deg]');
title('om plot');
grid on;


% th plot
nexttile
plot( tspan./3600, rad2deg(kep_matrix(:,6)), '-' )
xlabel('time [H]');
ylabel('f [deg]');
title('f plot');
grid on;

% potrebbe essere interessante uno studio sul deorbiting 
% movie of the orbit change