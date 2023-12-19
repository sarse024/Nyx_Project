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
e = 0.19760;
i = deg2rad(60); % deg
OM = deg2rad(270); %deg
w = deg2rad(45); %deg
f = deg2rad(230); %deg
%T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]

s0 = [a, e, i, OM, w, f];


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
%parameters.j2 = 0.00108263;

% Set ODE solver options (e.g., RelTol, AbsTol): options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Numerical integration of the equations of motion
[T, S ] = ode113( @(t,s) eq_motion_GAUSS( t, s, @(t,s) acc_pert_fun_RWS(t,s,parameters), parameters ), tspan, s0, options );

% Analyse and plot the results
figure()
tiledlayout(2,3);
title('Keplerian change by perturbation')

% a plot
nexttile
plot( tspan./3600, S(:,1), '-' )
xlabel('time [H]');
ylabel('a [km]');
title('a plot');
grid on;

% e plot
nexttile
plot( tspan./3600, S(:,2), '-' )
xlabel('time [H]');
ylabel('e [-]');
title('e plot');
grid on;

% i plot
nexttile
plot( tspan./3600, rad2deg(S(:,3)), '-' )
xlabel('time [H]');
ylabel('i [deg]');
title('i plot');
grid on;

% OM plot
nexttile
plot( tspan./3600, rad2deg(S(:,4)), '-' )
xlabel('time [H]');
ylabel('OM [deg]');
title('OM plot');
grid on;


% om plot
nexttile
plot( tspan./3600, rad2deg(S(:,5)), '-' )
xlabel('time [H]');
ylabel('om [deg]');
title('om plot');
grid on;


% th plot
nexttile
plot( tspan./3600, rad2deg(S(:,6)), '-' )
xlabel('time [H]');
ylabel('f [deg]');
title('f plot');
grid on;

% potrebbe essere interessante uno studio sul deorbiting 
% movie of the orbit change


