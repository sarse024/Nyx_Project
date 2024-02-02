clear all;
close all;
clc;

%%%% OM PROJECT -> Assignments 2: : Planetary Explorer Mission %%%%%

% Data from File
% Group ID:2336 
% a [10e4 km]: 0.6846
% e [-]: 0.0298
% i [deg]: 80.2068

% Repeating GT ratio k:m: 15:1
% Perturbations: J2 DRAG
% Parameters: cD = 2.1 A/M = 0.0043 m^2/kg

%nominal orbit data
a = 6846;         %[km]
e = 0.0298;       %[]
i = 80.2068;      %[deg]
k = 15;           %[]
m = 1;            %[]

%drag parameters
CD = 2.1;         %[]
AM = 0.0043;      %[m^2/kg]

%other keplerian elements (arbitrary)
OM = 285; %[deg] %dende lab OM = 284.67
om = 135; %[deg] % dende lab om = 135.52
th = 0.01; %[deg] % dende lab th = 297.23
kep0 = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th)];

%% DATA FOR BOTH PERTURBATION

%pert
mu = astroConstants(13);
rE = astroConstants(23);
j2 = astroConstants(9);
t_sid = 23*60*60 + 56*60 +4;
wE = 2*pi/t_sid;          %[rad/s] 

% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
N = 10000;
N_orbit = 2; 
N_filter1 = N/5;
N_filter2 = N/100;

T_period = 2*pi*sqrt( kep0(1)^3/mu); % Orbital period [s]

% Calculate and save initial condition
[r,v] = kep2carRAD(kep0, mu);
s0 = [r;v];

% Set your integration time span: tspan
tF = N_orbit*T_period; % 
tspan = linspace(0, tF, N);

% parameters
parameters.rE = rE;              %[Km]
parameters.wE = wE;              % rad/s;
parameters.mu = mu;
parameters.drag.CD = CD;         %[]
parameters.drag.AM = AM;      %[m^2/kg]

parameters.drag.rE = rE;         %[Km]
parameters.j2 = astroConstants(9);
parameters.kep = kep0;

% Set ODE solver options (e.g., RelTol, AbsTol): options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14, 'Events', @terminate);

%% CAR METHOD
% integration in CAR form  
tic
[T, S] = ode113( @(t,s) eq_motion_CAR( t, s, @(t,s) acc_pert_fun_CAR(t,s,parameters), parameters ), tspan, s0, options);
time_int_car = toc;

fprintf('\nIntegration time of car integration %4.2f s \n', time_int_car)

% Conversion of cartesian state vector matrix in keplerian 
kep_matrix = zeros(length(S), 6);

for i = 1:length(S) 
    kep = car2kepRAD(S(i,1:3), S(i,4:6), mu);
    kep_matrix(i,:) = kep;
end

%% Study influence of acceleration during orbit
acc = zeros(length(S),2);

% Calculation of acceleration part
for i = 1:length(S)
    [acc_sum, acc_drag, acc_j2] = acc_pert_fun_CAR(0, S(i,:), parameters);
    acc(i,:) = [norm(acc_drag), norm(acc_j2)];
end

% plot of acceleration effect
figure
semilogy(tspan./T_period, acc, 'LineWidth', 1.5)
hold on
grid on
legend('Drag', 'J2')
xlabel('Time [T]')
ylabel('Acceleration [m/s^2]', 'Interpreter', 'latex')
title('Module of acceleration during time')

%% GAUSS METHOD

% Set Initial condition of keplerian elements 
s0 = kep0;

% Numerical integration of the equations of motion
tic
[T, S] = ode113( @(t,s) eq_motion_GAUSS( t, s, @(t,s) acc_pert_fun_RWS(t,s,parameters), parameters ), tspan, s0, options );
time_int_gauss = toc;

fprintf('\nIntegration time of car integration %4.2f s \n', time_int_gauss)

%% Effect on single disturbation

% ONLY J2
parameters.drag.CD = 0;
parameters.j2 = j2;

[T, S_j2] = ode113( @(t,s) eq_motion_GAUSS( t, s, @(t,s) acc_pert_fun_RWS(t,s,parameters), parameters ), tspan, s0, options );

% ONLY DRAG
parameters.drag.CD = CD;
parameters.j2 = 0;

[T, S_drag] = ode113( @(t,s) eq_motion_GAUSS( t, s, @(t,s) acc_pert_fun_RWS(t,s,parameters), parameters ), tspan, s0, options );

% wrap angle between 0 and 2pi
S(:,3) = wrapToPi(S(:,3)); % inclination
S(:,4) = wrapTo2Pi(S(:,4)); % OM
S(:,5) = wrapTo2Pi(S(:,5)); % om

% wrap angle between 0 and 2pi
S_j2(:,3) = wrapToPi(S_j2(:,3)); % inclination
S_j2(:,4) = wrapTo2Pi(S_j2(:,4)); % OM
S_j2(:,5) = wrapTo2Pi(S_j2(:,5)); % om

% wrap angle between 0 and 2pi
S_drag(:,3) = wrapToPi(S_drag(:,3)); % inclination
S_drag(:,4) = wrapTo2Pi(S_drag(:,4)); % OM
S_drag(:,5) = wrapTo2Pi(S_drag(:,5)); % om

% unwrap th for car integration
kep_matrix(:,end) = unwrap(kep_matrix(:,end));

%% figure and comparison

%adjust tspan 
tspan = tspan./T_period;
%%
% plot all keplerian orbit variation and filter
figure();
tiledlayout(2,3);
% a plot
hold on;
grid on;
plot( tspan, kep_matrix(:,1), '-b', 'LineWidth', 2)
plot( tspan, S(:,1), '-r', 'LineWidth', 2)
plot( tspan, S_j2(:,1), '-g', 'LineWidth', 2)
plot( tspan, S_drag(:,1), '-c', 'LineWidth', 2)
xlabel('Time [T]');
ylabel('a [km]');
title('a plot');

legend('Cartesian', 'Gauss', 'J2 Effect', 'Drag Effect')

% Find indexes 
index1 = find(tspan > 0.875, 1, 'first');
index2 = find(tspan > 0.9, 1, 'first');

height = abs(S(index1,1) - min(S(index1:index2,1)));
rectangle('Position', [tspan(index1), S(index1,1) - height, tspan(index2) - tspan(index1), height], ...
	'EdgeColor', 'r', 'LineWidth', 1);

ax2 = axes('position',[.65 .175 .25 .25]);
box on % put box around new pair of axes
hold on
plot(tspan(index1:index2),kep_matrix(index1:index2,1), '-b', 'LineWidth', 2) % plot on new axes
plot(tspan(index1:index2),S(index1:index2,1), '-r', 'LineWidth', 2) % plot on new axes
plot(tspan(index1:index2),S_j2(index1:index2,1), '-g', 'LineWidth', 2) % plot on new axes
axis tight
%

figure
% e plot
hold on;
grid on;
plot( tspan, kep_matrix(:,2), '-b', 'LineWidth', 2)
plot( tspan, S(:,2), '-r', 'LineWidth', 2)
plot( tspan, S_j2(:,2), '-g', 'LineWidth', 2)
plot( tspan, S_drag(:,2), '-c', 'LineWidth', 2)
xlabel('Time [T]');
ylabel('e [-]');
title('e plot');
grid on;

legend('Cartesian', 'Gauss', 'J2 Effect', 'Drag Effect')

ax2 = axes('position',[.65 .175 .25 .25]);
box on % put box around new pair of axes
hold on
plot(tspan, S_drag(:,2), '-c', 'LineWidth', 2) % plot on new axes
axis tight

% i plot
figure
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,3)), '-b', 'LineWidth', 2)
plot( tspan, rad2deg(S(:,3)), '-r', 'LineWidth', 2)
plot( tspan, rad2deg(S_j2(:,3)), '-g', 'LineWidth', 2)
plot( tspan, rad2deg(S_drag(:,3)), '-c', 'LineWidth', 2)
xlabel('Time [T]');
ylabel('i [deg]');
title('i plot');
grid on;

legend('Cartesian', 'Gauss', 'J2 Effect', 'Drag Effect')

% Find indexes 
index1 = find(tspan > 1.1, 1, 'first');
index2 = find(tspan > 1.13, 1, 'first');

height = abs(rad2deg(S(index1,3)) - rad2deg(S(index2,3)));
rectangle('Position', [tspan(index1), rad2deg(S(index1,3)), tspan(index2) - tspan(index1), 0.0003], ...
	'EdgeColor', 'r', 'LineWidth', 1);

ax2 = axes('position',[.65 .175 .25 .25]);
box on % put box around new pair of axes
hold on
plot(tspan(index1:index2),rad2deg(kep_matrix(index1:index2,3)), '-b', 'LineWidth', 2) % plot on new axes
plot(tspan(index1:index2),rad2deg(S(index1:index2,3)), '-r', 'LineWidth', 2) % plot on new axes
plot(tspan(index1:index2),rad2deg(S_j2(index1:index2,3)), '-g', 'LineWidth', 2) % plot on new axes
axis tight

% OM plot
figure
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,4)), '-b', 'LineWidth', 2)
plot( tspan, rad2deg(S(:,4)), '-r', 'LineWidth', 2)
plot( tspan, rad2deg(S_j2(:,4)), '-g', 'LineWidth', 2)
plot( tspan, rad2deg(S_drag(:,4)), '-c', 'LineWidth', 2)
xlabel('Time [T]');
ylabel('\Omega [deg]');
title('\Omega plot');
grid on;

legend('Cartesian', 'Gauss', 'J2 Effect', 'Drag Effect')

% om plot
figure
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,5)), '-b', 'LineWidth', 2)
plot( tspan, rad2deg(S(:,5)), '-r', 'LineWidth', 2)
plot( tspan, rad2deg(S_j2(:,5)), '-g', 'LineWidth', 2)
plot( tspan, rad2deg(S_drag(:,5)), '-c', 'LineWidth', 2)
xlabel('Time [T]');
ylabel('\omega [deg]');
title('\omega plot');
grid on;

legend('Cartesian', 'Gauss', 'J2 Effect', 'Drag Effect')

ax2 = axes('position',[.65 .175 .25 .25]);
box on % put box around new pair of axes
hold on
plot(tspan, rad2deg(S_drag(:,5)), '-c', 'LineWidth', 2) % plot on new axes
axis tight

% th plot
figure
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,6)), '-b', 'LineWidth', 2)
plot( tspan, rad2deg(S(:,6)), '-r', 'LineWidth', 2)
plot( tspan, rad2deg(S_j2(:,6)), '-g', 'LineWidth', 2)
plot( tspan, rad2deg(S_drag(:,6)), '-c', 'LineWidth', 2)
xlabel('Time [T]');
ylabel('\theta [deg]');
title('\theta plot');

legend('Cartesian', 'Gauss', 'J2 Effect', 'Drag Effect')

% Find indexes 
index1 = find(tspan > 0.75, 1, 'first');
index2 = find(tspan > 0.95, 1, 'first');

height = abs(rad2deg(S(index1,6)) - rad2deg(S(index2,6)));
rectangle('Position', [tspan(index1), rad2deg(S(index1,6)), tspan(index2) - tspan(index1), height], ...
	'EdgeColor', 'r', 'LineWidth', 1);

ax2 = axes('position',[.65 .175 .25 .25]);
box on % put box around new pair of axes
hold on
plot(tspan(index1:index2),rad2deg(kep_matrix(index1:index2,6)), '-b', 'LineWidth', 2) % plot on new axes
plot(tspan(index1:index2),rad2deg(S(index1:index2,6)), '-r', 'LineWidth', 2) % plot on new axes
plot(tspan(index1:index2),rad2deg(S_j2(index1:index2,6)), '-g', 'LineWidth', 2) % plot on new axes
plot(tspan(index1:index2),rad2deg(S_drag(index1:index2,6)), '-c', 'LineWidth', 2) % plot on new axes
axis tight