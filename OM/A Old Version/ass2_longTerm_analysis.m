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
%CD = 3;
%AM = 10; % we see the effect of drag  

%other keplerian elements (arbitrary)
OM = 285; %[deg] -> più comodo per i calcoli  %dende lab OM = 284.67
om = 135; %[deg] % dende lab om = 135.52
th = 0.001; %[deg] % dende lab th = 297.23
kep0 = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th)];

%% DATA FOR BOTH PERTURBATION

%pert
mu = astroConstants(13);
rE = astroConstants(23);
j2 = astroConstants(9);
t_sid = 23*60*60 + 56*60 +4;
wE = 2*pi/t_sid;  %[rad/s] 

% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
N = 1e6;
N_orbit = 600; 
N_filter1 = N/10;
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
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);

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

%% Secular effect
% Secular variation of drag
a = kep0(1)*1000; %calculation are in m
e = kep0(2);
J2 = parameters.j2;
fun_a = @(E) rho_find(a.*(1-e*cos(E))/1000 - rE).*(1+e*cos(E)).^(3/2)./sqrt(1-e*cos(E));
delta_a_sec = - CD*AM*a^2*integral(fun_a, 0, 2*pi, "ArrayValued",true)/1000; %in Km
day_time = 24*60*60;

% secular variation of J2
n = sqrt(mu/a^3);
a = a/1000; %calculation are in Km

OM_dot_sec = -3/2*n*J2*rE^2/(a^2*(1-e^2)^2)*cos(deg2rad(i));
om_dot_sec = 3/2*n*J2*rE^2/(a^2*(1-e^2)^2)*(2 - 5/2*sin(deg2rad(i))^2);

delta_OM_sec = -3*pi*J2*rE^2/(a^2*(1-e^2)^2)*cos(deg2rad(i));
delta_om_sec = 3*pi*J2*rE^2/(a^2*(1-e^2)^2)*(2 - 5/2*sin(deg2rad(i))^2);

%% Study influence of acceleration during orbit
acc = zeros(length(S),2);

% Calculation of acceleration part
for i = 1:length(S)
    [acc_sum, acc_drag, acc_j2] = acc_pert_fun_CAR(0, S(i,:), parameters);
    acc(i,:) = [norm(acc_drag), norm(acc_j2)];
end

% plot of acceleration effect
figure
semilogy(tspan./T_period, acc)
hold on
grid on
semilogy(tspan./T_period, movmean(acc, N_filter1), 'LineWidth', 2)
legend('Drag', 'J2', 'Drag (filtered)', 'J2 (filtered)')
xlabel('time [T]')
ylabel('acceleration m/s^2')
title('Module of acceleration during time')


%% movie of change orbit
%{
% detect when the SC complete a full orbit
th_test = unwrap(kep_matrix(:,end)) -2*pi;

% Plot variation of orbit during time with colorbar
%{
%clean vector 
test = S(:,1:3);
orbit_num_plot = 125;
for i = 1:length(th_test)
    num_orbit = fix(th_test(i)./(2*pi));
    if( mod(num_orbit, orbit_num_plot) ~= 0 )
        test(i,:) = [0,0,0];
    end
end

figure
opt_earth.Units = 'Km';
%%Create two axes
planet3D('Earth', opt_earth);
hold on
grid on
scatter3(test(:,1), test(:,2), test(:,3), 5, fix(th_test./(2*pi)), 'filled');
%%Link them together
cb = colorbar();
clim([0,N_orbit]);
title(cb, 'Periods', 'Fontsize', 24)

xlabel('X [Km]', 'Fontsize', 24);
ylabel('Y [Km]', 'Fontsize', 24);
zlabel('Z [Km]', 'Fontsize', 24);
title('Orbit change');
%}
%}

%% GAUSS METHOD

% Set Initial condition of keplerian elements 
s0 = kep0;

% Set ODE solver options (e.g., RelTol, AbsTol): options
clear options;
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);

% Numerical integration of the equations of motion
tic
[T, S] = ode113( @(t,s) eq_motion_GAUSS( t, s, @(t,s) acc_pert_fun_RWS(t,s,parameters), parameters ), tspan, s0, options );
time_int_gauss = toc;

fprintf('\nIntegration time of car integration %4.2f s \n', time_int_gauss)

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

%adjust tspan 
tspan = tspan./T_period;

%% FILTRATION OF GAUSS METHOD
close all;
N_filter1 = N/10;
N_filter2 = N/3;
% filte
a_filter1 = movmean(S(:,1), N_filter1);
e_filter1 = movmean(S(:,2), N_filter1);
i_filter1 = rad2deg(wrapToPi(movmean(S(:,3), N_filter1)));
OM_filter1 = rad2deg(wrapTo2Pi(movmean(S(:,4), N_filter1)));
om_filter1 = rad2deg(wrapTo2Pi(movmean(S(:,5), N_filter1)));
th_filter1 = rad2deg(movmean(S(:,6), N_filter1));

% wrap angle between 0 and 2pi
S(:,3) = wrapToPi(S(:,3)); % inclination
S(:,4) = wrapTo2Pi(S(:,4)); % OM
S(:,5) = wrapTo2Pi(S(:,5)); % om

% unwrap th for car integration
kep_matrix(:,end) = unwrap(kep_matrix(:,end));

%% figure and comparison

% plot all keplerian orbit variation and filter

% a plot
figure
hold on;
grid on;
plot( tspan, S(:,1), '-r', 'LineWidth', 1.5)
plot( tspan, a_filter1, '-g', 'LineWidth', 2)
plot([0:N_orbit], kep0(1,1) + [0:N_orbit]*delta_a_sec, '-b', 'LineWidth', 2);
xlabel('time [T]');
ylabel('a [km]');
%title('a plot');

legend('Oscillation', 'Secular', 'Secular effect drag')

% e plot
figure
hold on;
grid on;
plot( tspan, S(:,2), '-r', 'LineWidth', 1.5)
plot( tspan, e_filter1, '-g', 'LineWidth', 2)
xlabel('time [T]');
ylabel('e [-]');
%title('e plot');
grid on;

legend('Oscillation', 'Secular')

% i plot
figure
hold on;
grid on;
plot( tspan, rad2deg(S(:,3)), '-r', 'LineWidth', 1.5)
plot( tspan, i_filter1, '-g', 'LineWidth', 2)
xlabel('time [T]');
ylabel('i [deg]');
%title('i plot');
grid on;

legend('Oscillation', 'Secular')

% OM plot
figure
hold on;
grid on;
plot( tspan, rad2deg(S(:,4)), '-r', 'LineWidth', 1.5)
plot( tspan, OM_filter1, '-g', 'LineWidth', 2)
plot([0:N_orbit], rad2deg(kep0(1,4) + [0:N_orbit]*delta_OM_sec), '-b', 'LineWidth', 2);
xlabel('Time [T]');
ylabel('\Omega [deg]');
%title('\Omega plot');
grid on;

legend('Oscillation', 'Secular', 'Secular effect J2')

% om plot
figure
hold on;
grid on;
plot( tspan,  rad2deg(S(:,5)), '-r', 'LineWidth', 2)
plot( tspan, om_filter1, '-g', 'LineWidth', 1.5)
plot([0:N_orbit], rad2deg(kep0(1,5) + [0:N_orbit]*delta_om_sec), '-b', 'LineWidth', 2);
xlabel('Time [day]');
ylabel('\omega [deg]');
%title('\omega plot');
grid on;

legend('Oscillation', 'Secular', 'Secular effect J2')
%%
% th plot
figure
hold on;
grid on;
plot( tspan, rad2deg(S(:,end)), '-r', 'LineWidth', 1)
plot( tspan, th_filter1, '-g' , 'LineWidth', 2)
xlabel('Time [T]');
ylabel('\theta [deg]');
%title('\theta plot');

legend('Oscillation', 'Secular')

%% comparison error
figure()
grid on;
err = abs((kep_matrix(:,1:end)-S(:,1:end)));

% normalize err
err(:, 1) = err(:,1)./kep0(1);
err(:, 2) = err(:,2)./kep0(2);
err(:, 3:end) = err(:,3:end)./(2*pi);
semilogy( tspan, err, '-', 'LineWidth', 1.5)
grid on
legend('a','e','i','\Omega','\omega', '\theta')
xlabel('Time [T]');
ylabel('err [-]');

%title('Error from GAUSS and CAR method')


