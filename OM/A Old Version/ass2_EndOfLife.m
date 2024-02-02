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
th = 0.001; %[deg] % dende lab th = 297.23
kep0 = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th)];

%% DATA FOR BOTH PERTURBATION

%pert
mu = astroConstants(13);
rE = astroConstants(23);
j2 = astroConstants(9);
t_sid = 23*60*60 + 56*60 +4;
wE = 2*pi/t_sid;          %[rad/s] 
day_time = 24*60*60;

% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
N = 1e4;
N_orbit = 5000; 
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
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14, 'Events', @terminate);

%% CAR METHOD
% integration in CAR form with event detection (if alt - 100 == 0 ode will
% stop)
tic
[T, S, te] = ode113( @(t,s) eq_motion_CAR( t, s, @(t,s) acc_pert_fun_CAR(t,s,parameters), parameters ), tspan, s0, options); % te is the time where happen the event 
time_int_car = toc;

% find index in tspan vector for plot
index_event = find(tspan > te, 1);
%% print some information
fprintf('\nIntegration time of car integration %4.2f s \n', time_int_car)
fprintf('\nReach 100 Km of altitude after %4.2f day \n', te/(24*60*60))

% Conversion of cartesian state vector matrix in keplerian 
kep_matrix = zeros(length(S), 6);

for i = 1:length(S) 
    kep = car2kepRAD(S(i,1:3), S(i,4:6), mu);
    kep_matrix(i,:) = kep;
end

%% Secular effect
% Secular variation
a = kep0(1)*1000;
e = kep0(2);
J2 = parameters.j2;
fun_a = @(E) rho_find(a.*(1-e*cos(E))/1000 - rE).*(1+e*cos(E)).^(3/2)./sqrt(1-e*cos(E));
delta_a_sec = - CD*AM*a^2*integral(fun_a, 0, 2*pi, "ArrayValued",true)/1000;

figure

yyaxis left
plot_a = plot(tspan(1:index_event)./day_time, kep_matrix(:,1));
hold on
plot_sec_a_effect = plot([0:N_orbit]*T_period./day_time, kep0(1) + [0:N_orbit]*delta_a_sec, 'LineWidth', 2);
grid on
ylabel('a [km]')
xlabel('Time [day]')

yyaxis right
plot_e = plot(tspan(1:index_event)./day_time, kep_matrix(:,2));
ylabel('e [-]')

legend([plot_sec_a_effect, plot_a, plot_e], 'Secular a', 'True semi-major axis', 'True eccentricity')


%% Analyse and plot the results in CAR
% Change of Periapsis and Pericenter
rp_CAR =  kep_matrix(:,1).*((ones(length(S),1) - kep_matrix(:,2)))- rE*ones(length(S),1);
%rp_CAR = movmean(rp_CAR, N_filter1);

%
ra_CAR =  kep_matrix(:,1).*((ones(length(S),1) + kep_matrix(:,2)))- rE*ones(length(S),1);
%ra_CAR = movmean(ra_CAR, N_filter1);

%plot Change of Periapsis and apocenter (for drag analysis)
figure()
grid on
plot(tspan(1:index_event)./(24*60*60), rp_CAR)
hold on
plot(tspan(1:index_event)./(24*60*60), ra_CAR)
xlabel('Time [day]');
ylabel('Altitude [Km]');
title('Altitude of Pericenter and Apocenter')
legend('radius of Pericentre', 'radius of Apocenter')