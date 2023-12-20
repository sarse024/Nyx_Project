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
OM = 50; %[deg] -> più comodo per i calcoli
om = 60; %[deg]
th = 0;
kep0 = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th)];

%% Codice Aaron ground TRACK E ROBA varia

%pert
mu = astroConstants(13);
rE = astroConstants(23);
j2 = astroConstants(9);
t_sid = 23*60*60 + 56*60 +4;
wE = 2*pi/t_sid;          %[rad/s] 
[r,v] = kep2car(a,e,i/180*pi,OM,om,th,astroConstants(13));    %[m;m/s]
r_v_vect = [r;v];         
plot_orbit_pert(mu,rE,j2,CD,AM,r_v_vect,1000,wE);

%unperturbed orbit plot
mu = astroConstants(13);  %[km^3/s^2]
[r,v]=kep2car(a,e,i/180*pi,OM,om,th,mu);
r_v_vect=[r;v];         %[m;m/s]
odefun=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
options=odeset('RelTol',1e-13,'AbsTol',1e-14);
[~,Y]=ode89(odefun,[0,2*pi*sqrt(a^3/mu)],r_v_vect,options);
figure()
Terra3d;
hold on;
plot3(Y(:,1),Y(:,2),Y(:,3),'-',LineWidth=3);
xlabel('X[km]');
ylabel('Y[km]');
zlabel('Z[km]');
axis equal;
grid on;

%% GROUND TRACKS

thG0=0;                 %[rad]
t_vect=0:24*3600;       %[s]
figure();
[alpha,delta,long,lat]=groundTrack(r_v_vect,kep0,thG0,t_vect,mu,wE,k,m);

%% DATA FOR BOTH PERTURBATION

% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
N = 100000;
N_orbit = 10; 
N_filter = N/10;

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
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

%% CAR METHOD
% integration in CAR form  
tic
[T, S] = ode113( @(t,s) eq_motion_CAR( t, s, @(t,s) acc_pert_fun_CAR(t,s,parameters), parameters ), tspan, s0, options);
time_int_car = toc;

fprintf('\nIntegration time of car integration %4.2f s \n', time_int_car)

%% Analyse and plot the results in CAR
figure()
Terra3d
plot3( S(:,1), S(:,2), S(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

% Conversion of cartesian state vector matrix in keplerian 
kep_matrix = zeros(length(S), 6);

for i = 1:length(S) 
    kep = car2kepRAD(S(i,1:3), S(i,4:6), mu);
    kep_matrix(i,:) = kep;
end

% Change of Periapsis and Pericenter
rp_CAR =  kep_matrix(:,1).*((ones(length(S),1) - kep_matrix(:,2)))- rE*ones(length(S),1);
rp_CAR = movmean(rp_CAR, N_filter);

%
ra_CAR =  kep_matrix(:,1).*((ones(length(S),1) + kep_matrix(:,2)))- rE*ones(length(S),1);
ra_CAR = movmean(ra_CAR, N_filter);

% Study influence of acceleration during orbit
acc = zeros(length(S),2);

% Calculation of acceleration part
for i = 1:length(S)
    [acc_sum, acc_drag, acc_j2] = acc_pert_fun_CAR(0, S(i,:), parameters);
    acc(i,:) = [norm(acc_drag), norm(acc_j2)];
end

% plot of acceleration effect
figure
semilogy(tspan, acc)
hold on
grid on
semilogy(tspan, movmean(acc, N_filter), 'LineWidth', 3)
legend('Drag', 'J2', 'Drag (filtered)', 'J2 (filtered)')
title('Module of acceleration during time')

%% movie of change orbit
%{
% detect when the SC complete a full orbit
a = ischange(rad2deg(kep_matrix(:,6)), 'Linear', 'MaxNumChanges', N_orbit); %understan when there is a big variation between two consecutive point
index = find(a,N_orbit); % find all position of the variatios point

figure
Terra3d
grid on;
title('CHANGE of orbit in time');

j = index(1);
k = index(2) + 2;
orbit0 = plot3(S(j:k,1), S(j:k,2), S(j:k,3), '-r', 'LineWidth',2);
orbit = plot3(1,1,1);

legend([orbit0, orbit], 'Initial Orbit', 'Actual Orbit');

for i = 2:length(index)-1 % i is number of orbit
    delete(orbit);
    j = index(i);
    k = index(i+1) + 2; %plus two to draw full orbit (whitout there is a discontinuity)
    orbit = plot3(S(j:k,1), S(j:k,2), S(j:k,3), '-b', 'LineWidth',2);
    str = append('Orbit ',num2str(i)); %concatenete two string
    legend([orbit0, orbit], 'Initial Orbit', str);

    title(append('CHANGE of orbit in time t = ', num2str(tspan(j)/(24*3600), '%2.1f'), ' day'));
    pause(0.1)
end
% plot last orbit
delete(orbit);
j = index(N_orbit - 1);
orbit = plot3(S(j:end,1), S(j:end,2), S(j:end,3), '-g', 'LineWidth',3);
legend([orbit0, orbit], 'Initial Orbit', 'Final Orbit');
%}

%% GAUSS METHOD

% Set Initial condition of keplerian elements 
s0 = kep0;

% Numerical integration of the equations of motion
tic
[T, S] = ode113( @(t,s) eq_motion_GAUSS( t, s, @(t,s) acc_pert_fun_RWS(t,s,parameters), parameters ), tspan, s0, options );
time_int_gauss = toc;

fprintf('\nIntegration time of car integration %4.2f s \n', time_int_gauss)

%% FILTRATION OF GAUSS METHOD

% filter
a_filter = movmean(S(:,1), N_filter);
e_filter = movmean(S(:,2), N_filter);
i_filter = rad2deg(wrapToPi(movmean(S(:,3), N_filter)));
OM_filter = rad2deg(wrapTo2Pi(movmean(S(:,4), N_filter)));
om_filter = rad2deg(wrapTo2Pi(movmean(S(:,5), N_filter)));
th_filter = rad2deg(movmean(S(:,6), N_filter));

% wrap angle between 0 and 2pi
S(:,3) = wrapToPi(S(:,3)); % inclination
S(:,4) = wrapTo2Pi(S(:,4)); % OM
S(:,5) = wrapTo2Pi(S(:,5)); % om

% unwrap th for car integration
kep_matrix(:,end) = unwrap(kep_matrix(:,end)) -2*pi;

%% figure and comparison

%adjust tspan 
tspan = tspan./T_period;

% plot all keplerian orbit variation and filter
figure();
tiledlayout(2,3);
title('Keplerian change by perturbation')

% a plot
nexttile
hold on;
grid on;
plot( tspan, kep_matrix(:,1), '-b' )
plot( tspan, S(:,1), '-r' )
plot( tspan, a_filter, '-y' )
xlabel('time [T]');
ylabel('a [km]');
title('a plot');

legend('Cartesian', 'Gauss', 'Filtered')

% e plot
nexttile
hold on;
grid on;
plot( tspan, kep_matrix(:,2), '-b' )
plot( tspan, S(:,2), '-r' )
plot( tspan, e_filter, '-y' )
xlabel('time [T]');
ylabel('e [-]');
title('e plot');
grid on;

legend('Cartesian', 'Gauss', 'Filtered')

% i plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,3)), '-' )
plot( tspan, rad2deg(S(:,3)), '-r' )
plot( tspan, i_filter, '-y' )
xlabel('time [T]');
ylabel('i [deg]');
title('i plot');
grid on;

legend('Cartesian', 'Gauss', 'Filtered')

% OM plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,4)), '-b' )
plot( tspan, rad2deg(S(:,4)), '-r' )
plot( tspan, OM_filter, '-y' )
xlabel('time [T]');
ylabel('OM [deg]');
title('OM plot');
grid on;

legend('Cartesian', 'Gauss', 'Filtered')

% om plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,5)), '-b' )
plot( tspan, rad2deg(S(:,5)), '-r' )
plot( tspan, om_filter, '-y' )
xlabel('time [T]');
ylabel('om [deg]');
title('om plot');
grid on;

legend('Cartesian', 'Gauss', 'Filtered')

% th plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,6)), '-' )
plot( tspan, rad2deg(S(:,6)), '-r' )
plot( tspan, th_filter, '-y' )
xlabel('time [T]');
ylabel('f [deg]');
title('f plot');

legend('Cartesian', 'Gauss', 'Filtered')

%% comparison error
figure()
grid on;
err = abs((kep_matrix(:,1:end)-S(:,1:end)));
semilogy( tspan, err, '-' )
grid on
legend('a','e','i','OM','om', 'th')
xlabel('time [T]');
ylabel('err [-]');
title('Error from GAUSS and CAR method')


%% plot Change of Periapsis and apocenter (for drag analysis)
figure()
grid on
plot(tspan, rp_CAR)
hold on
plot(tspan, ra_CAR)
xlabel('time [T]');
ylabel('Altitude [Km]');
title('Altitude of Pericenter and Apocenter (filtered)')
legend('radius of Pericentre', 'radius of Apocenter')


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
plot(tspan, a_filter, '-y')
xlabel('time [T]');
ylabel('a [km]');
title('a plot');
legend('Cartesian', 'Gauss', 'Secular (filtered)')

nexttile
hold on;
grid on;
plot( tspan(1:index), kep_matrix(1:index,1), '-b' )
plot( tspan(1:index), S(1:index,1), '-r' )
xlabel('time [T]');
ylabel('a [km]');
title('10 orbit study');
legend('Cartesian', 'Gauss')


% cose da fare:
%{
- Creare il movie va inserito exportgif 
- Provare ad implementare car usando quella di aaron
- Capire altri filtri e aggiungerli -> da vedere con vale
- Cercare quache dato reale
- Studio sul lungo periodo
%}

%%

%{
% verifica accelerazioni sistema di riferimento
kep0
a_CAR = acc_pert_fun_CAR(5, s0, parameters)

a_RWS = acc_pert_fun_RWS(5, kep0, parameters)

om = kep0(6) + kep0(5);

OM = kep0(4);

i = kep0(3);

Rot_mat=[cos(om)*cos(OM)-sin(om)*sin(OM)*cos(i), ...  
    -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM), sin(i)*sin(OM); sin(OM)*cos(om)+...
    sin(om)*cos(OM)*cos(i),-sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM), -sin(i)*cos(OM);...
    sin(om)*sin(i), cos(om)*sin(i), cos(i)];

test = Rot_mat'*a_CAR

err = a_RWS - test
%}



