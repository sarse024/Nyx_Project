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
OM = 284.67; %[deg] -> più comodo per i calcoli  %dende lab OM = 284.67
om = 135.52; %[deg] % dende lab om = 135.52
th = 297.23; %[deg] % dende lab th = 297.23
kep0 = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th)];

%% DATA FOR BOTH PERTURBATION

%pert
mu = astroConstants(13);
rE = astroConstants(23);
j2 = astroConstants(9);
t_sid = 23*60*60 + 56*60 +4;
wE = 2*pi/t_sid;          %[rad/s] 

% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
N = 1000;
N_orbit = 2; 
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
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

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
% Secular variation
a = kep0(1)*1000;
e = kep0(2);
J2 = parameters.j2;
fun_a = @(E) rho_find(a.*(1-e*cos(E))/1000 - rE).*(1+e*cos(E)).^(3/2)./sqrt(1-e*cos(E));
delta_a_sec = - CD*AM*a^2*integral(fun_a, 0, 2*pi, "ArrayValued",true)/1000;
%a = a - a_sec;

figure

yyaxis left
hold on
plot([0:N_orbit], kep0(1) + [0:N_orbit]*delta_a_sec, 'LineWidth', 2);
hold on;    
plot(tspan./T_period, kep_matrix(:,1))
grid on
ylabel('a [km]')

yyaxis right
plot(tspan./T_period, kep_matrix(:,2))
ylabel('e [-]')

legend('De Orbiting data')

n = sqrt(mu/a^3);
a = a/1000;

OM_dot_sec = -3/2*n*J2*rE^2/(a^2*(1-e^2)^2)*cos(deg2rad(i));
om_dot_sec = 3/2*n*J2*rE^2/(a^2*(1-e^2)^2)*(2 - 5/2*sin(deg2rad(i))^2);

delta_OM_sec = -3*pi*J2*rE^2/(a^2*(1-e^2)^2)*cos(deg2rad(i));
delta_om_sec = 3*pi*J2*rE^2/(a^2*(1-e^2)^2)*(2 - 5/2*sin(deg2rad(i))^2);





%% Analyse and plot the results in CAR
figure()
Terra3d
plot3( S(:,1), S(:,2), S(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;



% Change of Periapsis and Pericenter
rp_CAR =  kep_matrix(:,1).*((ones(length(S),1) - kep_matrix(:,2)))- rE*ones(length(S),1);
rp_CAR = movmean(rp_CAR, N_filter1);

%
ra_CAR =  kep_matrix(:,1).*((ones(length(S),1) + kep_matrix(:,2)))- rE*ones(length(S),1);
ra_CAR = movmean(ra_CAR, N_filter1);

%plot Change of Periapsis and apocenter (for drag analysis)
figure()
grid on
plot(tspan./(24*60*60), rp_CAR)
hold on
plot(tspan./(24*60*60), ra_CAR)
xlabel('time [Day]');
ylabel('Altitude [Km]');
title('Altitude of Pericenter and Apocenter (filtered)')
legend('radius of Pericentre', 'radius of Apocenter')

%% Study influence of acceleration during orbit
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
semilogy(tspan, movmean(acc, N_filter1), 'LineWidth', 3)
legend('Drag', 'J2', 'Drag (filtered)', 'J2 (filtered)')
title('Module of acceleration during time')




%% movie of change orbit
%
% detect when the SC complete a full orbit
th_test = unwrap(kep_matrix(:,end)) -2*pi;

%{
% count how much point there are in that orbit
index = histcounts(fix(th_test./(2*pi)), 1:N_orbit);

%figure plot
figure('WindowState', 'maximized')
opt_earth.Units = 'Km';
planet3D('Earth', opt_earth);
hold on
grid on

title('Movie of changing orbit');

k = index(1) + 2; 
pointer = k - 2; % save the position of index of last point of a orbit
orbit0 = plot3(S(1:k,1), S(1:k,2), S(1:k,3), '-r', 'LineWidth',2);
orbit = plot3(1,1,1);
ra_plot = plot3(1,1,1);
rp_plot = plot3(1,1,1);

legend([orbit0, orbit], 'Initial Orbit', 'Actual Orbit', 'Location', 'southeast');

for i = 2:length(index)-1 % i is number of orbit
    delete(orbit);
    delete(ra_plot);
    delete(rp_plot);
    k = pointer + index(i+1) + 2; %plus two to draw full orbit (whitout there is a discontinuity)
    orbit = plot3(S(pointer:k,1), S(pointer:k,2), S(pointer:k,3), '-b', 'LineWidth',2);

    r = vecnorm(S(pointer:k,1:3),2,2);
    [ra, index_ra] = max(r);
    [rp, index_rp] = min(r);

    index_ra = index_ra + pointer;
    index_rp = index_rp + pointer;

    ra_plot = plot3(S(index_ra,1), S(index_ra,2), S(index_ra,3), 'om', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
    rp_plot = plot3(S(index_rp,1), S(index_rp,2), S(index_rp,3), 'oc', 'MarkerSize', 8, 'MarkerFaceColor', 'c');

    str1 = append('Orbit ',num2str(i)); %concatenete two string
    str2 = append('Apocenter h = ',num2str(ra - rE, '%3.f'), ' Km'); %concatenete two string
    str3 = append('Pericenter h = ',num2str(rp - rE, '%3.f'), ' Km'); %concatenete two string
    
    legend([orbit0, orbit, ra_plot, rp_plot], 'Initial Orbit', str1, str2, str3, 'Location', 'southeast');

    title(append('CHANGE of orbit in time t = ', num2str(tspan(pointer)/(24*3600), '%2.1f'), ' day'));
    pointer = k - 2;
    pause(0.1)
end
% plot last orbit
delete(orbit);
delete(ra_plot);
delete(rp_plot);
orbit = plot3(S(pointer:end,1), S(pointer:end,2), S(pointer:end,3), '-b', 'LineWidth',3);

r = vecnorm(S(pointer:end,1:3),2,2);
[ra, index_ra] = max(r);
[rp, index_rp] = min(r);

index_ra = index_ra + pointer;
index_rp = index_rp + pointer;

ra_plot = plot3(S(index_ra,1), S(index_ra,2), S(index_ra,3), 'om', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
rp_plot = plot3(S(index_rp,1), S(index_rp,2), S(index_rp,3), 'oc', 'MarkerSize', 8, 'MarkerFaceColor', 'c');

str1 = append('Orbit ',num2str(i)); %concatenete two string
str2 = append('Apocenter h = ',num2str(ra - rE, '%3.f'), ' Km'); %concatenete two string
str3 = append('Pericenter h = ',num2str(rp - rE, '%3.f'), ' Km'); %concatenete two string
    
legend([orbit0, orbit, ra_plot, rp_plot], 'Initial Orbit', str1, str2, str3, 'Location', 'southeast');

title(append('CHANGE of orbit in time t = ', num2str(tspan(end)/(24*3600), '%2.1f'), ' day'));
%}

%%
%clean vector 
test = S(:,1:3);
orbit_num_plot = 25;
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
title(cb, 'Periods')

xlabel('X [Km]');
ylabel('Y [Km]');
zlabel('Z [Km]');
title('Orbit change');

%% GAUSS METHOD

% Set Initial condition of keplerian elements 
s0 = kep0;

% Numerical integration of the equations of motion
tic
[T, S] = ode45( @(t,s) eq_motion_GAUSS( t, s, @(t,s) acc_pert_fun_RWS(t,s,parameters), parameters ), tspan, s0, options );
time_int_gauss = toc;

fprintf('\nIntegration time of car integration %4.2f s \n', time_int_gauss)

%% FILTRATION OF GAUSS METHOD

% filter
a_filter1 = movmean(S(:,1), N_filter1);
e_filter1 = movmean(S(:,2), N_filter1);
i_filter1 = rad2deg(wrapToPi(movmean(S(:,3), N_filter1)));
OM_filter1 = rad2deg(wrapTo2Pi(movmean(S(:,4), N_filter1)));
om_filter1 = rad2deg(wrapTo2Pi(movmean(S(:,5), N_filter1)));
th_filter1 = rad2deg(movmean(S(:,6), N_filter1));

a_filter2 = movmean(S(:,1), N_filter2);
e_filter2 = movmean(S(:,2), N_filter2);
i_filter2 = rad2deg(wrapToPi(movmean(S(:,3), N_filter2)));
OM_filter2 = rad2deg(wrapTo2Pi(movmean(S(:,4), N_filter2)));
om_filter2 = rad2deg(wrapTo2Pi(movmean(S(:,5), N_filter2)));
th_filter2 = rad2deg(movmean(S(:,6), N_filter2));

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
plot( tspan, a_filter1, '-y' )
plot( tspan, a_filter2, '-c' )
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
plot( tspan, e_filter1, '-y' )
plot( tspan, e_filter2, '-c' )
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
plot( tspan, i_filter1, '-y' )
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
plot( tspan, OM_filter1, '-y' )
xlabel('time [T]');
ylabel('OM [deg]');
title('OM plot');
grid on;

legend('Cartesian', 'Gauss', 'Filtered')

% om plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,5)), '-b', 'LineWidth', 1.5)
plot( tspan, rad2deg(S(:,5)), '-r', 'LineWidth', 1.5)
plot( tspan, om_filter1, '-y', 'LineWidth', 2)
plot( tspan, om_filter2, '-c', 'LineWidth', 2)
xlabel('time [T]');
ylabel('om [deg]');
title('om plot');
grid on;

legend('Cartesian', 'Gauss', 'Filtered 1', 'Filtered2')

% th plot
nexttile
hold on;
grid on;
plot( tspan, rad2deg(kep_matrix(:,6)), '-' )
plot( tspan, rad2deg(S(:,6)), '-r' )
plot( tspan, th_filter1, '-y' )
xlabel('time [T]');
ylabel('f [deg]');
title('f plot');

legend('Cartesian', 'Gauss', 'Filtered')

%% comparison error
figure()
grid on;
err = abs((kep_matrix(:,1:end)-S(:,1:end)));

% normalize err
err(:, 1) = err(:,1)./kep0(1);
err(:, 2) = err(:,2)./kep0(2);
err(:, 3:end) = err(:,3:end)./(2*pi);
semilogy( tspan, err, '-' )
grid on
legend('a','e','i','OM','om', 'th')
xlabel('time [T]');
ylabel('err [-]');
title('Error from GAUSS and CAR method')





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



Possible True data: 
39369 	DANDE DEB (LAB)	2013-055AC	DEBRIS	US	2013-09-29	AFWTR		94.25	inclin = 80.95 	675	290	MEDIUM
TLE: 1 39369U 13055AC  23353.47960201  .00102165  13047-5  10990-2 0  9995
     2 39369  80.9520 250.3014 0280319  37.9799 324.0837 15.27977448536733

DANDE DEB (LAB)                 SATELLITE PROGETTO
          6859          a                6846
     0.0273             e                0.0298
         80.95          i                80.2068      

%}
