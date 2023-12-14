clear all;
close all;
clc;

%%%% OM PROJECT -> Assignments 1: : Interplanetary Explorer Mission %%%%%

% Data from File
% Group ID:2336 
% Departure: Saturn 
% Flyby: Jupiter 
% Arrival: Mars 
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058

% Constant DATA
N = 1000; % can be change  
Saturn = 6; % for ibody in uplanet
Jupiter = 5;
Mars = 4;
mu_S = astroConstants(4);

% time conversion
t0_departure = date2mjd2000([2028 01 01 0 0 0]);
tf_arrival = date2mjd2000([2058 01 01 0 0 0]);
dt_window = tf_arrival - t0_departure;
tspan = linspace(t0_departure, tf_arrival, N);

%% test figure and plot
% -> in questo modo poco pratico (non va bene
% -> Non capisco perchè vuole così tanti hold on

% figure set up
figure();
hold on;

% set common options for plot planet 
opts.Units = 'Km'; % Units - (char array) 'AU', 'ft', 'km', 'm', 'mi', or 'nmi'

% Plot Jupiter
kepJ = uplanet(t0_departure, Jupiter);
[rJ, vJ] = kep2carRAD(kepJ, mu_S);
opts.Position = rJ; % Km
planet3D('Jupiter', opts);
hold on;
plot_Jupiter = plot3(rJ(1), rJ(2), rJ(3), 'ok', 'MarkerSize', 10);
hold on;

% Plot Saturn
kepS = uplanet(t0_departure, Saturn);
[rS, vS] = kep2carRAD(kepS, mu_S);
opts.Position = rS;
planet3D('Saturn', opts);
hold on;
plot_Saturn = plot3(rS(1), rS(2), rS(3), 'ob', 'MarkerSize', 10);
hold on;

% Plot Mars
kepM = uplanet(t0_departure, Mars);
[rM, vM] = kep2carRAD(kepM, mu_S);
opts.Position = rM;
planet3D('Mars', opts);
hold on;
plot_Mars = plot3(rM(1), rM(2), rM(3), 'or', 'MarkerSize', 10);
hold on;
grid on;

% Plot Sun
opts.Position = [0 0 0];
plot_Sun = plot3(0, 0, 0, 'oy', 'MarkerSize', 10);
planet3D('Sun', opts);

% adjust figure layout
grid on
axis('auto xy')
zlim([-2 2].*1e8)

xlabel('X [Km]');
ylabel('Y [Km]');
zlabel('Z [Km]');

title('Solar system')

% Legend
legend([plot_Sun, plot_Mars, plot_Jupiter, plot_Saturn], 'Sun', 'Mars', 'Jupiter', 'Saturn')


%% first search of optimal DV


%GRID SEARCH OVER THE THREEE DEGREES OF FREEDOM
    %for each departure time
        %for each GA time
            %for each arrival time
%Compute and store DV(dep,ga,arr)


% deciding the TIMESPAN for departure(1), GA time(2) and arrival(3)

%al momento i tempi li ho scelti a caso
tspan_1 = linspace(t0_departure, t0_departure + 365*20, 1000);
tspan_2 = linspace(t0_departure+ 365, t0_departure + 365*30, 1000);
tspan_3 = linspace(tf_arrival-1, tf_arrival, 1);


% al momento DV è una matrice 2d ma deve essere 3d alla fine
DV = zeros(length(tspan_1), length(tspan_2)); %length(tspan));

for i = 1:length(tspan_1)
    flag1 = tspan_1(i);
    [kep_1,ksun] = uplanet(flag1, Saturn);
    [r_1, v_1] = kep2carRAD(kep_1, ksun);
    for j = 1:length(tspan_2)
        flag2 = tspan_2(j);
        [kep_2,ksun] = uplanet(flag2, Jupiter);
        [r_2, v_2] = kep2carRAD(kep_2, ksun);
        delta_T = flag2 - flag1;
        [~, ~, ~, ~, VI,VF,TPAR,THETA] = lambertMR(r_1, r_2, delta_T, ksun);
        delta_v1 = norm(VI' - v_1);     % watchout for the dimensions IF ONE IS ROW VECTOR AND THE OTHER IS COLUMN VECTOR 
        if delta_v1 > 1e3
            delta_v1 = 1e3;
        end
        DV(i, j) = delta_v1;
        for k = 1:length(tspan_3)
            flag3 = tspan_3(k);
            [kep_3,ksun] = uplanet(flag3, Mars);
            [r_3, v_3] = kep2carRAD(kep_3, ksun);
        end
    end
end

%%
figure
contour(tspan_1, tspan_2, DV)
figure
surf(tspan_1, tspan_2, DV, 'EdgeColor', 'none');
%clabel(C, h)
colorbar
xlabel('time of departure')
ylabel('time of arrival')


