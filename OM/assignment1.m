clear all;
close all;
clc;
starter_pack

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
mu_Jupiter = astroConstants(15);

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
tspan_1 = linspace(t0_departure, tf_arrival, 100);
tspan_2 = linspace(t0_departure, tf_arrival, 100);
tspan_3 = linspace(t0_departure, tf_arrival, 100);

% al momento DV è una matrice 2d ma deve essere 3d alla fine
DV = zeros(length(tspan_1), length(tspan_2)); %length(tspan));

for i = 1:length(tspan_1)

    % Saturn position and velocity at chosen time
    flag1 = tspan_1(i);
    [kep_1,ksun] = uplanet(flag1, Saturn);
    [r_1, v_1] = kep2carRAD(kep_1, ksun);

    for j = 1:length(tspan_2)

        % Jupiter position and velocity at chosen time
        flag2 = tspan_2(j);
        if flag2<=flag1
            continue
        end

        [kep_2,ksun] = uplanet(flag2, Jupiter);
        [r_2, v_2] = kep2carRAD(kep_2, ksun);
        %first interplanetary arc
        delta_T1 = flag2 - flag1;
        delta_T1 = delta_T1*3600*24;
        [~, ~, ~, ~, VI,V_minus,TPAR,THETA] = lambertMR(r_1, r_2, delta_T1, ksun);
        delta_v_dep = norm(VI' - v_1);     % watchout for the dimensions IF ONE IS ROW VECTOR AND THE OTHER IS COLUMN VECTOR

        for k = 1:length(tspan_3)

            % Mars position and velocity at chosen time
            flag3 = tspan_3(k);
            if flag3<=flag2
                continue
            end
            [kep_3,ksun] = uplanet(flag3, Mars);
            [r_3, v_3] = kep2carRAD(kep_3, ksun);
            %second interplanetary arc
            delta_T2 = flag3 - flag2;
            delta_T2 = delta_T2*3600*24;
            [~, ~, ~, ~, V_plus,VF,TPAR,THETA] = lambertMR(r_2, r_3, delta_T2, ksun);
            delta_v_arr = norm(v_3 - VF');

            % ga assist
            d = norm(r_2);     
            V_p = v_2;      % velocity of the planet where the gravity assist happens
            v_asym_minus = V_minus' - V_p;
            v_asym_plus = V_plus' - V_p;
            % in order to determine rp we need to solve an implicit function
            fun = @(x) lab4_fun(x,  v_asym_minus, v_asym_plus, mu_Jupiter);
            r_p = fzero(fun, 10);
            % calculation of delta_vp
            eps_hyp_plus = norm(v_asym_plus)^2/2;
            vp_plus = sqrt(2*eps_hyp_plus + mu_Jupiter/r_p);
            eps_hyp_minus = norm(v_asym_minus)^2/2;
            vp_minus = sqrt(2*eps_hyp_minus + mu_Jupiter/r_p);
            delta_vp = abs(vp_plus-vp_minus);
            DV(i, j, k) = delta_v_dep +delta_v_arr; %delta_vp;
        end
    end
end

%%
figure
contour(tspan_1, tspan_2, DV(1, :, :))
figure
surf(tspan_1, tspan_2, DV(1, :, :), 'EdgeColor', 'none');
%clabel(C, h)
colorbar
xlabel('time of departure')
ylabel('time of arrival')
grid on

%% ciao
for i= 1:size(DV, 1)
    for j= 1:size(DV, 2)
        for k= 1:size(DV, 3)
            if DV(i, j, k) == 0
                DV(i, j, k) = 100;
            end
        end
    end
end


[dv_min, index] = min(DV, [], 'all');
[ii, jj, kk] = ind2sub(size(DV), index);
t_dep = tspan_1(ii);
t_ga = tspan_2(jj);
t_arr = tspan_3(kk);

dv12 = interplanetary_transfer_cost(t_dep, t_ga, Saturn, Jupiter, 1);
dv23 = interplanetary_transfer_cost(t_ga, t_arr, Jupiter, Mars, 2);
dvv = dv12 + dv23;

%%
% ga assist
[kep_1,ksun] = uplanet(t_dep, Saturn);
[r_1, v_1] = kep2carRAD(kep_1, ksun);
[kep_2,ksun] = uplanet(t_ga, Jupiter);
[r_2, v_2] = kep2carRAD(kep_2, ksun);
[kep_3,ksun] = uplanet(t_arr, Mars);
[r_3, v_3] = kep2carRAD(kep_3, ksun);

d = norm(r_2);     
V_p = v_2;      % velocity of the planet where the gravity assist happens

[~, ~, ~, ~, VI,V_minus,TPAR,THETA] = lambertMR(r_1, r_2, (t_ga-t_dep), ksun);
[~, ~, ~, ~, V_plus,VF,TPAR,THETA] = lambertMR(r_2, r_3, (t_arr-t_ga), ksun);
v_asym_minus = V_minus' - V_p;
v_asym_plus = V_plus' - V_p;
% in order to determine rp we need to solve an implicit function
fun = @(x) lab4_fun(x,  v_asym_minus, v_asym_plus, mu_Jupiter);
r_p = fzero(fun, 10);
% calculation of delta_vp
eps_hyp_plus = norm(v_asym_plus)^2/2;
vp_plus = sqrt(2*eps_hyp_plus + mu_Jupiter/r_p);
eps_hyp_minus = norm(v_asym_minus)^2/2;
vp_minus = sqrt(2*eps_hyp_minus + mu_Jupiter/r_p);
delta_vp = abs(vp_plus-vp_minus);

%%
plot_propagated(r_1, v_1, 0, ksun)

