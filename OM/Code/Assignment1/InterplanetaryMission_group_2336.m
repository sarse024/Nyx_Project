%% OM PROJECT -> Assignments 1: : Interplanetary Explorer Mission 

% Data from File
% Group ID:2336 
% Departure: Saturn 
% Flyby: Jupiter 
% Arrival: Mars 
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058

% constant DATA
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

% SOI of Jupiter during flyby
r_p_Jupiter = 778.479*1e6;
r_soi = r_p_Jupiter*(mu_Jupiter/mu_S)^(2/5);

%% PRELIMINAR STUDIES
% TRANSFER COST FROM SATURN TO JUPITER OVER THE TIME WINDOW

% time span of departure
tspan_1 = linspace(t0_departure, tf_arrival, 200);

% time span of the time of flight for the transfer
tspan_2 = linspace(0, tf_arrival-t0_departure, 200);

% preallocation of DV and single impulses
DV12 = ones(length(tspan_1), length(tspan_2))*NaN;
delta_T1_mat = ones(length(tspan_1), length(tspan_2))*NaN;
tspan_dep = NaT(length(tspan_1), length(tspan_2));
tspan_arr = NaT(length(tspan_1), length(tspan_2));

for i = 1:length(tspan_1)

    % Saturn position and velocity at chosen time
    flag1 = tspan_1(i);
    [kep_1,ksun] = uplanet(flag1, Saturn);
    [r_1, v_1] = kep2carRAD(kep_1, ksun);

    for j = 1:length(tspan_2)

        % Jupiter position and velocity at chosen time
        flag2 = flag1 + tspan_2(j);
        if flag2>tf_arrival
            continue
        end
        [kep_2,ksun] = uplanet(flag2, Jupiter);
        [r_2, ~] = kep2carRAD(kep_2, ksun);

        % creating tspan matrix
        tspan_dep(i, j) = datetime(mjd20002date(flag1));
        tspan_arr(i, j) = datetime(mjd20002date(flag2));
        
        % first interplanetary arc
        delta_T1_mat(i, j) = tspan_2(j);
        [~, ~, ~, ~, VI,~,~,~] = lambertMR(r_1, r_2, tspan_2(j)*3600*24, ksun);
        DV12(i, j) = norm(VI' - v_1);  
        if DV12(i, j)> 11.95
            DV12(i,j) = NaN;
        end

    end
end

% DV of first Lambert's arc
figure
contour(datenum(tspan_dep), datenum(tspan_arr), DV12, 0:0.5:15)
xticks(datenum(tspan_dep(1)):1465:datenum(tspan_dep(end, 1)))
ax.Ytick = 8;
datetick('x', 'mmmyyyy', 'keeplimits', 'keepticks')
datetick('y', 'mmmyyyy', 'keeplimits', 'keepticks')
colorbar
clim([min(DV12(:)), max(DV12(:))])
grid on
grid minor
hold on
contour(datenum(tspan_dep), datenum(tspan_arr), delta_T1_mat, 1000:1000:6000, 'k', 'ShowText','on')
grid on
hold off
xlabel('Saturn departure date')
ylabel('Jupiter arrival date')
title('Porkchop plot for Saturn-Jupiter transfer')

% It can be graphically observed that the main area of interest which
% shows a low cost has a departure time from Saturn that goes from 2028 
% to 2040, so the time window for the study can be further decreased to
% ignore the other departure dates.
% Also the maximum deltaT it takes is aroud 6000 days or 16 years and the
% minimum is 2000 days or 5 years


%% TRANSFER COST FROM JUPITER TO MARS OVER THE TIME WINDOW

% time span of departure
tspan_1 = linspace(t0_departure, tf_arrival, 300);
% time span of the time of flight for the transfer
tspan_2 = linspace(0, tf_arrival-t0_departure, 300);

% preallocation of DV and single impulses
DV12 = ones(length(tspan_1), length(tspan_2))*NaN;
delta_T1_mat = ones(length(tspan_1), length(tspan_2))*NaN;
tspan_dep = NaT(length(tspan_1), length(tspan_2));
tspan_arr = NaT(length(tspan_1), length(tspan_2));

for i = 1:length(tspan_1)

    % Saturn position and velocity at chosen time
    flag1 = tspan_1(i);
    [kep_1,ksun] = uplanet(flag1, Jupiter);
    [r_1, ~] = kep2carRAD(kep_1, ksun);

    for j = 1:length(tspan_2)

        % Jupiter position and velocity at chosen time
        flag2 = flag1 + tspan_2(j);
        if flag2>tf_arrival
            continue
        end
        [kep_2,ksun] = uplanet(flag2, Mars);
        [r_2, v_2] = kep2carRAD(kep_2, ksun);
        
        % creating tspan matrix
        tspan_dep(i, j) = datetime(mjd20002date(flag1));
        tspan_arr(i, j) = datetime(mjd20002date(flag2));
        
        % first interplanetary arc
        delta_T1_mat(i, j) = tspan_2(j);
        [~, ~, ~, ~, ~,V_minus,~,~] = lambertMR(r_1, r_2, tspan_2(j)*3600*24, ksun);
        DV12(i, j) = norm(v_2- V_minus'); 
        if DV12(i, j)> 11.95
            DV12(i,j) = NaN;
        end

    end
end

% DV of first Lambert's arc
figure
contour(datenum(tspan_dep), datenum(tspan_arr), DV12, 5:1:20)
xticks(datenum(tspan_dep(1)):1465:datenum(tspan_dep(end, 1)))
yticks(datenum(tspan_dep(1)):1465:datenum(tspan_dep(end, 1)))
datetick('x', 'mmmyyyy', 'keeplimits', 'keepticks')
datetick('y', 'mmmyyyy', 'keeplimits', 'keepticks')
colorbar
clim([min(DV12(:)), max(DV12(:))])
grid on
grid minor
hold on
contour(datenum(tspan_dep), datenum(tspan_arr), delta_T1_mat, 500:500:3000, 'k', 'ShowText','on')
hold off
xlabel('Jupiter departure date')
ylabel('Mars arrival date')
title('Porkchop plot for Jupiter-Mars transfer')

% in this case we can observe that there is not a favorable time to start
% the second Lamber arc but it's interesting to see that the minimum are
% spaced uniformly with a periodicity that coincides with the mutual
% synodic period of Jupiter and Mars and the maximum time it must take to
% have an acceptable cost is 3000 days, around 8.2 years and minimum 500
% days or 1.4 years

%% First search of optimal DV

% GRID SEARCH OVER THE THREEE DEGREES OF FREEDOM
    % for each departure time
        % for each TOF of first Lambert's arc
            %f or each for each TOF of second Lambert's arc
% Compute and store DV(dep,ga,arr)
% TIMESPANS are chosen thanks to the preliminary studies they can be cut down

% one step every month (12 years = 144 months), for tspan_2 we consider 
% 16-5 years --> 132 months (5 years and 16 years are the minimum and maximum
% time it can take for the arc to be acceptable in cost), for tspan_3 
% maximum 3000 days and minimum 500 days with one step every 10 days
tspan_1 = linspace(t0_departure, t0_departure+ 12*365.25, 144); 
tspan_2 = linspace(5*365.25, 16*365.25, 132);        
tspan_3 = linspace(500, 3000, 250);                 

% calculation of optimal Dv and its corresponding times
[DV1, dv_min1, index_val1, t_dep_1, deltat1_1, deltat2_1] = grid_search(tspan_1, tspan_2, tspan_3);
t_ga_1 = t_dep_1 +deltat1_1;
t_arr_1 = t_dep_1 +deltat1_1 + deltat2_1;

%% Second search of optimal Dv

% step every week in an interval of + 6 months, except for Deltat2 which is
% every day (for Mars position, t_span=1 week is too much)
tspan_1 = linspace(t_dep_1-364/2, t_dep_1+364/2, 53);  

% using deltaT instead of specific time
tspan_2 = linspace(deltat1_1-364/2, deltat1_1+364/2, 53);                                                                        
tspan_3 = linspace(deltat2_1-364/2, deltat2_1+364/2, 365);

% calculation of optimal Dv and its corresponding times
[DV2, dv_min2, ~, t_dep_2, deltat1_2, deltat2_2] = grid_search(tspan_1, tspan_2, tspan_3);
t_ga_2 = t_dep_2 +deltat1_2;
t_arr_2 = t_dep_2 +deltat1_2 + deltat2_2;

%% Third search of optimal Dv

% step every day in an interval of +- 1 month (two months span)
tspan_1 = linspace(t_dep_2-31, t_dep_2+31, 63); 

% using deltaT instead of specific time
tspan_2 = linspace(deltat1_2-31, deltat1_2+31, 63);                                                                        
tspan_3 = linspace(deltat2_2-31, deltat2_2+31, 63);    

% calculation of optimal Dv and its corresponding times
[DV3, dv_min3, index_val2, t_dep_3, deltat1_3, deltat2_3] = grid_search(tspan_1, tspan_2, tspan_3);
t_ga_3 = t_dep_3 +deltat1_3;
t_arr_3 = t_dep_3 +deltat1_3 + deltat2_3;

%% Plot of the two interplanetry legs

% The time of departure, fly-by and arrival are the ones calculated with
% the last research grid refinement
t_dep = t_dep_3;
t_ga = t_ga_3;
t_arr = t_arr_3;

% Position of the planets in best case scenario
[kep_1,ksun] = uplanet(t_dep, Saturn);
[r_1, ~] = kep2carRAD(kep_1, ksun);
[kep_2,ksun] = uplanet(t_ga, Jupiter);
[r_2, v_2] = kep2carRAD(kep_2, ksun);
[kep_3,ksun] = uplanet(t_arr, Mars);
[r_3, v_3] = kep2carRAD(kep_3, ksun);

% Plot of the various planets and the Sun
opts.Units = 'm'; 
opts.Position = [r_1(1), r_1(2), r_1(3)];
planet3D('Saturn', opts);  
hold on;
grid on;
opts.Units = 'm';
opts.Position = [r_2(1), r_2(2), r_2(3)];
planet3D('Jupiter', opts);  
opts.Units = 'mars'; 
opts.Position = [r_3(1), r_3(2), r_3(3)];
planet3D('Mars', opts);  
opts.Units = 'sun';
opts.Position = [0 0 0]; 
planet3D('Sun', opts); 
% Plot of the Lambert's arcs
dv12 = interplanetary_transfer_cost(t_dep, t_ga, Saturn, Jupiter, 1);
dv23 = interplanetary_transfer_cost(t_ga, t_arr, Jupiter, Mars, 2); 
title('Road to Mars')
legend('', '', '', '', 'Orbit of Saturn', '' , 'First Lambert arc', 'Orbit of Jupiter', 'Orbit of Mars', 'Second Lambert arc')
hold off
%% Gravity assist plot

% velocity and position of flyby planet(Jupiter)
d = norm(r_2);     
V_p = v_2;      

% determination of V_minus, V_plus, v_asym_minus and v_asym_plus necessary 
% to determine the radius of pericentre of the flyby
[~, ~, ~, ~, VI,V_minus,~,~] = lambertMR(r_1, r_2, (t_ga-t_dep)*3600*24, ksun);
[~, ~, ~, ~, V_plus,VF,~,~] = lambertMR(r_2, r_3, (t_arr-t_ga)*3600*24, ksun);
v_asym_minus = V_minus' - V_p;
v_asym_plus = V_plus' - V_p;
delta_vp = norm(v_asym_plus-v_asym_minus);

% in order to determine rp we need to solve an implicit function
delta=acos(dot(v_asym_minus,v_asym_plus)/(norm(v_asym_minus)*norm(v_asym_plus)));
fun_delta=@(r_p) delta-(asin(1/(1+(r_p*norm(v_asym_plus)^2)/mu_Jupiter))+asin(1/(1+(r_p*norm(v_asym_minus)^2)/mu_Jupiter)));
options = optimset('TolFun',1e-15,'Display','off');
r_per = fsolve(fun_delta, 71000, options);

% calculation of delta_vp at pericentre
eps_hyp_plus = norm(v_asym_plus)^2/2;
vp_plus = sqrt(2*(eps_hyp_plus + mu_Jupiter/r_per));
eps_hyp_minus = norm(v_asym_minus)^2/2;
vp_minus = sqrt(2*(eps_hyp_minus + mu_Jupiter/r_per));
delta_vp_powered = norm(vp_plus-vp_minus);
ratio_dv = delta_vp_powered/delta_vp;

% definition of u, versor perpendicular to the fly-by plane
u = cross(v_asym_minus, v_asym_plus)/ norm(cross(v_asym_minus, v_asym_plus));

% incoming hyperbolic arc parameters
e_minus =  1 + norm(r_per)*(norm(v_asym_minus)^2)/mu_Jupiter;
delta_minus = 2*asin(1/e_minus);
a_minus = r_per/(1-e_minus);

% outcoming hyperboli arc parameters
e_plus =  1 + norm(r_per)*(norm(v_asym_plus)^2)/mu_Jupiter;
delta_plus = 2*asin(1/e_plus);
a_plus = r_per/(1-e_plus);

% direction of r_pericentre, thanks to rodrigues rotation of v_asym_plus around u
d1 = -rodrigues_rot(v_asym_plus/norm(v_asym_plus), (pi/2-(delta_plus/2)), u);
d1 = d1/norm(d1);
r_p_vect = r_per*d1;

% direction of v_p_plus and v_p_minus
d2 = -cross(u, d1);  
vp_minus_v = vp_minus*(-d2);
vp_plus_v = vp_plus*d2;
r0 = [r_p_vect, r_p_vect];
v0 = [vp_minus_v, vp_plus_v];

% plot the planet
opts.Units = 'Km'; % Units - (char array) 'AU', 'ft', 'km', 'm', 'mi', or 'nmi'
figure
planet3D('Jupiter', opts);  
hold on;

% plot the hyperbolic arcs
plot_propagated(r0, v0, [600000 600000], mu_Jupiter)

% plot of SOI influence, too large to plot
%plot_propagated([r_soi; 0; 0], [0; sqrt(mu_Jupiter/r_soi); 0], 0, mu_Jupiter, 1)

% plot the asymptotes
c1 = (-a_minus+r_per)*d1;
c2 = (-a_plus+r_per)*d1;
quiver3(0, 0, 0, d1(1)*3500000, d1(2)*3500000, d1(3)*3500000)
quiver3(0,0,0, V_p(1)*100000, V_p(2)*100000, V_p(3)*100000)
legend('', 'Arriving hyperbola', 'Leaving hyperbola', 'Direction of pericentre', 'Planet velocity')
title('Flyby on Jupiter')
axis equal;
grid on;
hold off

%% Time in Jupiter SOI estimation

% time from entrance to pericentre
n_minus = sqrt(mu_Jupiter/(-a_minus)^3);
h_minus = cross(r_p_vect, vp_minus_v);
p_minus = norm(h_minus)^2/mu_Jupiter;
theta = acos((p_minus/r_soi-1)/e_minus);
calc = sqrt(-(1+e_minus)/(1-e_minus));
F = 2*atanh(tan(theta/2)/calc);
t_minus = (e_minus*sinh(F)-F)/n_minus;

% time from pericentre to exit
n_plus = sqrt(mu_Jupiter/(-a_plus)^3);
h_plus = cross(r_p_vect, vp_plus_v);
p_plus = norm(h_plus)^2/mu_Jupiter;
theta = acos((p_plus/r_soi-1)/e_plus);
calc = sqrt(-(1+e_plus)/(1-e_plus));
F = 2*atanh(tan(theta/2)/calc);
t_plus = (e_plus*sinh(F)-F)/n_plus;

% total time from the entrance to the exit of the SOI( in seconds)
t_SOI=t_minus+t_plus;
t_SOI_days = t_SOI/(3600*24);
%% Direct interplanetary trajectory fom Saturn to Mars

% t_span definition
imax=3000;
jmax=1500;
% we know, after we ran the code, the time to arrive to Mars with the minimum Dv is about 5 years,
% so we can eliminate at least the last 6 years of possible departure from Saturn to
% increase the velocity of the code itself, avoiding useless departure
% windows
t_span1=linspace(date2mjd2000([2028,1,1,0,0,0]),date2mjd2000([2053,1,1,0,0,0]),imax); 

% pre-allocation
r_sat=zeros(3,length(t_span1));
v_sat=zeros(3,length(t_span1));
r_mar=zeros(3,length(t_span1),jmax);
v_mar=zeros(3,length(t_span1),jmax);
kep_sat=zeros(6,length(t_span1));
kep_mar=zeros(6,length(t_span1),jmax);
Dv_mar_jup=inf*ones(imax,jmax);
t_mat1=zeros(imax,jmax);
v1_tr1=zeros(3,length(t_span1),jmax);
v2_tr1=zeros(3,length(t_span1),jmax);

for i=1:length(t_span1)
    for j=1:jmax

        % using ephemerides 
        [kep_sat(:,i),ksun]=uplanet(t_span1(i),6);    
        %Mars position every 12h to be precise
        [kep_mar(:,i,j),~]=uplanet(t_span1(i)+5*365+j/2,4); 
        [r_sat(:,i),v_sat(:,i)]=kep2car(kep_sat(1,i),kep_sat(2,i),kep_sat(3,i),kep_sat(4,i),kep_sat(5,i),kep_sat(6,i),ksun);
        [r_mar(:,i,j),v_mar(:,i,j)]=kep2car(kep_mar(1,i,j),kep_mar(2,i,j),kep_mar(3,i,j),kep_mar(4,i,j),kep_mar(5,i,j),kep_mar(6,i,j),ksun);

        % time for each transfer in seconds
        t_mat1(i,j)=(5*365+j/2)*3600*24;    

        % Lambert function 
        [~,~,~,~,v1_tr1(:,i,j),v2_tr1(:,i,j),~,~]=lambertMR(r_sat(:,i),r_mar(:,i,j),t_mat1(i,j),ksun,0,0,0,0);
        Dv_mar_jup(i,j)=norm(v_sat(:,i)-v1_tr1(:,i,j))+norm(v_mar(:,i,j)-v2_tr1(:,i,j));

    end
end

% minimum Dv and its index calculation
dv_min_noflyby = min(min(Dv_mar_jup));
[v,loc] = min(Dv_mar_jup(:));
[i,j] = ind2sub(size(Dv_mar_jup),loc);

%% Plot of the direct transfer
% Position of the planets in best case scenario
t_dep_noflyby = t_span1(i);
tof = t_mat1(i, j)/(3600*24);
[kep_1,ksun] = uplanet(t_dep_noflyby, Saturn);
[r_1, v_1] = kep2carRAD(kep_1, ksun);
[kep_2,ksun] = uplanet(t_dep_noflyby + +tof, Mars);
[r_2, v_2] = kep2carRAD(kep_2, ksun);

% Plot of the various planets and the Sun
opts.Units = 'm'; 
opts.Position = [r_1(1), r_1(2), r_1(3)];
planet3D('Saturn', opts);  
hold on;
grid on;
opts.Units = 'mars'; 
opts.Position = [r_2(1), r_2(2), r_2(3)];
planet3D('Mars', opts);  
opts.Units = 'sun';
opts.Position = [0 0 0];
planet3D('Sun', opts); 
% Plot of the Lambert's arcs
dv_check = interplanetary_transfer_cost(t_dep_noflyby, t_dep_noflyby+tof, Saturn, Mars);
legend('', '', '', 'Orbit of Saturn', 'Orbit of Mars', 'First Lambert arc')
title('Direct transfer')