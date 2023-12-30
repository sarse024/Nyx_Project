function [alpha,delta,long,lat]=ground_Track(s0,kep0,thG0,parameters,k,m)
% Computation and plot of the ground track for unperturbed and perturbed 
% orbits, with different time spans and keplerian elements

% INPUT: 
% s0[6x1]                             Cartesian initial state of the satellite 
%                                     (rx, ry, rz, vx, vy, vz) [km, km/s]
% kep0[1x6]                           Keplerian elements (a, e, i, OM, om, th) 
%                                     [km, -, rad, rad, rad, rad, rad]
% thG0[1]                             Longitude of Greenwich meridian at initial time
% parameters[1x1 struct]              Set of constant parameters (rE, om_E,
%                                     mu, J2, CD, AM) [km, rad/s, km^3/s^2,
%                                     -, -, m^2/kg]
% k[1]                                Number of revolutions of the
%                                     satellite [-]
% m[1]                                Rotation of the Earth [-]
%
% OUTPUT:
% alpha[length(tspan)x1]              Right ascention [deg]
% delta[length(tspan)x1]              Declination [deg]
% long[length(tspan)x1]               Longitude [deg]
% lat[length(tspan)x1]                Latitude [deg]
%
% CONTRIBUTORS:
% Valentina D’Annunzio
% Mirko Mascaretti
% Edoardo Nicolucci Balocco
% Samuele Orsenigo
%
% 2023-12-28: First version

a = kep0(1);
e = kep0(2);
in = kep0(3);
OM = kep0(4);
om = kep0(5);
th = kep0(6);

mu = parameters.mu;
wE = parameters.wE;
T_period = 2*pi*sqrt(a^3/mu);
period = linspace(0, T_period, 5000);

%% UNPERTURBED ORBIT GROUND TRACK - 1 ORBIT

%plot map of the Earth
figure()
tiledlayout(2,2)
nexttile
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

%r and v at t0
r0 = s0(1:3);
v0 = s0(4:6);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0 = [r0;v0];
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,period,y0,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right rappresentable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'c',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Unperturbed orbit ground track - One orbit')
legend('Ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% UNPERTURBED ORBIT GROUND TRACK - 1 DAY

one_day = linspace(0, 24*60*60, 30000);

%plot map of the Earth
nexttile
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,one_day,y0,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right rappresentable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'c',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Unperturbed orbit ground track - One day')
legend('Ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% PERTURBED ORBIT GROUND TRACK - 1 ORBIT

%plot map of the Earth
nexttile
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
[T, Y] = ode113( @(t,y) eq_motion_CAR( t, y, @(t,y) acc_pert_fun_CAR(t,y,parameters), parameters ), period, y0, options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right representable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'c',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Perturbed orbit ground track - One orbit')
legend('Ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% PERTURBED ORBIT GROUND TRACK - 1 DAY

%plot map of the Earth
nexttile
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
[T, Y] = ode113( @(t,y) eq_motion_CAR( t, y, @(t,y) acc_pert_fun_CAR(t,y,parameters), parameters ), one_day, y0, options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right representable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'c',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Perturbed orbit ground track - One day')
legend('Ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% UNPERTURBED ORBIT - ONE ORBIT - DIFFERENT RAAN

OM = deg2rad(100);

%plot map of the Earth
figure()
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

[r0,v0]=kep2car(a, e, in, OM, om, th, mu);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0 = [r0;v0];
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,period,y0,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right rappresentable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'r',LineStyle='none',Marker='.');
plot(long(1),lat(1),'r',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'r',Marker='diamond',LineWidth=3,LineStyle='none');

%%

OM = deg2rad(200);

[r0,v0]=kep2car(a, e, in, OM, om, th, mu);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0 = [r0;v0];
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,period,y0,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right rappresentable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'g',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'g',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Different \Omega')
legend('\Omega = 100\circ','Starting point','Ending point','\Omega = 200\circ','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% UNPERTURBED ORBIT GROUND TRACK - ONE ORBIT - DIFFERENT ARGUMENT OF PERICENTER
a = kep0(1);
e = kep0(2);
in = kep0(3);
OM = kep0(4);
om = deg2rad(100);
th = kep0(6);

%plot map of the Earth
figure()
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

[r0,v0]=kep2car(a, e, in, OM, om, th, mu);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0 = [r0;v0];
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,period,y0,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right rappresentable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'r',LineStyle='none',Marker='.');
plot(long(1),lat(1),'r',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'r',Marker='diamond',LineWidth=3,LineStyle='none');

%%

om = deg2rad(200);

[r0,v0]=kep2car(a, e, in, OM, om, th, mu);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0 = [r0;v0];
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,period,y0,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right rappresentable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'g',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'g',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Different \omega')
legend('\omega = 100\circ','Starting point','Ending point','\omega = 200\circ','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% UNPERTURBED ORBIT GROUND TRACK - ONE ORBIT - DIFFERENT TRUE ANOMALY

a = kep0(1);
e = kep0(2);
in = kep0(3);
OM = kep0(4);
om = kep0(5);
th = deg2rad(0);

%plot map of the Earth
figure()
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

[r0,v0]=kep2car(a, e, in, OM, om, th, mu);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0 = [r0;v0];
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,period,y0,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right rappresentable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'r',LineStyle='none',Marker='.');
plot(long(1),lat(1),'r',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'r',Marker='diamond',LineWidth=3,LineStyle='none');

%%

th = deg2rad(100);

[r0,v0]=kep2car(a, e, in, OM, om, th, mu);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0 = [r0;v0];
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,period,y0,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right rappresentable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'g',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'g',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Different \theta')
legend('\theta = 0\circ','Starting point','Ending point','\theta = 100\circ','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% UNPERTURBED ORBIT REPEATING GROUND TRACK

a = kep0(1);
e = kep0(2);
in = kep0(3);
OM = kep0(4);
om = kep0(5);
th = kep0(6);

n = wE*k/m;               %[rad/s]
a_rgt = (mu/n^2)^(1/3);   %[km]
[r_rgt,v_rgt]=kep2car(a_rgt, e, in, OM, om, th, mu);
T_rep = (2*pi/wE)*(m/k);
tspan = linspace(0, k*T_rep, 50000);

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0_rgt=[r_rgt;v_rgt];
odefun=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y]=ode89(odefun,tspan,y0_rgt,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right representable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot map of the Earth
figure()
tiledlayout(2,1)
nexttile
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'c',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Unperturbed orbit repeating ground track')
legend('Ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% PERTURBED ORBIT REPEATING GROUND TRACK

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
[T, Y] = ode113( @(t,y) eq_motion_CAR( t, y, @(t,y) acc_pert_fun_CAR(t,y,parameters), parameters ), tspan, y0_rgt, options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitude and latitude calculation
for i=1:size(Y,1)
    r=norm(Y(i,1:3));
    lat(i)=asin(Y(i,3)./r);
    
    alpha(i)=atan2(Y(i,2),Y(i,1));
    thG(i)=thG0+wE*T(i);
    long(i)=alpha(i)-thG(i);

    %transport longitude and latitude in the right representable values
    while long(i)<-pi || long(i)>pi
        if long(i)<-pi 
            long(i)=long(i)+2*pi;
        else
            long(i)=long(i)-2*pi;
        end
    end
end

%convert radiants to degrees
long=long./pi*180;   
lat=-lat./pi*180;   %- because the y axis has negative numbers on the top

delta=lat;          %because same equatorial plane of reference

%plot map of the Earth
nexttile
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

%plot the ground track
hold on
plot(long,lat,'g',LineStyle='none',Marker='.');
plot(long(1),lat(1),'c',Marker='o',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='diamond',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Perturbed orbit repeating ground track')
legend('Ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');