function [alpha,delta,long,lat]=groundTrack(r_v_vect,kep0,thG0,parameters,k,m)

%works well but y-axis go from 90 to -90

a = kep0(1);
e = kep0(2);
in = kep0(3);
OM = kep0(4);
om = kep0(5);
th = kep0(6);

mu = parameters.mu;
wE = parameters.wE;
period = 2*pi*sqrt(a^3/mu);
period = linspace(0, period, 1000);

% UNPERTURBED ORBIT GROUND TRACK - 1 ORBIT

%plot map of the Earth
figure()
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

%r and v at t0
r0 = r_v_vect(1:3);
v0 = r_v_vect(4:6);

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
plot(long(1),lat(1),'m',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'k',Marker='*',LineWidth=3,LineStyle='none');

%% UNPERTURBED ORBIT REPEATING GROUND TRACK

n = wE*k/m;               %[rad/s]
a_rgt = (mu/n^2)^(1/3);   %[km]
[r_rgt,v_rgt]=kep2car(a_rgt, e, in, OM, om, th, mu);

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0_rgt=[r_rgt;v_rgt];
odefun=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y]=ode89(odefun,period,y0_rgt,options);

lat = zeros(size(Y,1),1);
alpha = zeros(size(Y,1),1);
thG = zeros(size(Y,1),1);
long = zeros(size(Y,1),1);

%longitud and latitude calculation
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
plot(long(1),lat(1),'y',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='*',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Unperturbed orbit ground tracks - One orbit')
legend('Ground track','Starting point','Ending point','Repeating ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% UNPERTURBED ORBIT GROUND TRACK - 1 DAY
one_day = linspace(0, 24*60*60, 30000);

%plot map of the Earth
figure()
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
plot(long,lat,'r',LineStyle='none',Marker='.');
plot(long(1),lat(1),'m',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'k',Marker='*',LineWidth=3,LineStyle='none');

%% UNPERTURBED ORBIT REPEATING GROUND TRACK - 1 DAY

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
odefun=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y]=ode89(odefun,one_day,y0_rgt,options);

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
plot(long(1),lat(1),'y',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='*',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Unperturbed orbit ground tracks - One day')
legend('Ground track','Starting point','Ending point','Repeating ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% PERTURBED ORBIT GROUND TRACK - 1 ORBIT

%plot map of the Earth
figure()
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
plot(long,lat,'r',LineStyle='none',Marker='.');
plot(long(1),lat(1),'m',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'k',Marker='*',LineWidth=3,LineStyle='none');

%% PERTURBED ORBIT REPEATING GROUND TRACK - 1 ORBIT

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
[T, Y] = ode113( @(t,y) eq_motion_CAR( t, y, @(t,y) acc_pert_fun_CAR(t,y,parameters), parameters ), period, y0_rgt, options);

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
plot(long(1),lat(1),'y',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='*',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Perturbed orbit ground tracks - One orbit')
legend('Ground track','Starting point','Ending point','Repeating ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% PERTURBED ORBIT GROUND TRACK - 1 DAY

%plot map of the Earth
figure()
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
plot(long,lat,'r',LineStyle='none',Marker='.');
plot(long(1),lat(1),'m',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'k',Marker='*',LineWidth=3,LineStyle='none');

%% PERTURBED ORBIT REPEATING GROUND TRACK - 1 DAY

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
[T, Y] = ode113( @(t,y) eq_motion_CAR( t, y, @(t,y) acc_pert_fun_CAR(t,y,parameters), parameters ), one_day, y0_rgt, options);

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
plot(long(1),lat(1),'y',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'c',Marker='*',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
title('Perturbed orbit ground tracks - One day')
legend('Ground track','Starting point','Ending point','Repeating ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');