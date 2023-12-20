function [alpha,delta,long,lat]=groundTrack(r_v_vect,kep0,thG0,t_vect,mu,wE,k,m)
%works well but y-axis go from 90 to -90

a = kep0(1);
e = kep0(2);
in = kep0(3);
OM = kep0(4);
om = kep0(5);
th = kep0(6);

% UNPERTURBED ORBIT GROUND TRACK

%plot map of the Earth
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);

%r and v at t0
r0 = r_v_vect(1:3);
v0 = r_v_vect(4:6);

%orbit and r,v calculation in cartesian coordinates
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0 = [r0;v0];
odefun = @(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y] = ode89(odefun,t_vect,y0,options);

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
plot(long,lat,'r',LineStyle='none',Marker='.');
plot(long(1),lat(1),'b',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'k',Marker='*',LineWidth=3,LineStyle='none');

%% UNPERTURBED ORBIT REPEATING GROUND TRACK

n = wE*k/m;               %[rad/s]
a_rgt = (mu/n^2)^(1/3);   %[km]
[r,v]=kep2car(a_rgt, e, in, OM, om, th);
r_v_vect=[r;v];         %[m;m/s]

%r and v at t0
r0=r_v_vect(1:3);
v0=r_v_vect(4:6);

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0=[r0;v0];
odefun=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
[T,Y]=ode89(odefun,t_vect,y0,options);

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
title('Unperturbed orbit ground tracks')
legend('Ground track','Starting point','Ending point','Repeating ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%% PERTURBED ORBIT GROUND TRACK

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
wE = [0,0,wE]';
rE = astroConstants(23);
j2 = astroConstants(9);
CD = 2.1;
AM = 0.0043;
num_orb = 1000;
odefun_pert=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)-0.5*rho_find(norm(y(1:3))-rE)*norm(y(4:6)-cross(wE,y(1:3)))*CD*AM*(y(4:6)-cross(wE,y(1:3)))+3/2*(j2*mu*rE^2)/(norm(y(1:3))^4)*[y(1)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-1);y(2)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-1);y(3)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-3)]]; 
[T,Y]=ode89(odefun_pert,linspace(0,num_orb*2*pi*sqrt(a^3/mu),10000),y0,options);

%longitude and latitude calculation
t_sid = 23*60*60 + 56*60 +4;
wE = 2*pi/t_sid;

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
plot(long(1),lat(1),'b',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'k',Marker='*',LineWidth=3,LineStyle='none');

%% PERTURBED ORBIT REPEATING GROUND TRACK

n = wE*k/m;               %[rad/s]
a_rgt = (mu/n^2)^(1/3);   %[km]
[r,v]=kep2car(a_rgt, e, in, OM, om, th);
r_v_vect=[r;v];         %[m;m/s]

%r and v at t0
r0=r_v_vect(1:3);
v0=r_v_vect(4:6);

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0=[r0;v0];
wE = [0,0,wE]';
odefun_pert=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)-0.5*rho_find(norm(y(1:3))-rE)*norm(y(4:6)-cross(wE,y(1:3)))*CD*AM*(y(4:6)-cross(wE,y(1:3)))+3/2*(j2*mu*rE^2)/(norm(y(1:3))^4)*[y(1)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-1);y(2)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-1);y(3)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-3)]]; 
[T,Y]=ode89(odefun_pert,linspace(0,num_orb*2*pi*sqrt(a^3/mu),10000),y0,options);

%longitude and latitude calculation
t_sid = 23*60*60 + 56*60 +4;
wE = 2*pi/t_sid;
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
title('Perturbed orbit ground tracks')
legend('Ground track','Starting point','Ending point','Repeating ground track','Starting point','Ending point');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');