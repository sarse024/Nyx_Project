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
a=6846;         %[km]
e=0.0298;       %[]
i=80.2068;      %[deg]
k=15;           %[]
m=1;            %[]

%drag parameters
CD=2.1;         %[]
AM=0.0043;      %[m^2/kg]

%other keplerian elements (arbitrary)
OM=0;
om=0;
th=0;

%pert
mu=astroConstants(13);
rE=astroConstants(23);
j2=0.00108263;
[r,v]=kep2car(a,e,i/180*pi,OM,om,th,astroConstants(13));    %[m;m/s]
r_v_vect=[r;v];         
plot_orbit_pert(mu,rE,j2,CD,AM,r_v_vect,1000);
%%
%unperturbed orbit plot
mu=astroConstants(13);  %[km^3/s^2]
[r,v]=kep2car(a,e,i/180*pi,OM,om,th,mu);
r_v_vect=[r;v];         %[m;m/s]
odefun=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)]; 
options=odeset('RelTol',1e-13,'AbsTol',1e-14);
[T,Y]=ode89(odefun,[0,2*pi*sqrt(a^3/mu)],r_v_vect,options);
figure()
Terra3d;
hold on;
plot3(Y(:,1),Y(:,2),Y(:,3),'-',LineWidth=3);
xlabel('X[km]');
ylabel('Y[km]');
zlabel('Z[km]');
axis equal;
grid on;

%unperturbed orbit ground track
wE=2*pi/86164;          %[rad/s]
[r,v]=kep2car(a,e,i/180*pi,OM,om,th,mu);
r_v_vect=[r;v];         %[m;m/s]
thG0=0;                 %[rad]
t_vect=0:24*3600;       %[s]
figure();
[alpha,delta,long,lat]=groundTrack(r_v_vect,thG0,t_vect,mu,wE);

%if repeating ground track
n=wE*k/m;               %[rad/s]
a_rgt=(mu/n^2)^(1/3);   %[km]
[r,v]=kep2car(a_rgt,e,i/180*pi,OM,om,th,mu);
r_v_vect_rgt=[r;v];     %[m;m/s]
thG0=0;                 %[rad]
t_vect=0:86164;         %[s]
figure();
[alpha_rgt,delta_rgt,long_rgt,lat_rgt]=groundTrack(r_v_vect_rgt,thG0,t_vect,mu,wE);


