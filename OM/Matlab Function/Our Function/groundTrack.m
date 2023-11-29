function [alpha,delta,long,lat]=groundTrack(r_v_vect,thG0,t_vect,mu,wE)
%works well but y-axis go from 90 to -90

%plot map of the Earth
A=imread("EarthTexture.jpg");
image([-180 180],[-90 90],A);
hold on;

%r and v at t0
r0=r_v_vect(1:3);
v0=r_v_vect(4:6);

%orbit and r,v calculation in cartesian coordinates
options=odeset('RelTol',1e-13,'AbsTol',1e-14); %options must be specified
y0=[r0;v0];
h0=cross(y0(1:3),y0(4:6));
e0=cross((1/mu*y0(4:6)),h0)-y0(1:3)/norm(y0(1:3));
p=norm(h0)^2/mu;
a=p/(1-norm(e0)^2);
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
plot(long,lat,'r',LineStyle='none',Marker='.');
plot(long(1),lat(1),'b',Marker='*',LineWidth=3,LineStyle='none');
plot(long(end),lat(end),'g',Marker='*',LineWidth=3,LineStyle='none');
axis([-180 180 -90 90]);
legend('track ground','starting point','ending point',fontsize=15);
xlabel('longitude');
ylabel('latitude');


return
