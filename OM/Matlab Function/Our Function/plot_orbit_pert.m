function plot_orbit_pert(mu,rE,j2,CD,AM,y0,num_orb)

%usare patch al posto di plot per colormap

%procedure
h0=cross(y0(1:3),y0(4:6));
e0=cross((1/mu*y0(4:6)),h0)-y0(1:3)/norm(y0(1:3));
p=norm(h0)^2/mu;
a=p/(1-norm(e0)^2);
odefun_pert=@(t,y) [y(4:6);-mu/norm(y(1:3))^3*y(1:3)-0.5*rho_find(norm(y(1:3))-rE)*norm(y(4:6))*CD*AM*y(4:6)]; 
options=odeset('RelTol',1e-13,'AbsTol',1e-14);
%+3/2*j2*mu*rE^2/norm(y(1:3))^4*[y(1)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-1);y(2)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-1);y(3)/norm(y(1:3))*(5*y(3)^2/norm(y(1:3))^2-3)]
%plot
Terra3d;
xlabel('X[km]');
ylabel('Y[km]');
zlabel('Z[km]');
axis equal;
grid on;
title('perturbated orbit (oblateness+drag effect)');
[~,Y]=ode89(odefun_pert,linspace(0,num_orb*2*pi*sqrt(a^3/mu),10000),y0,options);
plot3(Y(:,1),Y(:,2),Y(:,3),'-');
 


