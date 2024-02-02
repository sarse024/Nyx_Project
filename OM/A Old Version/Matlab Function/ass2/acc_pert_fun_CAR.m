function [acc_pert_vec, acc_drag, acc_j2] = acc_pert_fun_CAR(t, s, parameters )
% Evaluate the perturbing accelerations in a given
% reference frame (e.g., TNH, RTH, ECI, etc.)

% extract general parameters
rE = parameters.rE;              %[Km]
wE = [0, 0, parameters.wE]';              % rad/s;
mu = parameters.mu;
kep = parameters.kep;  %vactor of keplerian
e = kep(2); %eccentricity
i = kep(3); %true anomaly [rad]
OM = kep(4); %true anomaly [rad]
om = kep(5); %true anomaly [rad]
th = kep(6); %true anomaly [rad]


% DRAG MODELLING
CD = parameters.drag.CD;         %[]
AM = parameters.drag.AM;         %[m^2/kg]


% position of satellite and velocity
r_sc = s(1:3);            %[Km] position fo spacecraft
h = norm(r_sc) - rE;      %[Km] altitude of spacecraft 
vel_sc = s(4:6);          %[Km/s] velocity of spacecraft

% relative velocity between sc and air (neglecting wind)
vel_relative = vel_sc - cross(wE,r_sc);    %DA CAMBIARE WE è UN VETTORE STRANO
uv = vel_relative/norm(vel_relative);

% calculate acc due to drag in tangential normal refernce frame (TNH)
acc_drag = - 0.5*rho_find(h)*(norm(vel_relative)*1000)^2*CD*AM*uv; 

% J2 modelling in Cartesian rf
j2 = parameters.j2;

acc_j2 = [r_sc(1)/norm(r_sc)*(5*r_sc(3)^2/norm(r_sc)^2-1);
          r_sc(2)/norm(r_sc)*(5*r_sc(3)^2/norm(r_sc)^2-1);
          r_sc(3)/norm(r_sc)*(5*r_sc(3)^2/norm(r_sc)^2-3)];
acc_j2 = 3/2*j2*mu*rE^2/norm(r_sc)^4*acc_j2; 

% sum of all perturbations in Cartesian rf
 acc_pert_vec = acc_j2 + acc_drag/1000;

end