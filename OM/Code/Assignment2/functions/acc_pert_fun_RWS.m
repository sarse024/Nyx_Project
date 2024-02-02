function acc_pert_vec = acc_pert_fun_RWS(~, s, parameters)
% acc_pert_fun_RWS.m - Calculation total acceleration due to perturbation
% in RWS frame. Perturbance take in consideration are J2 and atmospheric
% drag
%
% PROTOTYPE:
% acc_pert_vec = acc_pert_fun_RWS(t, s, parameters)
%
% DESCRIPTION:
% Calculation total acceleration due to perturbation
% in RWS frame. Perturbance take in consideration are J2 and atmospheric
% drag. Useful to Gauss Integration.
%
% INPUT:
% t [1x1] Time                      [s]
% s [6x1] Keplerian vector:
%  - a  [1x1] Semi-major axis       [km]
%  - e  [1x1] Eccentricity          [-]
%  - i  [1x1] Inclination           [rad]
%  - OM [1x1] RAAN                  [rad]
%  - om [1x1] Pericentre anomaly    [rad]
%  - th [1x1] True anomaly          [rad]
% parameters: Struct with information 
%   - rE        [1x1]   Radius of Earth                                 [Km]   
%   - wE        [1x1]   Angular velocity of Earth                       [rad/s]
%   - mu        [1x1]   Gravitational parameter                         [km^3/s^2]
%   - drag.CD   [1x1]   Drag coefficient                                [-]
%   - drag.AM   [1x1]   A/m ratio between cross area A and mass of S/C  [m^2/kg]
%   - j2        [1x1]   Second zonal harmonic potential                 [-]
% 
% OUTPUT:
% acc_pert_vec  [3x1] Total acceleration        [Km/s^2]
%
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascheretti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

% extract general parameters
rE = parameters.rE;              %[Km]
wE = [0, 0, parameters.wE]';     % rad/s;
mu = parameters.mu;
a = s(1);
e = s(2);
i = s(3);
OM = s(4);
om = s(5);
th = s(6);

p = a*(1-e^2);
h = sqrt(p*mu);

% position of satellite and velocity
[r_sc, vel_sc] = kep2car(a,e,i,OM,om,th,mu);
altitude = norm(r_sc) - rE;      %[Km] altitude of spacecraft 

% DRAG MODELLING
CD = parameters.drag.CD;         %[]
AM = parameters.drag.AM;         %[m^2/kg]

% relative velocity between sc and air (neglecting wind)
vel_relative = vel_sc - cross(wE,r_sc);  

vel_relative = [norm(vel_relative), 0, 0]'; %change velocity in TNH
uv = vel_relative/norm(vel_relative); % TNH velocity versor

% calculate acc due to drag in tangential normal refernce frame (TNH)
acc_drag = - 0.5*rho_find(altitude)*(norm(vel_relative)*1000)^2*CD*AM*uv; 

% TNH to Radial trasversal out of plane frame (RSW)
v = norm(vel_sc);
Rot_mat = h/(p*v) * [e*sin(th), -(1+e*cos(th)), 0; (1+e*cos(th)), e*sin(th), 0; 0, 0, 1]; %from tnh to RSW
acc_drag = Rot_mat*acc_drag;

% J2 modelling in RSW rf
j2 = parameters.j2;

acc_j2 = [1-3*sin(i)^2*sin(th+om)^2;
          sin(i)^2*sin(2*(th+om));
          sin(2*i)*sin(th+om)];
acc_j2 = -3/2*j2*mu*rE^2/norm(r_sc)^4*acc_j2; 

% sum of all perturbations in RWS rf
acc_pert_vec = acc_j2 + acc_drag/1000;

end