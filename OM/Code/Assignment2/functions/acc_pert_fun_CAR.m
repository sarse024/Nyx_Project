function [acc_pert_vec, acc_drag, acc_j2] = acc_pert_fun_CAR(~, s, parameters )
% acc_pert_fun_CAR.m - Calculation total acceleration due to perturbation
% in ECI frame. Perturbance take in consideration are J2 and atmospheric
% drag
%
% PROTOTYPE:
% [acc_pert_vec, acc_drag, acc_j2] = acc_pert_fun_CAR(t, s, parameters)
%
% DESCRIPTION:
% Calculation total acceleration due to perturbation
% in ECI frame. Perturbance take in consideration are J2 and atmospheric
% drag. Useful to Cartesian Integration.
%
% INPUT:
% t [1x1] Time                                  [s]
% s [6x1] State vector:
%    - r     [3x1] Position vector                [km]
%    - v     [3x1] Velocity vector                [km/s]
% parameters
%
% OUTPUT:
% acc_pert_vec  [3x1] Total acceleration        [Km/s^2]
% acc_drag      [3x1] Acceleration due to drag  [Km/s^2]
% acc_j2        [3x1] Acceleration due to J2    [Km/s^2]
% parameters: Struct with information 
%   - rE        [1x1]   Radius of Earth                                 [Km]   
%   - wE        [1x1]   Angular velocity of Earth                       [rad/s]
%   - mu        [1x1]   Gravitational parameter                         [km^3/s^2]
%   - drag.CD   [1x1]   Drag coefficient                                [-]
%   - drag.AM   [1x1]   A/m ratio between cross area A and mass of S/C  [m^2/kg]
%   - j2        [1x1]   Second zonal harmonic potential                 [-]    
%
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascheretti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

% extract general parameters
rE = parameters.rE;              %[Km]
wE = [0, 0, parameters.wE]';              % rad/s;
mu = parameters.mu;

% DRAG MODELLING
CD = parameters.drag.CD;         %[]
AM = parameters.drag.AM;         %[m^2/kg]

% position of satellite and velocity
r_sc = s(1:3);            %[Km] position fo spacecraft
h = norm(r_sc) - rE;      %[Km] altitude of spacecraft 
vel_sc = s(4:6);          %[Km/s] velocity of spacecraft

% relative velocity between sc and air (neglecting wind)
vel_relative = vel_sc - cross(wE,r_sc); % ECI velocity
uv = vel_relative/norm(vel_relative);   % ECI velocity versor

% calculate acc due to drag in tangential normal reference frame (ECI)
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