function [r,v] = kep2carRAD(kep, mu)
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% radians.
%
% INPUT:
% kep  [1x6] that contains:
%  - a  [1x1] Semi-major axis       [km]
%  - e  [1x1] Eccentricity          [-]
%  - i  [1x1] Inclination           [rad]
%  - OM [1x1] RAAN                  [rad]
%  - om [1x1] Pericentre anomaly    [rad]
%  - th [1x1] True anomaly          [rad]
% mu   [1x1] Gravitational parameter  [km^3/s^2]
%
% OUTPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
%
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascharetti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

if nargin == 6
    mu = 398600.433;
end

a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
th = kep(6);

toll = 1e-6;

% check the convention
if (abs(i) <= toll && abs(OM) >= toll)
    OM = 0; %avoid singularity
    fprintf('\n\nWarning: update OM = 0 to follow the convention for i = 0\n\n')
end
if (abs(e) <= toll && abs(om) >= toll)
    om = 0; %avoid singularity
    fprintf('\n\nWarning: update om = 0 to follow the convention for e = 0\n\n')
end

p=a*(1-e^2);    

%modulo vettore posizione
r_mod=p/(1+e*cos(th));      

%vettori posizione e velocità nel sistema PF
r_PF=r_mod*[cos(th),sin(th),0];            
v_PF=sqrt(mu/p)*[-sin(th),e+cos(th),0];

%matrice di rotazione da PF a ECI
R_trasp=[cos(om)*cos(OM)-sin(om)*sin(OM)*cos(i), ...  
    -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM), sin(i)*sin(OM); sin(OM)*cos(om)+...
    sin(om)*cos(OM)*cos(i),-sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM), -sin(i)*cos(OM);...
    sin(om)*sin(i), cos(om)*sin(i), cos(i)];

%vettori posizione e velocità nel sistema ECI
r=R_trasp*r_PF';         
v=R_trasp*v_PF';

end