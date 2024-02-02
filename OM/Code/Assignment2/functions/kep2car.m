function [r, v] = kep2car(a, e, i, OM, om, th, mu)
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates.
% 
% PROTOTYPE:
%  [r, v] = kep2car(a, e, i, OM, om, th, mu)
% 
% INPUTS:
%   a           [1x1]   Semi-major axis                  [km]
%   e           [1x1]   Eccentricity                     [-]
%   i           [1x1]   Inclination                      [rad]
%   OM          [1x1]   RAAN                             [rad]
%   om          [1x1]   Pericenter anomaly               [rad]
%   th          [1x1]   True anomaly                     [rad]
%   mu          [1x1]   Gravitational parameter          [km^3/s^2]
%
% OUTPUTS:
%   r           [3x1]   Position vector                  [km]
%   v           [3x1]   Velocity vector                  [km/s]
% 
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascharetti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

if nargin == 6
    mu = 398600.433;
end

p = a * (1 - e^2);
r_norm = p / (1 + e * cos(th));

% Position and velocity in perifocal system
r_PF = r_norm * [cos(th) sin(th) 0]';
v_PF = sqrt(mu/p) * [-sin(th) (e + cos(th)) 0]';

% Rotational matrices PF --> ECI

% 1. Rotation of ω around the k'' axis
R3_om = [ cos(om) sin(om) 0;
         -sin(om) cos(om) 0;
                0       0 1];
% 2. Rotation of i around the i' axis
R1_i = [1      0      0;
        0  cos(i) sin(i);
        0 -sin(i) cos(i)];
% 3. Rotation of Ω around the k axis
R3_OM = [ cos(OM) sin(OM) 0;
         -sin(OM) cos(OM) 0;
                0       0 1];

T_ECI_PF = R3_om * R1_i * R3_OM;
T_PF_ECI = T_ECI_PF';

r = T_PF_ECI * r_PF;
v = T_PF_ECI * v_PF;

return