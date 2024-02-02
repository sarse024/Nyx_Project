function ds = eq_motion_GAUSS(t, s, ~, parameters)
% eq_motion_GAUSS.m - Calculation of derivative functions in RWS frame to integrate with
% ode. Perturbance take in consideration are J2 and atmospheric
% drag
%
% PROTOTYPE:
% ds = eq_motion_GAUSS(t, s, acc_per_fun_RWS, parameters)
%
% DESCRIPTION:
% Calculation of derivative functions in RWS frame to integrate with
% ode. Perturbance take in consideration are J2 and atmospheric
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
% acc_per_fun_RWS [3x1] Calculate acceleration of disturbances  [Km/s^2]
% parameters: is a struct with information 
%   - rE        [1x1]   Radius of Earth                                 [Km]   
%   - wE        [1x1]   Angular velocity of Earth                       [rad/s]
%   - mu        [1x1]   Gravitational parameter                         [km^3/s^2]
%   - drag.CD   [1x1]   Drag coefficient                                [-]
%   - drag.AM   [1x1]   A/m ratio between cross area A and mass of S/C  [m^2/kg]
%   - j2        [1x1]   Second zonal harmonic potential                 [-]
%
% OUTPUT:
% ds [6x1] Derivative Keplerian variable function to integrate with Gauss
%
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascheretti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

mu = parameters.mu;

a = s(1);
e = s(2);
i = s(3);
%OM = s(4);
om = s(5);
th = s(6);

p = a*(1-e^2);
h = sqrt(p*mu);
r = p/(1+e*cos(th));

% Evaluate the perturbing accelerations
acc_pert_vec = acc_pert_fun_RWS(t, s, parameters); %give acc in RWS form

acc_r = acc_pert_vec(1);
acc_s = acc_pert_vec(2);
acc_w = acc_pert_vec(3);

% Evaluate the equations of motion (Cartesian or Keplerian),
ds = [2*a^2/h*(e*sin(th)*acc_r + p/r*acc_s);
      1/h*(p*sin(th)*acc_r + ((p+r)*cos(th) + r*e)*acc_s);
      (r*cos(th+om))/h*acc_w;
      (r*sin(th+om))/(h*sin(i))*acc_w;
      1/(h*e)*(-p*cos(th)*acc_r + (p+r)*sin(th)*acc_s) - (r*sin(th+om)*cos(i))/(h*sin(i))*acc_w;
      h/r^2 + 1/(e*h)*(p*cos(th)*acc_r - (p+r)*sin(th)*acc_s)];

end