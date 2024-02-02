function ds = eq_motion_CAR( t, s, ~, parameters)
% eq_motion_CAR.m - Calculation of derivative functions in ECI frame to integrate with
% ode. Perturbance take in consideration are J2 and atmospheric
% drag
%
% PROTOTYPE:
% ds = eq_motion_CAR( t, s, acc_per_fun_CAR, parameters)
%
% DESCRIPTION:
% Calculation of derivative functions in ECI frame to integrate with
% ode. Perturbance take in consideration are J2 and atmospheric
% drag. Useful to Cartesian Integration.
%
% INPUT:
% t [1x1] time                      [s]
% s [6x1] state vector:
%  - a  [1x1] Semi-major axis       [km]
%  - e  [1x1] Eccentricity          [-]
%  - i  [1x1] Inclination           [rad]
%  - OM [1x1] RAAN                  [rad]
%  - om [1x1] Pericentre anomaly    [rad]
%  - th [1x1] True anomaly          [rad]
% acc_per_fun_CAR [3x1] Calculate acceleration of disturbances  [Km/s^2]
% parameters: Struct with information
%   - rE        [1x1]   Radius of Earth                                 [Km]   
%   - wE        [1x1]   Angular velocity of Earth                       [rad/s]
%   - mu        [1x1]   Gravitational parameter                         [km^3/s^2]
%   - drag.CD   [1x1]   Drag coefficient                                [-]
%   - drag.AM   [1x1]   A/m ratio between cross area A and mass of S/C  [m^2/kg]
%   - j2        [1x1]   Second zonal harmonic potential                 [-]
%
% OUTPUT:
% ds [6x1] Derivative state vector function to integrate
%
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascheretti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

% Estract usefull parameters
mu = parameters.mu;

% Evaluate the perturbing accelerations
acc_pert_vec = acc_pert_fun_CAR( t, s, parameters);

% Evaluate the equations of motion (Cartesian or Keplerian),
ds = [s(4:6);-mu/norm(s(1:3))^3*s(1:3) + acc_pert_vec];

end