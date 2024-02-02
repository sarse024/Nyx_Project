function R = rotationMatrix(i,OM,om)
% rotationMatrix.m - Function to calculate the rotation matrix to work on the orbital plane.
%
% PROTOTYPE:
% R = rotationMatrix(i,OM,om)
%
% DESCRIPTION:
% Function to calculate the rotation matrix to work on the orbital plane. Angles in input in degree
%
% INPUT:
% i     [1x1]   Inclination         [deg]
% OM    [1x1]   RAAN                [deg]
% om    [1x1]   Pericentre anomaly  [deg]
%
% OUTPUT:
% R     [3x3]   Rotation Matrix     [-]
%
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascharetti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

% conversione deg to rad
i = deg2rad(i);
OM = deg2rad(OM);
om = deg2rad(om);

% Calulate the rotation matrix
R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R3_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

R = R3_om*R1_i*R3_OM;

end