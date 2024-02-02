function v_rot=rodrigues_rot(v,delta,u)
%
%PROTOTYPE:
%v_rot=rodrigues_rot(v,delta,u)
%
% DESCRIPTION:
% rodrigues_rot performs a counter-clockwise rotation of an angle delta 
% around a vector u, given an initial vector v %
%
%INPUT:
%v[3x1]             Initial vector
%delta[rad]         Rotation angle
%u[3x1]             Vector around which rotation is performed 
%
%OUTPUT:
%v_rot[3x1]         Rotated vector
%
%CONTRIBUTORS:
%Edoardo Nicolucci Balocco
%Mirko Mascaretti
%Valentina D'Annunzio
%Samuele Orsenigo

cos_delta = cos(delta);
sin_delta = sin(delta);

v_rot = v*cos_delta + cross(u, v)*sin_delta + u*dot(u, v)*(1-cos_delta);
end
