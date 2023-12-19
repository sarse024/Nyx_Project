function kep =car2kepRAD(r,v,mu)
% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a, e, i, OM, om, th] = car2kep(r, v, mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in
% radians.
%
% INPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% kep[6x1] 
%       a [1x1] Semi-major axis [km]
%       e [1x1] Eccentricity [-]
%       i [1x1] Inclination [rad]
%       OM [1x1] RAAN [rad]
%       om [1x1] Pericentre anomaly [rad]
%       th [1x1] True anomaly [rad]

% toll
toll = 1e-6;

r_norm=norm(r); %norm position vector
v_norm=norm(v); %norm velocity vector

%angular momentum and his vector
h_vect=cross(r,v);
h_norm=norm(h_vect);

%inclination
i=acos(h_vect(end)/h_norm);

%eccentricity vector and eccentricity
e_vect=(1/mu)*((v_norm^2-mu/r_norm)*r-dot(r,v)*v);
e=norm(e_vect);

%specific mechanic energy and major semiaxis
eps=0.5*v_norm^2-mu/r_norm;
a=-mu/(2*eps);

%nodes line and his module
k=[0 0 1];
n_vect = cross(k,h_vect);
n_norm = norm(n_vect);

%RAAN (Ascensione retta del nodo ascendente)
if abs(i) <= toll
    OM = 0; %avoid singularity
    n_vect = [1, 0, 0]';
    n_norm = 1;
elseif n_vect(2)>=0
    OM = acos(n_vect(1)/n_norm);
else
    OM = 2*pi - acos(n_vect(1)/n_norm);
end
%pericentre anomaly
e_norm = norm(e_vect);

if(abs(e) <= toll)
    om = 0;
    e = 0;
    e_vect = n_vect;
elseif e_vect(end)>=0
    om=acos(dot(n_vect,e_vect)/(n_norm*e_norm));
else
    om=2*pi-acos(dot(n_vect,e_vect)/(n_norm*e_norm));
end

%radius velocity
vr=dot(r,v)/r_norm;

%true anomaly

if vr>=0
    th = acos(dot(e_vect,r)/(e_norm*r_norm));
else
    th=2*pi-acos(dot(e_vect,r)/(e_norm*r_norm));
end

% create keplerian array
kep = [a,e,i,OM,om,th];
end