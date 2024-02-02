function [Bx,By,Bz] = msph2inert(Br,Bt,Bp,LST,lat)
% Inputs
% Br    B in radial direction
% Bt    B in theta direction
% Bp    B in phi direction
% LST   Local sidereal time of location (in degrees)(right ascension)
% lat   Latitude measured positive north from equator (in degrees), it is
%       the same as the local declination
%
% Outputs
% Bx    B in x-direction  | Magnetic field strength (B)
% By    B in y-direction  | in geocentric inertial coordinates
% Bz    B in z-direction  |

% Angle conversion to radians
lat=lat*pi/180;  LST=LST*pi/180;
lat = pi/2-lat;

% Coordinate transformation
Bx = (Br*sin(lat)+Bt*cos(lat))*cos(LST) - Bp*sin(LST);
By = (Br*sin(lat)+Bt*cos(lat))*sin(LST) + Bp*cos(LST);
Bz = (Br*cos(lat)-Bt*sin(lat));


