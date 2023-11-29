function [X,Y,Z] = msph2cart(Br,Bt,Bp)
% Inputs - magnetic field strength (B) in local tangential coordinates
% Br    B in radial direction
% Bt    B in theta direction
% Bp    B in phi direction
% Outputs - magnetic field strength (B) in local tangential coordinates
% X     B in "north" direction
% Y     B in "east" direction
% Z     B in down direction
% Variables
% e     epsilon, the term used to correct for the oblateness of the
%       Earth.  In all programs online, e=0


e = 0.0*pi/180;
X = -Bt*cos(e) - Br*sin(e);
Y = Bp;
Z = Bt*sin(e) - Br*cos(e);