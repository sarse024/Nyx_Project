clear
clc

% Reaction wheels - pyramid configuration
a = 1/sqrt(3);
A = [-a a a -a; -a -a a a; a a a a];
b = sqrt(3)/4;
A_ = [-b -b b; b -b b; b b b; -b b b];
h_r_0 = [0 0 0 0]';
M_c = [1 0 0]';

%%
% Satellite data
Ix = 0.0700; %[kg*m^2]
Iy = 0.0550; %[kg*m^2]
Iz = 0.025; %[kg*m^2]
Ir = 0; %[kg*m^2]
I = [Ix 0 0; 0 Iy 0; 0 0 Iz];
I_inv = inv(I);
omega_0 = [0 0 0]'; %[rad/s]
omega_r0 = 0;
M = 0;