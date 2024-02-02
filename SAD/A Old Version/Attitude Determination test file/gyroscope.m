close all
clear
clc

t0 = 0;
tf = 200;
Ix = 0.0700; %[kg*m^2]
Iy = 0.0550; %[kg*m^2]
Iz = 0.025; %[kg*m^2]
Ir = 0; %[kg*m^2]
I = [Ix 0 0; 0 Iy 0; 0 0 Iz];
I_inv = inv(I);
omega_0 = [0.45 0.52 0.55]'; %[rad/s]
omega_r0 = 0;
rot = [0 0 Ir*omega_r0]';
M = 0;
q_0 = [0 0 0 1]';
bias = 0.05; %[deg/h]
bias = bias*(pi/180)*(1/3600); %[rad/s]
bias = bias*[1 1 1]';
ts = 0.1;
ARW = 0.007; %[deg/sqrt(h)]
ARW = ARW*(pi/180)*(1/60)*(1/sqrt(ts)); %[rad/sqrt(h)]
ARW = ARW*[1 1 1]';