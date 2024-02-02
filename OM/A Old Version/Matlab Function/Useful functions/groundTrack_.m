function [alpha, delta, lon, lat] = groundTrack_(Y,theta_G_0,omega_E,tspan)

theta_G = zeros(size(Y,1),1);
for i = 1:size(Y,1)
    theta_G(i) = theta_G_0 + omega_E*(tspan(i));
end

%%
% Declination
delta = zeros(size(Y,1),1);
for i = 1:size(Y,1)
    delta(i) = asin(Y(i,3)/vecnorm(Y(i,1:3),2,2));
end

%%
% Right ascention
alpha = zeros(size(Y,1),1);
for i = 1:size(Y,1)
    if Y(i,2)/vecnorm(Y(i,1:3),2,2) > 0
        alpha(i) = acos(Y(i,1)/(vecnorm(Y(i,1:3),2,2)*cos(delta(i))));
    else
        alpha(i) = 2*pi-acos(Y(i,1)/(vecnorm(Y(i,1:3),2,2)*cos(delta(i))));
    end
end
%alpha(i) = atan2(Y(i,2),Y(i,1));
%%
% Longitude
lon = alpha-theta_G;

%% 
% Latitude 
lat = delta;

%%
% Conversion to degrees
delta = rad2deg(delta);
alpha = rad2deg(alpha);
lon = rad2deg(lon);
lat = rad2deg(lat);
lon = wrapTo180(lon);
end

