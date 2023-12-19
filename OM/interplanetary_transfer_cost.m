function [dv, delta_T] = interplanetary_transfer_cost(t_dep, t_arr, dep_planet, arr_planet, ind)
%
% INPUTS
% dep_planet    Planet where the transfer starts
%
% arr_planet    Planet where the transfer ends
%
% t_dep     time when the transfer starts
% 
% t_arr     time when the transfer ends
%
% ind       ind establishes if the considered transfer requires to
%           calculate deltaV due to two impulses(at the initial and final
%           time) or just one in the case of a gravity assist as starting
%           point or arrival point. 
%           If ind = 0 it will calculate both deltav 
%           If inf = 1 it will calculate only starting impulse 
%           If ind = 2 it will calculate only ending impulse
%           If it there is no input it is assumed 0
%           
%   
% OUTPUT:
% dv        Magnitude of the variation of velocity required for the transfer 
%           to happen
%
% delta_T   Time required for the satellite to travel through the transfer
%           orbit
% the transfer is condidered as interplanetary so the focus will always be
% the Sun


if nargin == 4
    ind = 0;
end


% Physical parameters
ksun = astroConstants(4);  % Sun's gravitational parameter [km^3/s^2]

% t1 = date2mjd2000(t_dep);
% t2 = date2mjd2000(t_arr);
t1 = t_dep;
t2 = t_arr;

[kep_1, ~] = uplanet(t1, dep_planet);
[r_1, v_1] = kep2car(kep_1(1), kep_1(2), kep_1(3), kep_1(4), kep_1(5), kep_1(6), ksun);
[kep_2, ~] = uplanet(t2, arr_planet);
[r_2, v_2] = kep2car(kep_2(1), kep_2(2), kep_2(3), kep_2(4), kep_2(5), kep_2(6), ksun);

delta_T = t2 - t1;
delta_T = delta_T*3600*24;

[~, ~, ~, ~, VI,VF, ~, ~] = lambertMR(r_1, r_2, delta_T, ksun);
if ind == 2
    delta_v1 = 0;
else
    delta_v1 = norm(VI' - v_1);
end

if ind == 1
    delta_v2 = 0;
else
    delta_v2 = norm(v_2 - VF');     
end  
dv = delta_v1 + delta_v2;

r_0 = [r_1, r_2, r_1];
v_0 = [v_1, v_2, VI'];
plot_propagated(r_0, v_0, [0, 0, delta_T], ksun)


