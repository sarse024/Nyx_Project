function [dv_matrix, x, y, data_min] = Porkchop_plot(planet1, time_departure_window, planet2, time_arrival_window, N, M)
% Porkchop_plot.m - Calculate the best interplanetary leg between two
% planet using ephemerides, lambert and optimisation tool. 
% Gives back the data to make porkchop contour plot.
%
% PROTOTYPE:
% [dv_matrix, x, y, data_min] = Porkchop_plot(planet1, time_departure_window, planet2, time_arrival_window, N)
%
% DESCRIPTION:
% DA FAREEEE!!! -> capire bene se ha senso restituire i dati in data_min
%
%  INPUT :
%  planet1[1]                  Integer number identifying the celestial body (< 11)
%                                    1:   Mercury
%                                    2:   Venus
%                                    3:   Earth
%                                    4:   Mars
%                                    5:   Jupiter
%                                    6:   Saturn
%                                    7:   Uranus
%                                    8:   Neptune
%                                    9:   Pluto
%                                    10:  Sun
%  time_departure_window[2]     Vector of initial and final time for
%                               departure from planet1. Time in MJD2000;
%  planet2[1]                   Same of planet1 
%  time_arrival_window[2]       Vector of initial and final time for
%                               arrival to planet2. Time in MJD2000;	
%   N[1]                        Number of point evaluate in window
%                               departure
%   M[1]                        Number of point evaluate in window
%                               arrival
%
%  OUTPUT:
%   dv_matrix[NxM]              Matrix of velocity cost evaluete in for
%                               each departure and arrival time
%   x                           Mesh grid axial for conter plot
%   y                           Mesh grid vertical for conter plot
%   data_min                    Structure that contain data of optimal orbit
%       -> cost[1]              Min velocity transfer cost [Km/s]
%       -> t_dep[1]             Time of departure [date]
%       -> t_arr[1]             Time of arrival [date]
%	    -> vel_SC[3x4]    	    Velocity of S/C [v_p1', Vt1', Vt2', v_p2']
%       -> dt[1]                Time of transfer [s]
%       -> lamb_error[1]        Lambert algorithm error


% departure window Earth
t1_dep = time_departure_window(1);
t2_dep = time_departure_window(2);
windows_dep = linspace(t1_dep, t2_dep, N);

% departure window Mars
t1_arr = time_arrival_window(1);
t2_arr = time_arrival_window(2);
windows_arr = linspace(t1_arr, t2_arr, M);

% create output for contur plot
[x,y] = meshgrid(windows_dep, windows_arr);
dv_matrix = zeros(length(windows_arr), length(windows_dep));

for i = 1:length(windows_dep) %columns -> departure from planet1
    for j = 1:length(windows_arr) %rows -> arrive at Mars
        
        dv_matrix(j,i) = transfer_cost_interplanetary_leg(planet1, windows_dep(i), planet2, windows_arr(j));
    end
end

% search minimum in dv_matrix
minimun = min(min(dv_matrix));
[i,j]=find(dv_matrix==minimun);

% ---- Refine the solution using Matlab-s fminunc ----

% define function to minimize
fun = @(t) transfer_cost_interplanetary_leg(planet1,t(1),planet2,t(2));

% initial time guess condition for fminunc using the minimum of contour
% plot
x0 = [windows_dep(i), windows_arr(j)];

% options setup for fminunc
options = optimoptions('fminunc','Algorithm','quasi-newton');
%options.Display = 'iter'; %enable to see iteration

% call fminunc
[t_sol, fval, exitflag] = fminunc(fun,x0,options);

if(exitflag == 0) %check condition
    error('\n\nFminuc does not converge');
else
    % conversion time from mjd2000 t0 date for simpler visualisation
    t_min_departure = mjd20002date(t_sol(1));
    t_min_arrival = mjd20002date(t_sol(2));
    
    % calculation of transfer orbit optimal 
    [kep_p1_dep, mu_S] = uplanet(t_sol(1), planet1); %radius of departure at p1
    kep_p2_dep = uplanet(t_sol(1), planet2); %radius of departure at p2
    kep_p2_arr = uplanet(t_sol(2), planet2); %radius of arrival p2
    
    [r_p1_dep, v_p1_dep] = kep2carRAD(kep_p1_dep, mu_S);
    [r_p2_dep, v_p2_dep] = kep2carRAD(kep_p2_dep, mu_S);
    [r_p2_arr, v_p2_arr] = kep2carRAD(kep_p2_arr, mu_S);
    
    % time of transfer
    dt = t_sol(2) - t_sol(1);
    dt = dt*(24*60*60); %conversion from mjd2000 to sec
            
    %call lambert
    [~,~,~,ERROR,Vt1, Vt2] = lambertMR(r_p1_dep, r_p2_arr, dt, mu_S);
    
    % print output
    fprintf('\n\n\n -------- RESULT OF MISSION DESIGN ----------\n')
    fprintf('Departure from Earth -> %4.0f/%2.0f/%2.0f  %2.0f h:%2.0f m:%2.2f s\n', t_min_departure(1), t_min_departure(2), t_min_departure(3), t_min_departure(4), t_min_departure(5), t_min_departure(6))
    fprintf('Arrival to Mars -> %4.0f/%2.0f/%2.0f  %2.0f h:%2.0f m:%2.2f s\n', t_min_arrival(1), t_min_arrival(2), t_min_arrival(3), t_min_arrival(4), t_min_arrival(5), t_min_arrival(6))
    fprintf('Minimum cost of transfer: %2.5f [km/s]\n', fval)
    
    % Save data of min in a structure
    data_min.cost = fval;
    data_min.t_dep = t_sol(1);
    data_min.t_arr = t_sol(2);
    data_min.vel_SC = [v_p1_dep, Vt1', Vt2', v_p2_arr];
    data_min.dt = dt;
    data_min.lamb_error = ERROR;
end
end