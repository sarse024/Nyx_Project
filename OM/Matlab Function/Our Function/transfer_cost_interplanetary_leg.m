function dv = transfer_cost_interplanetary_leg(planet1, departureTime, planet2, arrivalTime)
% ransfer_cost_interplanetary_leg.m - Calculate the velocity cost for interplanetary leg between two
% planet using ephemerides and lambert.
%
% PROTOTYPE:
% [dv_matrix, x, y, data_min] = Porkchop_plot(planet1, time_departure_window, planet2, time_arrival_window, N)
%
% DESCRIPTION:
% DA FAREEEE!!!
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
%  departureTime[1]             Time for departure from planet1. Time in MJD2000;
%  planet2[1]                   Same of planet1 
%  arrivalTime[1]               Time for arrival to planet2. Time in MJD2000;	
%   N[1]                        Number of point evaluate in window
%                               departure
%   M[1]                        Number of point evaluate in window
%                               arrival
%
%  OUTPUT:
%   dv                          Velocity cost of interplanetary leg [km/s]

% State vector of planet1
[kepE, mu_S] = uplanet(departureTime, planet1);
[R1, V1] = mykep2car(kepE(1), kepE(2), rad2deg(kepE(3)), rad2deg(kepE(4)), rad2deg(kepE(5)), rad2deg(kepE(6)), mu_S);

% State vector of planet2
kepM = uplanet(arrivalTime, planet2);
[R2, V2] = mykep2car(kepM(1), kepM(2), rad2deg(kepM(3)), rad2deg(kepM(4)), rad2deg(kepM(5)), rad2deg(kepM(6)), mu_S);

% time of transfer
dt = arrivalTime - departureTime;
dt = dt*(24*60*60); %conversion from mjd2000 to sec
    
%call lambert
[At,P,E,ERROR,Vt1,Vt2,TPAR,THETA] = lambertMR(R1, R2, dt, mu_S);

% fill the matrix
dv1 = norm(Vt1'-V1);
dv2 = norm(V2-Vt2');
dv = dv1 + dv2;
end