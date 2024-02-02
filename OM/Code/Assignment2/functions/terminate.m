function [lookfor, stop, direction] = terminate(~,s)
% terminate.m - A function useful for interrupting ODE integration. 
% It stops the calculation after the spacecraft reaches an altitude of 100 km.
%
% PROTOTYPE:
% [lookfor, stop, direction] = terminate(t,s)
%
% INPUT:
% t [1x1] Time                        [s]
% s [6x1] State vector:
%    - r     [3x1] Position vector    [km]
%    - v     [3x1] Velocity vector    [km/s]
%
% OUTPUT:
% lookfor   [1x1]  % The value that we want to be zero
% stop      [1x1]  % Stop the integration when lookfor = 0
% direction [1x1]  % Direction to find the lookfor = 0 condition (above -1,
% 1 under)
%
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascheretti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

lookfor = norm(s(1:3)) - 100 - 6371; % = 0 when altitude = 100 km
stop = 1; % 1 means terminate at lookfor = 0; Otherwise 0
direction = -1; % -1 means zero crossing is from above

end