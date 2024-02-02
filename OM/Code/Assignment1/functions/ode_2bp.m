function dy = ode_2bp( ~, y, mu )
%ode_2bp ODE system for the two-body problem (Keplerian motion) %
% PROTOTYPE
% dy=ode_2bp(t,y,mu)
%
% INPUT:
% t[1]
% y[6x1]
% mu[1]
%
% OUTPUT:
% dy[6x1]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez %
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
       % Position and velocity
       r = y(1:3);
       v = y(4:6);
       % Distance from the primary
       rnorm = norm(r);
       % Set the derivatives of the state
dy=[ v ; (-mu/rnorm^3)*r ];
end