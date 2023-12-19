function plot_propagated(r_0, v_0, t_f, mu)
%
% INPUT
%
% r_0       r_0 is the initial position from which the orbit is propagated,
% [3xn]     in case of multiple orbits each column represent the initial
%           position of one orbit[r1, r2, ..]
%
% v_0       v_0 is the initial velocity from which the orbit is propagated,
% [3xn]     in case of multiple orbits each column represent the initial
%           velocity of one orbit
%
% t_f       final time until which the orbit is propagated, if it is equal
% [n]       to 0 then it will plot an entire period of revolution
% 
% mu        Gravitational constant of the celestial body placed in one of
%           the two foci
%
%
%
% OUTPUT:
% no output just plottinggggg

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Initial condition
for i = 1:size(r_0, 2)
    R_0 = r_0(:, i);
    V_0 = v_0(:, i);
    Y_0 = [ R_0; V_0];
    T_f = t_f(i);
    
    % Eccentricity vector and eccentricity
    e_vect = 1/mu * ((norm(V_0)^2 - mu/norm(R_0)) * R_0 - dot(R_0, V_0) * V_0);
    e = norm(e_vect);
    if e > 1
        %fprintf('the satellite is travelling in an hyperbolic motion\n')
    end


    % Set time span
    if T_f == 0 
        a = 1 / (2/norm(R_0) - dot(V_0,V_0)/mu); % Semi-major axis [km] 
        T = 2*pi*sqrt( a^3/mu ); % Orbital period [1/s]
        t_span = linspace( 0, T, 10000 );
    else 
        t_span = linspace( 0, T_f, 10000 );
    end
     
    % Perform the integration
    [ ~, Y ] = ode113( @(t,y) ode_2bp(t, y, mu), t_span, Y_0, options );
    % Plot the results
    %orbits
    plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]'); title('Two-body problem orbit');
    axis equal;
    grid on; 
    hold on;
    % maneuver points
    plot3(R_0(1), R_0(2), R_0(3), 'o', 'MarkerSize', 10)

end
