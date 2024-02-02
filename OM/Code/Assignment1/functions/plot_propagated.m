function plot_propagated(r_0, v_0, t_f, mu, line)
%
%PROTOTYPE:
%plot_propagated(r_0, v_0, t_f, mu)
%
% DESCRIPTION:
% plot_propagated does the plot of an orbit starting from a given position
% vector and a given velocity vector, integrating the velocity for a t_f
% time
%
% INPUT:
%
% r_0[3xn]      r_0 is the initial position from which the orbit is propagated,
%               in case of multiple orbits each column represent the initial
%               position of one orbit[r1, r2, ..]
%
% v_0[3xn]       v_0 is the initial velocity from which the orbit is propagated,
%               in case of multiple orbits each column represent the initial
%               velocity of one orbit
%
% t_f[n]        final time until which the orbit is propagated, if it is equal
%               to 0 then it will plot an entire period of revolution
% 
% mu            Gravitational constant of the celestial body placed in one of
%               the two focii
%
% line          If line = 0 than the plot will have a normal line, if line= 1 
%               than the plotline will be dashed. If there is no input line
%               is assumed to be equal to 0
%
% OUTPUT:
% none

if nargin == 4
    line = 0;
end

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
    
    % Plot the results:
    %orbits
    if line == 0
        plot3( Y(:,1), Y(:,2), Y(:,3), '-', 'LineWidth', 1);
    else
        plot3( Y(:,1), Y(:,2), Y(:,3), '--', 'LineWidth', 1);
    end

    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    axis equal;
    grid on; 
    hold on;
    % maneuver points
    plot3(R_0(1), R_0(2), R_0(3), 'o', 'MarkerSize', 5, 'HandleVisibility','off')

end
