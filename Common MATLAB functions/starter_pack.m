function starter_pack
set(groot, 'defaultTextInterpreter', 'latex')
% linee di costruzione su x e y
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultAxesFontWeight', 'bold')
%setting the colors used in graphs
set(groot, 'defaultFigureColormap', turbo(256));
%setting the font type
set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
% face transparency
set(groot, 'defaultSurfaceEdgeAlpha', 0.3);
set(groot, 'defaultLineLineWidth', 1.6);
%color of the background
set(groot, 'defaultFigureColor', [0; 0; 0]);
set(groot, 'defaultAxesColor', 'none');
set(groot, 'defaultAxesFontSize', 20);
% color and transparency of the axes and grid
set(groot, 'defaultAxesXcolor', 'w')
set(groot, 'defaultAxesYcolor', 'w')
set(groot, 'defaultAxesZcolor', 'w')
set(groot, 'defaultAxesGridcolor', 'w')
set(groot, 'defaultAxesGridAlpha', 0.5)
set(groot, 'defaultAxesMinorGridColor', 'w')
set(groot, 'defaultAxesMinorGridAlpha', 0.5)

addpath(genpath('/Users/edoardonicolucci/Documents/GitHub/Nyx_Project'))
end

