function setGraphicPlot
    set(groot, 'defaultTextInterpreter', 'tex')
    set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
    set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
    set(groot, 'defaultLegendInterpreter', 'tex')
    set(groot, 'defaultAxesTickLabelInterpreter', 'tex')
    set(groot, 'defaultAxesFontWeight', 'bold')
    set(groot, 'defaultFigurePosition', [470, 360, 700, 430])
    set(groot, 'defaultFigureColormap', turbo(256));
    set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
    set(groot, 'defaultSurfaceEdgeAlpha', 0.3);
    set(groot, 'defaultLineLineWidth', 1.6);
    set(groot, 'defaultFigureColor', [1; 1; 1]);
    set(groot, 'defaultAxesColor', 'white');
    set(groot, 'defaultAxesFontSize', 20);
end