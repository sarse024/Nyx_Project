function Terra3d(Rt)
% plotta la terra in 3d con il suo raggio nell origine degli assi
% basta richiamare la funzione

if nargin==0
    Rt=6378; % Km
end

load('topo.mat','topo','topomap1');
% whos topo topomap1
contour(0:359,-89:90,topo,[0 0],'b')
axis equal
% box on
% set(gca,'XLim',[0 360],'YLim',[-90 90], ...
%     'XTick',[0 60 120 180 240 300 360], ...
%     'Ytick',[-90 -60 -30 0 30 60 90]);
hold on
image([0 360],[-90 90],topo,'CDataMapping', 'scaled');
colormap(topomap1)
[x,y,z] = sphere(100);
% cla reset
% axis square off
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
surface(x*Rt,y*Rt,z*Rt,props);
light('position',[1 1 1]);
light('position',[-1.5 0.5 -0.5], 'color', [.6 .2 .2]);

% quiver3(0,0,0,0,0,5*Rt,'--k')