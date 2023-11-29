function [X,Y,Z]=plotOrbit(kepEl,mu,delta_th,step_th)

% kepEl è un vettore contenente i 6 parametri orbitali. In ordine sono:
a=kepEl(1);
e=kepEl(2);
i=kepEl(3);         
OM=kepEl(4);
om=kepEl(5);
th=kepEl(6);

% calcolo il semilato retto
p=a*(1-e^2);      

% inizializzo i vettori contenenti le coordinate dei punti dell'orbita
X=zeros(1,fix(delta_th/step_th));
Y=zeros(1,fix(delta_th/step_th));
Z=zeros(1,fix(delta_th/step_th));

% definisco n come l'anomalia vera durante l'avanzamento e 
% la inizializzo con il valore di th
n=th;               

% definisco j come contatore dei vari punti dell'orbita all'interno del
% ciclo while
j=0;            

% definisco la matrice di rotazione dal sistema perifocale a quello
% geocentrico equatoriale inerziale
R_trasp=[cos(om)*cos(OM)-sin(om)*sin(OM)*cos(i), ...  
    -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM), sin(i)*sin(OM); sin(OM)*cos(om)+...
    sin(om)*cos(OM)*cos(i),-sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM), -sin(i)*cos(OM);...
    sin(om)*sin(i), cos(om)*sin(i), cos(i)];

% ciclo while per l'avanzamento sull'orbita
while (n-th<delta_th)  

% calcolo il vettor posizione nei due sistemi
r_mod=p/(1+e*cos(n));
r_PF=r_mod*[cos(n),sin(n),0];

r=R_trasp*r_PF';

% aumento il contatore per poter calcolare poi le coordinate dei punti
% all'iterazione successiva
j=j+1;

% calcolo le coordinate dei punti nello spazio ad ogni iterazione
X(j)=r(1);                     
Y(j)=r(2);
Z(j)=r(3);

% incremento l'anomalia per passare al punto successivo
n=n+step_th;
end                     

%plot dell'orbita
figure();
Terra3d;
hold on;
xlabel('X[km]');
ylabel('Y[km]');
zlabel('Z[km]');
axis equal;
grid on;
plot3(X,Y,Z,LineWidth=3);

return




