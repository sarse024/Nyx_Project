function [r,v]=kep2car(a,e,i,OM,om,th,mu)
%angles in radiant

%semilato retto
p=a*(1-e^2);    

%modulo vettore posizione
r_mod=p/(1+e*cos(th));      

%vettori posizione e velocità nel sistema PF
r_PF=r_mod*[cos(th),sin(th),0];            
v_PF=sqrt(mu/p)*[-sin(th),e+cos(th),0];

%matrice di rotazione da PF a ECI
R_trasp=[cos(om)*cos(OM)-sin(om)*sin(OM)*cos(i), ...  
    -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM), sin(i)*sin(OM); sin(OM)*cos(om)+...
    sin(om)*cos(OM)*cos(i),-sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM), -sin(i)*cos(OM);...
    sin(om)*sin(i), cos(om)*sin(i), cos(i)];

%vettori posizione e velocità nel sistema ECI
r=R_trasp*r_PF';         
v=R_trasp*v_PF';


