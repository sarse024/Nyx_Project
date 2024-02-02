function [a,e,i,OM,om,th]=car2kep(r,v,mu) %ciao

%definisco h momento della quantità di moto specifica
h=cross(r,v);     

%definisco i moduli dei tre parametri cartesiani
r_mod=norm(r);     
v_mod=norm(v);
h_mod=norm(h);

%inclinazione dell'orbita
i=acos(h(3)/h_mod);      

%eccentricità
e_vect=1/mu*((v_mod^2 - mu/r_mod)*r - dot(r,v)*v);    
e=norm(e_vect);

%energia specifica
E=1/2*v_mod^2 - mu/r_mod;   

%semiasse maggiore
a=-mu/(2*E);       

%asse nodale
N=cross([0 0 1],h);   
N_mod=norm(N);

%ascensione retta del nodo ascendente
if i==0
    OM=0;       %arbitrary OM=0 if inclination is 0
else
if N(2)>=0               
    OM=acos(N(1)/N_mod);
else
    OM=2*pi - acos(N(1)/N_mod);
end
end

%anomalia del pericentro
if e==0
    om=0;       %arbitrary om=0 if eccentricity is 0
else
if e_vect(3)>=0                             
    om=acos(dot(N,e_vect)/(N_mod*e));
else 
    om=2*pi - acos(dot(N,e_vect)/(N_mod*e));
end
end

%velocità radiale
v_rad=dot(r,v)/r_mod;     

%anomalia vera
if v_rad>=0                           
    th=acos(dot(e_vect,r)/(e*r_mod));
else
    th=2*pi - acos(dot(e_vect,r)/(e*r_mod));
end




