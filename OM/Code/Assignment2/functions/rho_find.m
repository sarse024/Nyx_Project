function rho = rho_find(altitude)
% rho_find.m - Calculation of air density at given altitude
%
% PROTOTYPE:
% rho = rho_find(altitude)
%
% DESCRIPTION:
% Function that calculate air density at input altitude relative to the
% Earth's surface. Implement an exponential law.
%
% INPUT:
% altitude  [1x1] Altitude of the spacecraft relative to the Earth's surface    [Km]
%
% OUTPUT:
% rho       [1x1] Air density at input altitude     [Km/m^3]
%
% AUTHORS:
%   Valentina D'Annunzio      
%   Mirko Mascheretti 
%   Nicolucci Balocco Edoardo
%   Samuele Orsenigo

%data matrix based on r
rho_mat=[0,1.225,7.249;25,3.899e-2,6.349;30,1.774e-2,6.682;40,3.972e-3,7.554;50,1.057e-3,8.382;...
         60,3.206e-4,7.714;70,8.77e-5,6.549;80,1.905e-5,5.799;90,3.396e-6,5.382;100,5.297e-7,5.877;...
         110,9.661e-8,7.263;120,2.438e-8,9.473;130,8.484e-9,12.636;140,3.845e-9,16.149;150,2.07e-9,22.523;...
         180,5.464e-10,29.74;200,2.789e-10,37.105;250,7.248e-11,45.546;300,2.418e-11,53.628;350,9.158e-12,53.628;...
         400,3.725e-12,58.515;450,1.585e-12,60.828;500,6.967e-13,63.822;600,1.454e-13,71.835;700,3.614e-14,88.667;...
         800,1.170e-14,124.64;900,5.245e-15,181.05;1000,3.019e-15,268];

%finding proper parametres into rho_mat
if altitude<25
    h0=rho_mat(1,1);
    rho0=rho_mat(1,2);
    H=rho_mat(1,3);
elseif  altitude>=25 && altitude<30
    h0=rho_mat(2,1);
    rho0=rho_mat(2,2);
    H=rho_mat(2,3);
elseif altitude>=30 && altitude<150
    i=fix(altitude/10);
    h0=rho_mat(i,1);
    rho0=rho_mat(i,2);
    H=rho_mat(i,3);
elseif altitude>=150 && altitude<180
    h0=rho_mat(15,1);
    rho0=rho_mat(15,2);
    H=rho_mat(15,3);
elseif altitude>=180 && altitude<200
    h0=rho_mat(16,1);
    rho0=rho_mat(16,2);
    H=rho_mat(16,3);
elseif altitude>=200 && altitude<500
    i=fix(altitude/50)+13;
    h0=rho_mat(i,1);
    rho0=rho_mat(i,2);
    H=rho_mat(i,3);
elseif altitude>=500 && altitude<1000
    i=fix(altitude/100)+18;
    h0=rho_mat(i,1);
    rho0=rho_mat(i,2);
    H=rho_mat(i,3);
else 
    h0=rho_mat(28,1);
    rho0=rho_mat(28,2);
    H=rho_mat(28,3);
end

%rho formula
rho_fun=@(h) rho0*exp(-(h-h0)/H);

%finding the correct rho
rho=rho_fun(altitude);

return






