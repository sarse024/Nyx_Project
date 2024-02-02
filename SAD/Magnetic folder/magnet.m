function [Br,Bt,Bp] = magnet(r,theta,phi,days, N, A1, A2)
% Inputs
% r         Geocentric radius
% theta     Latitude measured in degrees positive from equator  
% phi       Longitude measured in degrees positive east from Greenwich
% days      Decimal days since January 1, 2020
% N         Order of expansion of the magnetic field model
%
% Outputs - magnetic field strength in local tangential coordinates
% Br    B in radial direction
% Bt    B in theta direction
% Bp    B in phi direction
%
% Checks to see if located at either pole to avoid singularities
if (theta>-0.00000001 && theta<0.00000001)
    theta=0.00000001;
elseif(theta<180.00000001 && theta>179.99999999)
    theta=179.99999999;
end
% The angles must be converted from degrees into radians
theta=(90-theta)*pi/180;
phi = phi*pi/180;
a=6371.2;  % Reference radius used in IGRF
days = (days - 7304.5);

if days < 0
    error('the time chosen must be after January 1, 2020')
end

% This section of the code simply reads in the g and h Schmidt
% quasi-normalized coefficients
gn = A1(:, 1); gm =A1(:, 2);
gvali = A1(:, 3); gsvi = A1(:, 4);

hn = A2(:, 1); hm =A2(:, 2);
hvali = A2(:, 3); hsvi = A2(:, 4);


g=zeros(N,N+1);
h=zeros(N,N+1);

for x=1:sum(2:N+1)
    g(gn(x),gm(x)+1) = gvali(x) + gsvi(x)*days/365;
    h(hn(x),hm(x)+1) = hvali(x) + hsvi(x)*days/365;
end
% Initialize each of the variables
% Br        B in the radial driection
% Bt        B in the theta direction
% Bp        B in the phi direction
% P         The associated Legendre polynomial evaluated at cos(theta)
%           The nomenclature for the recursive values generally follows
%           the form P10 = P(n-1,m-0)
% dP        The partial derivative of P with respect to theta
Br=0;  Bt=0;  Bp=0;
P11=1;  P10=P11;
dP11=0;  dP10=dP11;
P20 = 0; dP20 = 0;
for m=0:N
    for n=1:N
        if m<=n
            % Calculate Legendre polynomials and derivatives recursively
            if n==m
                P2 = sin(theta)*P11;
                dP2 = sin(theta)*dP11 + cos(theta)*P11;
                P11=P2;   P10=P11;   P20=0;
                dP11=dP2; dP10=dP11; dP20=0;
            elseif n==1
                P2 = cos(theta)*P10;
                dP2 = cos(theta)*dP10 - sin(theta)*P10;
                P20=P10;   P10=P2;
                dP20=dP10; dP10=dP2;
            else
                K = ((n-1)^2-m^2)/((2*n-1)*(2*n-3));
                P2 = cos(theta)*P10 - K*P20;
                dP2 = cos(theta)*dP10 - sin(theta)*P10 - K*dP20;
                P20=P10;   P10=P2;
                dP20=dP10; dP10=dP2;
            end
            % Calculate Br, Bt, and Bp
            Br = Br + (a/r)^(n+2)*(n+1)*...
                ((g(n,m+1)*cos(m*phi) + h(n,m+1)*sin(m*phi))*P2);
            Bt = Bt + (a/r)^(n+2)*...
                ((g(n,m+1)*cos(m*phi) + h(n,m+1)*sin(m*phi))*dP2);
            Bp = Bp + (a/r)^(n+2)*...
                (m*(-g(n,m+1)*sin(m*phi) + h(n,m+1)*cos(m*phi))* P2);
        end         
    end
end

Bt = -Bt;
Bp = -Bp/sin(theta);
