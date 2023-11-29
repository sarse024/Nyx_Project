function [E,t_span,i]=kepler_equation(a,e,delta_t,mu,t0,Npoints,tol)

%return a vector of E (in degrees) calculated in N points 
%if in input we have E0 instead of t0, transform it: t0=1/n*(E0-e*sin(E0))
%t_span vector of time in which we evaluate E, used to plot the results 

if nargin<7
    tol=1e-6;   %default tolerance
end

if nargin<6
    Npoints=10;   %default number of points
end

% time data and mean velocity
tf=delta_t+t0;
t=linspace(t0,tf,Npoints);
t_span=linspace(t0,tf,Npoints);
T=2*pi*sqrt(a^3/mu);
n=sqrt(mu/a^3);

%verify if more than 1 orbit
for i=1:length(t)
    num_orb=0;
    if t(i)>T 
        num_orb=fix(t(i)/T);
        t(i)=mod(t(i),T);
    end

    %initial guess to have fzero faster
    Eguess=n*t(i)+(e*sin(n*t(i)))/(1-sin(n*t(i)+e)+sin(n*t(i)));
    
    E(i)=fzero(@(E) E-e*sin(E)-n*t(i),Eguess,optimset('TolX',tol));
    E(i)=E(i)/2/pi*360+num_orb*360;
end



