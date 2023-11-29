function [delta_t,Ei,Ef]=timeOfFlight(a,e,thi,thf)

%definisco mu
mu=398600.433;                          

 %definisco e calcolo E
Ei=2*atan(sqrt((1-e)/(1+e))*tan(thi/2));           
Ef=2*atan(sqrt((1-e)/(1+e))*tan(thf/2));

%porto Ei,Ef>0
if Ei<0
    Ei=Ei+2*pi;
end
if Ef<0
    Ef=Ef+2*pi;
end

% if Ei<0                 %calcolo Ei>0
%     if cos(thi)>=0
%         Ei=Ei+2*pi;
%     else
%         Ei=Ei+pi;
%     end
% end
% 
% if Ef<0                 %calcolo Ef>0
%     if cos(thf)>=0
%         Ef=Ef+2*pi;
%     else
%         Ef=Ef+pi;
%     end
% end

%calcolo M e delta_M
Mi=Ei-e*sin(Ei);                
Mf=Ef-e*sin(Ef);
delta_M=Mf-Mi;

%calcolo delta_t
if delta_M>=0                       
    delta_t=delta_M/sqrt(mu/a^3);  
else
    delta_t=delta_M/sqrt(mu/a^3) + 2*pi*sqrt(a^3/mu); 
end




