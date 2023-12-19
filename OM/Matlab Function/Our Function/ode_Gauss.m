function dkep = ode_Gauss(~,kep,acc,rE,p,h,j2,mu,y,wE)


% Keplerian elements
da = 2*(kep(1)^2)/h*(kep(2)*sin(kep(6))*acc(1)+p/norm(y(1:3))*acc(2));
de = 1/h*(p*sin(kep(6))*acc(1)+((p+norm(y(1:3)))*cos(kep(6))*norm(y(1:3))*kep(2))*acc(2));
di = (norm(y(1:3))*cos(kep(6)+kep(5)))/h*acc(3);
dOM = (norm(y(1:3))*sin(kep(6)+kep(5)))/(h*sin(kep(3)))*acc(3);
dom = 1/(h*kep(2))*(-p*cos(kep(6))*acc(1)+(p+norm(y(1:3)))*sin(kep(6))*acc(2))-(norm(y(1:3))*sin(kep(6)+kep(5))+cos(kep(3))/(h*sin(kep(3)))*acc(3));
dth = h/norm(y(1:3))^2+1/(kep(2)*h)*(p*cos(kep(6))*acc(1)-(p+norm(y(1:3))*sin(kep(6))*acc(2)));

dkep = [da; de; di; dOM; dom; dth];
end