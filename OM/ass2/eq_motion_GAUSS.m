function ds = eq_motion_GAUSS(t, s, ~, parameters)
%our function must include J2 and DRAG resistance in RWS
mu = parameters.mu;

a = s(1);
e = s(2);
i = s(3);
OM = s(4);
om = s(5);
th = s(6);

p = a*(1-e^2);
h = sqrt(p*mu);
r = p/(1+e*cos(th));

% Evaluate the perturbing accelerations
acc_pert_vec = acc_pert_fun_RWS(t, s, parameters); %give acc in RWS form

acc_r = acc_pert_vec(1);
acc_s = acc_pert_vec(2);
acc_w = acc_pert_vec(3);

% Evaluate the equations of motion (Cartesian or Keplerian),
ds = [2*a^2/h*(e*sin(th)*acc_r + p/r*acc_s);
      1/h*(p*sin(th)*acc_r + ((p+r)*cos(th) + r*e)*acc_s);
      (r*cos(th+om))/h*acc_w;
      (r*sin(th+om))/(h*sin(i))*acc_w;
      1/(h*e)*(-p*cos(th)*acc_r + (p+r)*sin(th)*acc_s) - (r*sin(th+om)*cos(i))/(h*sin(i))*acc_w;
      h/r^2 + 1/(e*h)*(p*cos(th)*acc_r - (p+r)*sin(th)*acc_s)];

% function of t, s, acc_pert_vec, and parameters


end