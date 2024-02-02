function ds = eq_motion_CAR( t, s, acc_per_fun_CAR, parameters)
%our function must include J2 and DRAG resistance

parameters.kep = car2kepRAD(s(1:3), s(4:6), parameters.mu);
mu = parameters.mu;
% Evaluate the perturbing accelerations
acc_pert_vec = acc_pert_fun_CAR( t, s, parameters);

% Evaluate the equations of motion (Cartesian or Keplerian),
ds = [s(4:6);-mu/norm(s(1:3))^3*s(1:3) + acc_pert_vec];

% calcualte altitude 
alt = norm(s(1:3)) - parameters.rE;
end