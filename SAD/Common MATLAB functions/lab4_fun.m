function zero = lab4_fun(r_p, v_asym_minus, v_asym_plus, mu)
e_minus =  1 + norm(r_p)*(norm(v_asym_minus)^2)/mu;
delta_minus = 2*asin(1/e_minus);

e_plus = 1 + norm(r_p)*(norm(v_asym_plus)^2)/mu;
delta_plus = 2*asin(1/e_plus);

delta = acos(dot(v_asym_minus, v_asym_plus)/(norm(v_asym_minus*norm(v_asym_plus))));
zero = delta_minus/2 + delta_plus/2-delta;

end