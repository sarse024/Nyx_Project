function [DV, dv_min, index_val, t_dep, deltat1, deltat2] = grid_search(tspan_1, tspan_2, tspan_3)
%
%PROTOTYPE:
%[DV, dv_min, index_val, t_dep, deltat1, deltat2] = grid_search(tspan_1, tspan_2, tspan_3)
%
% DESCRIPTION:
% grid_search performs a search of the minimum cost for a transfer with
% flyby on the second planet. In our specific case the planets are already
% specified(Saturn, Jupiter and Mars)
%
% INPUTS
% tspan_1 [1xl]     Time span of the departure time
%
% tspan_2 [1xm]     Time span of the time of flight for the first Lambert's
%                   arc
%
% tspan_3 [1xn]     Time span of the time of flight for the second Lambert's
%                   arc
%
% OUTPUT:
% DV                Matrix containing all the possible deltav for the
%                   manoeuvre with dimension [l x m x n]
%
% dv_min [1]        Minimum value contained in the matrix DV
%
% index_val [1x3]   Vector containing the indexes of dv_min
%
% t_dep             Optimal time of departure
%
% deltat_1          Optimal time of flight for first Lambert's arc
%
% deltat_1          Optimal time of flight for second Lambert's arc
%
%CONTRIBUTORS:
%Edoardo Nicolucci Balocco
%Mirko Mascaretti
%Valentina D'Annunzio
%Samuele Orsenigo

Saturn = 6; 
Jupiter = 5;
Mars = 4;
mu_Jupiter = astroConstants(15);
% Preallocation of DV and single impulses
DV = ones(length(tspan_1), length(tspan_2), length(tspan_3))*NaN;
DV12 = ones(length(tspan_1), length(tspan_2))*NaN;
DV23 = ones(length(tspan_2), length(tspan_3))*NaN;
r_p_mat = ones(length(tspan_1), length(tspan_2), length(tspan_3))*NaN;
delta_T1_mat = ones(length(tspan_1), length(tspan_2))*NaN;
delta_T2_mat = ones(length(tspan_1), length(tspan_2))*NaN;

for i = 1:length(tspan_1)

    % Saturn position and velocity at chosen time
    flag1 = tspan_1(i);
    [kep_1,ksun] = uplanet(flag1, Saturn);
    [r_1, v_1] = kep2carRAD(kep_1, ksun);

    for j = 1:length(tspan_2)

        % Jupiter position and velocity at chosen time
        flag2 = tspan_2(j) + flag1;
        [kep_2,ksun] = uplanet(flag2, Jupiter);
        [r_2, v_2] = kep2carRAD(kep_2, ksun);
        %first interplanetary arc
        delta_T1_mat(i, j) = tspan_2(j);
        [~, ~, ~, ~, VI,V_minus,~,~] = lambertMR(r_1, r_2, tspan_2(j)*3600*24, ksun);
        DV12(i, j) = norm(VI' - v_1);
        if DV12(i, j)> 4
            DV12(i,j) = NaN;        % to reduce the total number of iterations when DV12 is too high the iteration is skipped
            continue
        end

        for k = 1:length(tspan_3)

            % Mars position and velocity at chosen time
            flag3 = tspan_3(k) + flag2;
            [kep_3,ksun] = uplanet(flag3, Mars);
            [r_3, v_3] = kep2carRAD(kep_3, ksun);
            %second interplanetary arc
            delta_T2_mat(j, k) = tspan_3(k);
            [~, ~, ~, ~, V_plus,VF,~,~] = lambertMR(r_2, r_3, tspan_3(k)*3600*24, ksun);
            DV23(j, k) = norm(v_3 - VF');
            if DV23(j, k)> 8
                DV23(j, k) = NaN;       % to reduce the total number of iterations when DV23 is too high the iteration is skipped
                continue
            end

            % ga assist    
            V_p = v_2;      % velocity of the planet where the gravity assist happens
            v_asym_minus = V_minus' - V_p;
            v_asym_plus = V_plus' - V_p;
            % in order to determine rp we need to solve an implicit function
            delta_true=acos(dot(v_asym_minus,v_asym_plus)/(norm(v_asym_minus)*norm(v_asym_plus)));
            r_per = 71000;
            delta = asin(1/(1+(r_per*norm(v_asym_plus)^2)/mu_Jupiter))+asin(1/(1+(r_per*norm(v_asym_minus)^2)/mu_Jupiter));
            if delta_true>delta
                continue            % if delta_true is bigger than delta than the satellite would crush so the iteration is skipped
            end

            fun_delta=@(r_per) delta_true-(asin(1/(1+(r_per*norm(v_asym_plus)^2)/mu_Jupiter))+asin(1/(1+(r_per*norm(v_asym_minus)^2)/mu_Jupiter)));
            options = optimset('TolFun',1e-15, 'Display', 'off');
            r_p_mat(i, j, k) = fsolve(fun_delta, 71000, options);

            % calculation of delta_vp
            eps_hyp_plus = norm(v_asym_plus)^2/2;
            vp_plus = sqrt(2*(eps_hyp_plus + mu_Jupiter/r_p_mat(i, j, k)));
            eps_hyp_minus = norm(v_asym_minus)^2/2;
            vp_minus = sqrt(2*(eps_hyp_minus + mu_Jupiter/r_p_mat(i, j, k)));
            delta_vp = norm(vp_plus-vp_minus);
            DV(i, j, k) = DV12(i, j) + DV23(j, k) + delta_vp;
        end
    end
end

%% Study case of minimum DV transfer
% estraction of minimum 
% value from DV matrix
[dv_min, index] = min(DV, [], 'all');
[ii, jj, kk] = ind2sub(size(DV), index); 
t_dep = tspan_1(ii);
deltat1 = tspan_2(jj);
deltat2 = tspan_3(kk);
index_val = [ii, jj, kk];