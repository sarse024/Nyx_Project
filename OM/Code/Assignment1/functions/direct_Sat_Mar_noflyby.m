function dv_min_noflyby=direct_Sat_Mar_noflyby(~)
%direct_Sat_Mar_noflyby searches the best deltaV in the case of direct interplanetary trajectory fom Saturn to Mars
%
%PROTOTYPE:
%dv_min_noflyby=direct_Sat_Mar_noflyby(~)
%
%INPUT:
%none
%
%OUTPUT:
%none
%
%CONTRIBUTORS:
%Edoardo Nicolucci Balocco
%Mirko Mascaretti
%Valentina D'Annunzio
%Samuele Orsenigo


% t_span definition
imax=3000;
jmax=1500;
% we know, after we ran the code, the time to arrive to Mars with the minimum Dv is about 7 years,
% so we can eliminate at least the last 6 years of possible departure from Saturn to
% increase the velocity of the code itself, avoiding useless departure
% windows
t_span1=linspace(date2mjd2000([2028,1,1,0,0,0]),date2mjd2000([2052,1,1,0,0,0]),imax); 

% pre-allocation
r_sat=zeros(3,length(t_span1));
v_sat=zeros(3,length(t_span1));
r_mar=zeros(3,length(t_span1),jmax);
v_mar=zeros(3,length(t_span1),jmax);
kep_sat=zeros(6,length(t_span1));
kep_mar=zeros(6,length(t_span1),jmax);
Dv_mar_jup=inf*ones(imax,jmax);
t_mat1=zeros(imax,jmax);
v1_tr1=zeros(3,length(t_span1),jmax);
v2_tr1=zeros(3,length(t_span1),jmax);


for i=1:length(t_span1)
    for j=1:jmax

        % using ephemerides 
        [kep_sat(:,i),ksun]=uplanet(t_span1(i),6);    
        %Mars position every 12h to be precise
        [kep_mar(:,i,j),~]=uplanet(t_span1(i)+5*365+j/2,4); 
        [r_sat(:,i),v_sat(:,i)]=kep2car(kep_sat(1,i),kep_sat(2,i),kep_sat(3,i),kep_sat(4,i),kep_sat(5,i),kep_sat(6,i),ksun);
        [r_mar(:,i,j),v_mar(:,i,j)]=kep2car(kep_mar(1,i,j),kep_mar(2,i,j),kep_mar(3,i,j),kep_mar(4,i,j),kep_mar(5,i,j),kep_mar(6,i,j),ksun);

        % time for each transfer in seconds
        t_mat1(i,j)=(5*365+j/2)*3600*24;    

        % Lambert function 
        [~,~,~,~,v1_tr1(:,i,j),v2_tr1(:,i,j),~,~]=lambertMR(r_sat(:,i),r_mar(:,i,j),t_mat1(i,j),ksun,0,0,0,0);
        Dv_mar_jup(i,j)=norm(v_sat(:,i)-v1_tr1(:,i,j))+norm(v_mar(:,i,j)-v2_tr1(:,i,j));

    end
end

% minimum Dv and its index calculation
dv_min_noflyby = min(min(Dv_mar_jup));
[~,loc] = min(Dv_mar_jup(:));
[~,~] = ind2sub(size(Dv_mar_jup),loc);

