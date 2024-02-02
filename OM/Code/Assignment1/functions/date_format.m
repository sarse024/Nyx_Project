function Date_span=date_format(t_span)
%
%PROTOTYPE:
%Date_span = date_format(t_span)
%
% DESCRIPTION:
% date_format returns a vector in datetime format given a vector of MJD2000
%
% INPUT:
% t_span[nx1]     Time span in mjd2000 
%
% OUTPUT:
% Date_span[nx1]  Time span in datetime format
%
%CONTRIBUTORS:
%Edoardo Nicolucci Balocco
%Mirko Mascaretti
%Valentina D'Annunzio
%Samuele Orsenigo

dates = zeros(6, length(t_span));
for i = 1:length(t_span)
    dates(1:6, i) = mjd20002date(t_span(i));
end

Date_span = datetime(dates(1, :), dates(2, :), dates(3, :));
return
