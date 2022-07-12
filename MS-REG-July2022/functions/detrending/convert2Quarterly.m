%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File name: convert2Quarterly.m       %
% Purpose: Constructs dates for CSVs.  %
% Date of Last Edit: 11/23/2018        %
% Original Author: Dawson Miller       %
% Author of Recent Edits: ___Name___   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data_avg] = convert2Quarterly(data_input,date_start,date_end)
%{
    Documentation:
    
        data_input: 
                n-by-1 vector of numerical values

        date_start:
            character array in form of YYYYFMM or YYYYFM, where F is A, Q, or M
            for annual, quarterly, or monthly conversions, respectively

        date_end:
            character array in same format as date_start

        data_avg:
            simple average of data within Year-Quarter groups
%}

dates = [str2num(date_start(1:4)):(1/12):(str2num(date_end(1:4))+1)]';
dates = dates(str2num(date_start(6:end)):end-1);
dates = dates(1:(end - (12 - str2num(date_end(6:end)))));
months = round(12*(dates - floor(dates))) + 1;
quarters = ceil(months/3);
groups = findgroups(floor(dates),quarters);
dates = [floor(dates),months,quarters,groups]; 

current_group = groups(1);
current_quarter = 1;
months_in_group = zeros(length(unique(groups)),1);
months_in_group(1) = 1;
data = zeros(length(unique(groups)),1);
data(1) = data_input(1);

for ii = 2:length(groups)
   if current_group == groups(ii)
      %continue summation
      data(current_quarter) = sum(data(current_quarter) + data_input(ii));
      months_in_group(current_quarter) = months_in_group(current_quarter) + 1;
   else
      %reset current groups
      current_quarter = current_quarter + 1;
      months_in_group(current_quarter) = 0;
      current_group = groups(ii);
      
      data(current_quarter) = sum(data(current_quarter) + data_input(ii));
      months_in_group(current_quarter) = months_in_group(current_quarter) + 1;
   end
end

data_avg = data./months_in_group;

end