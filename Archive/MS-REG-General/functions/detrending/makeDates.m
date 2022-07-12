%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File name: makeDates.m			   %
% Purpose: Constructs dates for CSVs.  %
% Arguments: (start_year,last_year,start_quarter,last_quarter,frequency)%
% Date of Last Edit: 11/23/2018        %
% Original Author: Dawson Miller       %
% Author of Recent Edits: ___Name___   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
	Description:
        Input parameters:
	
		start_year: 	first year of date_series
		last_year: 		last year of date_series
		start_quarter: 	first quarter of date_series
		last_quarter: 	last quarter of date_series
		frequency:		monthly (M), quartery (Q), semi-annual (S), or annual (A)
%}

function [date_series] = makeDates(start_year,last_year,start_quarter,last_quarter,frequency)
	
	time_raw = [start_year:0.25:(last_year + (last_quarter - 1)*0.25)];
	time_raw = time_raw(start_quarter:end);
	qtr_raw  = (time_raw - floor(time_raw))*4 + 1;
	
	switch frequency 
	
		case 'Q'
			year_tmp = strsplit(num2str(floor(time_raw)), ' ');
			qtr_tmp = strsplit(num2str(qtr_raw), ' ');
			freq = cell(1,length(year_tmp));
			freq(:) = {'q'};
			
			date_series = strcat(year_tmp,freq,qtr_tmp);
		
		case 'S'
			year_tmp = time_raw(find(mod(qtr_raw,2) == 0));
			year_tmp = strsplit(num2str(floor(year_tmp)), ' ');
			qtr_tmp = (time_raw(find(mod(qtr_raw,2) == 0)))/2;
			qtr_tmp = strsplit(num2str(floor(qtr_tmp)), ' ');
			freq = cell(1:length(year_tmp));
			freq(:) = {'s'};
			
			date_series = strcat(year_tmp,freq,qtr_tmp);
			
		case 'A'
			year_tmp = time_raw(find(mod(qtr_raw,4) == 0));
			year_tmp = strsplit(num2str(floor(year_tmp)), ' ');
			freq = cell(1:length(year_tmp));
			freq(:) = {'a'};
			
			date_series = strcat(year_tmp,freq);
            
        case 'M'
            time_raw_month = [start_year:(1/12):(last_year + (1/12)*(last_quarter*3 - 1))];
            time_raw_month = time_raw_month((start_quarter*3):end);
            month_raw  = round((time_raw_month - floor(time_raw_month))*12 + 1);
            
            year_tmp = strsplit(num2str(round(floor(time_raw_month))), ' ');
			month_tmp = strsplit(num2str(month_raw), ' ');
			freq = cell(1,length(year_tmp));
			freq(:) = {'m'};
			
			date_series = strcat(year_tmp,freq,month_tmp);
			
		otherwise
			display('Frequency input wrong')
			display('')
			display('M for Monthly \n Q for quarterly \n S for semi-annual \n A for annual')
	end
end