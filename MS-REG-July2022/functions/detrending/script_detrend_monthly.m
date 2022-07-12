% This script detrends monthly series 
%
% IMPORTANT: No NaNs can exist in middle of series.
%
% Written by: Pablo Cuba-Borda, Federal Reserve Board
% Last modified: June 21, 2022. Washington, DC
%==========================================================================

clear; clc; close all;

% Bandwidths for detrending frequency: 40 quarters = 10 years
bw_bw = [40 40];

% Read monthly data
data = readtable('../../data/11302020b.xlsx');

% Transform data to quarterly frequency
dataQ = groupsummary(data,"Dates","quarter","mean");

% Define series to be transformed
series_names = {'FEDFUNDS','PCEPILFE_CCA'};


% Loop over series
for jj=1:length(series_names)
   
    % Monthly series
    use_data = data.(series_names{jj});
    
    % Quarterly series
    use_data_q= dataQ.(['mean_' series_names{jj}]);
   
    % Censor data to 4SD before detrending
    % Avoids outliers to influece trend

    [~,data_q_censored]   = removeOutliers(use_data_q,4,1);    
    
    % Get trends and interpolate on monthly data from 1971m4 onwards
    [m_data_trend, m_data_trend_full]   = fn_local_detrend(data_q_censored,1,series_names{jj},bw_bw(1));
    
    % Detrend variables
    detrended.(series_names{jj}).raw  = use_data  - m_data_trend;
    
    % Collect trend
    trend.(series_names{jj})  = m_data_trend;

    % Add to data table
    data.([series_names{jj} '_TREND']) = trend.(series_names{jj});

end

% Save table to excel file

writetable(data,'../../data/11302020c.xlsx')
