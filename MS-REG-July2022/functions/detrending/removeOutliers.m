%% Header
% This function removes outliers (+/- SD), where SD is a positive integer.
%
% Last modified: 10 Jan 2019
% By: Dawson Miller
%==========================================================================

function [data_no_outliers, data_censored] = removeOutliers(data_in,SD,censor)
    data_no_outliers = data_in;
    data_censored    = data_in;

    mean_data = mean(data_in,'omitnan');
    std_data = std(data_in,'omitnan');

    for ii = 1:max(size(data_in))
         if abs(data_in(ii) - mean_data)/std_data > SD
             if censor
                data_censored(ii) = sign(data_no_outliers(ii))*SD*std_data;
             else
                data_no_outliers(ii) = NaN;
             end
         end
    end

end