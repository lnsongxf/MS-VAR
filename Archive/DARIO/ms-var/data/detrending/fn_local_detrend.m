function [m_data_trend,m_data_trend_full] = fn_local_detrend(data_with_trend,monthly,figtitle,bw_bw)
% Construct Local Mean in GDP Growth Rate

% bw_bw is the Bi-Weight Parameter for local demeaning: set in quarters
% Stock and Watson baseline is 100 for 25 years


% Parameters for ses 
nma = 4;

% Compute Local Means
data_place = find(~isnan(data_with_trend(1:end,1)));
data_with_trend = (data_with_trend(data_place(1):end,1));
[tmp_data,~] = bw_trend_se(data_with_trend,bw_bw(1),nma);

for ii=1:3*size(data_with_trend)
    if ii<3
        m_data_trend(ii,1) = tmp_data(1);
    elseif rem(ii,3)==0
        m_data_trend(ii,1) = tmp_data(ii/3);
    elseif rem(ii,3)==1
        m_data_trend(ii,1) = tmp_data(floor(ii/3)) + (1/3)*(tmp_data(ceil(ii/3))-tmp_data(floor(ii/3)));
    elseif rem(ii,3)==2
        m_data_trend(ii,1) = tmp_data(floor(ii/3)) + (2/3)*(tmp_data(ceil(ii/3))-tmp_data(floor(ii/3)));
    end
end

m_data_trend_full = m_data_trend;

% Fill trends for quarterly data
if ~monthly
    
    for ii=1:3*size(data_with_trend)
        if ii<3
            m_data_trend(ii,1) = NaN;
        elseif rem(ii,3)==0
            m_data_trend(ii,1) = tmp_data(ii/3);
        elseif rem(ii,3)==1
            m_data_trend(ii,1) = NaN;
        elseif rem(ii,3)==2
            m_data_trend(ii,1) = NaN;
        end
    end
    
    for ii=1:3*size(data_with_trend)
        if ii<3
            m_data_trend_full(ii,1) = tmp_data(1);
        elseif rem(ii,3)==0
            m_data_trend_full(ii,1) = tmp_data(ii/3);
        elseif rem(ii,3)==1
            m_data_trend_full(ii,1) = tmp_data(floor(ii/3)) + (1/3)*(tmp_data(ceil(ii/3))-tmp_data(floor(ii/3)));
        elseif rem(ii,3)==2
            m_data_trend_full(ii,1) = tmp_data(floor(ii/3)) + (2/3)*(tmp_data(ceil(ii/3))-tmp_data(floor(ii/3)));
        end
    end
    
end

%% Save Trends
m_data_trend     = [NaN((data_place(1) - 1)*3,1); m_data_trend];
m_data_trend_full = [NaN((data_place(1) - 1)*3,1); m_data_trend_full];

%% Adjust trends for partial quarter data
trend_values.data  = find(~isnan(data_with_trend(:,1)));

trend_values.trend  = find(~isnan(m_data_trend));
trend_values_full.trend = find(~isnan(m_data_trend_full));


trend_last.trend  = trend_values.trend(end);
trend_full_last.trend = trend_values_full.trend(end);

if(isnan(m_data_trend(end)))
    m_data_trend(trend_last.trend+1:end) = m_data_trend(trend_last.trend);
end

if(isnan(m_data_trend_full(end)))
    m_data_trend_full(trend_full_last.trend:end) = m_data_trend_full(trend_full_last.trend);
end

%% Plots
figure;
plot([data_with_trend,tmp_data])
title(figtitle)

