% Construct Local Mean in GDP Growth Rate

%@ Bi-Weight Parameter for local demeaning @
bw_bw = [100 100 100 100]; 

%@ Parameters for ses @
nma = 4; 

% Compute Local Means
rs_place = find(~isnan(data_with_trend(1:end,1)));
rs_with_trend = (data_with_trend(rs_place(1):end,1));
[tmp_rs,se_rs] = bw_trend_se(rs_with_trend,bw_bw(1),nma);

% disp 'ip_with_trend'
ip_place = find(~isnan(data_with_trend(1:end,2)));
ip_with_trend = (data_with_trend(ip_place(1):end,2)); %Start at first number 
[tmp_ip,se_ip] = bw_trend_se(ip_with_trend,bw_bw(2),nma);

pmi_place = find(~isnan(data_with_trend(1:end,3))); 
pmi_with_trend = (data_with_trend(pmi_place(1):end,3));
[tmp_pmi,se_pmi] = bw_trend_se(pmi_with_trend,bw_bw(3),nma);

gdp_place = find(~isnan(data_with_trend(1:end,4)));
gdp_with_trend = (data_with_trend(gdp_place(1):end,4));
[tmp_gdp,se_gdp] = bw_trend_se(gdp_with_trend,bw_bw(4),nma);


for ii=1:3*size(rs_with_trend)
    if ii<3
        m_rs_trend(ii,1) = tmp_rs(1);
    elseif rem(ii,3)==0
        m_rs_trend(ii,1) = tmp_rs(ii/3);
    elseif rem(ii,3)==1
        m_rs_trend(ii,1) = tmp_rs(floor(ii/3)) + (1/3)*(tmp_rs(ceil(ii/3))-tmp_rs(floor(ii/3)));
    elseif rem(ii,3)==2
        m_rs_trend(ii,1) = tmp_rs(floor(ii/3)) + (2/3)*(tmp_rs(ceil(ii/3))-tmp_rs(floor(ii/3)));
    end
end

for ii=1:3*size(ip_with_trend)
    if ii<3
        m_ip_trend(ii,1) = tmp_ip(1);
    elseif rem(ii,3)==0
        m_ip_trend(ii,1) = tmp_ip(ii/3);
    elseif rem(ii,3)==1
        m_ip_trend(ii,1) = tmp_ip(floor(ii/3)) + (1/3)*(tmp_ip(ceil(ii/3))-tmp_ip(floor(ii/3)));
    elseif rem(ii,3)==2
        m_ip_trend(ii,1) = tmp_ip(floor(ii/3)) + (2/3)*(tmp_ip(ceil(ii/3))-tmp_ip(floor(ii/3)));
    end
end

for ii=1:3*size(pmi_with_trend)
    if ii<3
        m_pmi_trend(ii,1) = tmp_pmi(1);
    elseif rem(ii,3)==0
        m_pmi_trend(ii,1) = tmp_pmi(ii/3);
    elseif rem(ii,3)==1
        m_pmi_trend(ii,1) = tmp_pmi(floor(ii/3)) + (1/3)*(tmp_pmi(ceil(ii/3))-tmp_pmi(floor(ii/3)));
    elseif rem(ii,3)==2
        m_pmi_trend(ii,1) = tmp_pmi(floor(ii/3)) + (2/3)*(tmp_pmi(ceil(ii/3))-tmp_pmi(floor(ii/3)));
    end
end
%  
 for ii=1:3*size(gdp_with_trend)
     if ii<3
        m_gdp_trend(ii,1) = NaN;
     elseif rem(ii,3)==0
         m_gdp_trend(ii,1) = tmp_gdp(ii/3);
     elseif rem(ii,3)==1
         m_gdp_trend(ii,1) = NaN;
     elseif rem(ii,3)==2
         m_gdp_trend(ii,1) = NaN;
     end
 end
 
 for ii=1:3*size(gdp_with_trend)
     if ii<3
        m_gdp_trend_full(ii,1)  = tmp_gdp(1);
     elseif rem(ii,3)==0
         m_gdp_trend_full(ii,1) = tmp_gdp(ii/3);
     elseif rem(ii,3)==1
         m_gdp_trend_full(ii,1) = tmp_gdp(floor(ii/3)) + (1/3)*(tmp_gdp(ceil(ii/3))-tmp_gdp(floor(ii/3)));
     elseif rem(ii,3)==2
         m_gdp_trend_full(ii,1) = tmp_gdp(floor(ii/3)) + (2/3)*(tmp_gdp(ceil(ii/3))-tmp_gdp(floor(ii/3)));
     end
 end
 %% Save Trends
 m_rs_trend  = [NaN((rs_place(1) - 1)*3,1); m_rs_trend];
 m_ip_trend  = [NaN((ip_place(1) - 1)*3,1); m_ip_trend];
 m_pmi_trend = [NaN((pmi_place(1) - 1)*3,1); m_pmi_trend];
 m_gdp_trend = [NaN((gdp_place(1) - 1)*3,1); m_gdp_trend];
 m_gdp_trend_full = [NaN((gdp_place(1) - 1)*3,1); m_gdp_trend_full];
 
 %% Adjust trends for partial quarter data
 rs_values.data  = find(~isnan(data_with_trend(:,1)));
 ip_values.data  = find(~isnan(data_with_trend(:,2)));
 pmi_values.data = find(~isnan(data_with_trend(:,3)));
 gdp_values.data = find(~isnan(data_with_trend(:,4)));
 
 rs_values.trend  = find(~isnan(m_rs_trend));
 ip_values.trend  = find(~isnan(m_ip_trend));
 pmi_values.trend = find(~isnan(m_pmi_trend));
 gdp_values.trend = find(~isnan(m_gdp_trend));
 gdp_values_full.trend = find(~isnan(m_gdp_trend_full));
 

 rs_last.trend  = rs_values.trend(end);
 ip_last.trend  = ip_values.trend(end);
 pmi_last.trend = pmi_values.trend(end);
 gdp_last.trend = gdp_values.trend(end); 
 gdp_full_last.trend = gdp_values_full.trend(end);
 
 if(isnan(m_rs_trend(end)))
     m_rs_trend(rs_last.trend+1:end) = m_rs_trend(rs_last.trend);
 end
 
 if(isnan(m_ip_trend(end)))
     m_ip_trend(ip_last.trend+1:end) = m_ip_trend(ip_last.trend);
 end
 
 if(isnan(m_pmi_trend(end)))
     m_pmi_trend(pmi_last.trend+1:end) = m_pmi_trend(pmi_last.trend);
 end
 
 if(isnan(m_gdp_trend(end)))
     m_gdp_trend(gdp_last.trend+3:(3):end) = m_gdp_trend(gdp_last.trend);
 end
 
 if(isnan(m_gdp_trend_full(end)))
     m_gdp_trend_full(gdp_full_last.trend:end) = m_gdp_trend_full(gdp_full_last.trend);
 end

%% Plots
figure;
subplot(2,3,1)
plot([rs_with_trend,tmp_rs])
title('Retail Sales')

subplot(2,3,2)
plot([ip_with_trend,tmp_ip])
title('Industrial Production')

subplot(2,3,3)
plot([pmi_with_trend,tmp_pmi])
title('PMI Composite Index')

subplot(2,3,4)
plot([gdp_with_trend,tmp_gdp])
title('GDP')
