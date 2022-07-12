%% Define Models

if s==1
    xx = [FF MF lagmatrix(FF,1:p) lagmatrix(MF,1:p)];
    spec_name = 'ff_mf';
    var_names = {'Financial Factor', 'Macroeconomic Factor'};
    var_shortnames = {'ff','mf'};
elseif s==2
    xx = [US_ECB_FF MF lagmatrix(US_ECB_FF,1:p) lagmatrix(MF,1:p)];
    spec_name = 'ecbff_mf';
    var_names = {'Financial Factor (ECB)', 'Macroeconomic Factor'};
    var_shortnames = {'ecbff','mf'};
elseif s==3
    xx = [KCFSI MF lagmatrix(KCFSI,1:p) lagmatrix(MF,1:p)];
    spec_name = 'kansasff_mf';
    var_names = {'Financial Factor (Kansas Fed)', 'Macroeconomic Factor'};
    var_shortnames = {'kansasff','mf'};
elseif s==4
    xx = [NFCI MF lagmatrix(NFCI,1:p) lagmatrix(MF,1:p)];
    spec_name = 'nfci_mf';
    var_names = {'NFCI', 'Macroeconomic Factor'};
    var_shortnames = {'nfci','mf'};
elseif s==5
    xx = [NFCI_L MF lagmatrix(NFCI_L,1:p) lagmatrix(MF,1:p)];
    spec_name = 'nfcil_mf';
    var_names = {'NFCI (Leverage)', 'Macroeconomic Factor'};
    var_shortnames = {'nfcil','mf'};
elseif s==6
    xx = [NFCI_R MF lagmatrix(NFCI_R,1:p) lagmatrix(MF,1:p)];
    spec_name = 'nfcir_mf';
    var_names = {'NFCI (Risk)', 'Macroeconomic Factor'};
    var_shortnames = {'nfcir','mf'};
elseif s==7
    xx = [NFCI_C MF lagmatrix(NFCI_C,1:p) lagmatrix(MF,1:p)];
    spec_name = 'nfcic_mf';
    var_names = {'NFCI (Credit)', 'Macroeconomic Factor'};
    var_shortnames = {'nfcic','mf'};
elseif s==8
    xx = [FF FF2 MF lagmatrix(FF,1:p) lagmatrix(FF2,1:p) lagmatrix(MF,1:p)];
    spec_name = 'ff_ff2_mf';
    var_names = {'Financial Factor', 'Financial Factor (2nd)', 'Macroeconomic Factor'};
    var_shortnames = {'ff','ff2','mf'};
elseif s==9
    xx = [FF GDP_GROWTH lagmatrix(FF,1:p) lagmatrix(GDP_GROWTH,1:p)];
    spec_name = 'ff_gdpg';
    var_names = {'Financial Factor', 'GDP Growth'};
    var_shortnames = {'ff','gdpg'};
elseif s==10
    xx = [NFCI GDP_GROWTH lagmatrix(NFCI,1:p) lagmatrix(GDP_GROWTH,1:p)];
    spec_name = 'nfci_gdpg';
    var_names = {'NFCI', 'GDP Growth'};
    var_shortnames = {'nfci','gdpg'};
elseif s==11
    xx = [MF lagmatrix(MF,1:p)];
    spec_name = 'mf';
    var_names = {'Macroeconomic Factor'};
    var_shortnames = {'mf'};
elseif s==12
    xx = [GDP_GROWTH lagmatrix(GDP_GROWTH,1:p)];
    spec_name = 'gdpg';
    var_names = {'GDP Growth'};
    var_shortnames = {'gdpg'};
elseif s==13
    xx = [FF ADS_Index_031720 lagmatrix(FF,1:p) lagmatrix(ADS_Index_031720,1:p)];
    spec_name = 'ff_phillyads';
    var_names = {'Financial Factor', 'Philly Fed ADS Index (03/17/20 Vintage)'};
    var_shortnames = {'ff','phillyads'};
elseif s==14
    xx = [FF MF lagmatrix(FF,1:p) lagmatrix(MF,1:p)];
    spec_name = 'noff_mf';
    var_names = {'Financial Factor ($f_t=0$)', 'Macroeconomic Factor'};
    var_shortnames = {'noff','mf'};
elseif s==15
    xx = [FF MF lagmatrix(FF,1:p) lagmatrix(MF,1:p)];
    spec_name = 'ff_nomf';
    var_names = {'Financial Factor', 'Macroeconomic Factor (${mf}_t=0$)'};
    var_shortnames = {'ff','nomf'};
elseif s==16
    xx = [FF lagmatrix(FF,1:p)];
    spec_name = 'ff';
    var_names = {'Financial Factor'};
    var_shortnames = {'mf'};
elseif s==17
    xx = [NFCI_ABG GDP_GROWTH lagmatrix(NFCI_ABG,1:p) lagmatrix(GDP_GROWTH,1:p)];
    spec_name = 'nfci_gdpg_abg';
    var_names = {'NFCI (ABG)', 'GDP Growth (ABG)'};
    var_shortnames = {'nfci','gdpg'};
    Qmatch = [0.05 0.25 0.75 0.95]; % quantiles used to match the skewed t-distribution
    quantiles_dist = 0.05:0.05:0.95; % quantiles to construct distribution
    quantiles      = 0.05:0.45:0.95; % quantiles to plot slopes
elseif s==18
    %xx = [FF MF [mean(diff(FF));diff(FF)] [mean(diff(MF));diff(MF)]];
    FFG = diff(FF)./FF(1:end-1);
    MFG = diff(MF)./MF(1:end-1);
    xx = [FF MF [mean(FFG);FFG] [mean(MFG);MFG]];
    %xx = [FF MF GDP_GROWTH lagmatrix(FF,1:p) lagmatrix(MF,1:p) lagmatrix(GDP_GROWTH,1:p)];
    %xx = [FF MF lagmatrix(FF,1:p) lagmatrix(MF,1:p)];
    spec_name = 'ff_mf_dff_dmf';
    var_names = {'Financial Factor', 'Macroeconomic Factor','$\Delta$ Financial Factor', '$\Delta$ Macroeconomic Factor'};
    var_shortnames = {'ff','mf','dff','dmf'};
elseif s==19
    xx = [FF MF1 lagmatrix(FF,1:p) lagmatrix(MF1,1:p)];
    spec_name = 'ff_mf1';
    var_names = {'Financial Factor', 'Macroeconomic Factor without IP'};
    var_shortnames = {'ff','mf1'};
elseif s==20
    xx = [FF MF2 lagmatrix(FF,1:p) lagmatrix(MF2,1:p)];
    spec_name = 'ff_mf2';
    var_names = {'Financial Factor', 'Macroeconomic Factor without PMI-NEO'};
    var_shortnames = {'ff','mf2'};
elseif s==21
    xx = [FF MF3 lagmatrix(FF,1:p) lagmatrix(MF3,1:p)];
    spec_name = 'ff_mf3';
    var_names = {'Financial Factor', 'Macroeconomic Factor without RS'};
    var_shortnames = {'ff','mf3'};
elseif s==22
    xx = [FF MF4 lagmatrix(FF,1:p) lagmatrix(MF4,1:p)];
    spec_name = 'ff_mf4';
    var_names = {'Financial Factor', 'Macroeconomic Factor without GDP'};
    var_shortnames = {'ff','mf4'};
elseif s==23
    xx = [FF MF5 lagmatrix(FF,1:p) lagmatrix(MF5,1:p)];
    spec_name = 'ff_mf5';
    var_names = {'Financial Factor', 'Macroeconomic Factor without UCLAIMS'};
    var_shortnames = {'ff','mf5'};
elseif s==24
    ADS = data_other.ADS;    
    xx = [FF ADS lagmatrix(FF,1:p) lagmatrix(ADS,1:p)];
    spec_name = 'ff_ads';
    var_names = {'Financial Factor', 'ADS Index'};
    var_shortnames = {'ff','ads'};
end

%% Housekeeping

slopes = (1:length(var_shortnames));
x = xx(1+p:end,:); 
T = length(x);
y = GDP_GROWTH(1+p:end); 

%% Take care of trends
if dt_gdp==0
    y_dfm_trend = GDPG_TREND_DFM(1+p:end);
    Yh_trend = filter(ones(1,hh)/hh, 1, y_dfm_trend);
    Yh_trend(1:hh-1) = NaN;
    yh_trend = Yh_trend(hh+1:end); 
    y = y-y_dfm_trend;
elseif dt_gdp==1
    y_trend = zeros(size(y));
    y_trend(1:dt_window) = mean(y(1:dt_window));
    for jj=dt_window+1:T
        y_trend(jj) = mean(y(jj-dt_window+1:jj));
    end
    Yh_trend = filter(ones(1,hh)/hh, 1, y_trend);
    Yh_trend(1:hh-1) = NaN;
    yh_trend = Yh_trend(hh+1:end); 
    y = y-y_trend;
elseif dt_gdp==2
    Yh_trend = zeros(T,1);
    Yh_trend(1:hh-1) = NaN;
    yh_trend = Yh_trend(hh+1:end); 
    %y = y;
end

yh_trend_fut = [yh_trend;yh_trend(end)*(ones(hh,1))];
