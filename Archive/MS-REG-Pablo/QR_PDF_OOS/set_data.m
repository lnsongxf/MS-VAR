%% Collect Data

% Dependent Variable
if gdpmeas(g)==1
    if country(c)==1
        GDP_GROWTH     = data_dfm.GDPG_US;
        if dt_gdp==0
            GDPG_TREND_DFM = data_dfm.TREND_US;
        end
    elseif country(c)==2
        GDP_GROWTH     = data_dfm.GDPG_FN;
        if dt_gdp==0
            GDPG_TREND_DFM = data_dfm.TREND_FN;
        end
    end
    gdp_name = 'GDP Growth';
    gdp_filename = 'GDP';
elseif gdpmeas(g)==2 && frequency==4 && country(c)==1
    GDP_GROWTH  = data_other.A191RL1Q225SBEA;
    if dt_gdp==0
        GDPG_TREND_DFM = data_dfm.TREND_US;
    end
    gdp_name = 'GDP Growth';
    gdp_filename = 'GDP_BEA';  
else
    error('Specification not available!')
end

% DFM Regressors
if country(c)==1
    MF    = data_dfm.MF_US;
elseif country(c)==2
    MF    = data_dfm.MF_FN;
end
FF	      = data_dfm.FF;
FF2       = data_dfm.FF2;
if specs(ss)>18 && specs(ss)<24
    MF1	      = data_dfm.MF1;
    MF2	      = data_dfm.MF2;
    MF3	      = data_dfm.MF3;
    MF4	      = data_dfm.MF4;
    MF5	      = data_dfm.MF5;
end

% Other Regressors
US_ECB_FF = data_other.US_ECB_FF;
NFCI      = data_other.NFCI;
NFCI_L    = data_other.NFCI_L;
NFCI_R    = data_other.NFCI_R;
NFCI_C    = data_other.NFCI_C;
KCFSI     = data_other.KCFSI;
ADS_Index_031720 = data_other.ADS_Index_031720;

% Overwrite GDP Growth with ABG vintage for Specification 17
if specs(ss)==17
    if frequency==4 && detrend(d)~=0
       NFCI_ABG  = data_other.NFCI_ABG;
       GDP_GROWTH = data_other.GDPG_ABG;
    else
        error('Specification not available!')
    end
end

