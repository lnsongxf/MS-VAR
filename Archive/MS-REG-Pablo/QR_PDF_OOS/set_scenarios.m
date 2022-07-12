%% Collect Vintages

% Data and Name
if scenario(sc)==1
    vin_data = 'MacroRisk_May28_US_only.xlsx';
    scen_name = 'Baseline';
    scen_name_plot = 'Baseline';
elseif scenario(sc)==2
    vin_data = 'MacroRisk_April2_US_only.xlsx';
    scen_name = 'Apr2';
    scen_name_plot = 'April 2';
elseif scenario(sc)==3
    vin_data = 'MacroRisk_March13_US_only.xlsx';
    scen_name = 'Mar13';
    scen_name_plot = 'March 13';
elseif scenario(sc)==4
    vin_data = 'MacroRisk_MayTB.xlsx';
    scen_name = 'MayTB';
    scen_name_plot = 'May TB';
elseif scenario(sc)==5
    vin_data = 'MacroRisk_FOMCRemarks_May.xlsx';
    scen_name = 'MayTB_FOMCRemarks';
    scen_name_plot = 'May TB FOMC Remarks';
elseif scenario(sc)==6
    vin_data = 'MacroRisk_FOMCRemarks_June.xlsx';
    scen_name = 'JuneTB_FOMCRemarks';
    scen_name_plot = 'June TB FOMC Remarks';
elseif scenario(sc)==7
    vin_data = 'MacroRisk_Francesca_Briefing.xlsx';
    scen_name = 'August_24_Briefing';
    scen_name_plot = 'August 24 Briefing';
elseif scenario(sc)==8
    vin_data = 'MacroRisk_Francesca_Briefing_2020Q2GDP_Uncensored.xlsx';
    scen_name = 'August_24_Briefing_2020Q2GDP_Uncensored';
    scen_name_plot = 'August 24 Briefing (2020Q2 GDP Uncensored)';
elseif scenario(sc)==9
    vin_data = 'MacroRisk_Francesca_Briefing_Uncensored.xlsx';
    scen_name = 'August_24_Briefing_Uncensored';
    scen_name_plot = 'August 24 Briefing (Uncensored)';
elseif scenario(sc)==10
    vin_data = 'MacroRisk_Francesca_Briefing_old.xlsx';
    scen_name = 'August_24_Briefing_Uncensored_Old';
    scen_name_plot = 'August 24 Briefing (Old)';
elseif scenario(sc)==11
    vin_data = 'MacroRisk_Francesca_Briefing_Uncensored_08182020.xlsx';
    scen_name = 'August_24_Briefing_Uncensored_08182020';
    scen_name_plot = 'August 24 Briefing (Uncensored)';
elseif scenario(sc)==12
    vin_data = 'MacroRisk_Francesca_Briefing_Final.xlsx';
    scen_name = 'August_24_Briefing_Final';
    scen_name_plot = 'August 24 Briefing';
elseif scenario(sc)==13
    vin_data = 'MacroRisk_FS_Briefing.xlsx';
    scen_name = 'FS_Briefing';
    scen_name_plot = 'FS Briefing';
% elseif scenario(sc)==14
%     vin_data = 'MacroRisk_AssessingMacroRisk.xlsx';
%     scen_name = 'AssessingMacroRisk';
%     scen_name_plot = 'Assessing Macroeconomic Risk';
end

% Select Sheets
if frequency==12
    if sample(sam)==1 || sample(sam)==3 || sample(sam)==5 || sample(sam)==7
        dfm_sheet = 'DFM_73_Monthly';
    else
        dfm_sheet = 'DFM_86_Monthly';       
    end
    data_sheet = 'Data_Monthly';
elseif frequency==4
    if sample(sam)==1 || sample(sam)==3 || sample(sam)==5|| sample(sam)==7
        dfm_sheet = 'DFM_73';
    else
        dfm_sheet = 'DFM_86';
    end
    data_sheet = 'Data';
end

% Load DFM-Estimated GDP, Macroeconomic Factor and Financial Factor
data_dfm_org = readtable(vin_data,'Sheet',dfm_sheet);
% Load Other Variables
data_other_org = readtable(vin_data,'Sheet',data_sheet);


% Country Names
if country(c)==1
    country_filename = 'US'; % used in filename
    country_name = 'United States'; % used in plotting
    country_shortname = 'U.S.'; % used in plotting
elseif country(c)==2
    country_filename = 'Foreign_TW';
    country_name = 'Foreign Economy';
    country_shortname = 'Foreign'; % used in plotting
end