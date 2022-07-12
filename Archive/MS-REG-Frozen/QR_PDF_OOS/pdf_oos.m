%% PDF - Out-of-Sample
% Author: Francesca Loria
% This Version: May 2020

clear; clc; close all;
addpath(genpath('scripts'));
set(0,'DefaultLegendAutoUpdate','off');

% colors and markers
ColorSet = [([0,0,1]); ([0,0,0,]);([34,139,34]./255);([255,165,0]./255);([128,0,128]./255);([255,255,0]./255);([255,100,0]./255)];
linS = {'-','--','-.',':','-*','-o'};

%% User Input

% Preferences
saveit = 0; % 1 to save plots

% Distribution Preferences (only relevant if compute_dist=1)
delta  = 0.1; % stepsize of grid to evaluate skewed t-density
% Qmatch = [0.10 0.25 0.75 0.9]; % quantiles used to match the skewed t-distribution
% quantiles_dist = 0.1:0.05:0.9; % quantiles to construct distribution
% quantiles      = 0.1:0.4:0.9; % quantiles to plot slopes
Qmatch = [0.1 0.25 0.75 0.9]; % quantiles used to match the skewed t-distribution
% October 5th, 2020: Added 50th quantile to target quantiles in QuantilesInterpolation.m
quantiles_dist = 0.05:0.05:0.95; % quantiles to construct distribution
quantiles      = 0.05:0.45:0.95; % quantiles to plot slopes

% Specifications
% Note: Except for the number of lags and the frequency, specifications can 
% be passed on as vectors if user is interested in running multiple ones
p         =  0; % number of lags (in months/quarters)
frequency =  12; % 12 = monthly, 4 = quarterly (only for US)
sample    = [7]; 
% 1 = 1973:Jan-2020:Apr, 2 = 1986:Jan-2020:Apr,
% 3 = 1973:Jan-2007:Dec, 4 = 1986:Jan-2007:Dec
% 5 = 1973:Jan-2020:Mar, 6 = 1986:Jan-2020:Mar
% 7 = 1973:Jan-2020:May, 8 = 1986:Jan-2020:May
country   = [1]; % 1 = US, 2 = Foreign (TW)
gdpmeas   = [1]; % 1 = Monthly GDP Growth, 2 = BEA GDP Growth (quarterly only)
specs     = [1 11]; % select which model(s) want to run
HH        = [frequency]; % delta ybar_{t,t+h}: average growth between t and t+h (in months/quarters)
scenario  = [1]; 
% 1 = May 28, 2 = April 2, 3 = March 13, 4 = May TB, 5 = May TB FOMC Remarks
% Detrending Preferences
detrend   = [0];  % 0 = already detrended (specify trend to be subtracted), 1 = detrend GDP, 2 = undetrended
dt_window = 10*frequency; % detrending window (in years) for option 1

% Forecast Preferences
% start date of out-of-sample evalutation for expanding window
if frequency==12
%     startOOS  = '2006-Aug';
%     GBfor = 2.2;
%     startOOS  = '2006-Sep';
%     GBfor = 1.4;
%     startOOS  = '2007-Jan';
%     GBfor = 1.5;
    startOOS  = '2007-Mar';
    GBfor = 0.7;
elseif frequency==4
    startOOS    = '2006-Q3';  
end

% Models
% 1  = FF & MF
% 2  = US_ECB_FF & MF (available from 1980-Jan/1980:Q1)
% 3  = KCFSI & MF (available from 1990-Feb/1990:Q2)
% 4  = NFCI & MF
% 5  = NFCI_L & MF
% 6  = NFCI_R & MF
% 7  = NFCI_C & MF
% 8  = FF & FF2 & MF
% 9  = FF & GDP_GROWTH
% 10 = NFCI & GDP_GROWTH
% 11 = MF
% 12 = GDP_GROWTH
% 13 = FF & ADS_Index_031720
% 14 = FF (=0) & MF
% 15 = FF & MF (=0)
% 16 = FF
% 17 = NFCI + GDP Growth (ABG vintage, US quarterly and detrend=1/2 only)
% 18 = FF & MF & Delta FF & Delta MF

% Plotting
FontSize = 14;
FontName = 'Times New Roman';
figSize = [14 8]/2;
numticks = frequency;


%% Run code for different specifications

if frequency==12
    dataformat  = 'yyyy-mmm';
elseif frequency==4
    dataformat  = 'yyyy-qq';
end

% colors
color = [([0,0,1]); ([0,0,0,]);([34,139,34]./255);([255,165,0]./255);([128,0,128]./255);([255,255,0]./255);([255,100,0]./255)];

for d=1:length(detrend)
    
    dt_gdp = detrend(d);
    if dt_gdp==1
        dt='true';
    elseif dt_gdp==0
        dt='given';
    elseif dt_gdp==2
        dt='false';
    end    
    
    for c=1:length(country)
        
        if country(c)==1
            xmin = -14; 
            xmax = 8; 
        elseif country(c)==4
            xmin = -8; 
            xmax = 8; 
        end
        
        for sc=1:length(scenario)
            
            for g=1:length(gdpmeas)
                
                for sam=1:length(sample)
                    
                    %  Set Data and Dates from Scenarios
                    set_scenarios;
                   
                    
                    for j=1:length(HH) % horizon
                        
                        % set horizon
                        hh = HH(j);
                        if frequency==12
                            freq_subfolder= 'Figures_Monthly/';
                        elseif frequency==4
                            freq_subfolder= 'Figures_Quarterly/';
                        end
                        folder = ['../Text_Paper/' freq_subfolder country_filename '/horizon_' num2str(hh) '/PITs/'];
                        if exist(folder,'dir')==0
                            mkdir(folder)
                        end
                        
                        specnum = 0; % initialize specifications number
                        

                        for ss=1:length(specs)
                            
                            specnum = specnum+1;
                            
                            tic % start clock
                            

                            s = specs(ss);

                            % Initialize figure                        
                            if ss==1
                                fig = figure;
                                ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
                            end
                        
                            % Set Dates
                            set_dates;
                            
                            % Set Data
                            set_data;
                            
                            % Set Models with Different Regressors
                            set_models;
                            
                            disp(['Ran specification ' num2str(specnum) ' with ' spec_name ])
                            
                            %% Compute Distribution over Time
                            
                            
                            % Compute moving average for dependent variable
                            if hh==0
                                yh = y;
                            else
                                Yh = filter(ones(1,hh)/hh, 1, y);
                                Yh(1:hh-1) = NaN;
                                yh = Yh(hh+1:end);
                            end
                            Zh = x(1:end-hh,:);
                            Zh = [ones(T-hh,1),Zh];
                                
                            jt=1;
                            % Quantile Regression for Flexible t-Distribution
                            ws = find(datenum(startOOS,dataformat)==dates)-1;
                            period_pdf=datestr(dates(ws+1+hh),dataformat);
                            title_period_pdf=datestr(dates(ws+1+hh),'mmmm yyyy');
                            BQ1 = rq(Zh(1:jt+ws,:), yh(1:jt+ws,:), Qmatch(1));
                            BQ2 = rq(Zh(1:jt+ws,:), yh(1:jt+ws,:), Qmatch(2));
                            BQ3 = rq(Zh(1:jt+ws,:), yh(1:jt+ws,:), Qmatch(3));
                            BQ4 = rq(Zh(1:jt+ws,:), yh(1:jt+ws,:), Qmatch(4));

                            % Plot Slopes and Scatterplots
                            BQ = [BQ1 BQ2 BQ3 BQ4];
                            if s==14
                               BQ(2,:)=0;
                            elseif s==15
                               BQ(3,:)=0;
                            end
                            
                            % Construct Quantiles base on time t information to fit PDF                            
                            YQ_OOS = Zh(jt+ws+hh,:)*BQ;
                            YQ_OOS_unc = Zh(jt+ws+hh,:)*BQ;
                            for rr=1:size(YQ_OOS,2)
                                YQ_OOS_fut(:,rr) = YQ_OOS(:,rr)+yh_trend_fut(jt+ws+hh);
                                YQ_OOS_unc_fut(:,rr) = YQ_OOS_unc(:,rr)+yh_trend_fut(jt+ws+hh);
                            end


                            % Quantiles used to match the skewed t-distribution
                            YY = min(min(YQ_OOS_fut))-20:delta:max(max(YQ_OOS_fut))*2; % grid of points to evaluate skewed t-density
                            % Fit skewed t-distribution to conditional quantiles for each observation
                            ResMatch = Step2match(YQ_OOS_fut, YQ_OOS_unc_fut, Qmatch, YY,delta,Qmatch);

                            % Markers and Colors
                            color=ColorSet(ss,:);

                            % Plot PDF
                            auxx_PST = ResMatch.PST;
                            PST = ResMatch.PST;
                            CDF = ResMatch.CST;

                            % construct specification names
                            if length(var_names)>1
                                plname = sprintf('%s + %s',var_names{:});
                            else
                                plname = sprintf('%s',var_names{:});
                            end

                            realized = yh(jt+ws+hh) + yh_trend(jt+ws+hh);
                            
                            hold(ax1,'on')
                            plot(ax1,ResMatch.YY,PST,linS{ss},'LineWidth',2,'DisplayName',plname,'Color',color);
                            [~,i10] = min(abs(CDF-0.1));
                            [~,i25] = min(abs(CDF-0.25));
                            [~,imedian] = min(abs(CDF-0.5));
                            [~,imode] = max(PST);
                            ym(ss) = round(max(PST),2);
                            plot(ax1,ResMatch.YY(:,imedian),PST(:,imedian)','o','MarkerFaceColor',color,'MarkerEdgeColor','black','Markersize',10,'DisplayName',['Median = ' num2str(round(ResMatch.YY(:,imedian),1)) '\%']);
                            if ss==length(specs)
                                ymax = max(ym)*1.05;
                                plot(ax1,[realized realized],[0 ymax],'r','DisplayName',['Realization = ' num2str(round(realized,1)) '\%']);
                                plot(ax1,[GBfor GBfor],[0 ymax],'-.','Color',([34,139,34]./255),'DisplayName',['Staff Forecast = ' num2str(GBfor) '\%']);
                                set(ax1, 'YLim', [0, ymax])
                                set(ax1, 'XLim', [xmin, xmax])
                                set(ax1,'YTick',0:0.1:ymax)
                                set(ax1,'XTick',xmin:2:xmax)
                                set(ax1,'TickLabelInterpreter','Latex')
                                ylabel(ax1,'PDF','fontsize',12,'interpreter','Latex')
                                set(ax1, 'FontName', 'Times New Roman');
                                set(ax1, 'FontSize', FontSize);
                                set(ax1,'Layer','top')
                                hold off
                                hL=legend(ax1,'-DynamicLegend','Location','NorthWest','interpreter','Latex');
                                legend boxoff
                                %title(['Conditional Distributions of Average ',country_name,' GDP Growth over the Next ',num2str(hh),' ',funit,' in ' char(period_pdf)],'fontsize',FontSize,'interpreter','Latex')
                                title(char(title_period_pdf),'fontsize',20,'interpreter','Latex')
                                set(fig,'PaperOrientation','portrait');
                                set(fig, 'PaperSize', figSize);
                                set(fig, 'PaperUnits', 'inches');
                                set(fig, 'Units','inches');
                                set(fig, 'PaperPositionMode', 'auto');
                                set(fig, 'Position', [figSize(1)/5 figSize(2)/3 figSize(1) figSize(2)]);
                                tightfig;
                                if saveit==1
                                    savefig(fig,[folder 'PDF_OOS_' '___Period=' char(period_pdf) '___Country=' country_filename '___GDP=' gdp_filename '___Scenario=' scen_name '___Lags=' num2str(p) '___Detrended=' dt '___Sample=' startdate '_to_' enddate]);
                                    print('-dpdf',fig,[folder 'PDF_OOS_' '___Period=' char(period_pdf) '___Country=' country_filename '___GDP=' gdp_filename  '___Scenario=' scen_name '___Lags=' num2str(p) '___Detrended=' dt '___Sample=' startdate '_to_' enddate],'-bestfit');
                                    saveas(fig,sprintf('%s.png',[folder 'PDF_OOS_' '___Period=' char(period_pdf) '___Country=' country_filename '___GDP=' gdp_filename '___SpecWith=' spec_name '___Scenario=' scen_name '___Lags=' num2str(p) '___Detrended=' dt '___Sample=' startdate '_to_' enddate ]));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end



