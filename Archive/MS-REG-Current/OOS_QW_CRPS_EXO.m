% Understanding Growth-at-Risk: A Markov-Switching Approach
% Model with f_t and m_t endogenous
%
% This script creates the QW-CRPS test statistics and comparison tables using
% the out-of-sample (OOS) results computed by OOS_blender_counterfactual.m
%
% Comparison between endogenous and exogenous tp for direct model
%=====================================================================


%% housekeeping
clear; clc;
close all;

% Important paths
addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');         % RISE Toolbox
% addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_Tbx_Beta-master');         % RISE Toolbox
addpath(genpath('functions'));      % Main functions
addpath(genpath('scripts'));        % Main scripts
addpath(genpath('auxtools'));       % Tools for plotting and other
addpath(genpath('cbrewer'));        % Color mixing tools

% Load RISE
rise_startup()

% Date formats
inputformat   = 'yyyy-MMM';
dataformat    = 'yyyy-mmm';

%========================================
% USER INPUT
%========================================
% Select which forecast types to run: 1 = direct, 2 = iterated
myforecasts = [1 2];
opt.simtype        = 1;     % 1 = forecast, 0 = drawst
opt.hh             = 12;    % forecast horizon
opt.const          = 1;     % 1 = have a constant in transition probability
opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
opt.paramsUse      = 'mode';

% Select which models to run
% Model	  C	 A0	 A1	 S	GDP(-1) restriction
% 101     s	 n	 n	 n	Yes
% 102     s	 n	 n	 s	Yes
% 103     s	 s	 n	 n	Yes
% 104     s	 n	 s	 n	Yes
% 105     s	 s	 n	 s	Yes
% 106     s	 n	 s	 s	Yes
% 107     s	 s	 s	 s	Yes
% 108  s(3)  n   n  s(3) Yes
% 201     s	 n	 n	 n	No
% 202     s	 n	 n	 s	No
% 203     s	 s	 n	 n	No
% 204     s	 n	 s	 n	No
% 205     s	 s	 n	 s	No
% 206     s	 n	 s	 s	No
% 207     s	 s	 s	 s	No
% 208  s(3)  n   n  s(3) No
mymodels = [105];

% Select list of counterfactuals
mycounterfactuals = []; %[2 3];

% Simulation options
opt.nDraws       = 5000;    % Number of simulations for each time period
opt.nParamDraws  = 1000;    % Draws from posterior. has to be less than options.N;

% MCMC options
options = struct();
options.alpha = 0.234;
options.thin = 10;
options.burnin = 10^3;
options.N = 2*10^4;
options.nchain = 1;

% VAR configuration
% important: if nlags>1 the routine that pulls the parameters
% (scriptParams_FULL_companion.m) does not work for MCMC
nlags=1;
exog={};
constant=true;
panel=[];

% Declare first sample and last sample
% (needs to match msreg_2reg_endoprob_SVAR_FULL_OOS.m)
efirst = 'nov 01 1999';
efirstok = datestr(efirst,dataformat);
  
% Data vintage, sample and country selection
datafilename = '11302020';
sheetuse     = 'US_DFM';
start_date   = '1973-Jan';
end_date     = '2019-Nov';
startdb      = '1973M1';
enddb        = '2019M1';
start_plot = efirstok;
end_plot   = '2019-Mar'; % also ends evaluation sample


% 11302020_US_DFM_1973M12019M11_C1N0H12CT_A0T_A1_SIGT_restr_GDP_results.mat
ifi = 0;
for i=myforecasts
    ifi = ifi+1;
    opt.dir_or_it = i; % 1 = direct, 2 = iterated
    ind = 0;
    for imodel=mymodels
        for tp=[1 0]
            ind = ind+1;
            modelspec = imodel;

            % Model spectag
            spectag = fSpecLabel(modelspec);

            % Model spectag
            opt.transprob = tp;
            if opt.transprob==1
                tag ='';
            else
                tag ='_EXO';
            end

            % Determine simulation type
            if opt.simtype==1
                simtype = 'forecastst';
            else
                simtype = 'drawst';
            end

            %========================================
            % MODEL SPECIFIC VARIABLES
            %========================================
            foldertag = 'OOS';

            if opt.dir_or_it ==1
                performfolder = ['results/Perfomance/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
            elseif opt.dir_or_it ==2
                performfolder = ['results/Perfomance/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
            end
            
            % Model name and title
            % Note: iterated model posterior mode depends on opt.hh through the estimation sample
            modelname = [datafilename '_' sheetuse '_' startdb enddb ...
            '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];

            %% Load Results
            results{ifi,ind} = load([performfolder modelname '.mat'],'forcperf');
            results{ifi,ind}.performfolder = performfolder;
            results{ifi,ind}.modelname = modelname;
            results{ifi,ind}.spectag = spectag;

        end
    end
end


%% Dates

% Vector of dates for the full sample and for the available estimation sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

% Index of dates for plotting
dates = datenum((datetime(efirstok,'InputFormat',inputformat)+calmonths(nlags)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';
sd      = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates);
ed      = find(datenum((datetime(end_plot,'InputFormat',inputformat)))==dates);
sd_full = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates_full);
ed_full = find(datenum((datetime(end_plot,'InputFormat',inputformat)))==dates_full);

%% Gneiting and Ranjan (JBES, 2011)

% choose benchmark model
bm(1) = find(mymodels==105); % direct
bm(2) = find(mymodels==105); % iterated


Rows={'Uniform','Center','Tail','Right Tail','Left Tail','Right Tail Pure','Left Tail  Pure'};
for i=myforecasts
    for j=1:ind
        QW_CRPS(i,j,1,:) = cell2mat(results{i,j}.forcperf.uni_qwcrps_time(sd:ed))';
        QW_CRPS(i,j,2,:) = cell2mat(results{i,j}.forcperf.center_qwcrps_time(sd:ed))';
        QW_CRPS(i,j,3,:) = cell2mat(results{i,j}.forcperf.tails_qwcrps_time(sd:ed))';
        QW_CRPS(i,j,4,:) = cell2mat(results{i,j}.forcperf.righttail_qwcrps_time(sd:ed))';
        QW_CRPS(i,j,5,:) = cell2mat(results{i,j}.forcperf.lefttail_qwcrps_time(sd:ed))';
        QW_CRPS(i,j,6,:) = cell2mat(results{i,j}.forcperf.righttailpure_qwcrps_time(sd:ed))';
        QW_CRPS(i,j,7,:) = cell2mat(results{i,j}.forcperf.lefttailpure_qwcrps_time(sd:ed))';
    end
end


nparts = size(QW_CRPS,3); % Uniform, Center, Tails, Right Tail, Left Tail
tparts = size(QW_CRPS,4); % number of forecasts available
for i=myforecasts
    for j=1:ind
        for n = 1:nparts
            benchm = QW_CRPS(i,bm(i),n,:);
            [DMtstat(i,j,n),nanny] = dieboldmariano(squeeze(QW_CRPS(i,j,n,:)), squeeze(benchm), opt.hh + 1);
            DMpvalue(i,j,n)        = 2*(1-normcdf(abs(DMtstat(i,j,n))));
            relScore(i,j,n)        = mean(squeeze(QW_CRPS(i,j,n,nanny))) / mean(benchm);
            Score(i,j,n)           = mean(squeeze(QW_CRPS(i,j,n,:)));
        end
    end
end

%% Scores Across Time - Direct Model

% Figure options
figureTitleTag = [strrep(spectag,'_',' ') ', ' opt.paramsUse ', ' strrep(sheetuse,'_',' ') ', h=' num2str(opt.hh) ', ' startdb ':' enddb];
FontSize  = 16;
numticks  = 24;
figSize1   = [12 6];
linestyle = {'-','--',':'};
colors    = cbrewer('div', 'RdYlBu', 64);
colors2   = cbrewer('qual', 'Set1', 8);
left_color= [0 0 0];
right_color= colors(2,:);

for i=myforecasts
    opt.dir_or_it = i; % 1 = direct, 2 = iterated
    fig1 = figure;
    figSize1 = [12 6];
    left_scoret_endo = squeeze(QW_CRPS(i,1,7,:));
    left_scoret_exo = squeeze(QW_CRPS(i,2,7,:));
    plot(dates(sd:ed),left_scoret_endo(sd:ed),'Color','k','LineWidth', 2.5); 
    hold on
    plot(dates(sd:ed),left_scoret_exo(sd:ed),'Color',[0.6 0.1 0.1],'LineWidth', 2.5); 
    if opt.dir_or_it==1
        title(['Direct: ' figureTitleTag ],'Interpreter','Latex');
    else
        title(['Iterated: ' figureTitleTag ],'Interpreter','Latex');
    end
    axis tight
    ylabel('QW-CRPS Scores - Left Tail','fontsize',10,'interpreter','Latex')
    xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
    set(gca,'XTick',datenum(dates(sd:numticks:ed)))
    datetick('x','yyyy','keepticks')
    set(gca, 'XLim', [dates(sd), dates(ed)])
    hold off
    hBands=recessionplot;
    legend('Endogenous','Exogenous','Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', FontSize);
    set(gca,'Layer','top')
    set(gca,'TickLabelInterpreter','Latex')
    % Set figure appereance
    set(fig1,'PaperOrientation','portrait');
    set(fig1, 'PaperSize', figSize1);
    set(fig1, 'PaperUnits', 'inches');
    set(fig1, 'Units','inches');
    set(fig1, 'Position', [figSize1(1)/5 figSize1(2)/5 figSize1(1) figSize1(2)]);
    tightfig; 
end



%% Scores Across Models

dc = 2; % rounding
for i=myforecasts
    if i==1
        fprintf('Direct Model \n')
    else
        fprintf('Iterated Model \n')
    end
    fprintf(' \n')
    basetag = results{i,1}.spectag;
    T1 = array2table(round(squeeze(Score(i,:,:)),dc)','VariableNames',{'Endogenous','Exogenous'},'RowName',Rows);
    disp(T1);
end
    
%% Comparison Table Between Models
% If T-Stat is positive, baseline model favored

dc = 2; % rounding
for i=myforecasts
    if i==1
        fprintf('Direct Model \n')
    else
        fprintf('Iterated Model \n')
    end
    fprintf(' \n')
    basetag = results{i,bm(i)}.spectag;
    for j=1:ind
        T1 = array2table([round(squeeze(relScore(i,j,:)),dc) round(squeeze(DMtstat(i,j,:)),dc) round(squeeze(DMpvalue(i,j,:)),dc)],'VariableNames',{'Rel. Score','T-Stat','p-value'},'RowName',Rows);
        spectag = results{i,j}.spectag;
        T1 = table(T1,'VariableNames',{[basetag ' Endo vs. Exo']});
        if j~=bm(i)
            disp(T1);
        end
    end
end

%%
rise_exit
