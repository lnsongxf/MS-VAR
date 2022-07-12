% Understanding Growth-at-Risk: A Markov-Switching Approach
% Model with f_t and m_t endogenous
%
% This script creates the QW-CRPS test statistics and comparison tables using
% the out-of-sample (OOS) results computed by OOS_blender_counterfactual.m
%
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
opt.transprob      = 1;     % 1 = endogenous, otherwise exogenous
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
mymodels = [101:108];

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
        ind = ind+1;
        modelspec = imodel;

        % Model spectag
        spectag = fSpecLabel(modelspec);

        % Model spectag
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


%% 
% Vector of dates for the full sample and for the available estimation sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

% Index of dates for plotting
dates = datenum((datetime(efirstok,'InputFormat',inputformat)+calmonths(nlags)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';
sd      = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates);
ed      = find(datenum((datetime(end_plot,'InputFormat',inputformat)))==dates);
sd_full = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates_full);
ed_full = find(datenum((datetime(end_plot,'InputFormat',inputformat)))==dates_full);


% clc;
% logfolder    = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags)  '/' spectag '/'  ];
% diary ([log_folder 'reduced.txt'])

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
    for j=1:ind
        left_score = squeeze(QW_CRPS(i,j,7,:));
        if j==bm(i)
            plot(dates(sd:ed),left_score(sd:ed),'k-*','LineWidth', 2.5); 
        else
            mycolor = colors(j+5*j,:);            
            plot(dates(sd:ed),left_score(sd:ed),'Color',mycolor,'LineWidth', 2.5); 
        end
        hold on
        specname{j} = strrep(results{i,j}.spectag,'_',' ');
    end
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
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', FontSize);
    set(gca,'Layer','top')
    set(gca,'TickLabelInterpreter','Latex')
    legend(specname,'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
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
models = arrayfun(@num2str, mymodels, 'UniformOutput', 0);
for i=myforecasts
    if i==1
        fprintf('Direct Model \n')
    else
        fprintf('Iterated Model \n')
    end
    fprintf(' \n')
    basetag = results{i,1}.spectag;
    T1 = array2table(round(squeeze(Score(i,:,:)),dc)','VariableNames',models,'RowName',Rows);
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
        T1 = table(T1,'VariableNames',{[basetag ' vs. ' spectag]});
        if j~=bm(i)
            disp(T1);
        end
    end
end



%% Comparison Table Between Direct and Iterated Best Models
% If T-Stat is positive, iterated model favored

% Benchmark is Iterated Model
b = 2;
m = 1;

for j=1:ind
    for n = 1:nparts
        benchm_di = QW_CRPS(b,bm(b),n,:);
        [DMtstat_di(j,n),nanny] = dieboldmariano(squeeze(QW_CRPS(m,j,n,:)), squeeze(benchm_di), opt.hh + 1);
        DMpvalue_di(j,n)        = 2*(1-normcdf(abs(DMtstat_di(j,n))));
        relScore_di(j,n)        = mean(squeeze(QW_CRPS(m,j,n,nanny))) / mean(benchm_di);
        Score_di(j,n)           = mean(squeeze(QW_CRPS(m,j,n,:)));
    end
end


dc = 2; % rounding
basetag = results{1,bm(b)}.spectag;
for j=bm(m)
    T1 = array2table([round(squeeze(relScore_di(j,:)),dc)' round(squeeze(DMtstat_di(j,:)),dc)' round(squeeze(DMpvalue_di(j,:)),dc)'],'VariableNames',{'Rel. Score','T-Stat','p-value'},'RowName',Rows);
    spectagi = results{m,j}.spectag;
    spectagd = results{b,j}.spectag;
%     T1 = table(T1,'VariableNames',{['Iterated ' spectagi ' vs. Direct' spectagd]});
    T1 = table(T1,'VariableNames',{['Best Iterated vs. Best Direct']});
    disp(T1);
end

% diary off 

%%
rise_exit
