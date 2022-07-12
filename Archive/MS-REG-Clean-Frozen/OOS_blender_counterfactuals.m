% Understanding Growth-at-Risk: A Markov-Switching Approach
% Model with f_t and m_t endogenous
%
% This script blends the out-of-sample (OOS) results computed by
% msreg_2reg_endoprob_SVAR_FULL_OOS.m
%
% There are two estimation options
%  1. Direct: estimates the measurement equation at horizon hh.
%  2. Iterated: estimates the contemporanous measurement equation and
%     then iterates the predictive density through horizon hh.
%=====================================================================



%% housekeeping
clear; clc;
close all;

% Important paths
addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');  % RISE Toolbox
addpath(genpath('functions'));      % Main functions
addpath(genpath('scripts'));        % Main scripts
addpath(genpath('textools'));       % Tools to produce tables
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
opt.dir_or_it      = 1;     % 1 = direct, 2 = iterated
opt.transprob      = 1;     % 1 = endogenous, otherwise exogenous 
opt.simtype        = 1;     % 1 = forecast, 0 = drawst
opt.hh             = 12;    % forecast horizon
opt.const          = 1;     % 1 = have a constant in transition probability
opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
opt.paramsUse      = 'mode';
opt.saveit         = 0;     % 1 = save figures
opt.showfig        = 1;     % 1 = show figures
opt.selfig         = [1 1 1 1 1 1 1 1 1 1 1];   % 1 = vector of figures


% Figure selection options
% 1 = quantiles predictive distribution
% 2 = A0 parameters time evolution
% 3 = A1 parameters time evolution
% 4 = C parameters time evolution
% 5 = SIG parameters time evolution
% 6 = transition probability parameters time evolution
% 7 = predictive score
% 8 = PITs
% 9 = CDF of PITs
%10 = Counterfactual predictive distributions
%11 = Counterfactual density cuts

% Select which models to run
mymodels = [11];

% Select list of counterfactuals
mycounterfactuals = [2 3];

% Simulation options
opt.nDraws       = 50;    % Number of simulations for each time period
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
efirst = 'dec 01 1999';
esample = 'oct 01 2020';
maxwin = months(efirst,esample);
efirstok = datestr(efirst,dataformat);


for imodel=mymodels
    tic;
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
    
    % Data vintage, sample and country selection
    datafilename = '11302020';
    sheetuse     = 'US_DFM';
    start_date   = '1973-Jan';
    
    
    %========================================
    % MODEL SPECIFIC VARIABLES
    %========================================
    foldertag = 'OOS';
    
    if opt.dir_or_it ==1
        fcstype= 'direct';
        varlist={'FF','MF','GDPGH'}';
        
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' foldertag '/' ];
        IS_paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' ];
        figfolder   = ['results/Figures/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/' foldertag '/'];

        % Folder to store counterfactuals
        counterfolder = ['results/Counterfactuals/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' foldertag '/' ];

    elseif opt.dir_or_it ==2
        fcstype= 'iterated';
        varlist={'FF','MF','GDPG'}';
        
        % Folders Containing Parameters
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' foldertag '/' ];
        IS_paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/'];
        
        % Folder To store figures
        figfolder   =  ['results/Figures/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' foldertag '/'];
        
        % Folder to store counterfactuals
        counterfolder = ['results/Counterfactuals/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' foldertag '/' ];
    end
    
    %% Create Folders
    if exist(counterfolder,'dir')==0;  mkdir(counterfolder); end
    
    
    %% Read data
    rawdata = readtable(['data/' datafilename '.xlsx'],'Sheet',sheetuse,'ReadVariableNames',true);
    Dates    = rawdata.Dates;
    
    
    %% Collect in-sample quantiles
    if opt.hh ==1
        IS_modelname = [datafilename '_' sheetuse '_' '1973M1' '2020M9' ...
            '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];
    elseif opt.hh ==12
        IS_modelname = [datafilename '_' sheetuse '_' '1973M1' '2019M10' ...
            '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];
    end
    
    IS_results = load([IS_paramsfolder IS_modelname '_results.mat']);
    IS_results = IS_results.results;
    
    %% Collect data for plots
    
    end_date = datestr(esample,dataformat);
    
    % Collect Full Database
    [db, db_full,startdb,enddb,tex] = fLoadData(datafilename,sheetuse,start_date,end_date,opt.hh);
    
    % Collect MF, FF and trend (estimation sample)
    FF    = db.FF.data(nlags+1:end);
    MF    = db.MF.data(nlags+1:end);
    TR    = db.TRENDH.data(nlags+1:end);
    GDPG  = db.GDPG.data(nlags+1:end);
    GDPGH = db.GDPGH.data(nlags+1:end);
    
    % Collect MF, FF and trend (full sample)
    FF_full = db_full.FF.data;
    MF_full = db_full.MF.data;
    TR_full = db_full.TRENDH.data;
    GDPG_full  = db_full.GDPG.data;
    GDPGH_full = db_full.GDPGH.data;
    GDPGH_full_wt = GDPGH_full + TR_full;
    
    % Complete the trend with the last observation
    for jj=1:length(TR_full)
        if TR_full(jj)==0; TR_full(jj) = TR_full(jj-1); end
    end
    
    % Vector of dates for the full sample and for the available estimation sample
    dates_full = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';
    
    % Index of dates for plotting
    dates = datenum((datetime(efirstok,'InputFormat',inputformat)+calmonths(nlags)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';
    start_plot = efirstok;
    end_plot   = '2018-Oct';
    sd      = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates);
    ed      = find(datenum((datetime(end_plot,'InputFormat',inputformat)))==dates);
    sd_full = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates_full);
    ed_full = find(datenum((datetime(end_plot,'InputFormat',inputformat)))==dates_full);
    
    
    
    %% Blend OOS results
    
    for i=1:maxwin
        
        % Information for loading stored estimation results
        end_date = datestr(datenum(datetime(efirstok,'InputFormat',inputformat)+calmonths(i)),dataformat);
        startdb = [num2str(year(datetime(start_date,'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(start_date,'InputFormat','yyyy-MMM')))];
        e=find(Dates==datetime(end_date,'InputFormat','yyyy-MMM'));
        enddb   = [num2str(year(datetime(Dates(e-(opt.hh)),'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(Dates(e-(opt.hh)),'InputFormat','yyyy-MMM')))];
        myindex = find(datenum((datetime(end_date,'InputFormat',inputformat)))==dates_full);
        
        % Model name and title
        % Note: iterated model posterior mode depends on opt.hh through the estimation sample
        modelname = [datafilename '_' sheetuse '_' startdb enddb ...
            '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];
        
        x = load([paramsfolder modelname '.mat']);
        results.sv = x.results.sv;
               
        %% Collect parameters
        sv = results.sv;
        % Posterior mode values
        pmode=posterior_mode(sv);

        if strcmp(opt.paramsUse,'mode')
            params_in   = pmode;
            opt.nParamDraws = 1; %length(pmode.a12);
        elseif strcmp(opt.paramsUse,'mcmc')
            params_in   = pmcmc;
            opt.nDraws  = 10;
        end

        % Number of time-periods in simulation of GDP = number of periods
        % from start_date - end_date. Recall end_date is expanding
        opt.tperiods = length(dates_full);
        fprintf('\n Model estimated through end_date: %s \n',end_date);

        % Define target periods for density cuts
        target_dates = datenum(datetime(end_date,'InputFormat',inputformat));

        % Extract the index for the desired density plots
        opt.date_index = find(ismember(dates_full, target_dates));

        % Simulate model
        if opt.dir_or_it==1
            % Direct forecast
            [quantiles,pdensity] = fPredDensityDirectFull(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt,simtype);
        elseif opt.dir_or_it==2
            % Iterated forecast
            [quantiles, pdensity] = fPredDensityIteratedFull(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,simtype);
        end

        %==================
        % COLLECT RESULTS
        %==================
        
        % Parameters
        results.outparams = scriptParams_FULL_companion(sv,params_in,1);

        % Quantiles
        results.quantiles.dYsim_25  = quantiles.dYsim_25(opt.date_index);
        results.quantiles.dYsim_75  = quantiles.dYsim_75(opt.date_index);
        results.quantiles.dYsim_10  = quantiles.dYsim_10(opt.date_index);
        results.quantiles.dYsim_90  = quantiles.dYsim_90(opt.date_index);
        results.quantiles.mean      = quantiles.mean(opt.date_index);
        results.quantiles.st_t_mean = quantiles.st_t_mean(opt.date_index);

        % Predictive density at specified episodes
        [results.pdfi,results.xi]  = ksdensity(pdensity.realized);
        [results.pdfi_good,results.xi_good]  = ksdensity(pdensity.good);
        [results.pdfi_bad,results.xi_bad]  = ksdensity(pdensity.bad);
        [results.cxi,results.cdfi] = ksdensity(pdensity.realized,'Function','icdf');
        
        % Store results corresponding to period-i
        benchmark.pdfi{i}      = results.pdfi;                    % PDF ergodic
        benchmark.pdfi_good{i} = results.pdfi_good;               % PDF good regime
        benchmark.pdfi_bad{i}  = results.pdfi_bad;                % PDF bad regime
        benchmark.xi{i}        = results.xi;                      % support ofPDF ergodic
        benchmark.xi_good{i}   = results.xi_good;                 % support ofPDF good regime
        benchmark.xi_bad{i}    = results.xi_bad;                  % support ofPDF bad regime
        benchmark.cdfi{i}      = results.cdfi;                    % CDF ergodic
        benchmark.cxi{i}       = results.cxi;                     % support ofCDF ergodic
        benchmark.dYsim_25(i,1)  = results.quantiles.dYsim_25;    % 25th quantile
        benchmark.dYsim_75(i,1)  = results.quantiles.dYsim_75;    % 75th quantile
        benchmark.dYsim_10(i,1)  = results.quantiles.dYsim_10;    % 10th quantile
        benchmark.dYsim_90(i,1)  = results.quantiles.dYsim_90;    % 90th quantile
        benchmark.qmean(i,1)      = results.quantiles.mean;        % average
        benchmark.st_t_mean(i,1) = results.quantiles.st_t_mean;   % transition probability

        %===============
        % Map parameters
        %===============
        %         outparams{i} = results.outparams;
        A0_sync_1(i,:,:) = results.outparams.A0_sync_1;
        A0_sync_2(i,:,:) = results.outparams.A0_sync_2;
        C_sync_1(i,:,:) = results.outparams.C_sync_1;
        C_sync_2(i,:,:) = results.outparams.C_sync_2;
        A1_sync_1(i,:,:) = results.outparams.A1_sync_1;
        A1_sync_2(i,:,:) = results.outparams.A1_sync_2;
        SIG_sync_1(i,:,:) = results.outparams.SIG_sync_1;
        SIG_sync_2(i,:,:) = results.outparams.SIG_sync_2;
        D_sync_1(i,:,:) = results.outparams.D_sync_1;
        D_sync_2(i,:,:) = results.outparams.D_sync_2;
        B_sync_1(i,:,:) = results.outparams.B_sync_1;
        B_sync_2(i,:,:) = results.outparams.B_sync_2;
        O_sync_1(i,:,:) = results.outparams.O_sync_1;
        O_sync_2(i,:,:) = results.outparams.O_sync_2;
        a12(i,:) = results.outparams.a12;
        a21(i,:) = results.outparams.a21;
        if opt.transprob==1
            b12(i,:) = results.outparams.b12;
            b21(i,:) = results.outparams.b21;
            c12(i,:) = results.outparams.c12;
            c21(i,:) = results.outparams.c21;
        end
        
        %         % Do not use last hh observations since no data for GDPG
        %         ps(i,1)        = results.ps;        % predictive score
        %         pits(i,1)      = results.pits;      % PIT
        
        %=======================================
        % Predictive score at specified episodes
        %=======================================
        [~,results.xi_ps]  = min(abs(results.xi-GDPGH_full_wt(myindex,:)));
        if isempty(results.xi_ps) || results.xi_ps==1 || results.xi_ps==100
            ps(i,1) = 0;
        else
            ps(i,1) = results.pdfi(1,results.xi_ps);
        end
        
        %===========================
        % PITs at specified episodes
        %===========================
        [~,results.xi_pits]  = min(abs(results.cxi-GDPGH_full_wt(myindex,:)));
        if isempty(results.xi_pits) || results.xi_pits==1 || results.xi_pits==100
            pits(i,1) = 0;
        else
            pits(i,1) = results.cdfi(1,results.xi_pits);
        end

        %==================================
        % Counterfactual densities
        %==================================
        if ~isempty(mycounterfactuals)
           
            for j=1:length(mycounterfactuals)
                opt_temp = opt;
                opt_temp.counterfactual = mycounterfactuals(j);

                 % Simulate model
                if opt.dir_or_it==1
                    % Direct forecast
                    [quantiles_cf,pdensity_cf] = fPredDensityDirectFull(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt_temp,simtype);
                elseif opt.dir_or_it==2
                    % Iterated forecast
                    [quantiles_cf, pdensity_cf] = fPredDensityIteratedFull(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt_temp,simtype);
                end

            % Quantiles
            cfresults.quantiles.dYsim_25  = quantiles_cf.dYsim_25(opt.date_index);
            cfresults.quantiles.dYsim_75  = quantiles_cf.dYsim_75(opt.date_index);
            cfresults.quantiles.dYsim_10  = quantiles_cf.dYsim_10(opt.date_index);
            cfresults.quantiles.dYsim_90  = quantiles_cf.dYsim_90(opt.date_index);
            cfresults.quantiles.mean      = quantiles_cf.mean(opt.date_index);
            cfresults.quantiles.st_t_mean = quantiles_cf.st_t_mean(opt.date_index);

            % Predictive density at specified episodes
            [cfresults.pdfi,cfresults.xi]           = ksdensity(pdensity_cf.realized);
            [cfresults.pdfi_good,cfresults.xi_good] = ksdensity(pdensity_cf.good);
            [cfresults.pdfi_bad,cfresults.xi_bad]   = ksdensity(pdensity_cf.bad);
            [cfresults.cxi,cfresults.cdfi]          = ksdensity(pdensity_cf.realized,'Function','icdf');
            

            % Store results corresponding to period-i and counterfactual-j
            cf{j}.pdfi{i}       = cfresults.pdfi;              % PDF ergodic
            cf{j}.pdfi_good{i}  = cfresults.pdfi_good;         % PDF good regime
            cf{j}.pdfi_bad{i}   = cfresults.pdfi_bad;          % PDF bad regime
            cf{j}.xi{i}         = cfresults.xi;                % support ofPDF ergodic
            cf{j}.xi_good{i}    = cfresults.xi_good;           % support ofPDF good regime
            cf{j}.xi_bad{i}     = cfresults.xi_bad;            % support ofPDF bad regime
            cf{j}.cdfi{i}       = cfresults.cdfi;              % CDF ergodic
            cf{j}.cxi{i}        = cfresults.cxi;               % support ofCDF ergodic
            cf{j}.dYsim_25(i,1) = cfresults.quantiles.dYsim_25;    % 25th quantile
            cf{j}.dYsim_75(i,1) = cfresults.quantiles.dYsim_75;    % 75th quantile
            cf{j}.dYsim_10(i,1) = cfresults.quantiles.dYsim_10;    % 10th quantile
            cf{j}.dYsim_90(i,1) = cfresults.quantiles.dYsim_90;    % 90th quantile
            cf{j}.qmean(i,1)    = cfresults.quantiles.mean;        % average
            cf{j}.st_t_mean(i,1)= cfresults.quantiles.st_t_mean;   % transition probability
            end
        end
        
    end
       
    
    %----------------------------------------------------------
    %% Plots
    %----------------------------------------------------------
    
    if opt.showfig==1
        % Figure options
        figureTitleTag = [strrep(spectag,'_',' ') ', ' opt.paramsUse ', ' strrep(sheetuse,'_',' ') ', h=' num2str(opt.hh) ', ' startdb ':' enddb];
        FontSize  = 16;
        numticks  = 48;
        figSize   = [12 6];
        linestyle = {'-','--',':'};
        colors    = cbrewer('div', 'RdYlBu', 64);
        colors2   = cbrewer('qual', 'Set1', 8);
        left_color= [0 0 0];
        right_color= colors(2,:);
        
        %% Figure 1 - Quantiles and Transition Probabilities
        if opt.selfig(1)
            fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            subplot(211)
            hold on
            l1=plot(dates(sd:ed), benchmark.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th out-of-sample');
            l2=plot(dates(sd:ed), benchmark.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th out-of-sample');
            l3=plot(dates_full(sd_full:ed_full), IS_results.quantiles.dYsim_10(sd_full:ed_full),'-.','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th in-sample');
            l4=plot(dates_full(sd_full:ed_full), IS_results.quantiles.dYsim_90(sd_full:ed_full),'-.','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th in-sample');
            legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            l5=plot(dates_full(sd_full:ed_full),GDPGH_full_wt(sd_full:ed_full),'-.','Color','k','LineWidth', 2,'DisplayName','Realized Value');
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            if opt.dir_or_it==1
                title(['Quantiles Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            else
                title(['Quantiles Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            end
            axis tight
            numticks = 12;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            if opt.hh==12
                ylim([-10 8])
            else
                ylim([-30 20])
            end
            
            subplot(212)
            hold on
            l3=plot(dates(sd:ed), 1-benchmark.st_t_mean(sd:ed),'-','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Simulated $s(t)$');
            l4=plot(dates_full(sd_full:ed_full), 1-IS_results.quantiles.st_t_mean(sd_full:ed_full),'-.','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Simulated s(t) in-sample');
            legend([l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
            axis tight
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            ylim([0 1])
            
            % Set figure appereance
            set(fig1,'PaperOrientation','portrait');
            set(fig1, 'PaperSize', figSize);
            set(fig1, 'PaperUnits', 'inches');
            set(fig1, 'Units','inches');
            set(fig1, 'PaperPositionMode', 'auto');
            set(fig1, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;
        end
        
        %% Figure 2 - A0 matrix
        if opt.selfig(2)
            fig2=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            lw = 2;
            fs = 12;
            
            subplot(3,1,1)
            g = A0_sync_1(:,2,1);
            b = A0_sync_2(:,2,1);
            plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
            hold on
            plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
            hold off
            legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
            title('$A_O(2,1)$','FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            
            subplot(3,1,2)
            g = A0_sync_1(:,3,1);
            b = A0_sync_2(:,3,1);
            plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
            hold on
            plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
            hold off
            legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
            title('$A_O(3,1)$','FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            
            subplot(3,1,3)
            g = A0_sync_1(:,3,2);
            b = A0_sync_2(:,3,2);
            plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
            hold on
            plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
            hold off
            legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
            title('$A_O(3,2)$','FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            
            % Set figure appereance
            set(fig2,'PaperOrientation','portrait');
            set(fig2, 'PaperSize', figSize);
            set(fig2, 'PaperUnits', 'inches');
            set(fig2, 'Units','inches');
            set(fig2, 'PaperPositionMode', 'auto');
            set(fig2, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;
        end
        
        %% Figure 3 - A1 matrix
        if opt.selfig(3)
            fig3=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            lw = 2;
            fs = 12;
            
            count = 0;
            for k=1:size(A1_sync_1,2)
                for l=1:size(A1_sync_1,3)
                    count = count + 1;
                    subplot(size(A1_sync_1,2),size(A1_sync_1,3),count)
                    g = A1_sync_1(:,k,l);
                    b = A1_sync_2(:,k,l);
                    plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
                    hold on
                    plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
                    hold off
                    legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
                    title(['$A_1(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
                    axis tight
                    datetick('x','yyyy','keepticks')
                    set(gca, 'XLim', [dates(sd), dates(ed)])
                end
            end
            % Set figure appereance
            set(fig3,'PaperOrientation','portrait');
            set(fig3, 'PaperSize', figSize);
            set(fig3, 'PaperUnits', 'inches');
            set(fig3, 'Units','inches');
            set(fig3, 'PaperPositionMode', 'auto');
            set(fig3, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;
        end
        
        %% Figure 4 - C matrix
        if opt.selfig(4)
            fig4=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            lw = 2;
            fs = 12;
            
            count = 0;
            for k=1:size(C_sync_1,2)
                count = count + 1;
                subplot(size(C_sync_1,2),1,count)
                g = C_sync_1(:,k,1);
                b = C_sync_2(:,k,1);
                plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
                hold on
                plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
                hold off
                legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
                title(['$C(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
                axis tight
                datetick('x','yyyy','keepticks')
                set(gca, 'XLim', [dates(sd), dates(ed)])
            end
            % Set figure appereance
            set(fig4,'PaperOrientation','portrait');
            set(fig4, 'PaperSize', figSize);
            set(fig4, 'PaperUnits', 'inches');
            set(fig4, 'Units','inches');
            set(fig4, 'PaperPositionMode', 'auto');
            set(fig4, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;
        end
        
        %% Figure 5 - SIG matrix
        if opt.selfig(5)
            fig5=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            lw = 2;
            fs = 12;
            
            count = 0;
            for k=1:size(SIG_sync_1,2)
                count = count + 1;
                subplot(size(SIG_sync_1,2),1,count)
                g = SIG_sync_1(:,k,k);
                b = SIG_sync_2(:,k,k);
                plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
                hold on
                plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
                hold off
                legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
                title(['$\Sigma(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
                axis tight
                datetick('x','yyyy','keepticks')
                set(gca, 'XLim', [dates(sd), dates(ed)])
            end
            % Set figure appereance
            set(fig5,'PaperOrientation','portrait');
            set(fig5, 'PaperSize', figSize);
            set(fig5, 'PaperUnits', 'inches');
            set(fig5, 'Units','inches');
            set(fig5, 'PaperPositionMode', 'auto');
            set(fig5, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;
        end
        
        %% Figure 6 - transition probability parameters
        if opt.selfig(6)
            fig6=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            lw = 2;
            fs = 12;
            
            subplot(3,1,1)
            g = a12;
            b = a21;
            plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
            hold on
            plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
            hold off
            legend('good to bad','bad to good','box','off','FontSize',fs-2,'Interpreter','Latex')
            title(['$a(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            
            subplot(3,1,2)
            g = b12;
            b = b21;
            plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
            hold on
            plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
            hold off
            legend('good to bad','bad to good','box','off','FontSize',fs-2,'Interpreter','Latex')
            title(['$b(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            
            subplot(3,1,3)
            g = c12;
            b = c21;
            plot(dates(sd:ed),g(sd:ed),'LineWidth', lw)
            hold on
            plot(dates(sd:ed),b(sd:ed),'LineWidth', lw)
            hold off
            legend('good to bad','bad to good','box','off','FontSize',fs-2,'Interpreter','Latex')
            title(['$c(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            
            % Set figure appereance
            set(fig6,'PaperOrientation','portrait');
            set(fig6, 'PaperSize', figSize);
            set(fig6, 'PaperUnits', 'inches');
            set(fig6, 'Units','inches');
            set(fig6, 'PaperPositionMode', 'auto');
            set(fig6, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;
        end
        
        %% Figure 7 - Predictive Score
        if opt.selfig(7)
            fig7=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            lw = 2;
            fs = 12;
            
            % Do not use last hh observations since no data for GDPG
            plot(dates(sd:(ed-opt.hh)),ps(sd:(ed-opt.hh)),'LineWidth', lw)
            hold on
            plot(dates_full(sd_full:(ed_full-opt.hh)),IS_results.ps(sd_full:(ed_full-opt.hh)),'-.','LineWidth', lw)
            hold off
            legend('OOS','IS','Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            if opt.dir_or_it==1
                title(['Predictive Score Direct: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            else
                title(['Predictive Score Iterated: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            end
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            ylim([0 1])
            
            % Set figure appereance
            set(fig7,'PaperOrientation','portrait');
            set(fig7, 'PaperSize', figSize);
            set(fig7, 'PaperUnits', 'inches');
            set(fig7, 'Units','inches');
            set(fig7, 'PaperPositionMode', 'auto');
            set(fig7, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;
        end
        
        %% Figure 8 - PITs
        if opt.selfig(8)
            fig8=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            lw = 2;
            fs = 12;
            
            % Do not use last hh observations since no data for GDPG
            plot(dates(sd:(ed-opt.hh)),pits(sd:(ed-opt.hh)),'LineWidth', lw)
            hold on
            plot(dates_full(sd_full:(ed_full-opt.hh)),IS_results.pits(sd_full:(ed_full-opt.hh)),'-.','LineWidth', lw)
            hold off
            legend('OOS','IS','Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            if opt.dir_or_it==1
                title(['PITs Direct: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            else
                title(['PITs Iterated: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            end
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(ed)])
            ylim([0 1])
            
            % Set figure appereance
            set(fig8,'PaperOrientation','portrait');
            set(fig8, 'PaperSize', figSize);
            set(fig8, 'PaperUnits', 'inches');
            set(fig8, 'Units','inches');
            set(fig8, 'PaperPositionMode', 'auto');
            set(fig8, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;
        end
        
        
        %% Figure 9 - CDF of PITs and Rossi-Sekhposyan Test
        if opt.selfig(9)
            fig9=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            hold on
            lw = 2;
            fs = 10;
            figSize9 = [10 8]/1.5;
            
            % Do not use last hh observations since no data for GDPG
            pits = pits(sd:(ed-opt.hh));
            P = length(pits);
            el = round(P^(1/3));
            bootMC = 300;
            rvec = (0:0.001:1);
            
            % Critical Values Rossi-Sekhposyan Test
            if opt.hh == 1
                results1CS = rs_teststat(pits);
                result = results1CS(1,1);
            elseif opt.hh > 1
                table1 = CVfinalbootstrapInoue(el,bootMC,pits,rvec);
            end
            
            % Plot Uniformity Test
            h = cdfplot(pits);
            set(h,'Color','b','LineWidth',2.5);
            if opt.hh == 1
                plot(rvec,rvec + 1.34/sqrt(P),'b:','LineWidth',2);
                plot(rvec,rvec - 1.34/sqrt(P),'b:','LineWidth',2,'HandleVisibility','off');
            elseif opt.hh > 1
                plot(rvec,rvec + table1(2,1)/sqrt(P),'b:','LineWidth',2);
                plot(rvec,rvec - table1(2,1)/sqrt(P),'b:','LineWidth',2,'HandleVisibility','off');
            end
            plot(rvec,rvec,'r--','LineWidth',2);
            hL = legend('Empirical CDF','5\% Critical Value','Uniform CDF');
            set(hL,'interpreter','Latex','Orientation','Vertical','Location','SouthEast')
            legend boxoff
            xlim([0 1])
            ylim([0 1])
            if opt.dir_or_it==1
                title(['Direct: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            else
                title(['Iterated: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            end
            ylabel('CDF of PITs','interpreter','Latex')
            xlabel('$\tau$','interpreter','Latex')
            set(gca, 'FontName', 'Times New Roman');
            set(gca, 'FontSize', fs);
            set(gca,'Layer','top')
            hold off
            % Set figure appereance
            set(fig9,'PaperOrientation','portrait');
            set(fig9, 'PaperSize', figSize9);
            set(fig9, 'PaperUnits', 'inches');
            set(fig9, 'Units','inches');
            set(fig9, 'PaperPositionMode', 'auto');
            set(fig9, 'Position', [figSize9(1)/2 figSize9(2)/2 figSize9(1) figSize9(2)]);
            tightfig;
        end

        %% Figure 10 - Compare predictive density to counterfacual predictive density
        if opt.selfig(10)
            fig10=figure;clf;
            subplot(2,1,1)
            l1=plot(dates(sd:end),benchmark.dYsim_10(sd:end),'-','LineWidth',2,'Color',colors(10,:),'DisplayName','10th out-of-sample'); hold on;
            l2=plot(dates(sd:end),benchmark.dYsim_90(sd:end),'-','LineWidth',2,'Color',colors(50,:),'DisplayName','90th out-of-sample');
            l3=plot(dates(sd:end),cf{1}.dYsim_10(sd:end),'--','LineWidth',2,'Color',colors(10,:),'DisplayName','Couterfactual 1');
            l4=plot(dates(sd:end),cf{1}.dYsim_90(sd:end),'--','LineWidth',2,'Color',colors(50,:),'DisplayName','Couterfactual 1');
            legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            if opt.dir_or_it==1
                title(['OOS Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            else
                title(['OOS Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            end
            axis tight
            numticks = 36;
            set(gca,'XTick',datenum(dates(sd:numticks:(end))))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates(end)])

            subplot(2,1,2)
            l1=plot(dates(sd:end),benchmark.dYsim_10(sd:end),'-','LineWidth',2,'Color',colors(10,:),'DisplayName','10th out-of-sample'); hold on;
            l2=plot(dates(sd:end),benchmark.dYsim_90(sd:end),'-','LineWidth',2,'Color',colors(50,:),'DisplayName','90th out-of-sample');
            l3=plot(dates(sd:end),cf{2}.dYsim_10(sd:end),'--','LineWidth',2,'Color',colors(10,:),'DisplayName','Couterfactual 2');
            l4=plot(dates(sd:end),cf{2}.dYsim_90(sd:end),'--','LineWidth',2,'Color',colors(50,:),'DisplayName','Couterfactual 2');
            legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            if opt.dir_or_it==1
                title(['OOS Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            else
                title(['OOS Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            end
            axis tight
            numticks = 36;
            set(gca,'XTick',datenum(dates(sd:numticks:(end))))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(sd), dates((end))])

            % Set figure appereance
            set(fig10,'PaperOrientation','portrait');
            set(fig10, 'PaperSize', figSize);
            set(fig10, 'PaperUnits', 'inches');
            set(fig10, 'Units','inches');
            set(fig10, 'PaperPositionMode', 'auto');
            set(fig10, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;

        end

        %% Figure 11 - Compare predictive density to counterfacual predictive density
        if opt.selfig(11)
            count = 1;
            fig11 = figure;clf; 
            for ix = [92 99 242 244]
                subplot(2,2,count)
                l1=plot(benchmark.xi{ix},benchmark.pdfi{ix},'k','LineWidth',2,'DisplayName','Estimated'); hold on;
                l2=plot(cf{1}.xi{ix},cf{1}.pdfi{ix},'--','LineWidth',2,'Color',colors(10*1,:),'DisplayName',['No Macro Factor']);
                l3=plot(cf{2}.xi{ix},cf{2}.pdfi{ix},'--','LineWidth',2,'Color',colors(10*2,:),'DisplayName',['No Financial Factor']);
                title(datestr(datenum(datetime(efirstok,'InputFormat',inputformat)+calmonths(ix)),dataformat),'Interpreter','Latex')
                legend([l1 l2 l3],'Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
                set(gca,'children',flipud(get(gca,'children'))); xlim([-12 6]);
                count = count + 1;
            end
            axes( 'Position', [0, 0.95, 1, 0.05] ) ;
            set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
            %     if opt.dir_or_it==1        
            %         text( 0.5, 0.1, ['OOS Direct: ' figureTitleTag ], 'FontSize', 14', 'FontWeight', 'Bold', ...
            %             'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
            %     else
            %         text( 0.5, 0.2, ['OOS Iterated: ' figureTitleTag ], 'FontSize', 14', 'FontWeight', 'Bold', ...
            %             'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','Interpreter','Latex' ) ;        
            %     end
            % Set figure appereance
            set(fig11,'PaperOrientation','portrait');
            set(fig11, 'PaperSize', figSize);
            set(fig11, 'PaperUnits', 'inches');
            set(fig11, 'Units','inches');
            set(fig11, 'PaperPositionMode', 'auto');
            set(fig11, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
            tightfig;

        end

        %% Table 1
        
        % Do not consider last 12 observations since no data for GDPG
        stot = length(sd:(ed-opt.hh));
        snum_10 = sum(GDPGH_full_wt(sd_full:(ed_full-opt.hh))<benchmark.dYsim_10(sd:(ed-opt.hh)));
        vratio_10 = snum_10/stot*100;
        snum_25 = sum(GDPGH_full_wt(sd_full:(ed_full-opt.hh))<benchmark.dYsim_25(sd:(ed-opt.hh)));
        vratio_25 = snum_25/stot*100;
        snum_75 = sum(GDPGH_full_wt(sd_full:(ed_full-opt.hh))<benchmark.dYsim_75(sd:(ed-opt.hh)));
        vratio_75 = snum_75/stot*100;
        snum_90 = sum(GDPGH_full_wt(sd_full:(ed_full-opt.hh))<benchmark.dYsim_90(sd:(ed-opt.hh)));
        vratio_90 = snum_90/stot*100;
        
        snum_IS_10 = sum(GDPGH_full_wt(sd_full:(ed_full-opt.hh))<IS_results.quantiles.dYsim_10(sd_full:(ed_full-opt.hh)));
        vratio_IS_10 = snum_IS_10/stot*100;
        snum_IS_25 = sum(GDPGH_full_wt(sd_full:(ed_full-opt.hh))<IS_results.quantiles.dYsim_25(sd_full:(ed_full-opt.hh)));
        vratio_IS_25 = snum_IS_25/stot*100;
        snum_IS_75 = sum(GDPGH_full_wt(sd_full:(ed_full-opt.hh))<IS_results.quantiles.dYsim_75(sd_full:(ed_full-opt.hh)));
        vratio_IS_75 = snum_IS_75/stot*100;
        snum_IS_90 = sum(GDPGH_full_wt(sd_full:(ed_full-opt.hh))<IS_results.quantiles.dYsim_90(sd_full:(ed_full-opt.hh)));
        vratio_IS_90 = snum_IS_90/stot*100;
        
        % Display table
        dc = 1;
        T = table([round(vratio_10,dc);round(vratio_IS_10,dc)],...
            [round(vratio_25,dc);round(vratio_IS_25,dc)],...
            [round(vratio_75,dc);round(vratio_IS_75,dc)],...
            [round(vratio_90,dc);round(vratio_IS_90,dc)],'VariableNames',{'10th','25th','75th','95th'},'RowName',{'OOS','IS'});
        T = table(T,'VariableNames',{'Violation Ratios'});
        disp(T)
        
        %% Save Figures
        % 1 = quantiles predictive distribution
        % 2 = A0 parameters time evolution
        % 3 = A1 parameters time evolution
        % 4 = C parameters time evolution
        % 5 = SIG parameters time evolution
        % 6 = transition probability parameters time evolution
        % 7 = predictive score
        % 8 = PITs
        % 9 = CDF of PITs
        modelname = [datafilename '_' sheetuse '_' startdb enddb ...
            '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];

        if opt.saveit
            if opt.selfig(1)
                print('-dpdf',fig1,[figfolder 'Quantiles_' modelname ],'-bestfit');
            end
            if opt.selfig(2)
                print('-dpdf',fig2,[figfolder 'A0_' datestr(dates_full(opt.date_index(ii)),dataformat) '_' modelname ],'-bestfit');
            end
            if opt.selfig(3)
                print('-dpdf',fig3,[figfolder 'A1_' modelname ],'-bestfit');
            end
            if opt.selfig(4)
                print('-dpdf',fig4,[figfolder 'C_' modelname ],'-bestfit');
            end
            if opt.selfig(5)
                print('-dpdf',fig5,[figfolder 'SIG_' modelname ],'-bestfit');
            end
            if opt.selfig(6)
                print('-dpdf',fig6,[figfolder 'TransProb_' modelname ],'-bestfit');
            end
            if opt.selfig(7)
                print('-dpdf',fig7,[figfolder 'PredScore_' modelname ],'-bestfit');
            end
            if opt.selfig(8)
                print('-dpdf',fig8,[figfolder 'PITs_' modelname ],'-bestfit');
            end
            if opt.selfig(9)
                print('-dpdf',fig9,[figfolder 'PITsCDF_' modelname ],'-bestfit');
            end
        end
        
    end
end



%%
rise_exit
