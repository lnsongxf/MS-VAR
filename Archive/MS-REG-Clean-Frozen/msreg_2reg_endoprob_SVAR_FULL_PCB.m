% Understanding Growth-at-Risk: A Markov-Switching Approach
% Model with f_t and m_t endogenous
%
% This script produces posterior mode estimates and predictive
% densities at the posterior mode for any dataset.
%
% There are two estimation optionns
%  1. Direct: estimates the measurement equation at horizon hh.
%  2. Iterated: estimates the contemporanous measurement equation and
%     then iterates the predictive density through horizon hh.
%=====================================================================

%% housekeeping
clear; close all; clc; tic;

% Important paths
addpath('../RISE_toolbox');         % RISE Toolbox
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
opt.dir_or_it      = 2;     % 1 = direct, 2 = iterated
opt.hh             = 12;    % forecast horizon
opt.const          = 1;     % 1 = have a constant in transition probability
opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
opt.paramsUse      = 'mode';
opt.force_estimate = 1; % 1 = force to reestimate model

% Posterior mode options
opt.nDraws      = 5000;  % Number of simulations for eacht time period

% MCMC simulation options
opt.nParamDraws  = 1000;  % Draws from posterior. has to be less thatn options.N;

% Figure options
opt.saveit       = 1;     % 1 = save figures

% Select which models to run
% mymodels = [3,4,6,7];
mymodels = 4;

if opt.dir_or_it==1
    tag = [];
else
    tag = ['A1restricted'];
end

for imodel=mymodels
    
    modelspec = imodel;
    
    % Model spectag
    spectag = fSpecLabel(modelspec);
        
    % Data vintage, sample and country selection
    datafilename = '11302020';
    sheetuse     = 'US_DFM_TEST';
    start_date   = '1973-Jan';
    end_date     = '2020-Oct';
    
    % Define target periods for density cuts
    target_dates = [...%datenum('2020-Jan',dataformat),...
        %%datenum('2020-Feb',dataformat),...
        datenum('2007-Sep',dataformat)
        datenum('2008-Sep',dataformat)
        datenum('2020-Mar',dataformat)
        datenum('2020-Apr',dataformat)
        ];
    
    % Define dates for plots
    start_plot    = '1973-Feb';
    end_plot      = '2019-Oct';
    end_plot_full = '2020-Oct';
    
    % VAR configuration
    nlags=6;
    exog={};
    constant=true;
    panel=[];
    
    % MCMC options
    options=struct();
    options.alpha=0.234;
    options.thin= 10;
    options.burnin=10^3;
    options.N = 2*10^4;
    options.nchain=1;
    
    
    %========================================
    % MODEL SPECIFIC VARIABLES
    %========================================
    
    if opt.dir_or_it ==1
        fcstype= 'direct';
        varlist={'FF','MF','GDPGH'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/'];
        %texfolder   = ['results/Params/' opt.paramsUse '/Direct/Tex/'];
        figfolder   = ['results/Figures/' opt.paramsUse '/Full/Direct/' spectag '/'];
        
    elseif opt.dir_or_it ==2
        fcstype= 'iterated';
        varlist={'FF','MF','GDPG'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/'];
        %texfolder   = ['results/Params/' opt.paramsUse '/Iterated/Tex/'];
        figfolder   =  ['results/Figures/' opt.paramsUse '/Full/Iterated/' spectag '/'];
    end
    
    %========================================
    % CONFIGURE DATES AND PLOT OPTIONS
    %========================================
    
    % Vector of dates for the full sample
    dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';
    
    % Index of dates
    sd      = find(datenum(start_plot,dataformat)==dates_full);
    ed      = find(datenum(end_plot,dataformat)==dates_full);
    ed_full = find(datenum(end_plot_full,dataformat)==dates_full);
    
    % Vector of dates for plotting
    dates    = dates_full(sd:ed);
    opt.tperiods = length(dates_full);  % Number of time-periods in simulation of GDP
    
    % Figure options
    FontSize  = 16;
    numticks  = 48;
    figSize   = [12 6];
    linestyle = {'-','--',':'};
    colors    = cbrewer('div', 'RdYlBu', 64);
    colors2   = cbrewer('qual', 'Set1', 8);
    left_color= [0 0 0];
    right_color= colors(2,:);
    
    % Create folders to store results
    if exist(paramsfolder,'dir')==0;  mkdir(paramsfolder); end
    %if exist(texfolder,'dir')==0;  mkdir(texfolder); end
    if exist(figfolder,'dir')==0;  mkdir(figfolder); end
    
    % Diary file
    dfile = ([figfolder sheetuse '_' spectag '_params.txt']);
    if exist(dfile, 'file')
        delete(dfile);
    end
    
    %==========================================================================
    % LOAD DATA
    %==========================================================================
    
    
    [db, db_full,startdb,enddb,tex] = fLoadData(datafilename,sheetuse,start_date,end_date,opt.hh);
    
    % Collect MF, FF and trend (estimation sample)
    FF    = db.FF.data(nlags+1:end);
    MF    = db.MF.data(nlags+1:end);
    TR    = db.TRENDH.data(nlags+1:end);
    GDPG  = db.GDPG.data(nlags+1:end);
    GDPGH = db.GDPGH.data(nlags+1:end);
    
    % Collect MF, FF and trend (full sample)
    FF_full = db_full.FF.data(nlags+1:end);
    MF_full = db_full.MF.data(nlags+1:end);
    TR_full = db_full.TRENDH.data(nlags+1:end);
    GDPG_full  = db_full.GDPG.data(nlags+1:end);
    GDPGH_full = db_full.GDPGH.data(nlags+1:end);
    
    % Model name and title for posterior mode table
    % Note: iterated model posterior mode depends on opt.hh through the estimation sample
    
    modelname = [datafilename '_' sheetuse '_' startdb enddb ...
        '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];

    
    
    diary(dfile)
    fprintf('\n Model: %s \n', modelname)
    
    %==========================================================================
    % RISE ESTIMATION
    %==========================================================================
    
    if strcmp(opt.paramsUse,'mode')
        
        % Load results or estimate posterior mode
        if exist([paramsfolder modelname '.mat'],'file')==0 || opt.force_estimate==1
            % Estimate posterior mode
            run scriptEstimateMode_FULL.m;
            
            % Store solution structure
            save([paramsfolder modelname '.mat'],'sPMode')
        else
            % Load stored results
            load([paramsfolder modelname '.mat']);
            
            % Map objects:
            sv = sPMode.sv;
        end
        
    else
        
        % Load or estimate mcmc
        if exist([paramsfolder modelname '.mat'],'file')==0
            % Estimat posterior distribution
            run scriptEstimateMCMC.m
            
            % Store solution structure
            save([paramsfolder modelname '.mat'],'sPmcmc')
            
        else
            
            % Load stored results
            load([paramsfolder modelname '.mat']);
            
            % Map objects:
            sv = sPmcmc.sv;
            pmcmc = sPmcmc.params;
        end
        
    end
    
    %%
    %==========================================================================
    % Posterior Mode Analysis
    %==========================================================================
    
    % Posterior mode values
    pmode=posterior_mode(sv);
    
    % Map posterior mode coefficient estimates
    fnames = fieldnames(pmode);
    for jj=1:length(fnames)
        eval([fnames{jj} ' = pmode.' fnames{jj} ';']);
    end
    % Print Structural Form
    print_structural_form(sv)
    
    
    % Print Reduced Form
    print_solution(sv)
    
    %==========================================================================
    % Predictive Density at Posterior Mode or MCMC
    %==========================================================================
    
    if strcmp(opt.paramsUse,'mode')
        params_in   = pmode;
        opt.nParamDraws = 1; %length(pmode.a12);
    elseif strcmp(opt.paramsUse,'mcmc')
        params_in   = pmcmc;
        opt.nDraws      = 10;
    end

    %% Print A0, A1 and SIG matrices
    outparams = scriptParams_FULL(params_in,1,modelspec);
    
    %%
    fprintf('\n *** A0 Matrix: State = 1 *** \n')
    disp(outparams.A0_sync_1)    
    fprintf('\n *** A0 Matrix: State = 2 *** \n')
    disp(outparams.A0_sync_2)    
    
    fprintf('\n *** A1 Matrix: State = 1 *** \n')
    disp(outparams.A1_sync_1)    
    fprintf('\n *** A1 Matrix: State = 2 *** \n')
    disp(outparams.A1_sync_2)    
    

    fprintf('\n *** C Matrix: State = 1 *** \n')
    disp(outparams.C_sync_1)    
    fprintf('\n *** C Matrix: State = 2 *** \n')
    disp(outparams.C_sync_2)    

    fprintf('\n *** SIG Matrix: State = 1 *** \n')
    disp(outparams.SIG_sync_1)    
    fprintf('\n *** SIG Matrix: State = 2 *** \n')
    disp(outparams.SIG_sync_2)    
    %%
    % Extract the index for the desired density plots
    opt.date_index = find(ismember(dates_full, target_dates));
    
    if opt.dir_or_it==1
        % Direct forecast
        [quantiles,pdensity] = fPredDensityDirectFull(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt,modelspec);
        
    elseif opt.dir_or_it==2
        % Iterated forecast                
        [quantiles, pdensity] = fPredDensityIteratedFull(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,'forecastst',modelspec);
    end
    
    
    % Predictive density at specified episodes
    for ii=1:numel(opt.date_index)
        %[pdfi(ii,:),xi(ii,:)]  = ksdensity(pdensity(:,ii));
        [pdfi(ii,:),xi(ii,:)]  = ksdensity(pdensity.realized(:,ii));
        [pdfi_good(ii,:),xi_good(ii,:)]  = ksdensity(pdensity.good(:,ii));
        [pdfi_bad(ii,:),xi_bad(ii,:)]  = ksdensity(pdensity.bad(:,ii));
    end
    
    %==========================================================================
    % Collect Probabilities
    %==========================================================================
    [~,~,~,f]=filter(sv);
    % Smoothed probabilities
    p_reg1_smoothed = f.smoothed_regime_probabilities.regime1.data;
    p_reg2_smoothed = f.smoothed_regime_probabilities.regime2.data;
    % Filtered probabilities
    p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
    p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;
    if bad==1
        p_reg_bad_filtered = p_reg1_filtered;
        p_reg_bad_smoothed = p_reg1_smoothed;
    elseif bad==2
        p_reg_bad_filtered = p_reg2_filtered;
        p_reg_bad_smoothed = p_reg2_smoothed;
    end
    
    
    %%
    % ==========================================================================
    % Quantile Plots
    % ==========================================================================
    figureTitleTag = [strrep(spectag,'_',' ') ', ' opt.paramsUse ', ' strrep(sheetuse,'_',' ') ', h=' num2str(opt.hh) ', ' startdb ':' enddb];
    
    % Add additional objects to structure
    fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
    subplot(211)
    hold on
    l1=plot(dates_full(sd:ed), quantiles.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (hard coded)');
    l2=plot(dates_full(sd:ed), quantiles.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (hard coded)');
    legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
    set(gca,'children',flipud(get(gca,'children')))
    hold off
    ylabel('Percent','interpreter','Latex','fontsize',10)
    if opt.dir_or_it==1
        title(['Quantiles Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
    else
        title(['Quantiles Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
    end
    axis tight
    datetick('x','yyyy','keepticks')
    set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
    
    subplot(212)
    hold on
    l3=plot(dates_full(sd:ed), quantiles.st_t_mean(sd:ed),'-','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Simulated s(t)');
    l4=plot(dates_full(sd:ed), p_reg_bad_filtered(sd:ed),'--','Color',colors(20,:),'LineWidth', 2.5,'DisplayName','Filtered s(t)');
    legend([l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;    
    hold off
    ylabel('Percent','interpreter','Latex','fontsize',10)
    title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
    axis tight
    datetick('x','yyyy','keepticks')
    set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
    ylim([0 1])
    
    % Set figure appereance
    set(fig1,'PaperOrientation','portrait');
    set(fig1, 'PaperSize', figSize);
    set(fig1, 'PaperUnits', 'inches');
    set(fig1, 'Units','inches');
    set(fig1, 'PaperPositionMode', 'auto');
    set(fig1, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
    tightfig;
    
    
    %%
    %==========================================================================
    % Density plots
    %=======================================================;===================
    linetype= {'--','-','-.',':','-.'};
    for ii=1:numel(opt.date_index)
        
        fig2=figure;
        hold on
        plot(xi(ii,:),pdfi(ii,:),linetype{2},'Color',[0 0 0 ],'LineWidth', 3,'DisplayName',[datestr(dates_full(opt.date_index(ii)),dataformat)]); hold on;
        plot(xi_good(ii,:),pdfi_good(ii,:),linetype{4},'Color',colors(60,:),'LineWidth', 3,'DisplayName','Good regime'); hold on;
        plot(xi_bad(ii,:),pdfi_bad(ii,:),linetype{5},'Color',colors(1,:),'LineWidth', 3,'DisplayName','Bad regime'); hold on;
        if opt.dir_or_it==1
            title(['Direct: ' figureTitleTag ],'Interpreter','Latex');
        else
            title(['Iterated: ' figureTitleTag ],'Interpreter','Latex');
        end
        axis tight
        ylimits = ylim; xlimits = xlim;
%         text(xlimits(1)*0.865,0.8*ylimits(2),['Prob bad regime = ' num2str((1-quantiles.st_t_mean(opt.date_index(ii)))*100) '\%'],'FontSize',13,'Interpreter','Latex');
        plot(zeros(1,length(xi_bad(1,:))),zeros(1,length(pdfi_bad(1,:))),'*','DisplayName',['Prob bad regime = ' num2str(round((1-quantiles.st_t_mean(opt.date_index(ii)))*100,1)) '\%']);
        vline(0,'k--');
        text(GDPGH_full(opt.date_index(ii)),0.95*ylimits(2),'realized value','color','red','FontSize',9);
        legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
        ylabel('PDF','fontsize',10,'interpreter','Latex')
        xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
        set(gca, 'FontName', 'Times New Roman');
        set(gca, 'FontSize', FontSize);
        set(gca,'Layer','top')
        %set(gca,'XTick',unique(sort([0 round(GDPGH_full(opt.date_index(ii)),2) xlimits(1):2:xlimits(2)])))
        set(gca,'TickLabelInterpreter','Latex')
        axis tight
        %xlim([-18 8])
        set(fig2,'PaperOrientation','portrait');
        set(fig2, 'PaperSize', figSize);
        set(fig2, 'PaperUnits', 'inches');
        set(fig2, 'Units','inches');
        set(fig2, 'PaperPositionMode', 'auto');
        set(fig2, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
        tightfig;
        
        if opt.saveit
            print('-dpdf',fig2,[figfolder 'Densities' datestr(dates_full(opt.date_index(ii)),dataformat) '_' modelname ],'-bestfit');
            
        end
        
    end
    
    %% Regime specific fitted values
    [Resids,Fit]=residuals(sv);
    if opt.dir_or_it==1
        fit1 = Fit.GDPGH.data(:,1);
        fit2 = Fit.GDPGH.data(:,2);
    elseif opt.dir_or_it==2
        fit1 = Fit.GDPG.data(:,1);
        fit2 = Fit.GDPG.data(:,2);
    end

    fitgood = fit1;
    fitbad  = fit2;

    if bad==2
       fitgood = fit1;
       fitbad  = fit2;
    elseif bad==1
       fitgood = fit2;
       fitbad  = fit1;
     end
    
    fig3=figure;
    %yyaxis left
    hold on
    plot(fitgood,'Color',colors(1,:),'LineWidth', 1.5)
    plot(fitbad,'--','Color',colors(15,:),'LineWidth', 1.5)
    hold off
%%
    ylabel('Percent','interpreter','Latex','fontsize',10)
    leg = {'Good Regime','Bad Regime'};
    axis tight
    ax=gca;
    ax.XTick = datenum(dates(sd:numticks:ed));
    datetick('x','yyyy','keepticks')
    set(gca, 'XLim', [dates(sd), dates(ed)])
    legend(leg,'Orientation','Vertical','Location','South','interpreter','Latex')
    legend boxoff
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', FontSize);
    set(gca,'Layer','top')
    set(gca,'TickLabelInterpreter','Latex')
    hold off
    set(fig3,'PaperOrientation','portrait');
    set(fig3, 'PaperSize', figSize);
    set(fig3, 'PaperUnits', 'inches');
    set(fig3, 'Units','inches');
    set(fig3, 'PaperPositionMode', 'auto');
    set(fig3, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
    tightfig;
    
    %% Model Fit
    clearvars y_fit realized
    y_fit = fit1.*p_reg1_smoothed + fit2.*p_reg2_smoothed;
    
    realized = db.GDPG.data ;
    realized = realized(nlags+1:end);
    
    fig4=figure;
    hold on
    plot(dates(sd:ed),y_fit(sd:ed),'linewidth',2)
    plot(dates(sd:ed),realized(sd:ed),'k','linewidth',2)
    hold off
    ylabel('Percent','interpreter','Latex','fontsize',10)
    legend('Ergodic Fitted Value','Realized','interpreter','Latex')
    legend boxoff
    axis tight
    ax=gca;
    ax.XTick = datenum(dates(sd:numticks:ed));
    datetick('x','yyyy','keepticks')
    set(gca, 'XLim', [dates(sd), dates(ed)])
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', FontSize);
    set(gca,'Layer','top')
    set(gca,'TickLabelInterpreter','Latex')
    set(fig4,'PaperOrientation','portrait');
    set(fig4, 'PaperSize', figSize);
    set(fig4, 'PaperUnits', 'inches');
    set(fig4, 'Units','inches');
    set(fig4, 'PaperPositionMode', 'auto');
    set(fig4, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
    tightfig;
    
    
    %% Regime Probabilities
    fig5=figure;
    hold on
    plot(dates(sd:ed),p_reg_bad_smoothed(sd:ed),'linewidth',2)
    plot(dates(sd:ed),p_reg_bad_filtered(sd:ed),'k','linewidth',2)
    hold off
    ylabel('Probability','interpreter','Latex','fontsize',10)
    legend('Smoothed','Filtered','interpreter','Latex')
    legend boxoff
    axis tight
    ax=gca;
    ax.XTick = datenum(dates(sd:numticks:ed));
    datetick('x','yyyy','keepticks')
    set(gca, 'XLim', [dates(sd), dates(ed)])
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', FontSize);
    set(gca,'Layer','top')
    set(gca,'TickLabelInterpreter','Latex')
    set(fig5,'PaperOrientation','portrait');
    set(fig5, 'PaperSize', figSize);
    set(fig5, 'PaperUnits', 'inches');
    set(fig5, 'Units','inches');
    set(fig5, 'PaperPositionMode', 'auto');
    set(fig5, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
    tightfig;
    
    if opt.saveit
        print('-dpdf',fig1,[figfolder 'Quantiles_' modelname ],'-bestfit');
        print('-dpdf',fig3,[figfolder 'Fitted_' modelname ],'-bestfit');
        print('-dpdf',fig4,[figfolder 'Resids_' modelname ],'-bestfit');
        print('-dpdf',fig5,[figfolder 'ProbReg_' modelname ],'-bestfit');
    end
    
    diary off
    
end
%%
rise_exit
