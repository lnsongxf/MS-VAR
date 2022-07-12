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
clear; clc; tic;
close all; 

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
opt.dir_or_it      = 1;     % 1 = direct, 2 = iterated
opt.hh             = 12;    % forecast horizon
opt.const          = 1;     % 1 = have a constant in transition probability
opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
opt.paramsUse      = 'mode';
opt.force_estimate = 0;     % 1 = force to reestimate model
opt.mdd_estimate   = 0;     % 1 = force to estimate MDD if paramsUse=mcmc
opt.printsol       = 0;     % 1 = print model solution

% Simulation options
opt.nDraws       = 5000;    % Number of simulations for each time period
opt.nParamDraws  = 1000;    % Draws from posterior. has to be less than options.N;
opt.showfig      = 1;       % 1 = show figures
opt.saveit       = 1;       % 1 = save figures
opt.selfig       = [0 0 0 0 0 0 1 0 0];  % 1 = vector of figures

%Figure selection options
% 1 = quantiles predictive distribution
% 2 = densities
% 3 = regime specific fitted values
% 4 = ergodic fitted values
% 5 = regime probabilities
% 6 = predictive score over time
% 7 = PITs CDFs
% 8 = PITs over time
% 9 = Figures for slides (currently matching PhillyFed presentation)

% MCMC options
options = struct();
options.alpha = 0.234;
options.thin = 10;
options.burnin = 10^3;
options.N = 2*10^4;
% options.thin = 1;
% options.burnin = 0;
% options.N = 200;
options.nchain = 1;

% VAR configuration
% important: if nlags>1 the routine that pulls the parameters 
% (scriptParams_FULL_companion.m) does not work for MCMC
nlags=1; 
exog={};
constant=true;
panel=[];

% Select which models to run
% mymodels = [1 13];
% mymodels = [2 11];
% mymodels = [4 12];
mymodels = [12];

tag ='';

for imodel=mymodels
    
    modelspec = imodel;
    
    % Model spectag
    spectag = fSpecLabel(modelspec);
        
    % Data vintage, sample and country selection
    datafilename = '11302020';
    sheetuse     = 'US_DFM';
    start_date   = '1973-Jan';
    end_date     = '2020-Oct';
    
    % Define target periods for density cuts
    target_dates = [...%datenum('2020-Jan',dataformat),...
        %%datenum('2020-Feb',dataformat),...
%         datenum('2007-Aug',dataformat)
%         datenum('2007-Sep',dataformat)
%         datenum('2008-Sep',dataformat)
        datenum('2020-Mar',dataformat)
        datenum('2020-Apr',dataformat)
        ];
    
    % Define dates for plots using full dates vector
    start_plot    = start_date;
    end_plot      = '2019-Oct';

    % Define dates for plots using available estimation sample
    start_plot_var = '1974-Jan';
    end_plot_var   = '2019-Oct';

       
    %========================================
    % MODEL SPECIFIC VARIABLES
    %========================================
    
    if opt.dir_or_it ==1
        fcstype= 'direct';
        varlist={'FF','MF','GDPGH'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' ];
        texfolder   = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' spectag '/Tex/' ];
        logfolder   = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/Log/' ];
        figfolder   = ['results/Figures/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/'];
        figslides   =  ['results/Figures_Slides/' opt.paramsUse '/VAR' num2str(nlags) '/Direct/' spectag '/'];

    elseif opt.dir_or_it ==2
        fcstype= 'iterated';
        varlist={'FF','MF','GDPG'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/'];
        texfolder   = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' '/Tex/' ];
        logfolder   = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags)  '/' spectag '/' '/Log/' ];
        figfolder   =  ['results/Figures/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/'];
        figslides   =  ['results/Figures_Slides/' opt.paramsUse '/VAR' num2str(nlags) '/Iterated/' spectag '/'];

    end
    
    %========================================
    % CONFIGURE DATES AND PLOT OPTIONS
    %========================================
    
    % Vector of dates for the full sample and for the available estimation sample
    dates_full = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

    % Index of dates for plotting
    sd      = find(datenum((datetime(start_plot,'InputFormat',inputformat)))==dates_full);
    ed      = find(datenum(end_plot,dataformat)==dates_full);

    % Vector of dates for the available estimation sample
    sd_var = find(datenum((datetime(start_date,'InputFormat',inputformat)+calmonths(nlags)))==dates_full);
    ed_var  = find(datenum((datetime(end_date,'InputFormat',inputformat)-calmonths(opt.hh)))==dates_full);
    dates_var  = dates_full(sd_var:ed_var);

    % Vector of dates for plotting
    sd_varp  = find(datenum((datetime(start_plot_var,'InputFormat',inputformat)))==dates_var);
    ed_varp  = find(datenum((datetime(end_plot_var,'InputFormat',inputformat)))==dates_var);
    dates    = dates_var(sd_varp:ed_varp);

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
    if exist(texfolder,'dir')==0;  mkdir(texfolder); end
    if exist(logfolder,'dir')==0;  mkdir(logfolder); end
    if exist(figfolder,'dir')==0;  mkdir(figfolder); end
    if exist(figslides,'dir')==0;  mkdir(figslides); end

%%
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
    FF_full = db_full.FF.data;
    MF_full = db_full.MF.data;
    TR_full = db_full.TRENDH.data;
    GDPG_full  = db_full.GDPG.data;
    GDPGH_full = db_full.GDPGH.data;
    GDPGH_full_wt = GDPGH_full + TR_full;
    
    % Model name and title
    % Note: iterated model posterior mode depends on opt.hh through the estimation sample
    modelname = [datafilename '_' sheetuse '_' startdb enddb ...
        '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];    
 
    % Plot Data
%     if opt.showfig==1
%         fig0=figure;
%         plot(dates_full(sd:ed),[FF_full(sd:ed) MF_full(sd:ed) GDPG_full(sd:ed) GDPGH_full(sd:ed)],'LineWidth',1);
%         legend('FF','MF','GDPG','GDPGH')
%         axis tight
%         datetick('x','yyyy','keepticks')
%         set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
%         % Set figure appereance
%         set(fig0,'PaperOrientation','portrait');
%         set(fig0, 'PaperSize', figSize);
%         set(fig0, 'PaperUnits', 'inches');
%         set(fig0, 'Units','inches');
%         set(fig0, 'PaperPositionMode', 'auto');
%         set(fig0, 'Position', [figSize(1)/2 figSize(2)/5 figSize(1) figSize(2)]);
%         tightfig;
%     end
    
%%
    %==========================================================================
    % RISE ESTIMATION
    %==========================================================================
    
    % Diary file
    dfile = ([logfolder sheetuse '_' spectag '_params.txt']);
    if exist(dfile, 'file')
        delete(dfile);
    end
    
    diary(dfile)
    fprintf('\n Model: %s \n', modelname)
    
    %%%%%%%%
    % MODE %
    %%%%%%%%
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
        %%%%%%%%
        % MCMC %
        %%%%%%%%
       
        % Load or estimate mcmc
        if exist([paramsfolder modelname '.mat'],'file')==0 || opt.force_estimate==1
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
            results = sPmcmc.results;
            [ff,lb,ub,x0,vcov,self]=pull_objective(sv);
            
            
            % Collect Posterior Estimates and Confidence Bands
            pmode=posterior_mode(sv);

            pnames=fieldnames(pmode);

            a2tilde_to_a=sv.estim_.linres.a2tilde_to_a;
            
            params=[results.pop.x];

            params_M = prctile(params,50,2);
            params_UB = prctile(params,95,2);
            params_LB = prctile(params,5,2);

            print_structural_form(sv)

            F = fieldnames(pmode);
            C = num2cell(struct2cell(pmode),3);
            C = cellfun(@(c)[c{:}],C,'uni',0);
            params_pmode = vertcat(C{:});

            params_M_Struct=a2tilde_to_a(params_M);
            params_UB_Struct=a2tilde_to_a(params_UB);
            params_LB_Struct=a2tilde_to_a(params_LB);

            paramsAll = a2tilde_to_a(params);
            for ii=1:numel(F)
                eval(['pmcmc.' F{ii} '= paramsAll(ii,:) ;']) 
            end

        %% Collect Parameters
            paramsCI_mode = [params_pmode params_pmode-(params_M_Struct-params_LB_Struct)  params_pmode+(params_UB_Struct-params_M_Struct)];
            paramsCI_mcmc = [params_M_Struct params_LB_Struct  params_UB_Struct];

            %% Print Table of Estimates
            tableTitle = ['Posterior/Distribution:/' fcstype ',/' sheetuse ',/h/=/' num2str(opt.hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
            pmcmc_tex(pmode,paramsCI_mcmc,texfolder,modelname,tableTitle);

            tableTitleMode = ['Posterior/Mode/Distribution:/' fcstype ',/' sheetuse ',/h/=/' num2str(opt.hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
            pmcmc_tex(pmode,paramsCI_mode,texfolder,[modelname '_mode'],tableTitleMode);


        end
        
        diary off
        
        
        %% Marginal data density
        if opt.mdd_estimate==1
            mddfile = ([logfolder 'MDD_' modelname '.txt']);
            if exist(mddfile, 'file')
                delete(mddfile);
            end
            diary(mddfile)

            mddobj=mdd(results,ff,lb,ub);

            fprintf('Importance sampling::%0.2f\n',is(mddobj,[],mdd.global_options))

            fprintf('Reciprocal Importance sampling::%0.2f\n',ris(mddobj,[],mdd.global_options))

            fprintf('Meng and Wong''s bridge::%0.2f\n',bridge(mddobj,true,mdd.global_options))

            fprintf('Ulrich Mueller::%0.2f\n',mueller(mddobj,[],mdd.global_options))

            fprintf('Laplace::%0.2f\n',laplace(mddobj))

            fprintf('Sims, Waggoner and Zha::%0.2f\n',swz(mddobj,[],mdd.global_options))

            fprintf('Laplace MCMC::%0.2f\n',laplace_mcmc(mddobj))

            fprintf('Chib and Jeliazkov::%0.2f\n',cj(mddobj,[],mdd.global_options))

            diary off
        end
        
    end
    
    %%
    %======================================================================
    % Posterior Mode Analysis
    %====================================================================== 
    
    % Posterior mode values
    pmode=posterior_mode(sv);
    
    % Map posterior mode coefficient estimates
    fnames = fieldnames(pmode);
    for jj=1:length(fnames)
        eval([fnames{jj} ' = pmode.' fnames{jj} ';']);
    end

    
    %%
    %======================================================================
    % Predictive Density at Posterior Mode or MCMC
    %======================================================================
    if strcmp(opt.paramsUse,'mode')
        params_in   = pmode;
        opt.nParamDraws = 1; %length(pmode.a12);
    elseif strcmp(opt.paramsUse,'mcmc')
        params_in   = pmcmc;
        opt.nDraws  = 10;
    end
    
    
    %% Print A0, A1 and SIG matrices

    if opt.printsol
        % Print Structural Form
        print_structural_form(sv)

        % Print Reduced Form
        print_solution(sv)
    end
%     if nlags==1
%         outparams = scriptParams_FULL(params_in,1,modelspec);
%     else
%         outparams = scriptParams_FULL_companion(sv,params_in,1); % works for any lag and model
%     end
    outparams = scriptParams_FULL_companion(sv,params_in,1); % works for any lag and model

    if opt.printsol
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
    end
    
    
    %%
    %======================================================================
    % Construct Densities
    %======================================================================

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
    
    % Predictive density at all episodes
    for ii=1:(opt.tperiods)
        [pdfi_full(ii,:),xi_full(ii,:)]  = ksdensity(pdensity.full(:,ii));
        [cxi_full(ii,:),cdfi_full(ii,:)] = ksdensity(pdensity.full(:,ii),'Function','icdf');
    end
    
   
    % Predictive score at all episodes
    for ii=1:(opt.tperiods)
%         prob = trapz(xi_full(ii,:),pdfi_full(ii,:))
        [~,xi_ps]  = min(abs(xi_full(ii,:)-GDPGH_full_wt(ii,:)));
        if isempty(xi_ps) || xi_ps==1 || xi_ps==100
            ps(ii,:) = 0;
        else
            ps(ii,:) = pdfi_full(ii,xi_ps);
        end
    end
    
    % PITs at all episodes
    for ii=1:(opt.tperiods)
%         prob = trapz(xi_full(ii,:),pdfi_full(ii,:))
        [~,xi_pits]  = min(abs(cxi_full(ii,:)-GDPGH_full_wt(ii,:)));
        if isempty(xi_pits) || xi_pits==1 || xi_pits==100
            pits(ii,:) = 0;
        else
            pits(ii,:) = cdfi_full(ii,xi_pits);
        end
    end
    
%     %%
%     %======================================================================
%     % COVID-19 DENSITIES
%     %======================================================================
% 
%     DATA_MAR13 = readtable('data/MacroRisk_March13_US_only.xlsx','Sheet','DFM_73_Monthly');
%     DATA_APR2 = readtable('data/MacroRisk_April2_US_only.xlsx','Sheet','DFM_73_Monthly');
%     DATA_CURRENT = readtable('data/temp.xlsx','Sheet','US_DFM');
%      
%     [~, db_mar13,~,~,~] = fLoadData('MacroRisk_March13_US_only','DFM_73_Monthly',start_date,end_date,opt.hh);
%     [~, db_apr2,~,~,~] = fLoadData('MacroRisk_April2_US_only','DFM_73_Monthly',start_date,end_date,opt.hh);
% 
%     
%     %% Dates for density
%     covid_dates=[
%         datenum('2020-Mar',dataformat)
%         ];
%     
%      opt.date_index = find(ismember(dates_full, covid_dates));
%     
%      % March-13-2020
%     FF_COVID = db_mar13.FF.data; 
%     MF_COVID = db_mar13.MF.data;
%     TR_COVID = db_mar13.TREND.data;
%     GDPGH_COVID = db_mar13.GDPGH.data;
%     if opt.dir_or_it==1
%         % Direct forecast
%         [quantiles_mar,pdensity_covid_mar] = fPredDensityDitrectFull(sv,params_in,...
%         FF_COVID,MF_COVID,GDPGH_COVID,TR_COVID,opt,modelspec);
%     elseif opt.dir_or_it==2
%         % Iterated forecast                
%         [quantiles_mar, pdensity_covid_mar] = fPredDensityIteratedFull(sv,params_in,...
%             FF_COVID,MF_COVID,GDPGH_COVID,TR_COVID,opt,'forecastst',modelspec);
%     end
%     %%
%     % April-2-2020
%     FF_COVID = db_apr2.FF.data; 
%     MF_COVID = db_apr2.MF.data;
%     TR_COVID = db_apr2.TREND.data;
%     GDPGH_COVID = db_apr2.GDPGH.data;
% 
%     if opt.dir_or_it==1
%         % Direct forecast
%         [quantiles_apr,pdensity_covid_apr] = fPredDensityDirectFull(sv,params_in,...
%         FF_COVID,MF_COVID,GDPGH_COVID,TR_COVID,opt,modelspec);
%     elseif opt.dir_or_it==2
%         % Iterated forecast                
%         [quantiles_apr, pdensity_covid_apr] = fPredDensityIteratedFull(sv,params_in,...
%             FF_COVID,MF_COVID,GDPGH_COVID,TR_COVID,opt,'forecastst',modelspec);
%     end
%     
% %% Compute densities for MS model
% % results correspond to the posterior parameter draws
% 
% 
% [pdfi_mar,xi_mar]  = ksdensity(pdensity_covid_mar.realized);
% [pdfi_mar_good,xi_mar_good]  = ksdensity(pdensity_covid_mar.good);
% [pdfi_mar_bad,xi_mar_bad]  = ksdensity(pdensity_covid_mar.bad);
% preg_mar= 1-quantiles_mar.st_t_mean(567);
% preg_apr= 1-quantiles_apr.st_t_mean(567);
% 
% [pdfi_apr,xi_apr]  = ksdensity(pdensity_covid_apr.realized);
% [pdfi_apr_good,xi_apr_good]  = ksdensity(pdensity_covid_apr.good);
% [pdfi_apr_bad,xi_apr_bad]  = ksdensity(pdensity_covid_apr.bad);
% linetype= {'--','-','-.',':','-.'};
% 
% fig20=figure;
% hold on
% plot(xi_mar,pdfi_mar,linetype{2},'Color',[0 0 0 ],'LineWidth', 3,'DisplayName','Mar-2020');
% plot(xi_mar_bad,pdfi_mar_bad,linetype{1},'Color',[1 0 0 ],'LineWidth', 3,'DisplayName','Mar-2020 (bad)');
% plot(xi_mar_good,pdfi_mar_good,linetype{1},'Color',[0 0 1 ],'LineWidth', 3,'DisplayName','Mar-2020 (good)');
% 
% plot(xi_apr,pdfi_apr,linetype{2},'Color',colors(20,:),'LineWidth', 3,'DisplayName','Apr-2020');
% plot(xi_apr_bad,pdfi_apr_bad,linetype{4},'Color',[1 0 0 ],'LineWidth', 3,'DisplayName','Apr-2020 (bad)');
% plot(xi_apr_good,pdfi_apr_good,linetype{4},'Color',[0 0 1 ],'LineWidth', 3,'DisplayName','Apr-2020 (good)');
% 
% axis tight
% ylimits = ylim; xlimits = xlim;
% %plot(zeros(1,length(xi_mar)),zeros(1,length(pdfi_mar)),'*','DisplayName',['Prob bad regime = ' num2str(round((1-quantiles.st_t_mean(567))*100,1)) '\%']);
% vline(0,'k--');
% %vline(GDPGH_full_wt(567),'r:','Realized Value');
% legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
% ylabel('PDF','fontsize',10,'interpreter','Latex')
% xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% %set(gca,'XTick',unique(sort([0 round(GDPGH_full(opt.date_index(ii)),2) xlimits(1):2:xlimits(2)])))
% set(gca,'TickLabelInterpreter','Latex')
% axis tight
% %xlim([-18 8])
% set(fig20,'PaperOrientation','portrait');
% set(fig20, 'PaperSize', figSize);
% set(fig20, 'PaperUnits', 'inches');
% set(fig20, 'Units','inches');
% set(fig20, 'PaperPositionMode', 'auto');
% set(fig20, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
% tightfig;
% 
% %%
% fig21= figure
% plot(dates_full,1-quantiles_mar.st_t_mean,'LineWidth',2); hold on
% plot(dates_full,1-quantiles_apr.st_t_mean,'--','LineWidth',2)
% axis tight
% ylabel('Predictive Score','fontsize',10,'interpreter','Latex')
% xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
% datetick('x','yyyy','keepticks')
% %set(gca, 'XLim', [dates_full(sd), dates_full(end)])
% ylim([0 1])
% %hBands=recessionplot;
% set(hBands,'FaceColor','r','FaceAlpha',0.3)
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', FontSize);
% set(gca,'Layer','top')
% set(gca,'TickLabelInterpreter','Latex')
% % Set figure appereance
% set(fig21,'PaperOrientation','portrait');
% set(fig21, 'PaperSize', figSize);
% set(fig21, 'PaperUnits', 'inches');
% set(fig21, 'Units','inches');
% set(fig21, 'PaperPositionMode', 'auto');
% set(fig21, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
% % tightfig; 
%     %%
    %======================================================================
    % Collect Probabilities
    %======================================================================
    [~,~,~,f]=filter(sv);
    % Smoothed probabilities
    p_reg1_smoothed = f.smoothed_regime_probabilities.regime1.data;
    p_reg2_smoothed = f.smoothed_regime_probabilities.regime2.data;
    % Filtered probabilities
    p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
    p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;
    
    
    %%
    %======================================================================
    % Quantile Plots
    %======================================================================
    figureTitleTag = [strrep(spectag,'_',' ') ', ' opt.paramsUse ', ' strrep(sheetuse,'_',' ') ', h=' num2str(opt.hh) ', ' startdb ':' enddb];
    
    % Add additional objects to structure
    if opt.showfig==1
        if opt.selfig(1)
            fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            subplot(211)
            hold on
            l1=plot(dates_full(sd:ed), quantiles.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (hard coded)');
            l2=plot(dates_full(sd:ed), quantiles.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (hard coded)');
            legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            l3=plot(dates_full(sd:ed),GDPGH_full_wt(sd:ed),'-.','Color','k','LineWidth', 2,'DisplayName','Realized Value');
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            if opt.dir_or_it==1
                title(['Quantiles Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            else
                title(['Quantiles Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            end
            axis tight
            numticks = 48;
            set(gca,'XTick',datenum(dates(sd:numticks:ed)))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
        %     ylim([-10 20])

            subplot(212)
            hold on
            l3=plot(dates_full(sd:ed), 1-quantiles.st_t_mean(sd:ed),'-','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Simulated s(t)');
            l4=plot(dates, p_reg2_filtered(sd_varp:ed_varp),'--','Color',colors(20,:),'LineWidth', 2.5,'DisplayName','Filtered s(t)');
            legend([l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;    
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
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
    end

    %%
    %======================================================================
    % Predictive score plot
    %======================================================================
    if opt.showfig==1
        if opt.selfig(6)
            fig6=figure;
            figSize6 = [12 6];
            plot(dates_full(sd:ed),ps(sd:ed),'Color',[0.3 0.3 0.3],'LineWidth', 2.5); 
            if opt.dir_or_it==1
                title(['Direct: ' figureTitleTag ],'Interpreter','Latex');
            else
                title(['Iterated: ' figureTitleTag ],'Interpreter','Latex');
            end
            axis tight
            ylabel('Predictive Score','fontsize',10,'interpreter','Latex')
            xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
            ylim([0 1])
            hBands=recessionplot;
            set(hBands,'FaceColor','r','FaceAlpha',0.3)
            set(gca, 'FontName', 'Times New Roman');
            set(gca, 'FontSize', FontSize);
            set(gca,'Layer','top')
            set(gca,'TickLabelInterpreter','Latex')
            % Set figure appereance
            set(fig6,'PaperOrientation','portrait');
            set(fig6, 'PaperSize', figSize6);
            set(fig6, 'PaperUnits', 'inches');
            set(fig6, 'Units','inches');
            set(fig6, 'PaperPositionMode', 'auto');
            set(fig6, 'Position', [figSize6(1)/5 figSize6(2)/5 figSize6(1) figSize6(2)]);
            tightfig; 
        end
    end
    
    
    %%
    %======================================================================
    % Predictive scores table
    %======================================================================
    
    % Recessions
    rec_ini = {'01-Dec-1969','01-Nov-1973','01-Jan-1980','01-Jul-1981','01-Jul-1990','01-Mar-2001','01-Dec-2007'};
    rec_end = {'01-Nov-1970','01-Mar-1975','01-Jul-1980','01-Nov-1982','01-Mar-1991','01-Nov-2001','01-Jun-2009'};
    recdates =[datenum(cellstr(rec_ini(:)),'dd-mmm-yyyy'),datenum(cellstr(rec_end(:)),'dd-mmm-yyyy')];
%     recdates = load('Data_Recessions.mat','Recessions');
%     recdates = recdates.Recessions;
    [tf,~]=ismember(dates_full,recdates(:,1));
    recs=[1:length(dates_full)];
    recs=recs(tf);
    [tf,~]=ismember(dates_full,recdates(:,2));
    rece=[1:length(dates_full)];
    rece=rece(tf);
    for ii=1:length(recs)
%         psr = ps((recs(ii)-5):(rece(ii)+5));
        psr = ps(recs(ii):rece(ii));
    end
    
    % Expansions
    for ii=1:length(recs)
        if ii==1
            pse = ps(1:(recs(ii)-1));
        else
            pse = ps((rece(ii-1)+1):(recs(ii)+1));
        end
    end

    fprintf('Average Predictive Score: %0.2f\n Sample: %s and %s\n ',mean(ps),datestr(dates_full(1)),datestr(dates_full(end)))
    fprintf('Average Predictive Score in Recessions: %0.2f\n',mean(psr))
    fprintf('Average Predictive Score in Expansions: %0.2f\n',mean(pse))

    
    %%
    %======================================================================
    % PITs CDF Plot
    %======================================================================
    if opt.selfig(7)
        fig7 = figure;
        figSize7 = [10 8]/1.5;
        h=cdfplot(pits);
        hold on
        rvec = (0:0.001:1);
        plot(rvec,rvec,'r--','LineWidth',2);
        hL = legend('Empirical CDF','Uniform CDF');
        set(hL,'interpreter','Latex','Orientation','Vertical','Location','SouthEast')
        legend boxoff
        xlim([0 1])
        ylim([0 1])
        ylabel(['CDF of PITs'],'fontsize',10,'interpreter','Latex')
        xlabel('$\tau$','interpreter','Latex')
        set(gca, 'FontSize', FontSize);
        if opt.dir_or_it==1
            title(['Direct: ' figureTitleTag ],'FontSize',8,'Interpreter','Latex');
        else
            title(['Iterated: ' figureTitleTag ],'FontSize',8,'Interpreter','Latex');
        end
        set(gca, 'FontName', 'Times New Roman');
        set(gca,'Layer','top')
        hold off
        set(fig7,'PaperOrientation','portrait');
        set(fig7, 'PaperSize', figSize7);
        set(fig7, 'PaperUnits', 'inches');
        set(fig7, 'Units','inches');
        set(fig7, 'PaperPositionMode', 'auto');
        set(fig7, 'Position', [figSize7(1)/3 figSize7(2)/5 figSize7(1) figSize7(2)]);
        tightfig;
    end

    %%
    %======================================================================
    % PITs over Time Plot
    %======================================================================
    if opt.selfig(8)
        fig8 = figure;
        figSize8 = [12 6];
        plot(dates_full(sd:ed),pits(sd:ed),'Color',[0.3 0.3 0.3],'LineWidth', 2.5); 
        if opt.dir_or_it==1
            title(['Direct: ' figureTitleTag ],'Interpreter','Latex');
        else
            title(['Iterated: ' figureTitleTag ],'Interpreter','Latex');
        end
        axis tight
        ylabel('PITs','fontsize',10,'interpreter','Latex')
        xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
        datetick('x','yyyy','keepticks')
        set(gca, 'XLim', [dates_full(sd), dates_full(ed)])
        ylim([0 1])
        hBands=recessionplot;
        set(hBands,'FaceColor','r','FaceAlpha',0.3)
        set(gca, 'FontName', 'Times New Roman');
        set(gca, 'FontSize', FontSize);
        set(gca,'Layer','top')
        set(gca,'TickLabelInterpreter','Latex')
        % Set figure appereance
        set(fig8,'PaperOrientation','portrait');
        set(fig8, 'PaperSize', figSize8);
        set(fig8, 'PaperUnits', 'inches');
        set(fig8, 'Units','inches');
        set(fig8, 'Position', [figSize8(1)/5 figSize8(2)/5 figSize8(1) figSize8(2)]);
        tightfig; 
    end
    
    
    %%
    %======================================================================
    % Density plots
    %======================================================================
    if opt.showfig==1
        if opt.selfig(2)
            linetype= {'--','-','-.',':','-.'};
            for ii=1:numel(opt.date_index)

                fig2=figure;
                hold on
                plot(xi(ii,:),pdfi(ii,:),linetype{2},'Color',[0 0 0 ],'LineWidth', 3,'DisplayName',[datestr(dates_full(opt.date_index(ii)),dataformat)]); 
                plot(xi_good(ii,:),pdfi_good(ii,:),linetype{4},'Color',colors(60,:),'LineWidth', 3,'DisplayName','Good regime'); 
                plot(xi_bad(ii,:),pdfi_bad(ii,:),linetype{5},'Color',colors(1,:),'LineWidth', 3,'DisplayName','Bad regime'); 
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
        %         text(GDPGH_full(opt.date_index(ii)),0.95*ylimits(2),'realized value','color','red','FontSize',9);
                vline(GDPGH_full_wt(opt.date_index(ii)),'r:','Realized Value');
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
            end
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

    if opt.showfig==1
        if opt.selfig(3)
            fig3=figure;
            %yyaxis left
            hold on
            plot(dates, fitgood(sd_varp:ed_varp),'Color',colors(1,:),'LineWidth', 1.5)
            plot(dates, fitbad(sd_varp:ed_varp),'--','Color',colors(15,:),'LineWidth', 1.5)
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            leg = {'Good Regime','Bad Regime'};
            axis tight
            ax=gca;
            ax.XTick = datenum(dates(1:numticks:end));
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(1), dates(end)])
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
        end
    end
    
    %% Model Fit
    clearvars y_fit realized
    y_fit = fit1.*p_reg1_smoothed + fit2.*p_reg2_smoothed;
    
    realized = db.GDPG.data ;
    realized = realized(nlags+1:end);
    
    if opt.showfig==1
        if opt.selfig(4)
            fig4=figure;
            hold on
            plot(dates,y_fit(sd_varp:ed_varp),'linewidth',2)
            plot(dates,realized(sd_varp:ed_varp),'k','linewidth',2)
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            legend('Ergodic Fitted Value','Realized','interpreter','Latex')
            legend boxoff
            axis tight
            ax=gca;
            ax.XTick = datenum(dates(1:numticks:end));
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(1), dates(end)])
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
        end
    end
    
    %% Regime Probabilities
    if opt.showfig==1
        if opt.selfig(5)
            fig5=figure;
            hold on
            plot(dates,p_reg2_smoothed(sd_varp:ed_varp),'linewidth',2)
            plot(dates,p_reg2_filtered(sd_varp:ed_varp),'k','linewidth',2)
            hold off
            ylabel('Probability','interpreter','Latex','fontsize',10)
            legend('Smoothed','Filtered','interpreter','Latex')
            legend boxoff
            axis tight
            ax=gca;
            ax.XTick = datenum(dates(1:numticks:end));
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates(1), dates(end)])
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
        end
    end
    
    
    %%
    if opt.showfig && opt.saveit
        if opt.selfig(1)
            print('-dpdf',fig1,[figfolder 'Quantiles_' modelname ],'-bestfit');
        elseif opt.selfig(2)
            print('-dpdf',fig2,[figfolder 'Densities' datestr(dates_full(opt.date_index(ii)),dataformat) '_' modelname ],'-bestfit');
        elseif opt.selfig(3)
            print('-dpdf',fig3,[figfolder 'Fitted_' modelname ],'-bestfit');
        elseif opt.selfig(4)
            print('-dpdf',fig4,[figfolder 'Resids_' modelname ],'-bestfit');
        elseif opt.selfig(5)
            print('-dpdf',fig5,[figfolder 'ProbReg_' modelname ],'-bestfit');
        elseif opt.selfig(6)
            print('-dpdf',fig6,[figfolder 'PredScore_' modelname ],'-bestfit');
        elseif opt.selfig(7)
            print('-dpdf',fig7,[figfolder 'PITsCDF' modelname ],'-bestfit');
        elseif opt.selfig(8)
            print('-dpdf',fig8,[figfolder 'PITs' modelname ],'-bestfit');
        elseif opt.selfig(9)
            script_slides_PhillyFed;            
        end
    end   
    
end
%%
rise_exit
