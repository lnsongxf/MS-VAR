% Understanding Growth-at-Risk: A Markov-Switching Approach
% Model with f_t and m_t endogenous
%
% This script produces posterior mode estimates and predictive
% densities at the posterior mode for any dataset.
%
% There are two estimation options
%  1. Direct: estimates the measurement equation at horizon hh.
%  2. Iterated: estimates the contemporanous measurement equation and
%     then iterates the predictive density through horizon hh.
%=====================================================================

%  AUG 16, 2021 CHANGE: GDP is now t to t+11 instead of t+1 to t+12
%  in direct model

%% housekeeping
clear; clc; tic;
close all; 
%%
% Important paths
addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');         % RISE Toolbox
% addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_Tbx_Beta-master');         % RISE Toolbox
addpath(genpath('functions'));      % Main functions
addpath(genpath('scripts'));        % Main scripts
addpath(genpath('auxtools'));       % Tools for plotting and other
addpath(genpath('cbrewer'));        % Color mixing tools

% Load RISE
rise_startup()

%%
% Date formats
inputformat   = 'yyyy-MMM';
dataformat    = 'yyyy-mmm';

%========================================
% USER INPUT
%========================================
opt.dir_or_it      = 2;     % 1 = direct, 2 = iterated
opt.hh             = 12;    % forecast horizon
opt.optimizer      = 1;     % 1 = fminunc, otherwise fmincon
opt.transprob      = 1;     % 1 = endogenous, otherwise exogenous
opt.simul          = 0;     % 1 = RISE forecast, 0 = our simulation function
opt.simtype        = 1;     % 1 = forecast, 0 = drawst
opt.const          = 1;     % 1 = have a constant in transition probability
opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
opt.paramsUse      = 'mode';
opt.force_estimate = 0;     % 1 = force to reestimate model
opt.mdd_estimate   = 0;     % 1 = force to estimate MDD if paramsUse=mcmc
opt.printsol       = 1;     % 1 = print model solution
opt.saveit         = 0;     % 1 = save figures
opt.showfig        = 1;     % 1 = show figures
opt.selfig         = [1 1 zeros(1,8)];  % 1 = vector of figures
opt.covid          = 0;     % 1 = run Covid exercise

% Figure selection options
% 1 = quantiles predictive distribution
% 2 = densities
% 3 = regime specific fitted values
% 4 = ergodic fitted values
% 5 = regime probabilities
% 6 = predictive score over time
% 7 = PITs CDFs
% 8 = PITs over time
% 9 = Figures for slides (currently matching PhillyFed presentation)
% 10 = Covid Exercise

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
opt.nlags=1;
nlags=opt.nlags;
exog={};
constant=true;
panel=[];

% Select which models to run
% Model	  C	 A0	 A1	 S	GDP(-1) restriction
% 101     s	 n	 n	 n	Yes
% 102     s	 n	 n	 s	Yes
% 103     s	 s	 n	 n	Yes
% 104     s	 n	 s	 n	Yes
% 105     s	 s	 n	 s	Yes
% 106     s	 n	 s	 s	Yes
% 107     s	 s	 s	 s	Yes
% 108  s(3)  n   n s(3) Yes
% 201     s	 n	 n	 n	No
% 202     s	 n	 n	 s	No
% 203     s	 s	 n	 n	No
% 204     s	 n	 s	 n	No
% 205     s	 s	 n	 s	No
% 206     s	 n	 s	 s	No
% 207     s	 s	 s	 s	No
% 208  s(3)  n   n s(3) No
% mymodels = [101:108,205];
mymodels = [105];

if opt.transprob==1
    tag ='';
else
    tag ='_EXO';
end

for imodel=mymodels
    
    modelspec = imodel;
    
    % Model spectag
    spectag = fSpecLabel(modelspec);
        
    % Data vintage, estimation sample and country selection
    datafilename = '11302020';
    sheetuse     = 'US_DFM';
    start_date   = '1973-Jan';
    end_date     = '2020-Jan';    % Last MF/FF observation available for estimation
    
    % Define target periods for density cuts
    target_dates = [
        datenum('2007-Aug',dataformat)
        ];
    
    % Define dates for plots using full dates vector
    start_plot    = '1973-Jan';    
    end_plot      = '2020-Oct';

    % Define dates for plots using available estimation sample
    start_plot_var = '1973-Feb';
    end_plot_var   = '2019-Jan';

    % Define dates for evaluation 
    end_eval       = '2019-Jan';
    
    
    %========================================
    % MODEL SPECIFIC VARIABLES
    %========================================
    
    if opt.dir_or_it ==1
        fcstype= 'direct';
        varlist={'FF','MF','GDPGH'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' spectag '/'];
        texfolder    = ['results/Tex/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' spectag  '/'];
        logfolder    = ['results/Logs/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag  '/'];
        figfolder    = ['results/Figures/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/'];
        % figslides   =  ['results/Figures_Slides/' opt.paramsUse '/VAR' num2str(nlags) '/Direct/' spectag '/'];

    elseif opt.dir_or_it ==2
        fcstype= 'iterated';
        varlist={'FF','MF','GDPG'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/'];
        texfolder    = ['results/Tex/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag  '/'];
        logfolder    = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags)  '/' spectag '/'  ];
        figfolder    =  ['results/Figures/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/'];
        % figslides   =  ['results/Figures_Slides/' opt.paramsUse '/VAR' num2str(nlags) '/Iterated/' spectag '/'];

    end
    
    %========================================
    % CONFIGURE DATES AND PLOT OPTIONS
    %========================================
    
    % All dates
    dates_all = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_plot,'InputFormat',inputformat)))';
    
    % Vector of dates for the full sample    
    dates_full = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_plot,'InputFormat',inputformat)))';

    % Starting date for plots
    sd = find(datenum((datetime(start_plot,'InputFormat',inputformat)))==dates_full);   
    sd_varp = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates_full);   
        
    % End date for plots
    ed_full = find(datenum(end_plot,dataformat)==dates_full);
    ed_varp = find(datenum((datetime(end_date,'InputFormat',inputformat)-calmonths(opt.hh)))==dates_full); 
    ed_edate= find(datenum(end_date,dataformat)==dates_full); 
    ed_eval = find(datenum(end_eval,dataformat)==dates_full);        
    
    % Vector of dates for the available estimation sample 
    dates    = dates_full(sd_varp:ed_varp);

    % Number of time-periods in simulation of GDP
    opt.tperiods = length(dates_full(sd:ed_full));  
     
    % Create folders to store results
    if exist(paramsfolder,'dir')==0;  mkdir(paramsfolder); end
    if exist(texfolder,'dir')==0;  mkdir(texfolder); end
    if exist(logfolder,'dir')==0;  mkdir(logfolder); end
    if exist(figfolder,'dir')==0;  mkdir(figfolder); end
    % if exist(figslides,'dir')==0;  mkdir(figslides); end

%%
    %==========================================================================
    % LOAD DATA
    %==========================================================================    
    
    [db, db_full,startdb,enddb,tex] = fLoadData(datafilename,sheetuse,start_date,end_date,opt);
    
    % enddb returns the last period for which one can compute all the
    % elemennts in Y_t given that last available data is end_date. 
 
    
    % Collect MF, FF and trend (estimation sample, adjusted by VAR lags)
    FF    = db.FF.data;
    MF    = db.MF.data;
    TR    = db.TRENDH.data;
    GDPG  = db.GDPG.data;
    GDPGH = db.GDPGH.data;
    
    % Collect MF, FF and trend (full sample, adjusted by VAR lags)
    % These are the objects passed to the simulation routines
    FF_full = db_full.FF.data;
    MF_full = db_full.MF.data;
    
    % This is GDP(t:t+h-1)
    GDPG_full     = db_full.GDPG.data;
    TR_full       = db_full.TRENDH.data;
    GDPGH_full    = db_full.GDPGH.data;
    GDPGH_full_wt = GDPGH_full + TR_full;

    % To plot GDP(t+1:t+h) we need to offset GDPGH_full_wt by one period
    % when plotting
    GDPGH_full_wt    = [GDPGH_full_wt(2:end); NaN];
    
    % Complete the trend with the last observation
    for jj=1:length(TR_full)
        if TR_full(jj)==0; TR_full(jj) = TR_full(jj-1); end
    end
    
    % Model name and title
    % Note: iterated model posterior mode depends on opt.hh through the estimation sample
    d_temp = datestr(dates_full(ed_varp));
    enddbf = [num2str(year(d_temp)) 'M' num2str(month(d_temp))];

    modelname = [datafilename '_' sheetuse '_' startdb enddbf ...
        '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];    
 %%
    % Plot Data
    %{
    if opt.showfig==1
        FontSize  = 16;
        numticks  = 24;
        figSize   = [12 6];
        linestyle = {'-','--',':'};
        colors    = cbrewer('div', 'RdYlBu', 64);
        colors2   = cbrewer('qual', 'Set1', 8);
        left_color= [0 0 0];
        right_color= colors(2,:);
    
        fig0=figure;
        plot(dates_full(sd:ed_full),[GDPG_full GDPGH_full],'LineWidth',1);
        legend('GDPG','GDPGH')
        axis tight
        set(gca,'XTick',datenum(dates_full(sd:48:ed_full)))
        datetick('x','yyyy','keepticks')
        set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
        % Set figure appereance
        set(fig0,'PaperOrientation','portrait');
        set(fig0, 'PaperSize', figSize);
        set(fig0, 'PaperUnits', 'inches');
        set(fig0, 'Units','inches');
        set(fig0, 'PaperPositionMode', 'auto');
        set(fig0, 'Position', [figSize(1)/2 figSize(2)/5 figSize(1) figSize(2)]);
        tightfig;
    end
    %}
    
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
        
        
        fprintf('\n *** TP Params: State = 1 *** \n')
        disp([outparams.a12 outparams.b12 outparams.c12])    
        fprintf('\n *** TP Params: State = 2 *** \n')
        disp([outparams.a21 outparams.b21 outparams.c21])    
    end
    
    
    %%
    %======================================================================
    % Construct Densities
    %======================================================================

    % Extract the index for the desired density plots
    opt.date_index = find(ismember(dates_full, target_dates));
    
    % Determine simulation type
    if opt.simtype==1
        simtype = 'forecastst';
    else
        simtype = 'drawst';
    end

    % Simulate densities
    if opt.dir_or_it==1
        % Direct forecast
        [quantiles,pdensity] = fPredDensityDirectFull(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt,simtype);
        %[quantiles,pdensity] = fPredDensityDirectFull_RISE(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt,simtype,dates_full,db_full);
    elseif opt.dir_or_it==2
        % Iterated forecast      
        [quantiles, pdensity] = fPredDensityIteratedFull(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,simtype);
%         [quantiles, pdensity] = fPredDensityIteratedFull_RISE(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,simtype,dates_full,db_full);
    end
    %%
    % Predictive density and CDF at specified episodes
    for ii=1:numel(opt.date_index)
        %[pdfi(ii,:),xi(ii,:)]  = ksdensity(pdensity(:,ii));
        [pdfi(ii,:),xi(ii,:)]  = ksdensity(pdensity.realized(:,ii));
        [pdfi_good(ii,:),xi_good(ii,:)]  = ksdensity(pdensity.good(:,ii));
        [pdfi_bad(ii,:),xi_bad(ii,:)]  = ksdensity(pdensity.bad(:,ii));
        [cxi(ii,:),cdfi(ii,:)] = ksdensity(pdensity.realized,'Function','icdf');
   end
    
    % Predictive density for all quantile simulations
    if opt.simul==1
        for ii=(1+nlags):(opt.tperiods)
            [pdfi_full(ii,:),xi_full(ii,:)]  = ksdensity(pdensity.full(:,ii));
            [cxi_full(ii,:),cdfi_full(ii,:)] = ksdensity(pdensity.full(:,ii),'Function','icdf');
        end
    else
        for ii=1:(opt.tperiods)
            if isnan(mean(pdensity.full(:,ii)))
                pdfi_full(ii,:) = NaN(1,100);
                xi_full(ii,:) = NaN(1,100);
                cxi_full(ii,:) = NaN(1,100);
                cdfi_full(ii,:) = NaN(1,100);
            else
                [pdfi_full(ii,:),xi_full(ii,:)]  = ksdensity(pdensity.full(:,ii));
                [cxi_full(ii,:),cdfi_full(ii,:)] = ksdensity(pdensity.full(:,ii),'Function','icdf');
            end
        end
    end
    
   %%
    % Predictive score at all episodes
    % Do not use last hh observations since no data for GDPG
    for ii=1:ed_eval
%         prob = trapz(xi_full(ii,:),pdfi_full(ii,:))
        [~,xi_ps]  = min(abs(xi_full(ii,:)-GDPGH_full_wt(ii,:)));
        if isempty(xi_ps) || xi_ps==1 || xi_ps==100
            ps(ii,:) = 0;
        else
            ps(ii,:) = pdfi_full(ii,xi_ps);
        end
    end
    
    % PITs at all episodes
    % Do not use last hh observations since no data for GDPG
    
    for ii=1:ed_eval
%         prob = trapz(xi_full(ii,:),pdfi_full(ii,:))
        [~,xi_pits]  = min(abs(cxi_full(ii,:)-GDPGH_full_wt(ii,:)));
        if isempty(xi_pits) || xi_pits==1 || xi_pits==100
            pits(ii,:) = 0;
        else
            pits(ii,:) = cdfi_full(ii,xi_pits);
        end
    end
    
    results.outparams  = outparams; % parameters
    results.quantiles  = quantiles; % quantiles and transition probability
    results.pdfi       = pdfi;      % PDF ergodic
    results.pdfi_good  = pdfi_good; % PDF good regime
    results.pdfi_bad   = pdfi_bad;  % PDF bad regime
    results.xi         = xi;        % support of PDF ergodic
    results.xi_good    = xi_good;   % support of PDF good regime
    results.xi_bad     = xi_bad;    % support of PDF bad regime
    results.cxi        = cxi;       % CDF ergodic
    results.cdfi       = cdfi;      % support of CDF
    results.ps         = ps;        % predictive score
    results.pits       = pits;      % PIT
    
    save([paramsfolder modelname '_results.mat'],'results')

    %%
    %======================================================================
    % Collect Probabilities
    %======================================================================
    % [LogLik,Incr,retcode,filtering] = filter(sv);
    % Returns:
    %
    %    LogLik : log likelihood values of the VAR
    %
    %    Incr : period-by-period contributions to the likelihood
    %
    %    retcode : flag, 0 if no problem
    %
    %    filtering : structure containing filtration
    %%
    [~,~,~,f]=filter(sv);
    % Probabilities
    % filtered_variables refer to one-step ahead forecasts: a_{t|t-1}
    % updated_variables refer to the updates a_{t|t}
    % smoothed_variables refer to the final estimates conditional on all
    % available information a_t{t|n}    
    
    % Filtered probabilities
    p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
    p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;
    
    % Updated probabilities
    p_reg1_updated = [f.updated_regime_probabilities.regime1.data];
    p_reg2_updated = [f.updated_regime_probabilities.regime2.data];
    
    % Smoothed probabilities
    p_reg1_smoothed = [f.smoothed_regime_probabilities.regime1.data];
    p_reg2_smoothed = [f.smoothed_regime_probabilities.regime2.data];
    
    
    %%
    %======================================================================
    % Quantile Plots
    %======================================================================
    figureTitleTag = [strrep(spectag,'_',' ') ', ' opt.paramsUse ', ' strrep(sheetuse,'_',' ') ', h=' num2str(opt.hh) ', ' startdb ':' enddb];
    
    % Figure options
    FontSize  = 16;
    numticks  = 24;
    figSize   = [12 6];
    linestyle = {'-','--',':'};
    colors    = cbrewer('div', 'RdYlBu', 64);
    colors2   = cbrewer('qual', 'Set1', 8);
    left_color= [0 0 0];
    right_color= colors(2,:);

    % Add additional objects to structure
    if opt.showfig==1
        if opt.selfig(1)
            fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
            subplot(211)
            hold on
            l1=plot(dates_full(sd:ed_full), quantiles.dYsim_10(sd:ed_full),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th');
            l2=plot(dates_full(sd:ed_full), quantiles.dYsim_90(sd:ed_full),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th');
            legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            l3=plot(dates_full(sd:ed_varp),GDPGH_full_wt(sd:ed_varp),'-.','Color','k','LineWidth', 2,'DisplayName','Realized Value');
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            if opt.dir_or_it==1
                title(['Quantiles Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            else
                title(['Quantiles Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            end
            axis tight
            numticks = 24;
            set(gca,'XTick',datenum(dates_full(sd:numticks:ed_full)))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
        %     ylim([-10 20])

            subplot(212)
            hold on
            l3=plot(dates_full(sd:ed_full), 1-quantiles.st_t_mean(sd:ed_full),'-','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Simulated $s(t+1|t)$');
            l4=plot(dates_full(sd:ed_varp), p_reg2_filtered(sd:ed_varp),'--','Color',colors(20,:),'LineWidth', 2.5,'DisplayName','Filtered $s(t+1|t)$');
            l5=plot(dates_full(sd_varp:ed_varp), p_reg2_updated(sd_varp:ed_varp),'-.','Color',colors(30,:),'LineWidth', 2.5,'DisplayName','Updated $s(t|t)$');
            l6=plot(dates_full(sd_varp:ed_varp), p_reg2_smoothed(sd_varp:ed_varp),':','Color',colors(45,:),'LineWidth', 2.5,'DisplayName','Smoothed $s(t|T)$');
            l7=plot(dates_full(sd:ed_full), 1-quantiles.st_h_mean(sd:ed_full),'-','Color','k','LineWidth', 2.5,'DisplayName','Simulated $s(t+1,t+h|t)$');
            legend([l3 l4 l5 l6 l7],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;    
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
            axis tight
            set(gca,'XTick',datenum(dates_full(sd:numticks:ed_full)))
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
            plot(dates_full(sd:ed_eval),ps(sd:ed_eval),'Color',[0.3 0.3 0.3],'LineWidth', 2.5); 
            if opt.dir_or_it==1
                title(['Direct: ' figureTitleTag ],'Interpreter','Latex');
            else
                title(['Iterated: ' figureTitleTag ],'Interpreter','Latex');
            end
            axis tight
            ylabel('Predictive Score','fontsize',10,'interpreter','Latex')
            xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd), dates_full(ed_eval)])
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

    fprintf('Average Predictive Score: %0.2f\n Sample: %s and %s\n ',mean(ps(sd:ed_eval)),datestr(dates_full(1)),datestr(dates_full(ed_eval)))
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
        plot(dates_full(sd:ed_eval),pits(sd:ed_eval),'Color',[0.3 0.3 0.3],'LineWidth', 2.5); 
        if opt.dir_or_it==1
            title(['Direct: ' figureTitleTag ],'Interpreter','Latex');
        else
            title(['Iterated: ' figureTitleTag ],'Interpreter','Latex');
        end
        axis tight
        ylabel('PITs','fontsize',10,'interpreter','Latex')
        xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
        datetick('x','yyyy','keepticks')
        set(gca, 'XLim', [dates_full(sd), dates_full(ed_eval)])
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
                %plot(zeros(1,length(xi_bad(1,:))),zeros(1,length(pdfi_bad(1,:))),'*','DisplayName',['Prob bad regime = ' num2str(round((1-quantiles.st_t_mean(opt.date_index(ii)))*100,1)) '\%']);
                plot(zeros(1,length(xi_bad(1,:))),zeros(1,length(pdfi_bad(1,:))),'*','DisplayName',['Prob bad regime $(t+1,t+h|t)$ = ' num2str(round((1-quantiles.st_h_mean(opt.date_index(ii)))*100,1)) '\%']);
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
    
    numticks = 48;

    if opt.showfig==1
        if opt.selfig(3)
            fig3=figure;
            %yyaxis left
            hold on
            plot(dates_full(sd:ed_varp), fitgood(sd:ed_varp),'Color',colors(1,:),'LineWidth', 1.5)
            plot(dates_full(sd:ed_varp), fitbad(sd:ed_varp),'--','Color',colors(15,:),'LineWidth', 1.5)
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            leg = {'Good Regime','Bad Regime'};
            axis tight
            ax=gca;
            ax.XTick = datenum(dates_full(sd:numticks:ed_varp));
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd), dates_full(ed_varp)])
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
    
    numticks=48;
    
    if opt.showfig==1
        if opt.selfig(4)
            fig4=figure;
            hold on
            plot(dates_full(sd:ed_varp),y_fit(sd:ed_varp),'linewidth',2)
            plot(dates_full(sd:ed_varp),realized(sd:ed_varp),'k','linewidth',2)
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            legend('Ergodic Fitted Value','Realized','interpreter','Latex')
            legend boxoff
            axis tight
            ax=gca;
            ax.XTick = datenum(dates_full(sd:numticks:ed_varp));
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd), dates_full(ed_varp)])
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
            plot(dates_full(sd:ed_varp),p_reg2_smoothed(sd:ed_varp),'linewidth',2)
            plot(dates_full(sd:ed_varp),p_reg2_filtered(sd:ed_varp),'k','linewidth',2)
            hold off
            ylabel('Probability','interpreter','Latex','fontsize',10)
            legend('Smoothed $p(s_t|T)$','Filtered $p(s_{t}|t-1)$','interpreter','Latex')
            legend boxoff
            axis tight
            ylim([0 1.15]);
            ax=gca;
            ax.XTick = datenum(dates_full(sd:numticks:ed_varp));
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd), dates_full(ed_varp)])
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
    
    %% Save Figures
    
    if opt.showfig && opt.saveit
        if opt.selfig(1)
            print('-dpdf',fig1,[figfolder 'Quantiles_' modelname ],'-bestfit');
        end
        if opt.selfig(2)
            print('-dpdf',fig2,[figfolder 'Densities_' datestr(dates_full(opt.date_index(ii)),dataformat) '_' modelname ],'-bestfit');
        end
        if opt.selfig(3)
            print('-dpdf',fig3,[figfolder 'Fitted_' modelname ],'-bestfit');
        end
        if opt.selfig(4)
            print('-dpdf',fig4,[figfolder 'Resids_' modelname ],'-bestfit');
        end
        if opt.selfig(5)
            print('-dpdf',fig5,[figfolder 'ProbReg_' modelname ],'-bestfit');
        end
        if opt.selfig(6)
            print('-dpdf',fig6,[figfolder 'PredScore_' modelname ],'-bestfit');
        end
        if opt.selfig(7)
            print('-dpdf',fig7,[figfolder 'PITsCDF_' modelname ],'-bestfit');
        end
        if opt.selfig(8)
            print('-dpdf',fig8,[figfolder 'PITs_' modelname ],'-bestfit');
        end
        if opt.selfig(9)
            script_slides_PhillyFed;            
        end
    end   
    

    %% Covid Densities
    if opt.showfig && opt.covid
        covid_densities
    end
end

fprintf('Removing path from RISE toolbox');
rise_exit
