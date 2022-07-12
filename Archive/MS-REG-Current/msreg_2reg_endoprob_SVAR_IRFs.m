% Understanding Growth-at-Risk: A Markov-Switching Approach
% Model with f_t and m_t endogenous
%
% This script produces IRFS at posterior mode estimates for different model
% configurations of the Iterated Specification
%
% There are two estimation options
%  2. Iterated: estimates the contemporanous measurement equation and
%     then iterates the predictive density through horizon hh.
%=====================================================================

%  AUG 16, 2021 CHANGE: GDP is now t to t+11 instead of t+1 to t+12
%  in direct model

%% housekeeping
clear; clc; tic;
set(0,'defaulttextinterpreter','latex')
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
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
opt.dir_or_it      = 2;     % 1 = direct, 2 = iterated [ONLY WORKS FOR 2]
opt.hh             = 12;    % forecast horizon
opt.optimizer      = 1;     % 1 = fminunc, otherwise fmincon
opt.transprob      = 1;     % 1 = endogenous, otherwise exogenous
opt.simul          = 0;     % 1 = RISE forecast, 0 = our simulation function
opt.simtype        = 1;     % 1 = forecast, 0 = drawst
opt.const          = 1;     % 1 = have a constant in transition probability
opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
opt.paramsUse      = 'mode';
opt.force_estimate = 0;     % 1 = force to reestimate model
opt.printsol       = 1;     % 1 = print model solution
opt.saveit         = 1;     % 1 = save figures
opt.showfig        = 1;     % 1 = show figures
opt.shock          = 1;     % size of shock (in s.d.)
opt.which_shock    = 1;     % 1 = shock to the financial factor, 2 = shock to macro factor

% IRF options
opt_irf = opt;              % Initialize options
opt_irf.nDrawsU  = 1;       % number of draws of regime shocks per iteration
opt_irf.nDrawsEPS= 10000;   % Number of iterations

% Size of impulses
sd_shock.f = +3;            % +f = tighter financial conditions
sd_shock.m = -3;            % -m = weaker macroeconomic conditions
sd_shock.y = -3;            % -y = negative GDP shock

opt_irf.scenario=0;         % 0: shocks are analyzed separetely 1: shocks are combined


% MCMC options
options = struct();
options.alpha = 0.234;
options.thin = 10;
options.burnin = 10^3;
options.N = 2*10^4;
options.nchain = 1;
opt.nParamDraws  = 1000;    % Draws from posterior. has to be less than options.N;

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

if opt.which_shock ==1
    shock_tag = 'financial';
elseif opt.which_shock ==2
    shock_tag = 'macro';
end
if opt.shock >0
    shock_tag = [shock_tag '_pos_' num2str(abs(opt.shock))];
else
    shock_tag = [shock_tag '_neg_' num2str(abs(opt.shock))];
end

for imodel=mymodels
    
    modelspec = imodel;
    
    % Model spectag
    spectag = fSpecLabel(modelspec);
        
    % Data vintage, estimation sample and country selection
    datafilename = '11302020';
    sheetuse     = 'US_DFM';
    start_date   = '1973-Jan';
    end_date     = '2019-Dec';    % Last MF/FF observation available for estimation
    
    % Define target periods for density cuts
    target_dates = [
        % datenum('1994-Jan',dataformat)
        datenum('2008-Jun',dataformat)
        ];
    
    % Define dates for plots using full dates vector
    start_plot    = '1973-Feb';     
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
        figfolder    = ['results/Figures/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/IRFs/'];
        % figslides   =  ['results/Figures_Slides/' opt.paramsUse '/VAR' num2str(nlags) '/Direct/' spectag '/'];

    elseif opt.dir_or_it ==2
        fcstype= 'iterated';
        varlist={'FF','MF','GDPG'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/'];
        texfolder    = ['results/Tex/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag  '/'];
        logfolder    = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags)  '/' spectag '/'  ];
        figfolder    =  ['results/Figures/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/IRFs/'];
        % figslides   =  ['results/Figures_Slides/' opt.paramsUse '/VAR' num2str(nlags) '/Iterated/' spectag '/'];

    end
    
    %========================================
    % CONFIGURE DATES AND PLOT OPTIONS
    %========================================
    
    % All dates
    dates_all = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_plot,'InputFormat',inputformat)))';
    
    % Vector of dates for the full sample    
    dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_plot,'InputFormat',inputformat)))';

    % Starting date for plots
    sd      = find(datenum((datetime(start_plot,'InputFormat',inputformat)))==dates_full);   
        
    % End date for plots
    ed_full = find(datenum(end_plot,dataformat)==dates_full);
    ed_varp = find(datenum((datetime(end_date,'InputFormat',inputformat)-calmonths(opt.hh-1)))==dates_full); 
    ed_edate= find(datenum(end_date,dataformat)==dates_full); 
    ed_eval = find(datenum(end_eval,dataformat)==dates_full);        
    
    % Vector of dates for the available estimation sample 
    dates    = dates_full(sd:ed_varp);

    % Number of time-periods in simulation of GDP
    opt.tperiods = length(dates_full(sd:ed_full));  
     
    % Create folders to store results
    if exist(paramsfolder,'dir')==0;  mkdir(paramsfolder); end
    if exist(texfolder,'dir')==0;  mkdir(texfolder); end
    if exist(logfolder,'dir')==0;  mkdir(logfolder); end
    if exist(figfolder,'dir')==0;  mkdir(figfolder); end
    % if exist(figslides,'dir')==0;  mkdir(figslides); end


    %==========================================================================
    % LOAD DATA
    %==========================================================================    
    
    [db, db_full,startdb,enddb,tex] = fLoadData(datafilename,sheetuse,start_date,end_date,opt);
    
    % enddb returns the last period for which one can compute all the
    % elemennts in Y_t given that last available data is end_date. 
 
    
    % Collect MF, FF and trend (estimation sample, adjusted by VAR lags)
    FF    = db.FF.data(nlags+1:end);
    MF    = db.MF.data(nlags+1:end);
    TR    = db.TRENDH.data(nlags+1:end);
    GDPG  = db.GDPG.data(nlags+1:end);
    GDPGH = db.GDPGH.data(nlags+1:end);
    
    % Collect MF, FF and trend (full sample, adjusted by VAR lags)
    % These are the objects passed to the simulation routines
    FF_full = db_full.FF.data(nlags+1:end);
    MF_full = db_full.MF.data(nlags+1:end);
    
    % This is GDP(t:t+h-1)
    GDPG_full     = db_full.GDPG.data(nlags+1:end);
    TR_full       = db_full.TRENDH.data(nlags+1:end);
    GDPGH_full    = db_full.GDPGH.data(nlags+1:end);
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

    %==========================================================================
    % RISE ESTIMATION / LOAD EXISTING ESTIMATION
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
    
    
    % Map parameters to structure
    outparams = scriptParams_FULL_companion(sv,params_in,1); % works for any lag and model

    % Print A0, A1 and SIG matrices
    if opt.printsol
        % Print Structural Form
        print_structural_form(sv)

        % Print Reduced Form
        print_solution(sv)

        % Prints the estimated matrices to the screen
        fPrintSolution(outparams)
    end
    


    % Collect filtered probabilities
    [~,~,~,f]=filter(sv);

    % Filtered probabilities
    p_reg2_updated  = f.updated_regime_probabilities.regime2.data;

    
    % IRF at Posterior Mode
    s2data = ones(1,length(p_reg2_updated));
    s2data(p_reg2_updated>0.9)=2;

    % USER OPTIONS


    
    opt_irf.hh = opt.hh;
    hh = opt_irf.hh;
    %% SETUP
    %Baseline shocks shock matrices
    
    % Draw negative shocks
    pd = makedist('Normal');
    t  = truncate(pd,-4,0);   
    for jj = 1:opt_irf.nDrawsEPS
       shocks_baseline(:,:,jj) = random(t,3,hh);    
    end
    
    %
    %shocks_baseline = 0*randn(3,hh,opt_irf.nDrawsEPS);

    %Generate cell array of target dates and convert full date array to cell
    %array
    % target_dates_cell = {'01-Jan-1999','01-Jan-2000','01-Jan-2007','01-Jan-2008','01-Aug-2008','01-Jan-2009'};

    target_dates_cell = { '01-Jan-1999'};

    dates_full_cell = cellstr(datestr(dates_full));

    %Initialize init structure
    init_irf.s0 = [];
    sd_base = 0; % keep at zero


    %----------------------
    % IRF WITH FIXED REGIME
    %----------------------
    
    opt_fixed = opt_irf;
    opt_fixed.nDrawsEPS= 100;

    shocks_baseline_fixed = shocks_baseline(:,:,1:100);
    [~, pctl_fixed] = fAllDatesAlgorithm0_exante(dates_full_cell,target_dates_cell,FF_full,MF_full,GDPG_full,s2data,p_reg2_updated,init_irf,shocks_baseline_fixed,outparams,opt_fixed,sd_base,sd_shock);

    fields = fieldnames(pctl_fixed);

    % Plot IRF
    suptitle_use = ['Model ' num2str(imodel) ', $\delta=$',num2str(sd_shock.f),', $f_{t-1}=$',num2str(round(pctl_fixed.(fields{1}).f0,2)),', $m_{t-1}=$',num2str(round(pctl_fixed.(fields{1}).m0,2)), ', $y_{t-1}=$',num2str(round(pctl_fixed.(fields{1}).y0,2))];
    filename_use = './Figures/irfFixedRegime.pdf';

    suptitle_use = [fields{1}(1:end-2) ': ' '$\delta=$',num2str(sd_shock.f),', $f_{t-1}=$',num2str(round(pctl_fixed.(fields{1}).f0,2)),', $m_{t-1}=$',num2str(round(pctl_fixed.(fields{1}).m0,2)), ', $y_{t-1}=$',num2str(round(pctl_fixed.(fields{1}).y0,2))];
    filename_use = ['./irfFixedRegime' fields{1}(1:end-2)  '.pdf'];
    plot_fan_gdp_alg0exante_2x2(pctl_fixed,{'Good Regime' 'Bad Regime'}, hh,filename_use,suptitle_use);

    %-------------------------------------------------
    % IRF WITH REGIME CHANGE: NO PARAMETER UNCERTAINTY
    %-------------------------------------------------
    
    [irf, pctl] = fAllDatesAlgorithm9A_exante(dates_full_cell,target_dates_cell,FF_full,MF_full,GDPG_full,s2data,p_reg2_updated,init_irf,shocks_baseline,outparams,opt_irf,sd_base,sd_shock);
    fields = fieldnames(pctl);
    
    %% Plot IRF
    for i=1:length(fields)
        suptitle_use = ['Model ' num2str(imodel) ', ' fields{i} ': ' '$\delta=$',num2str(sd_shock.f),', $f_{t-1}=$',num2str(round(pctl.(fields{i}).f0,2)),', $m_{t-1}=$',num2str(round(pctl.(fields{i}).m0,2)), ', $y_{t-1}=$',num2str(round(pctl.(fields{i}).y0,2))];
        filename_use = ['irfRegimeSwitchNoEPS_' fields{i} '.pdf'];       
        plot_fan_gdp_alg9exante(pctl.(fields{i}),hh,filename_use,suptitle_use, pctl_fixed);
    end
%     
    %% SCENARIO
    
%     opt_irf.scenario=1;         % 0: shocks are analyzed separetely 1: shocks are combined
% 
%     % Impulses
%     sd_shock.f = +1;
%     sd_shock.m = -1;
%     sd_shock.y = 0;
%     opt.scenario=1;
% 
%     
%     [irf, pctl] = fAllDatesAlgorithm9A_exante(dates_full_cell,target_dates_cell,FF_full,MF_full,GDPG_full,s2data,p_reg2_updated,init_irf,shocks_baseline,outparams,opt_irf,sd_base,sd_shock);
%     fields = fieldnames(pctl);
%     
%     % Plot IRF
%     for i=1:length(fields)
%         suptitle_use = ['Model ' num2str(imodel) ', ' fields{i} ': ' '$\delta=$',num2str(sd_shock.f),', $f_{t-1}=$',num2str(round(pctl.(fields{i}).f0,2)),', $m_{t-1}=$',num2str(round(pctl.(fields{i}).m0,2)), ', $y_{t-1}=$',num2str(round(pctl.(fields{i}).y0,2))];
%         filename_use = ['./Figures/irfRegimeSwitch_' fields{i} '.pdf'];
%         plot_fan_gdp_alg9exante(pctl.(fields{i}),hh,filename_use,suptitle_use, pctl_fixed);
%     end


%     pctl_fixed.([fields{1} '_1'])
end

fprintf('Removing path from RISE toolbox');
% rise_exit
