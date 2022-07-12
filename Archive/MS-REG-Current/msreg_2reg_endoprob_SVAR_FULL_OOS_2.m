function msreg_2reg_endoprob_SVAR_FULL_OOS_2(i)
% i = 1;
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

    %  AUG 16, 2021 CHANGE: GDP is now t to t+11 instead of t+1 to t+12
    %  in direct model

    %% housekeeping

    % Important paths
    if isdeployed
        i = str2double(i); % Comes from the SLURM loop
    end

    % Date formats
    inputformat   = 'yyyy-MMM';
    dataformat    = 'yyyy-mmm';

    %========================================
    % USER INPUT
    %========================================
    opt.optimizer      = 1;     % 1 = fminunc, 0 or not specified = fmincon
    opt.transprob      = 1;     % 1 = endogenous, otherwise exogenous
    opt.simul          = 0;     % 1 = RISE forecast, 0 = our simulation function
    opt.simtype        = 1;     % 1 = forecast, 0 = drawst
    opt.const          = 1;     % 1 = have a constant in transition probability
    opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
    opt.paramsUse      = 'mode';
    opt.force_estimate = 1;     % 1 = force to reestimate model
    opt.printsol       = 1;     % 1 = print model solution

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
%     mymodels = [101:108];
    mymodels = [202];

    % Select which horizons to run
    myhorizons = [12];
    
    % Select which forecast types to run: 1 = direct, 2 = iterated
    myforecasts = [2];


    % Declare first sample and last sample
    efirst = 'nov 01 1999'; % A period before the first forecast.
    efirstok = datestr(efirst,dataformat);

% for i = 1:50    

    for imodel=mymodels
        tic;
        modelspec = imodel;

        for ihorizon = myhorizons
            opt.hh = ihorizon; % forecast horizon

            for iforecast=myforecasts
                opt.dir_or_it = iforecast; % 1 = direct, 2 = iterated
                
                if opt.dir_or_it == 1 && opt.hh ==1
                    disp('Skip direct model 1-step ahead');
                else

                    % Model spectag
                    if opt.transprob==1
                        tag ='';
                    else
                        tag ='_EXO';
                    end
                    spectag = fSpecLabel(modelspec);

                    % Data vintage, sample and country selection
                    datafilename = '11302020';
                    sheetuse     = 'US_DFM';
                    start_date   = '1973-Jan';
                    end_date     = datestr(datenum(datetime(efirstok,'InputFormat',inputformat)+calmonths(i)),dataformat);

                    % Define target periods for density cuts
                    target_dates = datenum(datetime(end_date,'InputFormat',inputformat));

                    %========================================
                    % MODEL SPECIFIC VARIABLES
                    %========================================

                    foldertag = 'OOS';

                    if opt.dir_or_it ==1
                        fcstype= 'direct';
                        varlist={'FF','MF','GDPGH'}';
                        % Folders
                        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' spectag '/' foldertag '/'];
                        IS_paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' spectag '/' ];

                    elseif opt.dir_or_it ==2
                        fcstype= 'iterated';
                        varlist={'FF','MF','GDPG'}';
                        % Folders
                        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' foldertag '/'];
                        IS_paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/'];

                    end

                    %========================================
                    % CONFIGURE DATES AND PLOT OPTIONS
                    %========================================

                    % Vector of dates for the full sample and for the available estimation sample
                    dates_full = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

                    % Vector of dates for the available estimation sample
                    sd_var = find(datenum((datetime(start_date,'InputFormat',inputformat)+calmonths(nlags)))==dates_full);
                    ed_var  = find(datenum((datetime(end_date,'InputFormat',inputformat)-calmonths(opt.hh)))==dates_full);
                    dates_var  = dates_full(sd_var:ed_var);

                    % Don't need this as simulationsn will be done only for
                    % the relevant period. Not the whole history. 
                    
                    % Number of time-periods in simulation of GDP
                    opt.tperiods = length(dates_full);

                    % Extract the index for the desired density plots
                    opt.date_index = find(ismember(dates_full, target_dates));

                    % Create folders to store results
                    if exist(paramsfolder,'dir')==0;  mkdir(paramsfolder); end

                    %%
                    %==========================================================================
                    % LOAD DATA
                    % startdb and enddb are passed to scriptEstimateMode_FULL
                    % enddbf is used for filenames to be consistent across
                    % iterated and direct models
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
                    GDPGH_full_wt = [GDPGH_full_wt(2:end); NaN];
                    
                    % Complete the trend with the last observation
                    for jj=1:length(TR_full)
                        if TR_full(jj)==0; TR_full(jj) = TR_full(jj-1); end
                    end
                                        
                    % Model name and title
                    % Note: iterated model posterior mode depends on opt.hh through the estimation sample
                    d_temp = datestr(dates_full(opt.tperiods-opt.hh+1));
                    enddbf = [num2str(year(d_temp)) 'M' num2str(month(d_temp))];
                    
                    modelname = [datafilename '_' sheetuse '_' startdb enddbf ...
                        '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];

                    IS_modelname = [datafilename '_' sheetuse '_' '1973M1' '2019M12' ...
                        '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];

                    out_print = ['Performing iteration ending on ' enddbf ' of model ' num2str(modelspec) ];
                    disp(out_print);

                    %% RISE ESTIMATION

                    %%%%%%%%
                    % MODE %
                    %%%%%%%%
                    if strcmp(opt.paramsUse,'mode')

                        % Load results or estimate posterior mode
                        if exist([paramsfolder modelname '_results.mat'],'file')==0 || opt.force_estimate==1

                            % Estimate posterior mode
                            run scriptEstimateMode_FULL.m;

%                             % Save parameters
%                             save([paramsfolder modelname '.mat'],'sPMode')

                        else
                            % Load stored results
                            load([paramsfolder modelname '_results.mat']);
                            sPMode = results.sPMode;

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

                        else   

                            % Load stored results
                            load([paramsfolder modelname '.mat']);
                            sPmcmc = results.sPmcmc;

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

                        end
                    end

                    %% Posterior Mode Analysis

                    % Posterior mode values
                    pmode=posterior_mode(sv);

                    % Map posterior mode coefficient estimates
                    fnames = fieldnames(pmode);
                    for jj=1:length(fnames)
                        eval([fnames{jj} ' = pmode.' fnames{jj} ';']);
                    end


                    %% Collect parameters

                    if strcmp(opt.paramsUse,'mode')
                        params_in   = pmode;
                        opt.nParamDraws = 1; %length(pmode.a12);
                    elseif strcmp(opt.paramsUse,'mcmc')
                        params_in   = pmcmc;
                        opt.nDraws  = 10;
                    end

                    outparams = scriptParams_FULL_companion(sv,params_in,1); % works for any lag and model

                    %% Construct Densities
                    
                    % Determine simulation type
                    if opt.simtype==1
                        simtype = 'forecastst';
                    else
                        simtype = 'drawst';
                    end
                    
                    
                    % Simulate model                    
                    if opt.dir_or_it==1
                        % Direct forecast
                        [quantiles,pdensity] = fPredDensityDirectFull(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt,simtype);
                    elseif opt.dir_or_it==2
                        % Iterated forecast
                        [quantiles, pdensity] = fPredDensityIteratedFull(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,simtype);
                    end

                    % Quantiles
                    quantiles.dYsim_25  = quantiles.dYsim_25(opt.date_index);
                    quantiles.dYsim_75  = quantiles.dYsim_75(opt.date_index);
                    quantiles.dYsim_10  = quantiles.dYsim_10(opt.date_index);
                    quantiles.dYsim_90  = quantiles.dYsim_90(opt.date_index);
                    quantiles.mean      = quantiles.mean(opt.date_index);
                    quantiles.st_t_mean = quantiles.st_t_mean(opt.date_index);
                    
%                     dYsim_25(i)  = quantiles.dYsim_25;
%                     dYsim_75(i)  = quantiles.dYsim_75;
%                     dYsim_10(i)  = quantiles.dYsim_10;
%                     dYsim_90(i)  = quantiles.dYsim_90;
%                     st_t_mean(i) = quantiles.st_t_mean;
                    
                    % Predictive density at specified episodes
                    [pdfi,xi]  = ksdensity(pdensity.realized);
                    [pdfi_good,xi_good]  = ksdensity(pdensity.good);
                    [pdfi_bad,xi_bad]  = ksdensity(pdensity.bad);
                    [cxi,cdfi] = ksdensity(pdensity.realized,'Function','icdf');

                    %% Forecast only the period we care about                    
%                     cdatesplit = strsplit(enddb,'M');                    
%                     cdate = datetime(str2num(cdatesplit{1}),str2num(cdatesplit{2}),1);
%                     usedate = find(ismember(dates_full, datenum(cdate)));
%                     opt.tperiods = opt.hh;
%                     opt.date_index = opt.hh;                   
%                     
%                     if opt.dir_or_it==1
%                         
%                         [~,~,~,f]=filter(sv);
%                         p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;
%                         
%                         % Direct forecast
%                         FF_in = FF_full(usedate:usedate+opt.hh-1);
%                         MF_in = MF_full(usedate:usedate+opt.hh-1);
%                         YY_in = GDPGH_full(usedate:usedate+opt.hh-1);
%                         TR_in = TR_full(usedate:usedate+opt.hh-1);
%                         PB_in = p_reg2_filtered(usedate);
%                         
%                         % Direct forecast
%                         [quantiles2,pdensity2] = fPredDensityDirectFullv2(sv,params_in,FF_in,MF_in,YY_in,TR_in,opt,simtype,PB_in);
%                         
%                         
%                     elseif opt.dir_or_it==2
%                         % Iterated forecast
%                         [quantiles, pdensity] = fPredDensityIteratedFull(sv,params_in,FF_full(usedate+1:usedate+opt.hh),MF_full(usedate+1:usedate+opt.hh),GDPG_full(usedate+1:usedate+opt.hh),TR_full(usedate+1:usedate+opt.hh),opt,simtype);
%                     end
% 
%                     dYsim2_25(i)  = quantiles2.dYsim_25;
%                     dYsim2_75(i)  = quantiles2.dYsim_75;
%                     dYsim2_10(i)  = quantiles2.dYsim_10;
%                     dYsim2_90(i)  = quantiles2.dYsim_90;
%                     st2_t_mean(i) = quantiles2.st_t_mean;

%%                    
%                     %% Stored Output
                    results.outparams  = outparams; % parameters
                    results.sv         = sv;        % sv
                    results.params_in  = params_in; % parameters (in)
                    results.quantiles  = quantiles; % quantiles and transition probability
                    results.pdfi       = pdfi;      % PDF ergodic
                    results.pdfi_good  = pdfi_good; % PDF good regime
                    results.pdfi_bad   = pdfi_bad;  % PDF bad regime
                    results.xi         = xi;        % support of PDF ergodic
                    results.xi_good    = xi_good;   % support of PDF good regime
                    results.xi_bad     = xi_bad;    % support of PDF bad regime
                    results.cxi        = cxi;       % CDF ergodic
                    results.cdfi       = cdfi;      % support of CDF
                    results.sPMode     = sPMode;    % mode

                    if strcmp(opt.paramsUse,'mcmc')
                        results.sPmcmc = sPmcmc;    % mode
                    end

                    save([paramsfolder modelname '_results.mat'],'results');
                    toc;

                    disp(['Ran Forecast Type ' fcstype ' for Model ' num2str(imodel) ' and Horizon ' num2str(ihorizon)]);
                end
            end
        end
    end
%  end
%%
% figure(1); clf;
% subplot(221)
% plot(dYsim_25); hold on;
% plot(dYsim2_25,'--')
% subplot(222)
% plot(dYsim_75); hold on;
% plot(dYsim2_75,'--')
% 
% subplot(223)
% plot(dYsim_10); hold on;
% plot(dYsim2_10,'--')
% 
% subplot(224)
% plot(dYsim_90); hold on;
% plot(dYsim2_90,'--')