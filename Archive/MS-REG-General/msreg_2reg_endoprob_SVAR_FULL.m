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


%% Housekeeping
clear; clc; tic;
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

%%
%========================================
% BASIC CONFIGURATION
%========================================
% Date formats
inputformat   = 'yyyy-MMM';
dataformat    = 'yyyy-mmm';

% Data vintage, estimation sample and country selection
opt.datafilename = '11302020c';
opt.sheetuse     = 'US_DFM';
opt.start_date   = '1973-Jan';
% To keep estimation sample for direct and iterated models keep this date 
% consistent with the end of observation in dataset - hh periods. 
opt.end_date     = '2019-Dec';          

% Name of variables
% This determines the ordering in the VAR and the number of variables
% opt.vars    = {{'FF','MF','GDPG'};...
%     {'FF','MF','GDPG','FEDFUNDS','PCEPILFE_CCA'}};
% opt.vars    = {{'FF','MF','GDPG'}};
opt.vars    = {{'FF','MF','GDPG','FEDFUNDS','PCEPILFE_CCA'}};

% Adjust labels based on varlist
opt.vlabels     = {'Financial Factor'...
                   'Macroeconomic Factor'...
                   'Real GDP Growth'...
                   'Federal Funds Rate'...
                   'Core PCE Inflation (CCA)'};

% Adjust labels based on varlist
opt.units     = {'Units of Measure'...
                   'Units of Measure'...
                   'Percent'...
                   'Percent'...
                   'Percent'};

               
%========================================
% MODEL SELECTION
%========================================

% Select which specifications to run
mymodels = [105];
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

% Select which forecasting type
myfor = [1];
% 1 = direct, 2 = iterated
% Note: GDP is average from t to t+11 in direct model

%========================================
% MORE OPTIONS
%========================================
opt.hh             = 12;     % forecast horizon
opt.nlags          = 1;      % number of lags
% important: if opt.nlags>1 the routine that pulls the parameters
% (scriptParams_FULL_companion.m) does not work for MCMC
opt.optimizer      = 1;      % 1 = fminunc, otherwise fmincon
opt.transprob      = 1;      % 1 = endogenous, otherwise exogenous
opt.simul          = 0;      % 1 = RISE forecast, 0 = our simulation function
opt.simtype        = 1;      % 1 = forecast, 0 = drawst
opt.const          = 1;      % 1 = have a constant in transition probability
opt.normal         = 0;      % 1 = use normal distribution, 0 = gamma distribution
opt.paramsUse      = 'mode'; % mode = posterior mode, mcmc = parameter uncertainty
opt.force_estimate = 0;      % 1 = force to reestimate model
opt.force_simulate = 0;      % 1 = force to simulate model
opt.mdd_estimate   = 0;      % 1 = force to estimate MDD if paramsUse=mcmc
opt.printsol       = 1;      % 1 = print model solution
opt.saveit         = 0;      % 1 = save figures
opt.showfig        = 1;      % 1 = show figures

%========================================
% IDENTIFICATION OPTIONS
%========================================
opt.identify = 1; % 1 = Choleski, 2 = sign restrictions, 3 = systematic MP

%========================================
% FIGURE OPTIONS
%========================================

% Define target periods for density cuts
target_dates = [
    datenum('2007-Aug',dataformat)
    ];

% Define dates for plots using full dates vector
% This will set the limit of the simulation
opt.start_plot    = '1973-Jan';
opt.end_plot      = '2020-Nov';

% Define dates for evaluation
opt.end_eval       = '2019-Jan';

% Figure selection options
opt.selfig         = [ones(1,3) zeros(1,4) 1 0];  % 1 = select figures
% 1 = quantiles predictive distribution
% 2 = transition probabilities
% 3 = densities
% 4 = regime specific fitted values
% 5 = ergodic fitted values
% 6 = regime probabilities
% 7 = predictive score over time (for GDP)
% 8 = PITs CDFs (for GDP)
% 9 = PITs over time (for GDP)


%========================================
% SIMULATION OPTIONS
%========================================

% Simulation options
opt.nDraws       = 5000;    % Number of simulations for each time period
opt.nParamDraws  = 1000;    % Draws from posterior. has to be less than options.N;

%========================================
% SAMPLER OPTIONS
%========================================

% MCMC options
options = struct();
options.alpha = 0.234;
options.thin = 10;
options.burnin = 10^3;
options.N = 2*10^4;
options.nchain = 1;

%% Get ready

% Tag for transition probability type
if opt.transprob==1
    tag ='';
else
    tag ='_EXO';
end


%% Run Main Code

for ifor = myfor
    
    opt.dir_or_it = ifor;
    % 1 = direct, 2 = iterated
    
    for ivar=1:size(opt.vars,1)
        
        % Extract relevant variable list for this iteration
        opt.varlist = opt.vars{ivar};
                        
        % Identify GDP position in data
        opt.idxg = find(ismember(opt.varlist, 'GDPG'));

        for imodel=mymodels
            
            modelspec = imodel;
            
            % Model spectag
            spectag = fSpecLabel(modelspec);            
            
            %========================================
            % MODEL SPECIFIC VARIABLES
            %========================================
            
            if opt.dir_or_it ==1
                fcstype= 'direct';
                % Folders
                paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(opt.nlags)  '/' spectag '/'];
                texfolder    = ['results/Tex/' opt.paramsUse '/Full/Direct/VAR' num2str(opt.nlags)  '/' spectag  '/'];
                logfolder    = ['results/Logs/' opt.paramsUse '/Full/Direct/VAR' num2str(opt.nlags) '/' spectag  '/'];
                figfolder    = ['results/Figures/' opt.paramsUse '/Full/Direct/VAR' num2str(opt.nlags) '/' spectag '/'];
                % figslides   =  ['results/Figures_Slides/' opt.paramsUse '/VAR' num2str(opt.nlags) '/Direct/' spectag '/'];
                
            elseif opt.dir_or_it ==2
                fcstype= 'iterated';
                % Folders
                paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(opt.nlags) '/' spectag '/'];
                texfolder    = ['results/Tex/' opt.paramsUse '/Full/Iterated/VAR' num2str(opt.nlags) '/' spectag  '/'];
                logfolder    = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(opt.nlags)  '/' spectag '/'  ];
                figfolder    =  ['results/Figures/' opt.paramsUse '/Full/Iterated/VAR' num2str(opt.nlags) '/' spectag '/'];
                % figslides   =  ['results/Figures_Slides/' opt.paramsUse '/VAR' num2str(opt.nlags) '/Iterated/' spectag '/'];
                
            end
            
            %========================================
            % CONFIGURE DATES AND PLOT OPTIONS
            %========================================
            
            % All dates
            dates_all = datenum((datetime(opt.start_date,'InputFormat',inputformat)):calmonths(1):(datetime(opt.end_plot,'InputFormat',inputformat)))';
            
            % Vector of dates for the full sample
            dates_full = datenum((datetime(opt.start_date,'InputFormat',inputformat)):calmonths(1):(datetime(opt.end_plot,'InputFormat',inputformat)))';
            
            % Starting date for plots
            sd = find(datenum((datetime(opt.start_plot,'InputFormat',inputformat)))==dates_full);
            sd_varp = find(datenum((datetime(opt.start_plot,'InputFormat',inputformat))+calmonths(opt.nlags))==dates_full);
            
            % End date for plots
            ed_full = find(datenum(opt.end_plot,dataformat)==dates_full);
            ed_varp = find(datenum((datetime(opt.end_date,'InputFormat',inputformat)-calmonths(opt.hh)))==dates_full);
            ed_edate= find(datenum(opt.end_date,dataformat)==dates_full);
            ed_eval = find(datenum(opt.end_eval,dataformat)==dates_full);
            
            % Vector of dates for the available estimation sample
            dates    = dates_full(sd_varp:ed_varp);
            
            % Number of time-periods in simulation of GDP
            opt.tperiods = length(dates_full(sd:ed_full));
            opt.nvars = length(opt.varlist);
            
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

            [db, db_full,startdb,enddb,tex] = fLoadData(opt);
            
            % enddb returns the last period for which one can compute all the
            % elemennts in Y_t given that last available data is opt.end_date.
            
            % Model name and title
            % Note: iterated model posterior mode depends on opt.hh through the estimation sample
            d_temp = datestr(dates_full(ed_varp));
            enddbf = [num2str(year(d_temp)) 'M' num2str(month(d_temp))];
            
            modelname = [num2str(opt.nvars) 'VAR_' opt.datafilename '_' opt.sheetuse '_' startdb enddbf ...
                '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];
                        
            % Collect realized values of variables (i.e., with trend)
            % Both contemporaneous and average over next hh periods
            for i=1:opt.nvars
                if opt.dir_or_it ==1     
                    if i==opt.idxg
                        % GDP is average over hh periods
                        whichreal = 'H_wt';
                    else
                        % Other variables are forecasted only one period ahead
                        whichreal = '_wt';
                    end
                else
                    % Iterated model always forecasts average over next hh periods
                    whichreal = 'H_wt';
                end
                eval(['realized(:,i)= db_full.' opt.varlist{i} whichreal '.data;']);
            end
            
            % GDP name
            if opt.dir_or_it ==1
                % direct
                opt.gname = 'GDPGH';
            elseif opt.dir_or_it ==2
                % iterated
                opt.gname = 'GDPG';
            end
            opt.varlist(strcmp('GDPG',opt.varlist)) = cellstr(opt.gname);
            
            % For direct model, need to lead GDPGH variable (t+1,t+hh)
            % since GDPGH is average from t to t+hh-1 in direct model
            if opt.dir_or_it ==1     
                realized(:,opt.idxg) = [realized(2:end,opt.idxg); NaN];
            end
            
            
            %%
            %==========================================================================
            % RISE ESTIMATION
            %==========================================================================
            
            % Diary file
            dfile = ([logfolder modelname '_params.txt']);
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
                    tableTitle = ['Posterior/Distribution:/' fcstype ',/' opt.sheetuse ',/h/=/' num2str(opt.hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
                    pmcmc_tex(pmode,paramsCI_mcmc,texfolder,modelname,tableTitle);
                    
                    tableTitleMode = ['Posterior/Mode/Distribution:/' fcstype ',/' opt.sheetuse ',/h/=/' num2str(opt.hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
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
            % Simulation
            %======================================================================
                
            % Extract the index for the desired density plots
            opt.date_index = find(ismember(dates_full, target_dates));
            
            % Determine simulation type
            if opt.simtype==1
                simtype = 'forecastst';
            else
                simtype = 'drawst';
            end
            
            if exist([paramsfolder modelname '_results.mat'],'file')==0 || opt.force_simulate==1
                
                % Simulate densities
                if opt.dir_or_it==1
                    % Direct forecast
                    [quantiles,pdensity] = fPredDensityDirectFull(sv,params_in,db_full,opt,simtype);
                elseif opt.dir_or_it==2
                    % Iterated forecast
                    [quantiles, pdensity] = fPredDensityIteratedFull(sv,params_in,db_full,opt,simtype);
                end
                                             
                % Predictive density and CDF at specified episodes
                for ii=1:numel(opt.date_index)
                    for nn=1:opt.nvars
                        [pdfi(ii,nn,:),xi(ii,nn,:)]  = ksdensity(pdensity.realized(:,nn,ii));
                        [pdfi_good(ii,nn,:),xi_good(ii,nn,:)]  = ksdensity(pdensity.good(:,nn,ii));
                        [pdfi_bad(ii,nn,:),xi_bad(ii,nn,:)]  = ksdensity(pdensity.bad(:,nn,ii));
                        [cxi(ii,nn,:),cdfi(ii,nn,:)] = ksdensity(pdensity.realized(:,nn,ii),'Function','icdf');
                    end
                end
                
                % Predictive density for all quantile simulations
                if opt.simul==1
                    for nn=1:opt.nvars
                        for ii=(1+opt.nlags):(opt.tperiods)
                            [pdfi_full(ii,nn,:),xi_full(ii,nn,:)]  = ksdensity(pdensity.full(:,nn,ii));
                            [cxi_full(ii,nn,:),cdfi_full(ii,nn,:)] = ksdensity(pdensity.full(:,nn,ii),'Function','icdf');
                        end
                    end
                else
                    for ii=1:(opt.tperiods)
                        for nn=1:opt.nvars
                            if isnan(mean(pdensity.full(:,ii)))
                                pdfi_full(ii,nn,:) = NaN(1,opt.nvars,100);
                                xi_full(ii,nn,:) = NaN(1,opt.nvars,100);
                                cxi_full(ii,nn,:) = NaN(1,opt.nvars,100);
                                cdfi_full(ii,nn,:) = NaN(1,opt.nvars,100);
                            else
                                [pdfi_full(ii,nn,:),xi_full(ii,nn,:)]  = ksdensity(pdensity.full(:,opt.nvars,ii));
                                [cxi_full(ii,nn,:),cdfi_full(ii,nn,:)] = ksdensity(pdensity.full(:,opt.nvars,ii),'Function','icdf');
                            end
                        end
                    end
                end
                
                % Now focus on GDP
                id = opt.idxg;

                % Predictive score for GDPGH at all episodes
                % Do not use last hh observations since no data for GDPGH
                for ii=1:ed_eval
                    [~,xi_ps]  = min(abs(xi_full(ii,id,:)-realized(ii,id)));
                    if isempty(xi_ps) || xi_ps==1 || xi_ps==100
                        ps(ii,:) = 0;
                    else
                        ps(ii,:) = pdfi_full(ii,id,xi_ps);
                    end
                end
                
                % PITs for GDPGH at all episodes
                % Do not use last hh observations since no data for GDPGH               
                for ii=1:ed_eval
                    [~,xi_pits]  = min(abs(cxi_full(ii,id,:)-realized(ii,id)));
                    if isempty(xi_pits) || xi_pits==1 || xi_pits==100
                        pits(ii,:) = 0;
                    else
                        pits(ii,:) = cdfi_full(ii,id,xi_pits);
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
            else
                load([paramsfolder modelname '_results.mat'],'results')
            end
            
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
            
            if opt.dir_or_it==1
                forctype = 'Direct';
            else
                forctype = 'Iterated';

            end
            figureTitleTag = [forctype ': ' strrep(spectag,'_',' ') ', ' opt.paramsUse ', ' strrep(opt.sheetuse,'_',' ') ', h=' num2str(opt.hh) ', ' startdb ':' enddb];
            clear temp
            
            % Figure options
            FontSize  = 16;
            TickSize  = 10;
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
                    for nn=1:opt.nvars
                        fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
                        subplot(211)
                        hold on
                        l1=plot(dates_full(sd:ed_full), results.quantiles.dYsim_10(sd:ed_full,nn),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th');
                        l2=plot(dates_full(sd:ed_full), results.quantiles.dYsim_90(sd:ed_full,nn),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th');
                        l3=plot(dates_full(sd:ed_full),realized(sd:ed_full,nn),'-.','Color','k','LineWidth', 2,'DisplayName','Realized Value');
                        set(gca,'children',flipud(get(gca,'children')))
                        hold off
                        ylabel(opt.units{nn},'fontsize',10,'Interpreter','Latex')
                        if opt.dir_or_it==1
                            if nn==opt.idxg
                                mytitle = ['\bf{Average ' opt.vlabels{nn} ' over Next ' num2str(opt.hh) ' Months}'];
                            else
                                mytitle = ['\bf{' opt.vlabels{nn} ' Next Month}'];
                            end
                        else
                            mytitle = ['\bf{Average ' opt.vlabels{nn} ' over Next ' num2str(opt.hh) ' Months}'];
                        end
                        title({['Quantiles ' figureTitleTag ] ; mytitle},'Interpreter','Latex','FontSize',16);
                        axis tight
                        set(gca,'XTick',datenum(dates_full(sd:numticks:ed_full)))
                        datetick('x','yyyy','keepticks')
                        set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
                        hBands=recessionplot;
                        set(hBands,'FaceColor','k','FaceAlpha',0.1)
                        legend([l1 l2 l3],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;

                        subplot(212)
                        hold on
                        l4=plot(dates_full(sd:ed_varp), p_reg2_filtered(sd:ed_varp),'--','Color',colors(20,:),'LineWidth', 2.5,'DisplayName','Filtered $s(t+1|t)$');
                        l5=plot(dates_full(sd_varp:ed_varp), p_reg2_updated(sd_varp:ed_varp),'-.','Color',colors(30,:),'LineWidth', 2.5,'DisplayName','Updated $s(t|t)$');
                        l6=plot(dates_full(sd_varp:ed_varp), p_reg2_smoothed(sd_varp:ed_varp),':','Color',colors(45,:),'LineWidth', 2.5,'DisplayName','Smoothed $s(t|T)$');
                        if opt.dir_or_it==2
                            l7=plot(dates_full(sd:ed_full), 1-results.quantiles.st_h_mean(sd:ed_full),'k-','LineWidth', 2.5,'DisplayName','Simulated $s(t+1,t+h|t)$');
                        else
                            l7=plot(dates_full(sd:ed_full), 1-results.quantiles.st_t_mean(sd:ed_full),'k-','LineWidth', 2.5,'DisplayName','Simulated $s(t+1|t)$');
                        end
                        hold off
                        title('\bf{Probability of Bad Regime}','FontSize',16','Interpreter','Latex');
                        axis tight
                        set(gca,'XTick',datenum(dates_full(sd:numticks:ed_full)))
                        datetick('x','yyyy','keepticks')
                        set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
                        ylim([0 1])
                        hBands=recessionplot;
                        set(hBands,'FaceColor','k','FaceAlpha',0.1)
                        legend([l4 l5 l6 l7],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;                   
                        % Set figure appereance
                        set(fig1,'PaperOrientation','portrait');
                        set(fig1, 'PaperSize', figSize);
                        set(fig1, 'PaperUnits', 'inches');
                        set(fig1, 'Units','inches');
                        set(fig1, 'PaperPositionMode', 'auto');
                        set(fig1, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
                        tightfig;
                        if opt.saveit
                            print('-dpdf',fig1,[figfolder 'Quantiles_' opt.varlist '_' modelname ],'-bestfit');
                        end
                    end
                end
            end
            
            
            %%
            %======================================================================
            % Transition probabilities
            %======================================================================
            idxf = find(ismember(opt.varlist, 'FF'));
            idxm = find(ismember(opt.varlist, 'MF'));

            if opt.const==1
                if opt.transprob==0
                    tp_nb = a12;
                    tp_bn = a21;
                elseif opt.transprob==1
                    if opt.normal==1
                        tp_nb=1./(1+exp(a12-b12*realized(:,idxf)-c12*realized(:,idxm)));
                        tp_bn=1./(1+exp(a21-b21*realized(:,idxf)-c21*realized(:,idxm)));
                    elseif opt.normal==0
                        % This is benchmark specification
                        tp_nb=1./(1+exp(a12-b12*realized(:,idxf)+c12*realized(:,idxm)));
                        tp_bn=1./(1+exp(a21+b21*realized(:,idxf)-c21*realized(:,idxm)));
                    end
                end
            elseif opt.const==0
                if opt.transprob==0
                    errordlg('opt.transprob=0 incompatible with opt.const=0');
                elseif opt.transprob==1
                    if opt.normal==1
                            tp_nb=1./(1+exp(-b12*realized(:,idxf)-c12*realized(:,idxm)));
                            tp_bn=1./(1+exp(-b21*realized(:,idxf)-c21*realized(:,idxm)));
                    elseif opt.normal==0
                            tp_nb=1./(1+exp(-b12*realized(:,idxf)+c12*realized(:,idxm)));
                            tp_bn=1./(1+exp(+b21*realized(:,idxf)-c21*realized(:,idxm)));
                    end
                end
            end            
            
            % Add additional objects to structure
            if opt.showfig==1
                if opt.selfig(2)
                        fig2=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
                        subplot(211)
                        plot(dates_full(sd:ed_full), tp_nb(sd:ed_full),'-','Color',[0.3 0.7 0.3],'LineWidth', 2,'DisplayName','10th');
                        title({figureTitleTag;'\bf{Probability of Switching from Normal to Bad Regime}'},'Interpreter','Latex')
                        axis tight
                        set(gca,'XTick',datenum(dates_full(sd:numticks:ed_full)))
                        datetick('x','yyyy','keepticks')
                        set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
                        hBands=recessionplot;
                        set(hBands,'FaceColor','k','FaceAlpha',0.1)
                        set(gca, 'FontName', 'Times New Roman');
                        set(gca, 'FontSize', TickSize);
                        set(gca,'Layer','top')
                        set(gca,'TickLabelInterpreter','Latex')

                        subplot(212)
                        plot(dates_full(sd:ed_full), tp_bn(sd:ed_full),'-','Color',[0.8 0.2 0.2],'LineWidth', 2,'DisplayName','10th');
                        title('\bf{Probability of Switching from Bad to Normal Regime}','interpreter','Latex')
                        axis tight
                        set(gca,'XTick',datenum(dates_full(sd:numticks:ed_full)))
                        datetick('x','yyyy','keepticks')
                        set(gca, 'XLim', [dates_full(sd), dates_full(ed_full)])
                        ylim([0 1])
                        hBands=recessionplot;
                        set(hBands,'FaceColor','k','FaceAlpha',0.1)
                        set(gca, 'FontName', 'Times New Roman');
                        set(gca, 'FontSize', TickSize);
                        set(gca,'Layer','top')
                        set(gca,'TickLabelInterpreter','Latex')
                        % Set figure appereance
                        set(fig2,'PaperOrientation','portrait');
                        set(fig2, 'PaperSize', figSize);
                        set(fig2, 'PaperUnits', 'inches');
                        set(fig2, 'Units','inches');
                        set(fig2, 'PaperPositionMode', 'auto');
                        set(fig2, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
                        tightfig;                   
                        if opt.saveit
                            print('-dpdf',fig2,[figfolder 'TransProb_' modelname ],'-bestfit');
                        end
                end
            end
   
            
            %%
            %======================================================================
            % Density plots
            %======================================================================
            if opt.showfig==1
                if opt.selfig(3)
                    linetype= {'--','-','-.',':','-.'};
                    for nn=1:opt.nvars
                        for ii=1:numel(opt.date_index)
                            fig3 = figure;
                            hold on
                            plot(squeeze(results.xi(ii,nn,:)),squeeze(results.pdfi(ii,nn,:)),linetype{2},'Color',[0 0 0 ],'LineWidth', 3,'DisplayName',[datestr(dates_full(opt.date_index(ii)),dataformat)]);
                            plot(squeeze(results.xi_good(ii,nn,:)),squeeze(results.pdfi_good(ii,nn,:)),linetype{4},'Color',colors(60,:),'LineWidth', 3,'DisplayName','Good regime');
                            plot(squeeze(results.xi_bad(ii,nn,:)),squeeze(results.pdfi_bad(ii,nn,:)),linetype{5},'Color',colors(1,:),'LineWidth', 3,'DisplayName','Bad regime');
                            if opt.dir_or_it==1
                                if nn==opt.idxg
                                    mytitle = ['\bf{Average ' opt.vlabels{nn} ' over Next ' num2str(opt.hh) ' Months}'];
                                else
                                    mytitle = ['\bf{' opt.vlabels{nn} ' Next Month}'];
                                end
                            else
                                mytitle = ['\bf{Average ' opt.vlabels{nn} ' over Next ' num2str(opt.hh) ' Months}'];
                            end
                            title({figureTitleTag;mytitle},'Interpreter','Latex');
                            if opt.dir_or_it==2
                                plot(zeros(1,100),zeros(1,100),'*','DisplayName',['Prob bad regime $(t+1,t+h|t)$ = ' num2str(round((1-results.quantiles.st_h_mean(opt.date_index(ii)))*100,1)) '\%']);
                            else
                                plot(zeros(1,100),zeros(1,100),'*','DisplayName',['Prob bad regime $(t+1,t+h|t)$ = ' num2str(round((1-results.quantiles.st_t_mean(opt.date_index(ii)))*100,1)) '\%']);
                            end
                            vline(0,'k--');
%                             vline(realized(opt.date_index(ii)),'r:','Realized Value');
                            legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
                            ylabel('PDF','fontsize',10,'interpreter','Latex')
                            xlabel(opt.units{nn},'fontsize',10,'Interpreter','Latex')
                            axis tight
                            set(gca, 'FontName', 'Times New Roman');
                            set(gca, 'FontSize', FontSize);
                            set(gca,'Layer','top')
                            set(gca,'TickLabelInterpreter','Latex')
                            % Set figure appereance
                            set(fig3,'PaperOrientation','portrait');
                            set(fig3, 'PaperSize', figSize);
                            set(fig3, 'PaperUnits', 'inches');
                            set(fig3, 'Units','inches');
                            set(fig3, 'PaperPositionMode', 'auto');
                            set(fig3, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
                            tightfig;
                            if opt.saveit
                                print('-dpdf',fig3,[figfolder 'Densities_' opt.varlist '_' datestr(dates_full(opt.date_index(ii)),dataformat) '_' modelname ],'-bestfit');
                            end
                        end
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
                if opt.selfig(4)
                    fig4=figure;
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
                    % Set figure appereance
                    set(fig4,'PaperOrientation','portrait');
                    set(fig4, 'PaperSize', figSize);
                    set(fig4, 'PaperUnits', 'inches');
                    set(fig4, 'Units','inches');
                    set(fig4, 'PaperPositionMode', 'auto');
                    set(fig4, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
                    tightfig;
                    if opt.saveit
                        print('-dpdf',fig4,[figfolder 'Fitted_' modelname ],'-bestfit');
                    end
                end
            end

            
            %% Model Fit
            clearvars y_fit realized
            y_fit = fit1.*p_reg1_smoothed + fit2.*p_reg2_smoothed;
            
            realized = db.GDPG.data ;
            realized = realized(opt.nlags+1:end);
            
            numticks=48;
            
            if opt.showfig==1
                if opt.selfig(5)
                    fig5=figure;
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
                    % Set figure appereance
                    set(fig5,'PaperOrientation','portrait');
                    set(fig5, 'PaperSize', figSize);
                    set(fig5, 'PaperUnits', 'inches');
                    set(fig5, 'Units','inches');
                    set(fig5, 'PaperPositionMode', 'auto');
                    set(fig5, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
                    tightfig;
                    if opt.saveit
                        print('-dpdf',fig5,[figfolder 'Resids_' modelname ],'-bestfit');
                    end
                end
            end

            
            %% Regime Probabilities
            if opt.showfig==1
                if opt.selfig(6)
                    fig6=figure;
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
                    % Set figure appereance
                    set(fig6,'PaperOrientation','portrait');
                    set(fig6, 'PaperSize', figSize);
                    set(fig6, 'PaperUnits', 'inches');
                    set(fig6, 'Units','inches');
                    set(fig6, 'PaperPositionMode', 'auto');
                    set(fig6, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
                    tightfig;
                    if opt.saveit
                        print('-dpdf',fig6,[figfolder 'ProbReg_' modelname ],'-bestfit');
                    end
                end
            end

            
            %%
            %======================================================================
            % Predictive score plot
            %======================================================================
            if opt.showfig==1
                if opt.selfig(7)
                    fig7=figure;
                    figSize7 = [12 6];
                    plot(dates_full(sd:ed_eval),results.ps(sd:ed_eval),'Color',[0.3 0.3 0.3],'LineWidth', 2.5);
                    title({figureTitleTag;'\bf{Predictive Score}'},'Interpreter','Latex');
                    axis tight
                    ylabel('Predictive Score','fontsize',10,'interpreter','Latex')
                    xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
                    datetick('x','yyyy','keepticks')
                    set(gca, 'XLim', [dates_full(sd), dates_full(ed_eval)])
                    ylim([0 1])
                    hBands=recessionplot;
                    set(hBands,'FaceColor','k','FaceAlpha',0.1)
                    set(gca, 'FontName', 'Times New Roman');
                    set(gca, 'FontSize', FontSize);
                    set(gca,'Layer','top')
                    set(gca,'TickLabelInterpreter','Latex')
                    % Set figure appereance
                    set(fig7,'PaperOrientation','portrait');
                    set(fig7, 'PaperSize', figSize7);
                    set(fig7, 'PaperUnits', 'inches');
                    set(fig7, 'Units','inches');
                    set(fig7, 'PaperPositionMode', 'auto');
                    set(fig7, 'Position', [figSize7(1)/5 figSize7(2)/5 figSize7(1) figSize7(2)]);
                    tightfig;
                    if opt.saveit
                       print('-dpdf',fig6,[figfolder 'PredScore_' modelname ],'-bestfit');
                    end
                end
            end
 
            
            %%
            %======================================================================
            % PITs CDF Plot
            %======================================================================
            if opt.selfig(8)
                fig8 = figure;
                figSize8 = [10 8]/1.5;
                h=cdfplot(results.pits);
                h.LineWidth = 2;
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
                title({figureTitleTag;'\bf{CDF of PITs}'},'FontSize',8,'Interpreter','Latex');
                hold off
                set(gca, 'FontName', 'Times New Roman');
                set(gca,'Layer','top')
                set(gca,'TickLabelInterpreter','Latex')
                % Set figure appereance
                set(fig8,'PaperOrientation','portrait');
                set(fig8, 'PaperSize', figSize8);
                set(fig8, 'PaperUnits', 'inches');
                set(fig8, 'Units','inches');
                set(fig8, 'PaperPositionMode', 'auto');
                set(fig8, 'Position', [figSize8(1)/3 figSize8(2)/5 figSize8(1) figSize8(2)]);
                tightfig;
                if opt.saveit
                    print('-dpdf',fig7,[figfolder 'PITsCDF_' modelname ],'-bestfit');
                end
            end
            
            
            %%
            %======================================================================
            % PITs over Time Plot
            %======================================================================
            if opt.selfig(9)
                fig9 = figure;
                figSize9 = [12 6];
                plot(dates_full(sd:ed_eval),results.pits(sd:ed_eval),'Color',[0.3 0.3 0.3],'LineWidth', 2.5);
                title({figureTitleTag;'\bf{Probability Integral Transform (PIT)}'},'FontSize',8,'Interpreter','Latex');
                axis tight
                ylabel('PITs','fontsize',10,'interpreter','Latex')
                xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
                datetick('x','yyyy','keepticks')
                set(gca, 'XLim', [dates_full(sd), dates_full(ed_eval)])
                ylim([0 1])
                hBands=recessionplot;
                set(hBands,'FaceColor','k','FaceAlpha',0.1)
                set(gca, 'FontName', 'Times New Roman');
                set(gca, 'FontSize', FontSize);
                set(gca,'Layer','top')
                set(gca,'TickLabelInterpreter','Latex')
                % Set figure appereance
                set(fig9,'PaperOrientation','portrait');
                set(fig9, 'PaperSize', figSize9);
                set(fig9, 'PaperUnits', 'inches');
                set(fig9, 'Units','inches');
                set(fig9, 'Position', [figSize9(1)/5 figSize9(2)/5 figSize9(1) figSize9(2)]);
                tightfig;
                if opt.saveit
                    print('-dpdf',fig9,[figfolder 'PITs_' modelname ],'-bestfit');
                end            
            end
            
        end
    end
end

fprintf('Removing path from RISE toolbox');
rise_exit
