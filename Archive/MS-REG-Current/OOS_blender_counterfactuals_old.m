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
opt.dir_or_it      = 2;     % 1 = direct, 2 = iterated
opt.transprob      = 1;     % 1 = endogenous, otherwise exogenous
opt.simtype        = 1;     % 1 = forecast, 0 = drawst
opt.hh             = 12;    % forecast horizon
opt.const          = 1;     % 1 = have a constant in transition probability
opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
opt.fastsim        = 1;     % 1 = fast simulation code, 0
opt.freshsim       = 1;     % 1 = run fresh simulations, 0 load existing if stored
opt.paramsUse      = 'mode';
opt.qw             = 0.05:0.05:0.95;            % Quantiles for QW CRPS test
opt.saveit         = 0;                         % 1 = save figures
opt.showfig        = 1;                         % 1 = show figures
opt.selfig         = [1 0 0 0 0 0 0 0 1 0 0];   % 1 = vector of figures


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
% mymodels = [101:108];
mymodels = [105,205]; 

% Select list of counterfactuals
mycounterfactuals = []; %[2 3]; %[2 3];

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
esample = 'dec 01 2019';
end_date_IS = datestr(esample,'yyyy-mmm');
maxwin = months(efirst,esample);
efirstok = datestr(efirst,dataformat);

% Dates for Plotting
start_plot = efirstok;
end_plot   = '2020-Nov';

% Dates for forecasting performance evaluation
end_eval   = '2018-Mar';


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
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' spectag '/' foldertag '/' ];
        IS_paramsfolder = ['results/Params/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags)  '/' spectag '/'];
        figfolder   = ['results/Figures/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/' foldertag '/'];
        
        % Folder to store simulations and counterfactuals
        simfolder = ['results/Simulations/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
        counterfolder = ['results/Counterfactuals/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
        performfolder = ['results/Perfomance/' opt.paramsUse '/Full/Direct/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
        
    elseif opt.dir_or_it ==2
        fcstype= 'iterated';
        varlist={'FF','MF','GDPG'}';
        
        % Folders containing parameters
        paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
        IS_paramsfolder = ['results/Params/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/'];
        
        % Folder to store figures
        figfolder   =  ['results/Figures/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' foldertag '/'];
        
        % Folder to store counterfactuals et al.
        simfolder = ['results/Simulations/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
        counterfolder = ['results/Counterfactuals/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
        performfolder = ['results/Perfomance/' opt.paramsUse '/Full/Iterated/VAR' num2str(nlags) '/' spectag '/' foldertag '/' ];
    end
    
    %% Create Folders
    if exist(figfolder,'dir')==0;  mkdir(figfolder); end
    if exist(simfolder,'dir')==0;  mkdir(simfolder); end
    if exist(counterfolder,'dir')==0;  mkdir(counterfolder); end
    if exist(performfolder,'dir')==0;  mkdir(performfolder); end
    
    
    %% Read data
    rawdata = readtable(['data/' datafilename '.xlsx'],'Sheet',sheetuse,'ReadVariableNames',true);
    Dates    = rawdata.Dates;
    
    e_IS=find(Dates==datetime(end_date_IS,'InputFormat','yyyy-MMM'));
    
    % Select enddb based on model
    % Makes sure we forecast from the same origin
    enddb_IS   = [num2str(year(datetime(Dates(e_IS-opt.hh+1),'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(Dates(e_IS-opt.hh+1),'InputFormat','yyyy-MMM')))];
    
    
    %% Collect in-sample quantiles
    if opt.hh ==1
        IS_modelname = [datafilename '_' sheetuse '_' '1973M1' enddb_IS ...
            '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];
    elseif opt.hh ==12
        IS_modelname = [datafilename '_' sheetuse '_' '1973M1' enddb_IS ...
            '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];
    end
    
    IS_results = load([IS_paramsfolder IS_modelname '_results.mat']);
    IS_results = IS_results.results;
    
    %% Collect data for plots
    
    end_date_IS = datestr(esample,dataformat);
    
    % Collect Full Database
    [db, db_full,startdb,enddb,tex] = fLoadData(datafilename,sheetuse,start_date,end_date_IS,opt);
    
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
    
    
    %%
    
    % Vector of dates for the full sample
    dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_plot,'InputFormat',inputformat)))';
    
    % Dates from the start of estimation sample
    % e.g: Dec-1999 to Dec-2019
    dates_oos = datenum((datetime(efirstok,'InputFormat',inputformat)+calmonths(nlags)):calmonths(1):(datetime(end_date_IS,'InputFormat',inputformat)))';
    
    % Indexes in OOS estimation sample
    sd_oos      = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates_oos);
    ed_oos      = find(datenum(end_date_IS,dataformat)==dates_oos);
    ed_oos_eval = find(datenum(end_eval,dataformat)==dates_oos);
    ed_oos_varp = find(datenum((datetime(end_date_IS,'InputFormat',inputformat)-calmonths(opt.hh-1)))==dates_oos);
    
    % Indexes in Full data
    sd_full      = find(datenum((datetime(start_plot,'InputFormat',inputformat))+calmonths(nlags))==dates_full);
    ed_full      = find(datenum((datetime(end_date_IS,'InputFormat',inputformat)))==dates_full);
    ed_full_eval = find(datenum(end_eval,dataformat)==dates_full);
    ed_full_varp = find(datenum((datetime(end_date_IS,'InputFormat',inputformat)-calmonths(opt.hh-1)))==dates_full);
    
    
    %% Blend OOS results
    
    for i=1:maxwin
        
        % Information for loading stored estimation results
        end_date = datestr(datenum(datetime(efirstok,'InputFormat',inputformat)+calmonths(i)),dataformat);
        startdb = [num2str(year(datetime(start_date,'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(start_date,'InputFormat','yyyy-MMM')))];
        e=find(Dates==datetime(end_date,'InputFormat','yyyy-MMM'));
        
        % Select enddb based on model
        % Makes sure we forecast from the same origin
        enddb   = [num2str(year(datetime(Dates(e-opt.hh+1),'InputFormat','yyyy-MMM'))) 'M' num2str(month(datetime(Dates(e-opt.hh+1),'InputFormat','yyyy-MMM')))];
        
        % Get index of end_date. This is the forecast origin.
        myindex = find(datenum((datetime(end_date,'InputFormat',inputformat)))==dates_full);
        
        % Model name and title
        % Note: iterated model posterior mode depends on opt.hh through the estimation sample
        modelname = [datafilename '_' sheetuse '_' startdb enddb ...
            '_' 'C' num2str(opt.const) 'N' num2str(opt.normal) 'H' num2str(opt.hh) spectag tag];
        
        %==================
        % Load parameters
        %==================
        x = load([paramsfolder modelname '_results.mat']);
        temp = x.results;
        
        % Collect parameters
        sv = temp.sv;
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
        fprintf('\n Model estimated through end_date: %s \n',end_date);
        
        % Define target periods for density cuts
        target_dates = datenum(datetime(end_date,'InputFormat',inputformat));
        
        % Simulate through the end of the estimation sample.
        opt.tperiods = find(ismember(dates_full, target_dates)); %length(dates_full);
        
        % Extract the index for the desired density plots
        opt.date_index = find(ismember(dates_full, target_dates));
        
        %=============================
        % Simulate Predictive Density
        %=============================
        
        % Name of file with stored simulations
        sim_file = [simfolder modelname '_results' '.mat'];
        
        
        if exist(sim_file,'file') && opt.freshsim == 0
            
            % Load stored predictive density simulations
            load(sim_file);
            
        else
            % Simulate predictive density
            
            if ~opt.fastsim
                % Simulate from 1973M1:t+h
                
                if opt.dir_or_it==1
                    % Direct forecast
                    [quantiles,pdensity] = fPredDensityDirectFull(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt,simtype);
                elseif opt.dir_or_it==2
                    % Iterated forecast
                    [quantiles, pdensity] = fPredDensityIteratedFull(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,simtype);
                end
                
                % Collect Quantiles into Results Structure
                results.quantiles.dYsim_25  = quantiles.dYsim_25(opt.date_index);
                results.quantiles.dYsim_75  = quantiles.dYsim_75(opt.date_index);
                results.quantiles.dYsim_10  = quantiles.dYsim_10(opt.date_index);
                results.quantiles.dYsim_90  = quantiles.dYsim_90(opt.date_index);
                results.quantiles.mean      = quantiles.mean(opt.date_index);
                results.quantiles.st_t_mean = quantiles.st_t_mean(opt.date_index);
                
                % Parameters
                results.outparams = scriptParams_FULL_companion(sv,params_in,1);
                
                % Predictive density at specified target-date
                [results.pdfi,results.xi]  = ksdensity(pdensity.realized);
                [results.pdfi_good,results.xi_good]  = ksdensity(pdensity.good);
                [results.pdfi_bad,results.xi_bad]  = ksdensity(pdensity.bad);
                [results.cxi,results.cdfi] = ksdensity(pdensity.realized,'Function','icdf');
                
                
            elseif opt.fastsim
                
                % Fast simulation: simulate from t:t+h
                cdatesplit = strsplit(enddb,'M');
                cdate = datetime(str2num(cdatesplit{1}),str2num(cdatesplit{2}),1);
                usedate = find(ismember(dates_full, datenum(cdate)));
                
                if opt.dir_or_it==1
                    
                    [~,~,~,f]=filter(sv);
                    p_reg2 = f.updated_regime_probabilities.regime2.data;
                    
                    %                   (a)
                    %                     FF_in = FF_full(usedate-opt.hh+1:usedate);
                    %                     MF_in = MF_full(usedate-opt.hh+1:usedate);
                    %                     YY_in = GDPGH_full(usedate-opt.hh+1:usedate);
                    %                     TR_in = TR_full(usedate-opt.hh+1:usedate);
                    %                     PB_in = p_reg2_filtered(usedate-opt.hh+1);
                    
                    %                   (b)
                    FF_in = FF_full(usedate:usedate+opt.hh-1);
                    MF_in = MF_full(usedate:usedate+opt.hh-1);
                    YY_in = GDPGH_full(usedate:usedate+opt.hh-1);
                    TR_in = TR_full(usedate:usedate+opt.hh-1);
                    PB_in = p_reg2(usedate);
                    
                    % Direct forecast
                    [quantiles,pdensity] = fPredDensityDirectFullFast(sv,params_in,FF_in,MF_in,YY_in,TR_in,opt,simtype,PB_in);
                    
                elseif opt.dir_or_it==2
                    
                    [~,~,~,f]=filter(sv);
                    p_reg2 = f.updated_regime_probabilities.regime2.data;
                    opt.date_t = myindex;
                    
                    FF_in = FF_full;
                    MF_in = MF_full;
                    YY_in = GDPGH_full;
                    TR_in = TR_full;
                    PB_in = p_reg2(usedate);
                    
                    
                    % Iterated forecast
                    [quantiles, pdensity] = fPredDensityIteratedFullFast(sv,params_in,FF_in,MF_in,YY_in,TR_in,opt,simtype,PB_in);
                    
                end
                
                % Collect Quantiles into Results structure
                results.quantiles.dYsim_25  = quantiles.dYsim_25;
                results.quantiles.dYsim_75  = quantiles.dYsim_75;
                results.quantiles.dYsim_10  = quantiles.dYsim_10;
                results.quantiles.dYsim_90  = quantiles.dYsim_90;
                results.quantiles.mean      = quantiles.mean;
                results.quantiles.st_t_mean = quantiles.st_t_mean;
                results.qwvec               = quantiles.qwvec;
                
                % Collect parameters
                results.outparams = scriptParams_FULL_companion(sv,params_in,1);
                
                % Predictive density at specified target-date
                [results.pdfi,results.xi]  = ksdensity(pdensity.realized);
                [results.pdfi_good,results.xi_good]  = ksdensity(pdensity.good);
                [results.pdfi_bad,results.xi_bad]  = ksdensity(pdensity.bad);
                [results.cxi,results.cdfi] = ksdensity(pdensity.realized,'Function','icdf');
                
                
            end
            
            %Save simulation results
            save([simfolder modelname '_results.mat'],'results');
            
        end
        
        %%
        %=============================
        % COLLECT RESULTS
        %=============================
        
        % Pass results corresponding to period-i to structure keeping all results
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
        
        %=============================================
        % Quantile-weighted CRPS at specified episodes
        % Comment: Would be preferrable to store quantiles for qw vector
        %=============================================
        
        qw = opt.qw;
        J = length(qw);
        for j=1:J
            a = qw(j);
            
            % Get quantile from cdf
            % [~,qi] = min(abs(results.cdfi-a));
            % qf = results.cxi(qi);
            
            % Quantile from simulated predictive density
            qf = results.qwvec(j);
            
            if  GDPGH_full_wt(myindex,:)<=qf
                qs = (1-a)*(qf - GDPGH_full_wt(myindex,:));
            else
                qs = -a*(qf - GDPGH_full_wt(myindex,:));
            end
            u_qwcrps(j)      = qs;              % Uniform
            c_qwcrps(j)      = (a*(1-a))*qs;    % Center
            t_qwcrps(j)      = ((2*a-1)^2)*qs;  % Tails
            rt_qwcrps(j)     = ((a)^2)*qs;      % Right Tail
            lt_qwcrps(j)     = ((1-a)^2)*qs;    % Left Tail
            if a<0.5
                rtpure_qwcrps(j) = 0*qs;      % Right Tail (above 0.5)
                ltpure_qwcrps(j) = ((0.5-a)^2)*qs;    % Left Tail (below 0.5)
            else
                rtpure_qwcrps(j) = ((a-0.5)^2)*qs;      % Right Tail (above 0.5)
                ltpure_qwcrps(j) = 0*qs;    % Left Tail (below 0.5)
            end
        end
        uni_qwcrps_time(i)           = mean(u_qwcrps);
        center_qwcrps_time(i)        = mean(c_qwcrps);
        tails_qwcrps_time(i)         = mean(t_qwcrps);
        righttail_qwcrps_time(i)     = mean(rt_qwcrps);
        lefttail_qwcrps_time(i)      = mean(lt_qwcrps);
        righttailpure_qwcrps_time(i) = mean(rtpure_qwcrps);
        lefttailpure_qwcrps_time(i)  = mean(ltpure_qwcrps);
        
        %%
        %=============================
        % COLLECT RESULTS
        %=============================
        
        % Pass results corresponding to period-i to structure keeping all results
        forcperf.ps{i}                        = ps(i);                        % Predictive Score (t)
        forcperf.pits{i}                      = pits(i);                      % PITs(t)
        forcperf.uni_qwcrps_time{i}           = uni_qwcrps_time(i);           % QW-CRPS(t) - Uniform
        forcperf.center_qwcrps_time{i}        = center_qwcrps_time(i);        % QW-CRPS(t) - Center
        forcperf.tails_qwcrps_time{i}         = tails_qwcrps_time(i);         % QW-CRPS(t) - Tail
        forcperf.righttail_qwcrps_time{i}     = righttail_qwcrps_time(i);     % QW-CRPS(t) - Right Tail
        forcperf.lefttail_qwcrps_time{i}      = lefttail_qwcrps_time(i);      % QW-CRPS(t) - Left Tail
        forcperf.righttailpure_qwcrps_time{i} = righttailpure_qwcrps_time(i); % QW-CRPS(t) - Right Tail (above 0.5)
        forcperf.lefttailpure_qwcrps_time{i}  = lefttailpure_qwcrps_time(i);  % QW-CRPS(t) - Left Tail  (below 0.5)
        
        
        
        
        % Store benchmark results
        save([performfolder modelname '.mat'],'forcperf');
        
        
        %%
        %==================================
        % Counterfactual densities
        %==================================
        if ~isempty(mycounterfactuals)
            
            % Simulate for all counterfactuals
            for j=1:length(mycounterfactuals)
                
                opt_temp = opt;
                opt_temp.counterfactual = mycounterfactuals(j);
                
                cf_file = [counterfolder modelname '_cf' num2str(mycounterfactuals(j)) '.mat'];
                
                if exist(cf_file,'file') && opt.freshsim == 0
                    
                    load(cf_file);
                    
                else
                    
                    % FAST SIMULATION
                    cdatesplit = strsplit(enddb,'M');
                    cdate = datetime(str2num(cdatesplit{1}),str2num(cdatesplit{2}),1);
                    usedate = find(ismember(dates_full, datenum(cdate)));
                    
                    % Simulate model
                    if opt.dir_or_it==1
                        
                        [~,~,~,f]=filter(sv);
                        p_reg2 = f.updated_regime_probabilities.regime2.data;
                        
                        % Direct forecast
                        % % (a)
                        % FF_in = FF_full(usedate-opt.hh+1:usedate);
                        % MF_in = MF_full(usedate-opt.hh+1:usedate);
                        % YY_in = GDPGH_full(usedate-opt.hh+1:usedate);
                        % TR_in = TR_full(usedate-opt.hh+1:usedate);
                        % PB_in = p_reg2_filtered(usedate-opt.hh+1);
                        
                        %                     (b)
                        FF_in = FF_full(usedate:usedate+opt.hh-1);
                        MF_in = MF_full(usedate:usedate+opt.hh-1);
                        YY_in = GDPGH_full(usedate:usedate+opt.hh-1);
                        TR_in = TR_full(usedate:usedate+opt.hh-1);
                        PB_in = p_reg2(usedate);
                        
                        % Direct forecast
                        [quantiles_cf,pdensity_cf] = fPredDensityDirectFullFast(sv,params_in,FF_in,MF_in,YY_in,TR_in,opt_temp,simtype,PB_in);
                        
                    elseif opt.dir_or_it==2
                        
                        [~,~,~,f]=filter(sv);
                        p_reg2 = f.updated_regime_probabilities.regime2.data;
                        opt.date_t = myindex;
                        
                        FF_in = FF_full;
                        MF_in = MF_full;
                        YY_in = GDPGH_full;
                        TR_in = TR_full;
                        PB_in = p_reg2(usedate);
                        
                        % Iterated forecast
                        [quantiles, pdensity] = fPredDensityIteratedFullFast(sv,params_in,FF_in,MF_in,YY_in,TR_in,opt_temp,simtype,PB_in);
                        
                    end
                    
                    % Quantiles
                    cfresults.quantiles.dYsim_25  = quantiles_cf.dYsim_25;
                    cfresults.quantiles.dYsim_75  = quantiles_cf.dYsim_75;
                    cfresults.quantiles.dYsim_10  = quantiles_cf.dYsim_10;
                    cfresults.quantiles.dYsim_90  = quantiles_cf.dYsim_90;
                    cfresults.quantiles.mean      = quantiles_cf.mean;
                    cfresults.quantiles.st_t_mean = quantiles_cf.st_t_mean;
                    
                    % Predictive density at specified episodes
                    [cfresults.pdfi,cfresults.xi]           = ksdensity(pdensity_cf.realized);
                    [cfresults.pdfi_good,cfresults.xi_good] = ksdensity(pdensity_cf.good);
                    [cfresults.pdfi_bad,cfresults.xi_bad]   = ksdensity(pdensity_cf.bad);
                    [cfresults.cxi,cfresults.cdfi]          = ksdensity(pdensity_cf.realized,'Function','icdf');
                    
                    save([counterfolder modelname '_cf' num2str(mycounterfactuals(j)) '.mat'],'cfresults');
                    
                end
                
                % Map results into object for plotting
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
    % Plots
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
            fig1=figure('Units','normalized','Position', [0,0,1.1,0.9]); clf;
            subplot(211)
            hold on
            l1=plot(dates_oos(sd_oos:ed_oos), benchmark.dYsim_10(sd_oos:ed_oos),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th out-of-sample');
            l2=plot(dates_oos(sd_oos:ed_oos), benchmark.dYsim_90(sd_oos:ed_oos),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th out-of-sample');
            l3=plot(dates_full(sd_full:ed_full), IS_results.quantiles.dYsim_10(sd_full:ed_full),'-.','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th in-sample');
            l4=plot(dates_full(sd_full:ed_full), IS_results.quantiles.dYsim_90(sd_full:ed_full),'-.','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th in-sample');
            %
            legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            l5=plot(dates_full(sd_full:ed_full_varp),GDPGH_full_wt(sd_full:ed_full_varp),'-.','Color','k','LineWidth', 2,'DisplayName','Realized Value');
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            
            if opt.dir_or_it==1
                title(['Quantiles Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            else
                title(['Quantiles Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            end
            axis tight
            numticks = 12;
            set(gca,'XTick',datenum(dates_full(sd_full:numticks:ed_full)))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_full(sd_full), dates_full(ed_full)])
            
            if opt.hh==12
                ylim([-10 8])
            else
                ylim([-30 20])
            end
            
            subplot(212)
            hold on
            l3=plot(dates_oos(sd_oos:ed_oos), 1-benchmark.st_t_mean(sd_oos:ed_oos),'-','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Simulated $s(t)$');
            l4=plot(dates_full(sd_full:ed_full), 1-IS_results.quantiles.st_t_mean(sd_full:ed_full),'-.','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Simulated s(t) in-sample');
            legend([l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            hold off
            ylabel('Percent','interpreter','Latex','fontsize',10)
            title('Probability of Bad Regime','FontSize',16','Interpreter','Latex');
            axis tight
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos)))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
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
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos)))
            lw = 2;
            fs = 12;
            
            subplot(3,1,1)
            g = A0_sync_1(:,2,1);
            b = A0_sync_2(:,2,1);
            plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
            hold on
            plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
            hold off
            legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
            title('$A_O(2,1)$','FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
            
            subplot(3,1,2)
            g = A0_sync_1(:,3,1);
            b = A0_sync_2(:,3,1);
            plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
            hold on
            plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
            hold off
            legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
            title('$A_O(3,1)$','FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
            
            subplot(3,1,3)
            g = A0_sync_1(:,3,2);
            b = A0_sync_2(:,3,2);
            plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
            hold on
            plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
            hold off
            legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
            title('$A_O(3,2)$','FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
            
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
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos)))
            lw = 2;
            fs = 12;
            
            count = 0;
            for k=1:size(A1_sync_1,2)
                for l=1:size(A1_sync_1,3)
                    count = count + 1;
                    subplot(size(A1_sync_1,2),size(A1_sync_1,3),count)
                    g = A1_sync_1(:,k,l);
                    b = A1_sync_2(:,k,l);
                    plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
                    hold on
                    plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
                    hold off
                    legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
                    title(['$A_1(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
                    axis tight
                    datetick('x','yyyy','keepticks')
                    set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
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
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos)))
            lw = 2;
            fs = 12;
            
            count = 0;
            for k=1:size(C_sync_1,2)
                count = count + 1;
                subplot(size(C_sync_1,2),1,count)
                g = C_sync_1(:,k,1);
                b = C_sync_2(:,k,1);
                plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
                hold on
                plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
                hold off
                legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
                title(['$C(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
                axis tight
                datetick('x','yyyy','keepticks')
                set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
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
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos)))
            lw = 2;
            fs = 12;
            
            count = 0;
            for k=1:size(SIG_sync_1,2)
                count = count + 1;
                subplot(size(SIG_sync_1,2),1,count)
                g = SIG_sync_1(:,k,k);
                b = SIG_sync_2(:,k,k);
                plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
                hold on
                plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
                hold off
                legend('Good','Bad','box','off','FontSize',fs-2,'Interpreter','Latex')
                title(['$\Sigma(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
                axis tight
                datetick('x','yyyy','keepticks')
                set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
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
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos)))
            lw = 2;
            fs = 12;
            
            subplot(3,1,1)
            g = a12;
            b = a21;
            plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
            hold on
            plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
            hold off
            legend('good to bad','bad to good','box','off','FontSize',fs-2,'Interpreter','Latex')
            title(['$a(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
            
            subplot(3,1,2)
            g = b12;
            b = b21;
            plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
            hold on
            plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
            hold off
            legend('good to bad','bad to good','box','off','FontSize',fs-2,'Interpreter','Latex')
            title(['$b(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
            
            subplot(3,1,3)
            g = c12;
            b = c21;
            plot(dates_oos(sd_oos:ed_oos),g(sd_oos:ed_oos),'LineWidth', lw)
            hold on
            plot(dates_oos(sd_oos:ed_oos),b(sd_oos:ed_oos),'LineWidth', lw)
            hold off
            legend('good to bad','bad to good','box','off','FontSize',fs-2,'Interpreter','Latex')
            title(['$c(' num2str(k) ',' num2str(l) ')$'],'FontSize',fs,'Interpreter','Latex')
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
            
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
            lw = 2;
            fs = 12;
            
            % Do not use last hh observations since no data for GDPG
            plot(dates_oos(sd_oos:ed_oos_eval),ps(sd_oos:ed_oos_eval),'LineWidth', lw)
            hold on
            plot(dates_full(sd_full:ed_full_eval),IS_results.ps(sd_full:ed_full_eval),'-.','LineWidth', lw)
            hold off
            if opt.dir_or_it==1
                title(['Predictive Score Direct: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            else
                title(['Predictive Score Iterated: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            end
            axis tight
            ylim([0 1])
            hBands=recessionplot;
            legend('OOS','IS','Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos_eval)))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos_eval)])
            
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
            lw = 2;
            fs = 12;
            
            % Do not use last hh observations since no data for GDPG
            plot(dates_oos(sd_oos:ed_oos_eval),pits(sd_oos:ed_oos_eval),'LineWidth', lw)
            hold on
            plot(dates_full(sd_full:ed_full_eval),IS_results.pits(sd_full:ed_full_eval),'-.','LineWidth', lw)
            hold off
            legend('OOS','IS','Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            if opt.dir_or_it==1
                title(['PITs Direct: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            else
                title(['PITs Iterated: ' figureTitleTag ],'FontSize',fs,'Interpreter','Latex');
            end
            axis tight
            datetick('x','yyyy','keepticks')
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos_eval)))
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos_eval)])
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
            pits = pits(sd_oos:ed_oos_eval);
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
            l1=plot(dates_oos(sd_oos:end),benchmark.dYsim_10(sd_oos:ed_oos),'-','LineWidth',2,'Color',colors(10,:),'DisplayName','10th out-of-sample'); hold on;
            l2=plot(dates_oos(sd_oos:end),benchmark.dYsim_90(sd_oos:ed_oos),'-','LineWidth',2,'Color',colors(50,:),'DisplayName','90th out-of-sample');
            l3=plot(dates_oos(sd_oos:end),cf{1}.dYsim_10(sd_oos:ed_oos),'--','LineWidth',2,'Color',colors(10,:),'DisplayName','Couterfactual 1');
            l4=plot(dates_oos(sd_oos:end),cf{1}.dYsim_90(sd_oos:ed_oos),'--','LineWidth',2,'Color',colors(50,:),'DisplayName','Couterfactual 1');
            legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            if opt.dir_or_it==1
                title(['OOS Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            else
                title(['OOS Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            end
            axis tight
            numticks = 36;
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:(end))))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(end)])
            
            subplot(2,1,2)
            l1=plot(dates_oos(sd_oos:end),benchmark.dYsim_10(sd_oos:ed_oos),'-','LineWidth',2,'Color',colors(10,:),'DisplayName','10th out-of-sample'); hold on;
            l2=plot(dates_oos(sd_oos:end),benchmark.dYsim_90(sd_oos:ed_oos),'-','LineWidth',2,'Color',colors(50,:),'DisplayName','90th out-of-sample');
            l3=plot(dates_oos(sd_oos:end),cf{2}.dYsim_10(sd_oos:ed_oos),'--','LineWidth',2,'Color',colors(10,:),'DisplayName','Couterfactual 2');
            l4=plot(dates_oos(sd_oos:end),cf{2}.dYsim_90(sd_oos:ed_oos),'--','LineWidth',2,'Color',colors(50,:),'DisplayName','Couterfactual 2');
            legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
            set(gca,'children',flipud(get(gca,'children')))
            if opt.dir_or_it==1
                title(['OOS Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            else
                title(['OOS Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
            end
            axis tight
            numticks = 36;
            set(gca,'XTick',datenum(dates_oos(sd_oos:numticks:ed_oos)))
            datetick('x','yyyy','keepticks')
            set(gca, 'XLim', [dates_oos(sd_oos), dates_oos(ed_oos)])
            
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
            for ix = [92 101 240 241]
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
        
        %% Table 1: Violation Ratios
        
        % Do not consider last 12 observations since no data for GDPG
        stot = length(sd_oos:ed_oos_eval);
        snum_10 = sum(GDPGH_full_wt(sd_full:ed_full_eval)<benchmark.dYsim_10(sd_oos:ed_oos_eval));
        vratio_10 = snum_10/stot*100;
        snum_25 = sum(GDPGH_full_wt(sd_full:ed_full_eval)<benchmark.dYsim_25(sd_oos:ed_oos_eval));
        vratio_25 = snum_25/stot*100;
        snum_75 = sum(GDPGH_full_wt(sd_full:ed_full_eval)<benchmark.dYsim_75(sd_oos:ed_oos_eval));
        vratio_75 = snum_75/stot*100;
        snum_90 = sum(GDPGH_full_wt(sd_full:ed_full_eval)<benchmark.dYsim_90(sd_oos:ed_oos_eval));
        vratio_90 = snum_90/stot*100;
        
        snum_IS_10 = sum(GDPGH_full_wt(sd_full:ed_full_eval)<IS_results.quantiles.dYsim_10(sd_full:ed_full_eval));
        vratio_IS_10 = snum_IS_10/stot*100;
        snum_IS_25 = sum(GDPGH_full_wt(sd_full:ed_full_eval)<IS_results.quantiles.dYsim_25(sd_full:ed_full_eval));
        vratio_IS_25 = snum_IS_25/stot*100;
        snum_IS_75 = sum(GDPGH_full_wt(sd_full:ed_full_eval)<IS_results.quantiles.dYsim_75(sd_full:ed_full_eval));
        vratio_IS_75 = snum_IS_75/stot*100;
        snum_IS_90 = sum(GDPGH_full_wt(sd_full:ed_full_eval)<IS_results.quantiles.dYsim_90(sd_full:ed_full_eval));
        vratio_IS_90 = snum_IS_90/stot*100;
        
        % Display table 1
        dc = 1;
        T = table([round(vratio_10,dc);round(vratio_IS_10,dc)],...
            [round(vratio_25,dc);round(vratio_IS_25,dc)],...
            [round(vratio_75,dc);round(vratio_IS_75,dc)],...
            [round(vratio_90,dc);round(vratio_IS_90,dc)],'VariableNames',{'10th','25th','75th','90th'},'RowName',{'OOS','IS'});
        T = table(T,'VariableNames',{'Violation Ratios'});
        disp(T)
        
        %% Table 2: Quantile-Weighted CRPS (Gneiting and Ranjan, JBES 2011)
        
        % Time-average of QW-CRPS
        uqwcrps = mean(uni_qwcrps_time);
        cqwcrps = mean(center_qwcrps_time);
        tqwcrps = mean(tails_qwcrps_time);
        rtqwcrps = mean(righttail_qwcrps_time);
        ltqwcrps = mean(lefttail_qwcrps_time);
        
        % Display table 2
        dc = 2;
        T2 = table([round(uqwcrps,dc);round(cqwcrps,dc);...
            round(tqwcrps,dc);round(rtqwcrps,dc);...
            round(ltqwcrps,dc)],'VariableNames',{'QW-CRPS'},'RowName',{'Uniform','Center','Tail','Right Tail','Left Tail'});
        T2 = table(T2,'VariableNames',{'Quantile-Weighted CRPS'});
        disp(T2)
        
        
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
