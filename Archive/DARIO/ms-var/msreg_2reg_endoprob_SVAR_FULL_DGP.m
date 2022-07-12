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
opt.dir_or_it      = 2;     % 1 = direct, 2 = iterated
opt.hh             = 12;    % forecast horizon
opt.const          = 1;     % 1 = have a constant in transition probability
opt.normal         = 0;     % 1 = use normal distribution, 0 = gamma distribution
opt.paramsUse      = 'mode';
opt.force_estimate = 1; % 1 = force to reestimate model

% Posterior mode options
opt.nDraws       = 2000;    % Number of simulations to produce quantiles
opt.nParamDraws  = 1;       % Draws from posterior. has to be less thatn options.N;
opt.nSamples     = 5;       % Number of simulated DGP samples
opt.nSim         = 1200;
opt.tperiods     = opt.nSim-opt.hh;
opt.burn         = 500;
opt.saveit       = 1;     % 1 = save figures


SEED = 1;


% Load estimated parameters from
modeldata  = '11302020_US_DFM_1973M12019M10_C1N0H12CT_A0T_A1_SIGT';
paramsdata = ['results/Params/' opt.paramsUse '/Full/Iterated/'];

modelspec  = 2;
imodel     = modelspec;

modeldgp    = 7;

tag = [''];


use_dgp = 0;

%=========================================
% Load stored results from data estimation
%=========================================
if use_dgp==0
    load([paramsdata modeldata '.mat']);
    % Map objects:
    sv_data    = sPMode.sv;
    pmode_data = posterior_mode(sv_data);
    params_in  = pmode_data;
elseif use_dgp==1
    %A0 Matrix (Lower triangular elements)
    dgp.a0_2_1_sync_1 = 0.15;
    dgp.a0_3_1_sync_1 = 0.10;       % response of GDP to FF
    dgp.a0_3_2_sync_1 = -0.2;       % response of GDP to MF
    
    dgp.a0_2_1_sync_2 = 0.25;
    dgp.a0_3_1_sync_2 = 0.3;        % response of GDP to FF
    dgp.a0_3_2_sync_2 = -0.5;       % response of GDP to MF
    
    % A1 Matrix S=1
    dgp.a1_1_1_sync_1 =  0.40;      % Persistence FF
    dgp.a1_1_2_sync_1 = -0.05;
    dgp.a1_1_3_sync_1 =  0;
    
    dgp.a1_2_1_sync_1 = -0.10;
    dgp.a1_2_2_sync_1 =  0.40;      % Persistence MF
    dgp.a1_2_3_sync_1 =  0.0;       % response of MF to GDPG(-1)
    
    dgp.a1_3_1_sync_1 = 0.01113;
    dgp.a1_3_2_sync_1 = 0.0298;
    dgp.a1_3_3_sync_1 = 0.10; %Persistence GDP
    
    % A1 Matrix S=2
    dgp.a1_1_1_sync_2 =  0.89; % Persistence FF
    dgp.a1_1_2_sync_2 = -0.05;
    dgp.a1_1_3_sync_2 =  0.0;
    
    dgp.a1_2_1_sync_2 = -0.10;
    dgp.a1_2_2_sync_2 =  0.70;  % Persistence MF
    dgp.a1_2_3_sync_2 =  0.0;   % Response of FF to GDPG(-1)
    
    dgp.a1_3_1_sync_2 = 0.0143;
    dgp.a1_3_2_sync_2 = 0.0508;
    dgp.a1_3_3_sync_2 = 0.20; %Persistence GDP
    
    % C Matrix
    dgp.c_1_1_sync_1 = 0;
    dgp.c_2_1_sync_1 = 0;
    dgp.c_3_1_sync_1 = 0.15;
    
    dgp.c_1_1_sync_2 = 0.05;
    dgp.c_2_1_sync_2 = -.05;
    dgp.c_3_1_sync_2 = -0.05;
    
    % SIG
    dgp.s_1_1_sync_1 = 1.0;
    dgp.s_1_1_sync_2 = 1.0;
    dgp.s_2_2_sync_1 = 1.0;
    dgp.s_2_2_sync_2 = 1.0;
    dgp.s_3_3_sync_1 = 0.5;
    dgp.s_3_3_sync_2 = 1.75;
    
    % Transition probabilities
    dgp.a12 = 3.8;
    dgp.a21 = 2.1;
    dgp.b12 = 0.35;
    
    dgp.b21 = 0.4;
    dgp.c12 = 0.1;
    dgp.c21 = 0.1;
    
    % ZERO AND NORMALIZATION RESTRICTIONS
    
    % Upper triangular portion of A0
    dgp.a0_1_2_sync_1 = 0;
    dgp.a0_1_2_sync_2 = 0;
    dgp.a0_1_3_sync_1 = 0;
    dgp.a0_1_3_sync_2 = 0;
    dgp.a0_2_3_sync_1 = 0;
    dgp.a0_2_3_sync_2 = 0;
    
    % Normalizations
    dgp.a0_1_1_sync_1 = 1;
    dgp.a0_1_1_sync_2 = 1;   % FF contemporaneous coefficients
    dgp.a0_2_2_sync_1 = 1;
    dgp.a0_2_2_sync_2 = 1;   % MF contemporaneous coefficients
    dgp.a0_3_3_sync_1 = 1;
    dgp.a0_3_3_sync_2 = 1;   % GDP contemporaneous coefficients
    
    
    %% Alternative DGPs
    
    if modeldgp==2
        % Transition probabilities
        dgp.a12 = 3.8;
        dgp.b12 = 0.35;
        dgp.c12 = 0.1;
        
        dgp.a21 = 1.0;
        dgp.b21 = 0.1;
        dgp.c21 = 0.1;
        
    elseif modeldgp==3
        % SIG
        dgp.s_1_1_sync_1 = 1.0;
        dgp.s_1_1_sync_2 = 1.0;
        dgp.s_2_2_sync_1 = 1.0;
        dgp.s_2_2_sync_2 = 1.0;
        dgp.s_3_3_sync_1 = 0.5;
        dgp.s_3_3_sync_2 = 1.0;
        
    elseif modeldgp==4
        % Transition probabilities
        % This calibration produces data-like transition probs:
        dgp.a12 = 3.8;
        dgp.b12 = 0.45;
        dgp.c12 = 0.45;
        
        dgp.a21 = 0.5;
        dgp.b21 = 0.1;
        dgp.c21 = 0.1;
        
    elseif modeldgp==5
        % Transition probabilities
        % Similar to gdp=4, with time variation in variances of MF and FF
        dgp.a12 = 3.8;
        dgp.b12 = 0.45;
        dgp.c12 = 0.45;
        
        dgp.a21 = 0.075;
        dgp.b21 = 0.05;
        dgp.c21 = 0.05;
        
        % Variances FF and MF
        dgp.s_1_1_sync_1 = 1.0;
        dgp.s_1_1_sync_2 = 2.5;
        dgp.s_2_2_sync_1 = 1.0;
        dgp.s_2_2_sync_2 = 2.5;
        
        % Variance of GDP
        dgp.s_3_3_sync_1 = 0.5;
        dgp.s_3_3_sync_2 = 0.75;
        
        
        % A0 Matrix S=1
        dgp.a0_2_1_sync_1 = 0.15;
        
        % A0 Matrix S=2
        dgp.a0_2_1_sync_2 = 0.35;
        
        % A1 Matrix S=1
        dgp.a1_1_1_sync_1 =  0.70;      % Persistence FF
        dgp.a1_2_2_sync_1 =  0.70;      % Persistence MF
        
        % A1 Matrix S=2
        dgp.a1_1_1_sync_2 =  0.25; % Persistence FF
        dgp.a1_2_2_sync_2 =  0.25;  % Persistence MF
        
        
    elseif modeldgp==6
        % Transition probabilities
        % Similar to dgp=5, with higher persistence of regime 1
        dgp.a12 = 5;
        dgp.b12 = 0.45;
        dgp.c12 = 0.45;
        
        dgp.a21 = 0.075;
        dgp.b21 = 0.05;
        dgp.c21 = 0.05;
        
        % Variances FF and MF
        dgp.s_1_1_sync_1 = 1.0;
        dgp.s_1_1_sync_2 = 2.5;
        dgp.s_2_2_sync_1 = 1.0;
        dgp.s_2_2_sync_2 = 2.5;
        
        % Variance of GDP
        dgp.s_3_3_sync_1 = 0.5;
        dgp.s_3_3_sync_2 = 0.75;
        
        
        % A0 Matrix S=1
        dgp.a0_2_1_sync_1 = 0.15;
        
        % A0 Matrix S=2
        dgp.a0_2_1_sync_2 = 0.35;
        
        % A1 Matrix S=1
        dgp.a1_1_1_sync_1 =  0.70;      % Persistence FF
        dgp.a1_2_2_sync_1 =  0.70;      % Persistence MF
        
        % A1 Matrix S=2
        dgp.a1_1_1_sync_2 =  0.25; % Persistence FF
        dgp.a1_2_2_sync_2 =  0.25;  % Persistence MF
        
        
        
    elseif modeldgp==7
        
        % Transition probabilities
        % Similar to gdp=4, with time variation in variances of MF and FF
        dgp.a12 = 3.8;
        dgp.b12 = 0.45;
        dgp.c12 = 0.45;
        
        dgp.a21 = 1.5;
        dgp.b21 = 0.15;
        dgp.c21 = 0.15;
        
        % Variances FF and MF and GDP
        dgp.s_1_1_sync_1 = 0.3344;
        dgp.s_2_2_sync_1 = 0.4657;
        dgp.s_3_3_sync_1 = 0.7525;
        
        
        
        dgp.s_2_2_sync_2 = 0.9390;
        dgp.s_1_1_sync_2 = 0.7403;
        dgp.s_3_3_sync_2 = 2.8033;
        
        
        
        %A0 Matrix (Lower triangular elements)
        dgp.a0_2_1_sync_1 = -0.1781;
        dgp.a0_3_1_sync_1 =  0.2901;       % response of GDP to FF
        dgp.a0_3_2_sync_1 = -0.1334;       % response of GDP to MF
        
        dgp.a0_2_1_sync_2 = -0.0434;
        dgp.a0_3_1_sync_2 = 0.2941;        % response of GDP to FF
        dgp.a0_3_2_sync_2 = -0.5167;       % response of GDP to MF
        
        
        
        % A1 Matrix S=1
        dgp.a1_1_1_sync_1 =  0.8713;      % Persistence FF
        dgp.a1_1_2_sync_1 =  0.0083;
        dgp.a1_1_3_sync_1 =  -0.0394;
        
        dgp.a1_2_1_sync_1 = -0.21;
        dgp.a1_2_2_sync_1 =  0.4289;      % Persistence MF
        dgp.a1_2_3_sync_1 =  -0.0394;       % response of MF to GDPG(-1)
        
        dgp.a1_3_1_sync_1 = 0.2220;
        dgp.a1_3_2_sync_1 = 0.0743;
        dgp.a1_3_3_sync_1 = 0;       %Persistence GDP
        
        % A1 Matrix S=2
        dgp.a1_1_1_sync_2 =  0.8713;      % Persistence FF
        dgp.a1_1_2_sync_2 =  0.0083;
        dgp.a1_1_3_sync_2 =  -0.0394;
        
        dgp.a1_2_1_sync_2 = -0.21;
        dgp.a1_2_2_sync_2 =  0.4289;      % Persistence MF
        dgp.a1_2_3_sync_2 =  -0.0394;       % response of MF to GDPG(-1)
        
        dgp.a1_3_1_sync_2 = 0.2220;
        dgp.a1_3_2_sync_2 = 0.0743;
        dgp.a1_3_3_sync_2 = 0;       %Persistence GDP
        
        % C Matrix
        dgp.c_1_1_sync_1 = -0.0651;
        dgp.c_2_1_sync_1 =  0.0219;
        dgp.c_3_1_sync_1 =  0.5101;
        
        dgp.c_1_1_sync_2 =  0.1069;
        dgp.c_2_1_sync_2 =  0.0497;
        dgp.c_3_1_sync_2 = -0.8603;
        
        
    else
        error('Wrong DGP selected');
    end
    
    params_in = dgp;

end


%% Print A0, A1 and SIG matrices

[dgpparams,bad] = scriptParams_FULL(params_in,1,2);

fprintf('\n**************************************\n');
fprintf('         DGP Parameters \n');
fprintf('**************************************\n');

fprintf('\n *** A0 Matrix: State = 1 *** \n')
disp(dgpparams.A0_sync_1)
fprintf('\n *** A0 Matrix: State = 2 *** \n')
disp(dgpparams.A0_sync_2)

fprintf('\n *** A1 Matrix: State = 1 *** \n')
disp(dgpparams.A1_sync_1)
fprintf('\n *** A1 Matrix: State = 2 *** \n')
disp(dgpparams.A1_sync_2)

fprintf('\n *** C Matrix: State = 1 *** \n')
disp(dgpparams.C_sync_1)
fprintf('\n *** C Matrix: State = 2 *** \n')
disp(dgpparams.C_sync_2)

fprintf('\n *** SIG Matrix: State = 1 *** \n')
disp(dgpparams.SIG_sync_1)
fprintf('\n *** SIG Matrix: State = 2 *** \n')
disp(dgpparams.SIG_sync_2)

fprintf('\n *** Transition probabilities: Good to Bad *** \n')
disp([params_in.a12 params_in.b12 params_in.c12]);

fprintf('\n *** Transition probabilities: Bad to Good *** \n')
disp([params_in.a21 params_in.b21 params_in.c21]);

fprintf('\n *** Unconditional Probability: Good to Bad *** \n')
disp(1/(1+exp(params_in.a12)));

fprintf('\n *** Unconditional Probability: Bad to Good *** \n')
disp(1/(1+exp(params_in.a21)));


%%
% Model spectag
spectag = fSpecLabel(modelspec);


% VAR configuration
nlags=1;
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
% Simulate DGP
%========================================
rng(123);
[Y_out_all, Y_1_all,Y_2_all, s_out_all,eta,p12_out_all,p21_out_all] = fSimulateDGPIteretedFull(dgpparams,opt);

date0 = datestr(367,dataformat);

opt.tperiods2 = opt.tperiods;

start_date   = datestr(datetime(date0,'InputFormat',inputformat)+calmonths(nlags),dataformat);
end_date     = datestr(datetime(date0,'InputFormat',inputformat)+calmonths(opt.tperiods2),dataformat);

start_plot   = datestr(datetime(date0,'InputFormat',inputformat)+calmonths(nlags),dataformat);
end_plot     = datestr(datetime(date0,'InputFormat',inputformat)+calmonths(opt.tperiods2-opt.hh),dataformat);
end_plot_full= datestr(datetime(date0,'InputFormat',inputformat)+calmonths(opt.tperiods2),dataformat);

%========================================
% CONFIGURE DATES AND PLOT OPTIONS
%========================================

% Vector of dates for the full sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

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

% Organize data
rawdata.Dates = ((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';
rawdata.FF = squeeze(Y_out_all(1,:,:))';
rawdata.MF = squeeze(Y_out_all(2,:,:))';
rawdata.GDPG = squeeze(Y_out_all(3,:,:))';
rawdata.S  =  squeeze(s_out_all)';
rawdata.p12 = squeeze(p12_out_all)';
rawdata.p21 = squeeze(p21_out_all)';
%% DATA USED FOR ESTIMATION
fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
subplot(311)
h=plot([rawdata.FF(SEED,:)' rawdata.MF(SEED,:)' rawdata.GDPG(SEED,:)']);
set(h,{'LineWidth'},{1;1;2})
set(h,{'LineStyle'},{'-';'--';'-'})
axis tight
xlim([1 550]);legend('FF','MF','GDPG','Location','Best'); legend boxoff
title(['Data'],'Interpreter','Latex','FontSize',16);

subplot(312)
plot([rawdata.S(SEED,:)' ],'LineWidth',2);
axis tight
xlim([1 opt.tperiods2]);legend('FF','MF','GDPG')
title('Regime: 1= good, 2=bad','Interpreter','Latex','FontSize',16); legend off;

subplot(313)
plot(rawdata.p12(SEED,:)','LineWidth',2); hold on
plot(rawdata.p21(SEED,:)','--','LineWidth',2);
axis tight
legend('p12','p21','Location','Best'); legend boxoff;
xlim([1 550]);ylim([0 1]);
title('Transition probabilities','Interpreter','Latex','FontSize',16);

set(fig1,'PaperOrientation','portrait');
set(fig1, 'PaperSize', figSize);
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'Units','inches');
set(fig1, 'PaperPositionMode', 'auto');
set(fig1, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig

%%
%==========================================================================
% LOAD DATA FOR ESTIMATION
%==========================================================================

%[db, db_full,startdb,enddb,tex] = fLoadData(datafilename,sheetuse,start_date,end_date,opt.hh);
[db, db_full,startdb,enddb, tex]= fLoadDataSim(rawdata,SEED,start_date,end_date,opt.hh);

% Collect MF, FF and trend (estimation sample)
FF    = db.FF.data(1:end);
MF    = db.MF.data(1:end);
TR    = db.TRENDH.data(1:end);
GDPG  = db.GDPG.data(1:end);
GDPGH = db.GDPGH.data(1:end);

% Collect MF, FF and trend (full sample)
FF_full = db_full.FF.data(1:end);
MF_full = db_full.MF.data(1:end);
TR_full = db_full.TRENDH.data(1:end);
GDPG_full  = db_full.GDPG.data(1:end);
GDPGH_full = db_full.GDPGH.data(1:end);


%%
figure(20);clf;
fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
subplot(211)
h=plot([rawdata.FF(SEED,:)' rawdata.MF(SEED,:)' rawdata.GDPG(SEED,:)']);
set(h,{'LineWidth'},{1;1;2})
set(h,{'LineStyle'},{'-';'--';'-'})
axis tight
xlim([1 1000]);ylim([-8 8]);legend('FF','MF','GDPG','Location','Best'); legend boxoff
title(['Data'],'Interpreter','Latex','FontSize',16);
subplot(212)
h=plot([FF_full MF_full GDPGH_full]);
set(h,{'LineWidth'},{1;1;2})
set(h,{'LineStyle'},{'-';'--';'-'})
axis tight
xlim([1 1000]);ylim([-8 8]);legend('FF','MF','GDPGH','Location','Best'); legend boxoff
title(['Data'],'Interpreter','Latex','FontSize',16);

%%
figure(30);clf
plot(Y_1_all(3,:,SEED)); hold on;
plot(Y_2_all(3,:,SEED),'--');
plot(Y_out_all(3,:,SEED),'k-');
legend('Regime 1','Regime 2','data');

%%


return
%xlim([1 1000]);legend('FF','MF','GDP')
%%

% Model name and title for posterior mode table
% Note: iterated model posterior mode depends on opt.hh through the estimation sample

datafilename = ['model' num2str(modelspec)];
sheetuse = ['SEED' num2str(SEED)];

modelname = ['DGP' num2str(modeldgp) '_EST' num2str(modelspec) '_' strrep(spectag,'_','')  '_' sheetuse '_H' num2str(opt.hh) tag];
%%
% Folders
paramsfolder = ['results/Params/' opt.paramsUse '/FullDGP/'];
figfolder   = ['results/Figures/' opt.paramsUse '/FullDGP/DGP' num2str(modeldgp) '/'];

% Create folders to store results
if exist(paramsfolder,'dir')==0;  mkdir(paramsfolder); end
if exist(figfolder,'dir')==0;  mkdir(figfolder); end
%%
% Diary file
dfile = ([figfolder sheetuse '_' spectag '_params.txt']);
if exist(dfile, 'file')
    delete(dfile);
end


% diary(dfile)
fprintf('\n Model: %s \n', modelname)

%==========================================================================
% RISE ESTIMATION
%==========================================================================

if strcmp(opt.paramsUse,'mode')
    
    % Load results or estimate posterior mode
    if exist([paramsfolder modelname '.mat'],'file')==0 || opt.force_estimate==1
        % Estimate posterior mode
        opt.dir_or_it = 2;
        varlist={'FF','MF','GDPG'}';
        
        run scriptEstimateMode_FULL_PCB.m;
        %%
        sv_iter = sv;
        params_iter   = posterior_mode(sv_iter);
        
        
        % Store solution structure
        %save([paramsfolder modelname '.mat'],'sPMode')
        clearvars varlist sv pmode
        
        
        %%        % Estimate posterior mode for the direct model
        opt.dir_or_it = 1;
        
        varlist={'FF','MF','GDPGH'}';
        
        
        run scriptEstimateMode_FULL_PCB.m;
        
        sv_dir = sv;
        params_dir   = posterior_mode(sv_dir);
        
        % Store solution structure
        %save([paramsfolder modelname '.mat'],'sPMode')
        
        clearvars varlist sv pmode
        
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
% Predictive Density at Posterior Mode or MCMC
%==========================================================================

opt.nParamDraws = 1; %length(pmode.a12);

%% Print A0, A1 and SIG matrices from ESTIMATION ITERATED MODEL
estimatedparams = scriptParams_FULL(params_iter,1,modelspec);

fprintf('\n**************************************\n');
fprintf('         Estimation Results \n');
fprintf('**************************************\n');

fprintf('\n *** A0 Matrix: State = 1 *** \n')
fprintf('\n     ESTIMATED                 |      DGP \n')
disp([estimatedparams.A0_sync_1 dgpparams.A0_sync_1])

fprintf('\n *** A0 Matrix: State = 2 *** \n')
fprintf('\n     ESTIMATED                 |      DGP \n')
disp([estimatedparams.A0_sync_2 dgpparams.A0_sync_2])

fprintf('\n *** A1 Matrix: State = 1 *** \n')
fprintf('\n     ESTIMATED                 |      DGP \n')
disp([estimatedparams.A1_sync_1 dgpparams.A1_sync_1])
fprintf('\n *** A1 Matrix: State = 2 *** \n')
fprintf('\n     ESTIMATED                 |      DGP \n')
disp([estimatedparams.A1_sync_2 dgpparams.A1_sync_2])

fprintf('\n *** C Matrix: State = 1 *** \n')
fprintf('\n  EST.     |    DGP \n')
disp([estimatedparams.C_sync_1 dgpparams.C_sync_1] )


fprintf('\n *** C Matrix: State = 2 *** \n')
fprintf('\n  EST.     |    DGP \n')
disp([estimatedparams.C_sync_2 dgpparams.C_sync_2] )

fprintf('\n *** SIG Matrix: State = 1 *** \n')
fprintf('\n     ESTIMATED                 |      DGP \n')
disp([estimatedparams.SIG_sync_1 dgpparams.SIG_sync_1])

fprintf('\n *** SIG Matrix: State = 2 *** \n')
fprintf('\n     ESTIMATED                 |      DGP \n')
disp([estimatedparams.SIG_sync_2 dgpparams.SIG_sync_2])

fprintf('\n *** Transition probabilities: Good to Bad *** \n')
fprintf('\n  EST.     |    DGP \n')
disp([estimatedparams.a12 estimatedparams.b12 estimatedparams.c12; params_in.a12 params_in.b12 params_in.c12]');

fprintf('\n *** Transition probabilities: Bad to Good *** \n')
fprintf('\n  EST.     |    DGP \n')
disp([estimatedparams.a21 estimatedparams.b21 estimatedparams.c21;params_in.a21 params_in.b21 params_in.c21]');

fprintf('\n *** Unconditional Probability: Good to Bad *** \n')
fprintf('\n  EST.     |    DGP \n')
disp([1/(1+exp(estimatedparams.a12)) 1/(1+exp(params_in.a12))]);

fprintf('\n *** Unconditional Probability: Bad to Good *** \n')
fprintf('\n  EST.     |    DGP \n')
disp([1/(1+exp(estimatedparams.a21)) 1/(1+exp(params_in.a21))]);


%% Plot Filtered States vs Simulated from DGP
[~,~,~,f]=filter(sv_iter);
% Smoothed probabilities
p_reg1_smoothed = f.smoothed_regime_probabilities.regime1.data;
p_reg2_smoothed = f.smoothed_regime_probabilities.regime2.data;
% Filtered probabilities
p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;


%%
% % Extract the index for the desired density plots
opt.date_index = []; %find(ismember(dates_full, target_dates));
opt.tperiods=500;
% Direct forecast
[quantilesD,pdensityD] = fPredDensityDirectFull(sv_dir,params_dir,FF_full,MF_full,GDPGH_full,TR_full,opt,modelspec);

% Iterated forecast
[quantilesI, pdensityI] = fPredDensityIteratedFull(sv_iter,params_iter,FF_full,MF_full,GDPG_full,TR_full,opt,'drawst',modelspec);

% Iterated forecast: DGP
[quantilesDGP, pdensityDGP] = fPredDensityIteratedFull(sv_iter,dgp,FF_full,MF_full,GDPG_full,TR_full,opt,'drawst',4);

%%

%==========================================================================
% Collect Probabilities
%==========================================================================
[~,~,~,f]=filter(sv_iter);
% Smoothed probabilities
p_reg1_smoothed = f.smoothed_regime_probabilities.regime1.data;
p_reg2_smoothed = f.smoothed_regime_probabilities.regime2.data;
% Filtered probabilities
p_reg1_filtered = f.filtered_regime_probabilities.regime1.data;
p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;


p_reg_bad_filtered = p_reg2_filtered;
p_reg_bad_smoothed = p_reg2_smoothed;


%%
% ==========================================================================
% Quantile Plots
% ==========================================================================
figureTitleTag = ['DGP: ' num2str(modeldgp) ', Estimation: ' num2str(modelspec) ', ' strrep(spectag,'_',' ')  ', h=' num2str(opt.hh) ];

% Add additional objects to structure
fig2=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
subplot(211)
hold on
l1=plot(quantilesD.dYsim_10,'--','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Iterated)');
l2=plot( quantilesD.dYsim_90,'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Iterated)');
l3=plot(quantilesDGP.dYsim_10,'-','Color',colors(1,:),'LineWidth', 2,'DisplayName','10th (DGP)');
l4=plot(quantilesDGP.dYsim_90,'-','Color',colors(60,:),'LineWidth', 2,'DisplayName','90th (DGP)');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title(['Quantiles Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
axis tight
xlim([1 opt.tperiods]); ylim([-8, 2]);

subplot(212)
hold on
l1=plot(quantilesI.dYsim_10,'--','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th (Direct)');
l2=plot( quantilesI.dYsim_90,'--','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th (Direct)');
l3=plot(quantilesDGP.dYsim_10,'-','Color',colors(1,:),'LineWidth', 2,'DisplayName','10th (DGP)');
l4=plot(quantilesDGP.dYsim_90,'-','Color',colors(60,:),'LineWidth', 2,'DisplayName','90th (DGP)');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title(['Quantiles Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
axis tight
xlim([1 opt.tperiods]);ylim([-8, 2]);
set(fig2,'PaperOrientation','portrait');
set(fig2, 'PaperSize', figSize);
set(fig2, 'PaperUnits', 'inches');
set(fig2, 'Units','inches');
set(fig2, 'PaperPositionMode', 'auto');
set(fig2, 'Position', [figSize(1)/5 figSize(2)/10 figSize(1) figSize(2)]);
tightfig
%% Regime specific fitted values
[~,Fit]=residuals(sv_iter);
fit1I = Fit.GDPG.data(:,1);
fit2I = Fit.GDPG.data(:,2);

[~,Fit]=residuals(sv_dir);
fit1D = Fit.GDPGH.data(:,1);
fit2D = Fit.GDPGH.data(:,2);



fig3=figure;
subplot(211)
hold on
plot( fit1I,'Color',colors(1,:),'LineWidth', 1.5)
plot( fit2I,'--','Color',colors(15,:),'LineWidth', 1.5)
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
leg = {'Regime 1','Regime 2'};
title(['Regime Fitted Values, Iterated: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
axis tight
ax=gca;
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

subplot(212)
hold on
plot( fit1D,'Color',colors(1,:),'LineWidth', 1.5)
plot( fit2D,'--','Color',colors(15,:),'LineWidth', 1.5)
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
leg = {'Regime 1','Regime 2'};
title(['Regime Fitted Values, Direct: ' figureTitleTag ],'Interpreter','Latex','FontSize',16);
axis tight
ax=gca;
%ax.XTick = datenum(dates(sd:numticks:ed));
%datetick('x','yyyy','keepticks')
%set(gca, 'XLim', [dates(sd), dates(ed)])
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
%%

fig4=figure;clf; hold on
l3=plot(1-quantilesDGP.st_t_mean,'-','Color',colors(1,:),'LineWidth', 2.5,'DisplayName','Filtered from DGP');
l4=plot(1-quantilesI.st_t_mean,'--','Color',colors(11,:),'LineWidth', 2.5,'DisplayName','Filtered from Iterated');
l5=plot(1-quantilesD.st_t_mean,'-.','Color',colors(21,:),'LineWidth', 2.5,'DisplayName','Simulated from Direct');
l6 = plot(rawdata.S(SEED,:)-1,':','Color',[0 0 0],'LineWidth', 1.5,'DisplayName','DGP');
legend([l3 l4 l5 l6],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title(['Probability of Bad Regime: ' figureTitleTag],'FontSize',16','Interpreter','Latex');
axis tight
ylim([0 1])
xlim([1 opt.tperiods]);
% Set figure appereance
set(fig4,'PaperOrientation','portrait');
set(fig4, 'PaperSize', figSize);
set(fig4, 'PaperUnits', 'inches');
set(fig4, 'Units','inches');
set(fig4, 'PaperPositionMode', 'auto');
set(fig4, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
tightfig;

%%
% if opt.saveit
print('-dpdf',fig1,[figfolder 'Data_' modelname ],'-bestfit');
print('-dpdf',fig2,[figfolder 'Quantiles_' modelname ],'-bestfit');
print('-dpdf',fig3,[figfolder 'Fitted_' modelname ],'-bestfit');
print('-dpdf',fig4,[figfolder 'ProbReg_' modelname ],'-bestfit');
% end
%
% diary off

%%
rise_exit
