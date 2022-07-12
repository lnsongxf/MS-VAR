%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous



%% NOTE:
% Different from the previous simulation, now MF reacts contemporaneously
% to FF

%% housekeeping
clear; close all; clc;

rng(777); % Fixed seed for replication purposes

% Options
saveit    = 1; % 1 = save graphs
save_data = 0; % 1 = save simulated data
const     = 1; % 1= have a constant in transition probability
normal    = 0; % 1 = use normal distribution, 0 = gamma distribution

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

sheetuse = 'Sim';

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

% Data vintage, sample and country selection
start_date   = '2100-Jan';
end_date     = '2199-Dec';

% Define target periods for density cuts
target_dates = [...%datenum('2020-Jan',dataformat),...
    %%datenum('2020-Feb',dataformat),...
    ];

% Define dates for plots
% Data vintage, sample and country selection
start_plot     = '2100-Feb';
end_plot       = '2198-Dec';
end_plot_full  = '2199-Dec';

% VAR configuration
nlags=1;
exog={};
constant=true;
panel=[];


% Vector of dates for the full sample
Dates = datenum((datetime(start_date,'InputFormat',inputformat)):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';


% Figure options
FontSize = 16;
numticks = 48;
figSize = [12 6];
linestyle = {'-','--',':'};
colors = cbrewer('div', 'RdYlBu', 64);
colors2 = cbrewer('qual', 'Set1', 8);
left_color = [0 0 0];
right_color = colors(2,:);

% FF = c_1_1             + a1_1_1*FF(-1) + a1_1_2*MF(-1)
% FF s.d. = s_1_1
% MF = c_2_1 + a0_2_1*FF + a1_2_1*FF(-1) + a1_2_2*MF(-1)
% MF s.d. = s_2_2
% GDP, normal = c_3_1_sync_1 + a0_3_1_sync_1*FF + a0_3_2_sync_1*MF
% GDP, normal s.d. = s_3_3_sync_1
% GDP, bad = c_3_1_sync_2 + a0_3_1_sync_2*FF + a0_3_2_sync_2*MF
% GDP, bad s.d. = s_3_3_sync_2
% 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))'
% 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))'

dgp.c_1_1 = 0; dgp.a0_1_1 = 1;     dgp.a0_1_2 = 0; dgp.a0_1_3 = 0; dgp.a1_1_1 =  0.75; dgp.a1_1_2 = -0.05; dgp.a1_1_3=0; dgp.s_1_1 = 1;
dgp.c_2_1 = 0; dgp.a0_2_1 = -0.15; dgp.a0_2_2 = 1; dgp.a0_2_3 = 0; dgp.a1_2_1 = -0.1; dgp.a1_2_2 = 0.85;  dgp.a1_2_3=0; dgp.s_2_2 = 1;
dgp.c_3_1_sync_1 =  0.5; dgp.a0_3_1_sync_1 = -0.1; dgp.a0_3_2_sync_1 = 0.1; dgp.a0_3_3_sync_1 = 0; dgp.s_3_3_sync_1 = 0.5;
dgp.c_3_1_sync_2 = -0.5; dgp.a0_3_1_sync_2 = -0.3; dgp.a0_3_2_sync_2 = 0.3; dgp.a0_3_3_sync_2 = 0; dgp.s_3_3_sync_2 = 2;
dgp.a12 = 2; dgp.a21 = 2; dgp.b12 = 0.5; dgp.b21 = 0.5; dgp.c12 = 0.5; dgp.c21 = 0.5;


%% SIM


% 2 = BAD, 1 = NORMAL
burn = 500; % 500 burn, 500 keep
tt=1200+burn; % time periods
rr = 600; % draws
f_draw = zeros(rr,tt);m_draw = zeros(rr,tt);
p12_draw = NaN(tt-burn,1); p21_draw = NaN(tt-burn,1);
y_mat = NaN(rr,tt-burn);

st_mat = NaN(rr,tt); st_mat(:,1) = ones(rr,1);

for dd=1:rr
    for jj=2:tt
        eta1 = randn(2,1); % financial and macro shocks
        f_draw(dd,jj) = dgp.c_1_1 +                          + dgp.a1_1_1*f_draw(dd,jj-1) + dgp.a1_1_2*m_draw(dd,jj-1) + dgp.s_1_1*eta1(1,1);
        m_draw(dd,jj) = dgp.c_2_1 + dgp.a0_2_1*f_draw(dd,jj) + dgp.a1_2_1*f_draw(dd,jj-1) + dgp.a1_2_2*m_draw(dd,jj-1) + dgp.s_2_2*eta1(2,1);
        
        p12 = 1/(1+exp(dgp.a12-dgp.b12*(f_draw(dd,jj))+dgp.c12*(m_draw(dd,jj)))); % 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))' %USE THIS!
        p21 = 1/(1+exp(dgp.a21+dgp.b21*(f_draw(dd,jj))-dgp.c21*(m_draw(dd,jj)))); % 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))' %USE THIS!
        
        p11 = 1 - p12; % probability of remaining in normal
        p22 = 1 - p21; % probability of remaining in bad
        
        
        udraw = rand(1);
        
        if st_mat(dd,jj-1) == 1 % started in good regime
            if udraw > p11
                st_temp = 2; % switch from good to bad
            else
                st_temp = 1; % don't switch and remain in good
            end
            
        else % start in bad regime
            if udraw > p22
                st_temp = 1; % switch from bad to good
            else
                st_temp = 2; % don't switch and remain in bad
            end
            
        end
        
        eta2 = randn(1,1); % GDP shock
        
        if st_temp==2
            y_draw = dgp.c_3_1_sync_2 + dgp.a0_3_1_sync_2*f_draw(dd,jj) + dgp.a0_3_2_sync_2*m_draw(dd,jj) + dgp.s_3_3_sync_2*eta2;
        else
            y_draw = dgp.c_3_1_sync_1 + dgp.a0_3_1_sync_1*f_draw(dd,jj) + dgp.a0_3_2_sync_1*m_draw(dd,jj) + dgp.s_3_3_sync_1*eta2;
        end
        
        %         if dd==rr % randomly assign a state for time jj
        %             ind_st = round(rand(1)*(length(st_temp))+0.5,0);
        %             st_draw(jj) = st_temp(ind_st);
        %         end
        st_mat(dd,jj) = st_temp;
        if jj>burn
            y_mat(dd,jj-burn) = y_draw;
            %             f_mat(dd,jj-burn) = f_draw(dd,jj); m_mat(dd,jj-burn) = m_draw(dd,jj);
            p12_draw(dd,jj-burn) = p12; p21_draw(dd,jj-burn) = p21;
        end
    end
end

st_mat = st_mat(:,burn+1:end);
f_mat =  f_draw(:,burn+1:end);
m_mat =  m_draw(:,burn+1:end);
% Make y_mat of GDP(t+1):
% y_mat = y_mat(:,2:end);
% f_mat = f_mat(:,1:end-1);
% m_mat = m_mat(:,1:end-1);
% st_mat = st_mat(:,1:end-1);
% p12_draw = p12_draw(1:end-1,:);
% p21_draw = p21_draw(1:end-1,:);


%% Quantiles

% 12-months ahead
y_mat_fut = NaN(size(y_mat,1),size(y_mat,2)-12);
for hh=1:size(y_mat,1)
    for ww=1:size(y_mat,2)-12
        y_mat_fut(hh,ww) = mean(y_mat(hh,ww+1:ww+12));
    end
end

% Compute percentiles
dY_25 = prctile(y_mat,25)'; dY_25_fut = prctile(y_mat_fut,25)';
dY_75 = prctile(y_mat,75)'; dY_75_fut = prctile(y_mat_fut,75)';
dY_10 = prctile(y_mat,10)'; dY_10_fut = prctile(y_mat_fut,10)';
dY_90 = prctile(y_mat,90)'; dY_90_fut = prctile(y_mat_fut,90)';

ds_50 = prctile(st_mat,50)';
ds_10 = prctile(st_mat,10)';
ds_90 = prctile(st_mat,90)';

figure; clf;
hold on
l1=plot(1:tt-burn, dY_10,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(1:tt-burn, dY_25,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(1:tt-burn, dY_75,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(1:tt-burn, dY_90,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Simulated GDP distribution in t','FontSize',16','Interpreter','Latex');
axis tight

figure; clf;
hold on
l1=plot(1:tt-burn-12, dY_10_fut,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(1:tt-burn-12, dY_25_fut,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(1:tt-burn-12, dY_75_fut,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(1:tt-burn-12, dY_90_fut,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');

legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Simulated GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
axis tight

% States
figure; clf;
hold on
plot(1:tt-burn, mean(st_mat),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','States');
set(gca,'children',flipud(get(gca,'children')))
axis tight
ylim([1 2]);
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Average states','FontSize',16','Interpreter','Latex');


figure; clf;
hold on
l1=plot(1:tt-burn, mean(p12_draw),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Good to bad');
l2=plot(1:tt-burn, mean(p21_draw),'Color',colors(15,:),'LineWidth', 1.5,'DisplayName','Bad to good');
legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
axis tight
ylim([0 1]);
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Switching probabilities','FontSize',16','Interpreter','Latex');

%% SAVE

Res_iterated.GDPG = y_mat; Res_iterated.FF = f_mat; Res_iterated.MF = m_mat;
Res_iterated.dY_25 = dY_25; Res_iterated.dY_75 = dY_75;
Res_iterated.dY_10 = dY_10; Res_iterated.dY_90 = dY_90;
Res_iterated.st_mat = st_mat;
Res_iterated.Dates = datetime(Dates,'ConvertFrom','datenum');

Res_direct.y_mat = y_mat_fut; Res_direct.FF = f_mat; Res_direct.MF = m_mat;
Res_direct.dY_25 = dY_25_fut; Res_direct.dY_75 = dY_75_fut;
Res_direct.dY_10 = dY_10_fut; Res_direct.dY_90 = dY_90_fut;


% FF = c_1_1 + a0_1_1*FF + a0_1_2*MF + a0_1_3*GDP + a1_1_1*FF(-1) + a1_1_2*MF(-1) + a1_1_3*GDP(-1)
% FF s.d. = s_1_1
% MF = c_2_1 + a0_2_1*FF + a0_2_2*MF + a0_2_3*GDP + a1_2_1*FF(-1) + a1_2_2*MF(-1) + a1_2_3*GDP(-1)
% MF s.d. = s_2_2
% GDP, normal = c_3_1_sync_1 + a0_3_1_sync_1*FF + a0_3_2_sync_1*MF + a0_3_3_sync_1*GDP
% GDP, normal s.d. = s_3_3_sync_1
% GDP, bad = c_3_1_sync_2 + a0_3_1_sync_2*FF + a0_3_2_sync_2*MF + a0_3_3_sync_2*GDP
% GDP, bad s.d. = s_3_3_sync_2
% 'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))'
% 'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))'

% dgp.c_1_1 = A_FF(2); dgp.a0_1_1 = 0; dgp.a0_1_2 = 0; dgp.a0_1_3 = 0; dgp.a1_1_1 = A_FF(1); dgp.a1_1_2 = 0; dgp.a1_1_3=0; dgp.s_1_1 = STD_FF;
% dgp.c_2_1 = A_MF(2); dgp.a0_2_1 = 0; dgp.a0_2_2 = 0; dgp.a0_2_3 = 0; dgp.a1_2_1 = 0; dgp.a1_2_2 = A_MF(1); dgp.a1_2_3=0; dgp.s_2_2 = STD_MF;
% dgp.c_3_1_sync_1 = mu_1; dgp.a0_3_1_sync_1 = beta_1; dgp.a0_3_2_sync_1 = gamma_1; dgp.a0_3_3_sync_1 = 0; dgp.s_3_3_sync_1 = sig_1;
% dgp.c_3_1_sync_2 = mu_2; dgp.a0_3_1_sync_2 = beta_2; dgp.a0_3_2_sync_2 = gamma_2; dgp.a0_3_3_sync_2 = 0; dgp.s_3_3_sync_2 = sig_2;
% dgp.a12 = a12; dgp.a21 = a21; dgp.b12 = b12; dgp.b21 = b21; dgp.c12 = c12; dgp.c21 = c21;
if save_data ==1
    save Res_sim Res_iterated Res_direct dgp
end

%% Take a random seed (555)

SEED = 555;



for imodel=mymodels
    
    modelspec = imodel;
    
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
    % MODEL SPECIFIC VARIABLES
    %========================================
    
    if opt.dir_or_it ==1
        fcstype= 'direct';
        varlist={'FF','MF','GDPGH'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Simul/Direct/'];
        %texfolder   = ['results/Params/' opt.paramsUse '/Direct/Tex/'];
        figfolder   = ['results/Figures/' opt.paramsUse '/Simul/Direct/' spectag '/'];
        
    elseif opt.dir_or_it ==2
        fcstype= 'iterated';
        varlist={'FF','MF','GDPG'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Simul/Iterated/'];
        %texfolder   = ['results/Params/' opt.paramsUse '/Iterated/Tex/'];
        figfolder   =  ['results/Figures/' opt.paramsUse '/Simul/Iterated/' spectag '/'];
    end
    
    %========================================
    % CONFIGURE DATES AND PLOT OPTIONS
    %========================================
    
    
    % Index of dates
    % Vector of dates for the full sample
    dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';
    
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
    dfile = ([figfolder 'Sim_' spectag '_params.txt']);
    if exist(dfile, 'file')
        delete(dfile);
    end
    
    %==========================================================================
    % LOAD DATA
    %==========================================================================
    
    
    [db, db_full,startdb,enddb,tex] = fLoadDataSim(Res_iterated,SEED,start_date,end_date,opt.hh);
    
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
    
    modelname = ['Sim_' startdb enddb ...
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
    [outparams,~] = scriptParams_FULL(params_in,1,modelspec);
    
    disp(outparams.C_sync_1)
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
        [quantiles,pdensity, badind] = fPredDensityDirectFull(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt,modelspec);
        
    elseif opt.dir_or_it==2
        % Iterated forecast
        [quantiles, pdensity, badind] = fPredDensityIteratedFull(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,'forecastst',modelspec);
    end
    
    % Check bad regime indicator at mode
    [~,bad] = scriptParams_FULL(pmode,1,modelspec);
    if bad==3
        fprintf('\n Model with inverse mean/variance correlation \n')
        return
    end
    
    
    %     % Predictive density at specified episodes
    %     for ii=1:numel(opt.date_index)
    %         %[pdfi(ii,:),xi(ii,:)]  = ksdensity(pdensity(:,ii));
    %         [pdfi(ii,:),xi(ii,:)]  = ksdensity(pdensity.realized(:,ii));
    %         [pdfi_good(ii,:),xi_good(ii,:)]  = ksdensity(pdensity.good(:,ii));
    %         [pdfi_bad(ii,:),xi_bad(ii,:)]  = ksdensity(pdensity.bad(:,ii));
    %     end
    
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
    %     linetype= {'--','-','-.',':','-.'};
    %     for ii=1:numel(opt.date_index)
    %
    %         fig2=figure;
    %         hold on
    %         plot(xi(ii,:),pdfi(ii,:),linetype{2},'Color',[0 0 0 ],'LineWidth', 3,'DisplayName',[datestr(dates_full(opt.date_index(ii)),dataformat)]); hold on;
    %         plot(xi_good(ii,:),pdfi_good(ii,:),linetype{4},'Color',colors(60,:),'LineWidth', 3,'DisplayName','Good regime'); hold on;
    %         plot(xi_bad(ii,:),pdfi_bad(ii,:),linetype{5},'Color',colors(1,:),'LineWidth', 3,'DisplayName','Bad regime'); hold on;
    %         if opt.dir_or_it==1
    %             title(['Direct: ' figureTitleTag ],'Interpreter','Latex');
    %         else
    %             title(['Iterated: ' figureTitleTag ],'Interpreter','Latex');
    %         end
    %         axis tight
    %         ylimits = ylim; xlimits = xlim;
    % %         text(xlimits(1)*0.865,0.8*ylimits(2),['Prob bad regime = ' num2str((1-quantiles.st_t_mean(opt.date_index(ii)))*100) '\%'],'FontSize',13,'Interpreter','Latex');
    %         plot(zeros(1,length(xi_bad(1,:))),zeros(1,length(pdfi_bad(1,:))),'*','DisplayName',['Prob bad regime = ' num2str(round((1-quantiles.st_t_mean(opt.date_index(ii)))*100,1)) '\%']);
    %         vline(0,'k--');
    %         text(GDPGH_full(opt.date_index(ii)),0.95*ylimits(2),'realized value','color','red','FontSize',9);
    %         legend('Orientation','Vertical','Location','NW','interpreter','Latex');legend boxoff;
    %         ylabel('PDF','fontsize',10,'interpreter','Latex')
    %         xlabel([num2str(opt.hh) '-Months-Ahead of Real Activity Measure'],'fontsize',10,'Interpreter','Latex')
    %         set(gca, 'FontName', 'Times New Roman');
    %         set(gca, 'FontSize', FontSize);
    %         set(gca,'Layer','top')
    %         %set(gca,'XTick',unique(sort([0 round(GDPGH_full(opt.date_index(ii)),2) xlimits(1):2:xlimits(2)])))
    %         set(gca,'TickLabelInterpreter','Latex')
    %         axis tight
    %         %xlim([-18 8])
    %         set(fig2,'PaperOrientation','portrait');
    %         set(fig2, 'PaperSize', figSize);
    %         set(fig2, 'PaperUnits', 'inches');
    %         set(fig2, 'Units','inches');
    %         set(fig2, 'PaperPositionMode', 'auto');
    %         set(fig2, 'Position', [figSize(1)/5 figSize(2)/5 figSize(1) figSize(2)]);
    %         tightfig;
    %
    %         if opt.saveit
    %             print('-dpdf',fig2,[figfolder 'Densities' datestr(dates_full(opt.date_index(ii)),dataformat) '_' modelname ],'-bestfit');
    %
    %         end
    %
    %     end
    
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
    plot(dates(sd:ed), fitgood(sd:ed),'Color',colors(1,:),'LineWidth', 1.5)
    plot(dates(sd:ed), fitbad(sd:ed),'--','Color',colors(15,:),'LineWidth', 1.5)
    hold off
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
