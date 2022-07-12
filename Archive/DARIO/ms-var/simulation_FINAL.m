%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous



%% NOTE:
% Different from the previous simulation, now MF reacts contemporaneously
% to FF

%% housekeeping
clear;
close all; clc;

rng(777); % Fixed seed for replication purposes

% % Options
% saveit    = 1; % 1 = save graphs
% const     = 1; % 1= have a constant in transition probability
% normal    = 0; % 1 = use normal distribution, 0 = gamma distribution

% Important paths
addpath('../RISE_toolbox');         % RISE Toolbox
addpath(genpath('functions'));      % Main functions
addpath(genpath('scripts'));        % Main scripts
addpath(genpath('textools'));       % Tools to produce tables
addpath(genpath('auxtools'));       % Tools for plotting and other
addpath(genpath('cbrewer'));        % Color mixing tools

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
opt.force_estimate = 0; % 1 = force to reestimate model
opt.save_data      = 0; % 1 = save simulated data
sim_alg = 2;            %[1] Algorithm 1, [2] Algorithm 2
quantiles_fig      = 1; % 1 = make quantile comparison figures
thrash =0;

% Posterior mode options
opt.nDraws      = 5000;  % Number of simulations for each time period

% MCMC simulation options
opt.nParamDraws  = 1000;  % Draws from posterior. has to be less than options.N;

% Figure options
opt.saveit       = 1;     % 1 = save figures

% Select which models to run
mymodels = [8];
% if opt.dir_or_it ==1
%     mymodels = 4;
% elseif opt.dir_or_it ==2 % Correct specification is model 4
%     mymodels = 4;
% end

if sim_alg==1
    tag = ['Alg1'];
else
    tag = ['Alg2'];
end

if thrash ==1
    if sim_alg==1
        tag = ['Alg1_thrash'];
    else
        tag = ['Alg2_thrash'];
    end
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

SEED = 555;

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


%% Define DGP

% FF, normal       = c_1_1_sync_1 + a0_1_2_sync_1*MF + a1_1_1_sync_1*FF(-1) + a1_1_2_sync_1*MF(-1) + a1_1_3_sync_1*GDP(-1)
% FF, bad          = c_1_1_sync_2 + a0_1_2_sync_2*FF + a1_1_1_sync_2*FF(-1) + a1_1_2_sync_2*MF(-1) + a1_1_3_sync_2*GDP(-1)
dgp.c_1_1_sync_1 = 0; dgp.a0_1_2_sync_1 = 0; dgp.a1_1_1_sync_1 =  0.75; dgp.a1_1_2_sync_1 = -0.05; dgp.a1_1_3_sync_1 = 0;
dgp.c_1_1_sync_2 = 0; dgp.a0_1_2_sync_2 = 0; dgp.a1_1_1_sync_2 =  0.75; dgp.a1_1_2_sync_2 = -0.05; dgp.a1_1_3_sync_2 = 0;
% FF, normal s.d.  = s_1_1_sync_1
% FF, bad s.d.     = s_1_1_sync_2
dgp.s_1_1_sync_1 = 1;
dgp.s_1_1_sync_2 = 1;

% MF, normal       = c_2_1_sync_1 + a0_2_1_sync_1*FF + a1_2_1_sync_1*FF(-1) + a1_2_2_sync_1*MF(-1) + a1_2_3_sync_1*GDP(-1)
% MF, bad          = c_2_1_sync_2 + a0_2_1_sync_2*FF + a1_2_1_sync_2*FF(-1) + a1_2_2_sync_2*MF(-1) + a1_2_3_sync_2*GDP(-1)
dgp.c_2_1_sync_1 = 0; dgp.a0_2_1_sync_1 = -0.15; dgp.a1_2_1_sync_1 = -0.10;  dgp.a1_2_2_sync_1 = 0.85;  dgp.a1_2_3_sync_1 = 0;
dgp.c_2_1_sync_2 = 0; dgp.a0_2_1_sync_2 = -0.15; dgp.a1_2_1_sync_2 = -0.10;  dgp.a1_2_2_sync_2 = 0.85;  dgp.a1_2_3_sync_2 = 0;
% MF, normal s.d.  = s_2_2_sync_1
% MF, bad s.d.    = s_2_2_sync_2
dgp.s_2_2_sync_1 = 1;
dgp.s_2_2_sync_2 = 1;

% GDP, normal      = c_3_1_sync_1 + a0_3_1_sync_1*FF + a0_3_2_sync_1*MF + a1_3_1_sync_1*FF(-1) + a1_3_2_sync_1*MF(-1) + dgp.a1_3_3_sync_1*GDP(-1)
% GDP, bad         = c_3_1_sync_2 + a0_3_1_sync_2*FF + a0_3_2_sync_2*MF + a1_3_1_sync_2*FF(-1) + a1_3_2_sync_2*MF(-1) + dgp.a1_3_3_sync_2*GDP(-1)
dgp.c_3_1_sync_1 =  1.5; dgp.a0_3_1_sync_1 = -0.1; dgp.a0_3_2_sync_1 = 0.1; dgp.a1_3_1_sync_1 = 0; dgp.a1_3_2_sync_1 = 0; dgp.a1_3_3_sync_1 = 0;
dgp.c_3_1_sync_2 = -1.5; dgp.a0_3_1_sync_2 = -0.3; dgp.a0_3_2_sync_2 = 0.3; dgp.a1_3_1_sync_2 = 0; dgp.a1_3_2_sync_2 = 0; dgp.a1_3_3_sync_2 = 0;
% dgp.c_3_1_sync_1 =  0.25; dgp.a0_3_1_sync_1 = -0.1; dgp.a0_3_2_sync_1 = 0.1; dgp.a1_3_1_sync_1 = 0; dgp.a1_3_2_sync_1 = 0; dgp.a1_3_3_sync_1 = 0;
% dgp.c_3_1_sync_2 = -0.25; dgp.a0_3_1_sync_2 = -0.3; dgp.a0_3_2_sync_2 = 0.3; dgp.a1_3_1_sync_2 = 0; dgp.a1_3_2_sync_2 = 0; dgp.a1_3_3_sync_2 = 0;
% GDP, normal s.d. = s_3_3_sync_1
% GDP, bad s.d.    = s_3_3_sync_2
dgp.s_3_3_sync_1 = 0.5;
dgp.s_3_3_sync_2 = 2;

% sync_tp_1_2      = 1/(1+exp(a12-b12*(FF)+c12*(MF)))
% sync_tp_2_1      = 1/(1+exp(a21+b21*(FF)-c21*(MF)))
dgp.a12 = 2;   dgp.a21 = 2;
dgp.b12 = 0.5; dgp.b21 = 0.5;
dgp.c12 = 0.5; dgp.c21 = 0.5;

% EXTRAS:
dgp.a0_1_3_sync_1 = 0; dgp.a0_1_3_sync_2 = 0;   % GDP always zero contemporaneous to FF
dgp.a0_2_3_sync_1 = 0; dgp.a0_2_3_sync_2 = 0;   % GDP always zero contemporaneous to MF
dgp.a0_1_1_sync_1 = 1; dgp.a0_1_1_sync_2 = 1;   % FF contemporaneous coefficients
dgp.a0_2_2_sync_1 = 1; dgp.a0_2_2_sync_2 = 1;   % MF contemporaneous coefficients
dgp.a0_3_3_sync_1 = 1; dgp.a0_3_3_sync_2 = 1;   % GDP contemporaneous coefficients


% Values that work:
% dgp.c_1_1 = 0; dgp.a0_1_1 =  1;    dgp.a0_1_2 = 0; dgp.a0_1_3 = 0; dgp.a1_1_1 =  0.75; dgp.a1_1_2 = -0.05; dgp.a1_1_3= 0; dgp.s_1_1 = 1;
% dgp.c_2_1 = 0; dgp.a0_2_1 = -0.15; dgp.a0_2_2 = 1; dgp.a0_2_3 = 0; dgp.a1_2_1 = -0.1;  dgp.a1_2_2 = 0.85;  dgp.a1_2_3= 0; dgp.s_2_2 = 1;
% dgp.c_3_1_sync_1 =  1.5; dgp.a0_3_1_sync_1 = -0.1; dgp.a0_3_2_sync_1 = 0.1; dgp.a0_3_3_sync_1 = 1; dgp.a1_3_1_sync_1 = 0; dgp.a1_3_2_sync_1 = 0; dgp.a1_3_3_sync_1 = 0; dgp.s_3_3_sync_1 = 0.5;
% dgp.c_3_1_sync_2 = -1.5; dgp.a0_3_1_sync_2 = -0.3; dgp.a0_3_2_sync_2 = 0.3; dgp.a0_3_3_sync_2 = 1; dgp.a1_3_1_sync_2 = 0; dgp.a1_3_2_sync_2 = 0; dgp.a1_3_3_sync_2 = 0; dgp.s_3_3_sync_2 = 2;
% dgp.a12 = 2; dgp.a21 = 2; dgp.b12 = 0.5; dgp.b21 = 0.5; dgp.c12 = 0.5; dgp.c21 = 0.5;

%% CReate regime comparison tables
comparison_table_1 = [dgp.c_1_1_sync_1,dgp.c_1_1_sync_1,dgp.c_3_1_sync_1;dgp.a0_1_1_sync_1,-dgp.a0_2_1_sync_1,-dgp.a0_3_1_sync_1;...
    dgp.a0_1_2_sync_1,dgp.a0_2_2_sync_1,-dgp.a0_3_2_sync_1;dgp.a0_1_3_sync_1,dgp.a0_2_3_sync_1,dgp.a0_3_3_sync_1;...
    dgp.a1_1_1_sync_1,dgp.a1_2_1_sync_1,dgp.a1_3_1_sync_1;dgp.a1_1_2_sync_1,dgp.a1_2_2_sync_1,dgp.a1_3_2_sync_1;...
    dgp.a1_1_3_sync_1,dgp.a1_2_3_sync_1,dgp.a1_3_3_sync_1;dgp.s_1_1_sync_1,dgp.s_2_2_sync_1,dgp.s_3_3_sync_1];
comparison_table_2 = [dgp.c_1_1_sync_2,dgp.c_1_1_sync_2,dgp.c_3_1_sync_2;dgp.a0_1_1_sync_2,-dgp.a0_2_1_sync_2,-dgp.a0_3_1_sync_2;...
    dgp.a0_1_2_sync_2,dgp.a0_2_2_sync_2,-dgp.a0_3_2_sync_2;dgp.a0_1_3_sync_2,dgp.a0_2_3_sync_2,dgp.a0_3_3_sync_2;...
    dgp.a1_1_1_sync_2,dgp.a1_2_1_sync_2,dgp.a1_3_1_sync_2;dgp.a1_1_2_sync_2,dgp.a1_2_2_sync_2,dgp.a1_3_2_sync_2;...
    dgp.a1_1_3_sync_2,dgp.a1_2_3_sync_2,dgp.a1_3_3_sync_2;dgp.s_1_1_sync_2,dgp.s_2_2_sync_2,dgp.s_3_3_sync_2];

%% SIM

% 2 = BAD, 1 = NORMAL
burn = 500; % 500 burn, 500 keep
opt.tt=1200+burn; % time periods
opt.rr = 600; % draws

if sim_alg ==2
    [f_draw,m_draw,y_draw,st_mat,p12_draw,p21_draw,eta] =  fSimulateData(dgp,opt);
elseif sim_alg ==1
    [f_draw,m_draw,y_draw,st_mat,p12_draw,p21_draw,eta] =  fSimulateData_test(dgp,opt);
end

p12_draw = p12_draw(:,burn+1:end); p21_draw = p21_draw(:,burn+1:end); st_mat = st_mat(:,burn+1:end);
f_mat  = f_draw(:,burn+1:end); m_mat  = m_draw(:,burn+1:end); y_mat  = y_draw(:,burn+1:end);
eta = eta(:,:,burn+1:end); eta_seed = squeeze(eta(:,SEED,:));

% f_mat  = f_draw; m_mat  = m_draw; y_mat  = y_draw;

% [f_draw_test,m_draw_test,y_draw_test,st_mat_test,p12_draw_test,p21_draw_test] =  fSimulateData_test(dgp,opt);

% p12_draw_test = p12_draw_test(:,burn+1:end); p21_draw_test = p21_draw_test(:,burn+1:end); st_mat_test = st_mat_test(:,burn+1:end);
% f_mat_test  = f_draw_test(:,burn+1:end); m_mat_test  = m_draw_test(:,burn+1:end); y_mat_test  = y_draw_test(:,burn+1:end);

% f_mat_test  = f_draw_test; m_mat_test  = m_draw_test; y_mat_test  = y_draw_test;


%% Quantiles

% 12-months ahead
y_mat_fut = NaN(size(y_mat,1),size(y_mat,2)-12);
for hh=1:size(y_mat,1)
    for ww=1:size(y_mat,2)-12
        y_mat_fut(hh,ww) = mean(y_mat(hh,ww+1:ww+12));
    end
end

% y_mat_fut_test = NaN(size(y_mat_test,1),size(y_mat_test,2)-12);
% for hh=1:size(y_mat_test,1)
%     for ww=1:size(y_mat_test,2)-12
%         y_mat_fut_test(hh,ww) = mean(y_mat_test(hh,ww+1:ww+12));
%     end
% end

% Compute percentiles
dY_25 = prctile(y_mat,25)'; dY_25_fut = prctile(y_mat_fut,25)';
dY_75 = prctile(y_mat,75)'; dY_75_fut = prctile(y_mat_fut,75)';
dY_10 = prctile(y_mat,10)'; dY_10_fut = prctile(y_mat_fut,10)';
dY_90 = prctile(y_mat,90)'; dY_90_fut = prctile(y_mat_fut,90)';

ds_50 = prctile(st_mat,50)';
ds_10 = prctile(st_mat,10)';
ds_90 = prctile(st_mat,90)';

%% Figures
figure; clf;
hold on
l1=plot(1:opt.tt-burn, dY_10,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(1:opt.tt-burn, dY_25,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(1:opt.tt-burn, dY_75,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(1:opt.tt-burn, dY_90,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');
legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Simulated GDP distribution in t','FontSize',16','Interpreter','Latex');
axis tight

figure; clf;
hold on
l1=plot(1:opt.tt-burn-12, dY_10_fut,'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','MS 10th');
l2=plot(1:opt.tt-burn-12, dY_25_fut,'Color',colors(15,:),'LineWidth', 2,'DisplayName','MS 25th');
l3=plot(1:opt.tt-burn-12, dY_75_fut,'Color',colors(45,:),'LineWidth', 2,'DisplayName','MS 75th');
l4=plot(1:opt.tt-burn-12, dY_90_fut,'Color',colors(55,:),'LineWidth', 1.5,'DisplayName','MS 90th');

legend([l1 l2 l3 l4],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Simulated GDP distribution 1-year ahead','FontSize',16','Interpreter','Latex');
axis tight

% States
figure; clf;
hold on
plot(1:opt.tt-burn, mean(st_mat),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','States');
set(gca,'children',flipud(get(gca,'children')))
axis tight
ylim([1 2]);
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Average states','FontSize',16','Interpreter','Latex');


figure; clf;
hold on
l1=plot(1:opt.tt-burn, mean(p12_draw),'Color',colors(5,:),'LineWidth', 1.5,'DisplayName','Good to bad');
l2=plot(1:opt.tt-burn, mean(p21_draw),'Color',colors(15,:),'LineWidth', 1.5,'DisplayName','Bad to good');
legend([l1 l2],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
set(gca,'children',flipud(get(gca,'children')))
axis tight
ylim([0 1]);
hold off
ylabel('Percent','interpreter','Latex','fontsize',10)
title('Switching probabilities','FontSize',16','Interpreter','Latex');

% figure;plot([y_mat(SEED,:)' f_mat(SEED,:)' m_mat(SEED,:)']);xlim([1 500])
% figure;plot([y_mat_fut(SEED,:)' y_mat_fut_test(SEED,:)']);
% figure;plot([y_mat(SEED,:)' y_mat_test(SEED,:)']);
%
% figure;plot([p12_draw(SEED,:)' p21_draw(SEED,:)']);


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
if opt.save_data ==1
    if sim_alg ==1
        save Res_sim_alg1 Res_iterated Res_direct dgp
    else
        save Res_sim_alg2 Res_iterated Res_direct dgp
    end
end

%% Load RISE
rise_startup()

%% Take a random seed (555)
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
        figfolder   = ['results/Figures/' opt.paramsUse '/Simul/Direct/' spectag '/'];
        
    elseif opt.dir_or_it ==2
        fcstype= 'iterated';
        varlist={'FF','MF','GDPG'}';
        % Folders
        paramsfolder = ['results/Params/' opt.paramsUse '/Simul/Iterated/'];
        figfolder   =  ['results/Figures/' opt.paramsUse '/Simul/Iterated/' spectag '/'];
    end
    quantfolder   = ['results/Quantiles/' opt.paramsUse '/Simul/'];
    
    
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
    if exist(quantfolder,'dir')==0;  mkdir(quantfolder); end
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
    outparams = scriptParams_FULL(params_in,1,modelspec);
    
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
        [quantiles,pdensity] = fPredDensityDirectFull(sv,params_in,FF_full,MF_full,GDPGH_full,TR_full,opt,modelspec);
        
    elseif opt.dir_or_it==2
        % Iterated forecast
        [quantiles, pdensity] = fPredDensityIteratedFull(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,'forecastst',modelspec);
        %         [quantiles, pdensity] = fPredDensityIteratedFull_debug(sv,params_in,FF_full,MF_full,GDPG_full,TR_full,opt,'forecastst',modelspec);
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
    %     if bad==1
    %         p_reg_bad_filtered = p_reg1_filtered;
    %         p_reg_bad_smoothed = p_reg1_smoothed;
    %     elseif bad==2
    p_reg_bad_filtered = p_reg2_filtered;
    p_reg_bad_smoothed = p_reg2_smoothed;
    %     end
    
    
    %%
    % ==========================================================================
    % Quantile Plots
    % ==========================================================================
    
    st_t_mean_temp = ones(size(quantiles.st_t_mean))-quantiles.st_t_mean;
    st_dgp = st_mat(SEED,:); st_dgp(st_dgp==2) = 0; st_dgp = ones(size(st_dgp))-st_dgp;
    
    figureTitleTag = [strrep(spectag,'_',' ') ', ' opt.paramsUse ', ' strrep(sheetuse,'_',' ') ', h=' num2str(opt.hh) ', ' startdb ':' enddb];
    
    % Add additional objects to structure
    fig1=figure('Units','normalized','Position', [0,0,1.2,0.9]); clf;
    subplot(211)
    hold on
    l1=plot(dates_full(sd:ed), quantiles.dYsim_10(sd:ed),'-','Color',colors(10,:),'LineWidth', 2,'DisplayName','10th');
    l2=plot(dates_full(sd:ed), quantiles.dYsim_90(sd:ed),'-','Color',colors(50,:),'LineWidth', 2,'DisplayName','90th');
    l9=plot(dates_full(sd:ed), quantiles.mean(sd:ed),'b-','LineWidth', 3,'DisplayName','Median');
    l10=plot(dates_full(sd:ed), y_mat_fut(SEED,sd:ed),'k-','LineWidth', 3,'DisplayName','dgp');
    legend([l1 l2 l9 l10],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
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
    l11=plot(dates_full(sd:ed), st_dgp(:,sd:ed),'k-','LineWidth', 2,'DisplayName','dgp');
    l3=plot(dates_full(sd:ed), st_t_mean_temp(sd:ed),'-','Color',colors(10,:),'LineWidth', 2.5,'DisplayName','Simulated s(t)');
    l4=plot(dates_full(sd:ed), p_reg_bad_filtered(sd:ed),'--','Color',colors(20,:),'LineWidth', 2.5,'DisplayName','Filtered s(t)');
    legend([l3 l4 l11],'Orientation','Vertical','Location','Best','interpreter','Latex');legend boxoff;
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
    set(gca, 'FontName', 'Times New Roman');    set(gca, 'FontSize', FontSize);
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

%% Quantiles comparison
if strcmp(opt.paramsUse,'mode')
    if sim_alg ==1
        if opt.dir_or_it ==1
            if exist([quantfolder modelname 'quantiles_direct_alg1.mat'],'file')==0 || opt.force_estimate==1
                quantiles_direct_alg1 = quantiles;
                save([quantfolder modelname 'quantiles_direct_alg1.mat'],'quantiles_direct_alg1')
            end
        elseif opt.dir_or_it ==2
            if exist([quantfolder modelname 'quantiles_iterated_alg1.mat'],'file')==0 || opt.force_estimate==1
                quantiles_iterated_alg1 = quantiles;
                save([quantfolder modelname 'quantiles_iterated_alg1.mat'],'quantiles_iterated_alg1')
            end
            if exist([paramsfolder modelname 'quantiles_iterated_dgp_alg1.mat'],'file')==0 || opt.force_estimate==1
                p_reg_good = st_mat(SEED,:);
                for rr=1:length(p_reg_good); if p_reg_good(rr) ==2; p_reg_good(rr) =0; end; end
                params_in_dgp = dgp;
                params_in_dgp.a0_1_2_sync_1 = (-1)*params_in_dgp.a0_1_2_sync_1;
                params_in_dgp.a0_1_2_sync_2 = (-1)*params_in_dgp.a0_1_2_sync_2;
                params_in_dgp.a0_2_1_sync_1 = (-1)*params_in_dgp.a0_2_1_sync_1;
                params_in_dgp.a0_2_1_sync_2 = (-1)*params_in_dgp.a0_2_1_sync_2;
                params_in_dgp.a0_3_1_sync_1 = (-1)*params_in_dgp.a0_3_1_sync_1;
                params_in_dgp.a0_3_1_sync_2 = (-1)*params_in_dgp.a0_3_1_sync_2;
                params_in_dgp.a0_3_2_sync_1 = (-1)*params_in_dgp.a0_3_2_sync_1;
                params_in_dgp.a0_3_2_sync_2 = (-1)*params_in_dgp.a0_3_2_sync_2;
                opt_dgp = opt; opt_dgp.nParamDraws = 1;
                [quantiles_iterated_dgp_alg1, pdensity_dgp] = fPredDensityIteratedFull_dgp(params_in_dgp,FF_full,MF_full,GDPG_full,TR_full,p_reg_good,opt_dgp,'forecastst',4);
                quantiles_iterated_dgp_alg1.GDPG = y_mat(SEED,:); quantiles_iterated_dgp_alg1.GDPGH = y_mat_fut(SEED,:);
                save([quantfolder modelname 'quantiles_iterated_dgp_alg1.mat'],'quantiles_iterated_dgp_alg1')
            end
        end
    elseif sim_alg ==2
        if opt.dir_or_it ==1
            if exist([quantfolder modelname 'quantiles_direct_alg2.mat'],'file')==0 || opt.force_estimate==1
                quantiles_direct_alg2 = quantiles;
                save([quantfolder modelname 'quantiles_direct_alg2.mat'],'quantiles_direct_alg2')
            end
        elseif opt.dir_or_it ==2
            if exist([quantfolder modelname 'quantiles_iterated_alg2.mat'],'file')==0 || opt.force_estimate==1
                quantiles_iterated_alg2 = quantiles;
                save([quantfolder modelname 'quantiles_iterated_alg2.mat'],'quantiles_iterated_alg2')
            end
            if exist([quantfolder modelname 'quantiles_iterated_dgp_alg2.mat'],'file')==0 || opt.force_estimate==1
                p_reg_good = st_mat(SEED,:);
                for rr=1:length(p_reg_good); if p_reg_good(rr) ==2; p_reg_good(rr) =0; end; end
                params_in_dgp = dgp;
                params_in_dgp.a0_1_2_sync_1 = (-1)*params_in_dgp.a0_1_2_sync_1;
                params_in_dgp.a0_1_2_sync_2 = (-1)*params_in_dgp.a0_1_2_sync_2;
                params_in_dgp.a0_2_1_sync_1 = (-1)*params_in_dgp.a0_2_1_sync_1;
                params_in_dgp.a0_2_1_sync_2 = (-1)*params_in_dgp.a0_2_1_sync_2;
                params_in_dgp.a0_3_1_sync_1 = (-1)*params_in_dgp.a0_3_1_sync_1;
                params_in_dgp.a0_3_1_sync_2 = (-1)*params_in_dgp.a0_3_1_sync_2;
                params_in_dgp.a0_3_2_sync_1 = (-1)*params_in_dgp.a0_3_2_sync_1;
                params_in_dgp.a0_3_2_sync_2 = (-1)*params_in_dgp.a0_3_2_sync_2;
                opt_dgp = opt; opt_dgp.nParamDraws = 1;
                [quantiles_iterated_dgp_alg2, pdensity_dgp] = fPredDensityIteratedFull_dgp(params_in_dgp,FF_full,MF_full,GDPG_full,TR_full,p_reg_good,opt_dgp,'forecastst',4);
                quantiles_iterated_dgp_alg2.GDPG = y_mat(SEED,:); quantiles_iterated_dgp_alg2.GDPGH = y_mat_fut(SEED,:);
                save([quantfolder modelname 'quantiles_iterated_dgp_alg2.mat'],'quantiles_iterated_dgp_alg2')
            end
        end
    end
    if quantiles_fig ==1
        if sim_alg ==2
            if exist([quantfolder modelname 'quantiles_direct_alg2.mat'],'file')==2 && exist([quantfolder modelname 'quantiles_iterated_alg2.mat'],'file')==2 && exist([quantfolder modelname 'quantiles_iterated_dgp_alg2.mat'],'file')==2
                run quantile_comparison_figures.m
            else
                fprintf('\n Run Direct and Iterated before constructing comparisons. \n')
            end
        else
            fprintf('\n Change to Algorithm 2. \n')
        end
    end
end

%%
% rise_exit
