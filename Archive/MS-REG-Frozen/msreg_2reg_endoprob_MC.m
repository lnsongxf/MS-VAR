%% Markov-Switching Real-Time Growth-at-Risk
% Model with f_t and m_t endogenous
% Author: Francesca Loria
% This Version: June 2020

%% housekeeping
clear; close all; clc; tic;

% Important paths
addpath('/if/prod-tfs/production/GAR/MS-VAR/RISE_toolbox');
addpath(genpath('scripts'));
addpath(genpath('cbrewer'));

% Options
dir_or_it = 1; % 1 = direct, 2 = iterated
saveit    = 0; % 1 = save graphs
save_mcmc = 1; % 1 = save posterior sampling results
const     = 1; % 1= have a constant in transition probability
normal    = 0; % 1 = use normal distribution, 0 = gamma distribution

% Data vintage, sample and country selection
start_date   = '2100-Jan';
end_date     = '2147-Sep';
end_date_db  = '2146-Sep';

% VAR configuration 
nlags=1;
exog={};
constant=true;
panel=[];

% Create date formats for plotting
inputformat = 'yyyy-MMM';
dataformat  = 'yyyy-mmm';
start_plot    = '2100-Feb';
end_plot      = end_date_db;
end_plot_full = end_date;

% Vector of dates for the full sample
dates_full = datenum((datetime(start_date,'InputFormat',inputformat))+calmonths(nlags):calmonths(1):(datetime(end_date,'InputFormat',inputformat)))';

% Index of dates
sd      = find(datenum(start_plot,dataformat)==dates_full);
ed      = find(datenum(end_plot,dataformat)==dates_full);
ed_full = find(datenum(end_plot_full,dataformat)==dates_full);

% Vector of dates for plotting
dates   = dates_full(sd:ed);

% Figure options
FontSize = 16;
numticks = 48;
figSize = [12 6];
linestyle = {'-','--',':'};
colors = cbrewer('div', 'RdYlBu', 64);
colors2 = cbrewer('qual', 'Set1', 8);
left_color = [0 0 0];
right_color = colors(2,:);

%% Load RISE
rise_startup()

%% Prepare Data

%load data
load('Res_sim');
Res.Res_iterated = Res_iterated; Res.Res_direct = Res_direct;

% Number of simulated paths
n_draws = size(Res.Res_iterated.y_mat,1);

% Pick a particular simulation
draw =10;

[db, db_full, tex] = data_sim(Res,draw,start_date,end_date_db);
nlags=1;
% Collect MF, FF and trend (estimation sample)
FF = db.FF.data(nlags+1:end);
MF = db.MF.data(nlags+1:end);

% Collect MF, FF and trend (full sample)
FF_full = db_full.FF.data(nlags+1:end);
MF_full = db_full.MF.data(nlags+1:end);

if dir_or_it ==1
    varlist={'FF','MF','GDPGH'}';
elseif dir_or_it ==2
    varlist={'FF','MF','GDPG'}';
end

%% set up Markov chains

switches = {'c(3)','a0(3)','s(3)'};
model = 'SVAR_StateDependent_2Regimes';

if const==1
    if normal==1
        prob_fct = {{

        % Format is [name_tp_1_2]
        'sync_tp_1_2=1/(1+exp(a12-b12*(FF)-c12*(MF)))'
        'sync_tp_2_1=1/(1+exp(a21-b21*(FF)-c21*(MF)))'
        }};
    elseif normal==0
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))'
            }};
    end
    prob_params = {{'a12','a21','b12','c12','b21','c21'}};
elseif const==0
    if normal==1
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)-c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(-b21*(FF)-c21*(MF)))'
            }};
    elseif normal==0
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)+c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(+b21*(FF)-c21*(MF)))'
            }};
    end
    prob_params = {{'b12','c12','b21','c21'}};
end

markov_chains=struct('name','sync',...
    'number_of_states',2,...
    'controlled_parameters',{switches},...   % these correspond to the parameters that are switching
    'endogenous_probabilities',prob_fct,...
    'probability_parameters',prob_params);


%% Create the VAR
nlags=1;

exog={};

constant=true;

panel=[];

sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% set up restrictions

% syntax is alag(eqtn,vname)
%-------------------------------
if dir_or_it ==1
lin_restr={
    % first equation or "financial factor" equation
    %----------------------------------
    'a0(1,GDPGH)=0'
    'a1(1,GDPGH)=0'
    %'a1(1,MF)=0'
    % second equation or "macroeconomic factor" equation
    %-----------------------------------
    'a0(2,GDPGH)=0'
    'a0(2,FF)=0'
    'a1(2,GDPGH)=0'
    % third equation or "GDP Growth" equation
    %----------------------------------
    'a1(3,GDPGH)=0'
    'a1(3,FF)=0'
    'a1(3,MF)=0'
    };
elseif dir_or_it ==2
    lin_restr={
    % first equation or "financial factor" equation
    %----------------------------------
    'a0(1,GDPG)=0'
    'a1(1,GDPG)=0'
    %'a1(1,MF)=0'
    % second equation or "macroeconomic factor" equation
    %-----------------------------------
    'a0(2,GDPG)=0'
    'a0(2,FF)=0'
    'a1(2,GDPG)=0'
    % third equation or "GDP Growth" equation
    %----------------------------------
    'a1(3,GDPG)=0'
    'a1(3,FF)=0'
    'a1(3,MF)=0'
    };
end

restrictions=[lin_restr];

%% set priors

% priors for the VAR coefficients
var_prior=svar.prior_template();
%--------------------------------
% Sims and Zha (1998) Normal-Wishart Prior Hyperparameters
% L1 : Overall tightness
% L2 : Cross-variable specific variance parameter
% L3 : Speed at which lags greater than 1 converge to zero
% L4 : tightness on deterministic/exogenous terms
% L5 : covariance dummies(omega)
% L6 : co-persistence (lambda)
% L7 : Own-persistence (mu)
var_prior.type='sz';
var_prior.L1=1;
var_prior.L5=1;
var_prior.L6=0; % cointegration
var_prior.L7=0; % unit root
% var_prior.coefprior=0.9; % impose 0.9 eigenvalue for both

% priors for the sync transition probabilities
%------------------------------------------------
switch_prior=struct();
if const==1
    if normal==1
        switch_prior.a12={0,0,2,'normal'};
        switch_prior.b12={0,0,2,'normal'};
        switch_prior.c12={0,0,2,'normal'};
        switch_prior.a21={0,0,2,'normal'};
        switch_prior.b21={0,0,2,'normal'};
        switch_prior.c21={0,0,2,'normal'};
    elseif normal==0
        switch_prior.a12={0.5,0.5,0.5,'normal'};
        switch_prior.a21={0.5,0.5,0.5,'normal'};
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
    end
elseif const==0
    if normal==1
        switch_prior.b12={-0.25,-0.25,0.1,'normal'};
        switch_prior.c12={0.25,0.25,0.1,'normal'};
        switch_prior.b21={0.25,0.25,0.1,'normal'};
        switch_prior.c21={-0.25,-0.25,0.1,'normal'};
    elseif normal==0
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
    end
end

prior=struct();

prior.var=var_prior;

prior.nonvar=switch_prior;

%% Find posterior mode

sv=estimate(sv0,db,{'2100M1','2146M09'},prior,restrictions,'fmincon');

%% Store estimates

pmode=posterior_mode(sv)


print_structural_form(sv)
% print_solution(sv)

