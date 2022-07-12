%% Caldara and Herbst Proxy SVAR in Structural Form
% Author: Francesca Loria
% This Version: March 2022

%% Housekeeping
clear
close all
% clc()

% addpath('../RISE_Tbx_Beta-master');
addpath('../RISE_toolbox');
addpath(genpath('scripts'));

%% Load RISE

rise_startup()

%% Load the data

load data

%% Select model and sample
model_type=0;
startit = '1994M1';
endit = '2007M6';
nlags=12;



%% Create the VAR

varlist = {'EFFR','LIPM','UNRATE','LPPI','BAA10YMOODY','MHF'};

exog={};

constant=true;

panel=[];

[restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains0(nlags);
% [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains(nlags);

sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% Set priors

% first we create a template structure
sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

% priors for the VAR coefficients
var_prior=svar.prior_template();
%--------------------------------
% Sims and Zha (1998) opt.normal-Wishart Prior Hyperparameters
% L1 : Overall tightness
% L2 : Cross-variable specific variance parameter
% L3 : Speed at which lags greater than 1 converge to zero
% L4 : tightness on deterministic/exogenous terms
% L5 : covariance dummies(omega)
% L6 : co-persistence (lambda)
% L7 : Own-persistence (mu)
var_prior.type='sz';
var_prior.L6=10; % cointegration (10)
var_prior.L7=10; % unit root  (10)
% var_prior.coefprior=0.9; % impose 0.9 eigenvalue for both

prior=struct('var',var_prior,'nonvar',switch_prior);

is_prior=true;

if ~is_prior
    
    prior=rmfield(prior,'var');
    
end

%% Find posterior mode

sv=estimate(sv0,db,{startit,endit},prior,restrictions,'fmincon');

%% Estimates

pmode=posterior_mode(sv)

%% Printing solution

%structural form
print_structural_form(sv)

% reduced form
print_solution(sv)

%% Impulse responses to MP Shock

N = length(varlist);

shock_names = [];
myirfs=irf(sv,shock_names,48);

snames=fieldnames(myirfs);

figure('name','Simple impulse response functions')

iter=0;

jter=0;

for ishk=1 % MP shock is ordered first
    
    shk=snames{ishk};
    
    for iv=1:N
        
        v=varlist{iv};
        
        iter=iter+1;
        
        subplot(N,1,iter)
        
        plot(100*myirfs.(shk).(v),'linewidth',2)
        
        title(v)
        
    end
    
end

