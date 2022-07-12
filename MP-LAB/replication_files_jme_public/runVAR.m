%% --------------------- constant-parameter model --------------------- %%

%% Housekeeping
clear; clc; tic;
close all;

% RISE Toolbox
addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');         

% Load RISE
rise_startup()

%% Create dataset

[db,varlist0]=createData;
varlist=fieldnames(varlist0);

%% Create the VAR

clc
% from JME code:
% nlag = 12; % number of lags
% nex  = 0;  % set equal to 1 if a constant is included; 0 otherwise
nlags=12;
exog={};
constant=false;
sv0=svar(varlist,exog,nlags,constant);

%% set up restrictions

% syntax is alag(eqtn,vname)
%-------------------------------
lin_restr={
    % Restriction 1
    % systematic component of MP 
    'a0(6,TR)=0'
    'a0(6,NBR)=0'
    };
nonlin_restr={
    % Restriction 2
    % systematic component of MP 
    'a0(6,GDP)<0'
    'a0(6,PGDP)<0'
    };

restrictions=[lin_restr;nonlin_restr];

%% set priors 

var_prior=svar.prior_template();
var_prior.type='sz';
prior=struct('var',var_prior);

%% Find posterior mode
clc
sv=sv0;
startdb = '1965M1';
enddb   = '2007M6';
if exist('sPMode.mat','file')==0
    sv=estimate(sv,db,{startdb,enddb},prior,restrictions);
    % sv=estimate(sv,db,{startdb,enddb},prior,restrictions,'fmincon');
    
    % Structure to store posterior mode results
    sPMode.prior        = prior;
    sPMode.restrictions = restrictions;
    sPMode.startdb      = startdb;
    sPMode.enddb        = enddb;
    sPMode.sv           = sv;
    sPMode.nlags        = nlags;
    sPMode.exog         = exog;
    sPMode.constant     = constant;
    
    % Store solution structure
    save('sPMode.mat','sPMode');
else
    load('sPMode.mat');
    sv = sPMode.sv;
end

%% Printing estimates
clc

% Printing structural form
print_structural_form(sv)

% % Printing reduced form
% print_solution(sv)

%% Impulse responses

myirfs=irf(sv,[],61);

snames=fieldnames(myirfs);

figure('name','Simple impulse response functions')
iter=0;
% focus on MP shock
shk=snames{end};
for iv=1:numel(varlist)        
    v=varlist{iv};        
    iter=iter+1;       
    subplot(3,2,iter)
    ax = gca;
    plot(myirfs.(shk).(v),'linewidth',2) 
    set(ax,'XTick',0:12:60)
    set(ax,'XTickLabel',0:5)
    xlabel('Years')
    title(v)
end
