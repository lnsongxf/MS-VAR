%% --------------------- constant-parameter model --------------------- %%

%% Housekeeping
clear; clc; tic;
close all;

% RISE Toolbox
addpath('/if/gms_research_gar/GAR/MS-VAR/RISE_toolbox');         

% Load RISE
rise_startup()

%% Opts
opt.mp_restriction = 1;      % 1 = restrictions of as in Arias et al. (2019), 2 = coeffs fixed at their estimations
opt.paramsUse      = 'mode'; % mode = posterior mode, mcmc = parameter uncertainty
opt.force_estimate = 1;      % 1 = force to reestimate model

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
if opt.mp_restriction ==1
    lin_restr={
        % Restriction 1
        % systematic component of MP
        'a0(6,TR)=0'
        'a0(6,NBR)=0'
        };
    nonlin_restr={
        % Restriction 2
        % systematic component of MP
        'a0(6,GDP)<0' % Negative, as it is on the LHS
        'a0(6,PGDP)<0' % Negative, as it is on the LHS
        };
elseif opt.mp_restriction ==2
    lin_restr={
        % Restriction 1
        % systematic component of MP
        'a0(6,TR)=0'
        'a0(6,NBR)=0'
        'a0(6,GDP)=-0.84' % Negative, as it is on the RHS
        'a0(6,PGDP)=-2.73' % Negative, as it is on the RHS
        };
    nonlin_restr={};
end

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
if strcmp(opt.paramsUse,'mode')
    if exist('sPMode.mat','file')==0 || opt.force_estimate==1
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
else
    %%%%%%%%
    % MCMC %
    %%%%%%%%
    
    % Load or estimate mcmc
    if exist('sPMCMC.mat','file')==0 || opt.force_estimate==1
        % Estimate posterior distribution
        run scriptEstimateMCMC.m
        
        % Store solution structure
        save('sPMCMC.mat','sPmcmc')
        
    else
        
        % Load stored results
        load('sPMCMC.mat');
        
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
        pmcmc_tex(pmode,paramsCI_mcmc,texfolder,modelname,'MCMC');
        
        pmcmc_tex(pmode,paramsCI_mode,texfolder,[modelname '_mode'],'Mode');
        
        
    end
        
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

ci=[16,50,68];

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
%     d=myirfs.(shk);
%     out=fanchart(d,ci);
%     plot_fanchart(out)
    set(ax,'XTick',0:12:60)
    set(ax,'XTickLabel',0:5)
    xlabel('Years')
    title(v)
end
