
%==========================================================================
% Posterior sampling
%==========================================================================

%========================================
% SAMPLER OPTIONS
%========================================

options = struct();
options.alpha = 0.234;
options.thin = 10;
% options.burnin = 10^3;
% options.N = 2*10^4;
options.burnin = 10^1;
options.N = 2*10^2;
options.nchain = 1;

% Estimate or load existing posterior mode
if exist('sPMode.mat','file')==0
    % Estimate posterior mode
    error('First run the code with paramsUse = mode')
else
    % Load stored results
    load('sPMode.mat');    
    % Map objects:
    sv = sPMode.sv;
    prior = sPMode.prior;
    restrictions = sPMode.restrictions;
end


% Get posterior mode
pmode=posterior_mode(sv);

pnames=fieldnames(pmode);

a2tilde_to_a=sv.estim_.linres.a2tilde_to_a;

[ff,lb,ub,x0,vcov,self]=pull_objective(sv);

% Lower bounds
lb(~isfinite(lb))=-500;
ub(~isfinite(ub))=500;

% Run mcmc sampler
results=mh_sampler(ff,lb,ub,options,x0,vcov);

% Collect Posterior Estimates and Confidence Bands

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

%% Structure with Posterior Model and MCMC parameters
sPmcmc = sPMode;
sPmcmc.options = options;
sPmcmc.params = pmcmc;
sPmcmc.results = results;

%% Print Table of Estimates
% tableTitle = ['Posterior/Distribution:/' fcstype ',/' sheetuse ',/h/=/' num2str(opt.hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
% pmcmc_tex(pmode,paramsCI_mcmc,texfolder,modelname,tableTitle);
% 
% tableTitleMode = ['Posterior/Mode/Distribution:/' fcstype ',/' sheetuse ',/h/=/' num2str(opt.hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
% pmcmc_tex(pmode,paramsCI_mode,texfolder,[modelname '_mode'],tableTitleMode);


