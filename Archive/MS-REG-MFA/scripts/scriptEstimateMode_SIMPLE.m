
% This script specified the model structure, switching parameters
% and estimates the posterior given MF and FF
%==================================================================

%% set up Markov chains
% Model specifications as described in SZ-VAR.pdf
if modelspec==1                    
    switches = {'c(2)','a0(2)','s(2)'};
elseif modelspec==2            
    switches = {'c','a0','s'};
elseif modelspec==3
    switches = {'c','a0','a1'};
elseif modelspec==4
    switches = {'c','a0','a1','s'};
elseif modelspec==5
    switches = {'c','a1','s'};
elseif modelspec==6
    switches = {'c','a0'};
elseif modelspec==7
    switches = {'c(2)','a0(2)','a1(2)'};
elseif modelspec==8
    switches = {'c(2)','a0(2)','a1(2)','s(2)'};
end


%% Transition probabilities
if opt.const==1
    if opt.normal==1
        if opt.dir_or_it ==1
        % Direct Forecast Model
        error('Direct model not allowed in this version')
        % Doesn't work because we cannot pass GDPG in a model that 
        % uses GDPGH and FF 
    elseif opt.dir_or_it ==2
         % Iterated Forecast Model
        prob_fct = {{
                'sync_tp_1_2=1/(1+exp(a12-b12*(FF)-c12*(GDPG)))'
                'sync_tp_2_1=1/(1+exp(a21-b21*(FF)-c21*(GDPG)))'
                }};       
        end
    elseif opt.normal==0
        % This is benchmark specification
        if opt.dir_or_it ==1
        error('Direct model not allowed in this version')
        % Doesn't work because we cannot pass GDPG in a model that 
        % uses GDPGH and FF 
        elseif opt.dir_or_it ==2
         % Iterated Forecast Model
            prob_fct = {{
                'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(GDPG)))'
                'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(GDPG)))'
                }};
        end
    end
    prob_params = {{'a12','a21','b12','c12','b21','c21'}};
elseif opt.const==0
    if opt.normal==1
        if dir_or_it ==1
        % Direct Forecast Model
        error('Direct model not allowed in this version')
        % Doesn't work because we cannot pass GDPG in a model that 
        % uses GDPGH and FF 
    elseif opt.dir_or_it ==2
         % Iterated Forecast Model
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)-c12*(GDPG)))'
            'sync_tp_2_1=1/(1+exp(-b21*(FF)-c21*(GDPG)))'
            }};
        end
    elseif opt.normal==0
        if opt.dir_or_it ==1
        % Direct Forecast Model
        error('Direct model not allowed in this version')
        % Doesn't work because we cannot pass GDPG in a model that 
        % uses GDPGH and FF 
    elseif opt.dir_or_it ==2
         % Iterated Forecast Model
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)+c12*(GDPG)))'
            'sync_tp_2_1=1/(1+exp(+b21*(FF)-c21*(GDPG)))'
            }};
        end
    end
    prob_params = {{'b12','c12','b21','c21'}};
end

markov_chains=struct('name','sync',...
    'number_of_states',2,...
    'controlled_parameters',{switches},...   % these correspond to the parameters that are switching
    'endogenous_probabilities',prob_fct,...
    'probability_parameters',prob_params);


%% Create the VAR
sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% set up restrictions

% syntax is alag(eqtn,vname)
%-------------------------------
if opt.dir_or_it ==1
    % Direct Forecast Model
    if nlags==1
        lin_restr={
            % first equation or "financial factor" equation
            %----------------------------------
            'a0(1,GDPGH)=0' % No contemporaneous effect of GDP
            'a1(1,GDPGH)=0' % No effects of GDP
            % second equation or "GDP Growth" equation
            %----------------------------------
            'a1(2,GDPGH)=0' % No effect of GDP (t-1)
            'a1(2,FF)=0'    % No effect of FF (t-1)     
            };
    elseif nlags==2
        lin_restr={
            % first equation or "financial factor" equation
            %----------------------------------
            'a0(1,GDPGH)=0' % No contemporaneous effect of GDP
            'a1(1,GDPGH)=0' % No effects of GDP
            'a2(1,GDPGH)=0' % No effect of GDP

            % second equation or "GDP Growth" equation
            %----------------------------------
            'a1(2,GDPGH)=0' % No effect of GDP (t-1)
            'a1(2,FF)=0'    % No effect of FF (t-1)     
            'a2(2,GDPGH)=0' % No effect of GDP (t-2)
            'a2(2,FF)=0'    % No effect of FF (t-2)
            };        
    end
elseif opt.dir_or_it ==2
    % Iterated Forecast Model
    if imodel==1 || imodel==5 || imodel==7  || imodel==8
       lin_restr={       
           % first equation or "financial factor" equation
            %----------------------------------
            'a0(1,GDPG)=0' % No contemporaneous effect of GDP        
            };
    else          
        numberOfStates=markov_chains.number_of_states;
        lin_restr=cell(0,1);
        for istate=1:numberOfStates

            mystate=int2str(istate);
                lin_restr=[lin_restr
                {
                % first equation or "FFR" equation:
                %----------------------------------
                ['a0(1,GDPG,sync,',mystate,')=0'] % No contemporaneous effect of GDP
                }];
        end
    end
end

restrictions=lin_restr;

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
if opt.const==1
    if opt.normal==1
        switch_prior.a12={0,0,2,'normal'};
        switch_prior.b12={0,0,2,'normal'};
        switch_prior.c12={0,0,2,'normal'};
        switch_prior.a21={0,0,2,'normal'};
        switch_prior.b21={0,0,2,'normal'};
        switch_prior.c21={0,0,2,'normal'};
    elseif opt.normal==0
        switch_prior.a12={0.5,0.5,0.5,'normal'};
        switch_prior.a21={0.5,0.5,0.5,'normal'};
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
    end
elseif opt.const==0
    if opt.normal==1
        switch_prior.b12={-0.25,-0.25,0.1,'normal'};
        switch_prior.c12={0.25,0.25,0.1,'normal'};
        switch_prior.b21={0.25,0.25,0.1,'normal'};
        switch_prior.c21={-0.25,-0.25,0.1,'normal'};
    elseif opt.normal==0
        switch_prior.b12={0.5,0.5,0.25,'gamma'};
        switch_prior.c12={0.5,0.5,0.25,'gamma'};
        switch_prior.b21={0.5,0.5,0.25,'gamma'};
        switch_prior.c21={0.5,0.5,0.25,'gamma'};
    end
end

prior=struct();

prior.var=var_prior;

prior.nonvar=switch_prior;

% Find posterior mode
sv=estimate(sv0,db,{startdb,enddb},prior,restrictions,'fmincon');

% Structure to store posterior mode results
sPMode.prior        = prior;
sPMode.restrictions = restrictions;
sPMode.startdb      = startdb;
sPMode.enddb        = enddb;
sPMode.sv           = sv;
sPmode.nlags        = nlags;
sPmode.exog         = exog;
SPmode.constant     = constant;
sPmode.panel        = panel;

%% Print posterior mode to csv an produce Latex Table
% Need to use / for spaces to parse the tex table in Latex
tableTitle = ['Posterior/Mode:/' fcstype ',/' sheetuse ',/h/=/' num2str(opt.hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
%pmode_tex(posterior_mode(sv),texfolder,modelname,tableTitle);

