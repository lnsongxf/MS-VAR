
% This script specified the model structure, switching parameters
% and estimates the posterior given MF and FF
%==================================================================

%% set up Markov chains

switches = {'c(3)','a0(3)','s(3)'};

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
if vartype ==0
    if dir_or_it ==1
        % Direct Forecast Model
        lin_restr={
            % first equation or "financial factor" equation
            %----------------------------------
            'a0(1,GDPGH)=0' % No contemporaneous effect of GDP
            'a1(1,GDPGH)=0' % No lagged effects of GDP
            %'a1(1,MF)=0'

            % second equation or "macroeconomic factor" equation
            %-----------------------------------
            'a0(2,GDPGH)=0' % No contemporaneous effect of GDP
            'a0(2,FF)=0'    % No contemporaneous effect of FF
            'a1(2,GDPGH)=0' % No lagged effect of GDP

            % third equation or "GDP Growth" equation
            %----------------------------------
            'a1(3,GDPGH)=0' % No lagged effect of GDP
            'a1(3,FF)=0'    % No lagged effect of FF
            'a1(3,MF)=0'    % No lagged effect of MF
            };
    elseif dir_or_it ==2
        % Iterated Forecast Model
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
else
    % Impose Cholesky only
    if dir_or_it ==1
        % Direct Forecast Model
        lin_restr={
            % first equation or "financial factor" equation
            %----------------------------------
            'a0(1,GDPGH)=0' % No contemporaneous effect of GDP
            'a0(1,MF)=0' % No contemporaneous effect of MF

            % second equation or "macroeconomic factor" equation
            %-----------------------------------
            'a0(2,GDPGH)=0' % No contemporaneous effect of GDP

            };
    elseif dir_or_it ==2
        % Iterated Forecast Model
        lin_restr={
            % first equation or "financial factor" equation
            %----------------------------------
            'a0(1,GDPG)=0'  % No contemporaneous effect of GDP
            'a0(1,MF)=0' % No contemporaneous effect of MF

            % second equation or "macroeconomic factor" equation
            %-----------------------------------
            'a0(2,GDPG)=0'

            };
    end    
end

restrictions = lin_restr;

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
tableTitle = ['Posterior/Mode:/' fcstype ',/' sheetuse ',/h/=/' num2str(hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
pmode_tex(posterior_mode(sv),texfolder,modelname,tableTitle);

