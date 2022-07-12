
% This script specified the model structure, switching parameters
% and estimates the posterior given MF and FF
%==================================================================

%% set up Markov chains
% Model specifications as described in SZ-VAR.pdf
if modelspec==1                    
    switches = {'c(3)','a0(3)','s(3)'};
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
    switches = {'c(3)','a0(3)','a1(3)'};
elseif modelspec==8
    switches = {'c(3)','a0(3)','a1(3)','s(3)'};
end


%% Transition probabilities
if opt.const==1
    if opt.normal==1
        prob_fct = {{
            
        % Format is [name_tp_1_2]
        'sync_tp_1_2=1/(1+exp(a12-b12*(FF)-c12*(MF)))'
        'sync_tp_2_1=1/(1+exp(a21-b21*(FF)-c21*(MF)))'
        }};
    elseif opt.normal==0
        % This is benchmark specification
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(a12-b12*(FF)+c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(a21+b21*(FF)-c21*(MF)))'
            }};
    end
    prob_params = {{'a12','a21','b12','c12','b21','c21'}};
elseif opt.const==0
    if opt.normal==1
        prob_fct = {{
            'sync_tp_1_2=1/(1+exp(-b12*(FF)-c12*(MF)))'
            'sync_tp_2_1=1/(1+exp(-b21*(FF)-c21*(MF)))'
            }};
    elseif opt.normal==0
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
sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% set up restrictions

% syntax is alag(eqtn,vname)
%-------------------------------
if opt.dir_or_it ==1
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Direct Forecast Model %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Restrictions for A0 when constant
    %----------------------------------
    a0_const = {
            'a0(1,GDPGH)=0' % No contemporaneous effect of GDP on FF
            'a0(2,GDPGH)=0' % No contemporaneous effect of GDP on mF
            'a0(2,FF)=0'    % No contemporaneous effect of FF on MF
%             'a0(1,MF)=0'    % No contemporaneous effect of MF on FF
            };
        

    % Restrictions for A1 when constant
    %----------------------------------
    a1_const = {
        %             'a1(1,GDPGH)=0' % No effects of GDP(t-1) on FF
        %             'a1(2,GDPGH)=0' % No effect of GDP(t-1) on MF
        %             'a1(3,FF)=0' % No effect of FF(t-1) on GDP
        %             'a1(3,MF)=0' % No effect of MF(t-1) on GDP
            'a1(3,GDPGH)=0' % No effect of GDP(t-1) on GDP
            };
        
    % Restrictions for A0 and A1 when switching
    %------------------------------------------
    numberOfStates=markov_chains.number_of_states;
    a0_markov=cell(0,1);
    a1_markov=cell(0,1);

    for istate=1:numberOfStates
            
            mystate=int2str(istate);
            
            a0_markov=[a0_markov
                {
                ['a0(1,GDPGH,sync,',mystate,')=0'] % No contemporaneous effect of GDP on FF
                ['a0(1,MF,sync,',mystate,')=0'] % No contemporaneous effect of MF on FF
                ['a0(2,GDPGH,sync,',mystate,')=0']  % No contemporaneous effect of GDP on MF        
%                 ['a0(2,FF,sync,',mystate,')=0']    % No contemporaneous effect of FF on MF      
                }];
            
            a1_markov=[a1_markov
                {
                ['a1(3,GDPGH,sync,',mystate,')=0']  % No effect of GDP(t-1) on GDP   
%                 ['a1(3,FF,sync,',mystate,')=0']    % No effect of FF(t-1) on GDP   
%                 ['a1(3,MF,sync,',mystate,')=0']    % No effect of MF(t-1) on GDP 
                }];

    end

    if imodel==1
        lin_restr=[a0_const;a1_const];
    elseif imodel==2
        lin_restr=[a0_markov;a1_const];
    elseif imodel==3 || imodel==4
        lin_restr=[a0_markov;a1_markov];
    elseif imodel==5
        lin_restr=[a0_const;a1_markov];
    elseif imodel==6
        lin_restr=[a0_markov;a1_const];
    elseif imodel==7
        error('Model 7 not meaningful in direct approach with a1(3,GDPGH)=0 restrictions in both states');
    elseif imodel==8
        error('Model 8 not meaningful in direct approach with a1(3,GDPGH)=0 restrictions in both states');
    end

elseif opt.dir_or_it ==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterated Forecast Model %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Restrictions for A0 when constant
    %----------------------------------
    a0_const = {
            'a0(1,GDPG)=0' % No contemporaneous effect of GDP on FF        
%             'a0(1,MF)=0'   % No contemporaneous effect of MF on FF       
            'a0(2,GDPG)=0'  % No contemporaneous effect of GDP on MF              
            'a0(2,FF)=0'  % No contemporaneous effect of FF on MF             
            };
        
    % Restrictions for A0 when switching        
    %----------------------------------
    numberOfStates=markov_chains.number_of_states;
    a0_markov=cell(0,1);
    a1_markov=cell(0,1);

    for istate=1:numberOfStates
            
            mystate=int2str(istate);
%             
%             a0_markov=[a0_markov
%                 {
%                 ['a0(1,GDPG,sync,',mystate,')=0'] % No contemporaneous effect of GDP on FF  
% %                 ['a0(1,MF,sync,',mystate,')=0'] % No contemporaneous effect of MF on FF  
%                 ['a0(2,GDPG,sync,',mystate,')=0']  % No contemporaneous effect of GDP on MF       
%                 ['a0(2,FF,sync,',mystate,')=0']  % No contemporaneous effect of FF on MF     
%                 }];

            % Cholesky on A0
            a0_markov=[a0_markov
                {
                ['a0(1,GDPG,sync,',mystate,')=0'] % No contemporaneous effect of GDP on FF  
                ['a0(1,MF,sync,',mystate,')=0'] % No contemporaneous effect of MF on FF  
                ['a0(2,GDPG,sync,',mystate,')=0']  % No contemporaneous effect of GDP on MF       
                %['a0(2,FF,sync,',mystate,')=0']  % No contemporaneous effect of FF on MF     
                }];
            

            a1_markov=[a1_markov
            {
            ['a1(3,GDPG,sync,',mystate,')=0']  % No effect of GDP(t-1) on GDP   
%             ['a1(3,FF,sync,',mystate,')=0']    % No effect of FF(t-1) on GDP   
%             ['a1(3,MF,sync,',mystate,')=0']    % No effect of MF(t-1) on GDP 
            }];

    end
       
    if imodel==1 || imodel==5 || imodel==7  || imodel==8
        lin_restr=a0_const;
    elseif imodel==3 || imodel==4
        lin_restr=[a0_markov;a1_markov];
    else          
        lin_restr=a0_markov;
    end
    
end

restrictions=lin_restr;

%% set priors

% priors for the VAR coefficients
var_prior=svar.prior_template();
%--------------------------------
% Sims and Zha (1998) normal-Wishart Prior Hyperparameters
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

