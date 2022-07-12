
% This script specified the model structure, switching parameters
% and estimates the posterior given MF and FF
%==================================================================

%% set up Markov chains
% Model specifications as described in SZ-VAR.pdf
if modelspec==1 || modelspec==13             
    switches = {'c(3)','a0(3)','s(3)'};
elseif modelspec==2 || modelspec==10 || modelspec==11           
    switches = {'c','a0','s'};
elseif modelspec==3
    switches = {'c','a0','a1'};
elseif modelspec==4 || modelspec==12
    switches = {'c','a0','a1','s'};
elseif modelspec==5
    switches = {'c','a1','s'};
elseif modelspec==6
    switches = {'c','a0'};
elseif modelspec==7
    switches = {'c(3)','a0(3)','a1(3)'};
elseif modelspec==8 || modelspec==9
    switches = {'c(3)','a0(3)','a1(3)','s(3)'};
end


%% Transition probabilities
% Format is [name_tp_1_2]

if opt.const==1
    if opt.transprob==0
        prob_fct = {{            
            'sync_tp_1_2=a12'
            'sync_tp_2_1=a21'
        }};
        prob_params = {{'a12','a21'}};
    elseif opt.transprob==1
       if opt.normal==1
                prob_fct = {{            
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
    end
elseif opt.const==0
    if opt.transprob==0
        errordlg('opt.transprob=0 incompatible with opt.const=0');
    elseif opt.transprob==1
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
    
        
    % Restrictions for constant elements of A0
    %-----------------------------------------
    a0_const = {
            'a0(1,GDPGH)=0' % No contemporaneous effect of GDP on FF
            'a0(1,MF)=0'    % No contemporaneous effect of MF on FF
%             'a0(2,FF)=0'    % No contemporaneous effect of FF on MF
            'a0(2,GDPGH)=0' % No contemporaneous effect of GDP on mF
            };
        

    % Restrictions for constant elements of A1
    %-----------------------------------------
   a1_const=cell(0,1);
   for ilag=1:nlags
        mylag=int2str(ilag);
        a1_const=[a1_const
             {
            ['a',mylag,'(1,GDPGH)=0'] % No effect of GDP(t-i) on FF
            ['a',mylag,'(2,GDPGH)=0'] % No effect of GDP(t-i) on MF
%             ['a',mylag,'(3,GDPGH)=0'] % No effect of GDP(t-i) on GDP
            }];
   end
    if modelspec==10
       for ilag=1:nlags
            mylag=int2str(ilag);
            a1_const=[a1_const
                 {
                ['a',mylag,'(3,FF)=0'] % No effect of FF(t-i) on GDP
                ['a',mylag,'(3,MF)=0'] % No effect of MF(t-i) on GDP
%                 ['a',mylag,'(3,GDPGH)=0'] % No effect of GDP(t-i) on GDP
                }];
       end
    elseif modelspec==11 || modelspec==13
       for ilag=1:nlags
            mylag=int2str(ilag);
            a1_const=[a1_const
                 {
                ['a',mylag,'(3,GDPGH)=0'] % No effect of GDP(t-i) on GDP
                }];
       end
    elseif modelspec==9
       for ilag=1:nlags
            mylag=int2str(ilag);
            a1_const=[a1_const
                 {
                ['a',mylag,'(3,GDPGH)=0'] % No effect of GDP(t-i) on GDP
                }];
       end
    else
         a1_const = {
                };
    end
        
    % Restrictions for switching elements of A0 and A1, C and SIG
    %------------------------------------------------------------
    numberOfStates=markov_chains.number_of_states;
    a0_markov=cell(0,1);
    a1_markov=cell(0,1);
    a1_markov9=cell(0,1);

    % Restrictions for c(s_t)
    %------------------------
    c_markov = {
        ['c(3,1,sync,1)>=c(3,1,sync,2)'] % constant in state 1 higher than in state 2 (1=good regime, 2=bad regime)
        };

    % Restrictions for Sigma(s_t)
    %----------------------------
    s_markov = {
        ['s(3,3,sync,1)<=s(3,3,sync,2)'] % volatility in state 1 lower than in state 2 (1=good regime, 2=bad regime)
        };

    for istate=1:numberOfStates
            
            mystate=int2str(istate);

            % Restrictions on a0(s_t)
            a0_markov=[a0_markov
                {
                ['a0(1,GDPGH,sync,',mystate,')=0'] % No contemporaneous effect of GDP on FF
                ['a0(1,MF,sync,',mystate,')=0'] % No contemporaneous effect of MF on FF
%                 ['a0(2,FF,sync,',mystate,')=0']    % No contemporaneous effect of FF on MF      
                ['a0(2,GDPGH,sync,',mystate,')=0']  % No contemporaneous effect of GDP on MF        
                }];
            
            
            % Restrictions on a1(s_t) in most models
            for ilag=1:nlags
                mylag=int2str(ilag);
                a1_markov=[a1_markov
                    {
                    ['a',mylag,'(1,GDPGH,sync,',mystate,')=0']    % No effect of GDP(t-i) on FF 
                    ['a',mylag,'(2,GDPGH,sync,',mystate,')=0']    % No effect of GDP(t-i) on MF 
%                     ['a',mylag,'(3,GDPGH,sync,',mystate,')=0']    % No effect of GDP(t-i) on GDP 
                    }];
            end   
            if modelspec==12
                for ilag=1:nlags
                    mylag=int2str(ilag);
                    a1_markov=[a1_markov
                        {
                        ['a',mylag,'(3,GDPGH,sync,',mystate,')=0']    % No effect of GDP(t-i) on GDP 
                        }];
                end                    
            elseif modelspec==9
                for ilag=1:nlags
                    mylag=int2str(ilag);
                    a1_markov=[a1_markov
                        {
                        ['a1(3,GDPGH,sync,',mystate,')=0']  % No effect of GDP(t-1) on GDP   
                        }];
                end                    
            else
                a1_markov=[a1_markov
                    {
                    }];
            end


    end

    % across all models, FF and MF do not load on past GDPGH
    if imodel==1 || imodel==13
        lin_restr=[c_markov;a0_const;a1_const;s_markov];
    elseif imodel==2|| imodel==10 || imodel==11  
         % 10 = like model 2 but lags of MF and FF does not enter GDP equation
         % 11 = like model 2 but GDP lag set to zero
        lin_restr=[c_markov;a0_markov;a1_const;s_markov];
    elseif imodel==3
        lin_restr=[c_markov;a0_markov;a1_markov];
    elseif imodel==4 || imodel==12  % like model 4 but GDP lag set to zero
        lin_restr=[c_markov;a0_markov;a1_markov;s_markov];
    elseif imodel==5
        lin_restr=[c_markov;a0_const;a1_markov;s_markov];
    elseif imodel==6
        lin_restr=[c_markov;a0_markov;a1_const];
    elseif imodel==7
        lin_restr=[c_markov;a0_const;a1_markov]; % a0_const since only a0(3) switching and restrictions affect a0(1) and a0(2)
    elseif imodel==8
        lin_restr=[c_markov;a0_const;a1_markov;s_markov]; % a0_const since only a0(3) switching and restrictions affect a0(1) and a0(2)
    elseif imodel==9
        lin_restr=[c_markov;a0_const;a1_const;a1_markov;s_markov]; % a0_const since only a0(3) switching and restrictions affect a0(1) and a0(2)
    end
   
elseif opt.dir_or_it ==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterated Forecast Model %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Restrictions for constant elements of A0
    %-----------------------------------------
    a0_const = {
            'a0(1,GDPG)=0' % No contemporaneous effect of GDP on FF        
            'a0(1,MF)=0'   % No contemporaneous effect of MF on FF       
            'a0(2,GDPG)=0'  % No contemporaneous effect of GDP on MF              
%             'a0(2,FF)=0'  % No contemporaneous effect of FF on MF             
            };
        
    % Restrictions for constant elements of A1
    %-----------------------------------------
   a1_const=cell(0,1);
    if modelspec==10
       for ilag=1:nlags
            mylag=int2str(ilag);
            a1_const=[a1_const
                 {
                ['a',mylag,'(3,FF)=0'] % No effect of FF(t-i) on GDP
                ['a',mylag,'(3,MF)=0'] % No effect of MF(t-i) on GDP
%                 ['a',mylag,'(3,GDPG)=0'] % No effect of GDP(t-i) on GDP
                }];
       end
    elseif modelspec==11 || modelspec==13
       for ilag=1:nlags
            mylag=int2str(ilag);
            a1_const=[a1_const
                 {
                ['a',mylag,'(3,GDPG)=0'] % No effect of GDP(t-i) on GDP
                }];
       end
    elseif modelspec==9
       for ilag=1:nlags
            mylag=int2str(ilag);
            a1_const=[a1_const
                 {
                ['a',mylag,'(3,GDPG)=0'] % No effect of GDP(t-i) on GDP
                }];
       end
    else
         a1_const = {
    %                     'a1(1,GDPG)=0' % No effects of GDP(t-1) on FF
    %                     'a1(2,GDPG)=0' % No effect of GDP(t-1) on MF
    %                     'a1(3,FF)=0' % No effect of FF(t-1) on GDP
    
%                          'a1(3,MF)=0' % No effect of MF(t-1) on GDP
%                         'a2(3,MF)=0' % No effect of MF(t-2) on GDP
%                         'a3(3,MF)=0' % No effect of MF(t-3) on GDP
                     
%                          'a1(3,GDPG)=0' % No effect of GDP(t-1) on GDP
%                           'a2(3,GDPG)=0' % No effect of GDP(t-2) on GDP
%                           'a3(3,GDPG)=0' % No effect of GDP(t-3) on GDP
                         
                };
   end

    
    % Restrictions for switching elements of A0 and A1
    %-------------------------------------------------
    numberOfStates=markov_chains.number_of_states;
    a0_markov=cell(0,1);
    a1_markov=cell(0,1);

    % Restrictions for c(s_t)
    %------------------------
    c_markov = {
        ['c(3,1,sync,1)>=c(3,1,sync,2)'] % constant in state 1 higher than in state 2 (1=good regime, 2=bad regime)
        };

    % Restrictions for Sigma(s_t)
    %----------------------------
    s_markov = {
        ['s(3,3,sync,1)<=s(3,3,sync,2)'] % volatility in state 1 lower than in state 2 (1=good regime, 2=bad regime)
        };

    for istate=1:numberOfStates
            
            mystate=int2str(istate);

            % Restrictions on a0(s_t)
            a0_markov=[a0_markov
                {
                ['a0(1,GDPG,sync,',mystate,')=0'] % No contemporaneous effect of GDP on FF  
                ['a0(2,GDPG,sync,',mystate,')=0']  % No contemporaneous effect of GDP on MF       
                ['a0(1,MF,sync,',mystate,')=0'] % No contemporaneous effect of MF on FF  
                %['a0(2,FF,sync,',mystate,')=0']  % No contemporaneous effect of FF on MF     
                }];

            % Restrictions on a1(s_t) in most models
            if modelspec==12
                for ilag=1:nlags
                    mylag=int2str(ilag);
                    a1_markov=[a1_markov
                        {
                        ['a',mylag,'(3,GDPG,sync,',mystate,')=0']    % No effect of GDP(t-i) on GDP 
                        }];
                end
            elseif modelspec==9
                for ilag=1:nlags
                    mylag=int2str(ilag);
                    a1_markov=[a1_markov
                        {
                        ['a1(3,GDPG,sync,',mystate,')=0']  % No effect of GDP(t-1) on GDP   
                        }];
                end                    
            else
                a1_markov=[a1_markov
                    {
    %                 ['a1(1,GDPG,sync,',mystate,')=0']    % No effect of GDP(t-1) on FF   
    %                 ['a1(2,GDPG,sync,',mystate,')=0']    % No effect of GDP(t-1) on MF 
    %                 ['a1(3,FF,sync,',mystate,')=0']    % No effect of FF(t-1) on GDP   
    %                 ['a1(3,MF,sync,',mystate,')=0']    % No effect of MF(t-1) on GDP 
    %                ['a1(3,GDPG,sync,',mystate,')=0']  % No effect of GDP(t-1) on GDP   
                    }];
            end

    end
       
   
    % across all models, FF and MF do not load on past GDPG
    if imodel==1 || imodel==13
        lin_restr=[c_markov;a0_const;a1_const;s_markov];
    elseif imodel==2 || imodel==10 || imodel==11  
         % 10 = like model 2 but lags of MF and FF does not enter GDP equation
         % 11 = like model 2 but GDP lag set to zero
        lin_restr=[c_markov;a0_markov;a1_const;s_markov];
    elseif imodel==3
        lin_restr=[c_markov;a0_markov;a1_markov];
    elseif imodel==4 || imodel==12  % like model 4 but GDP lag set to zero
        lin_restr=[c_markov;a0_markov;a1_markov;s_markov];
    elseif imodel==5
        lin_restr=[c_markov;c_markov;c_markov;a0_const;a1_markov;s_markov];
    elseif imodel==6
        lin_restr=[c_markov;a0_markov;a1_const];
    elseif imodel==7
        lin_restr=[c_markov;a0_const;a1_markov]; % a0_const since only a0(3) switching and restrictions affect a0(1) and a0(2)
    elseif imodel==8     
        lin_restr=[c_markov;a0_const;a1_markov;s_markov]; % a0_const since only a0(3) switching and restrictions affect a0(1) and a0(2)
    elseif imodel==9   
        lin_restr=[c_markov;a0_const;a1_const;a1_markov;s_markov]; % a0_const since only a0(3) switching and restrictions affect a0(1) and a0(2)
    end

    
end

restrictions=lin_restr;
disp(lin_restr)

%% set priors

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
var_prior.L1=1;
var_prior.L5=1;
var_prior.L6=0; % cointegration (10)
var_prior.L7=0; % unit root  (10)
% var_prior.coefprior=0.9; % impose 0.9 eigenvalue for both

% priors for the sync transition probabilities
%------------------------------------------------
% start, mean, standardDev, distrib;
switch_prior=struct();
if opt.const==1
    if opt.transprob==0
            switch_prior.a12={0.5,0.1,0.3,'beta'};
            switch_prior.a21={0.5,0.1,0.3,'beta'};
    elseif opt.transprob==1
        if opt.normal==1
            switch_prior.a12={0,0,2,'normal'};
            switch_prior.a21={0,0,2,'normal'};
            switch_prior.b12={0,0,2,'normal'};
            switch_prior.c12={0,0,2,'normal'};
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
    end
elseif opt.const==0
    if opt.transprob==0
        errordlg('opt.transprob=0 incompatible with opt.const=0');
    elseif opt.transprob==1
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
end

prior=struct();

prior.var=var_prior;

prior.nonvar=switch_prior;

% Find posterior mode
if isfield(opt,'optimizer')==0
    sv=estimate(sv0,db,{startdb,enddb},prior,restrictions,'fmincon');
else
    if opt.optimizer==1
        % to center estimation around previously estimated mode
        % https://github.com/jmaih/RISE_toolbox/issues/13
        % sv0=load([IS_paramsfolder IS_modelname '.mat']);
        % sv0=sv0.sPMode.sv;
        % sv=estimate(sv0,db,{startdb,enddb},prior,restrictions,'fmincon',false,'estim_start_from_mode', true);
        sv=estimate(sv0,db,{startdb,enddb},prior,restrictions,'fminunc',false);
    elseif opt.optimizer==0 || isfield(opt,'optimizer')==0
        sv=estimate(sv0,db,{startdb,enddb},prior,restrictions,'fmincon');
    end
end

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
% tableTitle = ['Posterior/Mode:/' fcstype ',/' sheetuse ',/h/=/' num2str(opt.hh) '/\\\\/' 'Sample:/' startdb '-' enddb];
% pmode_tex(posterior_mode(sv),texfolder,modelname,tableTitle);

