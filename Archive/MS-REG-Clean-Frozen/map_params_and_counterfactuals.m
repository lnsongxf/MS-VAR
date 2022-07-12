% Script to create counterfactual parameters for predictive density


if opt.counterfactual==1
    % Map parameters to matrices for simulation
    param = scriptParams_FULL_companion(sv,params_in,dd);
    
    % Counterfactual 1: Kill MF in VAR
    fprintf('\n Counterfactual 1: No MF in VAR \n');
    
    % Contemporaneous effect of mf on ff
    param.A0_sync_1(1,2) = 0;
    param.A0_sync_2(1,2) = 0;
    
    % Contemporaneous effect of mf on y
    param.A0_sync_1(3,2) = 0;
    param.A0_sync_1(3,2) = 0;
    
    % Lagged affect of mf on ff and y
    param.A1_sync_1(1,2) = 0;
    param.A1_sync_1(3,2) = 0;
    param.A1_sync_2(1,2) = 0;
    param.A1_sync_2(3,2) = 0;
    
    
elseif opt.counterfactual==2

    % Map parameters to matrices for simulation
    param = scriptParams_FULL_companion(sv,params_in,dd);

    % Counterfactual 2: Kill MF in VAR and in Transition Prob
    fprintf('\n Counterfactual 2: No MF in VAR and Transprob \n');
    
    
    % Contemporaneous effect of mf on y
    param.A0_sync_1(3,2) = 0;
    param.A0_sync_2(3,2) = 0;
    
    % Lagged affect of mf on ff
    param.A1_sync_1(1,2) = 0;
    param.A1_sync_2(1,2) = 0;

    % Lagged affect of mf on y
    param.A1_sync_1(3,2) = 0;
    param.A1_sync_2(3,2) = 0;
    
    % Effect of mf on transition probability
    param.c12=0;
    param.c21=0;
    
    
elseif opt.counterfactual==3

    % Map parameters to matrices for simulation
    param = scriptParams_FULL_companion(sv,params_in,dd);

    % Counterfactual 3: Kill FF in VAR and in Transition Prob
    fprintf('\n Counterfactual 3: No FF in VAR and Transprob \n');
    
%     % Contemporaneous effect of ff on mf
%     param.A0_sync_1(2,1) = 0;
%     param.A0_sync_2(2,1) = 0;
    
    % Contemporaneous effect of ff on y
    param.A0_sync_1(3,1) = 0;
    param.A0_sync_2(3,1) = 0;
    
    % Lagged affect of mf on ff and mf
    param.A1_sync_1(2,1) = 0;
    param.A1_sync_2(2,1) = 0;

    % Lagged affect of mf on ff and y
    param.A1_sync_1(3,1) = 0;
    param.A1_sync_2(3,1) = 0;
    
    % Effect of ff on transition probability
    param.b12=0;
    param.b21=0;
    

        
else
    error('counterfactual not specified!');
end