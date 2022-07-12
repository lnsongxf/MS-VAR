
function [dY_sim_1,dY_sim_2, dY_sim, dY_sim_weighted] = fPredictiveDist_2reg_endoprob(MF,FF,TR,myparams,pnames, p_reg1_filtered, p_reg2_filtered,prob_reg2_new, nDraws,normal,const)

%===================
% MAP PARAMETERS
%===================

% Transition prob parameters
a12 = myparams(strcmp('a12',pnames),1);
b12 = myparams(strcmp('b12',pnames),1);
c12 = myparams(strcmp('c12',pnames),1);
a21 = myparams(strcmp('a21',pnames),1);
b21 = myparams(strcmp('b21',pnames),1);
c21 = myparams(strcmp('c21',pnames),1);

% Regime 1 GDP equation
c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),1);
a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),1);
a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),1);
s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),1);

% Regime 2 GDP equat1on
c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),1);
a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),1);
a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),1);
s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),1);

% Draw shocks
shocks_sim = randn(length(prob_reg2_new),nDraws);


% DRAW SIMULATION FOR PERIOD-t
for nsim=1:nDraws
    
    for tt= 1:length(prob_reg2_new)
        %fprintf('\n Procession period %i out of %i ',tt,length(prob_reg2_new));
        % Compute the time varying transition probabilities
        
        
        % Get regime in t-12
        
        % Recall here the indicator st is a indicator of switch to the Bad regime
        % Hence, st=0 (Good regime) st=1 (Bad regime)
        
        if tt<12
            st_lag = 1; % initialize in bad regime (Depends on when the sample starts)
            
            % Set st
            st_use = 1;
            
        else
            % **** Draw the initial state from filtered probability of
            % regime-2 ***
            st_lag = prob_reg2_new(tt-11);
            
            % Get transition probabilities from t-11:tt
            if const==1
                if normal
                    % Transition probabilities when coefficients have normal prior
                    p12 = 1./(1+exp(a12-b12*(FF(tt-11:tt))-c12*(MF(tt-11:tt))));
                    p21 = 1./(1+exp(a21-b21*(FF(tt-11:tt))-c21*(MF(tt-11:tt))));
                else
                    % Transition probabilities when coefficients have gamma prior
                    p12 = 1./(1+exp(a12-b12.*(FF(tt-11:tt))+c12.*(MF(tt-11:tt))));
                    p21 = 1./(1+exp(a21+b21.*(FF(tt-11:tt))-c21.*(MF(tt-11:tt))));
                end
            elseif const==0
                if normal
                    % Transition probabilities when coefficients have normal prior
                    p12 = 1./(1+exp(1-b12*(FF(tt-11:tt))-c12*(MF(tt-11:tt))));
                    p21 = 1./(1+exp(1-b21*(FF(tt-11:tt))-c21*(MF(tt-11:tt))));
                else
                    % Transition probabilities when coefficients have gamma prior
                    p12 = 1./(1+exp(1-b12.*(FF(tt-11:tt))+c12.*(MF(tt-11:tt))));
                    p21 = 1./(1+exp(1+b21.*(FF(tt-11:tt))-c21.*(MF(tt-11:tt))));
                end
            end

            
            % Compute the probability of remaining in a given regime
            p11 = ones(12,1) - p12;
            p22 = ones(12,1) - p21;
            
            
            % Draw sunspot from t-12:t
            for tt2 = 1:12
                
                % Draw the Markov Chain for period t-12:t
                udraw = rand(1);
                
                if st_lag == 0 % started in good regime
                    if udraw > p11(tt2)
                        st(tt2,1) = 1; % switch from good to bad
                    else
                        st(tt2,1) = 0; % don't switch and remain in good
                    end
                    
                else % start in bad regime
                    
                    if udraw > p22(tt2)
                        st(tt2,1) = 0; % switch from bad to good
                    else
                        st(tt2,1) = 1; % don't switch and remain in bad
                    end
                    
                end
                
                st_lag = st(tt2,1);
            end
            
            % st = 0 (Good regime), st = 1 (Bad regime)
            st_use = st(tt2);
        end
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        % Good regime in period t
        dY_sim_1(tt,nsim) = (c_3_1_sync_1 - a0_3_1_sync_1*FF(tt)...
            - a0_3_2_sync_1*MF(tt) + s_3_3_sync_1*shocks_sim(tt,nsim))...
            + TR(tt);
        
        % Bad regime in period t
        dY_sim_2(tt,nsim) = (c_3_1_sync_2 - a0_3_1_sync_2*FF(tt)...
            - a0_3_2_sync_2*MF(tt)+ s_3_3_sync_2*shocks_sim(tt,nsim))...
            + TR(tt);
        
        dY_sim(tt,nsim) = dY_sim_1(tt,nsim)*IND_good + dY_sim_2(tt,nsim)*(1-IND_good);
        
        St_sim(tt,nsim) = IND_good;
        
        % Weight regime simulation using filtered probabilities
        dY_sim_weighted(tt,nsim) = dY_sim_1(tt,nsim)*p_reg1_filtered(tt) + dY_sim_2(tt,nsim)*p_reg2_filtered(tt);
        
    end
    
    %if mod(nsim,1000)==0
    %    fprintf('\n Simulation %i out of %i',nsim,nDraws);
    %end
end

