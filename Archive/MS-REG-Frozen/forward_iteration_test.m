%% Recursive updating of s(t) from simulations of Markov-chain (FORWARD ITERATION)

% Number of draws for simulting predictive density
nDraws = 20000;

% Set a 0.9 threshold to switch from good to bad regime
prob_reg2_new(p_reg2(sd:ed)>0.9) = 1;

% Draw shocks for simulation
rng(123);
shocks_sim = randn(12,length(prob_reg2_new),nDraws); %NOW WE NEED 12 OF IT

% Collect MF, FF and trend
FF = db.FF.data(nlags+1:end);
MF = db.MF.data(nlags+1:end);
TR = db.TREND.data(nlags+1:end);

% DRAW SIMULATION FOR PERIOD-t
for nsim=1:nDraws

    for tt= 1:length(prob_reg2_new)
        %fprintf('\n Procession period %i out of %i ',tt,length(prob_reg2_new));
        % Compute the time varying transition probabilities

        % Get regime in t-12

        % Recall here the indicator st is a indicator of switch to the Bad regime
        % Hence, st=0 (Good regime) st=1 (Bad regime)

%       WE DO NOT NEED TO DRAW THE STATE FOR tt<12, AS WE ALREADY KNOW IT
%         if tt<12
%             st_lag = 1; % initialize in bad regime (Depends on when the sample starts)
% 
%             % Set st
%             st_use = 1;

%         else
            % **** TO DO: Draw the initial state from filtered probability of
            
            % regime-2 ***
            st_lag = prob_reg2_new(tt-11);

            % Get transition probabilities from t-11:tt
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12*(FF(tt-11:tt))-pmode.c12*(MF(tt-11:tt))));
                p21 = 1./(1+exp(pmode.a21-pmode.b21*(FF(tt-11:tt))-pmode.c21*(MF(tt-11:tt))));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(pmode.a12-pmode.b12.*(FF(tt-11:tt))+pmode.c12.*(MF(tt-11:tt))));
                p21 = 1./(1+exp(pmode.a21+pmode.b21.*(FF(tt-11:tt))-pmode.c21.*(MF(tt-11:tt))));
            end

            % Compute the probability of remaining in a given regime
            p11 = ones(12,1) - p12;
            p22 = ones(12,1) - p21;


            % Draw regimes from t+1:t+12
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
%         end


        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        % Create FF and MF from t+1 to t+12
        FF_loop = [FF(tt);NaN(12,1)];
        MF_loop = [MF(tt);NaN(12,1)];
        for gg=1:12
            FF_loop(gg+1) = FF_loop(gg)*coeff;         % ADJUST FROM THE VAR COEFFICIENTS
            MF_loop(gg+1) = MF_loop(gg)*coeff;         % ADJUST FROM THE VAR COEFFICIENTS
        end
        FF_loop = FF_loop(2:end); MF_loop = MF_loop(2:end);
        
        % WE EITHER: ITERATE FORWARD THE 12 PERIODS, SUM THEM UP AND THEN COMBINE, OR
        % RECOMBINE EACH PERIOD WITH THE REGIME DRAW AND THEN SUM THEM UP

        % Good regime in period t+12
        % ITERATE GOOD REGIME FORWARD
        if tt<length(prob_reg2_new)-12 % NOTE: WE DO NOT HAVE A RULE FOR THE TREND; USE THE HISTORY UP TO THE LAST OBSERVATION, WHEN WE KEEP IT FIXED
            TR_temp = TR(tt+1:tt+12);
        else
            TR_temp = ones(12,1)*TR(tt);
        end
        dY_sim_temp = NaN(12,1);
        for rr=1:12
            dY_sim_temp(rr) = (pmode.c_3_1_sync_1 - pmode.a0_3_1_sync_1*FF_loop(rr)...
            - pmode.a0_3_2_sync_1*MF_loop(rr) + pmode.s_3_3_sync_1*shocks_sim(tt,nsim,rr))...
            + TR_temp(rr);
        end
        dY_sim_1(tt,nsim) = sum(dY_sim_temp)/12;

        % Bad regime in period t+12
        % ITERATE BAD REGIME FORWARD
        dY_sim_temp = NaN(12,1);
        for rr=1:12
            dY_sim_temp(rr) = (pmode.c_3_1_sync_2 - pmode.a0_3_1_sync_2*FF_loop(rr)...
            - pmode.a0_3_2_sync_2*MF_loop(rr)+ pmode.s_3_3_sync_2*shocks_sim(tt,nsim,rr))...
            + TR_temp(rr);
        end        
        dY_sim_2(tt,nsim) = sum(dY_sim_temp)/12;

        dY_sim(tt,nsim) = dY_sim_1(tt,nsim)*IND_good + dY_sim_2(tt,nsim)*(1-IND_good);

        St_sim(tt,nsim) = IND_good;

        % Weight regime simulation using filtered probabilities
        dY_sim_weighted(tt,nsim) = dY_sim_1(tt,nsim)*p_reg1_filtered(tt) + dY_sim_2(tt,nsim)*p_reg2_filtered(tt);

    end

    if mod(nsim,1000)==0
        fprintf('\n Draw %i out of %i',nsim,nDraws);
    end
end