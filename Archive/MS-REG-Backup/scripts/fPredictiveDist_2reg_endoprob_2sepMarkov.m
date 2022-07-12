
function [dY_sim_1,dY_sim_2,dY_sim_3,dY_sim_4, dY_sim, dY_sim_weighted] = fPredictiveDist_2reg_endoprob_2sepMarkov(MF,FF,TR,myparams,pnames,p_reg1_filtered,p_reg2_filtered,p_reg3_filtered,p_reg4_filtered,prob_reg34_new,prob_reg24_new,nDraws,normal,st_lag,st_lag_vol)

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

a12_vol = myparams(strcmp('a12_vol',pnames),1);
b12_vol = myparams(strcmp('b12_vol',pnames),1);
c12_vol = myparams(strcmp('c12_vol',pnames),1);
a21_vol = myparams(strcmp('a21_vol',pnames),1);
b21_vol = myparams(strcmp('b21_vol',pnames),1);
c21_vol = myparams(strcmp('c21_vol',pnames),1);
s_3_3_vol_1 = myparams(strcmp('s_3_3_vol_1',pnames),1);
s_3_3_vol_2 = myparams(strcmp('s_3_3_vol_2',pnames),1);

% Regime 1 GDP equation
c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),1);
a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),1);
a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),1);

% Regime 2 GDP equat1on
c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),1);
a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),1);
a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),1);

% Draw shocks
shocks_sim = randn(length(prob_reg34_new),nDraws);

for nsim=1:nDraws

    for tt= 1:length(prob_reg34_new)
        %fprintf('\n Procession period %i out of %i ',tt,length(prob_reg2_new));
        % Compute the time varying transition probabilities

        % Get regime in t-12

        % Recall here the indicator st is a indicator of switch to the Bad regime
        % Hence, st=0 (Good regime) st=1 (Bad regime)

        if tt<12
%             st_lag = 1; % initialize in bad regime (Depends on when the sample starts)
%             st_lag_vol = 1; % initialize in high vol regime (Depends on when the sample starts)

            % Set st
            st_use = 1;
            st_use_vol = 1;
            

        else
            % **** TO DO: Draw the initial state from filtered probability of
            
            % regime-2 ***
            st_lag = prob_reg34_new(tt-11);
            st_lag_vol =  prob_reg24_new(tt-11);

            % Get transition probabilities from t-11:tt
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(a12-b12*(FF(tt-11:tt))-c12*(MF(tt-11:tt))));
                p21 = 1./(1+exp(a21-b21*(FF(tt-11:tt))-c21*(MF(tt-11:tt))));
                p12_vol = 1./(1+exp(a12_vol-b12_vol*(FF(tt-11:tt))-c12_vol*(MF(tt-11:tt))));
                p21_vol = 1./(1+exp(a21_vol-b21_vol*(FF(tt-11:tt))-c21_vol*(MF(tt-11:tt))));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(a12-b12.*(FF(tt-11:tt))+c12.*(MF(tt-11:tt))));
                p21 = 1./(1+exp(a21+b21.*(FF(tt-11:tt))-c21.*(MF(tt-11:tt))));
                p12_vol = 1./(1+exp(a12_vol-b12_vol.*(FF(tt-11:tt))+c12_vol.*(MF(tt-11:tt))));
                p21_vol = 1./(1+exp(a21_vol+b21_vol.*(FF(tt-11:tt))-c21_vol.*(MF(tt-11:tt))));
            end

            % Compute the probability of remaining in a given regime
            p11 = ones(12,1) - p12;
            p22 = ones(12,1) - p21;
            p11_vol = ones(12,1) - p12_vol;
            p22_vol = ones(12,1) - p21_vol;            


            % Draw regimes from t-12:t 
            for tt2 = 1:12

                % Draw the Markov Chain for period t-12:t
                udraw = rand(1);
                udraw_vol = rand(1);

                if st_lag == 0 && st_lag_vol == 0 % started in high growth / low vol

                    if udraw > p11(tt2)
                        st(tt2,1) = 1; % switch from good to bad
                    else
                        st(tt2,1) = 0; % don't switch and remain in good
                    end
                    if udraw_vol > p11_vol(tt2)
                        st_vol(tt2,1) = 1;  % switch from low to high vol
                    else
                        st_vol(tt2,1) = 0;  % don't switch and remain in low vol
                    end
               elseif st_lag == 0 && st_lag_vol == 1 % started in high growth / low vol
                    if udraw > p11(tt2)
                        st(tt2,1) = 1; % switch from good to bad
                    else
                        st(tt2,1) = 0; % don't switch and remain in good
                    end
                    if udraw_vol > p22_vol(tt2)
                        st_vol(tt2,1) = 0;  % switch from high to low vol
                    else
                        st_vol(tt2,1) = 1;  % don't switch and remain in high vol
                    end   
               elseif st_lag == 1 && st_lag_vol == 0 % started in low growth / low vol
                    if udraw > p22(tt2)
                        st(tt2,1) = 0; % switch from bad to good
                    else
                        st(tt2,1) = 1; % don't switch and remain in bad
                    end
                    if udraw_vol > p11_vol(tt2)
                        st_vol(tt2,1) = 1;  % switch from low to high vol
                    else
                        st_vol(tt2,1) = 0;  % don't switch and remain in low vol
                    end
               elseif st_lag == 1 && st_lag_vol == 1 % started in low growth / high vol
                    if udraw > p22(tt2)
                        st(tt2,1) = 0; % switch from bad to good
                    else
                        st(tt2,1) = 1; % don't switch and remain in bad
                    end
                    if udraw_vol > p22_vol(tt2)
                        st_vol(tt2,1) = 0;  % switch from high to low vol
                    else
                        st_vol(tt2,1) = 1;  % don't switch and remain in high vol
                    end
               end

                st_lag = st(tt2,1);
                st_lag_vol = st_vol(tt2,1);
            end

            % st = 0 (Good regime), st = 1 (Bad regime)
            st_use = st(tt2);
            % st = 0 (low vol), st = 1 (high vol)
            st_use_vol = st_vol(tt2);
        end


        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        IND_good_vol = 1 - st_use_vol;
        
        IND_hg_lv = IND_good*IND_good_vol;
        IND_hg_hv = IND_good*st_use_vol;
        IND_lg_lv = st_use*IND_good_vol;
        IND_lg_hv = st_use*st_use_vol;


        % Good regime in period t / low vol
        dY_sim_1(tt,nsim) = (c_3_1_sync_1 - a0_3_1_sync_1*FF(tt)...
            - a0_3_2_sync_1*MF(tt) + s_3_3_vol_1*shocks_sim(tt,nsim))...
            + TR(tt);

        % Good regime in period t / high vol
        dY_sim_2(tt,nsim) = (c_3_1_sync_1 - a0_3_1_sync_1*FF(tt)...
            - a0_3_2_sync_1*MF(tt)+ s_3_3_vol_2*shocks_sim(tt,nsim))...
            + TR(tt);
        
         % Bad regime in period t / low vol       
        dY_sim_3(tt,nsim) = (c_3_1_sync_2 - a0_3_1_sync_2*FF(tt)...
            - a0_3_2_sync_2*MF(tt) + s_3_3_vol_1*shocks_sim(tt,nsim))...
            + TR(tt);

        % Bad regime in period t / high vol
        dY_sim_4(tt,nsim) = (c_3_1_sync_2 - a0_3_1_sync_2*FF(tt)...
            - a0_3_2_sync_2*MF(tt)+ s_3_3_vol_2*shocks_sim(tt,nsim))...
            + TR(tt);
        
        dY_sim(tt,nsim) = dY_sim_1(tt,nsim)*IND_hg_lv + dY_sim_2(tt,nsim)*IND_hg_hv + dY_sim_3(tt,nsim)*IND_lg_lv + dY_sim_4(tt,nsim)*IND_lg_hv;

        St_sim1(tt,nsim) = IND_hg_lv;
        St_sim2(tt,nsim) = IND_hg_hv;
        St_sim3(tt,nsim) = IND_lg_lv;
        St_sim4(tt,nsim) = IND_lg_hv;
        
        % Weight regime simulation using filtered probabilities
        dY_sim_weighted(tt,nsim) = dY_sim_1(tt,nsim)*p_reg1_filtered(tt) + dY_sim_2(tt,nsim)*p_reg2_filtered(tt) + dY_sim_3(tt,nsim)*p_reg3_filtered(tt) + dY_sim_4(tt,nsim)*p_reg4_filtered(tt);

    end

end
end

