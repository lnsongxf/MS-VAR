
function [dY_sim_1,dY_sim_2,dY_sim_3,dY_sim_4, dY_sim] = fPredictiveDist_2reg_endoprob_2sepMarkov_POINT(MF,FF,TR,myparams,pnames,p_reg1_filtered,p_reg2_filtered,p_reg3_filtered,p_reg4_filtered,prob_reg34_new,prob_reg24_new,nDraws,normal,st_lag,st_lag_vol,i)

%===================
% MAP PARAMETERS
%===================

% Transition prob parameters
a12 = myparams(strcmp('a12',pnames),i);
b12 = myparams(strcmp('b12',pnames),i);
c12 = myparams(strcmp('c12',pnames),i);
a21 = myparams(strcmp('a21',pnames),i);
b21 = myparams(strcmp('b21',pnames),i);
c21 = myparams(strcmp('c21',pnames),i);

a12_vol = myparams(strcmp('a12_vol',pnames),i);
b12_vol = myparams(strcmp('b12_vol',pnames),i);
c12_vol = myparams(strcmp('c12_vol',pnames),i);
a21_vol = myparams(strcmp('a21_vol',pnames),i);
b21_vol = myparams(strcmp('b21_vol',pnames),i);
c21_vol = myparams(strcmp('c21_vol',pnames),i);
s_3_3_vol_1 = myparams(strcmp('s_3_3_vol_1',pnames),i);
s_3_3_vol_2 = myparams(strcmp('s_3_3_vol_2',pnames),i);

% Regime 1 GDP equation
c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);

% Regime 2 GDP equat1on
c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);

% Draw shocks
shocks_sim = randn(1,nDraws);

for nsim=1:nDraws


           % **** TO DO: Draw the initial state from filtered probability of
            


            % Get transition probabilities from t-11:tt
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(a12-b12*(FF)-c12*(MF)));
                p21 = 1./(1+exp(a21-b21*(FF)-c21*(MF)));
                p12_vol = 1./(1+exp(a12_vol-b12_vol*(FF)-c12_vol*(MF)));
                p21_vol = 1./(1+exp(a21_vol-b21_vol*(FF)-c21_vol*(MF)));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(a12-b12.*(FF)+c12.*(MF)));
                p21 = 1./(1+exp(a21+b21.*(FF)-c21.*(MF)));
                p12_vol = 1./(1+exp(a12_vol-b12_vol.*(FF)+c12_vol.*(MF)));
                p21_vol = 1./(1+exp(a21_vol+b21_vol.*(FF)-c21_vol.*(MF)));
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
        


        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        IND_good_vol = 1 - st_use_vol;
        
        IND_hg_lv = IND_good*IND_good_vol;
        IND_hg_hv = IND_good*st_use_vol;
        IND_lg_lv = st_use*IND_good_vol;
        IND_lg_hv = st_use*st_use_vol;


        % Good regime in period t / low vol
        dY_sim_1(1,nsim) = (c_3_1_sync_1 - a0_3_1_sync_1*FF(tt2)...
            - a0_3_2_sync_1*MF(tt2) + s_3_3_vol_1*shocks_sim(1,nsim))...
            + TR(tt2);

        % Good regime in period t / high vol
        dY_sim_2(1,nsim) = (c_3_1_sync_1 - a0_3_1_sync_1*FF(tt2)...
            - a0_3_2_sync_1*MF(tt2)+ s_3_3_vol_2*shocks_sim(1,nsim))...
            + TR(tt2);
        
         % Bad regime in period t / low vol       
        dY_sim_3(1,nsim) = (c_3_1_sync_2 - a0_3_1_sync_2*FF(tt2)...
            - a0_3_2_sync_2*MF(tt2) + s_3_3_vol_1*shocks_sim(1,nsim))...
            + TR(tt2);

        % Bad regime in period t / high vol
        dY_sim_4(1,nsim) = (c_3_1_sync_2 - a0_3_1_sync_2*FF(tt2)...
            - a0_3_2_sync_2*MF(tt2)+ s_3_3_vol_2*shocks_sim(1,nsim))...
            + TR(tt2);
        
        dY_sim(1,nsim) = dY_sim_1(1,nsim)*IND_hg_lv + dY_sim_2(1,nsim)*IND_hg_hv + dY_sim_3(1,nsim)*IND_lg_lv + dY_sim_4(1,nsim)*IND_lg_hv;

        St_sim1(1,nsim) = IND_hg_lv;
        St_sim2(1,nsim) = IND_hg_hv;
        St_sim3(1,nsim) = IND_lg_lv;
        St_sim4(1,nsim) = IND_lg_hv;
end
        


end