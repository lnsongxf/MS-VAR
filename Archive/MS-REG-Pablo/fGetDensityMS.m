function out = fGetDensityMS(target_date,data,a2tilde_to_a,results,pnames,nDraws,nParamDraws)

inputformat = 'yyyy-MMM';

FF_full = data.FF_full;
MF_full = data.MF_full;
TR_full = data.TR_full;
dates      = data.dates;
dates_full = data.dates_full;
prob_reg2_new = data.prob_reg2_new;
const   = data.const;
normal  = data.normal;

reference_date = datestr(datenum((datetime(target_date,'InputFormat',inputformat))-calmonths(12)));
first_month    = datestr(datenum((datetime(reference_date))+calmonths(1)));
last_month     = datestr(datenum((datetime(reference_date))+calmonths(12)));


% Last period where filtered data is available
t_start = find(datenum(reference_date)==dates);

% Dates for MF and FF
sd_density = find(datenum(first_month)==dates_full);
ed_density = find(datenum(last_month)==dates_full);

% Read
FF_DENSITY = FF_full(sd_density:ed_density);
MF_DENSITY = MF_full(sd_density:ed_density);
TR_DENSITY = TR_full(sd_density:ed_density);

% Allocate matrices
dY_erg_1_matuse = NaN(1,nDraws,nParamDraws);
dY_erg_2_matuse = NaN(1,nDraws,nParamDraws);
dY_erg_matuse   = NaN(1,nDraws,nParamDraws);

% Start waitbar
updateWaitbar = waitbarParfor(nParamDraws, ['Simulating Distribution for ' target_date ' ...']);

% % Collect parameters
 params=[results.pop.x];
 
 myparams=a2tilde_to_a(params);

tic
for i=1:nParamDraws
    %sol = solve(sv,myparams(:,i));
    % sv.estim_
   
    % Compute predictive density for each period (t, nSims,nParamDraws)
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
    
    % Regime 1 GDP equation
    c_3_1_sync_1  = myparams(strcmp('c_3_1_sync_1',pnames),i);
    a0_3_1_sync_1 = myparams(strcmp('a0_3_1_sync_1',pnames),i);
    a0_3_2_sync_1 = myparams(strcmp('a0_3_2_sync_1',pnames),i);
    s_3_3_sync_1  = myparams(strcmp('s_3_3_sync_1',pnames),i);

    % Regime 2 GDP equat1on
    c_3_1_sync_2  = myparams(strcmp('c_3_1_sync_2',pnames),i);
    a0_3_1_sync_2 = myparams(strcmp('a0_3_1_sync_2',pnames),i);
    a0_3_2_sync_2 = myparams(strcmp('a0_3_2_sync_2',pnames),i);
    s_3_3_sync_2  = myparams(strcmp('s_3_3_sync_2',pnames),i);
    
    % Draw shocks for October 2020
    shocks_sim = randn(1,nDraws);
    
    
    % DRAW SIMULATION FOR PERIOD-t
    for nsim=1:nDraws
        
        % **** Draw the initial state in October 2019
        
        % regime-2 ***
        st_lag = prob_reg2_new(t_start);
        
        % Get transition probabilities from Nov-2019 to Oct-2020
        if const==1
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(a12-b12*(FF_DENSITY)-c12*(MF_DENSITY)));
                p21 = 1./(1+exp(a21-b21*(FF_DENSITY)-c21*(MF_DENSITY)));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(a12-b12.*(FF_DENSITY)+c12.*(MF_DENSITY)));
                p21 = 1./(1+exp(a21+b21.*(FF_DENSITY)-c21.*(MF_DENSITY)));
            end
        elseif const==0
            if normal
                % Transition probabilities when coefficients have normal prior
                p12 = 1./(1+exp(-b12*(FF_DENSITY)-c12*(MF_DENSITY)));
                p21 = 1./(1+exp(-b21*(FF_DENSITY)-c21*(MF_DENSITY)));
            else
                % Transition probabilities when coefficients have gamma prior
                p12 = 1./(1+exp(-b12.*(FF_DENSITY)+c12.*(MF_DENSITY)));
                p21 = 1./(1+exp(+b21.*(FF_DENSITY)-c21.*(MF_DENSITY)));
            end
        end
        
        % Compute the probability of remaining in a given regime
        p11 = ones(12,1) - p12;
        p22 = ones(12,1) - p21;
        
        
        % Draw Markov Chain from Nov-2019 to Oct-2020
        for tt2 = 1:12
            
            % Draw the Markov Chain for period t-12:t
            udraw = rand(1);
            
            if st_lag == 0 % started in good regime
                if udraw > p11(tt2)
                    st = 1; % switch from good to bad
                else
                    st = 0; % don't switch and remain in good
                end
                
            else % start in bad regime
                
                if udraw > p22(tt2)
                    st = 0; % switch from bad to good
                else
                    st = 1; % don't switch and remain in bad
                end
                
            end
            
            st_lag = st;
            
        end
        
        % st = 0 (Good regime), st = 1 (Bad regime)
        st_use = st;
        
        
        % Regime indicator in period t (IND_good = 1 == GOOD REGIME)
        IND_good = 1 - st_use;
        
        
        % Good regime in period t
        dY_erg_1_matuse(1,nsim,i) = (c_3_1_sync_1 - a0_3_1_sync_1*FF_DENSITY(tt2)...
            - a0_3_2_sync_1*MF_DENSITY(tt2) + s_3_3_sync_1*shocks_sim(1,nsim))...
            + TR_DENSITY(tt2);
        
        % Bad regime in period t
        dY_erg_2_matuse(1,nsim,i) = (c_3_1_sync_2 - a0_3_1_sync_2*FF_DENSITY(tt2)...
            - a0_3_2_sync_2*MF_DENSITY(tt2)+ s_3_3_sync_2*shocks_sim(1,nsim))...
            + TR_DENSITY(tt2);
        
        dY_erg_matuse(1,nsim,i) = dY_erg_1_matuse(1,nsim,i)*IND_good + dY_erg_2_matuse(1,nsim,i)*(1-IND_good);
        
        St_sim_temp(1,nsim) = IND_good;
        
    end
    
     St_sim_use(:,i) = mean(St_sim_temp);
    

    updateWaitbar(); %#ok<PFBNS>


end
time_erg = toc;

fprintf('\n Total time posterior simulations =  %4.4f (min) \n',time_erg/60);


% Average out simulations -> returns a TxmParamDraws matrix
out.y_erg_bar   = (reshape(dY_erg_matuse,1,nDraws*nParamDraws));
out.y_erg_bar_1 = (reshape(dY_erg_1_matuse,1,nDraws*nParamDraws));
out.y_erg_bar_2 = (reshape(dY_erg_2_matuse,1,nDraws*nParamDraws));

% Ergodic probability of regime 1;
out.p_reg1_sim =  mean(St_sim_use);
out.p_reg2_sim =  1-out.p_reg1_sim;
