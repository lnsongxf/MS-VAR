%------------------------------------------------------------------
% Predictive Density - Iterated simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%
%------------------------------------------------------------------
function [out, pden_mat] = fPredDensityIterated(trend_fit,p_reg2_filtered,param_in,FF, MF, nParamDraws,nDraws,tperiods,hh,dates_index,temp)


% Pre-allocate
y_mat_fut_ex_temp = NaN(nParamDraws,nDraws,tperiods); % This stores all realizations for the ITERATED
stemp_hh          = NaN(nParamDraws,nDraws,tperiods);
stemp_t           = NaN(nParamDraws,nDraws,tperiods);


% Loop over parameter draws
for dd=1:nParamDraws
    
    % Loop over repetitions
    for uu=1:nDraws
        
        % Loop over periods
        for tt=1:tperiods
            
            % Macro and financial factor in period t
            f_t = FF(tt);
            m_t = MF(tt);
            
            if strcmp(temp,'drawst')
                % Use filtered probability to determine current state
                pbad_t = p_reg2_filtered(tt);
                udraw = rand(1); % draw a coin that determines the regime
                if udraw < (1-pbad_t)
                    st_sim = 1;
                else
                    st_sim = 2;
                end
                
                % % **** delete FROM HERE ***
            elseif strcmp(temp,'forecastst')
                % Compute transition probabilities for s(t-1) to s(t)
                % Transition probabilities when coefficients have gamma prior
                if tt==1
                    pbad_t = p_reg2_filtered(tt);
                    
                    udraw = rand(1); % draw a coin that determines the regime
                    
                    if udraw < (1-pbad_t)
                        st_sim = 1;
                    else
                        st_sim = 2;
                    end
                    
                else
                    p12 = 1./(1+exp(param_in.a12(dd)-param_in.b12(dd)*FF(tt-1)+param_in.c12(dd)*MF(tt-1)));
                    p21 = 1./(1+exp(param_in.a21(dd)+param_in.b21(dd)*FF(tt-1)-param_in.c21(dd)*MF(tt-1)));
                    
                    % Transition probabilities t|t-1
                    p11 = 1 - p12; % probability of remaining in normal
                    p22 = 1 - p21; % probability of remaining in bad
                    
                end
                
                % Simulate s(t) conditional on m(t) and f(t)
                if tt>1
                    st_sim = simulate_st(st_lag,p11,p22,1);
                end
            end
 
            
            % Simulate factors and transtion probabilities
            % conditional on f_t and m_t from t+1 to t+h
            [f_temp,m_temp,p11_temp,p22_temp] = simulate_factors(f_t, m_t, param_in, hh,dd,'iter');
            
            % Simulate states: s(t+1)...s(t+h)
            s_temp = simulate_st(st_sim,p11_temp,p22_temp,hh);
            
            % Simulate GDP from t+1:t+h
            y_temp = simulate_gdp(s_temp,f_temp,m_temp,param_in,dd,'iter');
            
            % Average future GDP t+1:t+h
            y_mat_fut_ex_temp(dd,uu,tt) = mean(y_temp(1:hh,1)) + trend_fit(tt);
            
            % Collect state in period t + hh
            stemp_hh(dd,uu,tt) = s_temp(hh);
            
            % Collect state in period t
            stemp_t(dd,uu,tt) = st_sim;
            
            % This will be the st_lag for the next tt period iteration
            st_lag = st_sim;
            
        end
    end
end

% Reshape matrices (Row= DrawsxParams, Columns= periods)
y_mat_fut_ex = reshape(y_mat_fut_ex_temp,[nDraws*nParamDraws,tperiods]);
stmat_hh     = reshape(stemp_hh,[nDraws  *nParamDraws,tperiods]);
stmat_t      = reshape(stemp_t,[nDraws*  nParamDraws,tperiods]);

% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
stmat_hh(stmat_hh==2) = 0;
stmat_t(stmat_t==2) = 0;

out.st_h_mean = mean(stmat_hh)';    % This is in t+hh
out.st_t_mean = mean(stmat_t)';     % This is in t
out.dYsim_25  = prctile(y_mat_fut_ex,25)';
out.dYsim_75  = prctile(y_mat_fut_ex,75)';
out.dYsim_10  = prctile(y_mat_fut_ex,10)';
out.dYsim_90  = prctile(y_mat_fut_ex,90)';
out.mean     = mean(y_mat_fut_ex)';

% This is the predictive density matrix (nDraws*nParamDraws)xselected_periods
if ~isempty(dates_index)
    pden_mat = y_mat_fut_ex(:,dates_index);
else
    pden_mat = [];
end
