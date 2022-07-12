%------------------------------------------------------------------
% Predictive Density - Direct simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%
%------------------------------------------------------------------
function [out, pden_mat] = fPredDensityDirectSimple(trend_fit,p_reg2_filtered,param_in,FF, MF, opt)

tperiods = opt.tperiods;

% Pre-allocate output
y_mat_fut_direct_temp = NaN(opt.nParamDraws,opt.nDraws,tperiods); % This stores all realizations for the direct
st_direct_temp_t      = NaN(opt.nParamDraws,opt.nDraws,tperiods);
% Loop over parameter draws
for dd=1:opt.nParamDraws
    
    % Counter for rows storing parameter*repetition
    %             waitbar(dd/opt.nParamDraws,wb,'Calculating paths for the direct version')
    
    % Loop over repetitions
    parfor uu=1:opt.nDraws
        %st_sim = 1; % initialization of the simulated states
        st_sim = NaN(tperiods,1);
        
        % Loop over periods
        for tt=1:tperiods
            
            % Simulate s(t) conditional on realized path of m and f
            % from t-12 to t-1.
            if tt<opt.hh+1
                % In the first opt.hh we can construct s(t) from the filtered probabilities
                
                % This is the filtered probability for s(t)
                pbad_t = p_reg2_filtered(tt);
                
                % Use filtered probability to determine current state
                udraw = rand(1); % draw a coin that determines the regime
                if udraw < (1-pbad_t)
                    st_temp = 1;
                else
                    st_temp = 2;
                end
            else
                
                % s(t-opt.hh)
                st_lag = st_sim(tt-opt.hh);
                
                % Realized path of f and m from t-opt.hh,....,t-1
                f_temp = FF(tt-opt.hh:tt-1);
                m_temp = MF(tt-opt.hh:tt-1);
                
                % Compute transition probabilities for s(t-11),...,s(t)
                [p12,p21] = fTranstionProb(param_in(dd),f_temp,m_temp,opt);
                
                p11 = ones(opt.hh,1) - p12; % probability of remaining in normal
                p22 = ones(opt.hh,1) - p21; % probability of remaining in bad
                
                % This returns s(t-11)...s(t))
                st_sim_temp = simulate_st(st_lag,p11,p22,opt.hh);
                
                % This is s(t)
                st_temp = st_sim_temp(opt.hh);
                
            end
            
            st_sim(tt) = st_temp;
            
            % Macro and financial factor in period t
            f_t = FF(tt);
            m_t = MF(tt);
            
            % Simulate GDP
            y_temp = simulate_gdp(st_temp,f_t,m_t,param_in,dd,'direct')+trend_fit(tt);
            
            % Average future GDP t+1:t+h
            y_mat_fut_direct_temp(dd,uu,tt) = y_temp(1);
            
            % Collect t state
            st_direct_temp_t(dd,uu,tt)      = st_temp;
            
        end
    end
end

% Reshape matrices
y_mat_fut_direct = reshape(y_mat_fut_direct_temp,[opt.nDraws*opt.nParamDraws,tperiods]);
st_direct_t = reshape(st_direct_temp_t,[opt.nDraws*opt.nParamDraws,tperiods]);

% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
st_direct_t(st_direct_t==2) = 0;

out.st_t_mean= mean(st_direct_t)'; % This is in t
out.dYsim_25 = prctile(y_mat_fut_direct,25)';
out.dYsim_75 = prctile(y_mat_fut_direct,75)';
out.dYsim_10 = prctile(y_mat_fut_direct,10)';
out.dYsim_90 = prctile(y_mat_fut_direct,90)';
out.mean     = mean(y_mat_fut_direct)';

% This is the predictive density matrix (opt.nDraws*opt.nParamDraws)xtperiods
if ~isempty(opt.date_index)
    pden_mat = y_mat_fut_direct(:,opt.date_index);
else
    pden_mat = [];
end

% *** Need to add regime specific predictive densities ***