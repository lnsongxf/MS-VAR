%------------------------------------------------------------------
% Predictive Density - Iterated simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%
%------------------------------------------------------------------
function [out, pden_mat] = fPredDensityIteratedFull_dgp(param_in,FF,MF,GDP,trend_fit,p_reg_good,opt,temp,modelspec)

tperiods = opt.tperiods;

% Pre-allocate simulation matrices
y_mat_fut_ex_temp  = NaN(opt.nParamDraws,opt.nDraws,tperiods); 
y_mat_fut_ex_temp_good = NaN(opt.nParamDraws,opt.nDraws,tperiods); 
y_mat_fut_ex_temp_bad = NaN(opt.nParamDraws,opt.nDraws,tperiods); 

% Pre-allocate state matrices
% stemp_hh  = NaN(opt.nParamDraws,opt.nDraws,tperiods);
stemp_t   = NaN(opt.nParamDraws,opt.nDraws,tperiods);

% Collect filtered probabilities
pbad = 1-p_reg_good;            
            
% Loop over parameter draws
for dd=1:opt.nParamDraws

    % Check regime indicator and collect VAR matrices
    param = scriptParams_FULL(param_in,dd,modelspec);

    
    % Loop over repetitions
    parfor uu=1:opt.nDraws
        % Loop over periods
        if mod(uu,1000)==0
        fprintf('Current draw: %i of %i\n',uu,opt.nDraws);
        end
        
        for tt=1:tperiods
            
  
            % Macro and financial factor in period t
            f_t = FF(tt);
            m_t = MF(tt);
            gdp_t = GDP(tt);
            

            % Get current state
            % Here we have two options:
            % drawst = draw from filtered states every period
            % forecastst = recurively simulate st

            if strcmp(temp,'drawst')
                % Use filtered probability to determine current state
                udraw = rand(1); % draw a coin that determines the regime
                if udraw < (1-pbad(tt))
                    st_sim = 1;
                else
                    st_sim = 2; % bad regime
                end
                
            elseif strcmp(temp,'forecastst')
                % Compute transition probabilities for s(t-1) to s(t)
                % Transition probabilities when coefficients have gamma prior
                if tt==1
                    
                    udraw = rand(1); % draw a coin that determines the regime
                    
                    if udraw < (1-pbad(tt))
                        st_sim = 1;
                    else
                        st_sim = 2; % bad regime
                    end
                    
                else
                    [p12,p21] = fTranstionProb(param,MF(tt-1),FF(tt-1),opt);

                    % Transition probabilities t|t-1
                    p11 = 1 - p12; % probability of remaining in normal
                    p22 = 1 - p21; % probability of remaining in bad
                    
                end
                
                % Simulate s(t) conditional on m(t) and f(t)
                if tt>1
                    st_sim = simulate_st(st_lag,p11,p22,1);
                end
            end
            
            % Collect GDPG simulation
            [y_temp,y_temp1,y_temp2] = simulate_IteratedFull(st_sim,f_t,m_t,gdp_t,param,opt);    
            
            % Average future GDP t+1:t+h
            y_mat_fut_ex_temp(dd,uu,tt) = mean(y_temp(1:opt.hh,1)) + trend_fit(tt);
            y_mat_fut_ex_temp_good(dd,uu,tt) = mean(y_temp1(1:opt.hh,1)) + trend_fit(tt);
            y_mat_fut_ex_temp_bad(dd,uu,tt) = mean(y_temp2(1:opt.hh,1)) + trend_fit(tt);

            % % Collect state in period t + opt.hh
            % stemp_hh(dd,uu,tt) = s_temp(opt.hh);
            
            % % Collect state in period t
            stemp_t(dd,uu,tt) = st_sim;
            
            % % This will be the st_lag for the next tt period iteration 
            
            st_lag = st_sim;
            
        end
    end
end

% Reshape matrices (Row= DrawsxParams, Columns= periods)
y_mat_fut_ex  = reshape(y_mat_fut_ex_temp,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_good  = reshape(y_mat_fut_ex_temp_good,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_bad  = reshape(y_mat_fut_ex_temp_bad,[opt.nDraws*opt.nParamDraws,tperiods]);


% stmat_hh     = reshape(stemp_hh,[opt.nDraws  *opt.nParamDraws,tperiods]);
stmat_t      = reshape(stemp_t,[opt.nDraws*  opt.nParamDraws,tperiods]);

% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
% stmat_hh(stmat_hh==2) = 0;
stmat_t(stmat_t==2) = 0;

% out.st_h_mean = mean(stmat_hh)';    % This is in t+opt.hh
out.st_t_mean = mean(stmat_t)';       % This is in t
out.dYsim_25  = prctile(y_mat_fut_ex,25)';
out.dYsim_75  = prctile(y_mat_fut_ex,75)';
out.dYsim_10  = prctile(y_mat_fut_ex,10)';
out.dYsim_90  = prctile(y_mat_fut_ex,90)';
out.mean     = mean(y_mat_fut_ex)';

% This is the predictive density matrix (opt.nDraws*opt.nParamDraws)xselected_periods
if ~isempty(opt.date_index)
    pden_mat.realized = y_mat_fut_ex(:,opt.date_index);
    pden_mat.good     = y_mat_fut_ex_good(:,opt.date_index);
    pden_mat.bad      = y_mat_fut_ex_bad(:,opt.date_index);
else
    pden_mat = [];
end
