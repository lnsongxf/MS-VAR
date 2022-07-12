%------------------------------------------------------------------
% Predictive Density - Iterated simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%
%------------------------------------------------------------------
function [out, pden_mat] = fPredDensityIteratedFullFast(sv,params_in,FF,MF,GDP,trend_fit,opt,temp,pbad)

tperiods = 1;
nlags = sv.nlags;

% Pre-allocate simulation matrices
y_mat_fut_ex_temp  = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_good = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_bad = NaN(opt.nParamDraws,opt.nDraws,tperiods);

% Pre-allocate state matrices
stemp_hh  = NaN(opt.nParamDraws,opt.nDraws,tperiods);
stemp_t   = NaN(opt.nParamDraws,opt.nDraws,tperiods);

tt=1;

% Loop over parameter draws
for dd=1:opt.nParamDraws
    
    
    % Map parameters
    if isfield(opt,'counterfactual')
       
        map_params_and_counterfactuals;
        
        % Recompute reduced form matrices based on counterfactuals
        param.D_sync_1 = param.A0_sync_1\param.C_sync_1;
        param.D_sync_2 = param.A0_sync_2\param.C_sync_2;
        param.B_sync_1 = param.A0_sync_1\param.A1_sync_1;
        param.B_sync_2 = param.A0_sync_2\param.A1_sync_2;
        param.O_sync_1 = param.A0_sync_1\param.SIG_sync_1;
        param.O_sync_2 = param.A0_sync_2\param.SIG_sync_2;

        
    else
        
        % Use estimated parameters
        param = scriptParams_FULL_companion(sv,params_in,dd); % works for any lag and model
        
    end
    
    
    % Loop over repetitions
    parfor uu=1:opt.nDraws
%       for uu=1:opt.nDraws

        % Initialize st_lag.
        st_lag = 1;

        %  for uu=1:opt.nDraws
        % Loop over periods
        if mod(uu,1000)==0
            fprintf('Current draw: %i of %i\n',uu,opt.nDraws);
        end
        
        % opt.date_t is the last period used for estimation
        for tt= opt.date_t
            
            
            % Macro and financial factor in period t
            f_t = FF(tt);
            m_t = MF(tt);
            gdp_t = GDP(tt);
            
            if nlags>1
                f_t_comp = FF((tt-nlags+1):tt);
                m_t_comp = MF((tt-nlags+1):tt);
                gdp_t_comp = GDP((tt-nlags+1):tt);
            end
            
            
            % Get current state
            % Here we have two options:
            % drawst = draw from filtered states every period
            % forecastst = recursively simulate st
            
            if strcmp(temp,'drawst')
                % Use filtered probability to determine current state
                udraw = rand(1); % draw a coin that determines the regime
                
                if udraw < (1-pbad)
                    st_sim = 1;
                else
                    st_sim = 2; % bad regime
                end
                
            elseif strcmp(temp,'forecastst')
                % Compute transition probabilities for s(t-1) to s(t)
                % Transition probabilities when coefficients have gamma prior
                if tt==nlags
                    
                    udraw = rand(1); % draw a coin that determines the regime
                    
                    if udraw < (1-pbad)
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
                
                % Simulate s(t) conditional on m(t-1) and f(t-1)
                if tt>nlags
                    st_sim = simulate_st(st_lag,p11,p22,1);
                end
            end
            
            % GDPG simulation for t+1 to t+h
            % Reduced form model is of the form:
            % y(t+1) = a(s(t+1|t))+b(s(t+1|t))*f(t)+c(s(t+1|t))*m(t)+d(s(t+1|t))*y(t) + sig(s(t+1|t))*eps(t)

            if nlags==1
                [y_temp,y_temp1,y_temp2,s_out] = simulate_IteratedFull(st_sim,f_t,m_t,gdp_t,param,opt);
            else
                [y_temp,y_temp1,y_temp2,s_out] = simulate_IteratedFull_companion(st_sim,f_t_comp,m_t_comp,gdp_t_comp,param,opt);
            end
            
            % Average future GDP t+1:t+h
            y_mat_fut_ex_temp(dd,uu,1)      = mean(y_temp(1:opt.hh,1))  + trend_fit(tt);
            y_mat_fut_ex_temp_good(dd,uu,1) = mean(y_temp1(1:opt.hh,1)) + trend_fit(tt);
            y_mat_fut_ex_temp_bad(dd,uu,1)  = mean(y_temp2(1:opt.hh,1)) + trend_fit(tt);
            
            
            % Collect state in period t + opt.hh
            % stemp_hh(dd,uu,tt) = s_temp(opt.hh);
            s_out(s_out==2) = 0;
            stemp_hh(dd,uu,1) = mean(s_out);
            
            % Collect state in period t
            stemp_t(dd,uu,1) = st_sim;
            
            % This will be the st_lag for the next tt period iteration
            st_lag = st_sim;
            
        end
    end
end

% Reshape matrices (Row= DrawsxParams, Columns= periods)
y_mat_fut_ex       = reshape(y_mat_fut_ex_temp,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_good  = reshape(y_mat_fut_ex_temp_good,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_bad   = reshape(y_mat_fut_ex_temp_bad,[opt.nDraws*opt.nParamDraws,tperiods]);

% stmat_hh     = reshape(stemp_hh,[opt.nDraws  *opt.nParamDraws,tperiods]);
stmat_t        = reshape(stemp_t,[opt.nDraws*opt.nParamDraws,tperiods]);
stmat_hh        = reshape(stemp_hh,[opt.nDraws*opt.nParamDraws,tperiods]);

% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
stmat_t(stmat_t==2) = 0;


% Quantiles for QW CRPS 
qw = opt.qw;
qwvec = zeros(1,length(qw));
for count = 1:length(qw)
    qwvec(1,count) = prctile(y_mat_fut_ex,qw(count)*100)';
end


out.st_h_mean = mean(stmat_hh)';    % This is in t+opt.hh
out.st_t_mean = mean(stmat_t)';       % This is in t
out.dYsim_25  = prctile(y_mat_fut_ex,25)';
out.dYsim_75  = prctile(y_mat_fut_ex,75)';
out.dYsim_10  = prctile(y_mat_fut_ex,10)';
out.dYsim_90  = prctile(y_mat_fut_ex,90)';
out.mean      = mean(y_mat_fut_ex)';
out.qwvec     = qwvec;

% This is the predictive density matrix (opt.nDraws*opt.nParamDraws)xselected_periods
pden_mat.realized = y_mat_fut_ex(:,end);
pden_mat.good     = y_mat_fut_ex_good(:,end);
pden_mat.bad      = y_mat_fut_ex_bad(:,end);

