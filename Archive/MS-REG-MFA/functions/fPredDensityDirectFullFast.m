%------------------------------------------------------------------
% Predictive Density - Direct simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
% # FAST VERSION: Simulates only the required number of periods
%------------------------------------------------------------------
function [out,pden_mat] = fPredDensityDirectFullFast(sv,params_in,FF,MF,GDP,trend_fit,opt,typest,pbad)

% Number of periods for simulation
tperiods = 1;
nlags    = sv.nlags;

% Pre-allocate simulation matrices
y_mat_fut_ex_temp      = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_good = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_bad  = NaN(opt.nParamDraws,opt.nDraws,tperiods);

% Pre-allocate state matrices
stemp_t   = NaN(opt.nParamDraws,opt.nDraws,tperiods);


tt = 1;

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
    for uu=1:opt.nDraws
        % Loop over periods
        if mod(uu,1000)==0
            fprintf('Current draw: %i of %i\n',uu,opt.nDraws);
        end
               
        
        % In the first hh we can construct s(t) from the filtered probabilities
        
        % Use filtered probability to determine state in t-h+1
        udraw = rand(1); % draw a coin that determines the regime
        if udraw < (1-pbad)
            st_temp = 1;
        else
            st_temp = 2;
        end
        
        
        % This is s(t-opt.hh)
        st_lag = st_temp; %st_sim(tt-opt.hh);
        
        %---------------------
        % Forecast the state
        %---------------------
        
        % Realized path of f and m from t-opt.hh+1 to t
        f_temp = FF(1:opt.hh);
        m_temp = MF(1:opt.hh);
        
        % Compute transition probabilities for s(. |t-opt.hh+1),...,s(t+1|t)
        [p12,p21] = fTranstionProb(param,m_temp,f_temp,opt);
        
        p11 = ones(opt.hh,1) - p12; % probability of remaining in normal
        p22 = ones(opt.hh,1) - p21; % probability of remaining in bad
        
        % This returns s(t-opt.hh+1)...s(t+1)
        st_sim_temp = simulate_st(st_lag,p11,p22,opt.hh);
        
        % This is s(t+1|t)
        st_temp = st_sim_temp(opt.hh);
        
        % Macro and financial factor in period t
        f_t = FF(opt.hh);
        m_t = MF(opt.hh);
        gdp_t = GDP(end);
        trend_t = trend_fit(end);
        
        % Simulate GDP
        [y_temp,y_temp1,y_temp2] = simulate_DirectFull(st_temp,f_t,m_t,gdp_t,param);
        
        % Average future GDP t+1:t+h
        y_mat_fut_ex_temp(dd,uu,tt)     = y_temp  + trend_t;        
        y_mat_fut_ex_temp_good(dd,uu,tt)= y_temp1 + trend_t;
        y_mat_fut_ex_temp_bad(dd,uu,tt) = y_temp2 + trend_t;
        
        % Collect s(t+1|t) state
        stemp_t(dd,uu,tt) = st_temp;
        
    end
end


% Reshape matrices
y_mat_fut_ex       = reshape(y_mat_fut_ex_temp,[opt.nDraws*opt.nParamDraws,tperiods-nlags+1]);
y_mat_fut_ex_good  = reshape(y_mat_fut_ex_temp_good,[opt.nDraws*opt.nParamDraws,tperiods-nlags+1]);
y_mat_fut_ex_bad   = reshape(y_mat_fut_ex_temp_bad,[opt.nDraws*opt.nParamDraws,tperiods-nlags+1]);


stmat_t = reshape(stemp_t,[opt.nDraws*opt.nParamDraws,tperiods-nlags+1]);


% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
stmat_t(stmat_t==2) = 0;


st_t_mean = mean(stmat_t)'; % This is in t

dYsim_25  = prctile(y_mat_fut_ex,25)';
dYsim_75  = prctile(y_mat_fut_ex,75)';
dYsim_10  = prctile(y_mat_fut_ex,10)';
dYsim_90  = prctile(y_mat_fut_ex,90)';
dYmean    = mean(y_mat_fut_ex)';


% Quantiles for QW CRPS 
qw = opt.qw;
qwvec = zeros(1,length(qw));
for count = 1:length(qw)
    qwvec(1,count) = prctile(y_mat_fut_ex,qw(count)*100)';
end


out.st_t_mean = st_t_mean;
out.dYsim_25  = dYsim_25;
out.dYsim_75  = dYsim_75;
out.dYsim_10  = dYsim_10;
out.dYsim_90  = dYsim_90;
out.mean      = dYmean;
out.qwvec     = qwvec;

% This is the predictive density matrix (opt.nDraws*opt.nParamDraws)xselected_periods
pden_mat.realized = y_mat_fut_ex(:,end);
pden_mat.good     = y_mat_fut_ex_good(:,end);
pden_mat.bad      = y_mat_fut_ex_bad(:,end);
