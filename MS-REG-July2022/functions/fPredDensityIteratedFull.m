%------------------------------------------------------------------
% Predictive Density - Iterated simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%
%------------------------------------------------------------------
function [out, pden_mat] = fPredDensityIteratedFull(sv,params_in,db_full,opt,temp)

varlist = opt.varlist;
tperiods = opt.tperiods;
date_index = opt.date_index;
nvars = opt.nvars;
nlags = opt.nlags;
nParamDraws = opt.nParamDraws;
nDraws = opt.nDraws;

% Find index for financial factor and macro factor
idxf = find(ismember(varlist, 'FF'));
idxm = find(ismember(varlist, 'MF'));

% Collect data and trends for simulation
YY = zeros(tperiods,nvars);
thetrend = zeros(tperiods,nvars);
for jj=1:nvars
    eval(['YY(:,jj) = db_full.' opt.varlist{jj} '.data;']);
    eval(['thetrend(:,jj) = db_full.' 'TRENDH_' opt.varlist{jj} 'H'  '.data;']);
end

% Pre-allocate simulation matrices
y_mat_fut_ex_temp  = NaN(nParamDraws,nDraws,nvars,tperiods); 
y_mat_fut_ex_temp_good = NaN(nParamDraws,nDraws,nvars,tperiods); 
y_mat_fut_ex_temp_bad = NaN(nParamDraws,nDraws,nvars,tperiods);

% Pre-allocate state matrices
stemp_hh  = NaN(nParamDraws,nDraws,tperiods);
stemp_t   = NaN(nParamDraws,nDraws,tperiods);

% Collect filtered probabilities
[~,~,~,f]=filter(sv);
% Filtered probabilities
% p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;
p_reg2_updated  = f.updated_regime_probabilities.regime2.data;

% if compiled, shut down parfor
if isdeployed
    parforArg = 0;
else
    parforArg = Inf;
end


% Loop over parameter draws
for dd=1:nParamDraws
          
    % Use estimated parameters
    param = scriptParams_FULL_companion(sv,params_in,dd); % works for any lag and model
            
    % Filtered probability of regime 2
    pbad = p_reg2_updated;
    
    % Loop over repetitions
    parfor (uu=1:nDraws,parforArg)
%     for uu=1:nDraws

        % Initialize st_lag.
        st_lag = 1;

        % Loop over periods
        if mod(uu,1000)==0
            fprintf('Current draw: %i of %i\n',uu,nDraws);
        end
        
        
        for tt=nlags:tperiods            

            % Get current state
            % Here we have two options:
            % drawst = draw from filtered states every period
            % forecastst = recursively simulate st
            
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
                if tt==nlags
                    
                    udraw = rand(1); % draw a coin that determines the regime
                    
                    if udraw < (1-pbad(tt))
                        st_sim = 1;
                    else
                        st_sim = 2; % bad regime
                    end
                    
                else
                    % Realized path of f and m at t-1
                    f_lag = YY(tt-1,idxf);
                    m_lag = YY(tt-1,idxm);
                    
                    [p12,p21] = fTranstionProb(param,m_lag,f_lag,opt);
                    
                    % Transition probabilities t|t-1
                    p11 = 1 - p12; % probability of remaining in normal
                    p22 = 1 - p21; % probability of remaining in bad
                    
                end
                
                % Simulate s(t) conditional on m(t) and f(t)
                if tt>nlags
                    st_sim = simulate_st(st_lag,p11,p22,1);
                end
            end
            
            % Simulate Model
            if nlags==1
                [y_temp,y_temp1,y_temp2,s_out] = simulate_IteratedFull(st_sim,YY(tt,:)',param,opt);
            else
                YY_comp = YY((tt-nlags+1):tt,:);
               [y_temp,y_temp1,y_temp2,s_out] = simulate_IteratedFull_companion(st_sim,YY_comp,param,opt);
            end
            
            % Collect simulations
            y_mat_fut_ex_temp(dd,uu,:,tt) = mean(y_temp,2)' + thetrend(tt,:);
            y_mat_fut_ex_temp_good(dd,uu,:,tt) = mean(y_temp1,2)' + thetrend(tt,:);
            y_mat_fut_ex_temp_bad(dd,uu,:,tt)  = mean(y_temp2,2)' + thetrend(tt,:);

            
            % Collect state in period t + hh
            % stemp_hh(dd,uu,tt) = s_temp(hh);
            s_out(s_out==2) = 0;
            stemp_hh(dd,uu,tt) = mean(s_out);

            % Collect state in period t
            stemp_t(dd,uu,tt) = st_sim;
            
            % This will be the st_lag for the next tt period iteration
            st_lag = st_sim;
            
        end
    end
end

% Reshape matrices (Row= DrawsxParams, Columns= periods)
y_mat_fut_ex       = reshape(y_mat_fut_ex_temp,[nDraws*nParamDraws,nvars,tperiods]);
y_mat_fut_ex_good  = reshape(y_mat_fut_ex_temp_good,[nDraws*nParamDraws,nvars,tperiods]);
y_mat_fut_ex_bad   = reshape(y_mat_fut_ex_temp_bad,[nDraws*nParamDraws,nvars,tperiods]);

% stmat_hh     = reshape(stemp_hh,[nDraws  *nParamDraws,tperiods]);
stmat_t        = reshape(stemp_t,[nDraws*nParamDraws,tperiods]);
stmat_hh        = reshape(stemp_hh,[nDraws*nParamDraws,tperiods]);

% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
stmat_t(stmat_t==2) = 0;

out.st_h_mean = mean(stmat_hh)';    % This is in t+hh
out.st_t_mean = mean(stmat_t)';       % This is in t
out.dYsim_25  = squeeze(prctile(y_mat_fut_ex,25))';
out.dYsim_75  = squeeze(prctile(y_mat_fut_ex,75))';
out.dYsim_10  = squeeze(prctile(y_mat_fut_ex,10))';
out.dYsim_90  = squeeze(prctile(y_mat_fut_ex,90))';
out.mean      = squeeze(mean(y_mat_fut_ex))';

% This is the predictive density matrix (nDraws*nParamDraws)xselected_periods
if ~isempty(date_index)
    pden_mat.full     = y_mat_fut_ex;
    pden_mat.realized = y_mat_fut_ex(:,:,date_index);
    pden_mat.good     = y_mat_fut_ex_good(:,:,date_index);
    pden_mat.bad      = y_mat_fut_ex_bad(:,:,date_index);
else
    pden_mat = [];
end

