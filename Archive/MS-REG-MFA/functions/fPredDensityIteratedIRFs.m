%------------------------------------------------------------------
% Predictive Density - Iterated simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%
%------------------------------------------------------------------
function [out, pden_mat] = fPredDensityIteratedIRFs(sv,params_in,FF,MF,GDP,trend_fit,opt,temp)
hh = opt.hh;
tperiods = opt.tperiods;
nlags = sv.nlags;

% Pre-allocate simulation matrices
y_mat_fut_ex_temp  = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_good = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_bad = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_shock  = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_good_shock = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_fut_ex_temp_bad_shock = NaN(opt.nParamDraws,opt.nDraws,tperiods);
y_mat_IRF_temp  = NaN(opt.nParamDraws,opt.nDraws,tperiods,hh);
y_mat_IRF_good_temp = NaN(opt.nParamDraws,opt.nDraws,tperiods,hh);
y_mat_IRF_bad = NaN(opt.nParamDraws,opt.nDraws,tperiods,hh);

% Pre-allocate state matrices
stemp_hh  = NaN(opt.nParamDraws,opt.nDraws,tperiods);
stemp_hh_shock  = NaN(opt.nParamDraws,opt.nDraws,tperiods);
stemp_hh_IRF  = NaN(opt.nParamDraws,opt.nDraws,tperiods,hh);
stemp_hh_IRF_shock  = NaN(opt.nParamDraws,opt.nDraws,tperiods,hh);
stemp_t   = NaN(opt.nParamDraws,opt.nDraws,tperiods);

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
    
    % Filtered probability of regime 2
    pbad = p_reg2_updated;
    
    % Loop over repetitions
     parfor (uu=1:opt.nDraws,parforArg)
%       for uu=1:opt.nDraws

        % Initialize st_lag.
        st_lag = 1;

        %  for uu=1:opt.nDraws
        % Loop over periods
        if mod(uu,1000)==0
            fprintf('Current draw: %i of %i\n',uu,opt.nDraws);
        end
        
        for tt=nlags:tperiods
            
            
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
                    [p12,p21] = fTranstionProb(param,MF(tt-1),FF(tt-1),opt);
                    
                    % Transition probabilities t|t-1
                    p11 = 1 - p12; % probability of remaining in normal
                    p22 = 1 - p21; % probability of remaining in bad
                    
                end
                
                % Simulate s(t) conditional on m(t) and f(t)
                if tt>nlags
                    st_sim = simulate_st(st_lag,p11,p22,1);
                end
            end
            
            % Collect GDPG simulation
            if nlags==1
                [y_out,y_out_1,y_out_2,s_out,y_out_shock,y_out_1_shock,y_out_2_shock,s_out_shock] = simulate_IteratedFull_IRF(st_sim,f_t,m_t,gdp_t,param,opt);
            else
                errordlg('TO DO');
%                 [y_temp,y_temp1,y_temp2,s_out] = simulate_IteratedFull_companion(st_sim,f_t_comp,m_t_comp,gdp_t_comp,param,opt);
            end
            
            y_temp = y_out(:,end); y_temp1 = y_out_1(:,end); y_temp2 = y_out_2(:,end);
            y_temp_shock = y_out_shock(:,end); y_temp1_shock = y_out_1_shock(:,end); y_temp2_shock = y_out_2_shock(:,end);
            
            % Average future GDP t+1:t+h
            y_mat_fut_ex_temp(dd,uu,tt) = mean(y_temp(1:opt.hh,1)) + trend_fit(tt);
            y_mat_fut_ex_temp_good(dd,uu,tt) = mean(y_temp1(1:opt.hh,1)) + trend_fit(tt);
            y_mat_fut_ex_temp_bad(dd,uu,tt) = mean(y_temp2(1:opt.hh,1)) + trend_fit(tt);
            y_mat_fut_ex_temp_shock(dd,uu,tt) = mean(y_temp_shock(1:opt.hh,1)) + trend_fit(tt);
            y_mat_fut_ex_temp_good_shock(dd,uu,tt) = mean(y_temp1_shock(1:opt.hh,1)) + trend_fit(tt);
            y_mat_fut_ex_temp_bad_shock(dd,uu,tt) = mean(y_temp2_shock(1:opt.hh,1)) + trend_fit(tt);
            
            % IRFs
%             y_mat_IRF_temp(dd,uu,tt,:) = y_temp_shock(1:opt.hh,1);
            y_mat_IRF_temp(dd,uu,tt,:) = y_temp_shock(1:opt.hh,1) - y_temp(1:opt.hh,1);
            y_mat_IRF_good_temp(dd,uu,tt,:) = y_temp1_shock(1:opt.hh,1) - y_temp1(1:opt.hh,1);
            y_mat_IRF_bad_temp(dd,uu,tt,:) = y_temp2_shock(1:opt.hh,1) - y_temp2(1:opt.hh,1);
            % Collect state in period t + opt.hh
            % stemp_hh(dd,uu,tt) = s_temp(opt.hh);
            s_out(s_out==2) = 0;
            s_out_shock(s_out_shock==2) = 0;

            stemp_hh(dd,uu,tt) = mean(s_out);
            stemp_hh_shock(dd,uu,tt) = mean(s_out_shock);
            stemp_hh_IRF(dd,uu,tt,:) = s_out;
            stemp_hh_IRF_shock(dd,uu,tt,:) = s_out_shock;
%             stemp_hh(dd,uu,tt) = s_out(end);
%             stemp_hh_shock(dd,uu,tt) = s_out_shock(end);
            
            % Collect state in period t
            stemp_t(dd,uu,tt) = st_sim;
            
            % This will be the st_lag for the next tt period iteration
            st_lag = st_sim;
            
        end
    end
end

% Reshape matrices (Row= DrawsxParams, Columns= periods)
y_mat_fut_ex       = reshape(y_mat_fut_ex_temp,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_good  = reshape(y_mat_fut_ex_temp_good,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_bad   = reshape(y_mat_fut_ex_temp_bad,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_shock       = reshape(y_mat_fut_ex_temp_shock,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_good_shock  = reshape(y_mat_fut_ex_temp_good_shock,[opt.nDraws*opt.nParamDraws,tperiods]);
y_mat_fut_ex_bad_shock   = reshape(y_mat_fut_ex_temp_bad_shock,[opt.nDraws*opt.nParamDraws,tperiods]);
% stmat_hh     = reshape(stemp_hh,[opt.nDraws  *opt.nParamDraws,tperiods]);
stmat_t        = reshape(stemp_t,[opt.nDraws*opt.nParamDraws,tperiods]);
stmat_hh        = reshape(stemp_hh,[opt.nDraws*opt.nParamDraws,tperiods]);
stmat_hh_shock        = reshape(stemp_hh_shock,[opt.nDraws*opt.nParamDraws,tperiods]);
out.stmat_hh_IRF        = reshape(stemp_hh_IRF,[opt.nDraws*opt.nParamDraws,tperiods,hh]);
out.stmat_hh_IRF_shock        = reshape(stemp_hh_IRF_shock,[opt.nDraws*opt.nParamDraws,tperiods,hh]);

% Reshape IRFs (Row= DrawsxParams, Columns= periods, add=1:hh)
out.y_mat_IRF       = reshape(y_mat_IRF_temp,[opt.nDraws*opt.nParamDraws,tperiods,hh]);
out.y_mat_IRF_good  = reshape(y_mat_IRF_good_temp,[opt.nDraws*opt.nParamDraws,tperiods,hh]);
out.y_mat_IRF_bad   = reshape(y_mat_IRF_bad_temp,[opt.nDraws*opt.nParamDraws,tperiods,hh]);
% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
stmat_t(stmat_t==2) = 0;

out.st_h_mean = mean(stmat_hh)';    % This is in t+opt.hh
out.st_h_mean_shock = mean(stmat_hh_shock)';    % This is in t+opt.hh
out.st_t_mean = mean(stmat_t)';       % This is in t
out.dYsim_25  = prctile(y_mat_fut_ex,25)';
out.dYsim_75  = prctile(y_mat_fut_ex,75)';
out.dYsim_10  = prctile(y_mat_fut_ex,10)';
out.dYsim_90  = prctile(y_mat_fut_ex,90)';
out.mean_shock      = mean(y_mat_fut_ex_shock)';
out.dYsim_25_shock  = prctile(y_mat_fut_ex_shock,25)';
out.dYsim_75_shock  = prctile(y_mat_fut_ex_shock,75)';
out.dYsim_10_shock  = prctile(y_mat_fut_ex_shock,10)';
out.dYsim_90_shock  = prctile(y_mat_fut_ex_shock,90)';
out.mean_shock      = mean(y_mat_fut_ex_shock)';

% This is the predictive density matrix (opt.nDraws*opt.nParamDraws)xselected_periods
if ~isempty(opt.date_index)
    pden_mat.full     = y_mat_fut_ex;
    pden_mat.realized = y_mat_fut_ex(:,opt.date_index);
    pden_mat.good     = y_mat_fut_ex_good(:,opt.date_index);
    pden_mat.bad      = y_mat_fut_ex_bad(:,opt.date_index);
    pden_mat.full_shock     = y_mat_fut_ex_shock;
    pden_mat.realized_shock = y_mat_fut_ex_shock(:,opt.date_index);
    pden_mat.good_shock     = y_mat_fut_ex_good_shock(:,opt.date_index);
    pden_mat.bad_shock      = y_mat_fut_ex_bad_shock(:,opt.date_index);
else
    pden_mat = [];
end

