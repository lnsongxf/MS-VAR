%------------------------------------------------------------------
% Predictive Density - Direct simulation
% Generates the Predictive Density for a given parameter vector
% Parameter vector is a structure in which each field corresponds
% to a parameter draw
%
%
%------------------------------------------------------------------
function [out,pden_mat] = fPredDensityDirectFull(sv,params_in,YY,trend_fit,opt,temp)

varlist = opt.varlist;
gname = opt.gname;
tperiods = opt.tperiods;
date_index = opt.date_index;
nvars = opt.nvars;
nlags = opt.nlags;
hh = opt.hh;
nParamDraws = opt.nParamDraws;
nDraws = opt.nDraws;

% Find index for financial factor and macro factor
idxf = find(ismember(varlist, 'FF'));
idxm = find(ismember(varlist, 'MF'));
idxg = find(ismember(varlist, gname));

% Pre-allocate simulation matrices
y_mat_fut_ex_temp  = NaN(nParamDraws,nDraws,nvars,tperiods); 
y_mat_fut_ex_temp_good = NaN(nParamDraws,nDraws,nvars,tperiods); 
y_mat_fut_ex_temp_bad = NaN(nParamDraws,nDraws,nvars,tperiods); 

% Pre-allocate state matrices
stemp_t   = NaN(nParamDraws,nDraws,tperiods);

% Collect filtered probabilities bad regime (regime 2)
[~,~,~,f]=filter(sv);
% Compute p(st+1|st)
% p_reg2_filtered = f.filtered_regime_probabilities.regime2.data;

% Compute p(st=2)
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
        % Loop over periods
        if mod(uu,1000)==0
        fprintf('Current draw: %i of %i\n',uu,nDraws);
        end
        
        st_sim      = NaN(tperiods,1);
        st_filtered = NaN(tperiods,1);
        
        st_lag = 1;
        
        % Loop over periods
        for tt=1:tperiods
            
            % Simulate s(t) conditional on realized path of m and f
            % from t-hh+1 to t.
            if tt<=hh
                % In the first hh we can construct s(t) from the filtered probabilities
                
                % Use updated probability to determine current state
                udraw = rand(1); % draw a coin that determines the regime
                if udraw < (1-pbad(tt))
                    st_temp = 1;
                else
                    st_temp = 2;
                end
                
                % This is s(t) drown from p(st=2)
                st_sim(tt) = st_temp;
                                               
            else
                
                % Get current state
                % Here we have two options:
                % drawst = draw from filtered states every period
                % forecastst = recursively simulate st

                if strcmp(temp,'drawst')
                    % Use filtered probability to determine current state
                    udraw = rand(1); % draw a coin that determines the regime
                    if udraw < (1-pbad(tt))
                        st_temp = 1;
                    else
                        st_temp = 2; % bad regime
                    end                    

                elseif strcmp(temp,'forecastst')
                    
                    % s(t-hh)
                    st_lag = st_sim(tt-hh);

                    % Realized path of f and m from t-hh+1 to t
                    f_temp = YY(tt-hh+1:tt,idxf);
                    m_temp = YY(tt-hh+1:tt,idxm);

                    % Compute transition probabilities for s(t-hh+1),...,s(t)
                    [p12,p21] = fTranstionProb(param,m_temp,f_temp,opt);

                    p11 = ones(hh,1) - p12; % probability of remaining in normal
                    p22 = ones(hh,1) - p21; % probability of remaining in bad

                    % This returns s(t-hh+1)...s(t+1)
                    st_sim_temp = simulate_st(st_lag,p11,p22,hh);
                    
                    % This is s(t)
                    % This is the object we obtain in the VAR
                    st_sim(tt) = st_sim_temp(hh-1);

                    % This is s(t+1) [Changed by PCB 11/19]
                    % This is the object that enters in t+1:t+h forecast
                    st_temp = st_sim_temp(hh);

                end
            end
            
            % Simulate Model
            [y_temp,y_temp1,y_temp2] = simulate_DirectFull(st_temp,YY(tt,:)',param,opt);
            
            % Collect simulations
            for vv=1:nvars
                if vv==idxg
                y_mat_fut_ex_temp(dd,uu,vv,tt) = y_temp(vv)+trend_fit(tt);
                y_mat_fut_ex_temp_good(dd,uu,vv,tt) = y_temp1(vv)+trend_fit(tt);
                y_mat_fut_ex_temp_bad(dd,uu,vv,tt)  = y_temp2(vv)+trend_fit(tt);
                else
                y_mat_fut_ex_temp(dd,uu,vv,tt) = y_temp(vv);
                y_mat_fut_ex_temp_good(dd,uu,vv,tt) = y_temp1(vv);
                y_mat_fut_ex_temp_bad(dd,uu,vv,tt)  = y_temp2(vv);
                end
            end

            % Collect t state
            stemp_t(dd,uu,tt) = st_temp;
            
        end
    end
end
            
% Reshape matrices
y_mat_fut_ex      = reshape(y_mat_fut_ex_temp,[nDraws*nParamDraws,nvars,tperiods-nlags+1]);
y_mat_fut_ex_good = reshape(y_mat_fut_ex_temp_good,[nDraws*nParamDraws,nvars,tperiods-nlags+1]);
y_mat_fut_ex_bad  = reshape(y_mat_fut_ex_temp_bad,[nDraws*nParamDraws,nvars,tperiods-nlags+1]);


stmat_t = reshape(stemp_t,[nDraws*nParamDraws,tperiods-nlags+1]);


% Change st=2 (bad) to st = 0 for plotting
% Mean(st) represents probability of good regime
stmat_t(stmat_t==2) = 0;


out.st_t_mean= mean(stmat_t)'; % This is in t
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
