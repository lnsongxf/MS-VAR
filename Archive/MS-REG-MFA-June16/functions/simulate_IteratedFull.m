% simulate_IteratedFull.m
%
% Simulates the reduced form of the estimated structural VAR model.
%
% Structural form: A0*Y(t+1) = C + A1*Y(t) + SIG*EPS(t+1)
% Reduced form: Y(t+1) = D + B*Y(t) + OMEGA*EPS(t+1)
%
% Note that A0, C, A1, SIG, could be time-varying depending on the model specification
%--------------------------------------------------------------------------

function [y_out,y_out_1,y_out_2,s_out] = simulate_IteratedFull(st_temp,Ymat,param,opt)
         
% Number of periods
nsim = opt.hh;

% Matrix of shocks
EPS = randn(opt.nvars,nsim);

% Store GDP simulations
Ynext   = NaN(opt.nvars,nsim);
Ynext_1 = NaN(opt.nvars,nsim);
Ynext_2 = NaN(opt.nvars,nsim);

y_out   = NaN(opt.nvars,nsim);
y_out_1 = NaN(opt.nvars,nsim);
y_out_2 = NaN(opt.nvars,nsim);
s_out   = NaN(1,nsim);

% Find index for financial factor and macro factor
idxf = find(ismember(opt.varlist, 'FF'));
idxm = find(ismember(opt.varlist, 'MF'));
f_temp = Ymat(idxf,:);
m_temp = Ymat(idxm,:);

%% Collect reduced-form matrices

D_sync_1 = param.D_sync_1;
D_sync_2 = param.D_sync_2;
B_sync_1 = param.B_sync_1;
B_sync_2 = param.B_sync_2;
O_sync_1 = param.O_sync_1;
O_sync_2 = param.O_sync_2;

%% Simulate variables based on simulated path for st and factors m(tt) and f(tt)

for tt=1:nsim
    
    % Get current state and observed variables
    st  = st_temp;
    
    % Regime 2
    Ynext_2(:,tt) = D_sync_2 + B_sync_2*Ymat + O_sync_2*EPS(:,tt);
    
    % Regime 1
    Ynext_1(:,tt) = D_sync_1 + B_sync_1*Ymat + O_sync_1*EPS(:,tt);
    
    if st==2
        % Regime 2
        Ynext(:,tt) = Ynext_2(:,tt);
    else
        % Regime 1
        Ynext(:,tt) = Ynext_1(:,tt);
    end
    
    % Collect the simulation
    y_out(:,tt) = Ynext(:,tt);
    y_out_1(:,tt) = Ynext_1(:,tt);
    y_out_2(:,tt) = Ynext_2(:,tt);
    
    
    % Forecast s(t+1)
    if opt.transprob ==1
        [p12,p21] = fTranstionProb(param,m_temp,f_temp,opt);
    elseif opt.transprob ==0
        p12 = param.a12; 
        p21 = param.a21; 
    end
            
    % Transition probabilities t+1|t
    p11 = 1 - p12; % probability of remaining in normal
    p22 = 1 - p21; % probability of remaining in bad
    
    s_out(tt) = simulate_st(st,p11,p22,1);
    
    % Update variables
    Ymat = Ynext(:,tt);
    f_temp = Ymat(idxf,:);
    m_temp = Ymat(idxm,:);

    % Update state
    st_temp = s_out(tt);
    
end