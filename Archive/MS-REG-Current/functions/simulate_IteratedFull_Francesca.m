% simulate_IteratedFull.m
%
% Simulates the reduced form of the estimated structural VAR model.
%
% Structural form: A0*Y(t+1) = C + A1*Y(t) + SIG*EPS(t+1)
% Reduced form: Y(t+1) = D + B*Y(t) + OMEGA*EPS(t+1)
%
% Note that A0, C, A1, SIG, could be time-varying depending on the model specification
%--------------------------------------------------------------------------

function [y_out,y_out_1,y_out_2] = simulate_IteratedFull(st_temp,f_temp,m_temp,gdp_temp,param,opt)

nsim = opt.hh;

% Matrix of shocks
EPS = randn(3,nsim);

% Store GDP simulations
Ynext   = NaN(3,nsim);
Ynext_1 = NaN(3,nsim);
Ynext_2 = NaN(3,nsim);

y_out   = NaN(nsim,1);
y_out_1 = NaN(nsim,1);
y_out_2 = NaN(nsim,1);
s_out   = NaN(nsim,1);

Ymat0 = [f_temp;m_temp;gdp_temp];


%% Collect reduced-form matrices

D_sync_1 = param.D_sync_1;
D_sync_2 = param.D_sync_2;
B_sync_1 = param.B_sync_1;
B_sync_2 = param.B_sync_2;
O_sync_1 = param.O_sync_1;
O_sync_2 = param.O_sync_2;


%% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
% Using reduced form matrices

for tt=1:nsim
    
    % Get current state and observed variables
    Ymat = [f_temp;m_temp;gdp_temp];
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
    
    % Collect the GDP simulation
    y_out(tt) = Ynext(end,tt);
    y_out_1(tt) = Ynext_1(end,tt);
    y_out_2(tt) = Ynext_2(end,tt);
    
    
    % Forecast s(t+1)
    [p12,p21] = fTranstionProb(param,m_temp,f_temp,opt);
    
        
    % Transition probabilities t+1|t
    p11 = 1 - p12; % probability of remaining in normal
    p22 = 1 - p21; % probability of remaining in bad
    
    s_out(tt) = simulate_st(st,p11,p22,1);
    
    % Update states
    f_temp = Ynext(1,tt);
    m_temp = Ynext(2,tt);
    st_temp = s_out(tt);
    
end



%% Collect structural form matrices
A0_sync_1 = params.A0_sync_1;
A0_sync_2 = params.A0_sync_2;
C_sync_1  = params.C_sync_1;
C_sync_2  = params.C_sync_2;
A1_sync_1 = params.A1_sync_1;
A1_sync_2 = params.A1_sync_2;
SIG_sync_1 = params.SIG_sync_1;
SIG_sync_2 = params.SIG_sync_2;

%% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
% Using structural form matrices

for tt=1:nsim
    
    % Get current state and observed variables
    Ymat = [f_temp;m_temp;gdp_temp];
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
    
    % Collect the GDP simulation
    y_out(tt) = Ynext(end,tt);
    y_out_1(tt) = Ynext_1(end,tt);
    y_out_2(tt) = Ynext_2(end,tt);
    
    
    % Forecast s(t+1)
    [p12,p21] = fTranstionProb(param,m_temp,f_temp,opt);
    
        
    % Transition probabilities t+1|t
    p11 = 1 - p12; % probability of remaining in normal
    p22 = 1 - p21; % probability of remaining in bad
    
    s_out(tt) = simulate_st(st,p11,p22,1);
    
    % Update states
    f_temp = Ynext(1,tt);
    m_temp = Ynext(2,tt);
    st_temp = s_out(tt);
    
end

%Ynext = [Y0 Ynext];