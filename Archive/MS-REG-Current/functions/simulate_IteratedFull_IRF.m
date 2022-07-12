% simulate_IteratedFull.m
%
% Simulates the reduced form of the estimated structural VAR model.
%
% Structural form: A0*Y(t+1) = C + A1*Y(t) + SIG*EPS(t+1)
% Reduced form: Y(t+1) = D + B*Y(t) + OMEGA*EPS(t+1)
%
% Note that A0, C, A1, SIG, could be time-varying depending on the model specification
%--------------------------------------------------------------------------

function [y_out,y_out_1,y_out_2,s_out,y_out_shock,y_out_1_shock,y_out_2_shock,s_out_shock] = simulate_IteratedFull_IRF(st_temp,FF_in,MF_in,GDP_in,param,opt)

nsim = opt.hh;
SHOCK = opt.shock;

% Matrix of shocks
EPS = randn(3,nsim);
EPS_SHOCK = EPS;
EPS_SHOCK(opt.which_shock,1) = EPS_SHOCK(opt.which_shock,1)+SHOCK;

% Store GDP simulations
Ynext   = NaN(3,nsim);
Ynext_1 = NaN(3,nsim);
Ynext_2 = NaN(3,nsim);
Ynext_1_condden = NaN(3,nsim);
Ynext_2_condden = NaN(3,nsim);
Ynext_shock   = NaN(3,nsim);
Ynext_1_shock = NaN(3,nsim);
Ynext_2_shock = NaN(3,nsim);
Ynext_1_shock_condden = NaN(3,nsim);
Ynext_2_shock_condden = NaN(3,nsim);

y_out   = NaN(nsim,3);
y_out_1 = NaN(nsim,3);
y_out_2 = NaN(nsim,3);
s_out   = NaN(nsim,1);
y_out_shock   = NaN(nsim,3);
y_out_1_shock = NaN(nsim,3);
y_out_2_shock = NaN(nsim,3);
s_out_shock   = NaN(nsim,1);

%% Collect reduced-form matrices

D_sync_1 = param.D_sync_1;
D_sync_2 = param.D_sync_2;
B_sync_1 = param.B_sync_1;
B_sync_2 = param.B_sync_2;
O_sync_1 = param.O_sync_1;
O_sync_2 = param.O_sync_2;

% Initialize vector of observations
f_temp = FF_in(1); f_temp_shock = f_temp;
m_temp  = MF_in(1); m_temp_shock = m_temp;
gdp_temp= GDP_in(1); gdp_temp_shock = gdp_temp;

% Initialize the st_temp_shock
st_temp_shock = st_temp;
%% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
for tt=1:nsim

%     if length(FF_in)==nsim
%         f_temp  = FF_in(tt);
%         m_temp  = MF_in(tt);
%         gdp_temp= GDP_in(tt);
%     end
    
    % Get current state and observed variables
    Ymat = [f_temp;m_temp;gdp_temp];
    Ymat_shock = [f_temp_shock;m_temp_shock;gdp_temp_shock];
    st  = st_temp;
    st_shock  = st_temp_shock;
    
    
    % Regime 2
    Ynext_2(:,tt) = D_sync_2 + B_sync_2*Ymat + O_sync_2*EPS(:,tt);
    Ynext_2_shock(:,tt) = D_sync_2 + B_sync_2*Ymat_shock + O_sync_2*EPS_SHOCK(:,tt);
    
    % Regime 1
    Ynext_1(:,tt) = D_sync_1 + B_sync_1*Ymat + O_sync_1*EPS(:,tt);
    Ynext_1_shock(:,tt) = D_sync_1 + B_sync_1*Ymat_shock + O_sync_1*EPS_SHOCK(:,tt);
    
    % CONDITIONAL DENSITIES
    % Regime 2
    if tt ==1
        Ynext_2_condden(:,tt) = D_sync_2 + B_sync_2*Ymat + O_sync_2*EPS(:,tt);
        Ynext_2_shock_condden(:,tt) = D_sync_2 + B_sync_2*Ymat_shock + O_sync_2*EPS_SHOCK(:,tt);
        
        % Regime 1
        Ynext_1_condden(:,tt) = D_sync_1 + B_sync_1*Ymat + O_sync_1*EPS(:,tt);
        Ynext_1_shock_condden(:,tt) = D_sync_1 + B_sync_1*Ymat_shock + O_sync_1*EPS_SHOCK(:,tt);
    else
        Ynext_2_condden(:,tt) = D_sync_2 + B_sync_2*Ynext_2_condden(:,tt-1) + O_sync_2*EPS(:,tt);
        Ynext_2_shock_condden(:,tt) = D_sync_2 + B_sync_2*Ynext_2_shock_condden(:,tt-1) + O_sync_2*EPS_SHOCK(:,tt);
        
        % Regime 1
        Ynext_1_condden(:,tt) = D_sync_1 + B_sync_1*Ynext_1_condden(:,tt-1) + O_sync_1*EPS(:,tt);
        Ynext_1_shock_condden(:,tt) = D_sync_1 + B_sync_1*Ynext_1_shock_condden(:,tt-1) + O_sync_1*EPS_SHOCK(:,tt);
    end
    
    if st==2
        % Regime 2
        Ynext(:,tt) = Ynext_2(:,tt);
    elseif st ==1
        % Regime 1
        Ynext(:,tt) = Ynext_1(:,tt);        
    end
     if st_shock==2
        % Regime 2
        Ynext_shock(:,tt) = Ynext_2_shock(:,tt);
     elseif st_shock==1
        % Regime 1
        Ynext_shock(:,tt) = Ynext_1_shock(:,tt);
    end   
    % Collect all simulations
    y_out(tt,:) = Ynext(:,tt);
%     y_out_1(tt,:) = Ynext_1(:,tt);
%     y_out_2(tt,:) = Ynext_2(:,tt);
    y_out_1(tt,:) = Ynext_1_condden(:,tt);
    y_out_2(tt,:) = Ynext_2_condden(:,tt);
    
%     y_out_shock(tt,:) = Ynext_shock(:,tt)-Ynext(:,tt);
%     y_out_1_shock(tt,:) = Ynext_1_shock(:,tt)-Ynext_1(:,tt);
%     y_out_2_shock(tt,:) = Ynext_2_shock(:,tt)-Ynext_2(:,tt);
    y_out_shock(tt,:) = Ynext_shock(:,tt);
%     y_out_1_shock(tt,:) = Ynext_1_shock(:,tt);
%     y_out_2_shock(tt,:) = Ynext_2_shock(:,tt);
    y_out_1_shock(tt,:) = Ynext_1_shock_condden(:,tt);
    y_out_2_shock(tt,:) = Ynext_2_shock_condden(:,tt);
    
    
    % Forecast s(t+1)
    if opt.transprob ==1
        [p12,p21] = fTranstionProb(param,m_temp,f_temp,opt);
        [p12_shock,p21_shock] = fTranstionProb(param,m_temp_shock,f_temp_shock,opt);
    elseif opt.transprob ==0
        p12 = param.a12; p12_shock = p12;
        p21 = param.a21; p21_shock = p21;
    end
    
        
    % Transition probabilities t+1|t
    p11 = 1 - p12; % probability of remaining in normal
    p22 = 1 - p21; % probability of remaining in bad
    p11_shock = 1 - p12_shock; % probability of remaining in normal
    p22_shock = 1 - p21_shock; % probability of remaining in bad
    
%     s_out(tt) = simulate_st(st,p11,p22,1); % PUT THE SAME COIN FOR BASELINE AND SHOCKED
%     s_out_shock(tt) = simulate_st(st_shock,p11_shock,p22_shock,1);
    st_temp_s = simulate_st_shock(st,p11,p22,st_shock,p11_shock,p22_shock,1);
    s_out(tt) = st_temp_s(1);
    s_out_shock(tt) = st_temp_s(2);
    
%     if length(FF_in)==1
        % Update states
        f_temp = Ynext(1,tt);
        m_temp = Ynext(2,tt);
        gdp_temp = Ynext(3,tt);
        f_temp_shock = Ynext_shock(1,tt);
        m_temp_shock = Ynext_shock(2,tt);
        gdp_temp_shock = Ynext_shock(3,tt);
%     end
    st_temp = s_out(tt);
    st_temp_shock = s_out_shock(tt);
    
end

%Ynext = [Y0 Ynext];