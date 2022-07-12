function [Y_out, s_out,p12_out,p21_out] = fSimulateIterated(init_in, shocks_in, param_in, opt_in, varargin)
    
nsim = opt_in.hh;

% Matrix of shocks
EPS = shocks_in;

% Store GDP simulations
Ynext   = NaN(3,nsim);
Ynext_1 = NaN(3,nsim);
Ynext_2 = NaN(3,nsim);

s_out   = NaN(nsim,1);
%% Collect reduced-form matrices

D_sync_1 = param_in.D_sync_1;
D_sync_2 = param_in.D_sync_2;
B_sync_1 = param_in.B_sync_1;
B_sync_2 = param_in.B_sync_2;
O_sync_1 = param_in.O_sync_1;
O_sync_2 = param_in.O_sync_2;

% Initialize vector of observations
f_temp = init_in.f0;
m_temp  = init_in.m0;
gdp_temp = init_in.y0;
st_temp = init_in.s0;

if isempty(varargin)
    udraw = NaN(1,nsim);
else
    udraw = varargin{:};
end
   

%% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
for tt=1:nsim
    
    % Get current state and observed variables
    Ymat = [f_temp;m_temp;gdp_temp];                % This is Y(t-1)
    st   = st_temp;                                 % This is s(t-1)
    
    %---------------------
    % Forecast s(t)|s(t-1)
    %---------------------
    if opt_in.transprob ==1
        [p12,p21] = fTranstionProb(param_in,m_temp,f_temp,opt_in);
    elseif opt_in.transprob ==0
        p12 = param.a12; 
        p21 = param.a21; 
    end
    
        
    % Transition probabilities t|t-1
    p11 = 1 - p12; % probability of remaining in normal
    p22 = 1 - p21; % probability of remaining in bad
    
    % Draw s(t)
    s_out(tt) = simulate_st(st,p11,p22,1,udraw(tt));

    %------------------------------
    % Compute Conditional Outcomes
    %------------------------------
    
    % Regime 2
    Ynext_2(:,tt) = D_sync_2 + B_sync_2*Ymat + O_sync_2*EPS(:,tt);
    
    % Regime 1
    Ynext_1(:,tt) = D_sync_1 + B_sync_1*Ymat + O_sync_1*EPS(:,tt);
    
    
    %----------------------------------------
    % Select the outcome based on the regime
    %----------------------------------------
    
    if s_out(tt)==2
        % Regime 2
        Ynext(:,tt) = Ynext_2(:,tt);
    else
        % Regime 1
        Ynext(:,tt) = Ynext_1(:,tt);
    end
    
    % Collect the GDP simulation
    Y_out(:,tt) = Ynext(:,tt);
    
    % Update states
    f_temp   = Ynext(1,tt);     % This is f(t)
    m_temp   = Ynext(2,tt);     % This is m(t)
    gdp_temp = Ynext(3,tt);     % This is y(t)
    st_temp  = s_out(tt);       % This is s(t)
    
    % Store transition prob   
    p12_out(tt) = p12;          % This is p12(t)
    p21_out(tt) = p21;          % This is p21(t)
    
end    

end