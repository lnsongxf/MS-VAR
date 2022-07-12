%==========================================================================
%  This function simulates the reduced form VAR
%
%  Y(t) = D(s(t)) + B(s(t))*Y(t-1) + O(s(t))*eps(t)
%
% 
%  Where Y(0) = (I-B)\D
%
% p12 = p(s(t+1=1)|s(t=2)) = 1/(1+exp(a12 + b12*f(t) + c13*m(t))
% p21 = p(s(t+1=2)|s(t=1)) = 1/(1+exp(a21 + b21*f(t) + c21*m(t))
%==========================================================================
%
%------------------------------------------------------------------
function [Y_out_all, Y_1_all, Y_2_all, s_out_all,eta,p12_out_all,p21_out_all] = fSimulateDGPIteretedFull(param_in,opt)


tsim = opt.nSim+opt.burn;

% Pre-allocate simulation matrices
Y_out_all  = NaN(3,opt.nSim,opt.nSamples);
Y_out_1  = NaN(3,opt.nSim);
Y_out_2  = NaN(3,opt.nSim);

Y_1_all  = NaN(3,opt.nSim,opt.nSamples);
Y_2_all  = NaN(3,opt.nSim,opt.nSamples);

p12_out_all= NaN(opt.nSim,opt.nSamples);

p21_out_all= NaN(opt.nSim,opt.nSamples);
s_out_all  = NaN(opt.nSim,opt.nSamples);
    eta    = NaN(3,opt.nSim,opt.nSamples);


% Check regime indicator and collect VAR matrices
%param = fMapParamsDGP(param_in,dd);

% Reduced form parameters
D_sync_1 = param_in.D_sync_1;
D_sync_2 = param_in.D_sync_2;
B_sync_1 = param_in.B_sync_1;
B_sync_2 = param_in.B_sync_2;
O_sync_1 = param_in.O_sync_1;
O_sync_2 = param_in.O_sync_2;

p12 = 1/(1+exp(param_in.a12));
p11 = 1-p12;

% Loop over repetitions
for uu=1:opt.nSamples
    % Loop over periods
    if mod(uu,1000)==0
        fprintf('Current sample draw: %i of %i\n',uu,opt.nSamples);
    end
    
    % Draw initial state
    udraw0 = rand(1); % draw a coin that determines the regime
    
    if udraw0 < p11
        s_0 = 1; % good regime
    else
        s_0 = 2; % bad regime
    end
    
    % Initialize vector of observations
    if s_0 ==1
        Y_0 = (eye(3)-B_sync_1)\D_sync_1;
    else
        Y_0 = (eye(3)-B_sync_2)\D_sync_2;
    end
    
    %% Clean matrices
    Y_out   = NaN(3,tsim);
    s_out   = NaN(tsim,1);
    p12_out   = NaN(tsim,1);
    p21_out   = NaN(tsim,1);
    %% Draw shocks
    %EPS = randn(3,tsim);    
    %eta_in(:,1) = zeros(3,1);
    %% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
    % Observed Y0 = (m0,f0,gdp0,s0);
    
    for tt=1:tsim
        
        
        eta(:,tt,uu) = randn(3,1); %eta_in(:,tt); %randn(3,1);
        
        EPS = eta(:,tt,uu);
        
        % Simulate s(t) given observed information through t-1
        
        % Financial and macro factor at time t-1
        if tt==1
            f_lag = Y_0(1);
            m_lag = Y_0(2);
            Y_lag = Y_0;
            st_lag= s_0;
        else
            f_lag = Y_out(1,tt-1);
            m_lag = Y_out(2,tt-1);
            Y_lag = Y_out(:,tt-1);
            st_lag= s_out(tt-1);

        end
        
        % Get transition probabilities p(s(t)=i |s(t-1) = j)
        [p12_temp,p21_temp] = fTranstionProb(param_in,m_lag,f_lag,opt);
        
        % Transition probabilities t+1|t
        p11 = 1 - p12_temp; % probability of remaining in normal
        p22 = 1 - p21_temp; % probability of remaining in bad
        
        % Draw the state at time t
        st = simulate_st(st_lag,p11,p22,1);        
        
        % Regime 2
        Y_out_2(:,tt) = D_sync_2 + B_sync_2*Y_lag + O_sync_2*EPS;
        
        % Regime 1
        Y_out_1(:,tt) = D_sync_1 + B_sync_1*Y_lag + O_sync_1*EPS;
        
        
        
        
        if st==2
            % Regime 2
            Y_out(:,tt) = Y_out_2(:,tt);
        else
            % Regime 1
            Y_out(:,tt) = Y_out_1(:,tt);
        end

        % Store the state
        s_out (tt) = st;
        
        % p12
        p12_out(tt) = p12_temp;
        % p21
        p21_out(tt) = p21_temp;

        
    end
    %[Y_out,s_out] = simulate_IteratedFullDGP(s_0,Y_0,param,opt);
    
    % Change st=2 (bad) to st = 0 for plotting
    %s_out(s_out==2) = 0;
    
    Y_out_all(:,:,uu) = Y_out(:,opt.burn+1:end);    
    Y_1_all(:,:,uu) = Y_out_1(:,opt.burn+1:end);    
    Y_2_all(:,:,uu) = Y_out_2(:,opt.burn+1:end);    
    s_out_all(:,uu)   = s_out(opt.burn+1:end);
    
    % p12
    p12_out_all(:,uu) = p12_out(opt.burn+1:end);
    % p21
    p21_out_all(:,uu) = p21_out(opt.burn+1:end);

end



