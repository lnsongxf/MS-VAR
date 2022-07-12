% simulate_IteratedFull.m
%
% Simulates the reduced form of the estimated structural VAR model.
%
% This script also accomodates nlags>1 through companion form
%
% Structural form: A0*Y(t+1) = C + A1*Y(t) + SIG*EPS(t+1)
% Reduced form: Y(t+1) = D + B*Y(t) + OMEGA*EPS(t+1)
%
% Note that A0, C, A1, SIG, could be time-varying depending on the model specification
%--------------------------------------------------------------------------

function [y_out,y_out_1,y_out_2,s_out] = simulate_IteratedFull_companion(st_temp,Ymat,param,opt)

nsim = opt.hh;

% Matrix of shocks
EPS = randn(opt.nvars,nsim);

% Store GDP simulations
Ynext   = NaN(opt.nvars*opt.nlags,nsim);
Ynext_1 = NaN(opt.nvars*opt.nlags,nsim);
Ynext_2 = NaN(opt.nvars*opt.nlags,nsim);

y_out   = NaN(opt.nvars,nsim);
y_out_1 = NaN(opt.nvars,nsim);
y_out_2 = NaN(opt.nvars,nsim);
s_out   = NaN(1,nsim);

% Find index for financial factor and macro factor
idxf = find(ismember(opt.varlist, 'FF'));
idxm = find(ismember(opt.varlist, 'MF'));
f_temp = Ymat(:,idxf);
m_temp = Ymat(:,idxm);

%% Collect reduced-form matrices

D_sync_1 = param.D_sync_1;
D_sync_1_comp = [D_sync_1;zeros((opt.nlags-1)*opt.nvars,1)];
D_sync_2 = param.D_sync_2;
D_sync_2_comp = [D_sync_2;zeros((opt.nlags-1)*opt.nvars,1)];

% Form the companion, explicitly excluding deterministic terms
B_sync_1 = param.B_sync_1;
B_sync_2 = param.B_sync_2;
[B_sync_1_comp,~]=vartools.companion(B_sync_1,[],0);
[B_sync_2_comp,~]=vartools.companion(B_sync_2,[],0);

O_sync_1 = param.O_sync_1;
O_sync_2 = param.O_sync_2;

%% Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
for tt=1:nsim
    
    % Get current state
    Ymat = vec(Ymat');
    st  = st_temp;
    
    % Regime 2
    Ynext_2(:,tt) = D_sync_2_comp + B_sync_2_comp*Ymat + [O_sync_2*EPS(:,tt);zeros((opt.nlags-1)*opt.nvars,1)];
    
    % Regime 1
    Ynext_1(:,tt) = D_sync_1_comp + B_sync_1_comp*Ymat + [O_sync_1*EPS(:,tt);zeros((opt.nlags-1)*opt.nvars,1)];
    
    if st==2
        % Regime 2
        Ynext(:,tt) = Ynext_2(:,tt);
    else
        % Regime 1
        Ynext(:,tt) = Ynext_1(:,tt);
    end
    
    % Collect the simulation
    y_out(:,tt) = Ynext(1:opt.nvars,tt);
    y_out_1(:,tt) = Ynext_1(1:opt.nvars,tt);
    y_out_2(:,tt) = Ynext_2(1:opt.nvars,tt);
    
    
    % Forecast s(t+1)
    if opt.transprob ==1
        [p12,p21] = fTranstionProb(param,m_temp(end),f_temp(end),opt);
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
    f_temp = Ynext((idxf:opt.nvars:(opt.nlags*opt.nvars-2)),tt);
    m_temp = Ynext((idxm:opt.nvars:(opt.nlags*opt.nvars-1)),tt);
    
    % Update state
    st_temp = s_out(tt);
    
end