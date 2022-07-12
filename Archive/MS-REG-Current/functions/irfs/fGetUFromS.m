%This function backs out ONE POSSIBLE (randomly chosen from infinitely many
%possible paths) path of u given initial conditions, shocks, and a set path for s.

%Recall: if s is 1 and u is greater than p11, you switch, if s is 2 and u
%is greater than p22 you switch. otherwise, stay in same state
function [u_out, p11_out, p22_out] = fGetUFromS(init_in, shocks_in, spath_in, opt_in, param_in)

nsim = numel(spath_in); %opt_in.hh;
EPS = shocks_in;

%Generate parameter matrices
D_sync_1 = param_in.D_sync_1;
D_sync_2 = param_in.D_sync_2;
B_sync_1 = param_in.B_sync_1;
B_sync_2 = param_in.B_sync_2;
O_sync_1 = param_in.O_sync_1;
O_sync_2 = param_in.O_sync_2;

% Initialize vector of observations
f_tminus1 = init_in.f0; %For first pass of loop, t-1 is initial value
m_tminus1  = init_in.m0; %For first pass of loop, t-1 is initial value
gdp_tminus1 = init_in.y0; %For first pass of loop, t-1 is initial value
s_tminus1 = init_in.s0; %For first pass of loop, s in t-1 is equal to s0
Y_tminus1 = [f_tminus1; m_tminus1; gdp_tminus1]; %For first pass of loop, Y (defined below) in t-1 is equal to Y0

%Initialize vectors for u, p11, and p22
u_out = NaN(1,nsim);
p11_out = NaN(1,nsim);
p22_out = NaN(1,nsim);

for tt=1:nsim
    %Calculate p12 and p21 for time t based on f and m in t-1
    if opt_in.const==1
        if opt_in.transprob==0
            p12 = param_in.a12;
            p21 = param_in.a21;
        elseif opt_in.transprob==1
            if opt_in.normal==1
                p12 = 1./(1 + exp(param_in.a12-param_in.b12*f_tminus1-param_in.c12*m_tminus1));
                p21 = 1./(1 + exp(param_in.a21-param_in.b21*f_tminus1-param_in.c21*m_tminus1));
            elseif opt_in.normal==0
                p12 = 1./(1 + exp(param_in.a12-param_in.b12*f_tminus1+param_in.c12*m_tminus1));
                p21 = 1./(1 + exp(param_in.a21+param_in.b21*f_tminus1-param_in.c21*m_tminus1));
            end
        end
    elseif opt_in.const==0
        if opt_in.transprob==0
            errordlg('opt.transprob=0 incompatible with opt.const=0');
        elseif opt_in.trasnprob==1
            if opt_in.normal==1
                p12 = 1./(1 + exp(-param_in.b12*f_tminus1-param_in.c12*m_tminus1));
                p21 = 1./(1 + exp(-param_in.b21*f_tminus1-param_in.c21*m_tminus1));
            elseif opt_in.normal==0
                p12 = 1./(1 + exp(-param_in.b12*f_tminus1+param_in.c12*m_tminus1));
                p21 = 1./(1 + exp(+param_in.b21*f_tminus1-param_in.c21*m_tminus1));
            end
        end
    end
    p11 = 1 - p12; % probability of remaining in normal
    p22 = 1 - p21; % probability of remaining in bad
    
    %Find possible values for u and draw u.  Must consider 4 cases based on
    %possible combinations of st and st-1
    if spath_in(tt)==1 & s_tminus1==1 %staying in state 1, so u should be less than p11
        u_out(tt) = unifrnd(0,p11);
    elseif spath_in(tt)==1 & s_tminus1==2 %moving from 2 to 1, so should be greater than p22
        u_out(tt) = unifrnd(p22,1);
    elseif spath_in(tt)==2 & s_tminus1==1 %moving from 1 to 2, so should be greater than p11
        u_out(tt) = unifrnd(p11,1);
    elseif spath_in(tt)==2 & s_tminus1==2 %staying in state 2, so should be less than p22
        u_out(tt) = unifrnd(0,p22);
    end
        
    %Simulate Yt (consisting of ft, mt, and yt) based on Yt-1.  Two cases,
    %since the matrices used depend on value of s currently
    if spath_in(tt)==1
        Yt = D_sync_1 + B_sync_1 * Y_tminus1 + O_sync_1 * EPS(:,tt);
    elseif spath_in(tt)==2
        Yt = D_sync_2 + B_sync_2 * Y_tminus1 + O_sync_2 * EPS(:,tt);
    end
    
    %Set values for next period.  Note that gdp_tminus1 individually is not
    %needed
    f_tminus1 = Yt(1);
    m_tminus1 = Yt(2);
    s_tminus1 = spath_in(tt);
    Y_tminus1 = Yt;
    
    %Collect p11 and p22
    p11_out(tt) = p11;
    p22_out(tt) = p22;
end

end