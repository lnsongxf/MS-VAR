function [param,bad] = scriptParams_SIMPLE(param,dd,modelspec)



%% modelspec defines how many equations in the VAR have switching coefficients

if modelspec==1
    % modelspec = 1: switching only in GDP equation, constant A1
    % switches = {'c(1)','c(2)','a0(2)','s(2)'};
    % label = 'CT_A0t2_A1_SIGt2';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2(dd); param.a0_2_1_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd); param.a0_2_1_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1(dd); param.c_2_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1(dd); param.c_2_1_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1(dd) 0; 0 param.s_2_2_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0; 0 param.s_2_2_sync_2(dd)];
    
    % Fixed coefficients matrix
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1(dd) param.a1_2_2(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1(dd) param.a1_2_2(dd)];
    
    
elseif modelspec==2
    % modelspec = 2: switching in GDP and FF equation, constant A1
    % switches = {'c','a0','s'};
    % label = 'CT_A0T_A1_SIGT';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2_sync_1(dd); param.a0_2_1_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2_sync_2(dd); param.a0_2_1_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1_sync_1(dd) 0; 0 param.s_2_2_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1_sync_2(dd) 0; 0 param.s_2_2_sync_2(dd)];
    
    % Fixed coefficients matrix
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1(dd) param.a1_2_2(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1(dd) param.a1_2_2(dd)];
    
elseif modelspec==3
    % modelspec = 3: switching in GDP and FF equation, constant Sigma
    % switches = {'c','a0','a1'};
    % label = 'CT_A0T_A1T_SIG';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2_sync_1(dd); param.a0_2_1_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2_sync_2(dd); param.a0_2_1_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1_sync_1(dd) param.a1_1_2_sync_1(dd); param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1_sync_2(dd) param.a1_1_2_sync_2(dd); param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)];
    
    % Fixed coefficients matrix
    SIG_sync_1 = [param.s_1_1(dd) 0; 0 param.s_2_2(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0; 0 param.s_2_2(dd)];
    
elseif modelspec==4
    % modelspec = 4: switching in GDP and FF equation
    % switches = {'c','a0','a1','s'};
    % label = 'CT_A0T_A1T_SIGT';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2_sync_1(dd); param.a0_2_1_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2_sync_2(dd); param.a0_2_1_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1_sync_1(dd) param.a1_1_2_sync_1(dd); param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1_sync_2(dd) param.a1_1_2_sync_2(dd); param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1_sync_1(dd) 0; 0 param.s_2_2_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1_sync_2(dd) 0; 0 param.s_2_2_sync_2(dd)];
    
elseif modelspec==5
    % modelspec = 5: switching in GDP and FF equation, constant A0
    % switches = {'c','a1','s'};
    % label = 'CT_A0_A1T_SIGT';
    
    % Switching matrices
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1_sync_1(dd) param.a1_1_2_sync_1(dd); param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1_sync_2(dd) param.a1_1_2_sync_2(dd); param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1_sync_1(dd) 0; 0 param.s_2_2_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1_sync_2(dd) 0; 0 param.s_2_2_sync_2(dd)];
    
    % Fixed coefficients matrix
    A0_sync_1 = [1 param.a0_1_2(dd); param.a0_2_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd); param.a0_2_1(dd) 1];
    
elseif modelspec==6
    % modelspec = 6: switching in GDP and FF equation, constant A1 and Sigma
    % switches = {'c','a0'};
    % label = 'CT_A0T_A1_SIG';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2_sync_1(dd); param.a0_2_1_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2_sync_2(dd); param.a0_2_1_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd)];
    
    % Fixed coefficients matrix
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1(dd) param.a1_2_2(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1(dd) param.a1_2_2(dd)];
    
    SIG_sync_1 = [param.s_1_1(dd) 0; 0 param.s_2_2(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0; 0 param.s_2_2(dd)];
    
elseif modelspec==7
    % modelspec = 1: switching only in GDP equation, constant Sigma
    % switches = {'c(2)','a0(2)','a1(2)'};
    % label = 'Ct2_A0t2_A1t2_SIG';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2(dd); param.a0_2_1_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd); param.a0_2_1_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1(dd); param.c_2_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1(dd); param.c_2_1_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)];
    
    % Fixed coefficients matrix
    SIG_sync_1 = [param.s_1_1(dd) 0; 0 param.s_2_2(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0; 0 param.s_2_2(dd)];
    
elseif modelspec==8
    % modelspec = 1: switching only in GDP equation
    % switches = {'c(2)','a0(2)','a1(2)','s(2)'};
    % label = 'Ct2_A0t2_A1t2_SIGt2';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2(dd); param.a0_2_1_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd); param.a0_2_1_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1(dd); param.c_2_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1(dd); param.c_2_1_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1(dd) 0; 0 param.s_2_2_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0; 0 param.s_2_2_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd); param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)];
    
end

%% Reduced form matrices
D_sync_1 = A0_sync_1\C_sync_1;
D_sync_2 = A0_sync_2\C_sync_2;

B_sync_1 = A0_sync_1\A1_sync_1;
B_sync_2 = A0_sync_2\A1_sync_2;

O_sync_1 = A0_sync_1\SIG_sync_1;
O_sync_2 = A0_sync_2\SIG_sync_2;


%% Collect matrices

param.A0_sync_1 = A0_sync_1;
param.A0_sync_2 = A0_sync_2;
param.C_sync_1 = C_sync_1;
param.C_sync_2 = C_sync_2;
param.A1_sync_1 = A1_sync_1;
param.A1_sync_2 = A1_sync_2;
param.SIG_sync_1 = SIG_sync_1;
param.SIG_sync_2 = SIG_sync_2;
param.D_sync_1 = D_sync_1;
param.D_sync_2 = D_sync_2;
param.B_sync_1 = B_sync_1;
param.B_sync_2 = B_sync_2;
param.O_sync_1 = O_sync_1;
param.O_sync_2 = O_sync_2;


%% Check if reg =1 is high growth, low vol
check_good = 0;
check_wrong = 0;
check_flip = 0;


if D_sync_1(end,1)>=D_sync_2(end,1)
    if O_sync_1(end,end)<O_sync_2(end,end)
        check_good=1; % if mean is larger, and s.d. is smaller
    else
        check_wrong=1; % if mean is larger, and s.d. is larger
    end
else
    if O_sync_1(end,end)<O_sync_2(end,end)
        check_wrong=1; % if mean is smaller, and s.d. is smaller
    else
        check_flip=1; % if mean is smaller, and s.d. is larger
    end
end

if check_good==1
    bad = 2;
elseif check_flip==1
    bad = 1;
end


if check_wrong==1
    disp('wrong')
    return
end
