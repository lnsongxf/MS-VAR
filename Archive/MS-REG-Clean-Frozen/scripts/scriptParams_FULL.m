function [params] = scriptParams_FULL(param,dd,modelspec)


%% modelspec defines how many equations in the VAR have switching coefficients
if modelspec==1
    % modelspec = 1: switching only in GDP equation, constant A1
    % switches = {'c(3)','a0(3)','s(3)'};
    % label = 'CT_A0t2_A1_SIGt2';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1(dd); param.c_2_1(dd); param.c_3_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1(dd); param.c_2_1(dd); param.c_3_1_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3_sync_2(dd)];
    
    % Fixed coefficients matrix
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1(dd) param.a1_3_2(dd)  param.a1_3_3(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1(dd) param.a1_3_2(dd)  param.a1_3_3(dd)];
    
    
elseif modelspec==2 || modelspec==10 || modelspec==11
    % modelspec = 2: switching in GDP and FF equation, constant A1
    % switches = {'c','a0','s'};
    % label = 'CT_A0T_A1_SIGT';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2_sync_1(dd) param.a0_1_3_sync_1(dd); param.a0_2_1_sync_1(dd) 1 param.a0_2_3_sync_1(dd); param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2_sync_2(dd) param.a0_1_3_sync_2(dd); param.a0_2_1_sync_2(dd) 1 param.a0_2_3_sync_2(dd); param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd); param.c_3_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd); param.c_3_1_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1_sync_1(dd) 0 0; 0 param.s_2_2_sync_1(dd) 0; 0 0 param.s_3_3_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1_sync_2(dd) 0 0; 0 param.s_2_2_sync_2(dd) 0; 0 0 param.s_3_3_sync_2(dd)];
    
    % Fixed coefficients matrix
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1(dd) param.a1_3_2(dd)  param.a1_3_3(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1(dd) param.a1_3_2(dd)  param.a1_3_3(dd)];
    
elseif modelspec==3
    % modelspec = 3: switching in GDP and FF equation, constant Sigma
    % switches = {'c','a0','a1'};
    % label = 'CT_A0T_A1T_SIG';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2_sync_1(dd) param.a0_1_3_sync_1(dd); param.a0_2_1_sync_1(dd) 1 param.a0_2_3_sync_1(dd); param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2_sync_2(dd) param.a0_1_3_sync_2(dd); param.a0_2_1_sync_2(dd) 1 param.a0_2_3_sync_2(dd); param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd); param.c_3_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd); param.c_3_1_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1_sync_1(dd) param.a1_1_2_sync_1(dd) param.a1_1_3_sync_1(dd); param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)  param.a1_2_3_sync_1(dd); param.a1_3_1_sync_1(dd) param.a1_3_2_sync_1(dd)  param.a1_3_3_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1_sync_2(dd) param.a1_1_2_sync_2(dd) param.a1_1_3_sync_2(dd); param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)  param.a1_2_3_sync_2(dd); param.a1_3_1_sync_2(dd) param.a1_3_2_sync_2(dd)  param.a1_3_3_sync_2(dd)];
    
    % Fixed coefficients matrix
    SIG_sync_1 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3(dd)];
    
elseif modelspec==4 || modelspec==12
    % modelspec = 4: switching in GDP and FF equation
    % switches = {'c','a0','a1','s'};
    % label = 'CT_A0T_A1T_SIGT';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2_sync_1(dd) param.a0_1_3_sync_1(dd); 
                param.a0_2_1_sync_1(dd) 1 param.a0_2_3_sync_1(dd); 
                param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];
            
    A0_sync_2 = [1 param.a0_1_2_sync_2(dd) param.a0_1_3_sync_2(dd); 
                param.a0_2_1_sync_2(dd) 1 param.a0_2_3_sync_2(dd); 
                param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd); param.c_3_1_sync_1(dd)];
    
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd); param.c_3_1_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1_sync_1(dd) param.a1_1_2_sync_1(dd) param.a1_1_3_sync_1(dd); 
                param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)  param.a1_2_3_sync_1(dd); 
                param.a1_3_1_sync_1(dd) param.a1_3_2_sync_1(dd)  param.a1_3_3_sync_1(dd)];
            
    A1_sync_2 = [ param.a1_1_1_sync_2(dd) param.a1_1_2_sync_2(dd) param.a1_1_3_sync_2(dd); 
                param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)  param.a1_2_3_sync_2(dd); 
                param.a1_3_1_sync_2(dd) param.a1_3_2_sync_2(dd)  param.a1_3_3_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1_sync_1(dd) 0 0; 
                  0 param.s_2_2_sync_1(dd) 0; 
                  0 0 param.s_3_3_sync_1(dd)];
        
    SIG_sync_2 = [param.s_1_1_sync_2(dd) 0 0; 
                 0 param.s_2_2_sync_2(dd) 0; 
                 0 0 param.s_3_3_sync_2(dd)];
    
elseif modelspec==5
    % modelspec = 5: switching in GDP and FF equation, constant A0
    % switches = {'c','a1','s'};
    % label = 'CT_A0_A1T_SIGT';
    
    % Switching matrices
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd); param.c_3_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd); param.c_3_1_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1_sync_1(dd) param.a1_1_2_sync_1(dd) param.a1_1_3_sync_1(dd); param.a1_2_1_sync_1(dd) param.a1_2_2_sync_1(dd)  param.a1_2_3_sync_1(dd); param.a1_3_1_sync_1(dd) param.a1_3_2_sync_1(dd)  param.a1_3_3_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1_sync_2(dd) param.a1_1_2_sync_2(dd) param.a1_1_3_sync_2(dd); param.a1_2_1_sync_2(dd) param.a1_2_2_sync_2(dd)  param.a1_2_3_sync_2(dd); param.a1_3_1_sync_2(dd) param.a1_3_2_sync_2(dd)  param.a1_3_3_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1_sync_1(dd) 0 0; 0 param.s_2_2_sync_1(dd) 0; 0 0 param.s_3_3_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1_sync_2(dd) 0 0; 0 param.s_2_2_sync_2(dd) 0; 0 0 param.s_3_3_sync_2(dd)];
    
    % Fixed coefficients matrix
    A0_sync_1 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1(dd) param.a0_3_2(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1(dd) param.a0_3_2(dd) 1];
    
elseif modelspec==6
    % modelspec = 6: switching in GDP and FF equation, constant A1 and Sigma
    % switches = {'c','a0'};
    % label = 'CT_A0T_A1_SIG';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2_sync_1(dd) param.a0_1_3_sync_1(dd); param.a0_2_1_sync_1(dd) 1 param.a0_2_3_sync_1(dd); param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2_sync_2(dd) param.a0_1_3_sync_2(dd); param.a0_2_1_sync_2(dd) 1 param.a0_2_3_sync_2(dd); param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1_sync_1(dd); param.c_2_1_sync_1(dd); param.c_3_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1_sync_2(dd); param.c_2_1_sync_2(dd); param.c_3_1_sync_2(dd)];
    
    % Fixed coefficients matrix
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1(dd) param.a1_3_2(dd)  param.a1_3_3(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1(dd) param.a1_3_2(dd)  param.a1_3_3(dd)];
    
    SIG_sync_1 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3(dd)];
    
elseif modelspec==7
    % modelspec = 7: switching only in GDP equation, constant Sigma
    % switches = {'c(3)','a0(3)','a1(3)'};
    % label = 'Ct2_A0t2_A1t2_SIG';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1(dd); param.c_2_1(dd); param.c_3_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1(dd); param.c_2_1(dd); param.c_3_1_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1_sync_1(dd) param.a1_3_2_sync_1(dd)  param.a1_3_3_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1_sync_2(dd) param.a1_3_2_sync_2(dd)  param.a1_3_3_sync_2(dd)];
    
    % Fixed coefficients matrix
    SIG_sync_1 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3(dd)];
    
elseif modelspec==8
    % modelspec = 8: switching only in GDP equation
    % switches = {'c(3)','a0(3)','a1(3)','s(3)'};
    % label = 'Ct2_A0t2_A1t2_SIGt2';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1(dd); param.c_2_1(dd); param.c_3_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1(dd); param.c_2_1(dd); param.c_3_1_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1_sync_1(dd) param.a1_3_2_sync_1(dd)  param.a1_3_3_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1_sync_2(dd) param.a1_3_2_sync_2(dd)  param.a1_3_3_sync_2(dd)];

elseif modelspec==9
    % modelspec = 8: switching in GDP and FF equation
    % switches = {'c(3)','a0(3)','a1(3)','s(3)'};
    % label = 'Ct2_A0t2_A1t2_SIGt2_restr';
    
    % Switching matrices
    A0_sync_1 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1_sync_1(dd) param.a0_3_2_sync_1(dd) 1];
    A0_sync_2 = [1 param.a0_1_2(dd) param.a0_1_3(dd); param.a0_2_1(dd) 1 param.a0_2_3(dd); param.a0_3_1_sync_2(dd) param.a0_3_2_sync_2(dd) 1];
    
    C_sync_1  = [param.c_1_1(dd); param.c_2_1(dd); param.c_3_1_sync_1(dd)];
    C_sync_2  = [param.c_1_1(dd); param.c_2_1(dd); param.c_3_1_sync_2(dd)];
    
    SIG_sync_1 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3_sync_1(dd)];
    SIG_sync_2 = [param.s_1_1(dd) 0 0; 0 param.s_2_2(dd) 0; 0 0 param.s_3_3_sync_2(dd)];
    
    A1_sync_1 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1_sync_1(dd) param.a1_3_2_sync_1(dd)  param.a1_3_3_sync_1(dd)];
    A1_sync_2 = [ param.a1_1_1(dd) param.a1_1_2(dd) param.a1_1_3(dd); param.a1_2_1(dd) param.a1_2_2(dd)  param.a1_2_3(dd); param.a1_3_1_sync_2(dd) param.a1_3_2_sync_2(dd)  param.a1_3_3_sync_2(dd)];
end

%% Reduced form matrices
D_sync_1 = A0_sync_1\C_sync_1;
D_sync_2 = A0_sync_2\C_sync_2;

B_sync_1 = A0_sync_1\A1_sync_1;
B_sync_2 = A0_sync_2\A1_sync_2;

O_sync_1 = A0_sync_1\SIG_sync_1;
O_sync_2 = A0_sync_2\SIG_sync_2;

%% Collect matrices

params.A0_sync_1 = A0_sync_1;
params.A0_sync_2 = A0_sync_2;
params.C_sync_1 = C_sync_1;
params.C_sync_2 = C_sync_2;
params.A1_sync_1 = A1_sync_1;
params.A1_sync_2 = A1_sync_2;
params.SIG_sync_1 = SIG_sync_1;
params.SIG_sync_2 = SIG_sync_2;
params.D_sync_1 = D_sync_1;
params.D_sync_2 = D_sync_2;
params.B_sync_1 = B_sync_1;
params.B_sync_2 = B_sync_2;
params.O_sync_1 = O_sync_1;
params.O_sync_2 = O_sync_2;
params.a12 = param.a12(dd);
params.a21 = param.a21(dd);
params.b12 = param.b12(dd);
params.b21 = param.b21(dd);
params.c12 = param.c12(dd);
params.c21 = param.c21(dd);




