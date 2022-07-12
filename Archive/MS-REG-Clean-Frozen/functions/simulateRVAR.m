% simulateRVAR.m 
%
% Simulates the reduced form of the estimated structural VAR model. 
%
% Structural form: A0*Y(t+1) = C + A1*Y(t) + SIG*EPS(t+1)
% Reduced form: Y(t+1) = D + B*Y(t) + OMEGA*EPS(t+1)
%
% Note that A0, C, A1, SIG, could be time-varying depending on the model specification
%--------------------------------------------------------------------------

function y_out = simulateRVAR(st_temp,f_temp,m_temp,rf_params,dd,nsim)

   
    % Reduced form matrices
    D_sync_1 = rf_params.D_sync_1;
    D_sync_2 = rf_params.D_sync_2;
    
    B_sync_1 = rf_params.B_sync_1;
    B_sync_2 = rf_params.B_sync_2;
    
    O_sync_1 = rf_params.O_sync_1;
    O_sync_2 = rf_params.O_sync_2;
    

    % Matrix of shocks
    EPS = randn(2,nsim); 
    
    % Store GDP simulations
    Ynext = NaN(2,nsim);

    y_out = NaN(nsim,1);
    s_out = NaN(nsim,1);

    
    %Y0 = [f_temp;m_temp];
    
    % Simulate GDP based on simulated path for st and factors m(tt) and f(tt)
    for tt=1:nsim 
        
        % Get current state and observed variables
        Ymat = [f_temp;m_temp];
        st  = st_temp;

        if st==2
            % Bad regime
            Ynext(:,tt) = D_sync_2 + B_sync_2*Ymat + O_sync_2*EPS(:,tt);
        else
            % Good regime
            Ynext(:,tt) = D_sync_1 + B_sync_1*Ymat + O_sync_1*EPS(:,tt);
        end
            
        % Collect the GDP simulation
        y_out(tt) = Ynext(2,tt);

        % Forecast s(t+1)
        p12 = 1./(1+exp(param.a12(dd)-param.b12(dd)*f_temp+param.c12(dd)*m_temp));
        p21 = 1./(1+exp(param.a21(dd)+param.b21(dd)*f_temp-param.c21(dd)*m_temp));
        
        % Transition probabilities t+1|t
        p11 = 1 - p12; % probability of remaining in normal
        p22 = 1 - p21; % probability of remaining in bad

        udraw = rand(1);

        if st == 1 % started in good regime
            if udraw > p11
                s_out(tt) = 2; % switch from good to bad
            else
                s_out(tt) = 1; % don't switch and remain in good
            end
        else % start in bad regime
            if udraw > p22
                s_out(tt) = 1; % switch from bad to good
            else
                s_out(tt) = 2; % don't switch and remain in bad
            end
        end

        % Update states
        f_temp = Ynext(1,tt);
        m_temp = Ynext(2,tt);
        st_temp = s_out(tt);

    end
    
    %Ynext = [Y0 Ynext];