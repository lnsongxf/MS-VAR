% simulate_st.m 
%
% Simulate state for nsim periods starting from st_lag and given
% p11 and p22
%
%--------------------------------------------------------------------------


function s_out = simulate_st_shock(st_lag,p11,p22,st_shock,p11_shock,p22_shock,nsim)

s_out = NaN(nsim,2);

if length(p11)==nsim
    
    
    for tt=1:nsim

        udraw = rand(1);

        if st_lag == 1 % started in good regime
            if udraw > p11(tt)
                s_out(tt,1) = 2; % switch from good to bad
            else
                s_out(tt,1) = 1; % don't switch and remain in good
            end
        else % start in bad regime
            if udraw > p22(tt)
                s_out(tt,1) = 1; % switch from bad to good
            else
                s_out(tt,1) = 2; % don't switch and remain in bad
            end
        end

        if st_shock == 1 % started in good regime
            if udraw > p11_shock(tt)
                s_out(tt,2) = 2; % switch from good to bad
            else
                s_out(tt,2) = 1; % don't switch and remain in good
            end
        else % start in bad regime
            if udraw > p22_shock(tt)
                s_out(tt,2) = 1; % switch from bad to good
            else
                s_out(tt,2) = 2; % don't switch and remain in bad
            end
        end
        % Update lagged state
        st_lag = s_out(tt,1);
        st_shock = s_out(tt,2);
        
    end
    
else
    error('Dimension of input arguments do not match');
end

