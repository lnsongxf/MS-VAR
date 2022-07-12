% simulate_st.m 
%
% Simulate state for nsim periods starting from st_lag and given
% p11 and p22
%
%--------------------------------------------------------------------------


function s_out = simulate_st(st_lag,p11,p22,nsim,varargin)

   

s_out = NaN(nsim,1);

if length(p11)==nsim
    
    
    for tt=1:nsim
        if isnan(varargin{1}(tt)) %right hand side of this assertion should be number of nonoptional inputs for the code to work
            udraw = rand(1);
        else
            udraw = varargin{1}(tt);
        end

        if st_lag == 1 % started in good regime
            if udraw > p11(tt)
                s_out(tt) = 2; % switch from good to bad
            else
                s_out(tt) = 1; % don't switch and remain in good
            end
        else % start in bad regime
            if udraw > p22(tt)
                s_out(tt) = 1; % switch from bad to good
            else
                s_out(tt) = 2; % don't switch and remain in bad
            end
        end
        
        % Update lagged state
        st_lag = s_out(tt);
        
    end
    
else
    error('Dimension of input arguments do not match');
end

