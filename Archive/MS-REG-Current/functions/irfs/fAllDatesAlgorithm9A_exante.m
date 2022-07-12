function [irf_out, pctl_out] = fAllDatesAlgorithm9A_exante(datesfull_in, datestarget_in, FFmat_in, MFmat_in, GDPGmat_in, sdata_in, psdata_in, init_in, shocks_in, param_in, opt_in, sd_base_in, sd_shock_in)


% path of u for regime selection
rng(12553)

% We set upath_shocks to zero to make sure there is no regime change
upath_shocks = rand(opt_in.nDrawsEPS,opt_in.hh);

for d = 1:length(datesfull_in)
    if any(strcmp(datesfull_in{d},datestarget_in)) %check if the current date is one of the target dates
        
        % Get current date
        current_date = datetime(datesfull_in{d}, 'InputFormat', 'dd-MMM-yyyy');
        
        % Transform current date to YYYYMM format
        yearmonth = strcat('Y', num2str(year(datetime(current_date))), 'M', num2str(month(datetime(current_date))));
        
        %Set path for s and initial values based on current date index
        init_in.f0 = FFmat_in(d-1);
        init_in.m0 = MFmat_in(d-1);
        init_in.y0 = GDPGmat_in(d-1);
        
        %Get the probability of s(t) in the period of the shock:
        
        % For this experiment, want probability of bad state to be same as initial value for all periods
        pbad_t = psdata_in(d);
        pspath = [];
        % Fix the regime based on threshold
        if pbad_t<0.2
            init_in.s0 = 1;
        else
            init_in.s0 = 2;
        end
        
        fprintf('p(bad) = %4.4f \n',pbad_t);
%         pspath = pbad_t*ones(1,opt_in.hh); 
        
        
        %Run algorithm 9, ex-ante; 
        [irf_temp, pctl_temp] = fAlgorithm9A_exante(init_in,shocks_in,param_in,opt_in,pspath,sd_base_in,sd_shock_in,upath_shocks);
        
        irf_out.(yearmonth) = irf_temp;
        pctl_out.(yearmonth) = pctl_temp;

        %Add extra information to pctl object
        pctl_out.(yearmonth).f0 = init_in.f0;
        pctl_out.(yearmonth).m0 = init_in.m0;
        pctl_out.(yearmonth).y0 = init_in.y0;
        pctl_out.(yearmonth).s0 = init_in.s0;
        pctl_out.(yearmonth).pspath = pspath;
    else
        continue
    end
end

end