if compare ==2
    if exist([sim_folder 'results_iterated.mat'],'file')==2
        load([sim_folder 'results_iterated.mat'])
        if exist([sim_folder 'results_iterated_samp.mat'],'file')==2;load([sim_folder 'results_iterated_samp.mat']);end
        
    
        %% EXAMPLE DRAW
        
        st_mat = st_fit'; % Now use the fitted values
        y_mat = y_fit'; % Now use the fitted values
        
        nDraws      = 10;           % Number of simulations per parameter draw
        nParamDraws = samp.options.N;    % Number of parameter draws
        y_mat_fut_ex_temp = NaN(nDraws,nParamDraws,size(y_mat,2)); % This stores all realizations for the ITERATED
        st_it_temp = NaN(nDraws,nParamDraws,size(y_mat,2));
        st_it_temp_t = NaN(nDraws,nParamDraws,size(y_mat,2));
        h = waitbar(0,'Calculating paths for the iterated version');
        y_temp_par1 = y_mat(example_draw,:);
        ff_temp_par1 = Res_iterated.FF(example_draw,:);
        mf_temp_par1 = Res_iterated.MF(example_draw,:);
        st_temp_par1 = st_mat(example_draw,:);
        for jj=1:100
            waitbar(jj/size(y_mat,2),h,'Calculating paths for the iterated version')
            y_temp_par = y_temp_par1(jj);
            ff_temp_par = ff_temp_par1(jj);
            mf_temp_par = mf_temp_par1(jj);
            if jj==1
                st_temp_par = st_temp_par1(jj);
            else
                st_temp_par = st_temp_par1(jj-1);
            end
            for dd=1:nParamDraws
                for uu=1:nDraws
                    y_temp = [y_temp_par;NaN(12,1)];
                    f_temp = [ff_temp_par;NaN(12,1)];
                    m_temp = [mf_temp_par;NaN(12,1)];
                    s_temp = NaN(13,1); % Simulate st; vfist
                    
                    udraw = rand(1);
                    eta1 = randn(2,1); % financial and macro shocks
                    
                    p12 = 1/(1+exp(samp.a12(dd)-samp.b12(dd)*(f_temp(1))+samp.c12(dd)*(m_temp(1))));
                    p21 = 1/(1+exp(samp.a21(dd)+samp.b21(dd)*(f_temp(1))-samp.c21(dd)*(m_temp(1))));
                    p11 = 1 - p12; % probability of remaining in normal
                    p22 = 1 - p21; % probability of remaining in bad
                    
                    if st_temp_par == 1 % started in good regime
                        if udraw > p11
                            s_temp(1) = 2; % switch from good to bad
                        else
                            s_temp(1) = 1; % don't switch and remain in good
                        end
                    else % start in bad regime
                        if udraw > p22
                            s_temp(1) = 1; % switch from bad to good
                        else
                            s_temp(1) = 2; % don't switch and remain in bad
                        end
                    end
                    
                    for ee=2:13 % Construct states, mf and ff for every period from t+1 to t+12, check state at each point in time and construct GDP accordingly, then avergage GDP (t+1 to t+12)
                        
                        eta1 = randn(2,1); % financial and macro shocks
                        
                        m_temp(ee) = samp.c_2_1(dd) +                                   samp.a1_2_1(dd)*f_temp(ee-1) + samp.a1_2_2(dd)*m_temp(ee-1) + samp.s_2_2(dd)*eta1(2,1);
                        f_temp(ee) = samp.c_1_1(dd) + samp.a0_1_2(dd)*m_temp(ee)*(-1) + samp.a1_1_1(dd)*f_temp(ee-1) + samp.a1_1_2(dd)*m_temp(ee-1) + samp.s_1_1(dd)*eta1(1,1);
                        
                        p12 = 1/(1+exp(samp.a12(dd)-samp.b12(dd)*(f_temp(ee))+samp.c12(dd)*(m_temp(ee))));
                        p21 = 1/(1+exp(samp.a21(dd)+samp.b21(dd)*(f_temp(ee))-samp.c21(dd)*(m_temp(ee))));
                        
                        p11 = 1 - p12; % probability of remaining in normal
                        p22 = 1 - p21; % probability of remaining in bad
                        
                        udraw = rand(1);
                        
                        if s_temp(ee-1) == 1 % started in good regime
                            if udraw > p11
                                s_temp(ee) = 2; % switch from good to bad
                            else
                                s_temp(ee) = 1; % don't switch and remain in good
                            end
                        else % start in bad regime
                            if udraw > p22
                                s_temp(ee) = 1; % switch from bad to good
                            else
                                s_temp(ee) = 2; % don't switch and remain in bad
                            end
                            
                        end
                        
                        eta2 = randn(1,1); % GDP shock
                        
                        if s_temp(ee)==2
                            y_temp(ee) = samp.c_3_1_sync_2(dd)+samp.a0_3_1_sync_2(dd)*f_temp(ee)*(-1) +samp.a0_3_2_sync_2(dd)*m_temp(ee)*(-1) +samp.s_3_3_sync_2(dd)*eta2;
                        else
                            y_temp(ee) = samp.c_3_1_sync_1(dd)+samp.a0_3_1_sync_1(dd)*f_temp(ee)*(-1) +samp.a0_3_2_sync_1(dd)*m_temp(ee)*(-1) +samp.s_3_3_sync_1(dd)*eta2;
                        end
                        
                        if ee==13
                            y_mat_fut_ex_temp(uu,dd,jj) = mean(y_temp(2:end,1));
                            st_it_temp(uu,dd,jj) = s_temp(ee);
                            st_it_temp_t(uu,dd,jj) = s_temp(1);
                        end
                    end
                end
            end
        end
        close(h);
        y_mat_fut_ex = NaN(nDraws*nParamDraws,size(y_mat,2)); % This stores all realizations for the ITERATED
        st_it = NaN(nDraws*nParamDraws,size(y_mat,2));
        st_it_t = NaN(nDraws*nParamDraws,size(y_mat,2));
        
        st_it_temp_OR = st_it_temp; st_it_temp_t_OR = st_it_temp_t; % Just storing the originals
        
        for jj=1:size(y_mat,2)
            for dd=1:nParamDraws
                for uu=1:nDraws
                    if st_it_temp(uu,dd,jj) ==1 % good regime
                        st_it_temp(uu,dd,jj) = 1;
                    elseif st_it_temp(uu,dd,jj) ==2 % bad regime
                        st_it_temp(uu,dd,jj) = 0;
                    end
                    if st_it_temp_t(uu,dd,jj) ==1 % good regime
                        st_it_temp_t(uu,dd,jj) = 1;
                    elseif st_it_temp_t(uu,dd,jj) ==2 % bad regime
                        st_it_temp_t(uu,dd,jj) = 0;
                    end
                end
            end
        end
        
        for rr=1:nDraws
            y_mat_fut_ex(rr*nParamDraws - nParamDraws+1:rr*nParamDraws,:) = squeeze(y_mat_fut_ex_temp(rr,:,:));
            st_it(rr*nParamDraws - nParamDraws+1:rr*nParamDraws,:) = squeeze(st_it_temp(rr,:,:));
            st_it_t(rr*nParamDraws - nParamDraws+1:rr*nParamDraws,:) = squeeze(st_it_temp_t(rr,:,:));
        end
        st_it_50 =     mean(st_it)';
        st_it_t_50 =     mean(st_it_t)';
        dY_25_fut_ex = prctile(y_mat_fut_ex,25)'; dY_75_fut_ex = prctile(y_mat_fut_ex,75)';
        dY_10_fut_ex = prctile(y_mat_fut_ex,10)'; dY_90_fut_ex = prctile(y_mat_fut_ex,90)';
    end
end