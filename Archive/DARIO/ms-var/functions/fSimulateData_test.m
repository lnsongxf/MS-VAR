function [f_draw,m_draw,y_draw,st_mat,p12_draw,p21_draw,eta] =  fSimulateData_test(dgp,opt)
if dgp.a0_1_2_sync_1 ~=0 && dgp.a0_2_1_sync_1 ~=0
    fprint('Simultaneous equations in FF and MF (regime 1)');
    return
end
if dgp.a0_1_2_sync_2 ~=0 && dgp.a0_2_1_sync_2 ~=0
    fprint('Simultaneous equations in FF and MF (regime 2)');
    return
end
rr = opt.rr; tt = opt.tt;

f_draw = zeros(rr,tt);m_draw = zeros(rr,tt);
p12_draw = NaN(rr,tt); p21_draw = NaN(rr,tt);
y_draw = zeros(rr,tt);

st_mat = NaN(rr,tt); st_mat(:,1) = ones(rr,1);
eta = NaN(3,rr,tt);

for dd=1:rr
    for jj=2:tt
        eta(:,dd,jj) = randn(3,1); % financial, macro and gdp shocks
        if st_mat(dd,jj-1)==2
            if dgp.a0_1_2_sync_2 ==0 % Then, MF does not cause FF, and FF should be drawn first
                f_draw(dd,jj) = dgp.c_1_1_sync_2                                   + dgp.a0_1_2_sync_2*m_draw(dd,jj) + dgp.a1_1_1_sync_2*f_draw(dd,jj-1) + dgp.a1_1_2_sync_2*m_draw(dd,jj-1) + dgp.a1_1_3_sync_2*y_draw(dd,jj-1) + dgp.s_1_1_sync_2*eta(1,dd,jj);
                m_draw(dd,jj) = dgp.c_2_1_sync_2 + dgp.a0_2_1_sync_2*f_draw(dd,jj)                                   + dgp.a1_2_1_sync_2*f_draw(dd,jj-1) + dgp.a1_2_2_sync_2*m_draw(dd,jj-1) + dgp.a1_2_3_sync_2*y_draw(dd,jj-1) + dgp.s_2_2_sync_2*eta(2,dd,jj);
            else                     % Then, FF does not cause MF, and MF should be drawn first
                m_draw(dd,jj) = dgp.c_2_1_sync_2 + dgp.a0_2_1_sync_2*f_draw(dd,jj)                                   + dgp.a1_2_1_sync_2*f_draw(dd,jj-1) + dgp.a1_2_2_sync_2*m_draw(dd,jj-1) + dgp.a1_2_3_sync_2*y_draw(dd,jj-1) + dgp.s_2_2_sync_2*eta(2,dd,jj);
                f_draw(dd,jj) = dgp.c_1_1_sync_2                                   + dgp.a0_1_2_sync_2*m_draw(dd,jj) + dgp.a1_1_1_sync_2*f_draw(dd,jj-1) + dgp.a1_1_2_sync_2*m_draw(dd,jj-1) + dgp.a1_1_3_sync_2*y_draw(dd,jj-1) + dgp.s_1_1_sync_2*eta(1,dd,jj);
            end
        else
            if dgp.a0_1_2_sync_1 ==0 % Then, MF does not cause FF, and FF should be drawn first
                f_draw(dd,jj) = dgp.c_1_1_sync_1                                   + dgp.a0_1_2_sync_1*m_draw(dd,jj) + dgp.a1_1_1_sync_1*f_draw(dd,jj-1) + dgp.a1_1_2_sync_1*m_draw(dd,jj-1) + dgp.a1_1_3_sync_1*y_draw(dd,jj-1) + dgp.s_1_1_sync_1*eta(1,dd,jj);
                m_draw(dd,jj) = dgp.c_2_1_sync_1 + dgp.a0_2_1_sync_1*f_draw(dd,jj)                                   + dgp.a1_2_1_sync_1*f_draw(dd,jj-1) + dgp.a1_2_2_sync_1*m_draw(dd,jj-1) + dgp.a1_2_3_sync_1*y_draw(dd,jj-1) + dgp.s_2_2_sync_1*eta(2,dd,jj);
            else                     % Then, FF does not cause MF, and MF should be drawn first
                m_draw(dd,jj) = dgp.c_2_1_sync_1 + dgp.a0_2_1_sync_1*f_draw(dd,jj)                                   + dgp.a1_2_1_sync_1*f_draw(dd,jj-1) + dgp.a1_2_2_sync_1*m_draw(dd,jj-1) + dgp.a1_2_3_sync_1*y_draw(dd,jj-1) + dgp.s_2_2_sync_1*eta(2,dd,jj);
                f_draw(dd,jj) = dgp.c_1_1_sync_1                                   + dgp.a0_1_2_sync_1*m_draw(dd,jj) + dgp.a1_1_1_sync_1*f_draw(dd,jj-1) + dgp.a1_1_2_sync_1*m_draw(dd,jj-1) + dgp.a1_1_3_sync_1*y_draw(dd,jj-1) + dgp.s_1_1_sync_1*eta(1,dd,jj);
            end
        end              
        [p12,p21] = fTranstionProb(dgp,m_draw(dd,jj),f_draw(dd,jj),opt);
        
        p11 = 1 - p12; % probability of remaining in normal
        p22 = 1 - p21; % probability of remaining in bad
        
        p12_draw(dd,jj) = p12; p21_draw(dd,jj) = p21;
        
        st_mat(dd,jj) = simulate_st(st_mat(dd,jj-1),p11,p22,1);
        if st_mat(dd,jj)==2
                y_draw(dd,jj) = dgp.c_3_1_sync_2 + dgp.a0_3_1_sync_2*f_draw(dd,jj) + dgp.a0_3_2_sync_2*m_draw(dd,jj) + dgp.a1_3_1_sync_2*f_draw(dd,jj-1) + dgp.a1_3_2_sync_2*m_draw(dd,jj-1) + dgp.a1_3_3_sync_2*y_draw(dd,jj-1) + dgp.s_3_3_sync_2*eta(3,dd,jj);
        else
                y_draw(dd,jj) = dgp.c_3_1_sync_1 + dgp.a0_3_1_sync_1*f_draw(dd,jj) + dgp.a0_3_2_sync_1*m_draw(dd,jj) + dgp.a1_3_1_sync_1*f_draw(dd,jj-1) + dgp.a1_3_2_sync_1*m_draw(dd,jj-1) + dgp.a1_3_3_sync_1*y_draw(dd,jj-1) + dgp.s_3_3_sync_1*eta(3,dd,jj);
        end              
    end
end
end