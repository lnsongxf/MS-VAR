% Function to return model switches and model labels


function label = fSpecLabel(modelspec)
   
    if modelspec==1                    
            label = 'Ct2_A0t2_A1_SIGt2';
        elseif modelspec==2            
            label = 'CT_A0T_A1_SIGT';
        elseif modelspec==3
            label = 'CT_A0T_A1T_SIG';
        elseif modelspec==4
            label = 'CT_A0T_A1T_SIGT';
        elseif modelspec==5
            label = 'CT_A0_A1T_SIGT';
        elseif modelspec==6
            label = 'CT_A0T_A1_SIG';
        elseif modelspec==7
            label = 'Ct2_A0t2_A1t2_SIG';
        elseif modelspec==8
            label = 'Ct2_A0t2_A1t2_SIGt2';
        elseif modelspec==9
            label = 'Ct2_A0t2_A1t2_SIGt2_restr';
        elseif modelspec==10
            label = 'CT_A0T_A1_SIGT_restr_MF_FF';
        elseif modelspec==11
            label = 'CT_A0T_A1_SIGT_restr_GDP';
        elseif modelspec==12
            label = 'CT_A0T_A1T_SIGT_restr_GDP';
        elseif modelspec==13
            label = 'Ct2_A0t2_A1_SIGt2_restr_GDP';
    end
    

end