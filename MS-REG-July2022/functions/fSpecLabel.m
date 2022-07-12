% Function to return model switches and model labels


function label = fSpecLabel(modelspec)

    % Model specifications as described in SZ-VAR.pdf
    if modelspec==1                    
            label = 'Ct2_A0t2_A1_SIGt2';
        elseif modelspec==2            
            label = 'CT_A0T_A1_SIGT';
        elseif modelspec==3
            label = 'CT_A0T_A1T_SIG';
        elseif modelspec==4
            label = 'CT_A0T_A1T_SIGT';
        elseif modelspec==5  % equivalent to model 105
            label = 'CT_A0_A1T_SIGT';
        elseif modelspec==6  % equivalent to model 102
            label = 'CT_A0T_A1_SIG';
        elseif modelspec==7
            label = 'Ct2_A0t2_A1t2_SIG';
        elseif modelspec==8
            label = 'Ct2_A0t2_A1t2_SIGt2';
        elseif modelspec==9
            label = 'Ct2_A0t2_A1t2_SIGt2_restr';
        elseif modelspec==10
            label = 'CT_A0T_A1_SIGT_restr_MF_FF';
        elseif modelspec==11 % equivalent to model 104
            label = 'CT_A0T_A1_SIGT_restr_GDP';
        elseif modelspec==12 % equivalent to model 106
            label = 'CT_A0T_A1T_SIGT_restr_GDP';
        elseif modelspec==13
            label = 'Ct2_A0t2_A1_SIGt2_restr_GDP';
        % Models Used in the Paper
        elseif modelspec==101
            label = 'CT_A0_A1_SIG_restr_GDP';
        elseif modelspec==102
            label = 'CT_A0_A1_SIGT_restr_GDP';
        elseif modelspec==103
            label = 'CT_A0T_A1_SIG_restr_GDP';
        elseif modelspec==104
            label = 'CT_A0_A1T_SIG_restr_GDP';
        elseif modelspec==105
            label = 'CT_A0T_A1_SIGT_restr_GDP';
        elseif modelspec==106
            label = 'CT_A0_A1T_SIGT_restr_GDP';
        elseif modelspec==107
            label = 'CT_A0T_A1T_SIGT_restr_GDP';
        elseif modelspec==108                    
            label = 'Ct2_A0_A1_restr_SIGt2';
        elseif modelspec==201
            label = 'CT_A0_A1_SIG_GDP';
        elseif modelspec==202
            label = 'CT_A0_A1_SIGT_GDP';
        elseif modelspec==203
            label = 'CT_A0T_A1_SIG_GDP';
        elseif modelspec==204
            label = 'CT_A0_A1T_SIG_GDP';
        elseif modelspec==205
            label = 'CT_A0T_A1_SIGT_GDP';
        elseif modelspec==206
            label = 'CT_A0_A1T_SIGT_GDP';
        elseif modelspec==207
            label = 'CT_A0T_A1T_SIGT_GDP';
        elseif modelspec==208
            label = 'Ct2_A0_A1_SIGt2_GDP';
    end
end