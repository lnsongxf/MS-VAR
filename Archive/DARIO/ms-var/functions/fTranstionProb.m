function [p12,p21] = fTranstionProb(param,MF,FF,opt)

ww = length(MF); p12 = NaN(ww,1); p21 = NaN(ww,1);

a12=param.a12;
a21=param.a21;
b12=param.b12;
b21=param.b21;
c12=param.c12;
c21=param.c21;

for jj=1:ww
    if opt.const==1
        if opt.normal==1
            p12(jj) = 1./(1+exp(a12-b12*FF(jj)-c12*MF(jj)));
            p21(jj) = 1./(1+exp(a21-b21*FF(jj)-c21*MF(jj)));
        elseif opt.normal==0
            p12(jj) = 1./(1+exp(a12-b12*FF(jj)+c12*MF(jj)));
            p21(jj) = 1./(1+exp(a21+b21*FF(jj)-c21*MF(jj)));
        end
    elseif opt.const==0
        if opt.normal==1
            p12(jj) = 1./(1+exp(-b12*FF(jj)-c12*MF(jj)));
            p21(jj) = 1./(1+exp(-b21*FF(jj)-c21*MF(jj)));
        elseif opt.normal==0
            p12(jj) = 1./(1+exp(-b12*FF(jj)+c12*MF(jj)));
            p21(jj) = 1./(1+exp(+b21*FF(jj)-c21*MF(jj)));
        end
    end
end

end

