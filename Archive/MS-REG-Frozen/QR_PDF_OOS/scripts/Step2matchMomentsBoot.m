function Res = Step2matchMomentsBoot(YQ,QQ,YY,delta,Qmatch)
%% Step2match: Match skewed t-distributions to quantiles

% Input arguments:
% - YQ : Matrix (with number of columns equal to length(QQ)) containing the
%        conditional quantiles to which the skewed t-distributions will be
%        fitted.
% - QQ : Vector of numbers between 0 and 1 (exclusive) indicating the
%        quantiles to which the columns of YQ and YQunc correspond. QQ
%        should include the quantiles listed above that will be used to
%        match the skewed t-distributions.
%
% Output arguments:
% - Res : struct containing the following fields (for convenience let T
%         denote the number of rows of YQ).
%   - PST : T-by-length(YY) matrix where PST(t,:) contains the values of
%           the probability density function of the skewed t-distribution
%           fit to the conditional quantiles YQ(t,:), evaluated at each
%           point in YY.
%   - MeanST : Vector with length T containing the mean of the skewed
%              t-distribution fit to the conditional quantiles YQ(t,:).
%   - VarianceST : Vector with length T containing the variance of the
%                  skewed t-distribution fit to the conditional quantiles
%                  YQ(t,:).
%   - SkewnessST : Vector with length T containing the skewness of the
%                  skewed t-distribution fit to the conditional quantiles
%                  YQ(t,:).
%   - KurtosisST : Vector with length T containing the kurtosis of the
%                  skewed t-distribution fit to the conditional quantiles
%                  YQ(t,:).
%   - YY : Vector containing the grid of values for which the skewed
%          t-distribution PDF and CDF are evaluated.


%% Fit skewed t-distribution to conditional quantiles
% Initialize matrices to store output
PST            = NaN(size(YQ,1),length(YY));
MeanST_store  = NaN(size(YQ,1),1);
VarianceST     = NaN(size(YQ,1),1);
SkewnessST     = NaN(size(YQ,1),1);
KurtosisST     = NaN(size(YQ,1),1);


% Fit skewed t-distribution for each observation
parfor jt = 1:size(YQ,1)

    disp(['Now matching quantile function for observation ',...
    num2str(jt),' (out of a total of ',num2str(size(YQ,1)),')'])

    % Target quantiles, for observation jt
    qqTarg = YQ(jt,:);
    
    % Display progress message
    if mod(jt,12)==0
%         disp(['Now matching quantile function for observation ',...
%               num2str(jt),' (out of a total of ',num2str(size(YQ,1)),')'])
    end
    
    if ~any(isnan(qqTarg))
        
        %%% NOTE: here I no longer use previously estimated parameters as
        %%%       initial conditions
        %if jt<=5 || mod(jt,4)==0
            [lc,sc,sh,df]  = QuantilesInterpolation(qqTarg,QQ,Qmatch);
        %else
        %    [lc,sc,sh,df]  = QuantilesInterpolation(qqTarg,QQ,lc0,sc0,sh0,df0);
        %end
        
        %lc0=lc; sc0=sc; sh0=sh; df0=df;
        %%% END NOTE SECTION
        
        % Compute PDF and CDF for fitted distribution
        PST(jt,:)   = dskt(YY,lc,sc,sh,df);
        
        % Compute moments for fitted distribution
        MeanST = sum(YY .* PST(jt,:) * delta);
        MeanST_store(jt,:) = MeanST;
        VarianceST(jt,:) = sum(((YY - MeanST).^2) .* PST(jt,:) * delta);
        SkewnessST(jt,:) = sum(((YY - MeanST).^3) .* PST(jt,:) * delta) / (VarianceST(jt,:)^(3/2));
        KurtosisST(jt,:) = sum(((YY - MeanST).^4) .* PST(jt,:) * delta) / (VarianceST(jt,:)^2);
        
    end
end

%% Store all results in struct Res
Res.PST            = PST;
Res.MeanST         = MeanST_store;
Res.VarianceST     = VarianceST;
Res.SkewnessST     = SkewnessST;
Res.KurtosisST     = KurtosisST;
Res.YY             = YY;
