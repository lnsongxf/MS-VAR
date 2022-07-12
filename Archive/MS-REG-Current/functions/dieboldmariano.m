function [tstat, nanny] = dieboldmariano(error1, error2, nwlag)
% DIEBOLDMARIANO ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 09-May-2015 17:59:23 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 8.1.0.604 (R2013a) 
% FILENAME  : dieboldmariano.m 


if nargin < 3
    nwlag = 4;
end

deltaloss = error1 - error2;
nanny     = ~isnan(deltaloss);
deltaloss = deltaloss(nanny);
Nobs      = length(deltaloss);
reggae    = nwest(deltaloss, ones(Nobs,1), nwlag);

tstat = reggae.tstat(1);
