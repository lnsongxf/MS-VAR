function [restrictions,markov_chains,switch_priors]=create_restrictions_and_markov_chains0(nlags)
% create_restrictions_and_markov_chains0 -- creates restrictions for
% the constant-parameter SVAR model
%
% ::
%
%
%   [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains0(tpl)
%
% Args:
%
%    - **markov_chains** [empty|struct]: structure of previously defined
%    markov chains
%
% Returns:
%    :
%
%    - **lin_restr** [cell]: cell array of restrictions (see below).
%
%    - **nonlin_restr** [cell]: cell array of inequality restrictions
%
%    - **markov_chains** [struct]: modified markov chains
%
% Note:
%
%    - The syntax to construct a restriction
%      --> ai(eqtn)
%      --> ai(eqtn,vbl)
%      --> ai(eqtn,vbl,chain_name,state)
%      --> a(eqtn)
%      --> a(eqtn,vbl)
%      --> a(eqtn,vbl,chain_name,state)
%      - **eqtn** [integer]: integer
%      - **vbl** [integer|char]: integer or variable name
%      - **i** [integer]: lag
%      - **chain_name** [char]: name of the markov chain
%      - **state** [integer]: state number
%
%    - The lag coefficients are labelled a0, a1, a2,...,ak, for a model with k
%    lags. Obviously, a0 denotes the contemporaneous coefficients.
%
%    - The constant terms labelled c_1_1, c_2_2,...,c_n_n, for a model with n
%    endogenous variables.
%
%    - The standard deviations labelled s_1_1, s_2_2,...,s_n_n, for a
%    model with n endogenous variables.
%
% Example:
%
%    See also:

markov_chains=struct('name',{},...
'number_of_states',{},...
'controlled_parameters',{},...
'endogenous_probabilities',{},...
'probability_parameters',{});

switch_priors=struct();
    

% syntax is alag(eqtn,vname)
%-------------------------------
% varlist = {'EFFR','LIPM','UNRATE','LPPI','BAA10YMOODY','MHF'};

% create restriction on Atilde_+ matrix
% first N equations, "MHF" column
% this is the (1,2) block of Atilde_+
%----------------------------------
aplus_restr=cell(0,1);
for ilag=1:nlags
    mylag=int2str(ilag);
    aplus_restr=[aplus_restr
         {
        ['a',mylag,'(6,EFFR)=a',mylag,'(1,EFFR)*(-0.0189/0.0231)']
        ['a',mylag,'(6,LIPM)=a',mylag,'(1,LIPM)*(-0.0189/0.0231)']
        ['a',mylag,'(6,UNRATE)=a',mylag,'(1,UNRATE)*(-0.0189/0.0231)']
        ['a',mylag,'(6,LPPI)=a',mylag,'(1,LPPI)*(-0.0189/0.0231)']
        ['a',mylag,'(6,BAA10YMOODY)=a',mylag,'(1,BAA10YMOODY)*(-0.0189/0.0231)']
        }];
end
   
a0_restr={
    % high relevance prior: sigma_nu = 0.5*std(db.MHF) = 0.0231
    % 1/sigma_nu = 10.8033
    % posterior estimate for rho = b^2/(b^2+sigma_nu^2) under this prior is 0.4
    % and thus beta = sqrt(sigma_nu^2*rho/(1-rho)) = 0.0189
    % resulting IRFs are in Figure 5 of original paper
    's(6,6)=0.0231'
    % sixth equation or "MHF" equation
    % this is the (2,1) block of Atilde_0
    %----------------------------------
    'a0(6,EFFR)=a0(1,EFFR)*(-0.0189/0.0231)'
    'a0(6,LIPM)=a0(1,LIPM)*(-0.0189/0.0231)'
    'a0(6,UNRATE)=a0(1,UNRATE)*(-0.0189/0.0231)'
    'a0(6,LPPI)=a0(1,LPPI)*(-0.0189/0.0231)'
    'a0(6,BAA10YMOODY)=a0(1,BAA10YMOODY)*(-0.0189/0.0231)'
    % first N equations, "MHF" column
    % this is the (1,2) block of Atilde_0
    %----------------------------------
    'a0(1,MHF)=0'
    'a0(2,MHF)=0'
    'a0(3,MHF)=0'
    'a0(4,MHF)=0'
    'a0(5,MHF)=0'
    };
lin_restr = [a0_restr;aplus_restr];

nonlin_restr={
%     's(1,1)<1/s(6,6)'
    };

restrictions=[nonlin_restr;lin_restr];

end
