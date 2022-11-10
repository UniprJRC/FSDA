function alpha = PDbdp(bdp)
%PDbdp finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator
%
%
%<a href="matlab: docsearchFS('PDbdp')">Link to the help function</a>
%
%  Required input arguments:
%
%      bdp    : breakdown point. Scalar. Scalar defining breakdown point
%               (i.e a number between 0 and 0.5)
%               Data Types - single|double
%
%  Optional input arguments:
%
% Output:
%
%  alpha : Requested tuning constant. Scalar. Tuning constatnt of minimum
%           power divergence estimator associated to requested breakdown point
%
%
% See also: TBbdp, OPTbdp, HYPbdp, HAbdp
%
% References:
% 
%  Riani, M. Atkinson, A.C., Corbellini A. and Perrotta A. (2020), Robust
%  Regression with Density Power Divergence: Theory, Comparisons and Data
%  Analysis, submitted.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('PDbdp')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
%
% Examples:
%
%{
    % Find alpha given bdp.
    % The constant alpha associated to a breakdown point of 30% in regression is 
    % alpha=1.040816326530613
    alpha=PDbdp(0.3)
%}
%

%% Beginning of code
   
% Convergence condition is E(\rho) = sup(\rho(x)) bdp
%  where sup(\rho(x)) is 1
alpha=1./((1-bdp).^2)-1;

end
%FScategory:UTISTAT