function alpha = PDeff(eff)
%PDeff finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator
%
%
%
%<a href="matlab: docsearchFS('PDeff')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    eff:       required efficiency. Scalar.
%               Scalar which contains the required efficiency (of location
%               or scale estimator).
%               Generally eff=0.85, 0.9 or 0.95
%               Data Types - single|double
%
%  Optional input arguments:
%
%
%
% Output:
%
%  alpha : Requested tuning constant. Scalar. Tuning constant for minimum
%           power divergence estimator associated to requested value of
%           efficiency
%
% See also: TBeff, OPTeff, HYPeff, HAeff
%
% References:
%
%  Riani, M. Atkinson, A.C., Corbellini A. and Perrotta A. (2020), Robust
%  Regression with Density Power Divergence: Theory, Comparisons and Data
%  Analysis, Entropy, Vol. 22, 399. 
%  https://www.mdpi.com/1099-4300/22/4/399 
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('PDeff')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

%
% Examples:
%
%{
    % Find tuning constant alpha for a given efficiency.
    % The constant alpha associated to a nominal location efficiency of 95% in regression is
    % alpha= 0.224515798935881
    c=PDeff(0.95);
%}
%
%

%% Beginning of code

F=eff.^(2/3);
alpha=(1-F+sqrt(1-F))./F;

% Alternative more complicated way to find alpha
% F=eff^2;
% coeff=[F, 6*F, 15*F, 20*F-8, 15*F-12, 6*F-6, F-1];
% r=roots(coeff);
% r = r(imag(r)==0); % Save only the real roots
% alpha = r(r>0);
end
%FScategory:UTISTAT
