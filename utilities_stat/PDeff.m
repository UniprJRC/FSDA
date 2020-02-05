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
%    eff:       required efficienty. Scalar.
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
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
% Copyright 2008-2019.
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
    % Find c for a given efficiency.
    % The constant alpha associated to a nominal location efficiency of 95% in regression is
    % alpha= 0.224515798935881
    c=PDeff(0.95)
%}
%
%

%% Beginning of code

F=eff^2;
coeff=[F, 6*F, 15*F, 20*F-8, 15*F-12, 6*F-6, F-1];
r=roots(coeff);
r = r(imag(r)==0); % Save only the real roots
alpha = r(r>0);
end
%FScategory:UTISTAT
