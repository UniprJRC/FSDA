function psix=TBpsix(u,c)
%TBpsix computes psi function (derivative of rho function) times x for Tukey's biweight  
%
%<a href="matlab: docsearchFS('TBpsix')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        tuning parameters. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...) 
%
%  Optional input arguments:
%
%  Output:
%
%
%   psix :     n x 1 vector which contains the values of TB psi(u)*u
%                function associated to the residuals or Mahalanobis
%                distances for the n units of the sample.
%
%
% More About:
%
%
% Function TBpsix transforms vector u as follows 
% 
%  \[
% TBpsix(u)= \left\{
%    \begin{array}{cc}
%   u^2[1-(u/c)^2]^2   if |u/c| \leq 1 \\
%  0                     &  |u/c|>1   \\
% \end{array}
%    \right.
%  \]
% See equation (2.38) p. 29 of Maronna et al. (2006)
%
% Remark: Tukey's biweight  psi-function is almost linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is approximately constant over the linear region of \psi,
% so the points in that region tend to get equal weight.
%
%
% See also: HApsix, HYPpsix, OPTpsix
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('TBpsix')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Plot of psi function (derivative of rho function) times x for Tukey's biweight. 
    x=-6:0.01:6;
    psixTB=TBpsix(x,2);
    plot(x,psixTB)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi (x)$','Interpreter','Latex')

%}

%% Beginning of code

psix = (abs(u) < c) .* (u.^2) .* ( 1 - (u./c).^2 ).^2 ;
end
%FScategory:UTISTAT