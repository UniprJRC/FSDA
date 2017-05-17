function psix=HUpsix(u,c)
%HUpsix computes psi function (derivative of rho function) times x for Huber 
%
%<a href="matlab: docsearchFS('HUpsix')">Link to the help function</a>
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
%   psix :     psi(u) function multiplied by u. Vector. n-by-1 vector which contains the values of HU psi(u)*u
%                function associated to the residuals or Mahalanobis
%                distances for the n units of the sample.
%
%
% More About:
%
%
% Function HUpsix transforms vector u as follows 
% 
% \[
% HUpsix(u)= \left\{
%    \begin{array}{cc}
%  u^2                            &  \mbox{if  }  |u/c| \leq 1 \\
%  c \times \mbox{sign}(u) u      &  |u/c|>1   \\
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
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HUpsix')">Link to the help page for this function</a>
% Last modified 31-05-2016

% Examples:

%{

    % Plot of psi(x) function multiplied  by x.
    x=-6:0.01:6;
    psixHU=HUpsix(x,2);
    plot(x,psixHU)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi (x)$','Interpreter','Latex')

%}

%% Beginning of code

w = (abs(u)<=c);

psix=(u.^2).*w +(1-w).*(c*sign(u).*u);
end
%FScategory:UTISTAT