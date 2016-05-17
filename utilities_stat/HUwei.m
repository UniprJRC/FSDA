function wei = HUwei(u,c)
%HUwei computes weight function psi(u)/u for Huber   
%
%<a href="matlab: docsearchFS('HUwei')">Link to the help function</a>
%
%
% Required input arguments:
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
%    wei :       Weights vector. Vector. 
%                n x 1 vector containing the Huber weights
%                associated to the scaled residuals or Mahalanobis
%                distances for the n units of the sample.
%
% More About:
%
% Function HUwei transforms vector u as follows 
% \[
% HUpsi(u)= \left\{
%    \begin{array}{cc}
%  1                            &  \mbox{if  }  |u/c| \leq 1 \\
%  c \times \frac{\mbox{sign}(u)}{u}      &  |u/c|>1   \\
% \end{array}
%    \right.
%  \]
%
% See equation (2.32) p. 28 of Maronna et al. (2006)
%
% See also TBwei, HYPwei, HAwei, OPTwei
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
% Riani M., Cerioli A., Torti F. (2014). On consistency factors and
% efficiency of robust S-estimators, TEST, Volume 23, Issue 2, pp 356-387.
% DOI: 10.1007/s11749-014-0357-7
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HUwei')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % Plot Huber weights.
    x=-6:0.01:6;
    weiHU=HUwei(x,2);
    plot(x,weiHU)
    xlabel('x','Interpreter','Latex')
    ylabel('$W (x) =\psi(x)/x$','Interpreter','Latex')

%}

%% beginning of code
w = (abs(u)<=c);

wei=w +(1-w).*(c*sign(u)./u);


end
%FScategory:UTISTAT