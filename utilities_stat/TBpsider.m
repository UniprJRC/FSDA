function psider=TBpsider(u,c)
%TBpsider computes derivative of psi function (second derivative of rho function) for Tukey's biweight  
%
%<a href="matlab: docsearchFS('HUpsider')">Link to the help function</a>
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
%
%  Optional input arguments:
%
%  Output:
%
%
%   psider :     derivative of psi function. Vector. 
%                n x 1 vector which contains the values of the derivative
%                of the Tukey biweight psi function associated to the
%                residuals or Mahalanobis distances for the n units of the
%                sample.
%
% More About:
%
% Function TBpsider transforms vector x as follows 
% \[
% TBpsider(x)=
% \left\{
%    \begin{array}{cc}
%  1- (x/c)^2 * [6- 5 (x/c)^2]                          &  \mbox{if  }  |x/c|<=1 
%    0            &                      |x/c|>1 \\       
% \end{array}
%    \right.
%  \] 

% Remark: Tukey's biweight  psi-function is almost linear around $u = 0$ in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  $\psi (u)/u$ is approximately constant over the linear region of $\psi$,
% so the points in that region tend to get equal weight.
%
% See also HUpsider, HYPpsider, HApsider, OPTpsider
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HUpsider')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{
    % Plot the derivative of Tukey's biweght psi function.
    x=-6:0.01:6;
    psiHUder=HUpsider(x,2);
    plot(x,psiHUder)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x)$','Interpreter','Latex')

%}

%% Beginning of code
psider = (abs(u) < c) .* (1 - 6*(u/c).^2  + 5*(u/c).^4);
end
%FScategory:UTISTAT