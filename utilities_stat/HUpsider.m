function psider=HUpsider(u,c)
%HUpsider computes derivative of psi function (second derivative of rho function) for Huber 
%
%<a href="matlab: docsearchFS('HUpsider')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1 vector containing 
%               for the n units of the sample scaled residuals or
%               Mahalanobis distances
%    c :        tuning constant. Scalar. scalar greater than 0 which controls the robustness/efficiency of the estimator 
%               (beta in regression or mu in the location case ...) 
%
% Optional input arguments:
%
%  Output: 
%
%   psider :     derivative of psi function. Vector. 
%                n x 1 vector which contains the values of the derivative
%                of the Huber psi function associated to the
%                residuals or Mahalanobis distances for the n units of the
%                sample.
%
%
% More About:
%
% Function HUpsider transforms vector x as follows 
% HUpsider(x)=
%
% \[
% HUpsider(u)= \left\{
%    \begin{array}{cc}
%  1                            &  \mbox{if  }  |u/c| \leq 1 \\
%  0                            &  |u/c|>1   \\
% \end{array}
%    \right.
%  \]
%
%
% See also: TBpsider, HYPpsider, OPTpsider
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HUpsider')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Plot of derivative of psi function.
    x=-6:0.01:6;
    psiHUder=HUpsider(x,2);
    plot(x,psiHUder)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x)$','Interpreter','Latex')

%}

%% Beginning of code
psider=zeros(length(u),1);
psider(abs(u) < c)=1;
end
%FScategory:UTISTAT