function psiderAS=ASpsider(u,c)
%ASpsider computes derivative of psi function (second derivative of rho function) or Andrew's sine function  
%
%<a href="matlab: docsearchFS('ASpsider')">Link to the help function</a>
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
%   psiderAS :     derivative of psi function. Vector. 
%                n x 1 vector which contains the values of the derivative
%                of the Andrew's psi function associated to the
%                residuals or Mahalanobis distances for the n units of the
%                sample.
%
% More About:
%
% Function ASpsider transforms vector x as follows 
% \[
% ASpsi(u)= \left\{
%    \begin{array}{cc}
%   (1/c) cos(u/c)      & |u/c| \leq \pi \\
%  0                     &  |u/c|> \pi   \\
% \end{array}
%    \right.
%  \]
%
% Remark: Andrew's sine functionon is almost linear around $u = 0$ in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  $\psi (u)/u$ is approximately constant over the linear region of $\psi$,
% so the points in that region tend to get equal weight.
%
% See also HUpsider, HYPpsider, HApsider, OPTpsider
%
%
% References:
%
% Andrews, D.F., Bickel, P.J., Hampel, F.R., Huber, P.J., Rogers, W.H., and
% Tukey, J.W. (1972), "Robust Estimates of Location: Survey and Advances",
% Princeton Univ. Press, Princeton, NJ. [p. 203]
%
% Andrews, D. F. (1974). A Robust Method for Multiple Linear Regression,
% "Technometrics", V. 16, pp. 523-531, https://doi.org/10.1080/00401706.1974.10489233
%
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ASpsider')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %% Plot the derivative of Andrew's psi function (when bdp=0.5).
    x=-6:0.01:6;
    c=ASbdp(0.5,1);
    psiASder=ASpsider(x,c);
    plot(x,psiASder)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x)$','Interpreter','Latex')
%}

%% Beginning of code
psiderAS = zeros(size(u));
inds1 = abs(u) < c*pi;
psiderAS(inds1) = (1/c)*(cos(u(inds1)/c));
end
%FScategory:UTISTAT