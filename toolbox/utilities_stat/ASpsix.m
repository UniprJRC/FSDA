function psix=ASpsix(u,c)
%ASpsix computes psi function (derivative of rho function) times x for Andrew's sine function   
%
%<a href="matlab: docsearchFS('ASpsix')">Link to the help function</a>
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
%   psix :     n x 1 vector which contains the values of AS psi(u)*u
%                function associated to the residuals or Mahalanobis
%                distances for the n units of the sample.
%
%
% More About:
%
%
% Function ASpsix transforms vector u as follows 
% 
% \[
% ASpsix(u)= \left\{
%    \begin{array}{cc}
%   (u/c) sin(u/c)      & |u/c| \leq \pi \\
%  0                     &  |u/c|> \pi   \\
% \end{array}
%    \right.
%  \]
%
% Remark: Andrew's  psi-function is almost linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is approximately constant over the linear region of \psi,
% so the points in that region tend to get equal weight.
%
%
% See also: HApsix, HYPpsix, OPTpsix, TBpsix, PDpsix
%
% References:
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
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ASpsix')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Plot of psi function (derivative of rho function) times x for Andrews weight function. 
    x=-6:0.01:6;
    psixAS=ASpsix(x,2);
    plot(x,psixAS)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi (x)$','Interpreter','Latex')

%}

%% Beginning of code
psix = zeros(size(u));
inds1 = abs(u) < c*pi;
psix(inds1) = (u(inds1)/c).*(sin(u(inds1)/c));

end
%FScategory:UTISTAT