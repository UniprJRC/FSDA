function psiOPT=OPTpsi(u,c)
%OPTpsi computes psi function (derivative of rho function) for optimal weight function
%
%<a href="matlab: docsearchFS('OPTpsi')">Link to the help function</a>
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
%   psiOPT :     n x 1 vector which contains the values of optimal psi
%                function associated to the residuals or Mahalanobis
%                distances for the n units of the sample.
%
% More About:
%
% function OPTpsi transforms vector u as follows
%
%
% \[
%   OPTpsi(u,c) = \left\{
%    \begin{array}{cc}
%                  \frac{u}{3.25*c^2}          &                                               |u| \leq 2c   \\
%                
%     =   (1/3.25) \left( -1.944  \frac{u}{c^2} + 1.728  \frac{u^3}{c^4} - 0.312 \frac{u^5}{c^6} + 0.016 \frac{u^7}{c^8}   \right)  &  \qquad 2c \leq |u| \leq 3c \\
%                
%                   0                                    &                                  |u|>3c \\
% \end{array}
%    \right.
%  \]
% 
% Remark: Optimal psi-function is almost linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  $\psi(u)/u$ is approximately constant over the linear region of $\psi$,
% so the points in that region tend to get equal weight.
%
% See also: HYPpsi, HApsi, TBpsi
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
%<a href="matlab: docsearchFS('OPTpsi')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Plot of psi function (derivative of rho function) for optimal weight
    % function.
    x=-6:0.01:6;
    psiOPT=OPTpsi(x,1.2);
    plot(x,psiOPT)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi (x)$','Interpreter','Latex')

%}

%% Beginning of code

% Computes Standardized optimal psi function (first derivative of rho function)
% \rho'(x)

psiOPT = zeros(length(u),1);
absx=abs(u);

% r /(3.25c^2) if r <=2*c
inds1 = absx <= 2*c;
psiOPT(inds1) = u(inds1) / (3.25*c^2);

% 1/(3.25) * ( -1.944* (r/c)^2 .... +8*0.002 (r/c)^8 )    if    2c< r <3c
inds1=(absx > 2*c)&(absx <= 3*c);
x1 = u(inds1);
psiOPT(inds1) = (-1.944 * x1 / c^2 + 1.728 * x1.^3 / c^4 - 0.312 * x1.^5 / c^6 + 0.016 * x1.^7 / c^8) / 3.25;

% 0 if r >3*c

end
%FScategory:UTISTAT