function psiAS=ASpsi(u,c)
%ASpsi computes psi function (derivative of rho function) for Andrew's sine function  
%
%<a href="matlab: docsearchFS('ASpsi')">Link to the help function</a>
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
%   psiAS :      n x 1 vector which contains the Tukey's psi
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
% More About:
%
% Function ASpsi transforms vector u as follows 
%
% \[
% ASpsi(u)= \left\{
%    \begin{array}{cc}
%   \sin(u/c)             &   |u/c| \leq \pi \\
%  0                     &    |u/c|> \pi   \\
% \end{array}
%    \right.
%  \]
%
%
% Remark: Andrews's psi-function is almost linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is approximately constant over the linear region of \psi,
% so the points in that region tend to get equal weight.
%
% See also HYPpsi, HApsi, OPTpsi, PDpsi
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
%<a href="matlab: docsearchFS('ASpsi')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Plot of psi function for Andrew's sine function.
    close all
    x=-6:0.01:6;
    c=1.5;
    psiAS=ASpsi(x,c);
    plot(x,psiAS,'LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel('$\psi(u,2)$','Interpreter','Latex','FontSize',14)
    hold('on')
    ax=axis;
    line([-c*pi;-c*pi],[ax(3);0],'LineStyle',':','LineWidth',1)
    line([c*pi;c*pi],[ax(3);0],'LineStyle',':','LineWidth',1)
%}

%% Beginning of code
psiAS = zeros(size(u));
inds1 = abs(u) < c*pi;
psiAS(inds1) = (sin(u(inds1)/c));
end
%FScategory:UTISTAT