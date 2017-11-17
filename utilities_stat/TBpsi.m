function psi=TBpsi(u,c)
%TBpsi computes psi function (derivative of rho function) for Tukey's biweight  
%
%<a href="matlab: docsearchFS('TBpsi')">Link to the help function</a>
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
%   tbpsi :      n x 1 vector which contains the Tukey's psi
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
% More About:
%
% Function TBpsi transforms vector u as follows 
%
% \[
% TBpsi(u)= \left\{
%    \begin{array}{cc}
%  (c^2/6) u[1-(u/c)^2]^2    \mbox{if} |u/c| \leq 1 \\
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
% See also HYPpsi, HApsi, OPTpsi
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
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('TBpsi')">Link to the help page for this function</a>
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Plot of psi function.
    close all
    x=-6:0.01:6;
    c=2;
    psiTB=TBpsi(x,c);
    plot(x,psiTB,'LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel('$\psi(u,2)$','Interpreter','Latex','FontSize',14)
    hold('on')
    ax=axis;
    line([-c;-c],[ax(3);0],'LineStyle',':','LineWidth',1)
    line([c;c],[ax(3);0],'LineStyle',':','LineWidth',1)
%}

%% Beginning of code

psi = (abs(u) < c) .* u .* ( 1 - (u./c).^2 ).^2 ;
end
%FScategory:UTISTAT