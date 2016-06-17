function psi=HUpsi(u,c)
%HUpsi computes psi function (derivative of rho function) for Huber
%
%<a href="matlab: docsearchFS('HUpsi')">Link to the help function</a>
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
%   HUpsi :      n x 1 vector which contains the Huber's psi
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
% More About:
%
% Function HUpsi transforms vector u as follows
%
% \[
% HUpsi(u)= \left\{
%    \begin{array}{cc}
%  u                            &  \mbox{if  }  |u/c| \leq 1 \\
%  c \times \mbox{sign}(u)      &  |u/c|>1   \\
% \end{array}
%    \right.
%  \]
% See equation (2.38) p. 29 of Maronna et al. (2006)
%
% Remark: Tukey's biweight  psi-function is linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is constant over the linear region of \psi,
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
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HUpsi')">Link to the help page for this function</a>
% Last modified 31-05-2016

% Examples:

%{
    close all
    x=-6:0.01:6;
    c=1.345;
    psiHU=HUpsi(x,c);
    plot(x,psiHU,'LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel('$\psi (u,1.345)$','Interpreter','Latex','FontSize',14)
    text(-c,-c,'-c=-1.345')
    text(c,c+0.1,'c=1.345')
    hold('on')
    stem(c,c,'LineStyle',':','LineWidth',2)
    stem(-c,-c,'LineStyle',':','LineWidth',2)

%}

%% Beginning of code
w = (abs(u)<=c);

psi=u.*w +(1-w).*(c*sign(u));
end
%FScategory:UTISTAT