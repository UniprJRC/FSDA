function rhoTB = TBrho(u,c)
%TBrho computes (rho) function for Tukey biweight
%
%<a href="matlab: docsearchFS('tbrho')">Link to the help function</a>
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
%
%  Output:
%
%
%   rhoTB :      n x 1 vector which contains the Tukey's biweight rho
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
% More About:
%
%
% function TBrho transforms vector u as follows 
% \[
% TBrho(u)= \left\{
%    \begin{array}{cc}
%  (c^2/6)*\left\{ 1-[1-(u/c)^2]^3 \right\}  &  |u/c| \leq 1  \\
%  (c^2/6)                      &  |u/c| >1   \\
% \end{array}
%    \right.
%  \]
%  
% See equation (2.37) p. 29 of Maronna et al. (2006).
% Remark: equation (2.37) is written in standardized terms in such a way
% that $\rho(c)=1$, so it is the previous expression divided by $(c^2/6)$
%
% See also HYPrho, HArho, OPTrho
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
% Riani M., Cerioli A., Torti F. (2014). On consistency factors and
% efficiency of robust S-estimators, TEST, Volume 23, Issue 2, pp 356-387.
% DOI: 10.1007/s11749-014-0357-7
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('tbrho')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    x=-6:0.01:6;
    rhoTB=TBrho(x,2);
    plot(x,rhoTB)
    xlabel('x','Interpreter','Latex')
    ylabel('$\rho (x,2)$','Interpreter','Latex')
    text(x(1)-0.8,rhoTB(1),'c^2/6')
    text(x(end)+0.2,rhoTB(end),'c^2/6')
    title('$\rho (x,c)$','Interpreter','Latex')
%}

%% Beginning of code

w = (abs(u)<=c);
rhoTB = (u.^2/(2).*(1-(u.^2/(c^2))+(u.^4/(3*c^4)))).*w +(1-w)*(c^2/6);

end
