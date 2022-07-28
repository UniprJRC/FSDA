function rhoTB = TBrho(u,c)
%TBrho computes rho function for Tukey's biweight
%
%<a href="matlab: docsearchFS('TBrho')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. 
%               Vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        tuning parameter. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...) 
%
%  Optional input arguments:
%
%
%  Output:
%
%
%   rhoTB :      vector of length n which contains the Tukey's biweight rho
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
%  (c^2/6) \left\{ 1-[1-(u/c)^2]^3 \right\}  &  |u/c| \leq 1  \\
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
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
% Riani, M., Cerioli, A. and Torti, F. (2014), On consistency factors and
% efficiency of robust S-estimators, "TEST", Vol. 23, pp. 356-387.
% http://dx.doi.org/10.1007/s11749-014-0357-7
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('TBrho')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Plot of rho function.
    close all
    x=-6:0.01:6;
    rhoTB=TBrho(x,2);
    plot(x,rhoTB,'LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel('$\rho (u,2)$','Interpreter','Latex')
    text(x(1)-0.8,rhoTB(1),'c^2/6')
    text(x(end)+0.2,rhoTB(end),'c^2/6')
    title('$\rho (u,c)$','Interpreter','Latex')
    hold('on')
    c=2;
    stem(c,c^2/6,'LineStyle',':','LineWidth',1)
    stem(-c,c^2/6,'LineStyle',':','LineWidth',1)

%}

%{
    %% Compare two rho functions for 2 different values of c.  
    % In the first we fix the bdp (value of efficiency is automatically given),
    % while in the second we find the efficiency (the value of bdp is
    % automatically given)
    close all
    x=-6:0.01:6;
    lwd=2;
    hold('on')
    c=TBbdp(0.5,1);
    rhoTB=TBrho(x,c);
    rhoTB=rhoTB/max(rhoTB);
    plot(x,rhoTB,'LineStyle','-','LineWidth',lwd)
    c=TBeff(0.95,1);
    rhoTB=TBrho(x,c);
    rhoTB=rhoTB/max(rhoTB);
    plot(x,rhoTB,'LineStyle','-.','LineWidth',lwd)
    xlabel('$x$','Interpreter','Latex','FontSize',16)
    ylabel('TB. Normalized $\rho_c(x)$','Interpreter','Latex','FontSize',20)
    legend({'$c_{(bdp=0.5 \mapsto eff=0.29)}$', '$c_{(eff=0.95 \mapsto bdp=0.12)}$'},'Interpreter','Latex','Location','SouthEast','FontSize',16)
%}

%% Beginning of code
c=c(1); % MATLAB Ccoder instruction to enforce that c is a scalar
w = (abs(u)<=c);
rhoTB = (u.^2/(2).*(1-(u.^2/(c^2))+(u.^4/(3*c^4)))).*w +(1-w)*(c^2/6);

end
%FScategory:UTISTAT