function rhoHU = HUrho(u,c)
%HUrho computes rho function for Huber
%
%<a href="matlab: docsearchFS('HUrho')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
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
%   rhoHU :      n x 1 vector which contains the Huber rho
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
% More About:
%
%
% function HUrho transforms vector u as follows 
% \[
% HUrho(u)= \left\{
%    \begin{array}{cc}
%  (u^2/2)    &  |u/c| \leq 1  \\
%  c|u| -c^2/2               &  |u/c| >1   \\
% \end{array}
%    \right.
%  \]
%  
% See equation (2.27) p. 26 of Maronna et al. (2006).
% Remark: equation (2.26) is written in standardized terms in such a way
% that $\rho(c)=1$, so it is the previous expression multiplied by $2$
%
% See also TBrho, HYPrho, HArho, OPTrho
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
% Riani, M., Cerioli, A. and Torti, F. (2014), On consistency factors and
% efficiency of robust S-estimators, "TEST", Vol. 23, pp. 356-387.
% http://dx.doi.org/10.1007/s11749-014-0357-7
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HUrho')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Plot Huber rho function.
    close all
    x=-3:0.001:3;
    c=1.345;
    rhoHU=HUrho(x,c);
    plot(x,rhoHU,'LineStyle','-','LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel('$\rho (u,1.345)$','Interpreter','Latex')
    text(-c,0,'-c=-1.345')
    text(c,0,'c=1.345')
    hold('on')
    plot(x,x.^2/2,'LineStyle',':','LineWidth',1.5)
    stem(c,c^2/2)
    stem(-c,c^2/2)

%}

%{
    % Huber rhos function for two values of c.
    x=-6:0.01:6;
    c=1.345;
    rhoHU=HUrho(x,c);
    plot(x,rhoHU)
    xlabel('u','Interpreter','Latex')
    ylabel('$\rho (x,1.345)$','Interpreter','Latex')
    text(-c,0,'-c')
    text(c,0,'c')
    title('$\rho (u,c)$ with $c=1,345$ and $c=2$','Interpreter','Latex')
    hold('on')
    rhoHU=HUrho(x,2);
    plot(x,rhoHU,'-.')
%}

%% Beginning of code

w = (abs(u)<=c);
rhoHU = (u.^2/2).*w +(1-w).*(c*abs(u)-c^2/2);

end
%FScategory:UTISTAT