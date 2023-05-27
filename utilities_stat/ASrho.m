function rhoAS = ASrho(u,c)
%ASrho computes rho function for Andrew's sine function
%
%<a href="matlab: docsearchFS('ASrho')">Link to the help function</a>
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
%   rhoAS :      vector of length n which contains the Andrew's sine function rho
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
% More About:
%
%
% function ASrho transforms vector u as follows 
% \[
% ASrho(u)= \left\{
%    \begin{array}{cc}
%   c (1-\cos (u / c))                   &  |u/c| \leq \pi  \\
%   2c                                   &  |u/c| > \pi   \\
% \end{array}
%    \right.
%  \]
%  
%
% See also HYPrho, HArho, OPTrho, TBrho, PDrho
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
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ASrho')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Plot of rho function.
    close all
    x=-7:0.01:7;
    rhoAS=ASrho(x,2);
    plot(x,rhoAS,'LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel('$\rho (u,2)$','Interpreter','Latex')
    text(-c*pi-1,2*c-0.1,'2*c')
    text(+c*pi+0.5,2*c-0.1,'2*c')
    title('$\rho (u,c)$','Interpreter','Latex')
    hold('on')
    c=2;
    stem(c*pi,2*c,'LineStyle',':','LineWidth',1)
    stem(-c*pi,2*c,'LineStyle',':','LineWidth',1)
    text(c*pi-0.8,0.1,'c \pi')
    text(-c*pi+0.2,0.1,'-c \pi')
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
    c=ASbdp(0.5,1);
    rhoAS=ASrho(x,c);
    rhoAS=rhoAS/max(rhoAS);
    plot(x,rhoAS,'LineStyle','-','LineWidth',lwd)
    c=ASeff(0.95,1);
    rhoAS=ASrho(x,c);
    rhoAS=rhoAS/max(rhoAS);
    plot(x,rhoAS,'LineStyle','-.','LineWidth',lwd)
    xlabel('$x$','Interpreter','Latex','FontSize',16)
    ylabel('AS. Normalized $\rho_c(x)$','Interpreter','Latex','FontSize',20)
    legend({'$c_{(bdp=0.5 \mapsto eff=0.2856)}$', '$c_{(eff=0.95 \mapsto bdp= 0.1217)}$'},'Interpreter','Latex','Location','SouthEast','FontSize',16)
%}

%% Beginning of code
c=c(1); % MATLAB Ccoder instruction to enforce that c is a scalar
w = (abs(u)<=c*pi);
rhoAS = c*(1-cos(u/c)).*w +(1-w)*2*c;

end
%FScategory:UTISTAT