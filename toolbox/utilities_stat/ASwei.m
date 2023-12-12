function weiAS = ASwei(u,c)
%ASwei computes weight function psi(u)/u for Andrew's sine function  
%
%<a href="matlab: docsearchFS('ASwei')">Link to the help function</a>
%
%
% Required input arguments:
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
%    weiAS :     n x 1 vector which contains the Andrew's sine weights
%                associated to the scaled residuals or Mahalanobis
%                distances for the n units of the sample.
%
% More About:
%
% Function ASwei transforms vector u as follows 
%
% \[
% ASwei(u)= \left\{
%    \begin{array}{cc}
%   sin(u/c)/(u/c)      & |u/c| \leq \pi \\
%  0                     &  |u/c|>\pi   \\
% \end{array}
%    \right.
%  \]
%
% Remark: Andrews's psi-function is almost linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is approximately constant over the linear region of \psi,
% so the points in that region tend to get equal weight.
%
%
% See also: TBwei, HYPwei, HAwei, OPTwei
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
%<a href="matlab: docsearchFS('ASwei')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Weight function for Andrew's sine link.
    x=-6:0.01:6;
    c=1.5;
    weiAS=ASwei(x,c);
    plot(x,weiAS)
    xlabel('x','Interpreter','Latex')
    ylabel('$W (x) =\psi(x)/x$','Interpreter','Latex')
%}

%{
  close all
    x=-6:0.01:6;
    lwd=2;
    hold('on')
    bdp1=0.25;
    c1=ASbdp(bdp1);    
    weiAS=ASwei(x,c1);
    weiAS=weiAS/max(weiAS);
    plot(x,weiAS,'LineStyle','-','LineWidth',lwd)
    bdp2=0.01;
    c=ASbdp(bdp2);
    weiAS=ASwei(x,c);
    weiAS=weiAS/max(weiAS);
    plot(x,weiAS,'LineStyle','-.','LineWidth',lwd)
    xlabel('$x$','Interpreter','Latex','FontSize',16)
    ylabel('AS weight function $\psi_c/x$','Interpreter','Latex','FontSize',20)
    legend({['bdp=' num2str(bdp1,2)],  ['bdp=' num2str(bdp2,2)]},...
    'Interpreter','Latex','Location','SouthEast','FontSize',12)
%}


%{
    %% Compare six different weight functions in each of them eff is 95 per cent.
    % Initialize graphical parameters.
    FontSize=14;
    x=-6:0.01:6;
    ylim1=-0.05;
    ylim2=1.05;
    xlim1=min(x);
    xlim2=max(x);
    LineWidth=2;

    subplot(2,3,1)
    ceff095HU=HUeff(0.95,1);
    weiHU=HUwei(x,ceff095HU);
    plot(x,weiHU,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Huber','FontSize',FontSize)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])

    subplot(2,3,2)
    ceff095HA=HAeff(0.95,1);
    weiHA=HAwei(x,ceff095HA);
    plot(x,weiHA,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Hampel','FontSize',FontSize)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])


    subplot(2,3,3)
    ceff095TB=TBeff(0.95,1);
    weiTB=TBwei(x,ceff095TB);
    plot(x,weiTB,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Tukey biweight','FontSize',FontSize)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])

    subplot(2,3,4)
    ceff095HYP=HYPeff(0.95,1);
    ktuning=4.5;
    weiHYP=HYPwei(x,[ceff095HYP,ktuning]);
    plot(x,weiHYP,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Hyperbolic','FontSize',FontSize)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])

    subplot(2,3,5)
    ceff095PD=PDeff(0.95);
    weiPD=PDwei(x,ceff095PD);
    weiPD=weiPD/max(weiPD);
    plot(x,weiPD,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Power divergence','FontSize',FontSize)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])

    subplot(2,3,6)
    ceff095AS=ASeff(0.95);
    weiAS=ASwei(x,ceff095AS);
    weiAS=weiAS/max(weiAS);
    plot(x,weiAS,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Andrew''s sine link','FontSize',FontSize)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])


%}


%% Beginning of code
c=c(1); % MATLAB Ccoder instruction to enforce that alpha is a scalar
weiAS = zeros(size(u));
inds1 = abs(u) < c*pi;
weiAS(inds1) = (sin(u(inds1)/c)./ (u(inds1)/c));

end
%FScategory:UTISTAT