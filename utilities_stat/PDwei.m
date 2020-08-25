function w = PDwei(u,alpha)
%PDwei computes weight function psi(u)/u for  for minimum density power divergence estimator  
%
%<a href="matlab: docsearchFS('PDwei')">Link to the help function</a>
%
%
% Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    alpha :    tuning parameter. Scalar. Scalar in the interval (0,1] which
%               controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...). The
%               greater is alpha the greater is the bdp and smaller is the
%               efficiency.
%
%  Optional input arguments:
%
%  Output:
%
%    w :         n x 1 vector which contains the Tukey's biweight weights
%                associated to the scaled residuals or Mahalanobis
%                distances for the n units of the sample.
%
% More About:
%
% function PDwei transforms vector u as follows
% \[
% PDwei(u,alpha)=  \alpha \exp(-\alpha (u^2/2));
%      \]
%
%
%
%
% See also: TBwei, HYPwei, HAwei, OPTwei
%
% References:
%
%  Riani, M. Atkinson, A.C., Corbellini A. and Perrotta A. (2020), Robust
%  Regression with Density Power Divergence: Theory, Comparisons and Data
%  Analysis, submitted.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('PDwei')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Weight function for minimum density power divergence estimator.
    x=-6:0.01:6;
    alpha=0.01;
    weiPD=PDwei(x,alpha);
    plot(x,weiPD)
    xlabel('x','Interpreter','Latex')
    ylabel('$W (x) =\psi(x)/x$','Interpreter','Latex')

%}

%{
    % Comparing two weight functions for two values of bdp.
    close all
    x=-6:0.01:6;
    lwd=2;
    hold('on')
    bdp1=0.25;
    alpha1=PDbdp(bdp1);    
    weiPD=PDwei(x,alpha1);
    weiPD=weiPD/max(weiPD);
    plot(x,weiPD,'LineStyle','-','LineWidth',lwd)
    bdp2=0.01;
    c=PDbdp(bdp2);
    weiPD=PDwei(x,c);
    weiPD=weiPD/max(weiPD);
    plot(x,weiPD,'LineStyle','-.','LineWidth',lwd)
    xlabel('$x$','Interpreter','Latex','FontSize',16)
    ylabel('PD weight function $\psi_\alpha(x)/x$','Interpreter','Latex','FontSize',20)
    legend({['bdp=' num2str(bdp1,2)],  ['bdp=' num2str(bdp2,2)]},...
    'Interpreter','Latex','Location','SouthEast','FontSize',12)
%}


%{
    %% Compare five different weight functions in each of them eff is 95 per cent.
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
%}


%% Beginning of code

% normalized wights in such a way that when u=0 w=1
w = exp(- alpha *(u.^2/2));
% Unnormalized weights are
% w = alpha * exp(- alpha *(u.^2/2));

end
%FScategory:UTISTAT