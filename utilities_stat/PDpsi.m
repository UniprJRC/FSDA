function psi=PDpsi(u,alpha)
%PDpsi computes psi function (derivative of rho function) for minimum density power divergence estimator
%
%
%<a href="matlab: docsearchFS('PDpsi')">Link to the help function</a>
%
%
%  Required input arguments:
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
%
%   PDpsi :      n x 1 vector which contains the minimum density power
%                divergence (MDPD) psi function
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
% More About:
%
%
% function PDpsi transforms vector u as follows
%  \[
% PDpsi(u,alpha)=  \alpha u \exp(-\alpha (u^2/2));
%  \]
%
%
% See also TBpsi, HYPpsi, HApsi, OPTpsi, HUpsi
%
% References:
%
%  Riani, M. Atkinson, A.C., Corbellini A. and Perrotta A. (2020), Robust
%  Regression with Density Power Divergence: Theory, Comparisons and Data
%  Analysis, submitted.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('PDpsi')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Plot of psi function.
    close all
    x=-6:0.01:6;
    alpha=1;
    psiPD=PDpsi(x,alpha);
    plot(x,psiPD,'LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel(['$\psi(u,' num2str(alpha) ')$'],'Interpreter','Latex','FontSize',14)
    hold('on')
%}

%{
    %% Comparison of psi function for two values of alpha.
    close all
    hold('on')
    x=-6:0.01:6;
    alpha1=1;
    psiPD1=PDpsi(x,alpha1);
    plot(x,psiPD1,'LineWidth',2)
    bdp1=1-1/sqrt(1+alpha1);
    eff1=(sqrt(1+2*alpha1)/(1+alpha1))^3;
    alpha2=0.1;
    psiPD2=PDpsi(x,alpha2);
    plot(x,psiPD2,'LineWidth',2)
    bdp2=1-1/sqrt(1+alpha2);
    eff2=(sqrt(1+2*alpha2)/(1+alpha2))^3;

    xlabel('$u$','Interpreter','Latex')
    ylabel('$\psi(u,\alpha)$','Interpreter','Latex','FontSize',14)
    legend({['$\alpha=' num2str(alpha1) '\mapsto bdp=' num2str(bdp1,2) '\;  eff=' num2str(eff1,2) '$'], ...
        ['$\alpha=' num2str(alpha2) '\mapsto bdp=' num2str(bdp2,2) '\;  eff=' num2str(eff2,2) '$']},...
        'Interpreter','Latex','Location','SouthEast','FontSize',12)
%}

%% Beginning of code

psi = alpha * u.*exp(- alpha *(u.^2/2));

end
%FScategory:UTISTAT