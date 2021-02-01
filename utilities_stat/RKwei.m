function weiRK = RKwei(u,c, M)
%RKwei computes weight function psi(u)/u for Rocke (translated Tukey's) biweight
%
%<a href="matlab: docsearchFS('RKwei')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        tuning parameter. Scalar. Non negative scalar which
%               (together with the other optional parameter M) controls the
%               robustness/efficiency of the estimator (beta in regression
%               or mu in the location case ...).
%    M :        tuning parameter. Scalar. Scalar greater than 0 which
%               (together with the other optional parameter c) controls the
%               robustness/efficiency of the estimator (beta in regression
%               or mu in the location case ...).
%
%  Optional input arguments:
%
%
%  Output:
%
%
%   weiRK :      n x 1 vector which contains the Rocke weights (translated
%                Tukey's biweight) associated to the residuals or
%                Mahalanobis distances for the n units of the sample.
%
% More About:
%
%
% function RKwei transforms vector u as follows
% \[
% RKwei(u)= \left\{
%    \begin{array}{cc}
% 1 &  0\leq u \leq M  \\
%  \left(1-\left( \frac{u-M}{c} \right)^2 \right)^2 &  M < u \leq M+c \\
% 0   &             u > M+c \\
% \end{array}
%    \right.
%  \]
%
% See equation (2.18) p. 1332 of Rocke (1996).
%
% See also HYPrho, HArho, OPTrho, TBrho
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
% Rocke D.M. (1996), Robustness properties of S-estimators of multivariate
% location and shape in high dimension, "The Annals of Statistics", Vol. 24,
% pp. 1327-1345.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('RKwei')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit

% Examples:

%{
   %% Plot of wei function.
    close all
    % Find the values of c and M given bdp=0.4 and v=5 for ARP=0.01
    x=0:0.01:5;
    bdp=0.4;
    v=5;
    ARP=0.01;
    [c,M]=RKbdp(bdp,v,ARP);
    weiRK=RKwei(x,c,M);
    % weiRK=weiRK/max(weiRK);
    plot(x,weiRK,'LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel('$w (u,c,M)$','Interpreter','Latex')
    title('$w (u,c,M)$','Interpreter','Latex')
    hold('on')
    stem(M,1,'LineStyle',':','LineWidth',1)
    text(M,0,'M')
    stem(M+c,0,'LineStyle',':','LineWidth',1)
    text(M+c,0,'M+c')
%}

%{
    %% Compare Rocke weight functions for 3 different values of bdp.
    close all
    x=0:0.01:6;
    % Number of variables v is fixed to 5
    v=5;
    % ARP is fixed to 0.01
    ARP=0.01;
    lwd=2;
    hold('on')
    % Use bdp=0.3
    bdp=0.3;
    [c,M]=RKbdp(bdp,v,ARP);
    weiRK030=RKwei(x,c,M);
    weiRK030=weiRK030/max(weiRK030);
    plot(x,weiRK030,'LineStyle','-','LineWidth',lwd)

    % Use bdp=0.4
    bdp=0.4;
    [c,M]=RKbdp(bdp,v,ARP);
    weiRK040=RKwei(x,c,M);
    weiRK040=weiRK040/max(weiRK040);
    plot(x,weiRK040,'LineStyle','-.','LineWidth',lwd)

    % Use bdp=0.5
    bdp=0.5;
    [c,M]=RKbdp(bdp,v,ARP);
    weiRK050=RKwei(x,c,M);
    weiRK050=weiRK050/max(weiRK050);
    plot(x,weiRK050,'LineStyle','--','LineWidth',lwd)
    
    xlabel('$x$','Interpreter','Latex','FontSize',16)
    ylabel('RK. Normalized $w(x,c,M)$','Interpreter','Latex','FontSize',20)
    legend({'$bdp=0.3$', '$bdp=0.4$' '$bdp=0.5$'},'Interpreter','Latex','Location','SouthEast','FontSize',16)
%}



%% Beginning of code

weiRK=zeros(size(u));

uu1 = u <= M;
uu2 = u > M & u <= M+c;

if ~isempty(uu1)
    weiRK(uu1) =  1;
end

if ~isempty(uu2)
    u2=u(uu2);
    weiRK(uu2)=(1- ((u2-M)/c).^2 ).^2;
 end

end
%FScategory:UTISTAT