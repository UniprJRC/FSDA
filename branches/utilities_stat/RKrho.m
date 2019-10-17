function rhoRK = RKrho(u,c, M)
%RKrho computes rho function for Rocke (translated Tukey's) biweight
%
%<a href="matlab: docsearchFS('RKrho')">Link to the help function</a>
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
%   rhoRK :      n x 1 vector which contains the Rocke rho (translated
%                Tukey's biweight) rho associated to the residuals or
%                Mahalanobis distances for the n units of the sample.
%
% More About:
%
%
% function RKrho transforms vector u as follows
% \[
% RKrho(u)= \left\{
%    \begin{array}{cc}
% \frac{u^2}{2} &  0\leq u \leq M  \\
%  \frac{M^2}{2} -M^2\frac{M^4-5 M^2 c^2 + 15c^4}{30c^4} + u^2 \left( 0.5+ \frac{M^4}{2c^4}  -\frac{M^2}{c^2} \right) \\
%  +u^3 \left( \frac{4M}{3c^2} -\frac{4 M^3}{3c^4} \right) +u^4 \left( \frac{3M^2}{2c^4}- \frac{1}{2c^2} \right) \\
%  -4M \frac{u^5}{5c^4} + \frac{u^6}{6c^4} &         M < u \leq M+c \\
% \frac{M^2}{2} + \frac{c(5c+ 16M)}{30}   &             u > M+c \\
% \end{array}
%    \right.
%  \]
%
% See equation (2.20) p. 1333 of Rocke (1996).
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
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('RKrho')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit

% Examples:

%{
   %% Plot of rho function.
    close all
    % Find the values of c and M given bdp=0.4 and v=5 for ARP=0.01
    x=0:0.01:5;
    bdp=0.4;
    v=5;
    ARP=0.01;
    [c,M]=RKbdp(bdp,v,ARP);
    rhoRK=RKrho(x,c,M);
    % rhoRK=rhoRK/max(rhoRK);
    plot(x,rhoRK,'LineWidth',2)
    xlabel('$u$','Interpreter','Latex')
    ylabel('$\rho (u,c,M)$','Interpreter','Latex')
    title('$\rho (u,c,M)$','Interpreter','Latex')
    hold('on')
    stem(M,M^2/2,'LineStyle',':','LineWidth',1)
    text(M,0,'M')
    stem(M+c,M^2/2+c*(5*c+16*M)/30,'LineStyle',':','LineWidth',1)
        text(M+c,0,'M+c')

%}

%{
    %% Compare Rocke rho functions for 3 different values of bdp.
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
    rhoRK030=RKrho(x,c,M);
    rhoRK030=rhoRK030/max(rhoRK030);
    plot(x,rhoRK030,'LineStyle','-','LineWidth',lwd)

    % Use bdp=0.4
    bdp=0.4;
    [c,M]=RKbdp(bdp,v,ARP);
    rhoRK040=RKrho(x,c,M);
    rhoRK040=rhoRK040/max(rhoRK040);
    plot(x,rhoRK040,'LineStyle','-.','LineWidth',lwd)

    % Use bdp=0.5
    bdp=0.5;
    [c,M]=RKbdp(bdp,v,ARP);
    rhoRK050=RKrho(x,c,M);
    rhoRK050=rhoRK050/max(rhoRK050);
    plot(x,rhoRK050,'LineStyle','--','LineWidth',lwd)
    
    xlabel('$x$','Interpreter','Latex','FontSize',16)
    ylabel('RK. Normalized $\rho(x,c,M)$','Interpreter','Latex','FontSize',20)
    legend({'$bdp=0.3$', '$bdp=0.4$' '$bdp=0.5$'},'Interpreter','Latex','Location','SouthEast','FontSize',16)
%}

%% Beginning of code

rhoRK=zeros(size(u));

uu1 = u <= M;
uu2 = u > M & u <= M+c;
uu3 = u > M+c;

if ~isempty(uu1)
    rhoRK(uu1) =  u(uu1).^2/2;
end

if ~isempty(uu2)
    u2=u(uu2);
    rhoRK(uu2)=M^2/2-M^2*(M^4-5*M^2*c^2 + 15*c^4)/(30*c^4) + u2.^2*(0.5+M^4/(2*c^4)-M^2/(c^2)) ...
        +u2.^3*(4*M/(3*c^2) -4*M^3/(3*c^4)) +u2.^4*(3*M^2/(2*c^4)-1/(2*c^2)) ...
        -4*M*u2.^5/(5*c^4) +u2.^6 /(6*c^4);
end

if ~isempty(uu3)
    rhoRK(uu3)= M^2/2+c*(5*c+16*M)/30;
end

end
%FScategory:UTISTAT