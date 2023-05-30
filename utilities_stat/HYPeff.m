function [c,A,B,d] = HYPeff(eff, v, k, traceiter)
%HYPeff finds constant c which is associated to the requested efficiency for hyperbolic estimator
%
%
%<a href="matlab: docsearchFS('HYPeff')">Link to the help page for this function</a>
%
%  Required input arguments:
%
%    eff:       efficiency. Scalar.  Scalar which contains the required
%               efficiency (of location
%               or scale estimator).
%               Generally eff=0.85, 0.9 or 0.95
%    v :        number of response variables. Scalar. Number of variables of
%               the  dataset (for regression v=1)
%               UP TO NOW v=1 (JUST REGRESSION) TO DO FOR MULTIVARIATE
%               ANALYSIS
%
%  Optional input arguments:
%
%   k        : supremum of the change of variance curve. Scalar.
%              $\sup CVC(psi,x) x \in R$
%              Default value is k=4.5.
%                 Example - 5
%                 Data Types - double
%
%  traceiter : Level of display. Scalar.
%              If traceiter = 1 it is possible to monitor
%              how the value of the objective function B^2/A
%              gets closer to the target (eff) during the iterations
%                 Example - 'traceiter',0
%                 Data Types - double
%
% Output:
%
%  c    : parameter c of hyperbolic tangent estimator. Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%   A   : parameter A of hyperbolic tangent estimator. Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%   B   : parameter B of hyperbolic tangent estimator. Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%   d   : parameter d of hyperbolic tangent estimator. Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%
% More About:
%
%  \[
%   HYPpsi(u) =
% \left\{
%   \begin{array}{cc}
%  	 u &        |u| \leq  d \\
%                  \sqrt{A (k - 1)}  \tanh \left( \sqrt{(k - 1) B^2/A} (c -|u|)/2 \right) sign(u) &
% 		         	                 d \leq |u| <  c, \\
%                 0 &                      |u| \geq c.
% \end{array}
%    \right.
%  \]
%  	It is necessary to have $0 < A < B < 2 normcdf(c)-1- 2 c \times normpdf(c) <1$
%
%
% See also TBeff, HAeff, OPTeff
%
%
%
% References:
%
%
% Hampel, F.R., Rousseeuw, P.J. and  Ronchetti E. (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% "Journal of the American Statistical Association", Vol. 76,
% pp. 643-648. [HRR]
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('HYPeff')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:

%{
    % Find parameter c for fixed efficiency. 
    % Find value of c for a nominal efficiency of 0.9
    % when k=4.5 (default value).
    c=HYPeff(0.9,1);
    % In this case
    % c = 3.3139
%}

%{
    % Example of specifying k. 
    % Find value of c for a nominal efficiency of 0.9
    % when k=4
    ktuning=4;
    c=HYPeff(0.9,1,ktuning);
    % In this case
    % c =  3.5443
%}


%{
    % Find parameters for fixed efficiency and k. 
    % Find value of c, A, B, for a nominal efficiency of 0.8427
    % when k=4.5
 
    ktuning=4.5;
    [c,A,B,d]=HYPeff(0.8427,1,ktuning);
    % In this case
    % c = 3.000130564905703
    % A = 0.604298601602487
    % B = 0.713612241773758
    % d= 1.304379168746527
    % See also Table 2 of HRR p. 645
%}

%{
    % Example of use of option Find parameters for fixed efficiency and k. 
    % Find value of c, A, B, for a nominal efficiency of 0.8427
    % when k=4.5
 
    ktuning=4.5;
    traceiter =true;
    [c,A,B,d]=HYPeff(0.8427,1,ktuning,traceiter);
    % In this case
    % c = 3.000130564905703
    % A = 0.604298601602487
    % B = 0.713612241773758
    % d= 1.304379168746527
    % See also Table 2 of HRR p. 645
%}

%{
    % Example to show the issue of multiple solutions problem. 
    % 5 redescending psi functions are used and the Huber psi.
    load Income1;
    y=Income1{:,"HTOTVAL"};
    % Use contaminated income data 
    y=[y(1:20); 600000; 575000; 590000];
    
    y=y';
    mady=mad(y,1)/0.675;
    
    eff=0.95;
    TBc=TBeff(eff,1);
    HUc=HUeff(eff,1);
    HAc=HAeff(eff,1);
    HYPc=HYPeff(eff,1);
    OPTc=OPTeff(eff,1);
    PDc=PDeff(eff);
    
    
    mu=0:1000:700000;
    avePSI=zeros(length(mu),6);
    for i=1:length(mu)
        % aveTB(i,2)=mean(TBrho((y-mu(i))./mady,c));
    
        avePSI(i,1)=mean(HUpsi((y-mu(i))./mady,HUc));
        avePSI(i,2)=mean(HApsi((y-mu(i))./mady,HAc));
        avePSI(i,3)=mean(TBpsi((y-mu(i))./mady,TBc));
        avePSI(i,4)=mean(HYPpsi((y-mu(i))./mady,[HYPc,5]));
        avePSI(i,5)=mean(OPTpsi((y-mu(i))./mady,OPTc));
        avePSI(i,6)=mean(PDpsi((y-mu(i))./mady,PDc));
    
    end
    % Plotting part
    close
    Link={'Huber', 'Hampel', 'Tukey', 'Hyperbolic' 'Optimal' 'Power divergence'} ;
    for i=1:6
        subplot(2,3,i)
        plot(mu',avePSI(:,i),'LineWidth',2,'Color','k')
        hold('on')
        yline(0) %  line([min(mu);max(mu)],[0;0],'LineStyle',':')
        title(Link(i),'FontSize',14)
        xlabel('$\mu$','FontSize',14,'Interpreter','Latex')
        ylabel('$\overline \psi \left( \frac{ y -\mu}{\hat \sigma} \right)$','FontSize',14,'Interpreter','Latex')
    end
%}

%{
    % Compare the weight function for 6 different links.
    FontSize=14;
    FontSizetitl=12;
    x=-6:0.01:6;
    ylim1=-0.05;
    ylim2=1.05;
    xlim1=min(x);
    xlim2=max(x);
    LineWidth=2;
    
    subplot(2,3,1)
    ceff05HU=HUeff(0.95,1);
    weiHU=HUwei(x,ceff05HU);
    plot(x,weiHU,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Huber','FontSize',FontSizetitl)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])
    
    subplot(2,3,2)
    ceff095HA=HAeff(0.95,1);
    weiHA=HAwei(x,ceff095HA);
    plot(x,weiHA,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Hampel','FontSize',FontSizetitl)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])
    
    
    subplot(2,3,3)
    ceff095TB=TBeff(0.95,1);
    weiTB=TBwei(x,ceff095TB);
    plot(x,weiTB,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Tukey biweight','FontSize',FontSizetitl)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])
    
    subplot(2,3,4)
    ceff095HYP=HYPeff(0.95,1);
    ktuning=4.5;
    weiHYP=HYPwei(x,[ceff095HYP,ktuning]);
    plot(x,weiHYP,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Hyperbolic','FontSize',FontSizetitl)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])
    
    
    subplot(2,3,5)
    ceff095OPT=OPTeff(0.95,1);
    % ceff095OPT=ceff095OPT/3;
    weiOPT=OPTwei(x,ceff095OPT);
    weiOPT=weiOPT/max(weiOPT);
    plot(x,weiOPT,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Optimal','FontSize',FontSizetitl)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])
    
    subplot(2,3,6)
    ceff095PD=PDeff(0.95);
    weiPD=PDwei(x,ceff095PD);
    weiPD=weiPD/max(weiPD);
    plot(x,weiPD,'LineWidth',LineWidth)
    xlabel('$u$','Interpreter','Latex','FontSize',FontSize)
    title('Power divergence','FontSize',FontSizetitl)
    ylim([ylim1 ylim2])
    xlim([xlim1 xlim2])
%}

%% Beginning of code


if (nargin >2)
    if (k < 0) 
        error('FSDA:HYPeff:WrongK',[' Illegal choice of parameters in hyperbolic tangent estimator: ' ...
            num2str(k) ]')
    end
else
    k=4.5;
end


if nargin <4
    traceiter =0;
end

if v==1
    
    % c=4 starting value of the iteration
    c=4;
    step=2;
    empeff=0;
    
    eps=1e-8;
    iter=0;
    while abs(empeff-eff)>eps
        
        % Find parameters A, B and d using routine HYPck
        % disp(c)
        [A,B,d]=HYPck(c,k);
        
        iter=iter+1;
        
        if iter==100
            disp(['Effective tolerance in routine HYPeff=' num2str(abs(empeff-eff))])
            break
        end
        
        if traceiter==1
            disp([iter c empeff eff])
        end
        
        step=step/2;
        empeff=B^2/A;
        if empeff<eff
            c=c+step;
        elseif empeff>eff
            c=max(c-step,0.1);
        else
        end
    end
else
    error('FSDA:HYPeff:Wrongv','Not yet implemented for v>1')
end

end
%FScategory:UTISTAT