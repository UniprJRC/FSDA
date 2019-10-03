function psiHYP = HYPpsi(u, cktuning)
%HYPpsi computes psi function for hyperbolic tangent estimator
%
%<a href="matlab: docsearchFS('HYPpsi')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    cktuning :  tuning parameters. Vector of length 2 or of length 5.
%                cktuning specifies the value of the tuning
%                constant c (scalar greater than 0 which controls the
%                robustness/efficiency of the estimator)
%                and the prefixed value k (sup of the
%                change-of-variance sensitivity) and the values of
%                parameters A, B and d.
%                cktuning(1) = c;
%                cktuning(2) = k = supCVC(psi,x) x \in R;
%                cktuning(3)=A;
%                cktuning(4)=B;
%                cktuning(5)=d;
%                Remark - if length(cktuning)==2 values of A, B and d will be
%                computed automatically
%
%  Optional input arguments:
%
%
%
%  Output:
%
%
%   psiHYP :     hyperbolic psi function. Vector. n x 1 vector which
%                contains the values of hyperbolic psi
%                function associated to the residuals or Mahalanobis
%                distances for the n units of the sample
%
% More About:
%
% Function HYPpsi transforms vector u as follows
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
% See also TBpsi, HApsi, OPTpsi
%
%
% References:
%
% Hampel, F.R., Rousseeuw, P.J. and  Ronchetti E. (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% "Journal of the American Statistical Association", Vol. 76,
% pp. 643-648 [HRR]
%
% Riani, M., Cerioli, A., Atkinson, A.C. and Perrotta, D. (2014), Monitoring
% Robust Regression, "Electronic Journal of Statistics", Vol. 8, pp. 646-677.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HYPpsi')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{

    % Plot of hyperbolic psi function.
    % Obtain Figure 2 of  p. 645 of HRR.
    x=-9:0.1:9;
    ctuning=6;
    ktuning=4.5;
    psiHYP=HYPpsi(x,[ctuning,ktuning]);
    plot(x,psiHYP)
    xlabel('$u$','Interpreter','Latex')
    ylabel(' Hyperbolic $\psi(u,c=6,k=4.5) $','Interpreter','Latex')
    text(ctuning,-0.1,'c','FontSize',14)
    text(-ctuning,0.1,'-c','FontSize',14)
%}

%{
    % Compare psi function for two values of parameter k.
    close all
    x=-9:0.1:9;
    ctuning=6;
    ktuning=5;
    psiHYP=HYPpsi(x,[ctuning,ktuning]);
    plot(x,psiHYP,'color','b')
    text(6,1.5,['k=' num2str(ktuning)],'Color','b','FontSize',14)
    xlabel('$u$','Interpreter','Latex','FontSize',14)
    ylabel(' Hyperbolic $\psi(u,c=6,k) $','Interpreter','Latex','FontSize',14)
    text(ctuning,-0.1,'c','FontSize',14)
    text(-ctuning,0.1,'-c','FontSize',14)
    hold('on')
    ktuning=4;
    psiHYP=HYPpsi(x,[ctuning,ktuning]);
    plot(x,psiHYP,'color','k')
    text(2,1.3,['k=' num2str(ktuning)],'Color','k','FontSize',14)

%}


%{
    % Parameters associated to a value of bdp=1/2.
    c=2.158325031399727
    k=4;
    A=0.000162707412432;
    B=0.006991738279441
    d=0.016982948780061
    x=-8:0.001:8;
    psiHYP=HYPpsi(x,[c,k,A,B,d]);
    plot(x,psiHYP)
    xlabel('x','Interpreter','Latex')
    ylabel(' Hyperbolic $\psi(x) $','Interpreter','Latex')

%}

%% Beginning of code

c = cktuning(1);
k = cktuning(2);
if length(cktuning)>2
    
    A=cktuning(3);
    B=cktuning(4);
    d=cktuning(5);
    
    if ((A < 0) || (B < A) || (B>1))
        disp('A must be >=0')
        disp('B must be >=A')
        disp('B must be <=1')
        disp(['B=' num2str(B) ' and A=' num2str(A)])
        error('FSDA:HYPpsi:WrongAorB','Illegal choice of parameters in hyperbolic tangent estimator:')
    else
    end
    
else
    % Find parameters A, B and d using routine HYPck
    [A,B,d]=HYPck(c,k);
    
    % For example if c=4 and k=5
    %     A = 0.857044;
    %     B = 0.911135;
    %     d =1.803134;
    % see Table 2 of HRR
end


psiHYP = zeros(size(u));
absu=abs(u);

%  u,		   |u| <=d
psiHYP(absu<=d) = u(absu<=d);


%                d <= |u| < c,
% \sqrt(A * (k - 1)) * tanh(sqrt((k - 1) * B^2/A)*(c -|u|)/2) .* sign(u)
psiHYP(absu > d & absu <=c) = sqrt(A * (k - 1)) * tanh(sqrt((k - 1) * B^2/A)...
    *(c - absu(absu > d & absu <=c ))/2) .* sign(u(absu > d & absu <=c));

% 0,			              |u| >= c.

end
%FScategory:UTISTAT