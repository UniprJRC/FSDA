function psiHYPx = HYPpsix(u, cktuning)
%HYPpsix computes psi function for hyperbolic tangent estimator times x
%
%<a href="matlab: docsearchFS('hyppsix')">Link to the help function</a>
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
%  Output:
%
%
%   psiHYPx :    psi(u)*u function. Vector.
%                n x 1 vector which contains the values of hyperbolic
%                psi(u)*u function associated to the residuals or Mahalanobis
%                distances for the n units of the sample.
%
%
% More About:
%
% Function HYPpsix transforms vector $u$ as follows
% 
% \[
%  HYPpsix(u) = 
%  \left\{
%    \begin{array}{cc}
%  	 u^2 &        |u| \leq  d \\
% 
% \sqrt{A (k - 1)}  \tanh \left( \sqrt{(k - 1) B^2/A} (c -|u|)/2 \right) sign(u) u &
% 		         	                 d \leq |u| <  c, \\
%                0 &                      |u| \geq c.
% \end{array}
%    \right.
%  \]
%  	It is necessary to have $0 < A < B < 2 normcdf(c)-1- 2 c normpdf(c) <1$
%
%
% See also TBpsix, HApsix, OPTpsix
%
% References:
%
%
% Frank R. Hampel, Peter J. Rousseeuw and Elvezio Ronchetti (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% Journal of the American Statistical Association , Vol. 76, No. 375,
% pp. 643-648 (HRR)
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('hyppsi')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

    % plot of psi(x)*x for Hyperbolic estimator.
    x=-9:0.1:9;
    ctuning=6;
    ktuning=4.5;
    psiHYPx=HYPpsix(x,[ctuning,ktuning]);
    plot(x,psiHYPx)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi (x) \times x$','Interpreter','Latex')

%}

%% Beginning of code

c = cktuning(1);
k = cktuning(2);
if length(cktuning)>2

        A=cktuning(3);
        B=cktuning(4);
        d=cktuning(5);

    if ((A < 0) || (B < A) || (B>1)),
        error('FSDA:HYPpsix:WrongAorB',[' Illegal choice of parameters in hyperbolic tangent estimator: ' ...
            num2str(param) ]')
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

psiHYPx = zeros(size(u));
absu=abs(u);

%  u,		   |u| <=d
psiHYPx(absu<=d) = u(absu<=d).^2;


%                d <= |u| < c,
% \sqrt(A * (k - 1)) * tanh(sqrt((k - 1) * B^2/A)*(c -|u|)/2) .* sign(u)
psiHYPx(absu > d & absu <=c) = sqrt(A * (k - 1)) * tanh(sqrt((k - 1) * B^2/A)...
    *(c - absu(absu > d & absu <=c ))/2) .* sign(u(absu > d & absu <=c)).*(u(absu > d & absu <=c));

% 0,			              |u| >= c.

end
%FScategory:UTISTAT
