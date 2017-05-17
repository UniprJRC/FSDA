function rhoHYP = HYPrho(u, cktuning)
%HYPrho computes rho function  using hyperbolic tangent estimator
%
%<a href="matlab: docsearchFS('HYPrho')">Link to the help function</a>
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
%   rhoHYP :     rho function for hyperbolic tangent estimator. Vector. 
%                n x 1 vector which contains the hyperbolic rho
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
%
% More About:
%
% Hampel et al. (1981) have introduced a rho function which
% minimizes the asymptotic variance of the regression M-estimate, subject
% to a bound on the supremum of the Change of Variance Curve of the
% estimate. This leads to the Hyperbolic Tangent $\rho$
% function, which, for suitable constants $c$, $k$, $A$, $B$ and
% $d$, is defined as
% 
%
% \[
%  HYPrho(u) =
%  \left\{
%  \begin{array}{cc}
%   	 u^2/2 &	        |u| \leq d, \\
%   d^2/2 -2 \frac{A}{B} \log  \left\{ \cosh \left[ 0.5 \sqrt{ \frac{(k - 1)  B^2}{A} } (c - |u|) \right] \right\} & \\
%                          +2 \frac{A}{B}\log \left\{  \cosh \left[ 0.5\sqrt{\frac{(k - 1)  B^2}{A}}(c -d)\right] \right\} &  \\		        	
%                &                               d \leq |u| <  c, \\
%                 d^2/2 +2 \frac{A}{B} \log \left\{ \cosh \left[ 0.5 \sqrt{ \frac{(k - 1)  B^2}{A} }(c -d) \right] \right\}	 &
%                                                        |u| \geq c. \\
%   \end{array}
%   \right.
% \]
% where $0 < d < c$ is such that
% \[
% d = \sqrt{[A(k-1)]}\tanh [\frac{1}{2}\sqrt{\frac{(k-1)B^2}{A}}(c - d)],
% \]
% $A$ and $B$ satisfy suitable conditions, and $k$ is related to the bound
% in the Change of Variance Curve.
%
% More precisely, it is necessary to have $0 < A < B < 2 *normcdf(c)-1- 2*c*normpdf(c) <1$
%
% See also TBrho, HArho, OPTrho
%
% References:
%
% Hampel F.R., Rousseeuw P.J. and Ronchetti E. (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% Journal of the American Statistical Association , Vol. 76, No. 375,
% pp. 643-648 (HRR)
% Riani M., Cerioli A., Atkinson A.C., Perrotta D.  (2014). Monitoring
% Robust Regression. Electronic Journal of Statistics, Vol. 8 pp.  646-677
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HYPrho')">Link to the help page for this function</a>
% Last modified 31-05-2016
%
% Examples:

%{

    % Plot of rho function for hyperbolic tangent estimator.
    x=-8:0.001:8;
    ctuning=6;
    ktuning=4.5;
    rhoHYP=HYPrho(x,[ctuning,ktuning]);
    plot(x,rhoHYP)
    xlabel('x','Interpreter','Latex')
    ylabel(' Hyperbolic $\rho(x) $','Interpreter','Latex')

%}

%{
    % Parameters associated to a value of bdp=1/2.
    c=2.158325031399727
    k=4;
    A=0.000162707412432;
    B=0.006991738279441   
    d=0.016982948780061
    x=-8:0.001:8;
    rhoHYP=HYPrho(x,[c,k,A,B,d]);
    plot(x,rhoHYP)
    xlabel('x','Interpreter','Latex')
    ylabel(' Hyperbolic $\rho(x) $','Interpreter','Latex')

%}

%% Beginning of code

c = cktuning(1);
k=cktuning(2);

if length(cktuning)>2

        A=cktuning(3);
        B=cktuning(4);
        d=cktuning(5);

    if ((A < 0) || (B < A) || (B>1))
        error('FSDA:HYPrho:WrongAorB',[' Illegal choice of parameters in hyperbolic tangent estimator: ' ...
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

rhoHYP = zeros(size(u));
absu=abs(u);

%  u,		   |u| <=d
rhoHYP(absu<=d) = u(absu<=d).^2/2;


%                   d <= |u| < c,
rhoHYP(absu > d & absu <=c) =  d^2/2    ...
    -2*(A/B) * log(cosh(0.5*sqrt((k - 1) * B^2/A)...
    *(c - absu(absu > d & absu <=c ))))...
    +(2*A/B)*log(cosh(0.5*sqrt((k - 1) * B^2/A)*(c -d)));

%|u| >= c.
rhoHYP(absu > c) = d^2/2 +2*(A/B)*log(cosh(0.5*sqrt((k - 1) * B^2/A)*(c -d)));

end

%FScategory:UTISTAT
