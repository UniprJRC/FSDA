function w = HYPwei(u, cktuning)
%HYPwei computes weight function psi(u)/u for hyperbolic tangent estimator
%
%<a href="matlab: docsearchFS('HYPwei')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    cktuning :  tuning parameters. Vector of length 2 or of length 5.
%                cktuning specifies  specifies the value of the tuning
%                constant c (scalar greater than 0 which controls the
%                robustness/efficiency of the estimator)
%                and the prefixed value k (sup of the
%                change-of-variance sensitivity) and the values of
%                parameters A, B and d:
%                cktuning(1) = c;
%                cktuning(2) = k = supCVC(psi,x) x \in R;
%                cktuning(3)=A;
%                cktuning(4)=B;
%                cktuning(5)=d;
%                Remark: if length(cktuning)==2 values of A, B and d will be
%                computed automatically
%
%  Optional input arguments:
%
%  Output:
%
%    w :         hyperbolic weights. Vector.
%                n x 1 vector contains the hyperbolic weights associated to
%                the scaled residuals or Mahalanobis distances for the n
%                units of the sample
%
%
% More About:
%
% Function HYPwei transforms vector u as follows
%
% \[
% HYPwei(u) =
% \left\{
%    \begin{array}{cc}
% 	1 &	          |u| \leq d, \\
%		         \sqrt(A * (k - 1)) * tanh(sqrt((k - 1) * B^2/A)*(c-|u|)/2) .* sign(u)/u
%		         	 &                d \leq |u| <  c, \\
%		    0	&                 |u| \geq c. \\
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
% More precisely, it is necessary to have $0 < A < B < 2 *normcdf(c)-1- 2*c*normpdf(c) <1$%
% Remark: hyperbolic  psi-function is linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is approximately constant over the linear region of \psi,
% so the points in that region tend to get equal weight.
%
% See also TBwei, HAwei, OPTwei
%
% References:
%
%
% Hampel, F.R., Rousseeuw, P.J. and Ronchetti, E. (1981), 
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% "Journal of the American Statistical Association", Vol. 76, 
% pp. 643-648. [HRR]
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HYPwei')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:

%{
    % Weight function for hyperbolic tangent estimator.
    x=-6:0.01:6;
    ctuning=4;
    ktuning=4.5;
    weiHYP=HYPwei(x,[ctuning,ktuning]);
    plot(x,weiHYP)
    xlabel('x','Interpreter','Latex')
    ylabel('$W (x) =\psi(x)/x$','Interpreter','Latex')

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



w = zeros(size(u));
absu=abs(u);

%  u,		   |u| <=d
w(absu<=d) = 1;


%                d <= |u| < c,
% \sqrt(A * (k - 1)) * tanh(sqrt((k - 1) * B^2/A)*(c -|u|)/2) .* sign(u)
w(absu > d & absu <=c) = sqrt(A * (k - 1)) * tanh(sqrt((k - 1) * B^2/A)...
    *(c - absu(absu > d & absu <=c ))/2) .* sign(u(absu > d & absu <=c))./u(absu > d & absu <=c);

% 0,			              |u| >= c.

end
%FScategory:UTISTAT