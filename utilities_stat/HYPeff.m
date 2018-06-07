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
%                 Example - 'k',5
%                 Data Types - double
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
% pp. 643-648 [HRR]'
%
% Copyright 2008-2018.
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