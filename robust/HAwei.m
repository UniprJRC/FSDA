function w = HAwei(u, ctuning)
%HAwei computes weight function psi(u)/u using Hampel proposal
%
%<a href="matlab: docsearchFS('HAwei')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    ctuning :  tuning parameters. Scalar or Vector. Scalar or vector of
%               length 4 which specifies the value of the tuning
%                constant c (scalar greater than 0 which controls the
%                robustness/efficiency of the estimator)
%                and the prefixed values of paramters a, b, c
%                ctuning(1) = tuning constant which will multiply
%                parameters a, b and c of Hampel rho (psi) function
%                ctuning(2) = paramter a of Hampel rho (psi) function
%                ctuning(3) = paramter b of Hampel rho (psi) function
%                ctuning(4) = paramter c of Hampel rho (psi) function
%                Remark: if length(ctuning)==1 values of a, b and c will be
%                set to a=2*ctuning b=4*ctuning c=4*ctuning 
%                With these choices, if ctuning=1  the
%                resulting influence function is nearly identical to the
%                biweight with parameter 8.
%
%
% Optional input arguments:
%
%  Output:
%
%
%  Output:
%
%    w :         n x 1 vector which contains the Hampel weights associated
%                to the residuals or Mahalanobis distances for the n units
%                of the sample.
%
%
% More About:
%
% Function HAwei transforms vector u as follows
%
%
%  \[
%  HAwei(u)  = \left\{   
%  \begin{array}{cc}
%    1 & |u| <= a                                       \\
%    \frac{a}{|u|}   & a <= |u| < b                    \\
%    \frac{a}{|u|} *  \frac{c-|u|}{c-b},  & b <= |u| <  c \\
%    0 & |u| >= c 
%  \end{array} \right.
% \]
%
%           
%             where ctun=ctuning(1).
%                   $a$= ctun *ctuning(2).
%                   $b$= ctun *ctuning(3).
%                   $c$= ctun *ctuning(4).
%
%             The default is
%                   $a$= 2*ctun. 
%                   $b$= 4*ctun. 
%                   $c$= 8*ctun. 
%
%
% See also: TBwei, HYPwei, OPTwei
%
% References:
%
% D. C. Hoaglin, F. Mosteller, J. W. Tukey (1982), Understanding Robust and
% Exploratory Data Analysis Wiley, New York.
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HAwei')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

    %% Obtain Figure 11.15 (panel b) p. 382 of
    % Hoaglin et al. (1987)
    x=-8:0.01:8;
    weiHA=HAwei(x,[1 2 4 8]);
    plot(x,weiHA)
    xlabel('x','Interpreter','Latex')
    ylabel('$W (x) =\psi(x)/x$','Interpreter','Latex')

%}

%% Beginning of code


if length(ctuning)>1,
    
    if ((ctuning(2) < 0) || (ctuning(3) < ctuning(2)) || (ctuning(4) < ctuning(3)))
        error('FSDA:HAwei:WrongAbc',[' Illegal choice of parameters in Hampel: ' ...
            num2str(ctuning(2:4)) ]')
    end
    a =  ctuning(2)*ctuning(1);
    b =  ctuning(3)*ctuning(1);
    c =  ctuning(4)*ctuning(1);
else
    a = 2*ctuning;
    b = 4*ctuning;
    c = 8*ctuning;
end


w = ones(size(u));
absu=abs(u);

% a/|u|,		 a <= |u| < b,
w(absu >= a) = a*(abs(u(absu >= a))).^(-1);

%a/|u| * (c-|u|)/(c-b),	 b <= |u| <  c,
w(absu >= b) = w(absu >= b).*( c  -abs(u(absu >= b)) )/(c-b);

% 0,			 |u| >= c
w(absu > c) = 0;

end
