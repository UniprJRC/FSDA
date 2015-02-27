function w = HAwei(u, ctuning,varargin)
%HAwei computes weight function psi(u)/u using Hampel proposal
%
%<a href="matlab: docsearchFS('HAwei')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    ctuning :  scalar or vector of length 4 which specifies the value of the tuning
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
%                biweight with parameter 8.%
%
%  Output:
%
%    w :         n x 1 vector contains the Hampel weights associated to the residuals or
%                Mahalanobis distances for the n units of the sample
%
% Function HAwei transforms vector u as follows
%
% HAwei(u) = 	{ 1,			              |u| <= a,
%		        { a/|u|,		         a <= |u| < b,
%		        { a/|u|*(c-|u|)/(c-b),	 b <= |u| <  c,
%		        { 0,			              |u| >= c.
%
%             where a= ctuning(2) *ctuning(1)
%                   b= ctuning(3) *ctuning(1)
%                   c= ctuning(4) *ctuning(1)
%              
%             The default (if input ctuning is a scalar) is  
%                   a= 2*ctuning
%                   b= 4*ctuning
%                   c= 8*ctuning
%
%	It is necessary to have 0 <= a <= b <= c%
%
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
%<a href="matlab: docsearchFS('hawei')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

    % Obtain Figure 11.15 (panel b) p. 382 of
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
