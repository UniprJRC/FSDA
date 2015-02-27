function rhoHA = HArho(u, ctuning,varargin)
%HArho computes rho function  using Hampel proposal
%
%<a href="matlab: docsearchFS('harho')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    ctuning :  scalar or vector of length 4 which specifies the value of the tuning
%                constant c (scalar greater than 0 which controls the
%                robustness/efficiency of the estimator)
%                and the prefixed values of parameters a, b, c
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
%  Output:
%
%
%   rhoHA :     n x 1 vector which contains the Hampel rho
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample
%
% Function HArho transforms vector u as follows
%
% HArho(u) = 	{ u^2/2,			                               |u| <= a,
%		        { a*|u|-0.5 a^2,		                        a <= |u| < b,
%		        { ab-0.5a^2+0.5*(c-b)*a(1- ((c-|u|)/(c-B))^2 ),	  b <= |u| <  c,
%		        { ab-0.5a^2+0.5*(c-b)*a,			              |u| >= c.
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
%	It is necessary to have 0 <= a <= b <= c
%
%
% See also HYPrho, Tbrho, OPTrho
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
%<a href="matlab: docsearchFS('harho')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

    % Obtain Figure 11.10 p. 375 of
    % Hoaglin et al. (1987)
    x=-9:0.1:9;
    rhoHA=HArho(x,1);
    plot(x,rhoHA)
    xlabel('x','Interpreter','Latex')
    ylabel(' Hampel $\rho(x) $','Interpreter','Latex')

%}


%{

    % Hampel rho function using a redescending slope of -1/3. 
    x=-9:0.1:9;
    rhoHA=HArho(x,[1,1.5,3.5,8]);
    plot(x,rhoHA)
    xlabel('x','Interpreter','Latex')
    ylabel(' Hampel $\rho(x) $','Interpreter','Latex')

%}

%% Beginning of code


if length(ctuning)>1,
    
    if ((ctuning(2) < 0) || (ctuning(3) < ctuning(2)) || (ctuning(4) < ctuning(3)))
        error('FSDA:HArho:WrongAbc',[' Illegal choice of parameters in Hampel: ' ...
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

rhoHA = ones(size(u));
absu=abs(u);


% 0.5* u^2,		   |u| <=a 
rhoHA(absu<=a) = 0.5*u(absu<=a).^2; 


% a/|u|,		 a <= |u| < b,
rhoHA(absu > a & absu <=b) = a*(abs(u(absu > a & absu <=b))) -0.5* a^2 ;

% ab-0.5a^2+0.5*(c-b)*a(1- ((c-|u|)/(c-b))^2 ),	  b <= |u| <  c,
rhoHA(absu>b & absu <=c) = a*b-0.5*a^2+0.5*(c-b)*a*(1-  ( (c-abs(u(absu>b & absu <=c)))/(c-b)  ).^2);

% ab-0.5a^2+0.5*(c-b)*a,			              |u| >= c.
rhoHA(absu > c) = a*b-0.5*a^2+0.5*(c-b)*a ;

end
