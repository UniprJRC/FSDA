function psiHAx = HApsix(u, ctuning,varargin)
%HApsix computes psi function  using Hampel proposal times x
%
%<a href="matlab: docsearchFS('HApsix')">Link to the help function</a>
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
%                biweight with parameter 8.
%
% Function HApsix transforms vector u as follows
%
% HApsi(u) = 	{ u^2,			                               |u| <= a,
%		        { a*sign(u)*u,		                      a <= |u| < b,
%		        { a((c-|u|)/(c-b))*sign(u)*u,	                 b <= |u| <  c,
%		        { 0,			                                |u| >= c.
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
%<a href="matlab: docsearchFS('hapsix')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

    % Obtain bottom panel of Figure 11.10 p. 375 of
    % Hoaglin et al. (1987)
    x=-9:0.1:9;
    psiHAx=HApsix(x,1);
    plot(x,psiHAx)
    xlabel('x','Interpreter','Latex')
    ylabel(' Hampel $\psi(x) \times x $','Interpreter','Latex')

%}

%% Beginning of code

if length(ctuning)>1,
    
    if ((ctuning(2) < 0) || (ctuning(3) < ctuning(2)) || (ctuning(4) < ctuning(3)))
        error('FSDA:HApsix:WrongAbc',[' Illegal choice of parameters in Hampel: ' ...
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


psiHAx = zeros(size(u));
absu=abs(u);


% u^2,		   |u| <=a 
psiHAx(absu<=a) = u(absu<=a).^2; 


% a*sign(u),		                      a <= |u| < b,
psiHAx(absu > a & absu <=b) = a*sign(u(absu > a & absu <=b)).*u(absu > a & absu <=b);

%	a((c-|u|)/(c-b))*sign(u),	                 b <= |u| <  c,
psiHAx(absu>b & absu <=c) = a*(  (c-abs(u(absu>b & absu <=c)))/(c-b)  ).*sign(u(absu>b & absu <=c)).*u(absu > b & absu <=c);

% 0,			              |u| >= c.
end
