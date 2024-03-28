function [logfn]=logfactorial(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluates the logarithm of n! with absolute accuracy of 1e-12 or better, and
% relative accuracy of 1e-15 or better. The following codes employ Stirling's integration
% formula of Gamma function.
% Input parameters:
% n - the positive integer to be evaluated;
% Output parameters:
% logfn - the logrithm value of n!
% Note:
%(1)if n<171, the value of log(n!) can be computed with high accuracy by Matlab itself;
%(2) n doesn't have to be an integer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  (n <= 50 )
    logfn=log(factorial(n));
elseif n<=170
    logfn = gammaln(n+1);
else
    eps=1e-18; % Preset accuracy for numerical evalueated log(n!) from Stirling's formula
    N=1e5; % Number of integration points in the numerical integral
    Lm=(1/2/pi)*log(1+1/(4*eps)); % Upper bound of the integral determined by the 'eps'
    t=Lm/N:Lm/N:(1-1/N)*Lm;
    g=2*atan(t./(n+1))./(exp(2*pi.*t)-1); % Integrand function in Stirling's integration formula
    % Multiplication coefficients of Simpson's open formula of fourth-order
    coe=ones(N-1,1);
    coe(1,1)=55/24;
    coe(2,1)=-1/6;
    coe(3,1)=11/8;
    coe(N-3,1)=11/8;
    coe(N-2,1)=-1/6;
    coe(N-1,1)=55/24;
    logfn=0.5*log(2*pi)+(n+0.5)*log(n+1)-n-1+g*coe*(Lm/N);
end
