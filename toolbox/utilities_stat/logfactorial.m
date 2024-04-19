function [logfn] = logfactorial(n,method)
%logfactorial returns the logarithm of the factorial
%
%<a href="matlab: docsearchFS('logfactorial')">Link to the help function</a>
%
% Factorials grow very quickly and results exceed the limits of computer
% floating point arithmetic already for 200! It is therefore more
% convenient to work with the logarithms of factorials. $log(n!) = log(1) +
% log(2) + log(3) + \ldots + log(n)$ would do the work but it is too slow
% to iterate the various log computations. There are alternatives to this
% naive approach. The first is to express the factorial in terms of the
% gamma function: $n! = \Gamma(n+1)$. The second is to use Stirling's
% approximation: $log \Gamma(n) \approx (n – 1/2) \log(n) – n  + (1/2)
% log(2 \pi) + 1/(12 n) – 1/(360 n^3) + 1/(1260 n^5) – \ldots$ which uses
% Bernoulli's numbers.
%
%  Required input arguments:
%
%           n    : nonnegative (integer) value. Scalar. It doesn't have to 
%                  be necessarily an integer.
%                  Data Types - single|double
%
%  Optional input arguments:
%
%     method     : calculation method. Integer. 
%                  Example - 'method',3
%                  Data Types - single | double
%
%
%  Output:
%
%           logfn : the logarithm value of n!. Scalar. 
%                  Data Types - single|double
%
% See also factorial.m
%
% References:
%
% Li, Y.C. (2006), A note on an identity of the gamma function and Stirling
% formula, "Real Analysis Exchange", Vol. 32 (1), pp.267–271.
% Fog, A. (2008), Calculation Methods for Wallenius' Noncentral Hypergeometric
% Distribution, "Communications in Statistics - Simulation and
% Computation", Vol. 37, pp. 258-273.
% Andy Huang (2024). Log Factorial of Large Positive Numbers 
% (https://www.mathworks.com/matlabcentral/fileexchange/
% 33687-log-factorial-of-large-positive-numbers), 
% MATLAB Central File Exchange. Retrieved April 3, 2024.
%
% Copyright 2008-2024.
%
%<a href="matlab: docsearchFS('logfactorial')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%
%{
    % logfactorial with default method: result is OK.
    lf250   = logfactorial(250)
    % logfactorial forced with the naive log(factorial(n)): result is inf.
    lf250m1 = logfactorial(250,1)
%}

%{
    % A test to check the relative accuracy between the methods. 
    % The naive one is clearly Ok only for n<170.
    i=0;
    for n=10:1000
        i=i+1;
        o1(i) = logfactorial(n,1);
        o2(i) = logfactorial(n,2);
        o3(i) = logfactorial(n,3);
        o4(i) = logfactorial(n,4);
    end
    o1(o1==inf) = -1;
    subplot(3,1,1)
    plot(o1-o4,'+')
    title('m1 - m4 -- log(factorial(n)) accurate till n~170');
    subplot(3,1,2)
    plot(o2-o4,'o')
    title('m2 - m4');
    subplot(3,1,3)
    plot(o3-o4,'*')
    title('m3 - m4');
%}

if nargin<2
    if (n <= 50) && (n==round(n))
        % uses factorial function of MATLAB. In principle it should be
        % accurate until n<171, but we prefer using it for smaller n. 
        method = 1;
    elseif n <= 170
        % uses gammaln function. 
        method = 2;
    else
        % uses Stirling's approximation: accuracy is in the order of 1e-12
        method = 3;
    end
end

switch method
    case 1
        % obvious method. 
        logfn=log(factorial(n));
    case 2
        % approximation using gammaln. Note that the gamma function is 
        % defined for all complex numbers, other than non-positive integers.
        logfn = gammaln(n+1);
    case 3
        % Stirling's approximation
        eps = 1e-18;  % Preset accuracy for numerical evaluated log(n!) from Stirling's formula
        N   = 1e5;    % Number of integration points in the numerical integral
        Lm  = (1/2/pi)*log(1+1/(4*eps)); % Upper bound of the integral determined by the 'eps'
        t   = Lm/N:Lm/N:(1-1/N)*Lm;
        g   = 2*atan(t./(n+1))./(exp(2*pi.*t)-1); % Integrand function in Stirling's integration formula
        % Multiplication coefficients of Simpson's open formula of fourth-order
        coe         = ones(N-1,1);
        coe(1,1)    = 55/24;
        coe(2,1)    = -1/6;
        coe(3,1)    = 11/8;
        coe(N-3,1)  = 11/8;
        coe(N-2,1)  = -1/6;
        coe(N-1,1)  = 55/24;
        logfn       = 0.5*log(2*pi)+(n+0.5)*log(n+1)-n-1+g*coe*(Lm/N);

    case 4
        % Stirling's approximation to two orders
        D  = 1;
        C0 =  0.918938533204672722; % ln(sqrt(2*pi))
        C1 =  1/12;
        C3 = -1/360;
        C5 =  1/1260;
        C7 = -1/1680;
        if (n < 6)
            if (n == 0 || n == 1)
                logfn=0;
                return
            else
                while (n < 6)
                    D = D*n;
                end
            end
        end
        r     = 1. / n;  r2 = r*r;
        logfn = (n + 0.5)*log(n) - n + C0 + r*(C1 + r2*(C3 + r2*(C5 + r2*C7)));
        if D ~= 1
            logfn = logfn - log(D);
        end
end
end

%FScategory:UTISTAT