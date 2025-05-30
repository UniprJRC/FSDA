function x = tobitinv(p, mu, sigma, left, right)
%tobitinv computes the inverse of the tobit cumulative distribution function.
%
%<a href="matlab: docsearchFS('tobitinv')">Link to the help function</a>
%
%  Required input arguments:
%
%    p:         Probability at which the inverse of the cdf must be evaluated
%               $0 \leq p \leq 1$.
%               Scalar, vector or matrix 3D array of the same size of x and b.
%               A scalar input functions as a constant matrix of the same
%               size as the other input.
%               See the section called "More About:" for more details about
%               the inverse gamma distribution.
%               Data Types - single | double
%
%
%  Optional input arguments:
%
%    mu :       location parameter of the tobit distribution.
%               Scalar, vector or matrix 3D array of the same size of x and sigma, Lower, Upper.
%               A scalar input functions as a constant matrix of the same
%               size as the other input. Default value of mu is 0.
%               See "More About:" for details about the tobit
%               distribution.
%               Example - 'mu',10
%               Data Types - single | double
%
%    sigma :    scale parameter of the tobit distribution.
%               Scalar, vector or matrix 3D array of the same size of x and sigma, Lower, Upper.
%               A scalar input functions as a constant matrix of the same
%               size as the other input. Default value of sigma is 1
%               See "More About:" for details about the tobit
%               distribution.
%               Example - 'sigma',800
%               Data Types - single | double
%
%    left :     lower limit for the censored random variable. Scalar.
%               If set to -Inf, the random variable is assumed to be not
%               left-censored; default value of left is zero (classical
%               tobit model).
%               Example - 'left',1
%               Data Types - double
%
%    right :    right limit for the censored random variable. Scalar.
%               If set to Inf, the random variable is assumed to be not
%               right-censored; default value of left is Inf (classical
%               tobit model).
%               Example - 'right',800
%               Data Types - double%
%
%  Output:
%
%    x:         inverse CDF value. Scalar, vector or matrix or 3D array of the same size
%               of input arguments p, mu, sigma, left, right. $p=\int_0^x f_{tobit}(t | \mu, \sigma, left, right) dt$ is the
%               inverse of tobit cdf with parameters mu, sigma, left, right
%               for the corresponding
%               probabilities in p.
%
% More About:
%
%
% The cdf of the tobit distribution defined over the support
% $x \in R $ with location parameter $mu$, scale parameter $sigma$ and
% lower and upper censor values left and right
%  \[
%  F_{tobit}(x, \mu, \sigma, left, right)  =\int_0^x f(t) dt
%  \]
%
%
% See also: tobitpdf, tobitcdf, tobitrnd
%
% References:
%
% Greene, W.H. (2008), "Econometric Analysis, Sixth Edition", Prentice Hall, pp. 871-875.
%
% Tobin, J. (1958), Estimation of Relationships for Limited Dependent
% Variables, "Econometrica", 26, pp. 24-36.
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('tobitinv')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples:
%

%{
 % Using all default options.
 % In this case we assume that mu=0, sigma=1, left=0, right=Inf
    x=tobitinv(0.5)
%}

%{
 % mu is specified.
 % In this case we assume that sigma=1, left=0, right=Inf
    x=tobitinv(0.8,1)
%}

%{
 % mu sigma are specified.
 % In this case we assume that left=0, right=Inf
    x=tobitinv(0.8,1,0.4)
%}

%{
 % mu sigma and left are specified.
 % In this case we assume that right=Inf
    x=tobitinv(0.95,2.5,0.4,3)
%}


%{
    %% Check accuracy of results, monitoring $|x-F_{tobit}^{-1} (F_{tobit}(x))|$.
    a=100; 
    b=200;
    mu=(a+b).*0.6;
    sigma=40;
    x=tobitrnd(mu,sigma,a,b,100,1);
    Y=zeros(length(x),1);
    Ychk=Y;

    for i=1:length(x)
        Y(i)=x(i)-tobitinv(tobitcdf(x(i),mu,sigma,a,b),mu,sigma,a,b);
    end
    disp('Maximum deviation from 0');
    disp(max(max(abs(Y))));
%}


%% Beginning of code

if nargin<5
    right =  Inf;
end

if nargin<4
    left =  0;
end

if nargin<3
    sigma =  1;
end

if nargin<2
    mu =  0;
end

if nargin < 1
    error(message('FSDA:tobitinv:TooFewInputs'));
end


x=norminv(p,mu,sigma);
    xleft=x<left;
    
x(xleft)=left(xleft);
    
xright=x>right;

x(xright)=right(xright);
 
end
%FScategory:ProbDist