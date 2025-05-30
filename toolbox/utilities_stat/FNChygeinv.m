function x = FNChygeinv(p,M,K,n,odds, accuracy)
%FNChygeinv computes the inverse of the Fisher non central hypergeometric cumulative distribution function (cdf).
%
%   Since the underlying distribution is discrete, FNChygeinv returns the
%   smallest integer x, such that the hypergeometric cdf evaluated at x,
%   equals or exceeds p. The size of x is the common size of the input
%   arguments. A scalar input functions as a constant matrix of the same
%   size as the other inputs.
%
%<a href="matlab: docsearchFS('FNChygeinv')">Link to the help page for this function</a>
%
%
%  Required input arguments:
%
%
%           p    : input probabilitiies. Scalar or vector or matrix. You
%                  can think of p as the probability of observing x
%                  defective items (red balls) in n drawings without
%                  replacement from a group of M items where K are
%                  defective (the total numnber of red balls is K) and the
%                  ratio of the probability of observing a defect with that
%                  of observing a non defect is equal to odds.
%                  Data Types -  single|double
%           M    : Total number of balls in urn before sampling. Scalar.
%                  Data Types - single|double
%           K    : Initial number of red balls in the urn. Scalar.
%                  Data Types - single|double
%           n    : Total number of balls sampled. Scalar.
%                  Data Types - single|double
%           odds : Probability ratio of red over white balls. Scalar.
%                  Data Types - single|double
%
%  Optional input arguments:
%
%
%       accuracy : accuracy of the calculations. Scalar. The default value
%                  of accuracy is 1e-08.
%                  Example - 1e-06
%                  Data Types - single|double
%
%
%
%  Output:
%
%
%           x : Fisher' quantile values.  Quantiles corresponding to input probabilities. 
%                  The size of x is the common size of the input
%                  arguments. A scalar input functions as a constant matrix
%                  of the same size as the other inputs.
%
%
%
% See also: FNChygepdf, FNChygecdf, FNChygeinv, FNChygernd, WNChygepdf, WNChygecdf, WNChygeinv, WNChygernd
%
% References:
%
% Fog, A. (2008), Calculation Methods for Wallenius' Noncentral Hypergeometric
% Distribution, "Communications in Statistics - Simulation and
% Computation", Vol. 37, pp. 258-273.
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('FNChygeinv')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Compute the inverse of Fisher non central hypergeometric distribution.
    % M = total number of balls;
    M=80; 
    n=10;  % number of balls taken
    odds=3; % Prob. of red balls vs other color balls
    K=50; % Number of red balls in the urn
    % Compute quantile 0.3
    x030=FNChygeinv(0.3,M,K,n,odds);
    disp(x030)
%}

%{
    %% Compute quantiles 0.1, 0.2, ..., 0.9 from nverse of Fisher non central hypergeometric distribution.
    % M= total number of balls
    M=80;
    % n= number of balls taken
    n=10;  
    % odds = Ratio of Prob. of red balls vs other color balls    
    odds=3; 
    % K = number of red balls in the urn
    K=50; 
    % Compute quantiles 0.1, ..., 0.9
    xquant=FNChygeinv(0.1:0.1:0.9,M,K,n,odds);
    disp(xquant)
%}


%% Beginning of code
if nargin < 5
    error(message('FSDA:FNChygeinv:TooFewInputs'));
end

if nargin<6
    accuracy=1e-08;
end

 [errorcode, p, M, K, n, odds] = distchck(5,p,M,K,n,odds);
 
 if errorcode > 0
     error(message('FSDA:FNChygeinv:InputSizeMismatch'));
 end

% Initialize X to zero.
outType = internal.stats.dominantType(p,M,K,n,odds);
x = zeros(size(p),"like",outType);
x(isnan(p) | isnan(M) | isnan(K) | isnan(n)) = NaN;

%   Return NaN for values of the parameters outside their respective limits.
k1 = M < 0 | K < 0 | n < 0 | round(M) ~= M | round(K) ~= K ...
    | round(n) ~= n | n > M | K > M | p < 0 | p > 1 | isnan(p);
if any(k1(:))
    x(k1) = NaN;
end

cumdist = FNChygepdf(x,M,K,n,odds, accuracy);
count = zeros(size(p));

% Compare P to the hypergeometric distribution for each value of N.
while any(p(:) > cumdist(:)) && count(1) < max(n(:)) && count(1) < max(K(:))
    count = count + 1;
    idx = find(cumdist < p - eps(p));
    x(idx) = x(idx) + 1;
    cumdist(idx) = cumdist(idx) + FNChygepdf(count(idx),M(idx),K(idx),n(idx),odds(idx),accuracy);
end

% Check using hygecdf
y = FNChygecdf(x,M,K,n,odds);
ynew = zeros(size(y), "like", y);
xnew = x;
under = y<p & ~k1;
while any(under(:))
    ynew(under) = FNChygecdf(xnew(under)+1,M(under),K(under),n(under),odds(under),accuracy);
    xnew(under) = xnew(under)+1;
    under = under & ynew<p;
end
x = xnew;
ynew = zeros(size(y), "like", y);
xnew = x;
over = y>p & ~k1 & ~under;
while any(over(:))
    ynew(over) = FNChygecdf(xnew(over)-1,M(over),K(over),n(over),odds(over), accuracy);
    over = over & ynew>=p;
    xnew(over) = xnew(over)-1;
end
x = xnew;
end


%FScategory:ProbDist