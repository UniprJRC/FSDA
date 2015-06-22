function c = bc(n,k)
%bc returns the Binomial coefficient
%
%<a href="matlab: docsearchFS('bc')">Link to the help page for this function</a>
%
% Required input arguments:
%
%       n:  a non negative integer
%       k:  a non negative integer
%
% Output:
%
%  bc(n,k) = where n and k are non-negative integers returns n!/k!(n-k)!,
%            i.e. the coefficient of the x^k term in the polynomial
%            expansion of the binomial power (1 + x)^n. This is also the so
%            called choose function of n and k (nchoosek in MATLAB), i.e.
%            the number of k-element subsets (the k-combinations) of a set
%            of n objects. When a coefficient is large, results may be
%            inexact. The result is only accurate to 15 digits for
%            double-precision inputs in 32bits computers.
%
% See also: nchoosek
%
% Copyright 2008-2015.
% Matlab function bc has been adapted to this toolbox by FSDA team
%
%<a href="matlab: docsearchFS('bc')">Link to the help page for this function</a>
%
% Last modified 06-Feb-2015
%
% Examples:
%{
    %Number of pairs chosen among 6 elements (returns 15)
    bc(6,2)
%}

% Ensure computations in doubles.
n = double(n);
k = double(k);

if k > n/2, k = n-k; end

nums = (n-k+1):n;
dens = 1:k;
nums = nums./dens;
c = round(prod(nums));

