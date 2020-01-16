function c = bc(n,k)
%bc returns the Binomial coefficient
%
%<a href="matlab: docsearchFS('bc')">Link to the help page for this function</a>
%
% Required input arguments:
%
%       n:  Number of elements. Non negative integer. 
%           Data Types - single|double
%       k:  Items to choose from the set of n elements. Non negative integer.
%           Data Types - single|double
%
% Optional input arguments:
%
% Output:    
% 
%       c  : The binomial coefficient $n!/k!(n-k)!$. Integer. This is the 
%            coefficient of the $x^k$ term in the polynomial expansion of
%            the binomial power $(1 + x)^n$. This is also the so called
%            choose function of n and k (nchoosek in MATLAB), i.e. the
%            number of k-element subsets (the k-combinations) of a set of n
%            objects. When a coefficient is large, results may be inexact.
%            The result is only accurate to 15 digits for double-precision
%            inputs in 32bits computers.
%
% See also: nchoosek
%
% References:
%
%    Knuth, D. E. (1997). "The Art of Computer Programming", Volume 1:
%    Fundamental Algorithms, Third ed. Addison-Wesley. [pp. 52--74]. 
%
% Acknowledgements: 
%
% Matlab function bc has been adapted to this toolbox by FSDA team
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('bc')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %% Number of pairs chosen among 6 elements.
    bc(6,2)
%}

%% Beginning of code

% Ensure computations in doubles.
n = double(n);
k = double(k);

if k > n/2, k = n-k; end

nums = (n-k+1):n;
dens = 1:k;
nums = nums./dens;
c = round(prod(nums));

end
%FScategory:UTICOMB