function x = shuffling(x)
%shuffling does a random permutation of the elements of input vector
%
%<a href="matlab: docsearchFS('shuffling')">Link to the help function</a>
%
%  Required input arguments:
%
%  x :  A set of elements. Vector of length t. 
%       Data Types - single|double
%
%  Optional input arguments:
%
%  Output: 
%
%  x :  A permutation of the set of elements in x. Vector of length t. 
%       Data Types - single|double
%
% More About:
%
% If set $x$ has $t$ elements, the objective is to obtain each of the $t!$
% pemutations with equal probability, especially when $t$ is large. To
% achieve this goal, we use backward Knuth's shuffling, which is actually
% based on the Fisher-Yates shuffle. Once compiled, Knuth solution is more
% efficient than the natural MATLAB solution x(randperm(numel(x)).
%
% See also: randperm.
%
% References:
%
% Knuth, Donald E. (1969). The Art of Computer Programming volume 2,
% Seminumerical algorithms, Reading, MA: Addison-Wesley, pp. 124-125.
%
% Fisher, R.A.; Yates, F. (1948) [1938]. Statistical tables for biological,
% agricultural and medical research (3rd ed.). London, Oliver & Boyd. pp
% 26-27.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('shuffling')">Link to the help function</a>
%
% Last modified 06-Feb-2015

% Examples:
%{
    %% shuffling applied to a set of 20 elements.
    shuffling(1:20)
%}

%{
    % shuffling applied to a set of 20 elements, but using a parsimonious
    data % type; this is convenient if the vector is big.
    shuffling(int8(1:20))
%}

%{
    % check of the permutation produced by shuffling
    x = 1:200000;
    numel(unique(shuffling(x)))
%}

%{
    % To profile against randperm, which uses sort.
    for i=1:5000
        x = randi(100000,10000,1);
        nn=numel(x);
        shuffling(x);
        x(randperm(nn));
    end
%}

%% Beginning of code
%if iscolumn(x) ; x = x' ; end;
n = numel(x);
I = n:-1:1;
J = ceil(rand(1,n) .* I);

for i=n:-1:1 
    % Note 1: to generate random integers between 1 and i, one by one in
    % the loop, would be less efficient. Such numbers are thus generated
    % outside the loop. To check, uncomment line below and profile the code.
    % j = ceil(rand() * i);
    % Note 2: this would be even less efficient: 
    % j = randi([1,i]);
    % Note 3: the double access to the vector on indices takes less time
    % than the allocation of the index in a separate variable with this line
    % j = J(i);
    
    % now, exchange element x[j] with x[i]. In MATLAB it would be
    % natural to do the exchange with a matrix operation, i.e. with 
    % x([i,j]) = x([j,i]).
    % Note that this switch forces MATLAB to pre-allocate a copy of x,
    % which takes time and, of corse, more memory. Using the MATLAB
    % profiler we found that the switch requires about 72% additional
    % time with respect to the following solution:
    t       = x(J(i));
    x(J(i)) = x(i);
    x(i)    = t;
end

end

