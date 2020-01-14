function c = nchoosekFS(v,k)
%nchoosekFS returns the Binomial coefficient or matrix containing all combinations
%
%<a href="matlab: docsearchFS('nchoosekFS')">Link to the help function</a>
%
% Required input arguments:
%
%       v:  Vector of length n. Integer or array of non-negative integers. 
%           Data Types - single|double
%       k:  Items to choose from the set of n elements. Non negative integer.
%           Data Types - single|double
%
% Optional input arguments:
%
% Output:  
%
%     c:    scalar $v!/k!(v-k)!$ if $v$ and $k$ are non-negative integers
%           or matrix with $n!/k!(n-k)!$ rows and $k$ columns if $v$ is a
%           vector of length $n. Binomial coefficient(s) or all combinations.
%           Data Types - single|double
% 
% More About:
%
%   This function is similar to nchoosek of Statistics Toolbox but it is
%   much faster and makes a more efficient use of memory.
%
%   Returns the scalar $v!/k!(v-k)!$ if $v$ and $k$ are non-negative integers.
%   This is the number of combinations of $v$ things taken $k$ at a time. In
%   this case it makes use of function bc.
%
%   Produces a matrix with $n!/k!(n-k)!$ rows and $k$ columns if $v$ is a vector
%   of length $n%. Each row contains a combination of k elements taken
%   without repetitions among n. In this case function combsFS is used.
%
% See also: nchoosek, perms
%
% References:
%
%    Riordan, J. (1958), "An Introduction to Combinatorial Analysis", 
%    Wiley & Sons, New York.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('nchoosekFS')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:
%
%{
    %% Binomial coefficient(s) or all combinations.
    %  Profile generation of 2118760 combinations.
    v = 1:50; k = 4; 

    tic
    for i=1:10, nchoosekFS(v,k); end
    t_nchoosekFS = toc

    tic
    for i=1:10, nchoosek(v,k); end
    t_nchoosek = toc

    fprintf('nchoosekFS has been %5.2f times faster than nchoosek\n\n\n', t_nchoosek/t_nchoosekFS); 
    fprintf('Try now again using k=5: in a 32 bit computer\n');
    fprintf('nchoosekFS will require about the same time (in order of magnitude)\n');
    fprintf('while nchoosek will start swaping into virtual memory.\n'); 
%}

%% Beginning of code

if ~isscalar(k) || k < 0 || ~isreal(k) || k ~= round(k)
    error('FSDA:nchoosekFS:InvalidArg2',...
        'The second input has to be a non-negative integer.');
end

[m, n] = size(v);

if min(m,n) ~= 1
    error('FSDA:nchoosekFS:InvalidArg1',...
        'The first argument has to be a scalar or a vector.');
end

% the first argument is a scalar integer
if isscalar(v) && v >= 0
    if k > v
        error('FSDA:nchoosekFS:WrongK','K must be an integer between 0 and N.');
    end
    
    % if the first argument is a scalar, then, we only return the number of
    % combinations. Not the actual combinations.
    c=bc(v,k);
else
    % the first argument is a vector, generate actual combinations.
    c=combsFS(v,k);
end

end

%FScategory:UTICOMB