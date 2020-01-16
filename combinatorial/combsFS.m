function P = combsFS(v,m)
%combsFS is an iterative algorithm equivalent to the MATLAB combs.m
%
% It generates m-combinations without repetition taken in lexicographic
% order from the vector v.
%
% REMARK: the MATLAB function combs.m uses recursive calls and it is
% therefore very inefficient. Our iterative counterpart also makes better
% use of memory, first because it works iteratively, and then because we
% force computations in the lowest possible precision. This is not a
% limitation, because the algotithm first builds the matrix P of all
% m-combinations starting from the first n natural numbers, for which
% double precision is not at all needed. Then, if the input vector b is
% different from vector 1:v, then the desired P is simply obtained as P =
% v(P). Note also that we build the matrix P by going over colums rather
% than over lines. This is faster, as MATLAB indexes the elements of a
% matrix by column first.
%
%<a href="matlab: docsearchFS('combsFS')">Link to the help function</a>
%
%  Required input arguments:
%
%    v:         A vector with n elements. It contains the response variable.
%               It can be either a row or a column vector.
%               Data Types - single|double
%    m:         Scalar. It specifies the size of the combinations.
%               Data Types - single|double
%
% Optional input arguments:
%
% Output:
%
%     P:        m-combinations without repetition taken in lexicographic
%               order from the vector v. Matrix containing the
%               m-combinations in the rows.
%               Data Types - single|double
%
% See also: nchoosek
%
% References:
%
%    Knuth, D. E. (1997). "The Art of Computer Programming", Volume 1:
%    Fundamental Algorithms, Third ed. Addison-Wesley. [pp. 52--74]. 
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('combsFS')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:

%{
    %% combsFS used to generate all possible combinations of size 3 of elements 5, 8, 9, 10, 11.
    combsFS([5 8:11],3)
%}


%% Beginning of code
if nargin~=2
    error('FSDA:combsFS:WrongInputNum', 'Requires 2 input arguments.')
end

v = v(:).';     % Make sure v is a row vector.
n = length(v);  % Elements of v.

% set the *theoretical* precision based on the number of elements in v.
% Of course so big n values will never be used in practice. The precision
% will always be int8.
if n < 128
    precision = 'int8';
    v = int8(v);
elseif n < 32768
    precision = 'int16';
    v = int16(v);
else
    precision = 'int32';
    v = int32(v);
end

if(m > n)
    error('FSDA:combsFS:WrongInputNum', 'm > n !!');
elseif n == m
    P = v;
elseif m == 1
    P = v.';
elseif(m == 0)
    P=[];
else
    %The binomial coefficient (n choose m) can be computed using
    %prod(np1-m:n)/prod(1:m). For large number of combinations
    %our function 'bc' is better.
    bcn = bc(n,m);
    % initialise the matrix of all m-combinations
    P = zeros(bcn,m,precision);
    np1 = n + 1;  % do once here n+1 (needed in the internal loop)
    toRow = np1 - m;
    % set the first n+1-m rows of the last column
    P(1:toRow , m) = (m:n)';
    for i = m-1:-1:1                     % external loop over colums
        [s1 , s2] = deal(toRow);
        % set the first n+1-m rows block of rows of colum i
        P(1:toRow , i) = i;
        for j = i+1:i+n-m                 % internal loop
            s1      = s1*(np1+i-j-m)/(np1-j);
            fromRow = toRow + 1;
            toRow   = fromRow + s1 - 1;
            P(fromRow:toRow , i+1:m) = P(s2-s1+1:s2 , i+1:m);
            P(fromRow:toRow , i)     = j;
        end
    end
    % find the true P if the vector of elements in v is not the default 1:n
    if ~isequal(v,1:n) , P = v(P); end
end

end
%FScategory:UTICOMB
