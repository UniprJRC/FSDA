function [out]= genr8(n, distrib, s1, s2)
%genr8 returns a vector of pseudorandom number sequence.
%
%<a href="matlab: docsearchFS('genr8')">Link to the help function</a>
%
%
%    This function generates a sequence of pseudorandom numbers
%    using the 32-bit rnadom number generator from figure 3 of
%    the article by Pierre L'Ecuyer.
%    The cycle length is claimed to be 2.30584E+18.
%
%  Required input arguments:
%
%
%   n       :    Number of random numbers. Positive integer. The number of
%                random numbers to generate.
%                Data Types - integer value
%
%
%
%  Optional input arguments:
%
%   distrib :     Distribution. Scalar 0 or 1. Default is distrib=0 that
%                   indicates the use of the uniform distribution in [0 1],
%                   while distrib=1 indicates the use of the standard normal
%                   with 0 mean and standard deviation 1
%                   Data Types - integer value
%   s1:       First value used as the seed for the sequence. Scalar.
%                   Whenever needed, the user can initialize s1 to a value
%                   between 1 and 2147483562.
%                   Example - 4356123
%                   Data Types - integer value
%   s2:       Second value used as the seed for the sequence. Scalar.
%                   Whenever needed, s2 can be initialized to
%                   a value between 1 and 2147483398.
%                   Example - 123
%                   Data Types - integer value
%
%
%
% Output:
%
%         out: random sequence. Vector. n x 1 vector containing the
%               generated random sequence.
%
% More About:
%
%   The Mersenne Twister algorithm is the most used RNG among statistical
%   software, but most implementations contain slightly variations that
%   make difficult to obtain random numbers sequences that are congruent
%   across platforms. This simple yet powerful algorithm is a variation of
%   the well known Wichmann Hill RNG but boasts an impressive cycle length
%   of 2.30584E+18. Given the sheer but clear structure and the lack of
%   specific mathematical operators, the coding should be easily
%   reproducible across all platforms.
%
% See also: mtR
%
%
% References:
%
% L'Ecuyer, P. (1988), Efficient and Portable Combined Random Number Generators,
%    "Communications of the ACM", Vol. 31, pp. 742-751.
% Wichmann, B. and Hill, D. (1982), An Efficient and Portable Pseudo-random
%    Number Generator, "Applied Statistics", Vol. 31, pp. 188-190.
%
% Acknowledgements:
%
%    Original PASCAL version by Pierre L'Ecuyer.
%    Modifications by John Burkardt.
%    Further modifications by FSDA team
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('genr8')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit
%
% Examples:

%{
    % genr8 using all default arguments.
    out = genr8(1)
%}

%{
    % genr8 specifing the type of distribution.
    
    % a vector of 5 uniform distributed random numbers
    out = genr8(5, 0)
    % a vector of 10 normal distributed random numbers
    out = genr8(10, 1)

%}

%{
    % genr8 specifing the type of distribution and 1 seed.
    
    % vector of 5 uniform distributed random numbers and 1 seed
    out = genr8(5, 0, 12345)
    % vector of 50 normal distributed random numbers and 1 seed
     out = genr8(10, 1, 45678)

%}

%{
    % genr8 specifing the type of distribution and the 2 seeds.
    
    % vector of 5 uniform distributed random numbers with both seeds
    out = genr8(5, 0, 12345, 45678)
    % vector of 50 normal distributed random numbers with both seeds
     out = genr8(10, 1, 12345, 45678)

%}


%% Beginning of code
if nargin <1 || isempty(n)
    error('FSDA:gen_r8uni:MissingInput','n must be specified');
end

randvec=zeros(n,1);

persistent seed1;
persistent seed2;

if nargin < 2 || isempty(distrib)
    % user supplied only n
    distrib=0;
    s1=12345;
    s2=56789;
    
elseif nargin < 3
    % user supplied n and distrib
    if ( ~isempty(seed1) ||  ~isempty(seed2))
        % not the first run
        s1=seed1;
        s2=seed2;
    else
        % first run, default seeds
        s1=12345;
        s2=56789;
    end
elseif nargin < 4
    % user supplied n, distrib and only one seed:
    % the seed supplied is assigned to seed1
    if ( ~isempty(seed1) ||  ~isempty(seed2))
        % not the first run
        s2=seed2;
    else
        % first run, default seed2
        s2=56789;
    end
elseif  nargin < 5
    % user supplied n, distrib, s1, s2
    % nothing to do
end


for i = 1 : n
    [ randvec(i), s1, s2 ] = r8_uni ( s1, s2 );
end


if distrib==1
    randvec=norminv(randvec);
end

seed1=s1;
seed2=s2;
out=randvec;
end

function [ r, s1, s2 ] = r8_uni ( s1, s2 )

k = floor ( s1 / 53668 );
s1 = 40014 * ( s1 - k * 53668 ) - k * 12211;
if ( s1 < 0 )
    s1 = s1 + 2147483563;
end

k = floor ( s2 / 52774 );
s2 = 40692 * ( s2 - k * 52774 ) - k * 3791;
if ( s2 < 0 )
    s2 = s2 + 2147483399;
end

z = s1 - s2;
if ( z < 1 )
    z = z + 2147483563;
end

r = z / 2147483563.0;

end
%FScategory:UTISTAT