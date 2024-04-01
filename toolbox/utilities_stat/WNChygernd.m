function r = WNChygernd(M,K,n,odds,   mm,nn,oo)
%WNChygernd random arrays from the Wallenius non central hypergeometric distribution.
% 
%   returns an array of random numbers of size mm-by-nn-by-oo chosen from the 
%   Wallenius non central hypergeometric distribution with parameters M, K, n, and odds
%
%<a href="matlab: docsearchFS('WNChygernd')">Link to the help page for this function</a>
%
%  Required input arguments:
%
%           M    : Total number of balls in urn before sampling. Scalar.
%                  Data Types - single|double
%           K    : Initial number of red balls in the urn. Scalar.
%                  Data Types - single|double
%           n    : Total number of balls sampled. Scalar.
%                  Data Types - single|double
%        odds    : Probability ratio of red over white balls. Scalar.
%                  Data Types - single|double
%
% Optional input arguments:
%
%         mm    : Length of first dimension. Scalar. Number of rows of the
%                 array which contains the random numbers. 
%               Example - 3
%               Data Types - double
%         nn    : Length of second dimension. Scalar. Number of columns of
%                 the array which contains the random numbers.
%               Example - 2
%               Data Types - double
%         oo    : Length of third dimension. Scalar. Number of 3D slides of
%                 the array which contains the random numbers.
%               Example - 5
%               Data Types - double
%       accuracy : accuracy of the calculations. Scalar. The default value
%                  of accuracy is 1e-10.
%                  Data Types - single|double
%                  Example - 1e-06
%
% Output:
%
%         r    : Random numnbers. Array of random numbers from the
%                Wallenius non central hypergeometric distribution.
%                The size of rr is determined by the optional input
%                parameters mm, nn, oo.
%
%
% See also: WNChygepdf, WNChygecdf, WNChygeinv
%
% References:
%
% Fog, A. (2008), Calculation Methods for Wallenius' Noncentral Hypergeometric
% Distribution, "Communications in Statistics - Simulation and
% Computation", Vol. 37, pp. 258-273.
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('WNChygernd')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Generate a random number from  Wallenius non central hypergeometric distribution.
    M=80; % total number of balls
    n=10;  % number of balls taken
    odds=3; % Prob. of red balls vs other color balls
    K=50; % Number of red balls in the urn
    samplesize=100000; % Sample size to extract
    edges=(0:(n+1))';
    % Generate a random number from the distribution above
    x=WNChygernd(M,K,n,odds);
    disp(x)
%}

%{
    % A difficult example which needs to avoid underflow/overflow.
    M=8000; % total number of balls
    n=2000;  % number of balls taken
    odds=3; % Prob. of red balls vs other color balls
    K=5000; % Number of red balls in the urn
    % Generate an array of size 2x3x5 of random numbers from the distribution above
    mm=2;   nn=3; oo=5;
    rng(12345)
    X=WNChygernd(M,K,n,odds);
    disp(X)
%}


%{
    % Generate a 3D array of size mmxnnxoo of random number from  Wallenius non central hypergeometric distribution.
    M=800; % total number of balls
    n=200;  % number of balls taken
    odds=3; % Prob. of red balls vs other color balls
    K=500; % Number of red balls in the urn
    % Generate an array of size 2x3x5 of random numbers from the distribution above
    mm=2;   nn=3; oo=5;
    X=WNChygernd(M,K,n,odds, mm,nn,oo);
    disp(X)
%}

%{
    %% Comparison density and random numbers.
    close all
    M=100; % total number of balls
    n=10;  % number of balls taken
    odds=2; % Prob. of red balls vs other color balls
    K=50; % Number of red balls in the urn
    samplesize=100000; % Sample size to extract
    edges=(0:(n+1))';
    % Compute the density
    Wpdf=WNChygepdf(edges(1:end-1),M,K,n,odds);
    % Generate random numbers
    x=WNChygernd(M,K,n,odds,samplesize);
    % bar plot of theoretical and relative frequencies
    freqANDdens=[(histcounts(x,edges-0.5)/samplesize)' Wpdf];
    bar(edges(1:end-1),freqANDdens)
    legend(["Theoretical density" "Empirical relative frequency"],'Location','best')
    xlabel('Number of successes')
%}



%% Beginning of code
if nargin < 4
    error(message('FSDA:WNChygernd:TooFewInputs'));
end

if nargin<8
    accuracy=1e-08;
end

if nargin<7
    oo=1;
end
if nargin <6
    nn=1;
end
if nargin <5
    mm=1;
end



r=zeros(mm,nn,oo);
for i=1:numel(r)
    r(i)=aux.WalleniusNCHypergeometricrnd(n,K,M,odds, accuracy);
end




