function r = FNChygernd(M,K,n,odds,   mm,nn,oo)
%FNChygernd random arrays from the Fisher non central hypergeometric distribution.
% 
%   returns an array of random numbers of size mm-by-nn-by-oo chosen from the 
%   Fisher non central hypergeometric distribution with parameters M, K, n, and odds
%
%<a href="matlab: docsearchFS('FNChygernd')">Link to the help page for this function</a>
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
%                 array which contains the random numbers 
%               Example - 3
%               Data Types - double
%         nn    : Length of second dimension. Scalar. Number of columns of
%                 the array which contains the random numbers.
%               Example - 2
%               Data Types - double
%         oo    : Length of third dimension. Scalar. Number of 3D slides of
%                 the array which contains the random numbers
%               Example - 5
%               Data Types - double
%
% Output:
%
%         r    : Random numnbers. Array of random numbers from the
%                Fisher non central hypergeometric distribution.
%                The size of rr is determined by the optional input
%                parameters mm, nn, oo.
%
% See also: FNChygepdf, FNChygecdf, FNChygeinv, WNChygernd, WNChygepdf, WNChygecdf, WNChygeinv
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
%<a href="matlab: docsearchFS('FNChygernd')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Generate a random number from  Wallenius non central hypergeometric distribution.
    % M = total number of balls.
    M=80;
    % n = number of balls taken 
    n=10;  
    % odds= prob. of red balls vs other color balls.
    odds=3; 
    % Number of red balls in the urn.
    K=50; 
    % Generate a random number from the distribution above
    x=FNChygernd(M,K,n,odds);
    disp(x)
%}

%{
    % A difficult example which needs to avoid underflow/overflow.
    % M is total number of balls.
    M=8000; 
    n=2000;  % number of balls taken
    odds=3; % Prob. of red balls vs other color balls
    K=5000; % Number of red balls in the urn
    % Generate an array of size 2x3x5 of random numbers from the distribution above
    mm=2;   nn=3; oo=5;
    rng(12345)
    X=FNChygernd(M,K,n,odds);
    disp(X)
%}

%{
    % Generate a matrix of size mmxnn of random number from  Fisher non central hypergeometric distribution.
    % M is total number of balls.
    M=80; 
    n=10;  % number of balls taken
    odds=3; % Prob. of red balls vs other color balls
    K=50; % Number of red balls in the urn
    % Generate a matrix of size 3x5 of random numbers from the distribution above
    mm=3;
    nn=5;
    X=FNChygernd(M,K,n,odds, mm,nn);
    disp(X)
%}

%{
    % Generate a 3D array of size mmxnnxoo of random number from  Fisher non central hypergeometric distribution.
    % M is total number of balls.
    M=800;
    n=200;  % number of balls taken
    odds=3; % Prob. of red balls vs other color balls
    K=500; % Number of red balls in the urn
    % Generate an array of size 2x3x5 of random numbers from the distribution above
    mm=2;   nn=3; oo=5;
    X=FNChygernd(M,K,n,odds, mm,nn,oo);
    disp(X)
%}

%{
    %% Comparison density and relatived frequencies based on random numbers.
    close all
    % M is total number of balls.
    M=100; 
    n=10;  % number of balls taken
    odds=2; % Prob. of red balls vs other color balls
    K=50; % Number of red balls in the urn
    samplesize=100000; % Sample size to extract
    edges=(0:(n+1))';
    % Compute the density
    Fpdf=FNChygepdf(edges(1:end-1),M,K,n,odds);
    % Generate random numbers
    x=FNChygernd(M,K,n,odds,samplesize);
    % bar plot of theoretical and relative frequencies
    freqANDdens=[(histcounts(x,edges-0.5)/samplesize)' Fpdf];
    bar(edges(1:end-1),freqANDdens)
    legend(["Theoretical density" "Empirical relative frequency"],'Location','best')
    xlabel('Number of successes')
%}


%% Beginning of code

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
    r(i)=aux.FisherNCHypergeometricrnd(n,K,M,odds, accuracy);
end
end

%FScategory:ProbDist



