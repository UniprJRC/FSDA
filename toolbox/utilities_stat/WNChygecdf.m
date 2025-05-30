function [Wcdf] = WNChygecdf(x,M,K,n,odds, accuracy)
%WNChygecdf returns Wallenius' non-central hypergeometric cumulative distribution function
%
%<a href="matlab: docsearchFS('WNChygecdf')">Link to the help function</a>
%
% This function calls function  WalleniusNCHypergeometricpdf.m which is a
% translation into MATLAB of the corresponding C++ function of Fog (2008).
% The notation which is used in WNChygecdf and the order of the arguments
% is the one of MATLAB hyge. The notation which is used inside
% WalleniusNCHypergeometriccdf is the original one of Fog.
%
%  Required input arguments:
%
%           x    : Number of red balls sampled. Scalar, vector or matrix.
%                  Data Types - single|double
%           M    : Total number of balls in urn before sampling. Scalar, vector or matrix.
%                  Data Types - single|double
%           K    : Initial number of red balls in the urn. Scalar, vector or matrix.
%                  Data Types - single|double
%           n    : Total number of balls sampled. Scalar, vector or matrix.
%                  Data Types - single|double
%        odds    : Probability ratio of red over white balls. Scalar, vector or matrix.
%                  Data Types - single|double
%
%
%  Optional input arguments:
%
%       accuracy : accuracy of the calculations. Scalar. The default value
%                  of accuracy is 1e-10.
%                  Example - 1e-06
%                  Data Types - single|double
%
%  Output:
%
%           Wcdf : Wallenius' cdf values. Cumulative probability of drawing
%                  exactly x or less than x of a possible K items in n
%                  drawings without replacement from a group of M objects,
%                  when objects are from two weighted groups. The size of
%                  Wcdf is the common size of the input arguments. A scalar
%                  input functions as a constant matrix of the same size as
%                  the other inputs.
%
%
% See also WNChygepdf, WNChygeinv, WNChygernd, randsampleFS.m, subsets.m
%
% References:
%
% Fog, A. (2008), Calculation Methods for Wallenius' Noncentral Hypergeometric
% Distribution, "Communications in Statistics - Simulation and
% Computation", Vol. 37, pp. 258-273.
%
% Copyright 2008-2025.
%
%<a href="matlab: docsearchFS('WNChygecdf')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Cumulative probability of getting 0 to x successes in n weighted drawns without replacement.
    % Problem description.
    % We have 500 balls in the urn
    M  = 500;   
    %we extract 3 balls, one at a time, without replacement
    n  = 3;     
    %initially, in the urn we have 250 red and 250 white balls
    K  = M/2;   
    % red balls are ten times more likely to be extracted than the white balls
    odds  = 10;    
    % We compute the cumulative probability of getting 0, 1, 2 red balls 
    % (in drawing the 2 balls without replacement).
    x = 2;
    wcdf = WNChygecdf(x,M,K,n,odds);
    disp('See WNChygecdf;');
    disp(wcdf);
%}

%% Beginning of code

if nargin <6
    accuracy=1e-10;
end

Wcdf = zeros(size(x));
for j=1:numel(x)
    for jn = 0:x(j)
        Wcdf(j) = Wcdf(j) + WNChygepdf(jn,M(j),K(j),n(j),odds(j),accuracy);
    end
end

end

%FScategory:ProbDist

