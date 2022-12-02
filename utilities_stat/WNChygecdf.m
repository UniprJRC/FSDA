function [Wcdf] = WNChygecdf(X,N,K,M,W)
%WNChygecdf returns Wallenius' non-central hypergeometric Cumulative Distribution Function values.
%
%<a href="matlab: docsearchFS('WNChygecdf')">Link to the help function</a>
%
% This function is taken from the toolbox "Generation of Random Variates"
% (function wallen_cdf.m) created by Jim Huntley. See WNChygepdf for details.
%
%  Required input arguments:
%           X    : Number of red balls sampled. Scalar.
%                  Data Types - single|double
%           N    : Total number of balls sampled. Scalar.
%                  Data Types - single|double
%           K    : Initial number of red balls in the urn. Scalar.
%                  Data Types - single|double
%           M    : Total number of balls in urn before sampling. Scalar.
%                  Data Types - single|double
%           W    : Probability ratio of red over white balls. Scalar.
%                  Data Types - single|double
%
%  Optional input arguments:
%
%  Output:
%
%           Wcdf : Wallenius' cdf values. Cumulative probability of drawing exactly X of a
%                  possible K items in N drawings without replacement from a
%                  group of M objects, when objects are from two weighted groups.
%                  Array of numel(x) values.
%                  Data Types - single|double.
%
%
% See also WNChygepdf, randsampleFS.m, subsets.m
%
% References:
% Fog, A. (2008), Calculation Methods for Wallenius' Noncentral Hypergeometric
% Distribution, "Communications in Statistics - Simulation and
% Computation", Vol. 37, pp. 258-273.
%
% Copyright 2008-2023.
% FSDA team adapted and documented Jim Huntley's function wallen_pdf for
% illustration purposes only.
%
%<a href="matlab: docsearchFS('WNChygecdf')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Cumulative probability of getting 0 to p successes in p weighted drawns without replacement.
    % Problem description.
    % we have 500 balls in the urn
    M  = 500;   
    %we extract 3 balls, one at a time, without replacement
    N  = 3;     
    %initially, in the urn we have 250 red and 250 white balls
    K  = M/2;   
    %red balls are ten times the white balls
    W  = 10;    
    % We compute the cumulative probability of getting 3 red balls in drawing
    % 2 balls without replacement.
    x = 2;
    wcdf = WNChygecdf(x,N,K,M,W);
    disp('See WNChygecdf;');
    disp(wcdf);
%}

%% Beginning of code

Wcdf = WNChygepdf(0,N,K,M,W);
for jn = 1:X
    Wcdf = Wcdf + WNChygepdf(jn,N,K,M,W);
end

end

%FScategory:UTISTAT

