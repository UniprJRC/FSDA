function odds = mWNChygeodds(mu, m, n)
%mWNChygeodds returns Wallenius' multivariate non-central hypergeometric odds given means
%
%<a href="matlab: docsearchFS('mWNChygeodds')">Link to the help function</a>
%
% This function calls function  oddsMWNCHypergeo_from_mean which is a
% translation into MATLAB of the corresponding C++ function of Fog (2008).
% The notation which is used in mWNChygeodds and the order of the arguments
% is the one of MATLAB hyge. The notation which is used inside
% oddsMWNCHypergeo_from_mean is the original one of Fog.
%
% To illustrate the meaning of Wallenius and Fisher' function parameters, let's use
% the classical biased urn example, with $K$ red balls and  $M-K$ white balls,
% totalling $M$ balls. $n$ balls are drawn at random from the urn
% without replacement. Each red ball has the weight $\omega_{1}$, and
% each white ball has the weight $\omega_{2}$; the probability ratio of red over
% white balls is then given by $odds = \omega_{1} / \omega_{2}$. Note that
% the odds are fixed once and for all during the drawings.
% 
% If the balls are taken one by one, the probability (say $p_1$) that the
% first ball picked is red is equal to the weight fraction of red balls:
% \[
% p_1= \frac{K w_1}{K w_1 + (M-K) w_2} 
% \] 
% In the Wallenius distribution the probability that the second ball picked
% is red depends on whether the first ball was red or white. If the first
% ball was red then the above formula is used with $K$ reduced by one. If
% the first ball was not red then the above formula is used with $M-K$
% reduced by one. The number of red balls that we get in this experiment is
% a random variable with Wallenius' noncentral hypergeometric distribution.
%
% The important fact that distinguishes Wallenius' distribution is that
% there is competition between the balls. The probability that a particular
% ball is taken in a particular draw depends not only on its own weight,
% but also on the total weight of the competing balls that remain in the
% urn at that moment. And the weight of the competing balls depends on the
% outcomes of all preceding draws. In the Fisher model, the fates of the
% balls are independent and there is no dependence between draws. One may
% as well take all $n$ balls at the same time. Each ball has no "knowledge"
% of what happens to the other balls.
% More formally, if the total number $n$ of balls taken is not known before
% the experiment (i.e n is determined just after the experiment), then the
% conditional distribution of the number of taken red balls for given $n$
% is Fisher's noncentral hypergeometric distribution.
%
% These two distributions have important applications in evolutionary
% biology and population genetics. If animals of a particular species are
% competing for a limited food resource so that individual animals are
% dying one by one until there is enough food for the remaining animals,
% and if there are different variants of animals with different chances of
% finding food, then we can expect the animals that die to follow a
% Wallenius noncentral hypergeometric distribution. Fisher’s noncentral
% hypergeometric distribution may be used instead of Wallenius distribution
% in cases where the fates of the individual animals are independent and
% the total number of survivors is known or controlled as part of
% a simulation experiment (Fog 2008).
% 
% Wallenius  distribution is a general model of biased sampling. Fisher’s
% noncentral hypergeometric distribution is used for statistical tests on
% contingency tables where the marginals are fixed (McCullagh and Nelder,
% 1983).
% 
% The multivariate Fishers and Wallenius noncentral hypergeometric
% distribution, are referred to the case where the type of different balls
% is greater than 2 and each ball has a different probability. 
%
% The two distributions are both equal to the (multivariate) hypergeometric
% distribution when the odds ratio is 1.
%
%  Required input arguments:
%
%           mu   : Measured mean for each color. Vector. Length of vector mu is equal
%                  to the number of colors. 
%                  The sum of the elements of mu must be equal to n.
%                  Data Types - single|double
%           m    : Initial number of balls of each color in the urn.
%                  Vector. The length of vector M is the number of colors
%                  in the urn. M=sum(m) is the the total number of balls.
%                  If the number of colors is equal to 2, m(1)+m(2)
%                  corresponds to M and m(1) corresponds to K (number of
%                  red balls in the urn).
%                  Data Types - single|double
%           n    : Total number of balls sampled. Scalar.
%                  Scalar which defines the number of balls which are
%                  drawn from the urn.
%                  Note that sum(mu) must be equal n, however, following
%                  the original C++ implementation of Fog we allow
%                  abs(sum(mu)-n)/n<=0.1
%                  Data Types - single|double
%
%  Optional input arguments:
%
%
%  Output:
%
%          odds :  Wallenius' odds (weights) for each color. Vector.
%                  Vector containing the estimated odds for each color
%                  given the input means. The reference color has odds=1.
%
%
% See also: WNChygepdf, WNChygecdf, WNChygeinv, WNChygernd, FNChygepdf, FNChygecdf, FNChygeinv, FNChygernd
%
% References:
%
% Fog, A. (2008), Calculation Methods for Wallenius' Noncentral Hypergeometric
% Distribution, "Communications in Statistics - Simulation and
% Computation", Vol. 37, pp. 258-273.
%
% Copyright 2008-2025.
%
%<a href="matlab: docsearchFS('mWNChygeodds')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% Problem description.
    % We have 70 balls in a urn of 3 different colors.
    % Initially in the urn we have 20 red balls, 30 white and 20 black balls.
    m=[20,30,20];
    % The vector of mean values is 
    mu=[15, 30, 5];
    % The extracted sample is n=50
    n=50;
    disp('We have an urn with 3 groups of balls; 20 red, 30 white and 20 black.');
    disp('In total the urn contains 70 balls.');
    disp('What are the weight of getting the 3 balls (odds)?');
    disp('given the above values of mu')
    odds=mWNChygeodds(mu,m,n);
    disp(odds);
%}

%{
    % Find the odds given mu, m and n.
    % Vector of mu values
    mu=[8 46 29 17]; 
    % Number of balls in the urn for each color
    m=[30,60,50,20];
    n=sum(mu); % number of extracted balls
    odds=mWNChygeodds(mu, m, n);
    disp('Odds')
    disp(odds)
%}

%{
    % An example in which one of the odds is Inf. 
    % The vector of mean values is 
    mu=[15, 30, 5]; 
    % The number of balls inside the urn for each color 
    m=[20,30,20]; 
    % The extracted sample is n=50 
    n=sum(mu); 
    disp('We have an urn with 3 groups of balls; 20 red, 30 white and 20 black.');
    disp('In total the urn contains 70 balls.'); 
    disp('What are the weight of getting the 3 balls (odds)?'); 
    disp('given the above values of mu')
    odds=mWNChygeodds(mu,m,n); 
    disp(odds);
%}

%{
    %% An example in which  one of the odds is 0
    m=[20,30,20];
    % The vector of mean values is 
    mu=[6, 0, 5];
    % The extracted sample is n=10
    n=10;
    disp('We have an urn with 3 groups of balls; 20 red, 30 white and 20 black.');
    disp('In total the urn contains 70 balls.');
    disp('What are the weight of getting the 3 balls (odds)?');
    disp('given the above values of mu')
    odds=mWNChygeodds(mu,m,n);
    disp(odds);
%}

%{
    % An example where one odd is much bigger.
    m=[200,300,200,10];
    % The vector of mean values is 
    mu=[6, 1, 5, 9];
    % The extracted sample n
    n=sum(mu);
    odds=mWNChygeodds(mu,m,n);
    disp(odds);
%}

%{
    %% Comparison between Fisher and Wallenius odds.
    % An example where one odd is much greater than the others
    m=[200,300,200,10];
    % The vector of mean values is 
    mu=[6, 1, 5, 9];
    % The extracted sample n
    n=sum(mu);
    oddsF=mFNChygeodds(mu,m,n);
    oddsW=mWNChygeodds(mu,m,n);
    % Row names and column names
    rn="odds"+(1:length(mu));
    cn=["Fisher" "Wallenius"];
    oddsT=array2table([oddsF, oddsW],"RowNames",rn,"VariableNames",cn);
    disp(oddsT);
%}

%% Beginning of code

if nargin < 3
    error(message('FSDA:mWNChygeodds:TooFewInputs'));
end
    odds = aux.oddsMWNCHypergeo_from_mean(mu,m,n);
end

%FScategory:ProbDist
