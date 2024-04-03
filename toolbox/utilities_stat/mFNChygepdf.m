function Wpdf = mFNChygepdf(x, m, w, accuracy)
%mFNChygepdf returns Fisher multivariate non-central hypergeometric probability density function
%
%<a href="matlab: docsearchFS('mFNChygepdf')">Link to the help function</a>
%
% This function calls function  CMultiFisherNCHypergeometricpdf which is a
% translation into MATLAB of the corresponding C++ function of Fog (2008).
% The notation which is used in mFNChygepdf and the order of the arguments
% is the one of MATLAB hyge. The notation which is used inside
% CMultiFisherNCHypergeometricpdf is the original one of Fog.
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
% In the Wallenius distribution the
% probability that the second ball picked is red depends on whether the
% first ball was red or white. If the first ball was red then the above
% formula is used with $K$ reduced by one. If the first ball was not red
% then the above formula is used with $M-K$ reduced by one. The number of red
% balls that we get in this experiment is a random variable with Wallenius'
% noncentral hypergeometric distribution.
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
% More formally, if the total number $n$ of balls taken is known then the
% conditional distribution of the number of taken red balls for given $n$ is
% Fisher's noncentral hypergeometric distribution. 
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
%           x    : Number of balls of each color sampled. Evaluation point.
%                  Vector or matrix.
%                  Evaluation points specified as 1-by-length(m) numeric
%                  vector or an ntrials-by-length(m) numeric matrix where
%                  ntrials is a positive scalar integer defining the
%                  number of evaluation points.
%                  Data Types - single|double
%           m    : Initial number of balls of each color in the urn.
%                  Vector. The length of vector M is the number of colors
%                  in the urn. M=sum(m) is the the total number of balls.
%                  If the number of colors is equal to 2, m(1)+m(2)
%                  corresponds to M and m(1) corresponds to K (number of
%                  red balls in the urn).
%                  Data Types - single|double
%           w    : Weights for each color. Vector. The length of vector w
%                  must be equal to the length of vector m. If the number
%                  of colors is equal to 2, odds=w(1)/w(2). If all the
%                  elements of w are equal we have the central
%                  hypergometric multivariate distribution.
%                  Data Types - single|double
%
%  Optional input arguments:
%
%       accuracy : accuracy of the calculations. Scalar. The default value
%                  of accuracy is 1e-08.
%                  Data Types - single|double
%                  Example - 1e-06
%
%
%  Output:
%
%           Wpdf : Fisher multivariate pdf values. 
%                  Scalar or column vector of length(ntrials).
%                  Wpdf(1) is the probability of drawing exactly x(1,1)
%                  balls of a possible m(1) items, x(1,2) balls of a
%                  possible m(2) items,  .., x(1,end) balls of a possible
%                  m(end) items in n drawings (where n=sum(1,:)) without
%                  replacement from a group of sum(m) objects, when objects
%                  are from length(m) weighted groups. The size of Wpdf is
%                  a scalar if x is a row vector. The size of Wpdf is a
%                  column vector of length ntrials if x is a matrix of size
%                  ntrials-by-length(m).
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
% Copyright 2008-2024.
%
%<a href="matlab: docsearchFS('mFNChygepdf')">Link to the help function</a>
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
    % The weight of each color is 1, 2,5 and 1.8
    w=[1, 2.5, 1.8];
    % We want to know the probability of extracting without replacement exactly 8 red balls, 10
    % white balls and 6 black balls (the total number of balls taken from
    % the urn is 8+10+6=20);
    x=[8, 10, 6];
    % Function mWNChygepdf is called without specifying fourth optional argument and
    % therefore the precision is set to 1e-08. 
    wpdf = mFNChygepdf(x,m,w);
    disp('We have an urn with 3 groups of balls; 20 red, 30 white and 20 black.');
    disp('In total the urn contains 70 balls.');
    disp('But the weight of getting the 3 balls is 1, 2.5 and 1.8;');
    disp('The probability of getting 8 red balls, 10, white and 6 black balls');
    disp('in an extraction without replacement is ');
    disp(wpdf);

%}

%{
    % In presence of just two colors in the urn it is possible to use  FNChygepdf or mFNChygepdf.
    % We have M balls in the urn.
    M  = 200;
    % We have K red balls
    K  = 6;
    %we extract n balls, one at a time, without replacement
    n  = 5;
    % red balls have a probability 100 times greater to be extracted than white balls
    odds  = 100;
    % We compute the probability of getting 0, 1, 2 or 3, ..., n red balls 
    x = (0:n)';
    % Densities using call to WNChygepdf
    wpdf = FNChygepdf(x,M,K,n,odds);
    % Densities using call to mWNChygepdf
    XX=[x n-x];
    m=[K M-K]; w=[100 1];
    wpdfCHK=mFNChygepdf(XX,m,w);
    assert(isequal(wpdf,wpdfCHK),'Call to FNChygepdf produces different results from call to mFNChygepdf')
%}

%{
    %% Check that the sum of densities is 1.
    weights=[5,2,1]'; % define the weights for each color
    m=[12 25 18];   % Define urn composition
    numberBallsExtracted=3;     % Number of balls which are drawn.
    
    % Define all possible cases in which for numberBallsExtracted;
    numcolors = length(m);      % Length of each permutation
    % Create all possible permutations (with repetition) of numcolors
    C = cell(numcolors, 1);         % Preallocate a cell array
    x=0:numberBallsExtracted;
    [C{:}] = ndgrid(x);     % Create grids of values
    Y = cellfun(@(x){x(:)}, C); % Convert grids to column vectors
    Y = [Y{:}];   % Create matrix with all possible permutations with repetion
    % Extract from Y the rows whose sum is equal to numberBallsExtracted
    boo=sum(Y,2)==numberBallsExtracted;
    Ysel=Y(boo,:);
    % Find the probability of each row of matrix Ysel
    Wpdf=mFNChygepdf(Ysel,m,weights);
    % Make sure that the sum of all probabilities is to (up to a certain
    % tolerance)
    tol=1e-08;
    assert(abs(sum(Wpdf)-1)<tol,"FSDA:Sum of densities is not equal to 1")
    lab=cellstr(num2str(Ysel));
    labx=categorical(lab,lab);
    bar(labx,Wpdf)
    title(['weights=[' num2str(weights') '] m=[' num2str(m) '] n=' num2str(numberBallsExtracted)])
    ylabel('Multivariate Fisher non central density')
%}

%{
    % Comparison between Wallenius and Fisher.
    weights=[100,2,1]'; % define the weights for each color
    m=[12 25 18];   % Define urn composition
    numberBallsExtracted=5;     % Number of balls which are drawn.
    
    % Define all possible cases in which for numberBallsExtracted;
    numcolors = length(m);      % Length of each permutation
    % Create all possible permutations (with repetition) of numcolors
    C = cell(numcolors, 1);         % Preallocate a cell array
    x=0:numberBallsExtracted;
    [C{:}] = ndgrid(x);     % Create grids of values
    Y = cellfun(@(x){x(:)}, C); % Convert grids to column vectors
    Y = [Y{:}];   % Create matrix with all possible permutations with repetion
    % Extract from Y the rows whose sum is equal to numberBallsExtracted
    boo=sum(Y,2)==numberBallsExtracted;
    Ysel=Y(boo,:);
    % Find the probability of each row of matrix Ysel
    WpdfFisher=mFNChygepdf(Ysel,m,weights);
    WpdfWallenius=mWNChygepdf(Ysel,m,weights);
    
    % Make sure that the sum of all probabilities is to (up to a certain
    % tolerance)
    lab=cellstr(num2str(Ysel));
    labx=categorical(lab,lab);
    bar(labx,[WpdfFisher WpdfWallenius])
    title(['weights=[' num2str(weights') '] m=[' num2str(m) '] n=' num2str(numberBallsExtracted)])
    ylabel('Multivariate Fisher and Wallenius non central hypergeometric density')
    legend(["Fisher" "Wallenius"])
%}

%% Beginning of code

if nargin<4
    accuracy =  1e-08;
end

if nargin < 3
    error(message('FSDA:mWNChygepdf:TooFewInputs'));
end

Wpdf = zeros(size(x,1),1);

for i=1:size(x,1)
    Wpdf(i) = aux.CMultiFisherNCHypergeometricpdf(x(i,:),m,w,accuracy);
end

end

%FScategory:ProbDist