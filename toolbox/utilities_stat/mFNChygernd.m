function r = mFNChygernd(m, w, n,    ntrials, accuracy)
%mFNChygernd returns Fisher multivariate non-central hypergeometric random variate generation
%
%<a href="matlab: docsearchFS('mFNChygernd')">Link to the help function</a>
%
% This function calls function  CMultiFisherNCHypergeometricrnd which is a
% translation into MATLAB of the corresponding C++ function of Fog (2008).
% The notation which is used in mFNChygernd and the order of the arguments
% is the one of MATLAB hyge. The notation which is used inside
% CMultiFisherNCHypergeometricrnd is the original one of Fog.
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
%                  hypergeometric multivariate distribution.
%                  Data Types - single|double
%           n    : Total number of balls sampled. Scalar.
%                  Scalar which defines the number of balls which are
%                  drawn from the urn.
%                  Data Types - single|double
%
% Optional input arguments:
%
%    ntrials    : Number of random variates which have to extracted. If
%                 this argument is not specified just one random variate is
%                 generated.
%               Example - 3
%               Data Types - double
%       accuracy : accuracy of the calculations. Scalar. The default value
%                  of accuracy is 1e-10.
%                  Example - 1e-06
%                  Data Types - single|double
%
% Output:
%
%         r    : Random numbers. Vector or matrix. Array of random numbers of size
%                ntrials-by length(w).
%
%
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
%<a href="matlab: docsearchFS('mFNChygernd')">Link to the help function</a>
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
    % n=number of balls which are taken
    n=10; 
    % Generate a random variate from this distribution
    r = mFNChygernd(m,w, n);
    disp(r);

%}

%{
    %% Generate 5 random variates from a multivariate non central Wallenius distribution.
    % m = number of balls of each type in the urn.
    m=[20,30,20, 10];
    % w = vector containing the weight of each color 
    w=[1, 2.5, 1.8, 100];
    % n=number of balls which are taken
    n=10; 
    % Generate 5 random variates from this distribution
    ntrials=5;
    R = mFNChygernd(m,w, n, ntrials);
    % R is matrix of size ntrials-by-length(m)
    disp(R);
%}

%{
    %% Check agreement between theoretical and empirical relative frequencies.
    % Define the weights for each color.
    weights=[10,2,1]; 
    % Define urn composition
    m=[12 25 18];   
    % Number of balls which are drawn.
    numberBallsExtracted=7;     
    % number of random variates which are generated 
    ntrials=100000; 
    % Define all possible cases in which for numberBallsExtracted;
    numcolors = length(m);     
    % Create all possible permutations (with repetition) of numcolors
    % Preallocate a cell array
    C = cell(numcolors, 1);         
    x=0:numberBallsExtracted;
    % Create grids of values
    [C{:}] = ndgrid(x);     
    % Convert grids to column vectors
    Y = cellfun(@(x){x(:)}, C); 
    % Create matrix with all possible permutations with repetion
    Y = [Y{:}];   
    % Extract from Y the rows whose sum is equal to numberBallsExtracted
    boo=sum(Y,2)==numberBallsExtracted;
    Ysel=Y(boo,:);
    % Find the probability of each row of matrix Ysel
    Wpdf=mFNChygepdf(Ysel,m,weights);
    % Make sure that the sum of all probabilities is to (up to a certain
    % tolerance)
    tol=1e-08;
    assert(abs(sum(Wpdf)-1)<tol,"FSDA:Sum of densities is not equal to 1")
    Outcomes=cellstr((num2str(Ysel)));
    % Generate a matrix of ntrials random variates
    Wrnd=mFNChygernd(m,weights,numberBallsExtracted,ntrials);
    % Compute the frequency distribution (pivot table) of observed outcomes
    Rd=cellstr(num2str(Wrnd));
    Rdtable=table(Rd);
    OutcomesObs=pivot(Rdtable,'Rows','Rd');
    
    % Make sure that there is a matching between the rows of theoretical and
    % empirical frequencies
    [int,ia,ib]=intersect(Outcomes,OutcomesObs{:,1});
    % Define the table which will contain both theoretical (density values) and
    % empirical frequencies
    Freq=table('Size',[size(Ysel,1),2],'VariableTypes',{'double' 'double'});
    Freq.Properties.RowNames=Outcomes;
    Freq{:,1}=Wpdf;
    Freq{ia,2}=OutcomesObs{:,2}/ntrials;
    Freq.Properties.VariableNames=["Densities" "Relative frequencies"];
    % Create categorical object in order to label x axis
    OutcomesC=categorical(Outcomes,Outcomes);
    bar(OutcomesC,Freq{:,1:2})
    legend(["Theoretical probabilities" "Relative frequencies"])
%}

%% Beginning of code

if nargin < 3
    error(message('FSDA:mWNChygernd:TooFewInputs'));
end
if nargin<5
    accuracy=1e-08;
end

if nargin <4
    ntrials=1;
end


r=zeros(ntrials,length(m));
for i=1:ntrials
    r(i,:)=aux.CMultiFisherNCHypergeometricrnd(m,w, n, accuracy);
end

end

%FScategory:ProbDist