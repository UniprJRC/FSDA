function y = randsampleFS(n,k,method)
%randsampleFS generates a random sample of k elements from the integers 1 to n (k<=n)
%
%<a href="matlab: docsearchFS('randsampleFS')">Link to the help function</a>
%
%  Required input arguments:
%
%       n : A vector of numbers will be selected from the integers 1 to n.
%           Scalar, a positive integer.
%           Data Types - single|double
%      k  : The number of elements to be selected. Non negative integer.
%           Data Types - single|double
%
%  Optional input arguments:
%
%   method : Sampling methods. Scalar or vector.
%            Methods used to extract the subsets. See more about for
%            further details. Default is method = 1.
%            - Scalar, from 1 to 3 determining the (random sample without
%            replacement) method to be used.
%            - Vector of weights: in such a case, Weighted Sampling Without
%              Replacement is applied using that vector of weights.
%            Example - randsampleFS(100,10,2)
%            Data Types - single|double
%
%   Output:
%
%   y :     A column vector of k values sampled at random from the integers 1:n.
%           For method 1 and 2, the elements  extracted are unique; For
%           method 3, there is no guarantee that the elements extracted are
%           unique.
%           Data Types - single|double.
%
% More About:
%
%   if method=1 (default option) the program proceeds as follows:
%   if  $4*k >n$ the programs does a random permutation of the population
%   and returns the first nsel elements else if $4*k<=n$ (that is if the
%   desired sample is small compared to all combinations, the program
%   repeatedly samples with replacement until there are nsel unique values.
%
%   if method=2 Systematic sampling is used where the starting point is
%   random and the step is also random.
%
%   if method=3 random sampling based on the old but well known Linear
%   Congruential Generator (LCG) method is used. In this case there is no
%   guarantee to get unique numbers.
%
%   if method is a vector of n weights, then Weighted Sampling Without
%   Replacement is applied. Our implementation could be improved. 
%   The best algorithm for Weighted Sampling Without Replacement, mentioned
%   in the references, is applied by MATLAB function datasample, which is
%   unfortunately very slow for the large amount of time spent on options
%   checking.
%
%
% See also: randsample, datasample, shuffling
%
% References:
%
%   For Method 1. Fisher, R.A.; Yates, F. (1948) [1938]. Statistical tables
%   for biological, agricultural and medical research (3rd ed.). London,
%   Oliver & Boyd, pp. 26-27.
%
%   For Method 2. Cochran, William G. (1977). Sampling techniques (Third ed.). Wiley.
%
%   For Method 3. D. E. Knuth. The Art of Computer Programming, Volume 2: Seminumerical
%   Algorithms, Third Edition. Addison-Wesley, 1997. Section 3.2.1: The
%   Linear Congruential Method, pp. 10-26.
%
%   For Weighted Sampling Without Replacement. Wong, C. K. and M. C.
%   Easton. An Efficient Method for Weighted Sampling Without Replacement.
%   SIAM Journal of Computing 9(1), pp. 111?113, 1980.
%
% Copyright 2008-2016.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('randsampleFS')">Link to the help function</a>
%
% Last modified 31-05-2016
%
% Examples:

%{
    %% randsampleFS with default options.
    % default method (1) is used.
    randsampleFS(100,10)
%}

%{
    %% randsampleFS with optional argument set to method (2).
    method = 2;
    randsampleFS(100,10,method)
%}

%{
    % randsampleFS with optional arguments set to method (3).
    method = 3;
    % Here, being nsel so big wrt nsamp, it is likely to obtain repetitions.
    randsampleFS(100,10,method)
%}

%{
    % randsampleFS Weighted Sampling Without Replacement.
    % Extract 10 number from [-1000 -900] with normal distributed weights.
     population = -1000:1:-900; 
     k=10; 
     wgts=sort(random('gamma',0.3,2,101,1),'descend'); 
     n = numel(population);
     y = randsampleFS(n,k,wgts);
     sample  = population(y);
     
     plot(wgts,'.r')
     hold on;
     text(y,wgts(y),'X');
     title('Weight distribution with the extracted numbers superimposed')
%}

%% Beginning of code
if nargin<3
    method=1;
end

% Weighted Sampling Without Replacement
% This is done if the third argument is provided as a vector of weights.
if nargin == 3 && ~isscalar(method)
    wgts = method;
    method = 0;
end

switch method
    
    case 0
        
        % Weighted Sampling Without Replacement
        
        % The current implementation is intuitive but sub-optimal
        % It will be improved in next releases
          
        if k>n, error('k must be smaller than n'), end
        if length(wgts)~=n,error('the length of the weight vector must be n'),end
        
        x = 1:n;
        y = zeros(1,k);
        p = wgts(:)' / sum(wgts);
        for i=1:k
            %the following three lines are equivalent (but faster than) 
            %v(i)=randsample(x,1,true,wgts)
            % i.e. return in yi a weighted sample, of 1 element only,
            % taken *with replacement* from the set 1:n, using a vector of
            % positive weights (probabilities) p whose length is n.
            % The 'min' function is to avoid probabilities begger than 1 
            % due to accumulation of roundoff errors (this was actually
            % observed).
            edges  = min([0 cumsum(p)],1);
            [~ , yi] = histc(rand(1,1),edges);
            y(i) = yi;

            % Now there is one element less in the sample
            % (weighted sampling must be *without* replacement)
            p(x==yi)=0;
            % and the new probabilities must be re-normalised
            p = p / sum(p);
        end
        
    case 1
        
        if 4*k > n
            % If the sample is a reasonable fraction of all combinations,
            % just randomize the whole population and take the first nsel.
            % Note that function shuffling (see the 'utilities' folder)
            % randomises the combinations without calling function sort.
            rp = shuffling(1:n);
            y = rp(1:k);
        else
            % If the desired sample is small compared to all combinations,
            % it may be more convenient to repeatedly sample with
            % replacement until there are nsel unique values.
            mindiff = 0;
            while mindiff == 0
                % sample w/replacement
                y = randi(n, 1 , k);
                mindiff = min(diff(sort(y)));
            end
        end
        
    case 2
        % Systematic sampling method, Cochran (1977), third edition, Sampling
        % Techniques, Wiley.
        stepk=floor(n/k);
        startk=randi(n);
        y=startk:stepk:startk+stepk*k-1;
        
        logi=y>n;
        y(logi)=y(logi)-n;
        
    case 3
        
        % A Linear Congruential Generator (LCG) represents one of the
        % oldest and best-known pseudorandom number generator algorithms
        terna = 5;
        switch terna
            case 1
                % Triple of Leormonth ? Lewis
            case 2
                m = 2^31;
                a = 2^16;
                c = 0;
            case 3
                % Triple of Knuth
                m = 2^31;
                a = floor(pi * 10^8);
                c = 45380624;
            case 4
                % Triple of Goodman ? Miller
                m = 2^31-1;
                a = 7^5;
                c = 0;
            case 5
                % Triple of Lehmer
                m = 2^31 - 1;
                a = 16807 ;
                c = 0;
        end
        
        y = NaN(1,k);
        y(1,1) =  randi(n) ;
        for i = 1 : k-1
            y(1,i+1) = mod(a * y(i) + c , m);
        end
        y = ceil(y * n / m);
end

end

%FScategory:UTICOMB


