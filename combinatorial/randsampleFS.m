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
%   If method=1 (default option) the program proceeds as follows:
%   if  $4*k >n$ the programs does a random permutation of the population
%   and returns the first nsel elements else if $4*k<=n$ (that is if the
%   desired sample is small compared to all combinations, the program
%   repeatedly samples with replacement until there are nsel unique values.
%
%   If method=2 Systematic sampling is used where the starting point is
%   random and the step is also random.
%
%   If method=3 random sampling based on the old but well known Linear
%   Congruential Generator (LCG) method is used. In this case there is no
%   guarantee to get unique numbers.
%
%   If method is a vector of n weights, then Weighted Sampling Without
%   Replacement is applied. Our implementation follows Efraimidis and Spirakis
%   (2006). MATLAB function datasample follows Wong and  Easton (1980), which is
%   also quite fast; note however that function datasample may be very slow if
%   applied repetedly, for the large amount of time spent on options checking.
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
%   For Weighted Sampling Without Replacement: 
%   Efraimidis, P.S. and Spirakis, P.G. (2006). Weighted random sampling with a reservoir. 
%   Information Processing Letters, 97, 181-185, 2006;
%   Wong, C. K. and M. C. Easton. An Efficient Method for Weighted Sampling Without Replacement.
%   SIAM Journal of Computing 9(1), pp. 111-113, 1980.
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
    % Extract k=10 number in [-1000 -900] with gamma distributed weights.
     population = -1000:1:-900; 
     n = numel(population);
     wgts = sort(random('gamma',0.3,2,n,1),'descend'); 

     k=10; 
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
    wgts = method(:)';
    method = 0;
end

switch method
    
    case 0
        
        % Weighted Sampling Without Replacement
        
        if k>n, error('k must be smaller than n'), end
        if length(wgts)~=n,error('the length of the weight vector must be n'),end

        % To achieve the best trade-off between simplicity and efficiency, following Efraimidis and
        % Spirakis (2006) we generate a weighted random sample without replacement of size $k<n$ by:
        % 1. Drawing uniformly in [0,1] $n$ independent values $u_i$;
        % 2. Computing the keys $U_i = u_i^{1/wgts_i}$;
        % 3. Keeping the $k$ elements with largest $U_i$. 
        % REMARK 1: This approach does not require the weights normalization wgts = wgts / sum(wgts)
        % REMARK 2: We take the log to avoid numerical issues when (positive) weights are << 1
        u = rand(1,n);
        U = (1./wgts) .* log(u); % U = u.^(1./wgts) may give rise to numerical issues 
        [~,y] = sort(U,'descend');
        y     = y(1:k); 

        % The next implementation is very intuitive, but also very inefficient
        %{
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
        %}

        % This implementation is equivalent to the previous one, but a bit more efficient  
        %{
            x = 1:n;
            p = wgts(:)' / sum(wgts);
            y2 = zeros(1,k);
            for i=1:k
                yi = 1 + sum( rand() > cumsum(p) );
                y2(i) = yi;
                p(x==yi)=0;
            end
        %}        
        
        % This is a similar implementation for Weighted Sampling *With* Replacement:
        %{
            [~, y] = histc(rand(k,1),cumsum([0;p(:)./sum(p)]));
        %}

        
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
        
        % Systematic sampling method, Cochran (1977), third edition, Sampling Techniques, Wiley.
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


