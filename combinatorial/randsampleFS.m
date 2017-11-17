function y = randsampleFS(n,k,method,after2011b)
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
%            Methods used to extract the subsets. See more about for details.
%            Default is method = 0.
%            - Scalar from 0 to 3 determining the method used to extract
%             (without replacement) the random sample.
%            - Vector of weights: in such a case, a weighted sampling without
%              replacement algorithm is applied using that vector of weights.
%            Example - randsampleFS(100,10,2)
%            Data Types - single|double
%
%after2011b: MATLAB version flag. Logical. Indicates if the MATLAB version 
%            in use is later than R2012a (7.14). Used to speed up
%            computations in function subsets and, more in general, in
%            simulation experiments which use randsampleFS intensively.
%            Example - randsampleFS(100,10,2,true) or
%            randsampleFS(100,10,2,~verLessThan('MATLAB','7.14'))
%            Data Types - logical
%
%   Output:
%
%   y :     A column vector of k values sampled at random from the integers 1:n.
%           For methods 0, 1, 2 and weighted sampling the elements  extracted
%           are unique; For method 3 (included for historical reasons) there is
%           no guarantee that the elements extracted are unique.
%           Data Types - single|double.
%
% More About:
%
%   The method=0 uses MATLAB function randperm. In old MATLAB releases
%   randperm was slower than FSDA function shuffling, which is used in
%   method 1 (for example, in R2009a - MATLAB 7.8 - randperm was at least
%   50% slower).  
%
%   If method=1 the approach depends on the population and sample sizes:
%   - if $n < 1000$ and $k < n/(10 + 0.007n)$, that is if the population is
%     relatively small and the desired sample is small compared to the
%     population, we repeatedly sample with replacement until there are k
%     unique values;
%   - otherwise, we do a random permutation of the population and return  
%     the first k elements. 
%   The threshold $k < n/(10 + 0.007n)$ has been determined by simulation
%   under MATLAB R2016b. Before, the threshold was $n < 4*k$.
%
%   If method=2 systematic sampling is used, where the starting point is
%   random and the step is also random. 
%
%   If method=3 random sampling is based on the old but well known Linear
%   Congruential Generator (LCG) method. In this case there is no guarantee
%   to get unique numbers. The method is included for historical reasons.
%
%   If method is a vector of n weights, then Weighted Sampling Without
%   Replacement is applied. Our implementation follows Efraimidis and
%   Spirakis (2006). MATLAB function datasample follows Wong and  Easton
%   (1980), which is also quite fast; note however that function datasample
%   may be very slow if applied repetedly, for the large amount of time
%   spent on options checking.
%
%   Remark on computation performances. Method=2 (systematic sampling) is
%   by far the fastest for any practical population size $n$. For example,
%   for $n \approx 10^6$ method=2 is two orders of magniture faster than
%   method=1. With recent MATLAB releases (after R2011b) method = 0 (which
%   uses compiled MATLAB function randperm) has comparable performances, at
%   least for reasonably small $k$. In releases before 2012a, randperm was
%   considerably slow.
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
%   For Method 3. D. E. Knuth. (1997). The Art of Computer Programming, Volume 2: Seminumerical
%   Algorithms, Third Edition. Addison-Wesley, Section 3.2.1: The
%   Linear Congruential Method, pp. 10-26.
%
%   For Weighted Sampling Without Replacement: Efraimidis, P.S. and Spirakis, P.G. (2006). 
%   Weighted random sampling with a reservoir.
%   Information Processing Letters, 97, 181-185.
%   Wong, C. K. and M. C. Easton, (1980). An Efficient Method for Weighted Sampling Without Replacement.
%   SIAM Journal of Computing 9(1), pp. 111-113.
%
% Copyright 2008-2017.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('randsampleFS')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %% randsampleFS with default options.
    % default method (1) is used.
    randsampleFS(1000,10)
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

% randsampleFS needs to check the MATLAB version in use in order to:
% - decide the sempling method to use, and
% - use properly the optional parameter of randperm.
% In the first case  the release to check is R2012a, i.e. 7.14
% In the second case the release to check is R2011b, i.e. 7.13
% For the sake of computational efficiency, we just check the latest
% To pass the argument, use:
% after2012a = ~verLessThan('MATLAB','7.14');
if nargin < 4 
    after2011b = ~verLessThan('MATLAB','7.14');
end

% choose the default sampling method
if nargin < 3 || isempty(method)
    if after2011b
        method=0;
    else 
        % in older releases we use systematic sampling because randperm
        % function (used for method = 0) was extremely inefficient
        method=2;
    end
end

% Weighted Sampling Without Replacement
% This is done if the third argument is provided as a vector of weights.
if nargin >= 3 && ~isscalar(method)
    wgts = method(:)';
    method = 999;
end

switch method
    
    case 0
        
        % Extract a random sample of k integers between 1 and n.
        if after2011b % it is a logical: no need to specify before714 == true
            y = randperm(n,k);
        else
            % In Matlab 2011b and older, input parameter k was not supported.
            % This solution is clearly less efficient for k much smaller than n.
            rp = randperm(n);
            y = rp(1:k);
        end
        
    case 1
        
        if  n < 1000 && k < n/(10 + 0.007*n)
            
            % If the population size n is snall, say smaller than 1000, and the
            % desired sample size k is small compared to n, it may be convenient
            % to repeatedly sample with replacement until there are k unique
            % values.
            % REMARK: till December 2016 the threshold used to choose between
            % this or the shuffling approach was n < 4*k. The new threshold has
            % been determined by simulation (with MATLAB R2016b) for n in the
            % interval [100 1000].

            mindiff = 0;
            while mindiff == 0
                y = randi(n, 1 , k);
                mindiff = min(diff(sort(y)));
            end

        else
            
            % This is the FSDA alternative to MATLAB function randsample, to
            % randomize the whole population and take the first k elements. Note
            % that function shuffling (see the 'utilities' folder) randomises
            % the combinations without calling function sort.
            
            rp = shuffling(1:n);
            y = rp(1:k);
        end
        
    case 2
        
        % Systematic sampling method, Cochran (1977), third edition,
        % Sampling Techniques, Wiley.
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
                % Triple of Goodman - Miller
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
        
        
    case 999
        
        % Weighted Sampling Without Replacement
        
        if k>n,  error('FSDA:randsampleFS:WrongInpArgs','k must be smaller than n'), end
        if length(wgts)~=n,error('FSDA:randsampleFS:WrongInpArgs','The length of the weight vector must be n'),end
        
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
        
end

end

%FScategory:UTICOMB


