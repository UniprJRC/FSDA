function y = randsampleFS(n,k,method)
%randsampleFS generates a random sample of k elements from the integers 1 to n (k<=n)
%
%<a href="matlab: docsearch('randsampleFS')">Link to the help function</a>
%
%  Required input arguments:
%
%       n : scalar, a positive integer. 
%           a vector of numbers will be selected from the integers 1 to n.
%      k  : the number of elements to be selected.
%
%  Optional input arguments:
%
%   method = scalar which determines the method used to extract the numbers. Default
%           is method = 1.
% 
%   if method=1 (default option) the program proceeds as follows:
%           if 4*k > n the programs does a random permutation of the
%              population and returns the first nsel elements
%           else if 4*k<=n (that is if the desired sample is small
%              compared to all combinations, the program repeatedly samples
%              with replacement until there are nsel unique values.
%   if method=2 Systematic sampling is used where the starting point
%           is random and the step is also random
%   if method=3 random sampling based on the old but well known Linear
%           Congruential Generator (LCG) method is used. In this case there
%           is no guarantee to get unique numbers. 
%
%   Output: 
%
%   y = randsampleFS(n,k) returns y as a column vector of k
%   values sampled at random from the integers 1:n. 
%   For method 1 and 2, the elements  extracted are unique. 
%   For method 3, there is no guarantee that the elements extracted are
%   unique.
%   
% See also: randsample (of statistics toolbox) and shuffling (FSDA)
%
% References (one for each method, in order):
%
%   Fisher, R.A.; Yates, F. (1948) [1938]. Statistical tables for biological,
%   agricultural and medical research (3rd ed.). London: Oliver & Boyd. pp.
%   26–27.
%
%   Cochran, William G. (1977). Sampling techniques (Third ed.). Wiley.
%
%   D. E. Knuth. The Art of Computer Programming, Volume 2: Seminumerical
%   Algorithms, Third Edition. Addison-Wesley, 1997. Section 3.2.1: The
%   Linear Congruential Method, pp. 10–26.
%
%
%
% Copyright 2008-2011.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti 
%            and Vytis Kopustinskas (2009-2010)
%
%
%<a href="matlab: docsearch('randsampleFS')">Link to the help function</a>
% Last modified 15-Nov-2011

% Examples:

%{

  randsampleFS(100,10) % default method (1) is used

  randsampleFS(100,10,2) % default method (1) is used

  randsampleFS(100,10,3) % default method (1) is used. Here, being nsel so
                         % biG wrt nsamp, it is likely to obtain repetitions

%}
%

%% Beginning of code
if nargin<3;
    method=1;
end

switch method
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
                % Triple of Leormonth – Lewis
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
                % Triple of Goodman – Miller
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




