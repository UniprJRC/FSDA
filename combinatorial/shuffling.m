function x = shuffling(x)
%shuffling does a random permutation of the elements of input vector
%
%<a href="matlab: docsearchFS('shuffling')">Link to the help function</a>
%
%  Required input arguments:
%
%  x :  A set of elements. Vector of length t. 
%       Data Types - single|double
%
%  Optional input arguments:
%
%  Output: 
%
%  x :  A permutation of the set of elements in x. Vector of length t. 
%       Data Types - single|double
%
% More About:
%
% If set $x$ has $t$ elements, the objective is to obtain each of the $t!$
% pemutations with equal probability, especially when $t$ is large. To
% achieve this goal we use backward Knuth's shuffling, which is based on
% the Fisher-Yates shuffle.
%
% shuffling has been introduced as an alternative to MATLAB function
% randperm. Randperm makes a call to sort(rand(1,n)) and, overall, is
% slower than shuffling (for example, in R2009a shuffling was on average
% 25% faster). If compiled as mex file, shuffling becomes much more
% efficient than x(randperm(numel(x))) solution (it is about 60% faster for
% n=10^6). C code that can be used for this purpose is available at the
% http://it.mathworks.com/matlabcentral/fileexchange/27076-shuffle website
% as part of Jan Simon's Shuffle library.
%
% See also: randperm
%
% References:
%
% Knuth, Donald E. (1969). The Art of Computer Programming volume 2,
% Seminumerical algorithms, Reading, MA: Addison-Wesley, pp. 124-125.
%
% Fisher, R.A.; Yates, F. (1948) [1938]. Statistical tables for biological,
% agricultural and medical research (3rd ed.). London, Oliver & Boyd. pp
% 26-27.
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('shuffling')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:
%{
    %% shuffling applied to a set of 20 elements.
    shuffling(1:20)
%}

%{
    % shuffling applied with parsimonious data type.
    % shuffling applied to a set of 20 elements, but using a parsimonious
    % data type; this is convenient if the vector is big.
    shuffling(int8(1:20))
%}

%{
    % check of the permutation produced by shuffling.
    x = 1:200000;
    numel(unique(shuffling(x)))
%}

%{
    % Computation time of shuffling and randperm. 
    % An extensive test for various sample sizes.
    % REMARK: shuffling code is interpreted whereas randperm is compiled;
    % therefore, the comparison has to be done using tic-toc statements, 
    % as in the example below (the MATLAB profiler would over-estimate the 
    % shuffling time). 
    RefTime = 1;
    for n = [10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7]
       X = 1:n;

       % Estimate the number of profiling loops:
       iniTime = cputime;
       countLoop = 0;
       while cputime - iniTime < RefTime
          Xp = X(randperm(n));  
          clear('Xp');   % Necessary to suppress JIT acceleration
          countLoop = countLoop + 1;
       end
       nDigit = max(1, floor(log10(max(1, countLoop))) - 1);
       nLoop  = max(4, round(countLoop / 10 ^ nDigit) * 10 ^ nDigit);

       % monitor randperm
       tic;
       for i = 1:nLoop
          Xp = X(randperm(n));
          clear('Xp');
       end
       RandPermTime = toc;

       % monitor shuffling
       tic;
       for i = 1:nLoop
          Xp = shuffling(X);
          clear('Xp');
       end
       ShufflingTime = toc;

       % results
       fprintf('\n%d elements shuffled %d times: \n', n, nLoop);
       disp(['    - shuffling etime (seconds): ' num2str(ShufflingTime)]);
       disp(['    - randperm  etime (seconds): ' num2str(RandPermTime)]);
       fprintf('SHUFFLING TIME IS %.1f%% of RANDPERM TIME\n', 100.0 * ShufflingTime / RandPermTime);

    end
%}

%{
    % Computation time of shuffling and randperm.
    % Now the sample size is chosen at random, between 1 and 1000000.
    % Note that results can differ between MATLAB releases. See below.
    stoc = 0; rtoc = 0; loops = 100; n = zeros(100,1);
    for i=1:loops
        n(i) = randi(1000000 , 1);
        %n(i) = floor(1000000*abs(randn));
        x = randi(1000000 , n(i) , 1);
        nn=numel(x);

        st = tic;
        xperm1 = shuffling(x);
        stoc = stoc+toc(st);

        rt = tic;
        ix = randperm(nn);
        xperm2 = x(ix);
        rtoc = rtoc+toc(rt);

        clear('xperm1','ix','xperm2'); % Necessary to suppress JIT acceleration, for realistic times
    end
    disp(['shuffling etime (seconds): ' num2str(stoc)]);
    disp(['randperm  etime (seconds): ' num2str(rtoc)]);
    fprintf('==> SHUFFLING TIME IS %.1f%% of RANDPERM TIME\n', 100.0 * stoc / rtoc);

    % Results on R2016b
    %     shuffling etime (seconds): 4.5303
    %     randperm  etime (seconds): 5.3804
    %     ==> SHUFFLING TIME IS 84.2% of RANDPERM TIME

    % Results on R2012a
    % shuffling etime (seconds): 7.9629
    % randperm  etime (seconds): 4.9526
    % ==> SHUFFLING TIME IS 160.8% of RANDPERM TIME

    % Results on R2009a
    % shuffling etime (seconds): 8.4239
    % randperm  etime (seconds): 9.5947
    % ==> SHUFFLING TIME IS 87.8% of RANDPERM TIME

%}


%% Beginning of code
%if iscolumn(x) ; x = x' ; end;
n = numel(x);
I = n:-1:1;
J = ceil(rand(1,n) .* I);

for i=n:-1:1 
    % Note 1: to generate random integers between 1 and i, one by one in
    % the loop, would be less efficient. Such numbers are thus generated
    % outside the loop. To check, uncomment line below and profile the code.
    % j = ceil(rand() * i);
    % Note 2: this would be even less efficient: 
    % j = randi([1,i]);
    % Note 3: the double access to the vector on indices takes less time
    % than the allocation of the index in a separate variable with this line
    % j = J(i);
    
    % now, exchange element x[j] with x[i]. In MATLAB it would be
    % natural to do the exchange with a matrix operation, i.e. with 
    % x([i,j]) = x([j,i]).
    % Note that this switch forces MATLAB to pre-allocate a copy of x,
    % which takes time and, of corse, more memory. Using the MATLAB
    % profiler we found that the switch requires about 72% additional
    % time with respect to the following solution:
    t       = x(J(i));
    x(J(i)) = x(i);
    x(i)    = t;
end

end

%FScategory:UTICOMB
