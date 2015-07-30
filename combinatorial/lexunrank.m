function [kcomb,calls]=lexunrank(n,k,N,varargin)
%lexunrank gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations
%
%<a href="matlab: docsearchFS('lexunrank')">Link to the help function</a>
%
% Required input arguments:
%
%    n:  Number of elements. A non negative integer > k.
%        Data Types - single|double
%
%    k:  Items to choose from the set of n elements. A non negative integer.
%        Data Types - single|double
%
%    N:  Position N in the reverse co-lexicographic order of such
%        combinations. A non negative integer between 0 and bc(n,p)-1.
%        Data Types - single|double
%
% Optional input arguments:
%
%   the Pascal matrix as given by the MATLAB function pascal(n).
%   In applications where lexunrank is called many times, it is preferable
%   to compute the Pascal matrix once outside lexunrank, and pass it
%   to lexunrank as optional argument. Otherwise, the required binomial
%   coeffients are computed inside lexunrank using function bc and, when
%   possible, using the traditional recurrent formula.
%
% Output:
%
%   kcomb : Vector of length k. The $k$-combination of n elements at position
%           N. The position is relative to a reverse co-lexicographic order
%           of the combinations or, equivalently, of position bc(n,k)-N in
%           the lexicographic order of the same combinations.
%           Data Types - single|double
%
%   calls : Scalar. The number of binomial coefficients used to compute the
%           $k$-combination. 
%           Data Types - single|double
%
% More About:
%
%   REMARKS ON THE INPUT ARGUMENTS.
%
%   Input checks are intentionally avoided, as lexunrank is supposed to be
%   called many times, for sampling subsets. Thus, please ensure that:
%   - k < n;
%   - N is an integer between 0 and bc(n,p)-1.
%   It is possible to enable checks, by changing an internal "if" statement to 1.
%
%   REMARKS ON THE OUTPUT ARGUMENTS.
%
%   As $n$ increases, 'calls' becomes much smaller than 'ncomb'. This means
%   that lexunrank(n,k,N) is extremely convenient if you are interested in
%   one or several, but not all, $k$-combinations at given generation
%   order(s) N.
%
%   To generate all combinations in lexicographic order, it is more 
%   convenient using the FSDA function combsFS. The MATLAB function
%   with the same purpose, nchoosek(1:4,3), is much less efficient.
%
%   ON THE LEXICOGRAPHIC ORDERING.
%
%   lexunrank(n,k,N) gives the $k$-combination of n elements of position N
%   in the reverse co-lexicographic order of such combinations or,
%   equivalently, of position bc(n,k)-N in the lexicographic order of the
%   same combinations.
%   
%   Note that, in this implementation of the lexicographic unrank, N ranges
%   over the integers between 0 and bc(n,k)-1. For details see the
%   "combinatorial number system" discussed by Knuth (2005), pp. 5-6.
%
%   To clarify with an example the meaning of the different orders, while
%   the lexicographic order of the 2-combinations of 3 elements are:
%   
%   \[ 
%     \left( 
%        \begin{array}{ccc}
%           1  &   2  &   3     \\
%           1  &   2  &   4     \\
%           1  &   3  &   4     \\
%           2  &   3  &   4 
%        \end{array} 
%      \right)
%   \] 
%
%   the co-lexicographic order of the same combinations are
%   
%   \[ 
%     \left( 
%        \begin{array}{ccc}
%           3   &  2  &   1     \\
%           4   &  2  &   1     \\
%           4   &  3  &   1     \\
%           4   &  3  &   2
%        \end{array} 
%      \right)
%   \] 
%   
%   and the reverse co-lexicographic order of the original combinations are:
%
%   \[ 
%     \left( 
%        \begin{array}{ccc}
%           4   &  3  &   2     \\
%           4   &  3  &   1     \\
%           4   &  2  &   1     \\
%           3   &  2  &   1     
%        \end{array} 
%      \right)
%   \] 
%
%   The reasons for choosing a co-lexicographic unrank is that right-to-left 
%   array filling is much faster and elegant. The reverse is due to a similar 
%   motivation.
%
%
%   ALGORITMIC DETAILS.
%
% Given the totally ordered set $S=\{1,2,\ldots,n\}$, a $k$-combination is
% a subset $\{x_1, \ldots, x_k\}$ of $S$. Consider the $n$-lists of
% elements of the set $\{0,1\}$, i.e. the vertices of the hypercube $V_n$.
% Each $k$-combination $\{x_1,\ldots,x_k\}$ can be associated to the
% $n$-list having a 1 at position $x_1$, \ldots, $x_k$, and a 0 elsewhere.
%
% Example:
%   2-combinations of $\{1,2,3,4\}$: $\{1,2\}$, $\{1,3\}$, $\{1,4\}$,
%   $\{2,3\}$, $\{2,4\}$, $\{3,4\}$. Corresponding 4-lists of $\{0,1\}$:
%   $1100$,  $1010$,  $1001$,  $0110$, $0101$,  $0011$.
%
% The $n$-lists of $\{0,1\}$ containing $k$ times 1, and therefore
% equivalently the $k$-combinations of $n$-elements of $S$, can be
% generated in lexicographic order with an algorithm that builds the
% $k$-list of position $t+1$ using only the $k$-list of position $t$, and
% which stops without counting the combinations generated. For example, the
% MATLAB function NCHOOSEK(S,k), where $S$ is the row vector of length $n$
% of the elements of $S$, creates in lexicographic order a $k$ columns
% matrix whose rows consist of all possible combinations of the $n$
% elements of $S$ taken $k$ at a time. The number of such combinations,
% given by the binomial coefficient $n!/((n-k)! k!)$, can be also computed
% with the function NCHOOSEK by replacing the first argument, the row
% vector $S$, with the scalar $n$.
%
% Unfortunately the binomial coefficient increases rapidly with $n$, which
% makes the generation of all $k$-combinations computationally hard: with
% NCHOOSEK the task is impractical even for values just above 15. However,
% a lexicographic algorithm implements a one-to-one correspondence between
% the $k$-combinations and the generation order, i.e. the set of numbers $s
% = 1,\ldots,(n!/((n-k)!k!))$. This fact is used in our function to
% determine the $n$-list corresponding to the $k$-combination $\{x_1,
% \ldots, x_k\}$ which would be generated by the lexicographic algorithm at
% a given desired position $N$. This is useful in a number of applications
% which require one or several, but not all, $k$-combinations at given
% generation order(s).
%
% See also: combsFS, nchoosek, bc
%
% References:
%
%   Lehmer, D. H. (1964). The machine tools of combinatorics. In E. F.
%   Beckenbach (Ed.), Applied Combinatorial Mathematics, pp. 5--31. New York, Wiley.
%
%   Knuth, D. (2005). Generating All Combinations and Partitions. The Art of
%   Computer Programming, Vol. 4, Fascicle 3. Reading, Mass., Addison-Wesley.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('lexunrank')">Link to the help function</a>
%
% Last modified 06-Feb-2015
%
% Examples:

%{
        %% 7th 2 combination chosen among 5 element.
        n = 5; 
        k = 2; 
        N = 7;
        kcomb=lexunrank(n,k,N)
%}

%{
        %% number of binomial coefficient calls necessary to compute the 7th 2 combination chosen among 5 element.
        [~,calls]=lexunrank(n,k,N)
%}

%{
        %% 7th 2 combination chosen among 5 element, using the pascal matrix.
        [kcomb,calls]=lexunrank(n,k,N,pascal(n))
%}

%{
    % Additional example on the use of lexunrank.
    % Standard use.
    n = 4; p = 3;
    % number of p-combinations out of n
    n_bc = bc(n,p);
    % Pascal matrix
    pascalM=pascal(n);
    % n_bc is the Pascal cell in position (n-p+1,p+1)
    n_bc==pascalM(n-p+1,p+1)
    % all p-combinations in reverse-colex order generated by lexunrank
    % using a loop with rank integers ranging from 0 to bc(n,p)-1
    all_recolex = nan(n_bc,p);
    for N_lex = 0:n_bc-1
        all_recolex(N_lex+1,:) = lexunrank(n,p,N_lex);
    end
    all_recolex
%}

%{
    % Additional example on the use of lexunrank.
    % To change from reverse-colex to colex.
    all_colex = flipud(all_recolex)
    % and to change from colex to lex, it is sufficient this
    all_lex = fliplr(all_colex)

    % all p-combinations in lexi order generated using combsFS
    all_lex_combs = combsFS(1:n,p)

    % the combination at Lexi position N_lex=3 is generated by lexiunrank
    % in Colex position
    N_lex = 3; N_colex = n_bc - N_lex ;
%}

%{
    % Additional example on the use of lexunrank.
    % Use of lexunrank with pascal matrix
    kcomb = lexunrank(n,p,N_colex,pascal(n))
    % This is without Pascal matrix
    kcomb2 = lexunrank(n,p,N_colex)
    % Just as confirmation, the combination in the lexi order is
    all_lex_combs(N_lex,:)
%}

%% initialization

% REMARK: checks and unnecessary computations are intentionally avoided, as
% this function in FSDA is supposed to be called many times, for sampling
% subsets. To enable checks change the if statement to 1.

if 0
    if k >= n %#ok<UNRCH>
        error('FSDA:lexunrank:Wrongk','k must be < n');
    end
    
    % rank numbers range over the integers in [0 ncomb-1] (see combinatorial
    % number system in Knuth).
    if N >= bc(n,k) || N < 0  % note that bc(n,k) == pascalM(n-k+1,k+1)
        N=0;
    end
end

% initialise the row vector for the k-combination to unrank from N
kcomb = zeros(1,k);

% initialise the count of the calls to binomial coefficient values (via
% call to bc function or access to Pascal matrix cells)
calls = 0;

if isempty(varargin)
    
    %% call_bc OPTION:
    
    N_kk = N;
    pas_col = ones(n,1);
    seq = (1:n)';
    
    for kk = k:-1:1
        % The next 'if' statement builds the required part of column kk+1
        % of the pascal matrix, which is the argument of the 'find'
        % statement which follows.
        % This replaces the loop with repeated calling of bc:
        %         for x = kk:n-1
        %             if  bc(x,kk)> N_kk, break, end
        %         end
        if kk == k;
            for x2 = kk:n-1
                pas_col(x2+1) = pas_col(x2)*(x2+1)/(x2+1-kk);
            end
        else
            pas_col(kk+1:n) = pas_col(kk+1:n)*(kk+1)./(seq(kk+1:n)-kk);
        end
        x = find(pas_col(kk:end) > N_kk,1);
        
        if isempty(x)
            maxx=n-1;
            calls=calls+maxx-kk+1;
        else
            maxx = x+kk-2;
            calls=calls+maxx-kk+2;
        end
        
        kcomb(1,kk)=n-maxx;
        
        if maxx>=kk
            N_kk = N_kk - bc(maxx,kk);
            calls = calls+1;
        end
    end
    
else
    
    %% FAST OPTION:
    % binomial coefficients are taken from the pascal matrix rather than
    % computing them using bc. Of course this option is space greedy.

    pascalM=varargin{1};

    N_kk = N;
    
    for kk = k:-1:1
        
        x=find(pascalM(1:n-kk,kk+1) > N_kk , 1);
        
        if isempty(x) % || x1==n-kk
            maxx=n-1;
            calls=calls+maxx-kk+1;
        else
            maxx = x+kk-2;
            calls=calls+maxx-kk+2;
        end
        
        kcomb(1,kk)=n-maxx;
        
        if maxx>=kk
            N_kk = N_kk - pascalM(maxx-kk+1,kk+1);  % this is: N_kk - bc(maxx,kk)
            calls = calls+1;
        end
    end
    
end

end

%% This part is kept for history.

% FAST OPTION:
% The following loop with break and if statement have been replaced by the
% the use of a "find" instruction and subsequent "if" statement.
%
%         for x = kk:n-1
%             if (pascalM(x-kk+1,kk+1) > N_kk), break, end  % this is: bc(x,kk)>N_kk
%         end
%
%         if x==n-1
%             maxx = x;
%             calls=calls+maxx-kk+1;
%         else
%             maxx = x-1;
%             calls=calls+maxx-kk+2;
%         end

% call_bc OPTION:
% The following two if statements
%
%         if x==n;
%             x=x-1;
%         elseif x==1;
%             x=kk;
%         end
%
%         if x==n-1
%             maxx = x;
%             calls=calls+maxx-kk+1;
%         else maxx = x-1;
%             calls=calls+maxx-kk+2;
%         end
%
% have been replaced by this one:
%
%         if isempty(x)
%             maxx=n-1;
%             calls=calls+maxx-kk+1;
%         else
%             maxx = x+kk-2;
%             calls=calls+maxx-kk+2;
%         end

