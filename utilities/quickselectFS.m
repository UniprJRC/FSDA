function kE = quickselectFS(A,k,kiniindex)
%quickselectFS finds the k-th order statistic
%
%<a href="matlab: docsearchFS('quickselectfs')">Link to the help function</a>
%
% Required input arguments:
%
%   A:  vector contatining a set of n (distinct) numbers.
%   k:  an integer between 1 and n indicating the desired order statistic.
%
%
%
% Optional input arguments:
%
%  kiniindex: Index of an element in A. Scalar.
%      The index of an element in A that is supposed to be "close"
%      to the desired k-th order statistic. This information is used to
%      choose the pivot so that the chance to fall into the worst case
%      performance (O(n^2)) is minimized and the average case performance
%      is maximized.
%      Example - 'kiniindex',1 
%      Data Types - double
%
% Output:
%
% - kE = element in A that is larger than exactly k - 1 other elements of A.

% References:
%
%       Riani M., Perrotta D. and Cerioli (2013), The Forward Search for
%       Very Large Datasets, submitted
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('quickselectFS')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:
%
%{
    %% quickselectFS with all default options
    n=200;
    Y=randn(n,1);
    k=10;
    [out]=quickselectFS(Y,k);
    % Check the result
    sorY=sort(Y);
    disp(out-sorY(k))
%}

%% Beginning of code

% Initialise the two sentinels
left    = 1;
right   = numel(A);
 
% if we know that element in position kiniindex is "close" to the desired order
% statistic k, than swap A(k) and A(kiniindex).
if nargin>2
    Ak    = A(k);
    A(k)  = A(kiniindex);
    A(kiniindex) = Ak;
end
    
% pivot is chosen at fixed position k. 
pivotIndex = k;

position = -999;
while ((left < right) && (position ~= k))
    
    pivot = A(pivotIndex);
    
    % Swap right sentinel and pivot element
    A(k)     = A(right);
    A(right) = pivot;
    
    position = left;
    for i = left:right
        if(A(i)<pivot)
            % Swap A(i) with A(position)!
            % A([i,position])=A([position,i]) would be more elegant but slower
            Ai          = A(i);
            A(i)        = A(position);
            A(position) = Ai;
            
            position=position+1;
        end
    end
    
    % Swap A(right) with A(position)
    A(right)    = A(position);
    A(position) = pivot;
    
    if  position < k
        left  = position + 1;
    else % this is 'elseif pos > k' as pos == k cannot hold (see 'while')
        right = position - 1;
    end
    
    % Pivot: extension to random choice has to be studied.
    %pivotIndex = ceil(( left + right ) / 2);
end

kE=A(k);

end