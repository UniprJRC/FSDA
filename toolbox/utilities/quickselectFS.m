function [kE , varargout] = quickselectFS(A,k,kiniindex)
%quickselectFS finds the k-th order statistic
%
%<a href="matlab: docsearchFS('quickselectFS')">Link to the help function</a>
%
% quickselectFS is a linear algorithm equivalent to quickselect (Hoare's
% Find), which computes order statistics with an approach that avoids
% recursion and repeated calls to partitioning functions.
% REMARK: we also provide the mex counterpart, quickselectFSmex; see the
% last example to understand how it works. 

%
% Required input arguments:
%
%   A:  a set of unique numbers. Vector. Vector containing a set of n
%       (distinct) numbers.
%                 Data Types - double
%   k:  order statistic index. Scalar. An integer between 1 and n indicating 
%       the desired order statistic.
%                 Data Types - double
%
%
% Optional input arguments:
%
%  kiniindex: Index of an element in A. Scalar.
%      The index of an element in A that is supposed to be "close" to the
%      desired k-th order statistic. This information is used to choose the
%      pivot so that the chance to fall into the worst case performance
%      ($O(n^2)$) is minimized and the average case performance is
%      maximized.
%      Example - 'kiniindex',1 
%      Data Types - double
%
% Output:
%
% kE : k-th order statistic. Scalar. Element in A that is larger than
% exactly k - 1 other elements of A.
%
% Optional Output:
%    
%    Asor   : Partially sorted vector. Vector. Elements of input vector 
%             A(1:k) are sorted in ascending order. Remark: this option
%             implies the application of a sorting algorithm on part of the
%             array, with obvious implications on performance.
%
% See also:  FSMmmd
%
% References:
%
% Azzini, I., Perrotta, D. and Torti, F. (2023), ﻿A practically efficient
% fixed-pivot selection algorithm and its extensible MATLAB suite,
% "arXiv, stat.ME, eprint 2302.05705"
%
% Riani, M., Perrotta, D. and Cerioli, A. (2015), The Forward Search for
% Very Large Datasets, "Journal of Statistical Software"
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('quickselectFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:
%
%{
    %% quickselectFS with all default options.
    n=200;
    Y=randn(n,1);
    k=10;
    [out]=quickselectFS(Y,k);
    % Check the result
    sorY=sort(Y);
    disp([out,sorY(k)])
%}

%{
    %% quickselectFS with kiniindex supplied.
    n=200;
    Y=randn(n,1);
    k=10;
    % kiniindex is supplied 
    [out]=quickselectFS(Y,k,20);
    % Check the result
    sorY=sort(Y);
    disp([out,sorY(k)])
%}

%{
    %% quickselectFS with two output arguments.
    n=10;
    Y=randn(n,1);
    k=3;
    % kiniindex is supplied 
    [out,Ysor]=quickselectFS(Y,k);
    % Check the result
    disp([Y, Ysor])
%}

%{
    %% quickselectFS: worst case scenario: see circshift.
    n=10;
    Y=1:n;
    Y = circshift(Y,-1);
    k=n;
    out=quickselectFS(Y,k);
    disp(out);
%}

%{
    % Use the mex function quickselectFSmex.
    % REMARK 1: it is necessary to pass the number of data elements.
    % REMARK 2: if you are looking for the order statistic k, then you
    % should pass to the function k-1, as in C the array index starts from
    % zero. 
    % REMARK 3: it is necessary to pass a modified copy of the data  
    % as indicated in the example, as the function change the
    % order of the elements in the original array (array is passed by
    % reference).

    n=1000;
    Y=randn(n,1);
    k=100;

    [k_star,Ysor] = quickselectFS(Y,k);

    % The next line is necessary to break the link between Y and the
    % copy which will be passed by reference to quickselectFSmex.
    Y_copy = Y; Y_copy(end+1)=999; Y_copy(end)=[]; 

    outmex = quickselectFSmex(Y_copy,n,k-1);

    disp('  ');
    disp(['this is k_star      = ' num2str(k_star)]);
    disp(['this is k_star_mex  = ' num2str(k_star)]);

    % if zero, the sorted arrays are equal. 
    sum(Y_copy - Y)

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

% The original loop was:
% while ((left < right) && (position ~= k)) 
% The (left < right) condition reduces the number of iterations, but the
% gain is not compensated by the cost of the additional check. Therefore,
% we just check that position ~= k.
position = -999;
while (position ~= k)
    
    % Swap the right sentinel with the pivot
    pivot    = A(pivotIndex);
    A(k)     = A(right);
    A(right) = pivot;
    
    position = left;
    for i = left:right
        if(A(i)<pivot)
            % Swap A(i) with A(position)
            % A([i,position])=A([position,i]) is more elegant but slower
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
    else % --> 'elseif pos > k' as pos == k cannot occur (see 'while')
        right = position - 1;
    end
    
end

kE=A(k);

if nargout == 2
    A(1:k-1) = sort(A(1:k-1));
    varargout{1} = A;
end

end
%FScategory:UTIGEN
