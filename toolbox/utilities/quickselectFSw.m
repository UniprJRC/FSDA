function [kD , kW , kstar, varargout]  = quickselectFSw(D,W,p)
% quickselectFSw finds the 100*p-th weighted order statistic for $0<p<1$
%
%<a href="matlab: docsearchFS('quickselectFSw')">Link to the help function</a>
%
% quickselectFSw generalises the calculation of the weighted median to a
% generic percentile $100pth (0<p<1)$. It finds the p-th weighted order
% statistic in the elements of a data vector $D$ considering the associated
% weights $W$. More precisely, the algorithm finds the element $D_{k^{*}}$
% that make the following difference as small as possible:
% $\sum_{i=1}^{k^{*}_{p} - 1} w_{i} - \sum_{i=k^{*}_{p}+1}^{n} w_{i}$.
% The linear algorithm used for the computation extends quickselectFS.
% REMARK: we also provide the mex counterpart, quickselectFSwmex; see the
% last example to understand how it works. 
%
%  Required input arguments:
%
%   D: Input data. Vector. A vector containing a set of $n$ data elements.
%      Data Types - double
%
%   W: Weights.  Vector. The vector contains a set of n positive weights.
%      It is not necessary that the weights sum to 1. 
%      Data Types - double
%
% Optional input arguments:
%
%   p: 100pth percentile ($0<p<1$). Scalar. A number between 0 and 1
%      indicating the fraction of total weights that should be considered
%      in partitioning the associated input data. Default is p=0.5, leading
%      to the weighted median.
%      Example - p,0.2
%      Data Types - double
%
% Output:
%
% kD    : weighted order statistic in D. Scalar. Element kstar in vector D
%         ($D(k^{*})$), for which we have $\sum_{i=1}^{k^{*}-1} w_i<=p$
%         and $\sum_{j=k^{*}+1}^{n} w_j<=1-p$.
%
% kW    : weight associated to kD. Scalar. Element kstar in vector W. The
%         weight partition around kW is optimal in the sense that the sum
%         of the weights on its left is as close as possible to $p$ and
%         those on its right account for the remaining $1-p$.
%
% kstar : the index of kD in D and of kW in W. Scalar. It is the position
%         of the weighted order statistic in the internal data and weight 
%         vectors, that is $D(k^{*})$ and $W(k^{*})$.
%
%
% Optional Output:
%
%   Ds  : Output vector. Array. At the and of the computation the vector
%         [D(:) , W(:)] is appropriately partitioned around kstar. It is
%         returned ordered in the interval (1:kstar-1,:). Remark: this
%         operation sorts only part of the array, but it can still slow
%         down the algorithm.
%
% More About:
%
%   quickselectFSw builds on quickselectFS, the algorithm used to find
%   order statistics. Data and weights are collapsed to a single vector
%   D=[D(:),W(:)] (point 2). The core part of quickselectFS is called on D
%   assuming that the weighted percentile is in position k=ceil(n*p) (point
%   3). This step permutes vector D so that the elements in positions
%   (1:k-1,:) are smaller than the element at position k. At this point we
%   check if the array fulfills the weighted median condition by Bleich and
%   Overton (1985) (point 5). If so, we stop computation (point 6.). If
%   not, we apply again the core part of quickselectFS on either D(1:k,:)
%   with k+1, or D(k:n,:) with k-1. The iteration will remove or add
%   weights in order to approach the optimality condition (see points 7,8
%   and 9). The process is iterated till success.
%
% See also: quickselectFS
%
% References:
%
%   Bleich, C. and Overton, M.L. (1983), A linear-time algorithm for the
%   weighted median problem. Technical Report 75, New York University,
%   Courant Institute of Mathematical Sciences.
%
%   Azzini, I., Perrotta, D. and Torti, F. (2023), ﻿A practically efficient
%   fixed-pivot selection algorithm and its extensible MATLAB suite,
%   "arXiv, stat.ME, eprint 2302.05705"
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('quickselectFSw')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% quickselectFSw without optional parameter p, giving the weighted median.
    % The median is 3 but the weighted one is 4, corresponding to the weight 0.3.
    A = [1 2 3 4 5];
    W = [0.15 0.1 0.2 0.3 0.25];
    [kD, kW , kstar] = quickselectFSw(A,W);
%}

%{
    %% quickselectFSw with sum of weights taking any positive value.
    % The weights usually are forced to sum to 1. However, our algorithm
    % does not have this limitation.
    %
    % Case of sum of weights greater than 1.
    % This is equivalent to the previous example, but the weights have been
    % multiplied by 10. Therefore, results should be identical even if the
    % sum of the weights is not one. 
    A = [1 2 3 4 5];
    W = [0.15 0.1 0.2 0.3 0.25];
    [kD, kW , kstar] = quickselectFSw(A,W);
    W2 = [1.5 1 2 3 2.5];
    [kD2, kW2 , kstar2] = quickselectFSw(A,W2);
    %
    % Case of sum of weights positive, but smaller than 1
    W3 = [0.015 0.01 0.02 0.03 0.025];
    [kD3, kW3 , kstar3] = quickselectFSw(A,W3);
    W4 = [0.0015 0.001 0.002 0.003 0.0025];
    [kD4, kW4 , kstar4] = quickselectFSw(A,W4);
%}

%{
    %% quickselectFSw with negative weights. 
    % This is an atypical case, but may find interesting applications,
    % for example in the simplex algorithm for the choice of the pivot.
    % quickselectFSw cannot be applied directly with the negative weights, 
    % but a trivial change of the 
    % observations and weights can lead to the solution.
    %
    % The approach proposed in this example is taken from 
    % G. R. Arce, "A general weighted median filter structure admitting
    % negative weights," in IEEE Transactions on Signal Processing, vol. 46,
    % no. 12, pp. 3195-3205, Dec. 1998, doi: 10.1109/78.735296. 
    
    p    = 0.5;
    A    = [-2 2 -1 3 6];
    W    = [0.1 0.2 0.3 -0.2 0.1];
    
    % Arce(1998) step 1a: take the sign of weights
    sW   = sign(W);
    % Arce(1998) step 1b: take the absolute value of the weights
    absW = abs(W);
    % Arce(1998) step 1c: take the "W-signed" observations 
    sWA  = sW .* A;

    % Arce(1998) step 2: compute the weighted median of sWA with weights absW 
    [swm, swwm , kstar] = quickselectFSw(sWA,absW,p);

    % Step 3: retrieve the weighted median and corresponding weight
    %         with the right original sign
    %   3a: index of weighted median and corrisponding weight 
    kstar_AW = find(sWA==swm & absW==swwm);
    %   3b: final weighted median, with the initial sign 
    wm       = sW(kstar_AW) * swm;
    %   3c: the weight of the weighted median, with the initial sign
    wwm      = W(kstar_AW);

%}

%{
    % quickselectFSw with negative weights - example 2. 
    % This is another atypical case, with negative weights, taken from
    % Gonzalo R. Arce (2002) Recursive Weighted Median Filters Admitting
    % Negative Weights and Their Optimization, IEEE TRANSACTIONS ON SIGNAL 
    % PROCESSING, VOL. 48, NO. 3, MARCH 2000, pp 768-779.
    % Note that Arce (2002), being this an even sample, takes the mean 
    % of two middle values with a different approach using sorting.
    
    p    = 0.5;
    A    = [-2 2 -1 3 6 8];
    W    = [0.2 0.4 0.6 -0.4 0.2 0.2];
    
    % Arce(1998) step 1a: take the sign of weights
    sW   = sign(W);
    % Arce(1998) step 1b: take the absolute value of the weights
    absW = abs(W);
    % Arce(1998) step 1c: take the "W-signed" observations 
    sWA  = sW .* A;

    % Arce(1998) step 2: compute the weighted median of sWA with weights absW 
    [swm, swwm , kstar] = quickselectFSw(sWA,absW,p);

    % Step 3: retrieve the weighted median and corresponding weight
    %         with the right original sign
    %   3a: index of weighted median and corrisponding weight 
    kstar_AW = find(sWA==swm & absW==swwm);
    %   3b: final weighted median, with the initial sign 
    wm       = sW(kstar_AW) * swm;
    %   3c: the weight of the weighted median, with the initial sign
    wwm      = W(kstar_AW);
%}

%{
    % quickselectFSw for computing a generic weighted percentile.
    A=[1 2 3 4 5];
    W=[0.15 0.1 0.2 0.3 0.25];
    perc = 25;
    p    = perc/100;
    [kD, kW , kstar] = quickselectFSw(A,W,p);
%}

%{
    % quickselectFSw with the optional output array partially sorted.
    A=randperm(10);
    W=abs(randn(10,1));
    W=W/sum(W);
    perc = 30;
    p    = perc/100;
    [~, ~ , ~, As] = quickselectFSw(A,W,p);
    disp(As(:,1));
%}

%{
    %% quickselectFSw when the weights are all equal.
    % The weighted order statistic may logically reduce to the simple
    % order statistic, but pay attention to the different wy to treat even
    % and odd arrays.

    n=1000;
    D=randperm(n);
    W=ones(1,n)./n;
    p=0.2;

    % Very fast way to check if elements of vector are all equal
    allequal=true; e=2;
    while allequal && e<n
        if W(1)==W(e), e=e+1; else , allequal=false; end
    end

    if allequal
        % Option 1: apply the standard quickselectFS. Note that if the array
        % has an even number of elements, we need to apply quickselectFS
        % twice and take the average of the two contiguous order statistics.
        k = floor(p*n);
        kD = quickselectFS(D,k+1);
        if mod(n,2) == 0
            kD2 = quickselectFS(D,k);
            kD  = (kD+kD2)/2;
        end
    end

    % Option2: apply quickselectFSw. No need to check for equal weights,
    % but the weighted order statistic with percentile p may produce a
    % different result if the number of elements in the array is even.
    kDw = quickselectFSw(D,W,p);

    disp(['quickselectFS  -> ' num2str(kD)]);
    disp(['quickselectFSw -> ' num2str(kDw)]);
%}

%{
    %% Verify the definition of weighted order statistic.
    % Check the weighted median, that is p=0.5.
    % The weighted median is the value $x$ which minimizes the weighted
    % mean absolute deviation $\sum w_{i} | x_{i} - x |$.

    N  = 10000;
    D  = randperm(N);
    W  = abs(randn(1,N));
    W1 = W/sum(W);

    p = 0.5;
    [wm , ww , kstar, wdSort] = quickselectFSw(D,W,p);
    kstar_D = find(D==wm);
    kstar_W = find(W==ww);
    disp(['this is the weighted median  = ' num2str(wm)]);
    disp(['this is kstar                = ' num2str(kstar)]);
    disp(['this is the index of wm in D = ' num2str(kstar_D)]);
    disp(['this is the index of ww in W = ' num2str(kstar_W)]);

    % compute the objective function with the definition
    obj = sum(W.*abs(D-wm));

    disp('check if there are elements with smaller/equal objective function');
    check = false;
    for i = 1 : N
            obj_i = sum(W.*abs(D-D(i)));
            if (i ~= kstar_D) && (obj_i <= obj)
                disp(['this is i = ' num2str(i)  ' - this is D(i) = ' num2str(D(i))]);
                check = true;
            end
    end
    if check
        disp('optimal obj is not unique');
    else
        disp('optimal obj is unique');
    end
%}

%{
    % Use the mex function quickselectFSwmex.
    % REMARK 1: it is necessary to pass the number of data elements.
    % REMARK 2: it is necessary to pass a modified copy of the data and 
    % weight arrays as indicated in the example, as the function change the
    % order of the elements in the original arrays (variables are passed by
    % reference).

    N = 10000;
    D=randperm(N);
    W = abs(randn(N,1));
    W = W/sum(W);

    p = 0.5;
    [wm , ww , kstar, wdSort]  = quickselectFSw(D,W,p);

    % The next two lines are necessary to break the link between D and the
    % copy which will be passed by reference to quickselectFSwmex
    D_copy = D; D_copy(end+1)=999; D_copy(end)=[]; 
    W_copy = W; W_copy(end+1)=999; W_copy(end)=[]; 
    wm_mex = aux.quickselectFSwmex(D_copy,W_copy,0.5,N);

    disp('  ');
    disp(['this is wm      = ' num2str(wm)]);
    disp(['this is wm_mex  = ' num2str(wm_mex)]);

    % if zero, the sorted arrays are equal. 
    sum(W_copy    - wdSort(:,2))
    sum(D_copy(:) - wdSort(:,1))

%}

%% Beginning of code

% 0. Input data size, sentinels and position as in quickselectFS
n        = length(D);
left     = 1;
right    = n;
position = -1;

% 1. Default is to compute the weighted median
if nargin<3 || p<0 %|| p>1
    p=0.5;
end

% 2. Data values and weights collapse in a nx2 column vector and move in pairs
D = [D(:) , W(:)];

% 3a. The pivot k is set to get elements partially ordered in D(1:k,:)
k = ceil(n*p);

% 3b. p is modified to treat the general case when sum of weights > 0 
% but different than 1
sumW = sum(W);
p = p*sumW;

% 4. The external loop checks the condition on weights (point 6),  
%    which generalises the ideas in Bleich and Overton (1983)
BleichOverton = true;
while BleichOverton
    
    %% The internal loop is like in quickselectFS %%  
    while (position~=k)
        
        pivot      = D(k,:);
        D(k,:)     = D(right,:);
        D(right,:) = pivot;
        
        position   = left;
        for i=left:right
            if D(i,1)<pivot(1)
                
                for s=1:2
                    buffer = D(i,s);
                    D(i,s) = D(position,s);
                    D(position,s) = buffer;
                end
                             
                position = position+1;
                
            end
        end
        
        D(right,:)    = D(position,:);
        D(position,:) = pivot;
        
        if  (position < k)
            left  = position + 1;
        else
            right = position - 1;
        end
        
    end
    
    %% Checks on weights extends Bleich-Overton %%
    
    % 5. The algebra of Bleich-Overton is re-written to check the
    %    conditions on wheigts efficiently. When the conditions are
    %    met, D is partially ordered and the optimal solution is reached.
    
    Le=sum(D(1:k-1,2));
    if Le-p<=0 && p-Le-D(k,2)<=0
        % 6. The condition is met: stop computation.
        kD=D(k,1);
        kW=D(k,2);
        kstar=k;
        BleichOverton=false;
        
    else
        % 7. The conditions not met: go back to quickselectFS with new
        %    sentinels - (k,n) or (1,k) - and new order statistics - 
        %    k+1 or k-1 (see point 8 and 9.).
        
        if  D(k,2)<2*(p-Le)
            % 8. Need to add weight to reach the condition in point 5. 
            %    Add an element (weight and data) to the left part.
            k     = k+1;
            left  = k;
            right = n;
        else
            % 9. Here, we have that D(k,2)> 2*(p-Le). Need to remove an 
            % element (weight and data) from the right part 
            % in point 5.
            k     = k-1;
            left  = 1;
            right = k;
        end
    end
    position=-1;
end

if nargout == 4
    D(1:kstar-1,:)  = sortrows(D(1:kstar-1,:));
    varargout{1} = D;
end

%FScategory:UTIGEN
