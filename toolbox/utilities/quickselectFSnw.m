function [kD , kW , kstar, varargout]  = quickselectFSnw(D,W,p)
% quickselectFSnw extends quickselectFSw to negative weights
%
%<a href="matlab: docsearchFS('quickselectFSnw')">Link to the help function</a>
%
% quickselectFSnw introduces the approach of Arce (1998) in quickselectFSw,
% to cope with possible negative weights in array W. The input/output
% arguments are the same of quickselectFSw. If the weights are all
% positive, the overhead for the extra checks required by Arce (1998) is
% almost negligible, especially when the sample size is large (~10^6).   
%
%  Required input arguments:
%
%   D: Input data. Vector. A vector containing a set of $n$ data elements.
%      Data Types - double
%
%   W: Weights.  Vector. The vector contains a set of $n$ arbitrary weights.
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
% The approach by Arce (1998) consists in the application of a standard
% algorithm for the weighted median to a data vector consisting in the
% "W-signed" observations (sign(W) .* D) weighted by the absolute values of
% the original weights (abs(W)). Therefore, in principle, it is possible to
% use quickselectFSw after the simple transformation of the input arrays;
% then, the weighted median position in the original arrays is the same of
% the one found by quickselectFSw on the transformed ones (the arguments in
% MATLAB are passed by value, therefore the input arrays do not change).
% This is shown in the examples below. 
% 
% Function quickselectFSnw follows a more elegant approach: it introduces a
% boolean array representing the sign of the original weights and uses a
% cmputationally efficient way to keep track of their position during the
% swaps (it uses the 'xor' operator to check if the weights to swap have
% different sign, and the 'not' to emulate a swap).
%
% Given that, in general, the potential presence of negative weights
% depends on the applicaation and is therefore known in advance, we decided
% to embed the new general approach in a copy of quickselectFSw to not
% introduce in the original function operations that may slow down
% computation in the standard positive weights case.
%
% As for quickselectFSw, we also provide the mex counterpart starting from
% code written in FORTRAN in collaboration with Ian Barrodale. The mex file
% is quickselectFSnwmex.
%
% See also: quickselectFSw
%
% References:
%
%   Arce, G.R. (1998), A general weighted median filter structure admitting
%   negative weights, "IEEE Transactions on Signal Processing, doi:
%   10.1109/78.735296", vol. 46, no. 12, pp. 3195-3205.
%
%   Arce, G.R. (2002), Recursive Weighted Median Filters Admitting Negative
%   Weights and Their Optimization, IEEE Transactions on Signal Processing,
%   Vol. 48, No. 3, pp. 768-779.
%
%   Azzini, I., Perrotta, D. and Torti, F. (2023), ﻿A practically efficient
%   fixed-pivot selection algorithm and its extensible MATLAB suite,
%   "arXiv, stat.ME, eprint 2302.05705"
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('quickselectFSnw')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Ian dataset, 15 January 2023
    
    p    = 0.5;
    X    = [-0.3 ; 0.1; 2.1 ; 0.3; 0.233333333; 0.366666667; 0.1; 0.3; -0.871428571; -0.85 ];
    W    = [-0.25; 0.5; 0.25; 0.5; 0.75; 0.75; 1.25; 1.5; -1.75; -2 ];

    [swm, swwm , kstar, XWS] = quickselectFSnw(X,W,p);		
    	
    rows_in_X = find(X==swm & swwm==swwm);	
    	

%}


%{
    %% quickselectFSnw with negative weights. 
    %
    % Execution with separate calls. The example is taken from Arce (1998).
    
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

    % Same result obtained with quickselectFSnw

    [swm1, swwm1 , kstar1] = quickselectFSnw(A,W,p);

%}

%{
    % quickselectFSnw with negative weights - example 2. 
    % This is another case with negative weights taken from Arce (2002).
    % Note that the sample now is even, then Arce (2002) takes the mean of
    % two middle values with an approach that uses sorting, which is 
    % therefore less efficient than applying twice quickselectFSnw and
    % taking the mean of the two middle order statistics.
    
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

    % Same result obtained with quickselectFSnw

    [swm1, swwm1 , kstar1] = quickselectFSnw(A,W,p);
%}


%% Pre-processing to apply Arce (1998)

% Arce(1998) step 1a : take the sign of weights
sW = sign(W); 
% Arce(1998) step 1b : take the absolute value of the weights
W  = abs(W);
% Arce(1998) step 1c: take the "W-signed" observations 
D  = sW .* D;

% boolean array representing the position of the positive weights: it is
% used to keep track of their position during the swaps
sWb  = logical(sW==1);

% The next block applies quickselectFSw on the transformed data and weights

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
sumW = sum(W); tol = eps(sumW) * 100; %TO DO: ADOPT IN MEX
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

        if xor(sWb(k),sWb(right))  % swap the weights sign if different
            sWb(k)      = not(sWb(k));
            sWb(right)  = not(sWb(right));
        end
        
        position   = left;
        for i=left:right
            if D(i,1)<pivot(1)
                
                for s=1:2
                    buffer = D(i,s);
                    D(i,s) = D(position,s);
                    D(position,s) = buffer;
                end

                if xor(sWb(i),sWb(position))  % swap the weights sign if different
                    sWb(i)          = not(sWb(i));
                    sWb(position)   = not(sWb(position));
                end
                             
                position = position+1;
                
            end
        end
        
        D(right,:)    = D(position,:);
        D(position,:) = pivot;

        if xor(sWb(right),sWb(position))  % swap the weights sign if different
            sWb(right)      = not(sWb(right));
            sWb(position)   = not(sWb(position));
        end
        
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
    if Le-p<=tol && p-Le-D(k,2)<=tol 
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

%% Post-processing, following Arce (1998)
% Assign to the weighted order statistic and the corresponding weight their
% original sign

% Reconstruct original sign at position kstar -- Maps: 1->1, 0->-1
original_sign = 2*sWb(kstar) - 1;    
%original_sign = (sWb(kstar) - not(sWb(kstar))); % ori

%   3b: final order statistic, with the initial sign 
kD            = kD * original_sign;
%   3c: the weight of the order statistic, with the initial sign
kW            = kW * original_sign;

%% optional step to partially sort the array before the required order statistic 
if nargout == 4
    %D(1:kstar-1,:)  = sortrows(D(1:kstar-1,:));

    D = [D(:,1) .* (2*sWb - 1) , D(:,2) .* (2*sWb - 1)];
    %D = [D(:,1) .* (sWb - not(sWb)) , D(:,2) .* (sWb - not(sWb))]; % ori
    varargout{1} = D;
end

%FScategory:UTIGEN
