function [kD , kW , kstar, varargout]  = quickselectFSw_demo(D,W,p)
% quickselectFSw_demo illustrates the functioning of quickselectFSw
%
%<a href="matlab: docsearchFS('quickselectFSw_demo')">Link to the help function</a>
%
%  Required input arguments:
%
%   D: Input data. Vector. A vector containing a set of $n$ data elements.
%      Data Types - double
%
%   W: Weights.  Vector. The vector contains a set of n positive weights
%      summing to 1. The 1-sum requirement does not imply loss of
%      generality: it is just to avoid extra checks and the normalization
%      step inside the function.
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
%         of the weighted order statistic in the data and weight vectors,
%         that is $D(k^{*})$ and $W(k^{*})$.
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
% See also: quickselectFSw
%
% References:
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('quickselectFSw_demo')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% quickselectFSw without optional parameter p gives the weighted
    % median. The median is 3, but the weighted one is 4, corresponding to
    % the weight 0.3.
    A = [1 2 3 4 5];
    W = [0.15 0.1 0.2 0.3 0.25];
    i=randperm(5);
    A=A(i); W=W(i);
    [kD, kW , kstar] = quickselectFSw_demo(A,W);
%}

%{
    % quickselectFSw for computing a generic weighted percentile
    A=[1 2 3 4 5];
    W=[0.15 0.1 0.2 0.3 0.25];
    i=randperm(5);
    A=A(i); W=W(i);
    perc = 25;
    p    = perc/100;
    [kD, kW , kstar] = quickselectFSw_demo(A,W,p);
%}

%{
    % quickselectFSw when the weights are all equal.
    % The weighted order statistic may logically reduce to the simple
    % order statistic, but pay attention to the different wy to treat even
    % and odd arrays.

    n=10;
    D=randperm(n);
    W=ones(1,n)./n;
    p=0.2;

    kDw = quickselectFSw_demo(D,W,p);

%}

%% Beginning of code

% 0. Input data size, sentinels and position as in quickselectFS
n        = length(D);
left     = 1;
right    = n;
position = -1;

% 1. Default is to compute the weighted median
if nargin<3 || p<0 || p>1
    p=0.5;
end

% 2. Data values and weights collapse in a nx2 column vector and move in pairs
D = [D(:) , W(:)];

% 3. The pivot k is set to get elements partially ordered in D(1:k,:)
k = ceil(n*p);

% 4. The external loop checks the condition on weights (point 6),
%    which generalises the ideas in Bleich and Overton (1983)
BleichOverton = true;
hl = []; hr = []; hp = []; hk = [];
%figure; sbD = subplot(2,1,1); sbW = subplot(2,1,2);
while BleichOverton
    
    %PLOT ---------------------------
    %set(gcf,'CurrentAxes',sbD);
    hhh = cellplotFS(D,k);
    [hl,hr,~,hk]=frecce(left,right,hhh,hl,hr,hp,hk);
    %set(gcf,'CurrentAxes',sbW);
    %hhhW = cellplotFS(D(:,2),k);
    %PLOT ---------------------------
    
    %% The internal loop is like in quickselectFS %%
    while (position~=k)
        
        %PLOT ---------------------------
        hhh=cellplotFS(D,k,right);
        [hl,hr,hp,hk]=frecce(left,right,hhh,hl,hr,hp,hk);
        %PLOT ---------------------------
        
        pivot      = D(k,:);
        D(k,:)     = D(right,:);
        D(right,:) = pivot;
        
        position   = left;
        for i=left:right
            
            %PLOT ---------------------------
            hhh=cellplotFS(D,i,position);
            [hl,hr,hp,hk]=frecce(left,right,hhh,hl,hr,hp,hk);
            %PLOT ---------------------------
            
            if (D(i,1)<pivot(1,1))
                
                %buffer = D(i,:);
                %D(i,:) = D(position,:);
                %D(position,:) = buffer;
                % the individual swap is much faster than the block swap above
                for s=1:2
                    buffer = D(i,s);
                    D(i,s) = D(position,s);
                    D(position,s) = buffer;
                end
                
                position = position+1;
                
            end
        end
        
        %PLOT ---------------------------
        hhh=cellplotFS(D,right,position);
        [hl,hr,hp,hk]=frecce(left,right,hhh,hl,hr,hp,hk);
        %PLOT ---------------------------
        
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

%PLOT ---------------------------
hhh=cellplotFS(D,k);
frecce(left,right,hhh,hl,hr,hp,hk);
%PLOT ---------------------------

title(['Weight balance for P= ', num2str(p) ': left = ' num2str(sum(D(1:k-1,2))) ' -- right = ' num2str(sum(D(k+1:end,2)))],'Fontsize',16,'Interpreter','Latex');

        
%% sub-functions

    function hhh=cellplotFS(C,p1,p2)
        
        C=C'; %C=C(:)';
        [nn,vv] = size(C);
        hhh = cellplot(num2cell(round(C,3)));
        for iii=1:nn*vv
            hhh(2*iii).FaceColor='none';
        end
        for iii=1:vv
            hhh(vv*2+2*iii).FaceColor=[0.9 0.9 0.9];
        end
        
        if nargin==2
            hhh(2*p1).FaceColor=[0.7 0.7 0.7];
        else
            hhh(2*p1).FaceColor='green';
        end
        hhh(2*p1+1).FontWeight='bold';
        
        if nargin==3
            if p1==p2
                hhh(2*p2).FaceColor='yellow';
            else
                hhh(2*p2).FaceColor='cyan';
            end
            hhh(2*p2+1).FontWeight='bold';
        end
        
        pause(0.51);
    end

    function [hl,hr,hp,hk]=frecce(left,right,H,hhl,hhr,hhp,hhk)
        if nargin > 3
            delete([hhl,hhr,hhp]);
        end
        if nargin == 7 && ~isempty(hhk)
            delete(hhk);
        end
        nn = (size(H,1) - 1)/4;
        pos = gca; pos=pos.InnerPosition;
        posl = pos(1); sizs = pos(3)/nn;
        posL = posl+sizs*left-sizs/2;
        posR = posl+sizs*right-sizs/2;
        posk = posl+sizs*k-sizs/2;
        hk=annotation('textarrow',[posk posk],[0.85 0.80],'String','k','Interpreter','latex','Fontsize',18);
        if position>0
            posP = posl+sizs*position-sizs/2;
        end
        hl=annotation('textarrow',[posL posL],[0.2 0.25],'String','L','Interpreter','latex','Fontsize',18);
        hr=annotation('textarrow',[posR posR],[0.2 0.25],'String','R','Interpreter','latex','Fontsize',18);
        if position>0
            hp=annotation('textarrow',[posP posP],[0.10 0.25],'String','p','Interpreter','latex','Fontsize',18);
        else
            hp = [];
        end
    end

end
%FScategory:GUI
