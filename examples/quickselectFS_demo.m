function kE = quickselectFS_demo(A,k,kiniindex)
%quickselectFS_demo illustrates the functioning of quickselectFS
%
%<a href="matlab: docsearchFS('quickselectFS_demo')">Link to the help function</a>
%
% Required input arguments:
%
%   A:  a set of unique numbers. Vector. Vector containing a set of n (distinct) numbers.
%   k:  order statistic index. Scalar. An integer between 1 and n indicating the
%       desired order statistic.
%                 Data Types - double
%
% Optional input arguments:
%
%  kiniindex: Index of an element in A. Scalar.
%      The index of an element in A that is supposed to be "close"
%      to the desired k-th order statistic. This information is used to
%      choose the pivot so that the chance to fall into the worst case
%      performance ($O(n^2)$) is minimized and the average case performance
%      is maximized.
%      Example - 'kiniindex',1
%      Data Types - double
%
% Output:
%
% kE : k-th order statistic. Scalar. Element in A that is larger than exactly k - 1 other elements of A.
%
%
% See also:  quickselectFS
%
% References:
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('quickselectFS_demo')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-01-15 11:53:08 #$: Date of the last commit

% Examples:
%
%{
    %% quickselectFS_demo with all default options.
    rng('default');
    rng(12345);
    n=15;
    Y=1:n; Y=shuffling(Y);
    
    k=7;
    [out]=quickselectFS_demo(Y,k);
    % Check the result
    sorY=sort(Y);
    disp([out,sorY(k)])
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

%PLOT ---------------------------
hhh=cellplotFS(A,k);
hp = [];
[hl,hr]=frecce(left,right,hhh);
%PLOT ---------------------------

while (position ~= k)
    %while ((left < right) && (position ~= k))
    
    %PLOT ---------------------------
    hhh=cellplotFS(A,k,right);
    [hl,hr,hp]=frecce(left,right,hhh,hl,hr,hp);
    %PLOT ---------------------------
    
    % Swap right sentinel and pivot element
    pivot    = A(pivotIndex);
    A(k)     = A(right);
    A(right) = pivot;
    
    position = left;
    for i = left:right
        
        %PLOT ---------------------------
        hhh=cellplotFS(A,i,position);
        [hl,hr,hp]=frecce(left,right,hhh,hl,hr,hp);
        %PLOT ---------------------------
        
        if(A(i)<pivot)
            if i~=position
                % Swap A(i) with A(position)!
                % A([i,position])=A([position,i]) would be more elegant but slower
                Ai          = A(i);
                A(i)        = A(position);
                A(position) = Ai;
            end
            position=position+1;
        end
    end
    
    %PLOT ---------------------------
    hhh=cellplotFS(A,right,position);
    [hl,hr,hp]=frecce(left,right,hhh,hl,hr,hp);
    %PLOT ---------------------------
    
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

%PLOT ---------------------------
hhh=cellplotFS(A,k);
frecce(left,right,hhh,hl,hr,hp);
%PLOT ---------------------------


%% sub-functions

    function hhh=cellplotFS(C,p1,p2)
        
        C=C(:)';
        nn = size(C,2);
        hhh = cellplot(num2cell(C));
        for iii=1:nn
            hhh(2*iii).FaceColor='none';
        end
        %hhh(2*k).FaceColor=[0.7 0.7 0.7];
        
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

    function [hl,hr,hp]=frecce(left,right,H,hhl,hhr,hhp)
        if nargin > 3
            delete([hhl,hhr,hhp]);
        end
        nn = (size(H,1) - 1)/2;
        pos = gca; pos=pos.InnerPosition;
        posl = pos(1); sizs = pos(3)/nn;
        posL = posl+sizs*left-sizs/2;
        posR = posl+sizs*right-sizs/2;
        posk = posl+sizs*k-sizs/2;
        annotation('textarrow',[posk posk],[0.65 0.60],'String','k','Interpreter','latex','Fontsize',18);
        if position>0
            posP = posl+sizs*position-sizs/2;
        end
        hl=annotation('textarrow',[posL posL],[0.4 0.45],'String','L','Interpreter','latex','Fontsize',18);
        hr=annotation('textarrow',[posR posR],[0.4 0.45],'String','R','Interpreter','latex','Fontsize',18);
        if position>0
            hp=annotation('textarrow',[posP posP],[0.25 0.35],'String','p','Interpreter','latex','Fontsize',18);
        else
            hp = [];
        end
    end


end

%FScategory:UTIGEN
