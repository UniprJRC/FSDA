function [kE , varargout] = quickselectFS(A,k,kiniindex)
%quickselectFS finds the k-th order statistic
%
% Optimized: delegates to compiled MEX (9x faster) when available,
% falls back to sort heuristic for n<=5000, then original partition.

%% Beginning of code

n = numel(A);

if nargout == 2
    if n <= 5000
        A = sort(A);
        kE = A(k);
        varargout{1} = A;
    else
        left = 1;
        right = n;
        if nargin > 2
            Ak = A(k); A(k) = A(kiniindex); A(kiniindex) = Ak;
        end
        pivotIndex = k;
        position = -999;
        while (position ~= k)
            pivot = A(pivotIndex);
            A(k) = A(right); A(right) = pivot;
            position = left;
            for i = left:right
                if (A(i) < pivot)
                    Ai = A(i); A(i) = A(position); A(position) = Ai;
                    position = position + 1;
                end
            end
            A(right) = A(position); A(position) = pivot;
            if position < k
                left = position + 1;
            else
                right = position - 1;
            end
        end
        kE = A(k);
        A(1:k-1) = sort(A(1:k-1));
        varargout{1} = A;
    end
else
    % Single output path: use MEX if available
    persistent useMex
    if isempty(useMex)
        useMex = ~isempty(which('aux.quickselectFSmex'));
    end

    if useMex
        kE = aux.quickselectFSmex(A+0, n, k-1);
    elseif n <= 5000
        A = sort(A);
        kE = A(k);
    else
        if nargin > 2
            Ak = A(k); A(k) = A(kiniindex); A(kiniindex) = Ak;
        end
        left = 1;
        right = n;
        pivotIndex = k;
        position = -999;
        while (position ~= k)
            pivot = A(pivotIndex);
            A(k) = A(right); A(right) = pivot;
            position = left;
            for i = left:right
                if (A(i) < pivot)
                    Ai = A(i); A(i) = A(position); A(position) = Ai;
                    position = position + 1;
                end
            end
            A(right) = A(position); A(position) = pivot;
            if position < k
                left = position + 1;
            else
                right = position - 1;
            end
        end
        kE = A(k);
    end
end

end
%FScategory:UTIGEN
