function [y,i] = datasampleFS(x,k,wgts)

nargs = nargin;

dim = find(size(x)~=1, 1); % first non-singleton dimension
if isempty(dim), dim = 1; end
n = size(x,dim);

replace = false;

if nargs > 2
    
    if ~isempty(wgts)
        if ~isvector(wgts) || length(wgts) ~= n
            error('FSDA:datasample:InputSizeMismatch');
        else
            wgts = double(wgts);
            sumw = sum(wgts);
            if ~(sumw > 0) || ~all(wgts>=0) % catches NaNs
                error('FSDA:datasample:InvalidWeights');
            end
        end
    end
end

% Sample with replacement
if replace
    if n == 0
        if k == 0
            i = zeros(0,1);
        else
            error('FSDA:datasample:EmptyData');
        end
        
    elseif isempty(wgts) % unweighted sample
        
        i = randi(n,1,k);
        
    else % weighted sample
        
        p = wgts(:)' / sumw;
        edges = min([0 cumsum(p)],1); % protect against accumulated round-off
        edges(end) = 1; % get the upper edge exact
        [~, i] = histc(rand(1,k),edges);
    end
    
    % Sample without replacement
else
    if k > n
        error('FSDA:datasample:SampleTooLarge');
        
    elseif isempty(wgts) % unweighted sample
        
        i = randperm(n,k);
        
    else % weighted sample
        if sum(wgts>0) < k
            error('FSDA:datasample:TooFewPosWeights');
        end
        
        % REMARK: THIS IS THE CALL TO A MATLAB PRIVATE MEX
        i = wsworFS(wgts,k);
    end
    
end

% Use the index vector to sample from the data.
if ismatrix(x) % including vectors and including dataset or table
    if dim == 1
        y = x(i,:);
    elseif dim == 2
        y = x(:,i);
    else
        reps = [ones(1,dim-1) k];
        y = repmat(x,reps);
    end
else % N-D
    subs = repmat({':'},1,max(ndims(x),dim));
    subs{dim} = i;
    y = x(subs{:});
end
