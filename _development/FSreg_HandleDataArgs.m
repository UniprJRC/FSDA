function [X,y,haveDataset,otherArgs] = FSreg_HandleDataArgs(X,y,varargin) % [X, y | DS], ...
if isa(X,'dataset')
    X = dataset2table(X);
end

haveDataset = isa(X,'table');
if haveDataset
    % Have a dataset/table array, so no y.
    if nargin > 1
        % Put the second arg in with the rest
        otherArgs = [{y} varargin];
    else
        otherArgs = {};
    end
    y = []; % just put anything in y
elseif nargin < 2
    error(message('stats:classreg:regr:TermsRegression:MissingY'))
else
    if isrow(X)
        nx = length(X);
        if (isvector(y) && numel(y)==nx) || (size(y,1)==nx)
            X = X';
        end
    end
    isNumVarX = isnumeric(X) || islogical(X);
    isCatVecX = isa(X,'categorical') && isvector(X);
    if ~(isNumVarX || isCatVecX)
        error(message('stats:classreg:regr:FitObject:PredictorMatricesRequired'));
    end
    otherArgs = varargin;
end
end