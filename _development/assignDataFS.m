function model = assignDataFS(X,y,w,asCat,varNames,excl)
% This is a generic implementation suitable for column-oriented
% numeric/categorical data.  The data can be
%    * in one dataset/table array X containing matrix or vector variables (and y is []),
%    * in separate matrices X and y
% We don't check types, we leave that for the selectVariables method.
% varNames is required for data in X,y, ignored for data in a dataset/table array
if isa(X,'dataset')
    X = dataset2table(X);
end
haveDataset = isa(X,'table');

% Do some preliminary processing of the input data
if haveDataset && all(varfun(@ismatrix,X,'OutputFormat','uniform'))
    [nobs,nvars] = size(X);
    predLocs = 1:(nvars-1);
    respLoc = nvars;
elseif ismatrix(X) && ismatrix(y)
    % recognize vectors of either orientation, this assumes n > 1
    if isvector(X)
        X = X(:);  % force to a column
    end
    [nobs,p] = size(X); % p is number of predictors
    if ischar(X)
        p = 1;
    end
    if isvector(y)
        y = y(:);  % force to a column
    end
    if size(y,1) ~= nobs
        error(message('stats:classreg:regr:FitObject:PredictorResponseMismatch'));
    end
    nvars = p + 1;
    predLocs = 1:(nvars-1);
    respLoc = nvars;
else
    error(message('stats:classreg:regr:FitObject:MatricesRequired'));
end

% Process the exclusion vector
[excl,predLocs] = getDataVariableFS(excl,nobs,X,predLocs,respLoc,'exclude');
if islogical(excl)
    if numel(excl)~=nobs
        error(message('stats:classreg:regr:FitObject:BadExcludeLength',nobs));
    end
else
    if any(excl<0) || any(excl>nobs) || any(excl~=round(excl))
        error(message('stats:classreg:regr:FitObject:BadExcludeValues'));
    end
    tmp = excl;
    excl = false(nobs,1);
    excl(tmp) = true;
end
if all(excl)
    error(message('stats:classreg:regr:FitObject:AllExcluded'));
end

% Now compute variable ranges and other information
if haveDataset
    viName = X.Properties.VariableNames;
    viClass = varfun(@class,X,'OutputFormat','cell');
    
    % Provisionally choose the last variable as the response, and
    % the rest as predictors.
    viInModel = [true(1,nvars-1) false];
    
    viIsCategorical = varfun(@internal.stats.isDiscreteVar,X,'OutputFormat','uniform');
    
    % Never inside because AsCat is always empty
    if ~isempty(asCat)
        viIsCategorical = classreg.regr.FitObject.checkAsCat(viIsCategorical,asCat,nvars,true,viName);
    end
    viRange = cell(nvars,1);
    for i = 1:nvars
        vi = X.(viName{i});
        vir = getVarRangeFS(vi,viIsCategorical(i),excl);
        if iscellstr(vi) && isvector(vi) &&...
                (numel(unique(strtrim(vir)))< numel(vir))
            X.(viName{i}) = strtrim(X.(viName{i}));
            vir = getVarRangeFS(X.(viName{i}),viIsCategorical(i),excl);
        end
        viRange{i} = vir;
    end
    
    data = X;
    obsNames = X.Properties.RowNames;
elseif ismatrix(X) && ismatrix(y)
    % The Variables property puts the response last, match that
    viName = varNames;
    viClass = [repmat({class(X)},1,p) {class(y)}];
    
    % Provisionally choose all variables in X as predictors.
    viInModel = [true(1,p) false];
    
    viIsCategorical = [repmat(internal.stats.isDiscreteVar(X),1,p) internal.stats.isDiscreteVar(y)];
    if ~isempty(asCat)
        viIsCategorical = classreg.regr.FitObject.checkAsCat(viIsCategorical,asCat,nvars,false,viName);
    end
    viRange = cell(nvars,1);
    if ~any(viIsCategorical)
        % Handle special all-continuous case quickly
        viMax = max(X(~excl,:),[],1);
        viMin = min(X(~excl,:),[],1);
        temp = [viMin(:),viMax(:)];
        viRange(1:nvars-1) = mat2cell(temp,ones(size(temp,1),1),2);
    else
        for i = 1:(nvars-1)
            viRange{i} = getVarRangeFS(X(:,i),viIsCategorical(i),excl);
        end
    end
    viRange{end} = getVarRangeFS(y,viIsCategorical(end),excl);
    
    data = struct('X',{X},'y',{y});
    obsNames = {};
end

[w,predLocs] = getDataVariableFS(w,nobs,X,predLocs,respLoc,'weights');
if isempty(w)
    w = ones(nobs,1);
elseif any(w<0) || numel(w)~=nobs
    error(message('stats:classreg:regr:FitObject:BadWeightValues', nobs));
end

% We do not check for empty data, a class may want to allow models
% that are not fit to data, or may have some edge case behavior.
model.Data = data;
model.PredLocs = predLocs; % provisionally
model.RespLoc = respLoc;   % provisionally
model.VariableInfo = table(viClass(:),viRange(:),viInModel(:),viIsCategorical(:), ...
    'VariableNames',{'Class' 'Range' 'InModel' 'IsCategorical'}, ...
    'RowNames',viName);
model.ObservationInfo = table(w,excl,false(nobs,1),false(nobs,1), ...
    'VariableNames',{'Weights' 'Excluded' 'Missing' 'Subset'},...
    'RowNames',obsNames);
model.NumObservations_ = sum(~excl);
end

%% INNER FUNCTIONS

function [w,predLocs] = getDataVariableFS(w,~,X,predLocs,respLoc,vtype)
if isempty(w)
    return
end
if isa(X,'dataset') 
    X = dataset2table(X);
end
if isa(X,'table') && internal.stats.isString(w) % a dataset variable name
    [tf,wloc] = ismember(w,X.Properties.VariableNames);
    if ~tf
        error(message('stats:classreg:regr:FitObject:BadVariableName', vtype, w));
    end
    w = X.(w);
    predLocs = setdiff(predLocs,wloc);
    if wloc == respLoc
        respLoc = max(predLocs);
        predLocs = setdiff(predLocs,respLoc);
    end
end
if ~(isnumeric(w) || islogical(w)) || ~isvector(w) || ~isreal(w)
    error(message('stats:classreg:regr:FitObject:BadVariableValues', vtype));
end
w = w(:); % force to a column
end

function range = getVarRangeFS(v,asCat,excl)
v(excl,:) = [];
if asCat % create a list of the unique values
    if isa(v,'categorical')
        % For categorical classes, get the values actually present in the
        % data, not the set of possible values.
        range = unique(v(:));
        range = range(~isundefined(range)); % empty if NaNs
            
    else
        % For classes other than categorical, the list of unique values is
        % also the list of possible values.  But get it in the same order as
        % what grp2idx defines for each class.
        [~,~,range] = grp2idx(v); % leaves NaN or '' out of glevels
    end
    if ~ischar(range)
        range = range(:)'; % force a row
    end
elseif isnumeric(v) || islogical(v) % find min and max
    range = [min(v,[],1)  max(v,[],1)]; % ignores NaNs unless all NaNs
else
    range = NaN(1,2);
end
end