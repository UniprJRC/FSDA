function [varNames,predictorVars,responseVar] = ...
    getVarNamesFS(varNames,predictorVars,responseVar,nx)
if isempty(varNames)
    % Create varNames from predictorVars and responseVar, supplied or
    % default
    if ~isempty(predictorVars) && (iscell(predictorVars) || ischar(predictorVars))
        if iscell(predictorVars)
            predictorVars = predictorVars(:);
        else
            predictorVars = cellstr(predictorVars);
        end
        if length(predictorVars)~=nx || ...
                ~internal.stats.isStrings(predictorVars,true)
            error(message('stats:classreg:regr:FitObject:BadPredNames'));
        end
        pnames = predictorVars;
    else
        pnames = internal.stats.numberedNames('x',1:nx)';
        if ~isempty(responseVar) && internal.stats.isString(responseVar)
            pnames = genvarname(pnames,responseVar);
        end
        if isempty(predictorVars)
            predictorVars = pnames;
        else
            predictorVars = pnames(predictorVars);
        end
    end
    if isempty(responseVar)
        responseVar = genvarname('y',predictorVars);
    end
    varNames = [pnames; {responseVar}];
else
    if ~internal.stats.isStrings(varNames,true)
        error(message('stats:classreg:regr:FitObject:BadVarNames'));
    end
    
    % If varNames is given, figure out the others or make sure they are
    % consistent with one another
    if isempty(responseVar) && isempty(predictorVars)
        % Default is to use the last name for the response
        responseVar = varNames{end};
        predictorVars = varNames(1:end-1);
        return
    end
    
    % Response var must be a name
    if ~isempty(responseVar)
        [tf,rname] = internal.stats.isString(responseVar,true);
        if tf
            responseVar = rname;
            if ~ismember(responseVar,varNames)
                error(message('stats:classreg:regr:FitObject:MissingResponse'));
            end
        else
            error(message('stats:classreg:regr:FitObject:BadResponseVar'))
        end
    end
    
    % Predictor vars must be names or an index
    if ~isempty(predictorVars)
        [tf,pcell] = internal.stats.isStrings(predictorVars);
        if tf
            predictorVars = pcell;
            if ~all(ismember(predictorVars,varNames))
                error(message('stats:classreg:regr:FitObject:InconsistentNames'));
            end
        elseif isValidIndexVector(varNames,predictorVars)
            predictorVars = varNames(predictorVars);
        else
            error(message('stats:classreg:regr:FitObject:InconsistentNames'))
        end
    end
    
    % One may still be empty
    if isempty(predictorVars)
        predictorVars = setdiff(varNames,{responseVar});
    elseif isempty(responseVar)
        % If predictorVar is given, there should be just the response left
        responseVar = setdiff(varNames,predictorVars);
        if isscalar(responseVar)
            responseVar = responseVar{1};
        else
            error(message('stats:classreg:regr:FitObject:AmbiguousResponse'));
        end
    else
        if ~ismember({responseVar},varNames) || ...
                ~all(ismember(predictorVars,varNames)) || ...
                ismember({responseVar},predictorVars)
            error(message('stats:classreg:regr:FitObject:InconsistentNames'))
        end
    end
end
end