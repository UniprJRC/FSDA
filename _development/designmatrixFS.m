function [doubleX,terms,cols2vars,cols2terms,colnames,termnames,doubley,namey,catnums] = designmatrixFS(A,varargin)
%DESIGNMATRIX Construct regression design matrix from dataset/table array.
%   doubleX = DESIGNMATRIX(A) returns a design matrix doubleX created from variables in
%   the NOBS-by-NVARS dataset/table array A.  doubleX is a design matrix for a linear
%   regression model that includes an intercept and main effects for each of
%   the predictor variables in A.  DESIGNMATRIX treats the first through
%   next-to-last variables in A as predictor variables, and creates columns in
%   X corresponding to those variables.  DESIGNMATRIX treates the last
%   variable in A as a response variable, and does not create any
%   corresponding columns in X.
%
%   DESIGNMATRIX copies numeric variables in A directly to columns of X, and
%   converts variables in A that are categorical, logical, cell arrays of
%   strings, or char into one or more binary dummy variables and copies those
%   to columns of X.  Only numeric variables in A may have more than one
%   column, and all variables in A must have two dimensions.
%
%   X is NOBS-by-NCOLS, where NCOLS is the sum of the numbers of columns in
%   the numeric predictor variables in A, plus the number of dummy variables
%   created from the remaining predictor variables.  X is a double matrix
%   unless any of the predictor variables in A are single.
%
%   [X,TERMS] = DESIGNMATRIX(A) returns an NTERMS-by-NVARS model terms matrix
%   TERMS for the regression model.  Rows of TERMS correspond to terms in the
%   model, columns to variables in the model.  TERMS(I,J) contains the power
%   of the J-th variable in the I-th term.  A row of zeros corresponds to
%   the intercept term.
%
%   [X,TERMS,COLS2VARS] = DESIGNMATRIX(A) returns an NVARS-by-NCOLS matrix
%   VARS2COLS, whose (I,J)th element contains a 1 if the J-th column of X
%   depends on the I-th variable in A, and 0 otherwise.
%
%   [X,TERMS,COLS2VARS,COLS2TERMS,COLNAMES] = DESIGNMATRIX(A) returns a vector
%   COLS2TERMS with NCOLS elements specifying the term number for each column
%   of X, and a cell array COLNAMES containing the names of the columns in X.
%
%   [...] = DESIGNMATRIX(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control how DESIGNMATRIX creates
%   columns in X from variables in A.  Parameters are:
%
%      'Model'             Specifies the regression model.  The value may be
%                          a string, or a terms matrix as described above.
%                          Valid strings are 'constant', 'linear', 'interactions',
%                          'purequadratic', or 'quadratic', or a string of the
%                          form 'polyMN...', where [M N ...] is a sequence of P
%                          numeric digits representing a P-variate polynomial
%                          in which the model contains main and interaction
%                          terms for the 1st through M-th powers of the first
%                          variable, 1st though N-th powers of the second
%                          variable, etc.
%
%      'PredictorVars'     Specifies which variables in A to treat as predictors.
%                          Only these variables are copied/converted to columns
%                          in X. Default is 2:NVARS.
%
%      'ResponseVar'       Specifies which variable in A to treat as the
%                          response.  This variable is not copied/converted to
%                          columns in X.  Default is last column.
%
%      'CategoricalVars'   Specifies which variables in A to treat as categorical
%                          variables.  These variables are converted to one or
%                          more binary dummy variables before copying to
%                          columns in X, and linear terms are created for
%                          these variables.  Default is false for numeric
%                          variables, true for categorical, logical, cell
%                          arrays of strings or char.
%
%      'CategoricalLevels' A cell array with NCOLS elements. Each entry is
%                          typically the VariableInfo.Range value for that
%                          column, where VariableInfo is a FitObject
%                          property. Only the entries corresponding to
%                          categorical variables are used. Each such entry
%                          is a cell or categorical array listing all
%                          unique values of that categorical variable.
%
%      'Intercept'         A logical value indicating whether or not to include
%                          an intercept term in the model.  Default is true.
%
%      'DummyVarCoding'    A string specifying the coding to use for dummy
%                          variables created from categorical variables, or a
%                          cell array with NVARS elements containing a string
%                          for each categroical variable.  Valid coding schemes
%                          are 'reference' (coefficient for first category set
%                          to zero),'referencelast' (coefficient for last
%                          category set to zero), 'effects' (coefficients sum
%                          to zero with last as reference), 'effectsfirst'
%                          (coefficients sum to zero with first as
%                          reference), 'full' (one dummy variable for each 
%                          category), 'difference' (contrast from previous
%                          category), and 'backwardDifference' (contrast from
%                          next category). and 'ordinal' (+1 for this and
%                          later categories, -1 for earlier categories).
%                          Default is 'reference'. 
%
%   See also DATASET, PREDICTORMATRIX, DATASET/DOUBLE, DATASET/SINGLE.

%   Copyright 2010-2016 The MathWorks, Inc.

%% Beginning of code 

if isa(A,'dataset') 
    A=dataset2table(A);
end

[nobs,nvars] = size(A);

if isa(A,'table')  
    isDataset = true;
    % Treat categoricals, strings, and logicals as categorical by default.
    catVars = varfun(@(x) isa(x,'categorical') || iscellstr(x) || ischar(x) || islogical(x), A,'OutputFormat','uniform');
elseif (isfloat(A) || islogical(A)) && ismatrix(A)
    isDataset = false;
    % Treat everything as numeric by default.
    catVars = false(1,size(A,2));
elseif isa(A,'categorical')
    isDataset = false;
    catVars = true(1,size(A,2));
else
    error(message('stats:classreg:regr:modelutils:BadDataType'));
end

paramNames = { 'Model' 'VarNames' 'PredictorVars' 'ResponseVar' 'Intercept' ...
               'CategoricalVars' 'CategoricalLevels' 'DummyVarCoding'};
paramDflts = {'linear'        []     1:(nvars-1)         nvars        true  ...
                        catVars                  []  {'reference'}};
[model,varNames,predictorVars,responseVar,includeIntercept, ...
 treatAsCategorical,catLevels,dummyCoding,supplied] = ...
    internal.stats.parseArgs(paramNames, paramDflts, varargin{:});

if isDataset
    varNames = A.Properties.VariableNames;
elseif ~supplied.VarNames
    varNames = strcat({'X'},num2str((1:nvars)','%-d'))';
end

% Default is to treat everything except the response as predictors, unless
% there's a terms matrix.
if supplied.PredictorVars
    if isnumeric(predictorVars)
        if any(predictorVars < 1) || any(predictorVars > nvars)
            error(message('stats:classreg:regr:modelutils:BadPredVarsNumeric'));
        end
    elseif isDataset && internal.stats.isStrings(predictorVars)
        [~,predVarInds] = ismember(predictorVars,varNames);
        if any(predVarInds == 0)
            error(message('stats:classreg:regr:modelutils:BadPredVarsName',predictorVars(predVarInds(1))));
        end
        predictorVars = predVarInds;
    else
        if isDataset
            error(message('stats:classreg:regr:modelutils:BadPredVarsDataset'));
        else
            error(message('stats:classreg:regr:modelutils:BadPredVarsMatrix'));
        end
    end
    
    if ~supplied.ResponseVar
        % Treat the one variable not in the predictors as the response.
        responseVar = setdiff(1:nvars,predictorVars);
        if numel(responseVar) > 1 % none is OK
            error(message('stats:classreg:regr:modelutils:CannotDetermineResponse'));
        end
    end
end

% Default is to treat the first var as the response, unless there's a terms
% matrix.
if supplied.ResponseVar
    if isempty(responseVar)
        % OK
    elseif isscalar(responseVar) && isnumeric(responseVar)
        if (responseVar < 1) || (responseVar > nvars)
            error(message('stats:classreg:regr:modelutils:BadResponseNumeric'));
        end
    elseif isDataset && internal.stats.isString(responseVar)
        responseVarIdx = find(strcmp(responseVar,varNames));
        if isempty(responseVarIdx)
            error(message('stats:classreg:regr:modelutils:BadResponseName',responseVar));
        end
        responseVar = responseVarIdx;
    else
        if isDataset
            error(message('stats:classreg:regr:modelutils:BadResponseDataset'));
        else
            error(message('stats:classreg:regr:modelutils:BadResponseMatrix'));
        end
    end
    
    if ~supplied.PredictorVars
        % Treat everything except the response as predictors.
        predictorVars = setdiff(1:nvars,responseVar);
    end
end

% Default is to include an intercept as the first term, unless there's a terms
% matrix.
if ~isscalar(includeIntercept) || ~islogical(includeIntercept)
    error(message('stats:classreg:regr:modelutils:BadIntercept'));
end

% Default is to treat categoricals and cellstrs as categorical.
if ~islogical(treatAsCategorical) || numel(treatAsCategorical)~=nvars
    error(message('stats:classreg:regr:modelutils:BadCategorical'));
end

% Default is a linear model.
if internal.stats.isString(model)
    whichVars = false(1,nvars); whichVars(predictorVars) = true;
    terms = classreg.regr.modelutils.model2terms(model,whichVars,includeIntercept,treatAsCategorical);
    
elseif isnumeric(model) && ismatrix(model)
    if size(model,2) ~= nvars
        error(message('stats:classreg:regr:modelutils:BadModelMatrix'));
    end
    terms = model;
    interceptRow = find(sum(terms,2) == 0);
    if supplied.Intercept && includeIntercept && isempty(interceptRow)
        % Add an intercept if requested, even if the terms matrix didn't have one
        terms = [zeros(1,nvars); terms];
    elseif ~includeIntercept && ~isempty(interceptRow)
        error(message('stats:classreg:regr:modelutils:InterceptAmbiguous'));
    end
    predVarsIn = predictorVars;
    predictorVars = find(sum(terms,1) > 0);
    if supplied.PredictorVars && ~isempty(setxor(predictorVars,predVarsIn))
        error(message('stats:classreg:regr:modelutils:PredictorsAmbiguous'));
    end
    if ~supplied.ResponseVar
        responseVar = setdiff(1:nvars,predictorVars);
    end
    
else
    error(message('stats:classreg:regr:modelutils:ModelStringOrMatrix'));
end

if isDataset
    doubley=table2array(A(:,varNames(responseVar)));
    namey=varNames(responseVar);
else
      doubley=A(:,responseVar);
      namey='y';
end

% Default is to use reference (first) coding
if supplied.DummyVarCoding
    codingStrings = {'full' 'reference' 'referencelast' 'effects' 'effectsfirst' 'difference' 'backwardDifference' 'ordinal'};
    [haveStrs,dummyCoding] = internal.stats.isStrings(dummyCoding);
    if ~haveStrs
        error(message('stats:classreg:regr:modelutils:BadDummyCode'));
    else
        for i = 1:length(dummyCoding)
            if ~isempty(dummyCoding{i}) && ~any(strcmpi(dummyCoding{i},codingStrings))
                error(message('stats:classreg:regr:modelutils:UnrecognizedDummyCode',internal.stats.listStrings(codingStrings)));
            end
        end
    end
end
if isscalar(dummyCoding)
    dummyCoding = repmat(dummyCoding,1,nvars);
    dummyCoding(~treatAsCategorical) = {''};
end

if nargout>=9 && ~all(cellfun('isempty',dummyCoding) | strcmp('full',dummyCoding))
    error(message('stats:classreg:regr:modelutils:FullCodingRequired'));
end

if isDataset
    varClasses = varfun(@class,A,'InputVariables',predictorVars,'OutputFormat','cell');
    if any(strcmp('single',varClasses))
        outClass = 'single';
    else
        outClass = 'double';
    end
elseif isnumeric(A)
    outClass = class(A);
else % islogical(data)
    outClass = 'double';
end

% Take a stab at preallocating, assume 1 column for each continuous var and 4
% columns for each categorical.  But X can grow if we underestimate, and is
% truncated if we overestimate.
ncolsX = length(predictorVars) + 3*sum(treatAsCategorical(predictorVars));
doubleX = zeros(nobs,ncolsX,outClass);
colnames = cell(1,ncolsX);
cols2vars = false(nvars,ncolsX);
ncols = 0;

% Create main terms, including dummy vars.
nvarlevels = ones(1,nvars);
for j = predictorVars
    vnamej = varNames{j};
    if isDataset
        Xj = A.(vnamej);
    else
        Xj = A(:,j);
    end
    ncolsj = size(Xj,2);
    if ~ismatrix(Xj) || ncolsj==0
        error(message('stats:classreg:regr:modelutils:NotMatrices'));
    end
    
    % Convert categorical vars to dummy vars.
    if treatAsCategorical(j)
        if supplied.CategoricalLevels
            [Xj,dummynames] = dummyVars(dummyCoding{j},Xj,catLevels{j});
        else
            [Xj,dummynames] = dummyVars(dummyCoding{j},Xj);
        end
        ncolsj = size(Xj,2);
        colnamesj = strcat(vnamej,'_',dummynames);
        nvarlevels(j) = ncolsj;
        
    % Include continuous vars as is.
    elseif isnumeric(Xj) || islogical(Xj)
        colnamesj = {vnamej};
        % Expand labels for a matrix variable
        if ncolsj > 1
            cols = cellstr(num2str((1:ncolsj)','%-d'))';
            colnamesj = strcat(vnamej,'_', cols);
        end
    else
        if catVars(j)
            error(message('stats:classreg:regr:modelutils:CatNotContinuous'));
        else
            error(message('stats:classreg:regr:modelutils:BadVariableType'));
        end
    end
    
    % Add columns to X if necessary.
    if ncols+ncolsj > ncolsX
        predVarsLeft = predictorVars(predictorVars > j);
        ncolsX = ncols + ncolsj + length(predVarsLeft) + 3*sum(treatAsCategorical(predVarsLeft));
        doubleX(nobs,ncolsX) = 0;
        colnames{ncolsX} = '';
        cols2vars(nvars,ncolsX) = 0;
    end
    
    colsj = (ncols+1):(ncols+ncolsj);
    doubleX(:,colsj) = Xj;
    colnames(:,colsj) = colnamesj;
    ncols = ncols + ncolsj;
    
    % Keep track of which columns of X correspond to which vars in data.
    cols2vars(j,colsj) = true;
end

Xmain = doubleX(:,1:ncols);
colnamesMain = colnames(:,1:ncols);
cols2varsMain = cols2vars(:,1:ncols);

nterms = size(terms,1);
doubleX = zeros(nobs,0,outClass);
colnames = {};
cols2vars = zeros(nvars,0);
cols2terms = zeros(1,nterms);
termnames = {};
if nargout>=9
    catnums = cell(1,nterms);
end
for j = 1:nterms
    varsj = find(terms(j,:));
    if isempty(varsj)
        Xj = ones(nobs,1,outClass);
        colnames = [colnames '(Intercept)'];
        cols2vars = [cols2vars false(nvars,1)];
        termnames = [termnames '(Intercept)'];
    else
        colsj = cols2varsMain(varsj(1),:);
        Xj = Xmain(:,colsj);
        colnamesj = colnamesMain(colsj);
        termnamesj = varNames{varsj(1)};
        expon = terms(j,varsj(1));
        if expon > 1
            if size(Xj,2) == 1
                Xj = Xj.^expon;
                colnamesj = strcat(colnamesj,'^',num2str(expon));
                termnamesj = strcat(termnamesj,'^',num2str(expon));
            else
                Xj1 = Xj;
                colnamesj1 = colnamesj;
                for e = 2:expon
                    [rep1,rep2] = allpairs2(1:size(Xj1,2),1:size(Xj,2));
                    Xj = Xj1(:,rep1) .* Xj(:,rep2);
                    colnamesj = strcat(colnamesj1(rep1),'_',colnamesj(rep2));
                end
            end
        end
        for k = 2:length(varsj)
            colsjk = cols2varsMain(varsj(k),:);
            Xjk = Xmain(:,colsjk);
            colnamesjk = colnamesMain(colsjk);
            termnamesjk = varNames{varsj(k)};
            expon = terms(j,varsj(k));
            if expon > 1
                if size(Xjk,2) == 1
                    Xjk = Xjk.^expon;
                    colnamesjk = strcat(colnamesjk,'^',num2str(expon));
                    termnamesjk = strcat(termnamesjk,'^',num2str(expon));
                else
                    Xjk1 = Xjk;
                    colnamesjk1 = colnamesjk;
                    for e = 2:expon
                        [rep1,rep2] = allpairs2(1:size(Xjk1,2),1:size(Xjk,2));
                        Xjk = Xjk1(:,rep1) .* Xjk(:,rep2);
                        colnamesjk = strcat(colnamesjk1(rep1),'_',colnamesjk(rep2));
                    end
                end
            end
            [rep1,rep2] = allpairs2(1:size(Xj,2),1:size(Xjk,2));
            Xj = Xj(:,rep1) .* Xjk(:,rep2);
            colnamesj = strcat(colnamesj(rep1),':',colnamesjk(rep2));
            termnamesj = [termnamesj, ':', termnamesjk];
        end
        colnames = [colnames colnamesj];
        termnames = [termnames termnamesj];
    
        % Keep track of which columns of X correspond to which vars in data.
        cols2varsj = false(nvars,size(Xj,2));
        cols2varsj(varsj,:) = true;
        cols2vars = [cols2vars cols2varsj];
    end
    cols2terms(size(doubleX,2) + (1:size(Xj,2))) = j;
    doubleX = [doubleX Xj];
    if nargout>=9 && any(treatAsCategorical(varsj))
        termlevels = nvarlevels(varsj);
        termlevels = termlevels(treatAsCategorical(varsj));
        catnums{j} = fullfact(termlevels);
    end
end


%-----------------------------------------------------------------------------
function [rep1,rep2] = allpairs2(i,j)
[rep1,rep2] = ndgrid(i,j);
rep1 = rep1(:)';
rep2 = rep2(:)';


%-----------------------------------------------------------------------------
function [X,colnames] = dummyVars(method,group,glevels)

if nargin < 3 || isempty(glevels)
    [gidx,gn] = grp2idx(group);        
    if isa(group,'categorical')
        uidx = unique(gidx);
        if length(uidx) < length(gn)
            warning(message('stats:classreg:regr:modelutils:LevelsNotPresent'));
            [gidx,gn] = grp2idx(removecats(group));
        end
    end
else % use the specified set of possible values
    [~,gn,glevels] = grp2idx(glevels); % let grp2idx order them
    group = convertVar(group,glevels); % convert group to be compatible with glevels
    if ischar(group)
        [tf,gidx] = ismember(group, glevels,'rows');
    else
        [tf,gidx] = ismember(group,glevels);
    end
    gidx(tf==0) = NaN;
end
ng = length(gn);

switch lower(method)
case 'full' % all categories
    X0 = eye(ng);
    colnames = gn;
case 'reference' % reference is first category
    X0 = eye(ng); X0(:,1) = [];
    colnames = gn(2:end);
case 'referencelast' % reference is last category
    X0 = eye(ng,ng-1);
    colnames = gn(1:end-1);
case 'effects' % sum to zero, -1 for last
    X0 = eye(ng,ng-1);
    X0(end,:) = -1;
    colnames = gn(1:end-1);
case 'effectsfirst' % sum to zero, -1 for first
    X0 = [-ones(1,ng-1); eye(ng-1)];
    colnames = gn(2:end);
case 'difference' % relative to previous
    X0 = tril(ones(ng,ng)); X0(:,1) = [];
    colnames = strcat(gn(2:end),'_increment');    
case 'backwarddifference' % relative to following
    X0 = triu(ones(ng,ng)); X0(:,end) = [];
    colnames = strcat(gn(1:end-1),'_decrement'); 
case 'ordinal' 
    X0 = tril(ones(ng))-triu(ones(ng),1);
    X0(:,1) = [];
    colnames = gn(2:end); 
otherwise
    error(message('stats:classreg:regr:modelutils:UnknownCodingMethod',method));
end

k = isnan(gidx);
if any(k)
    X(k,1:size(X0,2)) = NaN;
    k = ~k;
    X(k,:) = X0(gidx(k),:);
else
    X = X0(gidx,:);
end


%-----------------------------------------------------------------------------
function a = convertVar(a,b)

if isa(a,'categorical')
    if ischar(b)|| iscell(b)
        a = cellstr(a);
    elseif isa(b,'categorical')
        [tf,loc] = ismember(a,b);
        a = matlab.internal.datatypes.defaultarrayLike(size(a),'Like',b);
        a(tf) = b(loc(tf));
    else
        error(message('stats:classreg:regr:modelutils:CannotConvert',class(a),class(b)));
    end
elseif isnumeric(a) || islogical(a)
    if isnumeric(b) || islogical(b)
        % OK
    else
        error(message('stats:classreg:regr:modelutils:CannotConvert',class(a),class(b)));
    end
elseif iscell(a) || ischar(a) || isstring(a)
    if ischar(b) || iscell(b) || isstring(b)
        % OK
    elseif isa(b,'categorical')
        [tf,loc] = ismember(cellstr(a),b);
        if iscell(a)
            a = matlab.internal.datatypes.defaultarrayLike(size(a),'Like',b);
        else
            a = matlab.internal.datatypes.defaultarrayLike([size(a,1),1],'Like',b);
        end
        a(tf) = b(loc(tf));
    else
        error(message('stats:classreg:regr:modelutils:CannotConvert',class(a),class(b)));
    end
else
    error(message('stats:classreg:regr:modelutils:BadPredictorType'));
end