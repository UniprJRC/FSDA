function out = FSreg_fit(X,varargin)
% Not intended to be called directly. Use FSreg instead.

[X,y,haveDataset,otherArgs] = FSreg_HandleDataArgs(X,varargin{:});

% VarNames are optional names for the X matrix and y vector.  A
% dataset/table defines its own list of names, so this is not accepted
% with a dataset/table.

% PredictorVars is an optional list of the subset of variables to
% actually use as predictors in the model, and is only needed with
% an alias.  A terms matrix or a formula string already defines
% which variables to use without that.  ResponseVar is an optional
% name that is not needed with a formula string.

% rankwarn is undocumented and is used during stepwise fitting

paramNames = {'Intercept' 'PredictorVars' 'ResponseVar' ...
    'Weights' 'Exclude' 'CategoricalVars' 'VarNames'...
    'DummyVarCoding' 'rankwarn' ...
    'Monitoring' 'Control' 'Family' 'Estimator'};

% The if below is just to understand what is modelDef
if isempty(otherArgs)
    % if otheerArgs is empty the user has not specified model string and
    % therefore the default model is main effects only.
    modelDef = 'linear';
else
    arg1 = otherArgs{1};
    if mod(length(otherArgs),2)==1
        % odd, model followed by pairs.
        % In this case the user has called
        % the procedure specifying model string ('linear' 'quadratic') and a
        % series of name/value pairs)
        modelDef = arg1;
        otherArgs(1) = [];
    elseif ischar(arg1) && ...
            any(strncmpi(arg1,paramNames,length(arg1)))
        % even, omitted model but included name/value pairs
        % In this case the user has called
        % the procedure without specfying model string ('linear' 'quadratic') with a
        % series of name/value pairs)
        
        
        % DOME REMARK: nella formulazione originale invece di essere ischar(arg1)
        % era internal.stats.isString(arg1)
        
        modelDef = 'linear';
    end
end

% The following lines are not necessary anymore
% paramDflts = {[] [] [] [] [] [] [] [] 'reference' true};
% [intercept,predictorVars,ResponseVar,weights,exclude, ...
%     asCatVar,varNames,robustOpts,dummyCoding,rankwarn,supplied] = ...
%     internal.stats.parseArgs(paramNames, paramDflts, otherArgs{:});


CategoricalVars='';
DummyVarCoding='referencelast';
ResponseVar='';
intercept=true;
% Remark: 'VarNames' is not applicable to variables in a table or dataset
% array, because those variables already have names.
VarNames='';

Family='homo';
Estimator='FS';
Monitoring=false;
Control='';

UserOptions=otherArgs(1:2:length(otherArgs));
if ~isempty(UserOptions)
    
    options=struct('ResponseVar',ResponseVar,'CategoricalVars',CategoricalVars, ...
        'DummyVarCoding',DummyVarCoding,'VarNames',VarNames,'intercept',intercept,...
        'Family',Family,'Estimator',Estimator,'Monitoring',Monitoring,'Control',Control);
    
    
    % Check if number of supplied options is valid
    if length(otherArgs) ~= 2*length(UserOptions)
        error('FSDA:FSreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:FSreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    for i=1:2:length(otherArgs)
        options.(otherArgs{i})=otherArgs{i+1};
    end
    
    ResponseVar=options.ResponseVar;
    CategoricalVars=options.CategoricalVars;
    DummyVarCoding=options.DummyVarCoding;
    intercept=options.intercept;
    VarNames=options.VarNames;
    Family=options.Family;
    Estimator=options.Estimator;
    Control=options.Control;
    Monitoring=options.Monitoring;
end

if haveDataset == true
    VarNames=X.Properties.VariableNames;
else
    % If we are not dealing with a dataset (table) then it is necessary to
    % create variable names for response and predictors
    
    if isempty(VarNames)
        
        nx=size(X,2);
        
        PredictorVars='';
        [VarNames,~,ResponseVar] = getVarNamesFS(VarNames,PredictorVars,ResponseVar,nx);
        % [VarNames,PredictorVars,ResponseVar] = getVarNamesFS(VarNames,PredictorVars,ResponseVar,nx);
        
    end
    X=array2table([X y]);
    X.Properties.VariableNames=VarNames;
    
end


% Now we call routins LinearFormula which is inside
% [ matlabroot \toolbox\stats\classreg\+classreg\+regr]  in order to create
% formula object
% In order to see the MATLAB help for this formula press F1 key after
% selecting the string: classreg.regr.LinearFormula

% VarNames is a cell array o strings containing the names of all variables present in the table
formula = classreg.regr.LinearFormula(modelDef,VarNames,'',intercept,'identity');
% formula which is used
disp(formula)
% The most important properties of formula are
% InModel           % which vars are in model as a predictor? (logical array)
% LinearPredictor   % the thing on the RHS of the '~'
% PredictorNames    % names of predictors (regressors) actually in model
% Terms             % NTerms-by-NVars model terms matrix TERMS for the regression model
%                   Rows of Terms correspond to terms in the model, columns
%                   to variables in the model.  Terms(i,j) contains the
%                   power of the i-th variable in the i-th term.  A row of
%                   zeros corresponds to the intercept term.
% For example if one supplies a quadratic model in Weight and Model_Year
%       and Weight and Model_Year are the first two variables of a set of 3 regressors,
%       Terms matrix is as follows
%                      0     0     0
%                      1     0     0
%                      0     1     0
%                      1     1     0
%                      2     0     0
%                      0     2     0
%
% TermNames         % names of the rows of Terms
%                    For the example given above TermNames is as follows
%                     '(Intercept)'
%                     'Weight'
%                     'Model_Year'
%                     'Weight:Model_Year'
%                     'Weight^2'
%                     'Model_Year^2'
%
% NTerms            % number of rows in Terms
% NVars             % number of columns in Terms (in the example above
%                     NVars  is 3)
% NPredictors       % number of predictors actually in model (in the
%                     example above NPredictors is 2)
% HasIntercept      % logical, true if the model contains the intercept term
% disp              % what is shown when one types disp(formula)
% char              % string associated with formula which is used
% VariableNames     % all the variable names in the table (dataset)
if isempty(ResponseVar)
    ResponseVar=formula.ResponseName;
end

weights='';
asCatVar='';
exclude='';

% If CategoricalVars is empty than there is an automatic detection of
% categorical variables
if isempty(CategoricalVars)
    mode1=assignDataFS(X,y,weights,asCatVar,formula.VariableNames,exclude);
    CategoricalVars=mode1.VariableInfo.IsCategorical;
end


% Given matrix formula.Terms procedure below computes doubleX
[doubleX,terms,cols2vars,cols2terms,NamesdoubleX,termnames,doubley,namey]=designmatrixFS(X,'Model',formula.Terms,...
    'ResponseVar',ResponseVar,...
    'CategoricalVars',CategoricalVars,...
    'DummyVarCoding',DummyVarCoding,'Intercept',true);


% define X
doubleX=doubleX(:,2:end);

% DAFARE
% Formula = classreg.regr.LinearFormula(modelDef,varNames,'ResponseVar',interc,'identity');

%% Extract the estimator the user has chosen
%
% if strcmp(estimator,'hetero')
% Now we need to transform structure control varargin cell
if isstruct(Control)
    
    fcontrol=  fieldnames(Control);
    a=[fcontrol struct2cell(Control)]';
    ControlNamePairs=a(:)';
    %   ControlNamePairs=ControlNamePairs{:};
    
else
    ControlNamePairs={};
end

disp(['Response is' namey])

if strcmp(Estimator,'FS')
    
    if strcmp(Family,'homo')
        if Monitoring == true
            out=FSReda(doubley,doubleX,ControlNamePairs{:});
        else
            out=FSR(doubley,doubleX,ControlNamePairs{:});
        end
        
        
    elseif strcmp(Family,'hetero')
        if Monitoring == true
            out=FSRH(doubley,doubleX,Z,ControlNamePairs{:});
        else
            out=FSRHeda(doubley,doubleX,Z,ControlNamePairs{:});
        end
        
    elseif strcmp(Family,'Bayes')
        out=FSRB(doubley,doubleX,Z,ControlNamePairs{:});
    end
    
elseif   strcmp(Estimator,'S')
    if Monitoring == true
        out=Sregeda(doubley,doubleX,ControlNamePairs{:});
    else
        out=Sreg(doubley,doubleX,ControlNamePairs{:});
    end
    
elseif  strcmp(Estimator,'MM')
    if Monitoring == true
        out=MMregeda(doubley,doubleX,ControlNamePairs{:});
    else
        out=MMreg(doubley,doubleX,ControlNamePairs{:});
    end
    
elseif strcmp(Estimator,'LXS')
    if Monitoring == true
        out=LXSeda(doubley,doubleX,ControlNamePairs{:});
    else
        out=LXS(doubley,doubleX,ControlNamePairs{:});
        
    end
    
else
    error('FSDA:FSreg_fit:WrongEstimator','Estimator not recognized')
end

out.NamesdoubleX=NamesdoubleX;
out.namey=namey;
out.doubleX=doubleX;

end