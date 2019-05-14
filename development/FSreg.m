function out = FSreg(TBL,varargin) % [y, X | TBL], modelDef, ...
%FSreg creates linear (Bayesian or heteroskedastic) regression model by fitting to data.
%
%<a href="matlab: docsearchFS('FSreg')">Link to the help function</a>
%
%  Required input arguments:
%
%
%   TBL  : Input data, specified as a table or dataset array or a matrix. 
%          Table, dataset or numeric matrix.
%          If TBL is a table or a dataset it contains the variables to be
%          used inside FSreg. Alternative if TBL is a numeric matrix it
%          refers to the matrix of explanatory variables. In this case the
%          second argument must be y (the vector containing the response).
%          In other words, FSRreg can be called using FSreg(TBL) or
%          FSreg(y,X).
%          If FSRreg is called as FSreg(y,X,'modelspec') or
%          FSRreg(TBL,'modelspec') modelspec must be a character vector
%          representing a formula in the form 'Y ~ terms', where the terms
%          are in Wilkinson Notation, or a  t-by-(p + 1) terms matrix. For
%          more information about modelspec please see the help of fitlm
%
%  Optional input arguments:
%
%   Intercept :  Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               else no constant term will be included.
%               Example - 'intercept',1 
%               Data Types - double
%
% CategoricalVars : Categorical variables. 
%            Cell array of character vectors | logical or numeric index
%            vector. Categorical variables in the fit, specified as the
%            comma-separated pair consisting of 'CategoricalVars' and
%            either a cell array of character vectors of the names of the
%            categorical variables in the table or dataset array tbl, or a
%            logical or numeric index vector indicating which columns are
%            categorical. If data is in a table or dataset array tbl, then
%            the default is to treat all categorical or logical variables,
%            character arrays, or cell arrays of character vectors as
%            categorical variables. If data is in matrix X, then the
%            default value of this name-value pair argument is an empty
%            matrix []. That is, no variable is categorical unless you
%            specify it.
%            For example, you can specify the observations 2 and 3 out of 6
%            as categorical using either 'CategoricalVars',[2,3] or
%            'CategoricalVars',logical([0 1 1 0 0 0])
%               Example - 'CategoricalVars',[2,3]
%               Data Types - single | double | logical
%
%PredictorVars : Predictor variables.
%                Cell array of character vectors | logical or numeric
%                index vector.
%               Predictor variables to use in the fit, specified as the
%               comma-separated pair consisting of 'PredictorVars' and
%               either a cell array of character vectors of the variable
%               names in the table or dataset array tbl, or a logical or
%               numeric index vector indicating which columns are predictor
%               variables. The character vectors should be among the names
%               in tbl, or the names you specify using the 'VarNames'
%               name-value pair argument. The default is all variables in
%               X, or all variables in tbl except for ResponseVar.
%               For example, you can specify the second and third variables
%               as the predictor variables using 'PredictorVars',[2,3] or 
%               'PredictorVars',logical([0 1 1 0 0 0])
%               Example - 'CategoricalVars',[1,3]
%               Data Types: single | double | logical | cell
%
% ResponseVar : Response variable. 
%               Last column in tbl (default) | character vector containing
%               variable name | logical or numeric index vector.
%               Response variable to use in the fit, specified as the
%               comma-separated pair consisting of 'ResponseVar' and either
%               a character vector containing the variable name in the
%               table or dataset array tbl, or a logical or numeric index
%               vector indicating which column is the response variable.
%               You typically need to use 'ResponseVar' when fitting a
%               table or dataset array tbl.
%               For example, you can specify the fourth variable, say
%               yield, as the response out of six variables, in one of the
%               following ways.
%               'ResponseVar','yield', 'ResponseVar',[4], 'ResponseVar',logical([0 0 0 1 0 0])
%               Example - 'CategoricalVars',[3]
%               Data Types: single | double | logical | char
%   VarNames : Names of variables in fit. 
%              {'x1','x2',...,'xn','y'} (default) | cell array of character
%              vectors. Names of variables in fit, specified as the
%              comma-separated pair consisting of 'VarNames' and a cell
%              array of character vectors including the names for the
%              columns of X first, and the name for the response variable y
%              last. Note that 'VarNames' is not applicable to variables in
%              a table or dataset array, because those variables already
%              have names.
%               Example - 'VarNames',{'Horsepower','Acceleration','Model_Year','MPG'}
%               Data Types: cell
% 
%  Estimator : string which defines the estimator or structure containing
%              options for the estimator. Character. Estimator is a string
%              which may contain one of the following values 'FS', 'LXS',
%              'S' or 'MM'.
%  Family  :   regression context. 
%              Character specifying the context of the regression model or
%              structure specifying the estimator.
%              If Family is a character possible values of Family are
%              'homoskedastic' ('homo') this is the default choice,
%              'heteroskedastic' ('hetero') or 'Bayesian' (bayes).
%              If family is a structure it is possible to specify in the
%              field estimator the estimator which is used.
%               Family specifies which function of FSDA is used.
%               Example - 'Family','homo'
%               Data Types: character
% 
% Control  :   structure containing optional arguments to pass to the
%              specific FSDA function. Structure.
%              Given that FSreg is a container for most of the functions
%              present in FSDA, structure control contains the optional
%              arguments of the associated subfunction which is called.
%              The fields of control vary depending on the estimator which
%              is used. 
%               Example - 'Control',Control
%               Data Types: structure
%
% Monitoring:  Statistics to monitor. 
%              Logical value. 
%              If Monitoring is true ot is a structure the eda (explorative
%              data analysis) of the estimator is used else the standard
%              version of the estimator is called.
%              For example if estimator if 'FS', and 'Family' is 'homo',
%              when monitoring is true FSReda function is called else
%              function FSR is called. If Monitoring is a structure it
%              enables to specify which statistics must be monitored. For
%              example Monitoring.tstat=true and Monitoring.beta=true calls
%              function FSReda and enables us to monitor beta coefficients
%              and t statistics during the fwd search
%
%

% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit




%% Beginning of code

out = FSreg_fit(TBL,varargin{:});
end
