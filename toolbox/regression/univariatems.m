function [Tsel, Texcl] = univariatems(y,X,varargin)
%univariatems performs preliminary univariate robust model selection in linear regression
%
%<a href="matlab: docsearchFS('univariatems')">Link to the help function</a>
%
% In the situations in which the number of potential explanatory variables
% is very large it is necessary to preliminary understand the subset of
% variables which surely must be excluded before running the proper
% variable selection procedure. This procedure estimates a univariate
% regression model (with intercept) between each column of X and the
% response. Just the variables that have an R2 greater than (or a
% $p$-value smaller than) a certain threshold in the univariate regressions
% are retained. The p-value or R2 threshold are based on robust univariate
% models (with intercept), but the unrobust models can be chosen using
% option thresh. Option fsr enables the user to select the preferred robust
% regression procedure.
%
%
% Required input arguments:
%
%       y:      Response variable. Vector or table with just one column. A
%               vector with n elements that contains the response variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                 Data Types - vector or table
%
%       X :     Predictor variables. Matrix or table. Data matrix of
%               explanatory variables (also called 'regressors') of
%               dimension (n-by-p). Note that, in this case, p can be greater
%               than n. Rows of X represent observations, and columns
%               represent variables. Missing values (NaN's) and infinite
%               values (Inf's) are allowed, since observations (rows) with
%               missing or infinite values will automatically be excluded
%               from the computations. If X is table for the variables,
%               which are defined as categorical function dummyvar is
%               called, and the p-value (R2) refers to the univariate model
%               with all the levels (minus 1) of the categorical variable.
%                 Data Types - matrix or table
%
% Optional input arguments:
%
%        fsr  : Method to use for robust regression. Boolean.
%               If fsr is true univariate robust is based on function FSR
%               (forward search) else it is based on function LXS (LTS or
%               LMS). The default value of fsr is true.
%                 Example - 'fsr',false
%                 Data Types - logical
%
%
%         h   : The number of observations that have determined the robust
%               regression estimator for each univariate regression. Scalar.
%               h is an integer greater or equal than [(n+p+1)/2], but
%               smaller than n.
%                 Example - 'h',round(n*0,75)
%                 Data Types - double
%
%
%         lms : Criterion to use to find the robust  beta
%                 in the univariate regressions. Scalar.
%                   If lms=1 (default) Least Median of Squares is
%               computed, else Least Trimmed of Squares is computed.
%                 Example - 'lms',1
%                 Data Types - double
%
%     nsamp   : Number of subsamples that will be extracted to find the
%                 robust estimator. Scalar.
%               If nsamp=0, all subsets will be extracted.
%               They will be (n choose p).
%               Remark. If the number of all possible subset is <1000, the
%               default is to extract all subsets otherwise just 1000.
%                 Example - 'nsamp',1000
%                 Data Types - double
%
%  PredictorNames : names of the explanatory variables. Cell array of
%               characters or string array.
%               Cell array of characters or string array of length p containing the
%               names of the explanatory variables.
%               If PredictorNames is missing and X is not a table,
%               the following sequence of strings will be
%               automatically created for labelling the column of matrix X
%               ("X1", "X2", "X3", ..., "Xp").
%               Example - 'PredictorNames',{'X1','X2'}
%               Data Types - cell array of characters or string array
%
%
%   thresh:  threshold that defines the variables to retain. Scalar or struct.
%               The default value of thresh is 0.10, that is all variables
%               which in univariate robust regression had a p-value smaller
%               than 0.10 are retained.
%               If thresh is a struct it is possible to specify whether the
%               threshold is based on (robust) p-values or (robust) R2.
%               If thresh is a struct it may contain one of the the
%               following fields:
%               thresh.pval = all variables that, in univariate regression,
%                   have a p-value smaller or equal to thresh.pval are
%                   selected.
%               thresh.pvalrob = all variables that, in univariate
%               regression,
%                   have a robust p-value smaller or equal to thresh.pvalr are
%                   selected.
%               thresh.R2 = all variables that, in univariate regression,
%                   have a R2 square greater or equal to thresh.R2 are
%                   selected.
%               thresh.R2rob = all variables that, in univariate
%               regression,
%                   have a robust R2 square greater or equal to thresh.R2r are
%                   selected.
%                   Note that, if thresh is a struct with both
%                   fields an error is produced because just one between
%                   thresh.pval, thresh.pvalrob, thresh.R2, thresh.R2rob
%                   must be present.
%               Example - 'thresh',0.10
%               Data Types - double
%
%  theoreticalSigns: theoretical sign that the beta from univariate
%               regression must have. Vector of length p.
%               1 denotes a positive sign for the corresponding variable,
%               while -1 denotes a negative sign for the corresponding
%               variable. 0 denotes that any sign is possible. For example,
%               if p is equal to 5, if theoreticalSigns=[1 1 0 -1 1] means
%               that variables 1, 2, and 5 must have a positive robust
%               estimate of the slope in the univariate regression, while
%               variable 4 must have a negative estimate for the robust
%               beta coefficient. Finally, variable 3 can have any sign.
%               Example - 'theoreticalSigns',[-1 1 1 0 -1]. 
%               If theoreticalSigns is empty or it is not specified, no filter
%               based on sign is applied.
%               Data Types - double
%
%
% Output:
%
%         Tsel:   Details of the variables that were
%                 important from univariate analysis. table. 
%                 The details of table Tsel are as follows.
%                 The rownames of this table contain the names of the
%                 variables.
%               1st col: the indexes of the important variables
%               2nd-4th col: estimates of beta p-value and R2 from non robust
%                   univariate regression
%               5th-7th col: estimates of beta p-value and R2 from robust
%                   univariate regression
%               8th col: number of units declared as outliers in robust
%                   univariate regression.
%               The rows of table Tsel are ordered in terms of variable
%               importance, in the sense that row 1 refers to the variable
%               with highest robust R2 (smallest robust p-value). Row 2
%               contains the variable with the second highest robust R2
%               (second smallest robust p-value)...
%               If no explanatory variable survives the criteria, Tsel is a
%               0×8 empty table.
%
%        Texcl:   Details of the variables that were
%                 not important from univariate analysis. table.
%                 The details of table Texcl are as follows. 
%               The rownames of this table contain the names of the
%               selected variables.
%               1st col: the indexes of the important variables;
%               2nd-4th col: estimates of beta p-value and R2 from non robust
%                   univariate regression;
%               5th-7th col: estimates of beta p-value and R2 from robust
%                   univariate regression;
%               8th col: number of units declared as outliers in robust
%                   univariate regression.
%               If no explanatory variable is excluded Texcl is a
%               0×8 empty table.
%
% See also stepwiselm,lasso
%
% References:
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('univariatems')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

%{
    % Call of univariatems with all default arguments. 
    n=30;
    out=simulateLM(n,'R2',0.8);
    y=out.y;
    % Just the first 3 variables are significant.
    X=[out.X randn(n,57)];
    Tsel=univariatems(y,X);
    disp('Table with the details of the vars that are significant in univ. regr.')
%}

%{
    % Call of univariatems with a personalized threshold for p-value.
    % Generate a regression model with 3 expl. variables.
    n=30;
    out=simulateLM(n,'R2',0.5);
    y=out.y;
    % Just the first 3 variables are significant.
    X=[out.X randn(n,7)];
    % The variables that are retained are those for which the robust p-value
    % of univariate regression is not greater than 0.2
    mypval=0.20;
    Tsel=univariatems(y,X,'thresh',mypval);
%}

%{
    % Call of univariatems with a personalized threshold based on R2rob.
    n=30;
    out=simulateLM(n,'R2',0.8,'nexpl',5);
    y=out.y;
    % Just the first 10 variables are significant.
    X=[out.X randn(n,20)];
    % The variables that are retained are those for which the robust R2 
    % of univariate regression is not smaller than 0.08.
    thresh=struct;
    thresh.R2rob=0.08;
    Tsel=univariatems(y,X,'thresh',thresh);
    disp(Tsel)
%}

%{
    %% Compare stepwiselm and lasso with and without preliminary calling univariatems.
    rng(1)
    n=30;
    out=simulateLM(n,'R2',0.8);
    y=out.y;
    % Just variables 6:8 are significant.
    X=[randn(n,5) out.X randn(n,52)];
    % Use a threshold of p-values based on non robust models.
    Tsel=univariatems(y,X);
    Xtab=array2table(X(:,Tsel.i),'VariableNames',Tsel.Properties.RowNames);
    mdl=stepwiselm(Xtab,y,'Upper','linear','Verbose',0);
    fprintf('<strong>Final chosen model after univariatems filter and then stepwiselm</strong>')
    disp(mdl)
    disp('<strong>Final chosen model after applying directly stepwiselm</strong>')
    mdltrad=stepwiselm(X,y,'Upper','linear','Verbose',0);
    disp(mdltrad)
    
    % Lasso part
    cvlassoParameter=5;
    [B,FitInfo] = lasso(Xtab{:,:},y,'CV',cvlassoParameter);
    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    seqp=1:size(Xtab,2);
    minMSEModelPredictors = seqp(B(:,idxLambdaMinMSE)~=0);
    mdl=fitlm(Xtab(:,minMSEModelPredictors),y);
    disp('<strong>Final chosen model after univariatems filter and then lasso</strong>')
    disp(mdl)
    disp('<strong>Final chosen model after applying directly lasso</strong>')
    [B,FitInfo] = lasso(X,y,'CV',cvlassoParameter);
    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    seqp=1:size(X,2);
    minMSEModelPredictors = seqp(B(:,idxLambdaMinMSE)~=0);
    mdl=fitlm(X(:,minMSEModelPredictors),y,'VarNames',["X"+minMSEModelPredictors "y"]);
    disp(mdl)
%}

%{
    % Example of use of option theoreticalSigns.
    rng(1)
    n=30;
    thresh=struct;
    thresh.pval=0.10;
    out=simulateLM(n,'R2',0.8);
    y=out.y;
    % Just the first 3 variables are significant.
    X=[out.X randn(n,57)];
    % Suppose that it is known that variable 11:20 must have a positive sign in
    % univariate regressions.
    theoreticalSigns=zeros(1,size(X,2));
    theoreticalSigns(11:20)=1;
    % and that variables 41:50 must have a negative sign in
    % univariate regressions.
    theoreticalSigns(41:50)=-1;
    % call to univariatems with option theoreticalSigns.
    [Tsel,Texcl]=univariatems(y,X,'theoreticalSigns',theoreticalSigns); 
%}

%{
    % Example of X as a table and a variable inside X is categorical.
    % This example is taken from https://stats.oarc.ucla.edu/r/dae/tobit-models/
    % Consider the situation in which we have a measure of academic aptitude
    % (scaled 200-800) which we want to model using reading and math test
    % scores, as well as, the type of program the student is enrolled in
    % (academic, general, or vocational). For this example, tobit regression
    % (see function regressCens) is more appropriate. We just use this
    % dataset to show the case in which one of the variables is categorical
    % and X is table.
    link="https://stats.idre.ucla.edu/stat/data/tobit.csv";
    XX=readtable(link,"ReadRowNames",true);
    XX.prog=categorical(XX.prog);
    % Define y and X
    y=XX(:,"apt");
    X=XX(:,["read", "math" "prog"]);
    % In this case both y and X are tables
    [Tsel,Texcl]=univariatems(y,X);
%}

%% Beginning of code

if istable(y)
    y=y{:,1};
end

[n,p]=size(X);

if istable(X)
    PredictorNames=X.Properties.VariableNames;
    tableX=true;
    catColumns = varfun(@iscategorical, X, 'OutputFormat', 'uniform');
else
    catColumns=false(1,p);
    PredictorNames="X"+(1:p);
    tableX=false;
end

hmin=floor(0.5*(n+2));
h=hmin;
nsamp=2000;
lms=0;
thresh=0.10;
theoreticalSigns=[];
fsr=true;

options=struct('fsr',fsr,'h',h,'nsamp',nsamp, ...
    'lms',lms,'PredictorNames',PredictorNames,'thresh',thresh,'theoreticalSigns',theoreticalSigns);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid.
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:univariatems:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options.
    aux.chkoptions(options,UserOptions)
end

if nargin<2
    error('FSDA:univariatems:missingInputs','response y or X is missing');
end

if nargin >2
    % Write in structure 'options' the options chosen by the user.
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    % Get options chosen by the user.
    lms=options.lms;
    nsamp=options.nsamp;
    h=options.h;
    theoreticalSigns=options.theoreticalSigns;
    fsr=options.fsr;

    if tableX == false
        PredictorNames=options.PredictorNames;
    else
        if ~isempty(options.PredictorNames)
            warning('FSDA:univariatems:toomanyInputs','labels names are supplied, but are ignored because X is a table and already contains the variable names');
        end
    end
    thresh=options.thresh;
end

if h < hmin
    error('FSDA:univariatems:Wrongh',['The LTS (LMS) must cover at least ' int2str(h) ' observations.'])
elseif h >= n
    error('FSDA:univariatems:Wrongh','h is greater or equal to the number of non-missings and non-infinites.')
end


R2bet=zeros(p,8);
for j=1:p
    if tableX==true
        if catColumns(j)==true
            Xcatj=X{:,j};
            Xj = dummyvar(Xcatj);
            Xj=Xj(:,2:end);
        else
            Xj=X{:,j};
        end
    else
        Xj=X(:,j);
    end

    if sum(~ismissing(Xj))>4
        outj=fitlm(Xj,y);
        if length(outj.Coefficients.Estimate)>1
            % Non robust analysis.
            Betaj=outj.Coefficients.Estimate(2);
            R2j=outj.Rsquared.Ordinary;
            pval=coefTest(outj);

            % Robust analysis.
            if fsr== false
                outjROB=LXS(y,Xj,'msg',0,'h',h,'nsamp',nsamp,'lms',lms,'conflev',1-0.01/n);
            else
                outjROB=FSR(y,Xj,'msg',0,'init',n*0.6,'nsamp',nsamp,'lms',lms,'plots',false);
                if isnan(outjROB.outliers)
                    outjROB.outliers=[];
                end
            end
            outjR=fitlm(Xj,y,'Exclude',outjROB.outliers);
            BetajR=outjR.Coefficients.Estimate(2);
            R2jR=outjR.Rsquared.Ordinary;
            pvalR=coefTest(outjR);

            R2bet(j,:)=[j Betaj R2j pval BetajR R2jR pvalR length(outjROB.outliers)];
        else
            R2bet(j,1)=j;
        end
    else
        R2bet(j,1)=j;
    end
end

R2betT=array2table(R2bet,"RowNames",PredictorNames,'VariableNames',{'i', 'b' 'R2' 'pval' 'brob' 'R2rob' 'pvalrob' 'nout'});

% Delete all the explanatory variables that have a sign which is not
% in agreement with theory.
if ~isempty(theoreticalSigns)
    booToKeep=sign(R2betT.brob).*theoreticalSigns(:)>=0;
    R2betT=R2betT(booToKeep,:);
end

% Now we select the variable according to the criteria.
if isstruct(thresh)
    if length(fieldnames(thresh))>1
        error('FSDA:univariatems:WrgInp','struct thresh must have just one fieldname')

    elseif isfield(thresh,"pval")
        boo=R2betT.pval<=thresh.pval;
        Tsel=R2betT(boo,:);
        Tsel=sortrows(Tsel,"pval","ascend");

    elseif isfield(thresh,"pvalrob")
        boo=R2betT.pvalrob<=thresh.pvalrob;
        Tsel=R2betT(boo,:);
        Tsel=sortrows(Tsel,"pvalrob","ascend");

    elseif isfield(thresh,"R2")
        boo=R2betT.R2>=thresh.R2;
        Tsel=R2betT(boo,:);
        Tsel=sortrows(Tsel,"R2","descend");

    elseif isfield(thresh,"R2rob")
        boo=R2betT.R2rob>=thresh.R2rob;
        Tsel=R2betT(boo,:);
        Tsel=sortrows(Tsel,"R2rob","descend");

    else
        error('FSDA:univariatems:WrgInp',"the field name of thresh must be 'pval' 'pvalrob' 'R2' or 'R2rob'")
    end

else
    boo=R2betT.pvalrob<=thresh;
    Tsel=R2betT(boo,:);
    Tsel=sortrows(Tsel,"pvalrob","ascend");
end

Texcl=R2betT(~boo,:);

end

%FScategory:REG-Regression
