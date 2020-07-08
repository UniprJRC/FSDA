function [y,q,X,n,k,rmv_obs,input_full] = chkinputTAR(y, q, X, intercept)
% chkinputTAR makes some input parameters checking.
%
%<a href="matlab: docsearchFS('chkinputTAR')">Link to the help function</a>
%
% Required input arguments:
%
% y:            Response variable. Vector.
%               A vector with n elements that contains the response
%               variable, possibly with missing values (NaN's) and
%               infinite values (Inf's).
% q:            Threshold variable. Vector.
%               A vector with n elements that contains the threshold
%               variable, possibly with missing values (NaN's) and
%               infinite values (Inf's).
% X :           Predictor variables. Matrix.
%               Data matrix of explanatory variables (also called
%               'regressors') of dimension (n x p-1), possibly with missing
%               values (NaN's) and infinite values (Inf's). Rows of X
%               represent observations, and columns represent variables.
%               X can be empty, but in this case intercept should be set equal to true.
% intercept:    Indicator for constant term, true to include or false to not include
%               the constant term in the model.
%
%
% Optional input arguments:
%
% Output:
%
% y:            Response without missing and infs. Vector. The new response variable, with observations (rows) with
%               missing or infinite values excluded.
% q:            Threshold variable without missing and infs. Vector. The new threshold variable,
%               with observations (rows) with missing or infinite values excluded.
% X:            Predictor variables without infs and missings. Matrix.
%               The new matrix of explanatory variables, with missing or
%               infinite values excluded.
% n:            Number of rows of X (observations). Scalar.  Number of
%               rows after listwise exclusion.
% k:            Number of columns of X (variables). Scalar.
%               Number of parameters to be estimated.
% rmv_obs:      Indices of removed observations/rows (because of missings or infs). Scalar vector.
% input_full:   A structure containing y, q and X after adjustements BUT with observations (rows) with
%               missing or infinite values included.
%
%
% More About:
%
% This routines preforms the following operations:
% 1) If y and q are row vectors they are transformed in column vectors;
% 2) Checks that X is a matrix that has not more than two dimensions;
% 3) Adds to matrix X a column of ones if option intercept is 1;
% 4) Checks dimension consistency of X, y and q;
% 5) Removes observations with missing or infinite values from X, y or q
% (listwise exclusion);
% 6) Checks if there are constant columns in matrix X. In other words, if
% Xj is a generic column of X (excluding the column which contains the
% intercept) it removes it if max(Xj)=min(Xj) and produces a warning.
% 7) Computes final values of n and k after previous operations;
% 8) Makes sure that n>=k;
% 9) Makes sure that new X is full rank.
%
% See also SETARX
%
% Copyright 2008-2020.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('chkinputTAR')">Link to the help function</a>
%
%$LastChangedDate:: 2020-06-09 17:02:50 $: Date of the last commit
%
%
%
% Example:

%{
 % Example 1: error - first input is not a vector.
    rng(10)
    n=200;
    k=3;
    X=randn(n,k);
    [out]=chkinputTAR(X);
%}

%{
% Example 2: error - y and q are equal.
    rng(10)
    n=200;
    k=3;
    X=randn(n,k);
    y=randn(n,1);
    [out] = chkinputTAR(y, y, X);
%}

%{
    % Example 3: error - no regressors in the model.
    rng(10)
    n=200;
    k=3;
    X=randn(n,k);
    y=randn(n,1);    [out] = chkinputTAR(y, q, [], 'false');
%}

%{
%% Example 4: this works.
    rng(10)
    n=200;
    k=3;
    X=randn(n,k);
    y=randn(n,1);
    const = 'true';
    [y,q,X,n,k,rmv_obs,input_full] = chkinputTAR(y, q, X, const);
%}

%{
    % Example 5: this works with a warning on extra constant columns.
    % Only the last one corresponding to the intercept is kept.
    rng(10)
    n=200;
    k=3;
    X=randn(n,k);
    y=randn(n,1);
    XX=[X ones(n,k)];
    [y,q,X,n,k,rmv_obs,input_full] = chkinputTAR(y, q, XX, const);
%}

%{
%% Example 6: this works and rows with missing observations are removed.
     rng(10)
    n=200;
    k=3;
    X=randn(n,k);
    y=randn(n,1);
    q(3) = NaN;
    X(20) = NaN;
    const = 'true';
    [y,q,X,n,k,rmv_obs,input_full] = chkinputTAR(y, q, X, const);
%}


%% Beginning of code

% The first argument which is passed is y.
if nargin<1 || isempty(y)==1
    error('FSDA:chkinputTAR:MissingInputs','Input vector y not specified.');
end

% y must be a vector.
[nn,kk] = size(y);
if min(nn,kk)>1
    error('FSDA:chkinputTAR:Wrongy','y must be a vector.');
elseif kk~=1
    % If y is a row vector it is transformed in a column vector.
    y=y';
end

% The second argument which is passed is q.
if nargin<2 || isempty(q)==1
    error('FSDA:chkinputTAR:MissingInputs','Input vector q not specified.');
end

% q must be a vector
[nn,kk] = size(q);
if min(nn,kk)>1
    error('FSDA:chkinputTAR:Wrongq','q must be a vector.');
elseif kk~=1
    % If q is a row vector it is transformed in a column vector.
    q=q';
end



% The dependent variable y and the threshold variable q cannot be equal.
if isequal(y,q)
    error('FSDA:chkinputTAR:yqEqual','Invalid input: y and q are equal.');
end



% The third argument which is passed is X.
if nargin<3
    error('FSDA:chkinputTAR:MissingInputs','Input matrix X not specified.');
end

% X must be a 2-dimensional array.
if ~isempty(X) && ~ismatrix(X)
    error('FSDA:chkinputTAR:WrongX','Invalid data set X.');
end



% The fourth argument which is passed is intercept.
if nargin<4 || isempty(intercept)==1
    error('FSDA:chkinputTAR:MissingInputs','Input indicator intercept not specified.');
    % If the value of intercept is equal to 1, add to matrix X the column of ones.
elseif isequal(intercept,1)
    if isempty(X)
        X = ones(size(y));
    else
        X = cat(2,X,ones(size(X,1),1));
    end
else
    if isempty(X)
        % At least one regressor is needed, so if input X is empty, the intercept should be 1.
        error('FSDA:chkinputTAR:NoRegressors','Missing regressors: Input matrix X is empty and intercept is false.');
    end
end



% y and X must be different.
if isequal(y,X)
    error('FSDA:chkinputTAR:XyEqual','Invalid input: y and X are equal.');
end



% Check dimension consistency of X, y and q.
na.X = ~isfinite(X*ones(size(X,2),1));
na.y = ~isfinite(y);
na.q = ~isfinite(q);
if size(na.X,1)~=size(na.y,1)
    error('FSDA:chkinputTAR:NxDiffNy','Number of observations in X and y not equal.');
end
if size(na.X,1)~=size(na.q,1)
    error('FSDA:chkinputTAR:NxDiffNq','Number of observations in X and q not equal.');
end

% Save X, y and q before removing missing and infinite values
input_full.X = X;
input_full.y = y;
input_full.q = q;

% Observations with missing or infinite values are removed from X, y and q.
ok = ~(na.X|na.y|na.q);
X = X(ok,:);
y = y(ok,:);
q = q(ok,:);
% rmv_obs = scalar vector of the indices of removed rows (observations).
rmv_obs = find(ok == 0);
if numel(rmv_obs)>0
    disp(['Warning: observations in positions ' num2str(rmv_obs.') ' have been removed!']);
end


% Now n is the new number of observations
n = length(y);


% constcols = scalar vector of the indices of possible constant columns.
constcols = find(max(X,[],1)-min(X,[],1) == 0);
if numel(constcols)>1
    X(:,constcols(1:end-1)) = [];
    input_full.X(:,constcols(1:end-1)) = [];
    disp(['Warning: columns ' num2str(constcols) ' are constant and just col ' num2str(constcols(end)) ' has been kept!']);
end


% k is the number of parameters to be estimated.
k = size(X,2);

if n < k
    error('FSDA:chkinputTAR:NSmallerK',['Need more observations than variables: n = ' ...
        int2str(size(X,1)) ' and k = ' int2str(size(X,2)) ]);
end

rk = rank(X);
if rk < k
    error('FSDA:chkinputTAR:NoFullRank','Matrix X is singular');
end

end
