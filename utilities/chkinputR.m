function [y,X,n,p] = chkinputR(y, X, nnargin, vvarargin)
%chkinputR makes some input parameters and user options checking in regression
%
% Required input arguments:
%
% y:            Response variable. Vector.
%               A vector with n elements that contains the response
%               variables, possibly with missing values (NaN's) and
%               infinite values (Inf's).
% X :           Predictor variables. Matrix.
%               Data matrix of explanatory variables (also called
%               'regressors') of dimension (n x p-1), possibly with missing
%               values (NaN's) and infinite values (Inf's). Rows of X
%               represent observations, and columns represent variables.
% nnargin:      nargin. Scalar. The number of input arguments specified for the caller
%               function.
% vvarargin:    nvarargin. Scalar. The variable length input argument list
%               specified for the
%               caller function.
%
%
%  Optional input arguments:
%
% Output:
%
% y:            response without missing and infs. Vector. The new response variable, with observations (rows) with
%               missing or infinite values excluded.
% X:            Predictor variables without infs and missings. Matrix.
%               The new matrix of explanatory variables, with missing or
%               infinite values excluded.
% n:            Number of rows of X (observations). Scalar.  Number of
%               rows after listwise exclusion.
% p:            Number of columns of X (variables). Scalar.
%               Number of parameters to be estimated.
%
%
% More About:
%
% This routines preforms the following operations:
% 1) If y is a row vector it is transformed in a column vector;
% 2) Checks that X is a 2-dimensional array;
% 3) Checks dimension consistency of X and y;
% 4) Removes observations with missing or infinite values from X or y
% (listwise exclusion);
% 5) Adds to matrix X a column of ones if option intercept is 1;
% 6) Checks if there are constant columns in matrix X. In other words, if
% Xj is a generic column of X (excluding the column which contains the
% intercept) it removes it if max(Xj)=min(Xj) and produces a warning.
% 7) Computes final values of n and p after previous operations;
% 8) Makes sure than n>=p;
% 9) Makes sure that new X is full rank
%
% See also chkinputRB
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
% Last modified 06-Feb-2015
%
% Example:
%{
%% example_producing_error
    %To examplify the behaviour of chkinputR, we call function FSR without a
    %compulsory parameter ('y').
    n=200;
    p=3;
    state1=123498;
    randn('state', state1);
    X=randn(n,p);
    [out]=FSR(X);
%}
optargin = size(vvarargin,2);
stdargin = nnargin - optargin;

if nnargin > stdargin
    
    % chklist is the vector containing the names of optional arguments
    chklist=vvarargin(1:2:length(vvarargin));
    
    % chkchk is the position of the option nocheck in vector chklist
    % chkchk = strmatch('nocheck',chklist,'exact');
    chkchk = find(strcmpi('nocheck',chklist));
    
    % chkint is the position of the option intercept in vector chklist
    % chkint = strmatch('intercept',chklist,'exact');
    chkint = find(strcmpi('intercept',chklist));
else
    
    %chkchk and chkint are empty if not specified by the user
    chkchk='';
    chkint='';
end

% If nocheck=1, then skip checks on y and X
if ~isempty(chkchk) && vvarargin{2*chkchk}==1
    [n,p]=size(X);
else
    
    % The first argument which is passed is y
    if nnargin<1 || isempty(y)==1
        error('FSDA:chkinputR:missingInputs','Input vector y not specified.');
    end
    
    [m,q]=size(y);
    if min(m,q)>1
        error('FSDA:chkinputR:Wrongy','y is not one-dimensional.');
    elseif q~=1
        
        % If y is a row vector it is transformed in a column vector
        y=y';
    end
    
    
    % The second argument which is passed is X
    if nnargin<2  || isempty(X)
        error('FSDA:chkinputR:missingInputs','Input matrix X not specified.');
        
        % X must be a 2-dimensional array
    elseif ~ismatrix(X)
        error('FSDA:chkinputR:WrongX','Invalid data set X.');
    end
    
    % Check dimension consistency of X and y
    na.X=~isfinite(X*ones(size(X,2),1));
    na.y=~isfinite(y);
    if size(na.X,1)~=size(na.y,1)
        error('FSDA:chkinputR:NxDiffNy','Number of observations in X and y not equal.');
    end
    
    % Observations with missing or infinite values are removed from X and y
    ok=~(na.X|na.y);
    X=X(ok,:);
    y=y(ok,:);
    
    % Now n is the new number of non missing observations
    n=length(y);
    
    
    % Now add to matrix X a column of ones for the intercept.
    if nnargin <= stdargin
        
        % If the user has not specified a value for the intercept than add
        % a column of ones.
        X = cat(2,ones(n,1),X);
    else
        
        % If a value for the intercept has not been specified or if this value is
        % equal to 1, add to matrix X the column of ones. The position of
        % the option intercept in chklist, which contains the optional is
        % given in chkint. chkint is empty if the option intercept is not
        % specified.
        if isempty(chkint) || vvarargin{2*chkint}==1
            X = cat(2,ones(n,1),X);
        end
    end;
    
    
    % constcols = scalar vector of the indices of possible constant columns.
    constcols = find(max(X,[],1)-min(X,[],1) == 0);
    if numel(constcols)>1
        X(:,constcols(2:end))=[];
        disp(['Warning: columns ' num2str(constcols) ' are constant and just col ' num2str(constcols(1)) ' has been kept!']);
    end
    
    
    % p is the number of parameters to be estimated
    p=size(X,2);
    
    if n < p
        error('FSDA:chkinputR:NSmallerP',['Need more observations than variables: n=' ...
            int2str(size(X,1)) ' and p=' int2str(size(X,2)) ]);
    end
    
    rk=rank(X);
    if rk < p
        error('FSDA:chkinputR:NoFullRank','Matrix X is singular');
    end
end

end

%FScategory:UTIGEN