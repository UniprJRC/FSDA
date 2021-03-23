function [X,n,p] = chkinputM(X, nnargin, vvarargin)
%chkinputM makes some input parameters and user options checking in multivariate analysis
%
% Required input arguments:
%
% X :          Input data. Matrix.
%               n x v data matrix; n observations and v variables. Rows of
%               X represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
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
% X:            Data matrix without missing and infs. Matrix.
%               The new matrix of variables, with missing or
%               infinite values excluded.
% n:            Number of rows of X (observations). Scalar.  Number of
%               rows after listwise exclusion.
% p:            Number of columns of X (variables). Scalar.
%               Number of variable in the input data matrix.
%
% See also
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Example:
%{
    %% example_producing_error
    %To examplify the behaviour of chkinputM, we call function FSM with a
    %X with more columns then rows.
    n=3;
    p=200;
    state1=123498;
    randn('state', state1);
    X=randn(n,p);
    [out]=FSM(X);
%}
%% Beginning of code
optargin = size(vvarargin,2);
stdargin = nnargin - optargin;

if nnargin > stdargin
    
    if coder.target('MATLAB')
        % chklist is the vector containing the names of optional arguments
        chklist=vvarargin(1:2:length(vvarargin));
    else
        chklist=cell(1,length(vvarargin)/2);
        ij=0;
        for i=1:2:(length(vvarargin)-1)
            ij=ij+1;
            chklist{ij}=vvarargin{i};
        end
    end
    % chkchk is the position of the option nocheck in vector chklist
    % chkchk = strmatch('nocheck',chklist,'exact');
    chkchk = find(strcmpi('nocheck',chklist));
else
    %chkchk and chkint are empty if not specified by the user
    chkchk='';
end

% If nocheck=1, then skip checks on y and X
if ~isempty(chkchk) && vvarargin{2*chkchk}==1
    [n,p]=size(X);
else
    
    % The second argument which is passed is X
    if nnargin<1 || isempty(X)
        error('FSDA:chkinputM:missingInputs','Input matrix X not specified.');
        
        % X must be a 2-dimensional array
    elseif ~ismatrix(X)
        error('FSDA:chkinputM:WrongX','Invalid data set X.');
    end
    
    % Check dimension consistency of X and y
    na.X=~isfinite(X*ones(size(X,2),1));
    
    % Observations with missing or infinite values are removed from X and y
    ok=~(na.X);
    X=X(ok,:);
    % Now n is the new number of non missing observations
    n=length(X);
    % constcols = scalar vector of the indices of possible constant columns.
    constcols = find(max(X,[],1)-min(X,[],1) == 0);
    if numel(constcols)>=1
        X(:,constcols(1:end))=[];
        if coder.target('MATLAB')
            disp(['Warning: columns ' num2str(constcols) ' are constant and just col ' num2str(constcols(1)) ' has been kept!']);
        else
            disp('Warning: columns constant and just on of them has been kept!');
        end
    end
    
    % p is the number of parameters to be estimated
    p=size(X,2);
    
    if n < p
        error('FSDA:chkinputM:NSmallerP',['Need more observations than variables: n=' ...
            int2str(size(X,1)) ' and p=' int2str(size(X,2)) ]);
    end
    
    rk = rank(X);
    if rk < p
        error('FSDA:chkinputM:NoFullRank','Matrix X is singular');
    end
end

end