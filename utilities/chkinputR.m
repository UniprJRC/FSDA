function [y,X,n,p] = chkinputR(y, X, nnargin, vvarargin)
%chkinputR makes some input parameters and user options checking
%
% Required input arguments:
%
% y:            A vector with n elements that contains the response
%               variables, possibly with missing values (NaN's) and
%               infinite values (Inf's).
% X :           Data matrix of explanatory variables (also called
%               'regressors') of dimension (n x p-1), possibly with missing
%               values (NaN's) and infinite values (Inf's). Rows of X
%               represent observations, and columns represent variables.
% nnargin:      The number of input arguments specified for the caller
%               function.
% vvarargin:    The variable length input argument list specified for the
%               caller function.
%
% Output:
%
% y:            The new response variable, with observations (rows) with
%               missing or infinite values excluded.
% X:            The new matrix of explanatory variables, with missing or
%               infinite values excluded.
% n:            Number of rows of X (observations).
% p:            Number of parameters to be estimated.
%
% See also
%
% Copyright 2008-2011.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti 
%            and Vytis Kopustinskas (2009-2010)
%
% Last modified 15-Nov-2011

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
        error('Not enough input arguments.');
    end
    
    [m,q]=size(y);
    if min(m,q)>1
        error('y is not one-dimensional.');
    elseif q~=1
        
        % If y is a row vector it is transformed in a column vector
        y=y';
    end
    
    % The second argument which is passed is X
    if nnargin<2  || isempty(X)
        error('Input matrix X not specified.');
        
    % X must be a 2-dimensional array
    elseif ~ismatrix(X)
        error('Invalid data set X.');
    end
    
    % Check dimension consistency of X and y
    na.X=~isfinite(X*ones(size(X,2),1));
    na.y=~isfinite(y);
    if size(na.X,1)~=size(na.y,1)
        error('Number of observations in X and y not equal.');
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
        error(['Need more observations than variables: n=' ...
            int2str(size(X,1)) ' and p=' int2str(size(X,2)) ]);
    end
    
    rk=rank(X);
    if rk < p
        error('Matrix X is singular');
    end
end

end
