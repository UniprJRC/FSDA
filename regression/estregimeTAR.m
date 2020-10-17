function [out] = estregimeTAR(y, X)
% estregimeTAR estimate a regression model with OLS in one of the regimes of a TAR model.
%
%<a href="matlab: docsearchFS('estregimeTAR')">Link to the help function</a>
%
% Required input arguments:
%
% y :   Response variable. Vector.
%       A vector with n elements that contains the response variable.
%       Missing values (NaN's) and infinite values (Inf's) are
%       not allowed and an error message is returned.
% X :   Predictor variables. Matrix.
%       Data matrix of explanatory variables (also called
%       'regressors') of dimension (n x p). Rows of X
%       represent observations, and columns represent variables.
%       If an intercept is included in X, it must be in the last column.
%       Missing values (NaN's) and infinite values (Inf's) are
%       not allowed and an error message is returned.
%
% Optional input arguments:
%
% Output:
%
% out : A structure containing the following fields
%
%       out.beta =  Estimated parameters of the regression model. Vector. See out.covar.
%         out.se =  Estimated heteroskedasticity-consistent (HC) standard errors. Vector.
%                   See section 'More about'.
%      out.covar =  Estimated variance-covariance matrix. Matrix. It is the
%                   heteroskedasticity-consistent (HC) covariance matrix. See section 'More about'.
%    out.sigma_2 =  Estimated residual variance. Scalar.
%       out.yhat =  Fitted values. Vector.
%        out.res =  Residuals of the regression model. Vector.
%        out.RSS =  Residual Sum of Squared. Scalar.
%        out.TSS =  Total Sum of Squared. Scalar.
%        out.R_2 =  R^2. Scalar.
%          out.n =  Number of observations entering in the estimation. Scalar.
%          out.k =  Number of regressors in the model left after the checks. It is the number of
%                   betas to be estimated by OLS. The betas corresponding to the removed columns of X
%                   will be set to 0 (see section 'More about'). Scalar.
%    out.rmv_col =  Indices of columns removed from X before the model estimation. Scalar vector.
%                   Columns containing only zeros are removed. Then, to avoid multicollinearity, in the
%                   case of presence of multiple non-zero constant columns, the code leave only the first
%                   constant column (see section 'More about').
% out.rk_warning =  Warning for skipped estimation. String. If the matrix X is singular after the
%                   adjustments, the OLS estimation is skipped, the parameters are set to NaN and a
%                   warning is produced.
%
% More About:
%
% This routine performs the following operations:
% 1) If y is a row vector it is transformed in a column vector.
% 2) Checks that X is a matrix that has not more than two dimensions.
% 3) Checks dimension consistency of X and y.
% 4) Checks for missing or infinite values in X and y.
% 5) Checks if there are constant columns in matrix X. Firstly, all
% the 0-columns are removed (the associated regressors cannot enter
% in the parameter estimation step). Then, to avoid multicollinearity,
% in the case of presence of multiple non-zero constant columns, the
% code leave only the first constant column. If using the SETARX function,
% this corresponds to the preference order:
%  (i) autoregressive variables,
%  (ii) exogenous variables,
%  (iii) dummies.
% The indices of all the removed columns are saved in out.rmv_col.
% 6) Computes final values of n and k after previous operations;
% 7) Makes sure that n>=k;
% 8) Checks the rank of X before OLS estimation: if X is singular despite the adjustments,
% the OLS estimation is skipped with a warning and the parameters values are set to NaNs.
% 9) Performs OLS estimation if matrix X is not singular.
% The estimation of the covariance matrix (and standard errors) of the OLS estimator
% $\hat{\beta} = (X^{\prime} X)^{-1}X^{\prime}y$ is made via the Huber-White "Sandwich" estimator:
% $$\mathrm{cov}(\boldsymbol{\beta})=(\mathbf{X}^{\prime} \mathbf{X})^{-1}\mathbf{X}^{\prime} \mathbf{\Sigma} \mathbf{X}(\mathbf{X}^{\prime}
% \mathbf{X})^{-1}.$$
% This is a heteroskedasticity-robust variance-covariance matrix, without imposing any assumption on
% the structure of heteroskedasticity. Assuming the indipendence of the regression errors $\mathbf{u}$,
% the adjusted estimation of matrix $\Sigma$ is:
% $$\hat{\mathbf{\Sigma}}= \frac{n}{n-k} \mathrm{diag}(\mathbf{u}^2).$$
% If there is no heteroskedasticity, the robust standard errors will become just conventional
% OLS standard errors. Thus, the robust standard errors are appropriate even under homoskedasticity.
% 10) Extends the vectors out.beta and out.se with the extendVEC function.
% The beta values, corresponding to the removed columns of X, are set to 0.
% The se values, corresponding to the removed columns of X, are set to NaN.
%
% See also extendVEC, SETARX
%
%
% References:
%
% Copyright 2008-2020.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('estregimeTAR')">Link to the help function</a>
%
%
%$LastChangedDate:: 2020-06-09 17:36:50 $: Date of the last commit
%
% Example:
%{
    %% Example 1: estregimeTAR with all default options.
    rng('default')
    rng(10)
    n=200;
    k=3;
    X=randn(n,k);
    y=randn(n,1);
        
    X1=[X [1:200]' [1:200]'];
    [out1] = estregimeTAR(y, X1);
%}

%{
    % Example 2: adjustments for constant columns. Only the first non-zero
    %   constant column is kept in the model estimation. Check beta and se values.
    n=200;
    k=3;
    X=randn(n,k);
    y=randn(n,1);   
    X2 = [zeros(200,1) ones(200,1) X repmat(2,200,1)];
    [out2] = estregimeTAR(y, X2);
%}


%% Beginning of code

% The first argument which is passed is y.
if nargin<1 || isempty(y)==1
    error('FSDA:estregimeTAR:MissingInputs','Input vector y not specified.');
end

% y must be a vector.
[nn,kk] = size(y);
if min(nn,kk)>1
    error('FSDA:estregimeTAR:Wrongy','y must be a vector.');
elseif kk~=1
    % If y is a row vector it is transformed in a column vector.
    y=y';
end


% The second argument which is passed is X.
if nargin<2 || isempty(X)==1
    error('FSDA:estregimeTAR:MissingInputs','Input matrix X not specified.');
end

% X must be a 2-dimensional array.
if ~isempty(X) && ~ismatrix(X)
    error('FSDA:estregimeTAR:WrongX','Invalid data set X.');
end


% y and X must be different.
if isequal(y,X)
    error('FSDA:estregimeTAR:XyEqual','Invalid input: y and X are equal.');
end


% Check dimension consistency of X and y.
na.X = ~isfinite(X*ones(size(X,2),1));
na.y = ~isfinite(y);
if size(na.X,1)~=size(na.y,1)
    error('FSDA:estregimeTAR:NxDiffNy','Number of observations in X and y not equal.');
end


% Check if observations with missing or infinite values are present.
bad = (na.X|na.y);
if sum(bad)>0
    error('FSDA:estregimeTAR:NaNInfFound','Missing or infinite values found.');
end



% rmv_col is the scalar vector of the indices of constant columns, possible to remove.
rmv_col = find(max(X,[],1)-min(X,[],1) == 0);

if numel(rmv_col)>0
    % Remove all the 0-columns (the associated regressors cannot enter in the parameter estimation step).
    % Leave in the matrix only the first column with constant values ~= 0 (all the others cause multicollinearity).
    % In the presence of multiple constant columns ~= 0, the code keeps the first one ...
    % following the order of importance of the regressors: as instance AR, Exogenous, ...
    % Dummy/Deterministic, Intercept.
    constkeep = find(X(1,rmv_col) ~= 0); % Find which of the constant columns are ~= 0.
    if ~isempty(constkeep)
        rmv_col(constkeep(1)) = []; % Keep only the first of the constant columns ~= 0 (remove it from the list of the columns to remove).
    end
    X(:,rmv_col) = []; % Remove all the other constant columns.
    % disp(['Warning: columns ' num2str(rmv_col) ' are constant and have been removed!']);
end


% n is the number of observations and k is the number of parameters to be estimated.
[n, k] = size(X);

if n < k
    error('FSDA:estregimeTAR:NSmallerK',['Need more observations than variables: n = ' ...
        int2str(n) ' and k = ' int2str(k) ]);
end

TSS = sum((y - mean(y)).^2); % Total Sum of Squares (variance of dependent variable).

rk = rank(X); % Check the rank of X
rk_warning = '';
if rk < k
    rk_warning = sprintf('Warning: Matrix X is singular despite the adjustments! OLS estimation is skipped.');
    beta = NaN; % Parameter.
    yhat = NaN; % Fitting.
    res = NaN; % Residuals.
    RSS = Inf; % Residual Sum of Squares.
    sigma_2 = NaN; % Residual variance.
    covar = NaN; % Covariance matrix.
    se = NaN; % Standard errors.
    R_2 = NaN; % R^2
else
    % Estimate OLS parameters
    beta = X\y; % Parameter estimation with OLS.
    yhat = X*beta; % Fitting.
    res = y - yhat; % Residuals.
    RSS = res'*res; % Residual Sum of Squares.
    sigma_2 = RSS/(n-k); % Residual variance.
    
    
    % Huber-White "Sandwich" estimator: heteroskedasticity-robust variance-covariance matrix.
    %XXinv = inv(X'*X);
    %covar = XXinv*X'*diag(res.*res.*(n/(n-k)))*X*XXinv;
    covar = (X'*X)\X'*diag(res.*res.*(n/(n-k)))*X/(X'*X);
    % Heteroskedasticity-robust standard errors. It works also in the case of homoskedasticity.
    se = sqrt(diag(covar));
    % covar_hom = XXinv.*sigma_2; % Covariance matrix. Assumption: homoskedasticity.
    % se_hom = sqrt(diag(covar_hom)); % Standard errors. Assumption: homoskedasticity.
    
    R_2 = 1 - RSS/TSS; % R^2.
    
    % Insert the missing parameter values, corresponding to the removed columns of X, with extendVEC()
    if numel(rmv_col)>0
        beta = extendVEC(beta, rmv_col, 0);
        se = extendVEC(se, rmv_col, NaN);
    end
    
end



% Structured outputs
out = struct;

out.beta = beta;
out.se = se;
out.covar = covar;
out.sigma_2 = sigma_2;
out.yhat = yhat;
out.res = res;
out.RSS = RSS;
out.TSS = TSS;
out.R_2 = R_2;
out.n = n;
out.k = k;
out.rmv_col = rmv_col;
out.rk_warning = rk_warning;

end
%FScategory:REG-Regression