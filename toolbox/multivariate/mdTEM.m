function out = mdTEM(Y, varargin)
% mdTEM  EM algorithm with trimming (TEM) for data with missing values.
%
%
%<a href="matlab: docsearchFS('mdTEM')">Link to the help function</a>
%
% The algorithm:
%  - At each iteration compute adjusted partial Mahalanobis distances;
%  - Rank them and set weights w_i = 1 for the lowest n*(1-alpha) rows,
%       else 0;
%  - Run E-step and M-step using these weights;
%  - Repeat until convergence or maxiter.
%
%
% Required input arguments:
%
% Y :           Input data. Matrix. n x p data matrix; n observations and v
%               variables possibly with missing values (NaN's). Rows of Y
%               represent observations, and columns represent variables.
%               Data Types - single | double
%
%  Optional input arguments:
%
%       alpha   : proportion to trim. Real number in the interval [0 0.5]
%                 or empty value.
%                 At each iteration compute adjusted partial Mahalanobis
%                 distance and set weights w_i = 1 for the lowest
%                 n*(1-alpha) rows. (e.g., 0.5 -> keep 50% with smallest
%                 distances). If alpha is empty the default value which is
%                 used is 0.5.
%                 Example - 'alpha',0.1
%                 Data Types - single | double
%
%       mus     : initial mean. p x 1 vector or empty double.
%                 Initial  mean vector. If empty (default), column nanmeans
%                 are used.
%                 Example - 'mus',[]
%                 Data Types - single | double
%
%           sigs:  initial covariance matrix.
%                  p x p matrix or empty double.
%                  Initial p x p covariance matrix.
%                  If empty, uses nan-cov
%                 Example - 'sigs',eye(p)
%                 Data Types - single | double
%
%        maxiter:  maximum number of iterations. Positive integer.
%                  The default value is 100
%                 Example - 'maxiter',50
%                 Data Types - single | double
%
%           tol :  tolerance for convergence. Positive real number.
%                  The default value of the tolerance is 1e-5
%                 Example - 'tol',1e-10
%                 Data Types - single | double
%
%   tol_sigma  :   Use tolerance for both mu sigs. Boolean .
%                  If true use both mu and sigma diffs (default true)
%                 Example - 'tol_sigma',false
%                 Data Types - logical
%
%   method : method used to rescale the distances. String scalar or char vector.
%            Possible values are.
%
%            'pri'      = principled EM rescaling (default),
%                         d2_partial + (p - pobs).
%
%            'expScale' = expectation scaling,
%                         d2_partial * (p / pobs).
%
%            'zMap'     = standardization mapping,
%                         p + sqrt(2*p) * ((d2_partial - pobs) ./ sqrt(2*pobs)).
%
%            'detMap'   = determinant-based rescaling,
%                         d2_partial * (p / pobs) * (g_full / g_obs).
%
%            'chiMap'   = chi-square quantile mapping. Use the cdf and
%                         inverse of the cdf of Chi2 distribution.
%
%            'betaMap'  = Beta quantile mapping. Use the cdf and
%                         inverse of the cdf of Beta distribution.
%            'impMD'    = MD on EM-imputed data.
%            Example - 'method','chiMap'
%            Data Types - string scalar | char vector
%
%   condmeanimp :  Also give the matrix of conditional mean imputed values. Boolean.
%                 if true structure out also contains the matrix of imputed values called Yimp.
%                 The default value of condmeanimp is false.
%                 Example - 'condmeanimp',true
%                 Data Types - logical
%
%   stochimp :     Also give the matrix of stochastic imputed values. Boolean.
%                 if true structure out also contains the matrix of imputed values called stochYimp.
%                 The default value of stochimp is false.
%                 Example - 'stochimp',true
%                 Data Types - logical
%
%   storeobj :    Compute value of the objective function in each iteration.
%                 Also give the matrix of stochastic imputed values. Boolean.
%                 if true structure out also contains the matrix of imputed values called stochYimp.
%                 The default value of stochimp is false.
%                 Example - 'stochimp',true
%                 Data Types - logical
%
%
%  Output:
%
%
%         out:   structure which contains the following fields
%              out.loc = final estimates of means
%              out.cov = final estimate of cov matrix
%              out.iter = number of iterations to convergence.
%              out.Yimp = empty value of matrix Y with imputed values
%                   (depending on input option condmeanimp)
%              out.stochYimp = empty value of matrix Y with imputed values
%                   (only if input option stochimp is true)
%              out.obj  = empty value or value of the objective function (trimmed sum of
%                   smallest MD) in each iteration  
%                   (only if input option storeobj is true)
%
%
% See also: mdEM, mdImputeCondMean.m, mdPartialMD.m, mdPartialMD2full
%
% References:
%
% Little, R. J. A., & Rubin, D. B. (2019). Statistical Analysis with
% Missing Data (3rd ed.). Hoboken, NJ: John Wiley & Sons.
% van Buuren, S. (2018). Flexible Imputation of Missing Data (2nd ed.).
% Boca Raton, FL: Chapman & Hall/CRC (Taylor & Francis Group).
% Templ, M. (2023). Visualization and Imputation of Missing Values: With
% Applications in R. Cham, Switzerland: Springer Nature.
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('mdTEM')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Call to mdTEM with all the default options.
    % True model (choose something correlated)
    p=5; n=200;
    A = randn(p);
    SigmaTrue = A'*A;
    D = diag(1 ./ sqrt(diag(SigmaTrue)));
    SigmaTrue = D * SigmaTrue * D;      % "correlation-like"
    muTrue = linspace(-1,1,p)';
    
    %  generate complete data
    Yfull = mvnrnd(muTrue', SigmaTrue, n);             % n x p
    missRate = 0.25;     % MCAR missing probability per entry
    missMask = rand(n,p) < missRate;
    Y=Yfull;
    Y(missMask) = NaN;
    out=mdTEM(Y);
    % Show true means and inputed means
    scatter(out.loc,muTrue)
    refline(1)
    xlabel('Imputed means')
    ylabel('True means')
%}

%{
    %% Example of use of option condmeanimp.
    % number of variables
    p = 15;                
    % number of observations
    n = 1000;            
    % target pairwise correlation (0<rho<1)
    rho = 0.9;            
    % Covariance matrix (unit variances)
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);      % upper-triangular such that Sigma = R'*R
    % Generate samples ~ N(0,Sigma)
    Yfull = randn(n,p) * R;   % Strong positive correlation between the vars
    missRate = 0.25;     % MCAR missing probability per entry
    missMask = rand(n,p) < missRate;
    Y=Yfull;
    Y(missMask) = NaN;
    % md with missing imputation
    out=mdTEM(Y,'condmeanimp',true);
    % Mahalanobis distances using original matrix
    d2Ori=mahalFS(Yfull,mean(Yfull),cov(Yfull));
    % Calculate the Mahalanobis distance for the imputed data
    d2Imp = mahalFS(out.Yimp, mean(out.Yimp), cov(out.Yimp));
    % Compare original with distances for the imputed data
    % Calculate the differences between original and imputed Mahalanobis distances
    scatter(d2Ori,d2Imp)
    % Add axis labels
    xlabel('Original Mahalanobis Distances');
    ylabel('Imputed Mahalanobis Distances');
    grid on
%}

%% Beginning of code
alpha=0.5;
mus=[];
sigs=[];
maxiter=100;
tol=1e-5;
tol_sigma=true;
condmeanimp=false;
stochimp=false;
storeobj=true;
method='pri';

if nargin>1
    options=struct('storeobj',storeobj,'alpha',alpha,'mus',mus,'sigs',sigs,'maxiter',maxiter,'tol',tol, ...
        'tol_sigma',tol_sigma,'condmeanimp',condmeanimp,'stochimp',stochimp,'method',method);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:mdTEM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        aux.chkoptions(options,UserOptions)
    end

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    alpha=options.alpha;
    mus=options.mus;
    sigs=options.sigs;
    maxiter=options.maxiter;
    tol=options.tol;
    tol_sigma=options.tol_sigma;
    condmeanimp=options.condmeanimp;
    stochimp=options.stochimp;
    method=string(options.method);
    storeobj=options.storeobj;
end

if storeobj==true
    obj=zeros(maxiter,1);
end
[n, p] = size(Y);

% initialize mus and sigs if not provided
if isempty(mus)
    mus = mean(Y,1,"omitmissing")';       % p x 1
end
if isempty(sigs)
    X0 = Y;
    for j = 1:p
        miss = isnan(X0(:,j));
        X0(miss,j) = mus(j);
    end
    sigs = cov(X0, 1);
end

dif = Inf;
iter = 0;

% number to keep:
keep_count = max(0, floor(n * (1 - alpha)));

while (dif > tol) && (iter < maxiter)
    iter = iter + 1;
    mus_old = mus;
    sigs_old = sigs;


    if method=="impMD"
        Yimp=mdImputeCondMean(Y, mus, sigs);
        % In this case compute Mahalanobis distances on imputed data
        d2_adj=mahalFS(Yimp,mus',sigs);

    elseif method=="detMap"
        [d2, poss] = mdPartialMD(Y, mus, sigs);
        d2_adj = mdPartialMD2full(d2, p, poss,'method',method,'Y',Y,'Sigma',sigs);

    else
        % Trimming step: compute adjusted partial Mahalanobis distances
        [d2, poss] = mdPartialMD(Y, mus, sigs);
        d2_adj = mdPartialMD2full(d2, p, poss,'method',method);
    end

    % rank and select the smallest n*(1-alpha)
    % We treat NaN distances as large (so they're trimmed)
    nan_mask = isnan(d2_adj);

    % find indices of smallest distances
    % create sorted index from available (non-NaN) adj distances
    [~, idx_sorted] = sort(d2_adj, 'ascend', 'MissingPlacement', 'last');
    keep_idx = idx_sorted(1:min(keep_count, sum(~nan_mask)));


    w = zeros(n,1);
    w(keep_idx) = 1;

    % E-step with weights w
    [T1, T2] = aux.NAcompute_expected_stats(Y, mus, sigs, w);
    % M-step
    [mus, sigs] = aux.NAmaximization_step(T1, T2, w);

    % Apply consistency factor
    mm=sum(w);
    a=chi2inv(mm/n,p);
    corr=(n./mm).*(chi2cdf(a,p+2));
    sigs=sigs/corr;

    if storeobj==true
        obj(iter)=sum(d2_adj(keep_idx))/corr;
    end

    % convergence check
    mu_diff = max(abs(mus(:) - mus_old(:)));
    sigma_diff = max(abs(sigs(:) - sigs_old(:)));
    if tol_sigma
        dif = max(mu_diff, sigma_diff);
    else
        dif = mu_diff;
    end
end

%% EM imputation of missing values (conditional means or stochastic imputation)
if condmeanimp ==true
    Yimp = mdImputeCondMean(Y, mus, sigs);
else
    Yimp=[];
end

if stochimp == true
    stochYimp = mdImputeStochastic(Y, mus, sigs);
else
    stochYimp=[];
end

if storeobj==true
    obj=obj(1:iter);
end

out.loc = mus;
out.cov = sigs;
out.iter = iter;
out.Yimp=Yimp;
out.stochYimp=stochYimp;
out.obj=obj;

end

%FScategory:MULT-MissingData




