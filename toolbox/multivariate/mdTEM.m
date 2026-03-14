function out = mdTEM(Y, varargin)
% mdTEM  EM algorithm with trimming (TEM) for data with missing values.
%
%
%<a href="matlab: docsearchFS('mdTEM')">Link to the help function</a>
%
% The algorithm:
%  - At each iteration compute adjusted partial Mahalanobis distances
%  - Rank them and set weights w_i = 1 for the lowest n*(1-alpha) rows, else 0
%  - Run E-step and M-step using these weights
%  - Repeat until convergence or maxiter
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
%  typeAdj : method to rescale. Positive integer.
%             Method which must be used to rescale the variables.
%             typeAdj=1 use of asymptotic Chi2 distribution.
%             typeAdj=2 use of exact Beta distribution.
%             typeAdj=3 use of expectation correction
%             typeAdj=4 use of standardization correction.
%                 Example - 'typeAdj',1
%                 Data Types - single | double
%
%   imputation :  Also give the matrix of imputed values. Boolean.
%                 if true structure out also contains the matrix of imputed values.
%                 The default value of imputation is false.
%                 Example - 'imputation',true
%                 Data Types - logical
%
%  Output:
%
%
%         out:   structure which contains the following fields
%              out.loc = final estimates of means
%              out.cov = final estimate of cov matrix
%              out.iter = number of iterations to convergence.
%              out.Yimp = empty value of matrix Y with imputed values
%                   (depending on input option imputation)
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

%% Beginning of code
alpha=0.5;
mus=[];
sigs=[];
maxiter=100;
tol=1e-5;
tol_sigma=true;
imputation=false;
typeAdj=2;

if nargin>1
 options=struct('alpha',alpha,'mus',mus,'sigs',sigs,'maxiter',maxiter,'tol',tol, ...
     'tol_sigma',tol_sigma,'imputation',false,'typeAdj',typeAdj);

 [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:emNA:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    imputation=options.imputation;
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

while (dif > tol) && (iter < maxiter)
    iter = iter + 1;
    mus_old = mus;
    sigs_old = sigs;

    % Trimming step: compute adjusted partial Mahalanobis distances
    [d2, poss] = NApartialMD(Y, mus, sigs);
    d2_adj = NApartialMD2full(d2, p, poss,'typeAdj',typeAdj);

    % rank and select the smallest n*(1-alpha)
    % We treat NaN distances as large (so they're trimmed)
    nan_mask = isnan(d2_adj);
   
    % number to keep:
    keep_count = max(0, floor(n * (1 - alpha)));
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

    % convergence check
    mu_diff = max(abs(mus(:) - mus_old(:)));
    sigma_diff = max(abs(sigs(:) - sigs_old(:)));
    if tol_sigma
        dif = max(mu_diff, sigma_diff);
    else
        dif = mu_diff;
    end
end

%% EM single imputation of missing values (conditional means)
if imputation ==true
    Yimp = NAimputeConditionalmean(Y, mus, sigs);
else
    Yimp=[];
end


out.loc = mus;
out.cov = sigs;
out.iter = iter;
out.Yimp=Yimp;

end

%FScategory:MULT-MissingData
