function out = mdEM(Y, varargin)
% mdEM  EM algorithm for data with missing values (no trimming).
%
%<a href="matlab: docsearchFS('mdEM')">Link to the help function</a>
%
% Required input arguments:
%
% Y :           Input data. Matrix. n x p data matrix; n observations and p
%               variables possibly with missing values (NaN's). Rows of Y
%               represent observations, and columns represent variables.
%               Data Types - single | double
%
%  Optional input arguments:
%
%       mus     : intial mean. p x 1 vector or empty double.
%                 Initial  mean vector. If empty (default), column nanmeans are used.
%                 Example - 'mus',[]
%                 Data Types - single | double
%           sigs:  initial covariance matrix.
%                  p x p matrix or empty double.
%                  Initial p x p covariance matrix.
%                  If empty, uses nan-cov
%                 Example - 'sigs',eye(p)
%                 Data Types - single | double
%        maxiter:  maximum number of iterations. Positive integer.
%                  The default value is 100
%                 Example - 'maxiter',50
%                 Data Types - single | double
%           tol :  tolerance for convergence. Positive real number.
%                  The default value of the tolerance is 1e-5
%                 Example - 'tol',1e-10
%                 Data Types - single | double
%   tol_sigma  :   Use tolerance for both mu sigs. Boolean .
%                  If true use both mu and sigma diffs (default true)
%                 Example - 'tol_sigma',false
%                 Data Types - logical
%   imputation :  Also give the matrix of inputed values. Boolean.
%                 if true structure out also contains the matrix of inputed values.
%                 The default value of imputation is false.
%                 Example - 'imputation',true
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
%                   (depending on input option imputation)
%
%  
% See also: mdTEM, mdImputeCondMean.m, mdPartialMD.m, mdPartialMD2full
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
%<a href="matlab: docsearchFS('mdEM')">Link to the help page for this function</a>
%edit 
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Call to emNA with all the default options.
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
    out=mdEM(Y);
    % Show true means and inputed means
    scatter(out.mus,muTrue)
    xlabel('Imputed means')
    ylabel('True means')
%}

%{
    %% Example of use of option imputation.
    p = 5;                % number of variables
    n = 10;             % number of observations
    rho = 0.9;            % target pairwise correlation (0<rho<1)
    
    % Covariance matrix (unit variances)
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);      % upper-triangular such that Sigma = R'*R
    % Generate samples ~ N(0,Sigma)
    Yfull = randn(n,p) * R;   % Strong positive correlation between the vars
    missRate = 0.25;     % MCAR missing probability per entry
    missMask = rand(n,p) < missRate;
    Y=Yfull;
    Y(missMask) = NaN;
    
    out=mdEM(Y,'imputation',true);
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
mus=[];
sigs=[];
maxiter=100;
tol=1e-5;
tol_sigma=true;
imputation=false;


if nargin>1
 options=struct('mus',mus,'sigs',sigs,'maxiter',maxiter,'tol',tol, ...
     'tol_sigma',tol_sigma,'imputation',false);

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
    % Initialization with means and initial values based on means ??
    % Replace missing values with means then compute sample cov
    % on the matrix which replaces missing with means of the respective variable
    X0 = Y;
    for j = 1:p
        miss = isnan(X0(:,j));
        X0(miss,j) = mus(j);
    end
    sigs = cov(X0, 1); % use normalization by n (1) to be consistent with sums, acceptable init
end

dif = Inf;
iter = 0;
w = ones(n,1);

while (dif > tol) && (iter < maxiter)
    iter = iter + 1;
    mus_old = mus;
    sigs_old = sigs;

    % E-step with equal weights (no trimming)
    [T1, T2] = aux.NAcompute_expected_stats(Y, mus, sigs, w);

    % M-step
    [mus, sigs] = aux.NAmaximization_step(T1, T2, w);

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

out=struct;
out.loc = mus;
out.cov = sigs;
out.iter = iter;
out.Yimp=Yimp;

end

