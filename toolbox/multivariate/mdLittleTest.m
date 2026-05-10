function out = mdLittleTest(Y, varargin)
%mdLittleTest Little's test for Missing Completely At Random (MCAR).
%
%<a href="matlab: docsearchFS('mdLittleTest')">Link to the help function</a>
%
%  Little's test assesses the null hypothesis that the missing-data
%  mechanism is Missing Completely At Random (MCAR). The test is based on
%  pattern-specific mean deviations from the global maximum likelihood
%  estimate of the mean vector, using the corresponding submatrices of the
%  global covariance matrix.
%
%  Required input arguments:
%
%    Y :           Input data. Matrix. n x p data matrix; n observations and
%                  p variables possibly containing missing values (NaN's).
%                  Rows of Y represent observations, and columns represent
%                  variables.
%                  Data Types - single | double
%
%  Optional input arguments:
%
%    emOut :       Structure containing EM estimates.
%                  If supplied, it must contain fields:
%                  emOut.mu   : p x 1 estimated mean vector
%                  emOut.Sigma: p x p estimated covariance matrix
%                  If empty, function mdEM is called internally.
%                  Default is [].
%                  Example - 'emOut',outEM
%                  Data Types - struct
%
%    maxiter :     Maximum number of iterations for mdEM.
%                  Used only if option 'emOut' is empty.
%                  Default is 200.
%                  Example - 'maxiter',100
%                  Data Types - single | double
%
%    tol :         Convergence tolerance for mdEM.
%                  Used only if option 'emOut' is empty.
%                  Default is 1e-7.
%                  Example - 'tol',1e-6
%                  Data Types - single | double
%
%    msg :         Display messages from mdEM.
%                  Used only if option 'emOut' is empty.
%                  Default is false.
%                  Example - 'msg',true
%                  Data Types - logical
%
%
%  Output:
%
%    out :         Structure containing the following fields:
%
%      out.stat        = Little's test statistic.
%      out.df          = Degrees of freedom of the chi-square reference
%                        distribution.
%      out.pvalue      = p-value of the test.
%      out.mu          = Global MLE of the mean vector.
%      out.Sigma       = Global MLE of the covariance matrix.
%      out.patterns    = R x p logical matrix. Each row is a distinct
%                        missingness pattern; true means observed entry.
%      out.npatterns   = Number of distinct missingness patterns.
%      out.patternInfo = Table with one row for each informative pattern.
%
%  More About:
%
%  Let r=1,...,R index the distinct missingness patterns. For pattern r,
%  let n_r be the number of units following that pattern, let p_r be the
%  number of observed variables in that pattern, and let \bar{y}_r be the
%  sample mean vector computed using the observed variables only. Let
%  \hat{\mu} and \hat{\Sigma} be the global maximum likelihood estimates
%  under multivariate normality, typically obtained by EM. Denote by
%  \hat{\mu}_r and \hat{\Sigma}_r the subvector and submatrix corresponding
%  to the variables observed in pattern r. Little's statistic is
%
%  \[
%    T = sum_r n_r * (\bar{y}_r - \hat{\mu}_r)' * \hat{\Sigma}_r^(-1) ...
%                * (\bar{y}_r - \hat{\mu}_r).
%  \]
%  Under the null hypothesis of MCAR, T is asymptotically distributed as a
%  chi-square random variable with degrees of freedom
%
%    df = sum_r p_r - p
%
%  where p is the total number of variables.
%
%  Rows with all variables missing do not contribute to the test statistic.
%
%  References:
%
%  Little, R. J. A. (1988), "A Test of Missing Completely at Random for
%  Multivariate Data with Missing Values", Journal of the American
%  Statistical Association, 83, pp. 1198-1202.
%
%  See also: mdEM, mdPartialMD
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdLittleTest')">Link to the help page for this function</a>
%

%{
    %% Example 1: Little's MCAR test with default options.
    % Generate a data matrix with missing values and run the test using the
    % internal EM estimates.

    rng(1);
    Y = randn(100,3);
    Y(rand(100,3)<0.15) = NaN;

    out = mdLittleTest(Y);
    disp(out)
%}

%{
    %% Example 2: Supply EM estimates externally.
    % First compute the EM estimates using mdEM, then pass them to
    % mdLittleTest.

    rng(2);
    Y = randn(150,4);
    Y(rand(150,4)<0.20) = NaN;

    outEM = mdEM(Y);
    out = mdLittleTest(Y,'emOut',outEM);

    disp(out.stat)
    disp(out.pvalue)
%}


%{
    %% Example 3: Inspect missingness patterns.
    % The output contains the distinct observed-data patterns and a table
    % with information about the informative patterns.

    rng(4);
    Y = randn(120,4);
    Y(1:30,1) = NaN;
    Y(31:60,2) = NaN;
    Y(61:90,[3 4]) = NaN;

    out = mdLittleTest(Y);

    disp(out.patterns)
    disp(out.patternInfo)
%}

%{
    % Example 4: Data with rows completely missing.
    % Rows with all variables missing are ignored in the computation of
    % Little's statistic.
    rng(5);
    Y = randn(80,3);
    Y(rand(80,3)<0.15) = NaN;
    Y(1:5,:) = NaN;

    out = mdLittleTest(Y);

    disp(out.npatterns)
    disp(out.stat)
%}


%% Beginning of code

% Input parameters checking
if nargin < 1
    error('FSDA:mdLittleTest:TooFewInputs', ...
        'At least one input argument is required.');
end

if ~ismatrix(Y) || ~isnumeric(Y)
    error('FSDA:mdLittleTest:WrongInput', ...
        'Input argument Y must be a numeric matrix.');
end

[p] = size(Y,2);

% Default options
maxiter   = 200;
tol       = 1e-7;
msg       = false;
emOut=[];

% Optional arguments
if ~isempty(varargin)

    options=struct('maxiter',maxiter,'tol',tol,'msg',msg,'emOut',emOut);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSMmmd:mdLittleTest','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        aux.chkoptions(options,UserOptions)
    end

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    maxiter=options.maxiter;

    tol=options.tol;
    msg=options.msg;
    emOut=options.emOut;



    % Basic checks on options
    if ~isscalar(maxiter) || maxiter <= 0
        error('FSDA:mdLittleTest:WrongMaxiter', ...
            'Option ''maxiter'' must be a positive scalar.');
    end

    if ~isscalar(tol) || tol <= 0
        error('FSDA:mdLittleTest:WrongTol', ...
            'Option ''tol'' must be a positive scalar.');
    end

    if ~islogical(msg) && ~(isnumeric(msg) && isscalar(msg))
        error('FSDA:mdLittleTest:WrongMsg', ...
            'Option ''msg'' must be a logical scalar.');
    end

end

if isempty(emOut)
    % Obtain global MLE estimates
    outEM = mdEM(Y, 'maxiter', maxiter, 'tol', tol);
else
    outEM=emOut;
end

muhat = outEM.loc;
Sigmahat = outEM.cov;

% Distinct missingness patterns
obsMatrix = ~isnan(Y);
[Patterns, ~, idxPatterns] = unique(obsMatrix, 'rows', 'stable');
% boo=sum(Patterns,2)>0;
% Patterns=Patterns(boo,:);

nPatterns = size(Patterns, 1);

% Initialize accumulators
stat = 0;
sumpr = 0;

patternID = zeros(nPatterns,1);
nrVec     = patternID;
prVec     = patternID;
termVec   = patternID;
condVec   = patternID;


for r = 1:nPatterns
    idxr = (idxPatterns == r);
    nr = sum(idxr);
    obsr = Patterns(r,:);
    pr = sum(obsr);

    patternID(r,1) = r;
    nrVec(r,1)     = nr;
    
    % Completely missing rows carry no information
    if pr == 0
        continue
    end


    Yr = Y(idxr, obsr);
    ybarr = mean(Yr, 1, 'omitnan')';

    mur = muhat(obsr);
    Sigmar = Sigmahat(obsr, obsr);

    % Numerical stabilization
    Sigmar = (Sigmar + Sigmar')/2;

    % Check positive definiteness
    [~,pdflag] = chol(Sigmar);
    if pdflag ~= 0
        error('FSDA:mdLittleTest:NonPDSubmatrix', ...
            ['The covariance submatrix for pattern %d is not positive ' ...
            'definite.'], r);
    end
    termr=nr*mahalFS(ybarr',mur',Sigmar);

    stat = stat + termr;
    sumpr = sumpr + pr;

    prVec(r,1)     = pr;
    termVec(r,1)   = termr;
    condVec(r,1)   = cond(Sigmar);
end

df = sumpr - p;

if df <= 0
    pvalue = NaN;
    warning('FSDA:mdLittleTest:NonPositiveDf', ...
        ['The computed degrees of freedom are not positive. ' ...
        'The p-value is returned as NaN.']);
else
    pvalue = 1 - chi2cdf(stat, df);
end


patternInfo = table( nrVec, prVec, termVec, condVec, ...
    'VariableNames', {'nPattern','pObserved','Contribution','CondSigma'}, ...
    'RowNames',string(patternID));

plots=true;
if plots==true
    scatter(patternInfo,"pObserved","Contribution")

end
% Store output
out = struct;
out.stat        = stat;
out.df          = df;
out.pvalue      = pvalue;
out.mu          = muhat;
out.Sigma       = Sigmahat;
out.patterns    = Patterns;
out.npatterns   = nPatterns;
out.patternInfo = patternInfo;

end