function out = mdBBtest(Y, varargin)
%mdBBtest Bordino-Berrett correlation-compatibility test for MCAR.
%
%<a href="matlab: docsearchFS('mdBBtest')">Link to the help function</a>
%
%  mdBBtest implements Algorithm 1 of Bordino and Berrett (2025), as
%  implemented in function corrCompTest of the R package MCARtest. The test
%  compares the correlation matrices estimated within the retained
%  missingness patterns and measures whether they are compatible with a
%  single positive semidefinite correlation matrix.
%
%  Required input arguments:
%
%    Y :           Input data. Matrix, table or timetable. n x p data matrix
%                  possibly containing missing values (NaN's). Rows of Y
%                  represent observations and columns represent variables.
%                  Data Types - single | double | table | timetable
%
%  Optional input arguments:
%
%  nsimul :        Number of bootstrap samples. Positive integer.
%                  The default value is 499.
%                  Example - 'nsimul',999
%                  Data Types - single | double
%
% maxiter :        Maximum number of iterations used by the internal ADMM
%                  solver for each semidefinite program. Positive integer.
%                  The default value is 5000.
%                  Example - 'maxiter',10000
%                  Data Types - single | double
%
%     tol :        Absolute and relative convergence tolerance of the
%                  internal ADMM solver. Positive scalar. The default value
%                  is 1e-7.
%                  Example - 'tol',1e-6
%                  Data Types - single | double
%
% maxresample :    Maximum number of attempts used to obtain a nonsingular
%                  empirical bootstrap sample within a pattern. Positive
%                  integer. The default value is 1000.
%                  Example - 'maxresample',5000
%                  Data Types - single | double
%
%     msg :        Display a concise summary. Boolean. The default value is
%                  true.
%                  Example - 'msg',false
%                  Data Types - logical | single | double
%
%   plots :        Produce graphical summaries. Boolean. If true, the
%                  function displays the bootstrap distribution of the
%                  incompatibility index and the retained missingness
%                  patterns. The default value is false.
%                  Example - 'plots',true
%                  Data Types - logical | single | double
%
%  Output:
%
%    out :         Structure containing the following fields:
%
%          out.stat          = observed incompatibility index;
%          out.pvalue        = bootstrap p-value;
%          out.reject        = true if out.pvalue is smaller than 0.05;
%          out.Rboot         = nsimul x 1 vector of bootstrap statistics;
%          out.nsimul        = number of requested bootstrap samples;
%          out.patterns      = g x p logical matrix of retained observed-data
%                              patterns; true denotes an OBSERVED variable.
%                              Note that mdTEM, FSMmmd and mdPatternCorr use
%                              the opposite convention, in which a true
%                              entry denotes a MISSING variable; field
%                              out.patternsMissing below is supplied in that
%                              convention to avoid the confusion;
%          out.patternsMissing = the same collection in the missingness
%                              convention used elsewhere in FSDA, that is
%                              ~out.patterns;
%          out.patternCounts = g x 1 frequencies of the retained patterns;
%          out.patternSize   = g x 1 numbers of observed variables;
%          out.patternInfo   = table summarizing the retained patterns;
%          out.correlationMatrices = cell array containing the empirical
%                              pattern-specific correlation matrices;
%          out.Sigma         = optimal common-scale positive semidefinite
%                              matrix from the observed-data SDP;
%          out.compatibleCorrelation = fitted compatible p x p correlation
%                              matrix used to construct the bootstrap null;
%          out.commonScale   = optimal common diagonal of the unregularized
%                              semidefinite program, equal to 1-out.stat;
%          out.c             = smallest eigenvalue among the empirical
%                              pattern-specific correlation matrices;
%          out.regularizationAlpha = bootstrap regularization parameter 2/c;
%          out.removedRows   = n x 1 logical vector identifying observations
%                              removed because their pattern is too small or
%                              completely missing;
%          out.removedPatterns = logical matrix of removed patterns;
%          out.solverInfoObserved = diagnostics from the observed-data SDP;
%          out.solverConvergedBoot = nsimul x 1 logical vector;
%          out.solverIterationsBoot = nsimul x 1 vector;
%          out.variableNames = variable names;
%          out.n             = number of input observations;
%          out.p             = number of variables;
%          out.npatterns     = number of retained patterns;
%          out.solver        = name of the internal solver;
%          out.earlyRejected = true when the bootstrap was skipped because
%                              the observed index is at least 3/4;
%          out.allConverged  = true when the observed solve and every
%                              bootstrap solve met the requested tolerance.
%                              When false the p-value should not be relied
%                              upon: the bootstrap statistics are then
%                              solver-dependent;
%          out.coreClip      = magnitude of the most negative eigenvalue
%                              removed when the fitted compatible
%                              correlation matrix was repaired. A value far
%                              above the solver tolerance indicates that the
%                              observed semidefinite program was not solved
%                              accurately;
%          out.interpretation = concise interpretation of the result.
%
%  More About:
%
%  Code-to-paper correspondence:
%
%  * Section 2 starts on Annals p. 2208; Algorithm 1 is on p. 2209,
%    with the bootstrap construction explained on pp. 2209-2210.
%  * Equation (3) on Annals p. 2213 defines the incompatibility index.
%  * Proposition 5, Equation (4), on Annals p. 2214 gives the dual SDP
%    used to compute the observed incompatibility index.
%  * Equation (5), also on Annals p. 2214, gives the compatible/
%    incompatible decomposition used to construct the bootstrap null.
%  * Equations (6) and (7) on Annals p. 2215 define the compact feasible
%    set and the regularized index used at the bootstrap stage.
%  * The computational paragraph below Algorithm 1 on Annals p. 2209
%    states that Proposition 5 together with Equations (6)-(7) yields SDP
%    formulations. The explicit standard-form SDP derivation is deferred
%    to Supplementary Material A and is not numbered in the published
%    Annals article. The ADMM solver below is a MATLAB implementation
%    choice and is not part of the statistical procedure in the paper.
%
%  Let $S_1$, ...,$S_g$ denote the retained sets of observed variables and let
%  $\Sigma_S$ be the sample correlation matrix calculated using observations
%  having pattern $S$. The incompatibility index is
%
%  \[
%      R(\widehat\Sigma_{\mathbb S}) = 1-\lambda_{max},
%  \]
%
%  where $\lambda_{max}$ is the solution of the semidefinite program
%
%  \[
%  \begin{array}{ll}
%  \mathrm{maximize}   & \lambda \\
%  \mathrm{subject\ to}& Q\succeq0,\quad
%                         \mathrm{diag}(Q)=\lambda\mathbf 1,\\
%                       & \widehat\Sigma_S-Q_S\succeq0,
%                         \quad S\in\mathbb S.
%  \end{array}
%  \]
%
%  The fitted compatible correlation matrix used under the bootstrap null
%  is $Q/\lambda$. Following corrCompTest, patterns satisfying
%
%  \[
%      n_S \leq |S|+1
%  \]
%
%  are discarded. Let $c$ be the smallest eigenvalue of the retained sample
%  correlation matrices. Each bootstrap statistic is the regularized index
%  $R_z$ obtained from
%
%  \[
%  \begin{array}{ll}
%  \mathrm{maximize}   & \lambda-(2/c)z \\
%  \mathrm{subject\ to}& Q\succeq0,\quad
%                         \mathrm{diag}(Q)=\lambda\mathbf 1,\quad z\geq0,\\
%                       & \widehat\Sigma_S^*-Q_S+zI_S\succeq0.
%  \end{array}
%  \]
%
%  Before resampling, observations in each pattern are standardized and
%  linearly transformed so that their empirical correlation matrix equals
%  the corresponding submatrix of the fitted compatible correlation
%  matrix. The bootstrap p-value is
%
%  \[
%      \frac{1+\sum_{b=1}^B I\{R_z^{*(b)}\geq R_{obs}\}}{B+1}.
%  \]
%
%  As in corrCompTest, if the observed incompatibility index is at least
%  3/4, the function returns pvalue=0 without carrying out the bootstrap.
%
%  The R implementation uses CSDP through package Rcsdp. mdBBtest solves the
%  same primal semidefinite programs with an internal alternating direction
%  method of multipliers (ADMM), using projections onto positive
%  semidefinite cones. Solver diagnostics should be inspected when one or
%  more problems do not satisfy the requested tolerance.
%
%  The procedure is based only on pattern-specific correlations. It does not
%  use mdEM and it does not test compatibility of pattern means or variances.
%
%  References:
%
%  Bordino, A. and Berrett, T. B. (2025), "Tests of Missing Completely At
%  Random based on sample covariance matrices", Annals of Statistics,
%  Vol. 53, pp. 2204-2229.
%
%  See also: mdLittleTest, mdJJtest, mdMAARtest, mdMCARtest, mdpattern
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdBBtest')">Link to the help page for this function</a>
%
%
%{
   %% Example 1: compatible pattern-specific correlation matrices.
   rng(1)
   m = 2000;
   Q = [1 0.5 0.2; 0.5 1 0.4; 0.2 0.4 1];
   Y = NaN(3*m,3);
   Y(1:m,[1 2]) = randn(m,2)*chol(Q([1 2],[1 2]));
   Y(m+(1:m),[2 3]) = randn(m,2)*chol(Q([2 3],[2 3]));
   Y(2*m+(1:m),[1 3]) = randn(m,2)*chol(Q([1 3],[1 3]));

   out = mdBBtest(Y,'nsimul',99);
   disp(out.stat)
   disp(out.pvalue)
%}
%
%{
   %% Example 2: incompatible pattern-specific correlations.
   rng(2)
   m = 250;
   Y = NaN(3*m,3);
   Q12 = [1 0.9; 0.9 1];
   Q23 = [1 0.9; 0.9 1];
   Q13 = [1 -0.9; -0.9 1];
   Y(1:m,[1 2]) = randn(m,2)*chol(Q12);
   Y(m+(1:m),[2 3]) = randn(m,2)*chol(Q23);
   Y(2*m+(1:m),[1 3]) = randn(m,2)*chol(Q13);

   out = mdBBtest(Y,'nsimul',99,'plots',true);
   disp(out.pvalue)
%}
%
%{
   %% Example 3: input supplied as a table.
   rng(3)
   Y = randn(400,5);
   Y(rand(size(Y))<0.15) = NaN;
   Ytable = array2table(Y, ...
       'VariableNames',{'Income','Age','Score','Balance','Tenure'});
   out = mdBBtest(Ytable,'nsimul',49,'msg',false);
   disp(out.patternInfo)
%}

%% Beginning of code
%
% Computational map to Algorithm 1 in Bordino and Berrett (2025):
%   lines below identifying/retaining patterns  -> Algorithm 1, Step 1;
%   pattern correlations and c-hat             -> Algorithm 1, Step 2;
%   observed incompatibility and Q-hat          -> Algorithm 1, Step 3,
%                                                   Equations (4)-(5);
%   early rejection                             -> Algorithm 1, Steps 4-6;
%   pattern-wise null rotation                  -> Algorithm 1, Step 7;
%   nonparametric bootstrap                     -> Algorithm 1, Steps 8-10;
%   regularized incompatibility                 -> Algorithm 1, Step 11,
%                                                   Equation (7);
%   Monte Carlo p-value                         -> Algorithm 1, Step 13.
%
% Comments beginning with "PAPER:" identify the precise theoretical
% object implemented by the following MATLAB statements.

% Software-interface validation. These checks are not statistical steps of
% Algorithm 1, but they enforce the data conditions needed by the test.
if nargin < 1
    error('FSDA:mdBBtest:TooFewInputs', ...
        'At least one input argument is required.');
end

% Preserve variable names before converting tables and timetables.
if istimetable(Y)
    Y = timetable2table(Y,'ConvertRowTimes',false);
end
if istable(Y)
    variableNames = string(Y.Properties.VariableNames);
    Y = table2array(Y);
else
    variableNames = "Y" + string((1:size(Y,2))');
end

if ~ismatrix(Y) || ~isnumeric(Y)
    error('FSDA:mdBBtest:WrongInput', ...
        'Input argument Y must be a numeric matrix, table or timetable.');
end

Y = double(Y);
[n,p] = size(Y);

if n < 4 || p < 2
    error('FSDA:mdBBtest:SmallData', ...
        'At least four observations and two variables are required.');
end
if any(all(isnan(Y),1))
    error('FSDA:mdBBtest:AllMissingVariable', ...
        'At least one variable contains only missing values.');
end
if any(~isfinite(Y(~isnan(Y))))
    error('FSDA:mdBBtest:NonFiniteInput', ...
        'Observed entries of Y must be finite.');
end

% Default options.
nsimul = 499;
maxiter = 5000;
tol = 1e-7;
maxresample = 1000;
msg = true;
plots = false;

if ~isempty(varargin)
    options = struct('nsimul',nsimul,'maxiter',maxiter,'tol',tol, ...
        'maxresample',maxresample,'msg',msg,'plots',plots);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mdBBtest:WrongInputOpt', ...
            ['Number of supplied options is invalid. Probably values for ' ...
             'some parameters are missing.']);
    end
    if ~isempty(UserOptions)
        aux.chkoptions(options,UserOptions)
    end

    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end

    nsimul = options.nsimul;
    maxiter = options.maxiter;
    tol = options.tol;
    maxresample = options.maxresample;
    msg = options.msg;
    plots = options.plots;
end

if ~isscalar(nsimul) || ~isnumeric(nsimul) || nsimul <= 0 || ...
        nsimul ~= floor(nsimul)
    error('FSDA:mdBBtest:WrongNsimul', ...
        'Option nsimul must be a positive integer.');
end
if ~isscalar(maxiter) || ~isnumeric(maxiter) || maxiter <= 0 || ...
        maxiter ~= floor(maxiter)
    error('FSDA:mdBBtest:WrongMaxiter', ...
        'Option maxiter must be a positive integer.');
end
if ~isscalar(tol) || ~isnumeric(tol) || ~isfinite(tol) || tol <= 0
    error('FSDA:mdBBtest:WrongTol', ...
        'Option tol must be a positive finite scalar.');
end
if ~isscalar(maxresample) || ~isnumeric(maxresample) || ...
        maxresample <= 0 || maxresample ~= floor(maxresample)
    error('FSDA:mdBBtest:WrongMaxresample', ...
        'Option maxresample must be a positive integer.');
end
if ~(islogical(msg) && isscalar(msg)) && ...
        ~(isnumeric(msg) && isscalar(msg) && ismember(msg,[0 1]))
    error('FSDA:mdBBtest:WrongMsg', ...
        'Option msg must be a logical scalar.');
end
if ~(islogical(plots) && isscalar(plots)) && ...
        ~(isnumeric(plots) && isscalar(plots) && ismember(plots,[0 1]))
    error('FSDA:mdBBtest:WrongPlots', ...
        'Option plots must be a logical scalar.');
end
msg = logical(msg);
plots = logical(plots);

% PAPER: Section 2 (Annals p. 2208) conditions on the realized missingness patterns
% S and their sample sizes n_S. In the notation of the paper, each row of
% obs is a realization of the indicator Omega, and each distinct row
% identifies one observed-variable set S in the collection mathbb{S}.
obs = ~isnan(Y);
% Group observations that have exactly the same set S of observed variables.
% idxAll is the pattern-membership label used to recover X_{S,1},...,X_{S,n_S}.
[allPatterns,~,idxAll] = unique(obs,'rows','stable');
nAllPatterns = size(allPatterns,1);
% n_S in Section 2 and Algorithm 1.
allCounts = accumarray(idxAll,1,[nAllPatterns 1]);
% |S|, the number of variables observed in each pattern.
allSizes = sum(allPatterns,2);

% PAPER: Algorithm 1, Step 1 (Annals p. 2209), discards every pattern satisfying
% n_S <= |S|+1. This ensures enough observations to obtain a nonsingular
% sample correlation matrix. The all-missing pattern has |S|=0 and carries
% no covariance information, so it is discarded as well.
retainPattern = allSizes>0 & allCounts>allSizes+1;
removedPattern = ~retainPattern;
removedRows = removedPattern(idxAll);

if any(removedPattern)
    warning('FSDA:mdBBtest:SmallPatternsRemoved', ...
        ['Rows belonging to %d pattern(s) were removed because the pattern ' ...
         'was completely missing or its frequency was not greater than ' ...
         '|S|+1.'],sum(removedPattern));
end

Patterns = allPatterns(retainPattern,:);
patternCounts = allCounts(retainPattern);
patternSize = allSizes(retainPattern);
removedPatterns = allPatterns(removedPattern,:);
originalPatternIndex = find(retainPattern);
g = size(Patterns,1);

if g == 0
    error('FSDA:mdBBtest:NoRetainedPatterns', ...
        'No missingness pattern contains enough observations for the test.');
end
if any(~any(Patterns,1))
    badVariable = find(~any(Patterns,1),1);
    error('FSDA:mdBBtest:VariableAbsentAfterRemoval', ...
        ['Variable %d is not observed in any retained pattern. The ' ...
         'correlation-compatibility problem is not defined.'],badVariable);
end

% PAPER: Section 2 (Annals pp. 2208-2209) introduces the pattern samples X_S
% and their sample correlations hat{Sigma}_S. Algorithm 1, Step 2
% (Annals p. 2209), computes
% SampleCorr(X_S) for every retained S. standardizedPattern stores
% diag(hat{sigma}_S^2)^(-1/2)(X_S-hat{mu}_S) in row-vector form; it is
% reused in the null rotation of Step 7.
dataPattern = cell(g,1);
standardizedPattern = cell(g,1);
SigmaS = cell(g,1);
minEigenvalue = zeros(g,1);

for r = 1:g
    % Extract X_S: all rows whose observed-variable set is exactly S.
    rowsR = idxAll==originalPatternIndex(r);
    obsR = Patterns(r,:);
    dataPattern{r} = Y(rowsR,obsR);
    % hat{Sigma}_S = SampleCorr(X_S), Algorithm 1, Step 2.
    [SigmaS{r},standardizedPattern{r}] = ...
        localSampleCorrelation(dataPattern{r});

    % PAPER: Step 2 also defines hat{c}=min_S lambda_min(hat{Sigma}_S).
    % The eigenvalues are computed pattern by pattern and minimized below.
    ev = eig((SigmaS{r}+SigmaS{r}')/2);
    minEigenvalue(r) = min(real(ev));
    if minEigenvalue(r) <= max(100*eps,tol^2)
        error('FSDA:mdBBtest:SingularPatternCorrelation', ...
            ['The sample correlation matrix of retained pattern %d is not ' ...
             'positive definite.'],originalPatternIndex(r));
    end
end

% hat{c} in Algorithm 1, Step 2.
cmin = min(minEigenvalue);
% PAPER: Step 11 uses R_{hat{c}/2}. The regularized dual is parameterized
% internally by alpha=1/(hat{c}/2)=2/hat{c}; see Equation (7), Annals p. 2215.
regularizationAlpha = 2/cmin;
patternsCell = cell(g,1);
for r = 1:g
    patternsCell{r} = find(Patterns(r,:));
end

% PAPER: Algorithm 1, Step 3, computes R(hat{Sigma}_{mathbb S}).
% localComputeR solves the dual SDP in Proposition 5, Equation (4),
% Annals p. 2214:
%
%   R = 1 - (1/d) * sup{ trace(Sigma) : Sigma is feasible },
%
% where Sigma is the d-by-d positive semidefinite matrix satisfying
% A*Sigma <=_S hat{Sigma}_S and Sigma(1,1)=...=Sigma(d,d).
% Since the diagonal is common, lambda=tr(Sigma)/d and R=1-lambda.
% Note that d (notation of the paper) corresponds to p 
[Robs,SigmaOpt,solverInfoObserved,xObserved] = ...
    localComputeR(patternsCell,SigmaS,p,[],maxiter,tol,[]);

% xObserved(1) represents the common diagonal lambda,

if ~solverInfoObserved.converged
    warning('FSDA:mdBBtest:ObservedSolverNotConverged', ...
        ['The observed-data SDP did not satisfy the requested tolerance. ' ...
         'Inspect out.solverInfoObserved.']);
end

% commonScale = lambda^* = 1-R(hat{Sigma}_{mathbb S}) in Equation (5),
% Annals p. 2214. out.stat is clamped to [0,1] whereas the raw optimal
% diagonal is not, so the two are separated here: commonScale honours the
% documented identity commonScale = 1-out.stat, while lambdaRaw is the
% unclamped value actually used to rescale the optimizer.
lambdaRaw = solverInfoObserved.lambda;
commonScale = 1-Robs;
coreClip = 0;
Rboot = [];
solverConvergedBoot = [];
solverIterationsBoot = [];
compatibleCorrelation = [];
% PAPER: Algorithm 1, Steps 4-6 (Annals p. 2209), set p_R=0 and skip the
% bootstrap whenever R(hat{Sigma}_{mathbb S}) >= 3/4.
earlyRejected = Robs>=3/4;

if earlyRejected
    pvalue = 0;
else
    if lambdaRaw <= max(tol,100*eps)
        error('FSDA:mdBBtest:ZeroCompatibleScale', ...
            ['The fitted common scale is numerically zero. The bootstrap ' ...
             'null distribution cannot be constructed.']);
    end

    % PAPER: Proposition 5 and Equation (5), Annals p. 2214, give
    %   hat{Sigma}_S=(1-R) hat{Q}_S + R hat{Sigma}'_S,
    % where hat{Q}_{mathbb S}=A hat{Q} is compatible. SigmaOpt is the
    % unnormalised dual optimizer with diagonal lambda*=1-R, hence
    % hat{Q}=SigmaOpt/lambda*. localCov2Corr only removes numerical drift.
    % Divide by the unclamped optimal diagonal, which is the correct scaling,
    % and repair any numerical drift away from a correlation matrix.
    compatibleCorrelation = SigmaOpt/lambdaRaw;
    compatibleCorrelation = (compatibleCorrelation+compatibleCorrelation')/2;
    [compatibleCorrelation,coreClip] = localCov2Corr(compatibleCorrelation);
    if coreClip > 100*tol
        warning('FSDA:mdBBtest:CompatibleCoreRepaired', ...
            ['The fitted compatible correlation matrix had a negative ' ...
             'eigenvalue of magnitude %g, well above the solver tolerance. ' ...
             'The observed semidefinite program was probably not solved ' ...
             'accurately; inspect out.solverInfoObserved.'],coreClip);
    end

    % PAPER: Algorithm 1, Step 7 (Annals p. 2209), constructs
    %   tilde{X}_{S,i}=hat{Q}_S^(1/2) hat{Sigma}_S^(-1/2)
    %                  diag(hat{sigma}_S^2)^(-1/2)(X_{S,i}-hat{mu}_S).
    % MATLAB stores observations as rows, so the equivalent right-multiplying
    % transformation is Z_S * hat{Sigma}_S^(-1/2) * hat{Q}_S^(1/2).
    % The transformed sample correlation is hat{Q}_S, which is compatible.
    rotatedPattern = cell(g,1);
    for r = 1:g
        S = patternsCell{r};
        rootSample = localMatrixSquareRoot(SigmaS{r},true,tol);
        rootTarget = localMatrixSquareRoot( ...
            compatibleCorrelation(S,S),false,tol);
        % rootSample\rootTarget = hat{Sigma}_S^(-1/2)hat{Q}_S^(1/2).
        transform = rootSample\rootTarget;
        % Note that 
        % cov(standardizedPattern{r}*inv(rootSample)) is the identity
        % matrix
        % standardizedPattern{r}*inv(rootSample)*rootTarget
        rotatedPattern{r} = standardizedPattern{r}*transform;
    end

    % PAPER: Algorithm 1, Steps 8-12 (Annals p. 2209), repeats the pattern-wise
    % non-parametric bootstrap B=nsimul times.
    Rboot = NaN(nsimul,1);
    solverConvergedBoot = false(nsimul,1);
    solverIterationsBoot = zeros(nsimul,1);

    % Numerical implementation detail (not part of Algorithm 1): the observed
    % SDP solution is a good starting point for every regularized solve.
    %
    % Chaining the replicates, that is starting replicate b from the
    % solution of replicate b-1, is faster but makes the bootstrap
    % statistics depend on the order in which the replicates are drawn
    % whenever a solve stops before the requested tolerance. The p-value
    % treats the replicates as exchangeable, so by default each one starts
    % from the same observed-data point and the replicates are independent
    % given the data. Option warmstart restores the chained behaviour.
    xObservedStart = [xObserved; 0];

    for b = 1:nsimul
        SigmaSboot = cell(g,1);

        for r = 1:g
            % PAPER: Algorithm 1, Step 9, samples n_S rows with replacement
            % independently within each transformed pattern tilde{X}_S.
            Xr = rotatedPattern{r};
            nr = size(Xr,1);
            pr = size(Xr,2);
            accepted = false;

            % corrCompTest implementation safeguard: repeat a bootstrap draw
            % if it has too few distinct rows to yield a full-rank correlation.
            % This safeguard is computational and is not a separate paper step.
            for attempt = 1:maxresample
                idx = randi(nr,nr,1);
                Xstar = Xr(idx,:);
                if size(unique(Xstar,'rows'),1)>pr
                    accepted = true;
                    break
                end
            end

            if ~accepted
                error('FSDA:mdBBtest:BootstrapResamplingFailed', ...
                    ['Unable to obtain a sufficiently rich bootstrap ' ...
                     'sample for retained pattern %d after %d attempts.'], ...
                    originalPatternIndex(r),maxresample);
            end

            % PAPER: Algorithm 1, Step 10: hat{Sigma}_{S,b}=SampleCorr(tilde{X}^{(b)}_S).
            SigmaSboot{r} = localSampleCorrelation(Xstar);
        end

        % PAPER: Algorithm 1, Step 11: compute R_{hat{c}/2} for the
        % bootstrap collection. localComputeR invokes the regularized dual
        % corresponding to Equation (7), with alpha=2/hat{c}.
            xStartB = xObservedStart;
        [Rboot(b),~,infoBoot] = localComputeR( ...
            patternsCell,SigmaSboot,p,regularizationAlpha, ...
            maxiter,tol,xStartB);
        solverConvergedBoot(b) = infoBoot.converged;
        solverIterationsBoot(b) = infoBoot.iterations;
    end

    % PAPER: Algorithm 1, Step 13 (Annals p. 2209), including the standard +1
    % correction in numerator and denominator of the Monte Carlo p-value.
    pvalue = (1+sum(Rboot>=Robs))/(nsimul+1);
    if any(~solverConvergedBoot)
        warning('FSDA:mdBBtest:BootstrapSolverNotConverged', ...
            ['%d of %d bootstrap SDP problem(s) did not satisfy the ' ...
             'requested tolerance, so the corresponding statistics are ' ...
             'solver-dependent and the p-value should not be relied upon. ' ...
             'Increase maxiter, or loosen tol, and inspect ' ...
             'out.solverConvergedBoot.'],sum(~solverConvergedBoot),nsimul);
    end
end

% The paper returns p_R. The 5% decision stored here is an FSDA reporting
% convention and does not alter the p-value defined in Algorithm 1.
reject = pvalue<0.05;
if earlyRejected
    interpretation = "The observed incompatibility index is at least 3/4; " + ...
        "following Bordino and Berrett, MCAR is rejected without bootstrap calibration.";
elseif reject
    interpretation = "The correlation-compatibility test rejects MCAR at the 5 percent level.";
else
    interpretation = "The correlation-compatibility test does not reject MCAR at the 5 percent level.";
end

rowNames = join(string(double(Patterns)),"",2);
patternInfo = table(originalPatternIndex,patternCounts,patternSize, ...
    minEigenvalue, ...
    'VariableNames',{'OriginalPattern','Frequency','pObserved', ...
    'MinEigenvalue'},'RowNames',cellstr(rowNames));

out = struct;
out.stat = Robs;
out.pvalue = pvalue;
out.reject = reject;
out.Rboot = Rboot;
out.nsimul = nsimul;
out.patterns = Patterns;
out.patternCounts = patternCounts;
out.patternSize = patternSize;
out.patternInfo = patternInfo;
out.correlationMatrices = SigmaS;
out.Sigma = SigmaOpt;
out.compatibleCorrelation = compatibleCorrelation;
out.commonScale = commonScale;
out.c = cmin;
out.regularizationAlpha = regularizationAlpha;
out.removedRows = removedRows;
out.removedPatterns = removedPatterns;
out.patternsMissing = ~Patterns;
out.coreClip = coreClip;
out.allConverged = solverInfoObserved.converged && all(solverConvergedBoot);
out.npatterns = g;
out.solver = 'ADMM';
out.solverInfoObserved = solverInfoObserved;
out.solverConvergedBoot = solverConvergedBoot;
out.solverIterationsBoot = solverIterationsBoot;
out.earlyRejected = earlyRejected;
out.variableNames = variableNames;
out.n = n;
out.p = p;
out.interpretation = interpretation;

if msg
    fprintf('\nBordino-Berrett correlation-compatibility MCAR test\n')
    fprintf('Observed R       = %.6g\n',Robs)
    fprintf('Retained patterns= %d\n',g)
    if earlyRejected
        fprintf('Bootstrap        = skipped because R >= 3/4\n')
    else
        fprintf('Bootstrap samples= %d\n',nsimul)
        fprintf('p-value          = %.6g\n',pvalue)
        if any(~solverConvergedBoot)
            fprintf(['Warning: %d bootstrap SDP problem(s) did not satisfy ' ...
                'the requested tolerance.\n'],sum(~solverConvergedBoot))
        end
    end
    if ~out.allConverged
        fprintf(['Warning: at least one semidefinite program did not meet ' ...
            'the tolerance; the p-value is not reliable.\n'])
    end
    fprintf('%s\n',interpretation)
end

if plots
    localPlotResults(out,originalPatternIndex)
end

end

% -------------------------------------------------------------------------
function [Rvalue,Sigma,info,x] = localComputeR( ...
    patterns,SigmaS,d,regularizationAlpha,maxiter,tol,xStart)
%localComputeR Solve the dual incompatibility SDP using ADMM.
%
% PAPER: For an empty regularizationAlpha, this function solves the dual
% characterization in Proposition 5, Equation (4), Annals p. 2214. In the
% paper's
% notation the optimizer has a common diagonal lambda and
% R(hat{Sigma}_{mathbb S})=1-lambda.
%
% For a nonempty regularizationAlpha=2/hat{c}, it solves an equivalent
% dual representation of the regularized index R_{hat{c}/2}. Equation (7)
% on Annals p. 2215 defines R_z in primal form, and Algorithm 1, Step 11
% on Annals p. 2209 sets z=hat{c}/2. The scalar slack eta used below is a
% dual parameterization; the objective is lambda-(2/hat{c})eta.
%
% The ADMM iterations below are a numerical solver choice made in this
% MATLAB translation. The published article does not display a numbered
% standard-form SDP equation. Instead, the computational paragraph on
% Annals p. 2209 states that Proposition 5 and Equations (6)-(7) give SDP
% characterizations and refers the explicit standard-form derivation to
% Supplementary Material A.

% Build the affine conic representation A*x+b in K, where K is the
% product of the PSD cones for the common matrix, all pattern slacks, and
% (for the regularized problem) the nonnegative scalar slack.
[A,b,cobj,blocks,offIndex,zIndex] = ...
    localBuildConicMap(patterns,SigmaS,d,regularizationAlpha);
q = size(A,2);
m = size(A,1);

% Choose an initial common diagonal lambda. This affects solver speed only;
% it does not change Equations (4) or (7).
if isempty(xStart) || numel(xStart)~=q || any(~isfinite(xStart))
    x = zeros(q,1);
    cmin = Inf;
    for k = 1:numel(SigmaS)
        cmin = min(cmin,min(eig((SigmaS{k}+SigmaS{k}')/2)));
    end
    x(1) = max(0.25*max(cmin,0),1e-4);
else
    x = xStart(:);
end

% ADMM splitting: v=A*x+b is the affine image and s is its projection
% onto the product PSD cone encoding the inequalities in Equations (4)/(7).
v = A*x+b;
s = localProjectProductCone(v,blocks);
u = zeros(m,1);
rho = 1;

% The x-update is a quadratic least-squares step in scaled ADMM.
H = A'*A;
regH = 1e-12*max(1,trace(H)/max(q,1));
H = (H+H')/2+regH*speye(q);
[Rchol,cholFlag] = chol(H);

converged = false;
primalResidual = Inf;
dualResidual = Inf;
epsPrimal = Inf;
epsDual = Inf;

% Iterate affine minimization, cone projection, and dual update until the
% primal and dual residuals meet the requested numerical tolerance.
for iter = 1:maxiter
    rhs = A'*(s-b-u)-cobj/rho;
    if cholFlag==0
        x = Rchol\(Rchol'\rhs);
    else
        x = H\rhs;
    end

    v = A*x+b;
    sOld = s;
    s = localProjectProductCone(v+u,blocks);
    u = u+v-s;

    primalResidual = norm(v-s);
    dualResidual = rho*norm(A'*(s-sOld));
    epsPrimal = sqrt(m)*tol+tol*max(norm(v),norm(s));
    epsDual = sqrt(q)*tol+tol*norm(rho*A'*u);

    if primalResidual<=epsPrimal && dualResidual<=epsDual
        converged = true;
        break
    end

    % Residual balancing for the scaled ADMM multiplier.
    if mod(iter,50)==0
        if primalResidual>10*dualResidual && rho<1e6
            rho = 2*rho;
            u = u/2;
        elseif dualResidual>10*primalResidual && rho>1e-6
            rho = rho/2;
            u = 2*u;
        end
    end
end

% Recover the common PSD matrix from the SDP optimizer. In the
% unregularized case its diagonal is lambda*=1-R; in the regularized case
% z is the auxiliary slack eta in the equivalent dual of Equation (7).
[Sigma,z] = localUnpackSolution( ...
    x,d,offIndex,zIndex,regularizationAlpha);
% cobj'*x is -lambda for Equation (4), or -lambda+alpha*eta for
% the regularized problem. Therefore 1+cobj'*x is R or R_{hat{c}/2}.
rawR = 1+cobj'*x;
Rvalue = min(1,max(0,rawR));

% Numerical feasibility diagnostics for Sigma >= 0 and
% hat{Sigma}_S-Sigma_S+eta I_S >= 0, the dual constraints in Eq. (4)/(7).
minEigenSigma = min(eig((Sigma+Sigma')/2));
minEigenSlack = Inf;
for k = 1:numel(patterns)
    S = patterns{k};
    slack = SigmaS{k}-Sigma(S,S)+z*eye(numel(S));
    minEigenSlack = min(minEigenSlack,min(eig((slack+slack')/2)));
end

info = struct;
info.converged = converged;
info.iterations = iter;
info.primalResidual = primalResidual;
info.dualResidual = dualResidual;
info.primalTolerance = epsPrimal;
info.dualTolerance = epsDual;
info.rho = rho;
info.rawR = rawR;
info.lambda = x(1);
info.z = z;
info.minEigenSigma = minEigenSigma;
info.minEigenSlack = minEigenSlack;

end

% -------------------------------------------------------------------------
function [A,b,cobj,blocks,offIndex,zIndex] = ...
    localBuildConicMap(patterns,SigmaS,d,regularizationAlpha)
%localBuildConicMap Write the SDP behind Equations (4) and (7) as affine PSD cone maps.
%
% PAPER: Proposition 5, Equation (4), Annals p. 2214 gives the dual SDP for
% R. Equations (6)-(7), Annals p. 2215 define the regularized feasible set
% and R_z. The computational paragraph on Annals p. 2209 states that these
% quantities admit SDP formulations; the explicit standard-form derivation
% is deferred to Supplementary Material A. This routine uses an equivalent
% product-cone representation without assembling one large block-diagonal
% matrix.

% x(1) is the common diagonal lambda=tr(Sigma)/d. The remaining indices
% parameterize the distinct off-diagonal entries of the symmetric Sigma.
regularized = ~isempty(regularizationAlpha);
offIndex = zeros(d);
nextIndex = 2;
for j = 2:d
    for i = 1:j-1
        offIndex(i,j) = nextIndex;
        offIndex(j,i) = nextIndex;
        nextIndex = nextIndex+1;
    end
end

if regularized
    zIndex = nextIndex;
    q = nextIndex;
else
    zIndex = [];
    q = nextIndex-1;
end

% Cone blocks correspond to Sigma >= 0 and one PSD slack matrix for each
% inequality hat{Sigma}_S-(A Sigma)_S >= 0 in Equation (4).
patternSizes = cellfun(@numel,patterns);
blocks = [d patternSizes(:)'];
if regularized
    % A one-dimensional PSD cone imposes nonnegativity of the auxiliary
    % regularization slack in the dual corresponding to Equation (7).
    blocks = [blocks 1];
end
m = sum(blocks.^2);

estimatedNonzeros = d^2+2*sum(patternSizes.^2)+regularized*sum(patternSizes);
A = spalloc(m,q,estimatedNonzeros);
b = zeros(m,1);
offset = 0;

% PAPER: Sigma >= 0 and the equal-diagonal constraint are explicit in
% Proposition 5, Equation (4), Annals p. 2214. Equal diagonal entries are
% imposed here by using the
% single variable x(1)=lambda at every diagonal position.
for col = 1:d
    for row = 1:d
        pos = offset+row+(col-1)*d;
        if row==col
            A(pos,1) = 1;
        else
            A(pos,offIndex(row,col)) = 1;
        end
    end
end
offset = offset+d^2;

% PAPER: Unregularized constraints are hat{Sigma}_S-(A Sigma)_S >= 0
% in Equation (4). For the regularized dual corresponding to Equation (7),
% the constraint becomes hat{Sigma}_S-(A Sigma)_S+eta I_S >= 0.
for k = 1:numel(patterns)
    S = patterns{k};
    sk = numel(S);
    b(offset+(1:sk^2)) = SigmaS{k}(:);

    for col = 1:sk
        for row = 1:sk
            pos = offset+row+(col-1)*sk;
            globalRow = S(row);
            globalCol = S(col);
            if globalRow==globalCol
                A(pos,1) = -1;
            else
                A(pos,offIndex(globalRow,globalCol)) = -1;
            end
            if regularized && row==col
                % x(zIndex)=alpha*eta, hence x(zIndex)/alpha=eta is
                % added to each diagonal of the pattern slack matrix.
                A(pos,zIndex) = 1/regularizationAlpha;
            end
        end
    end
    offset = offset+sk^2;
end

% The optimization variable in the final scalar block is t=alpha*eta,
% where alpha=2/hat{c}=1/(hat{c}/2). This is the reciprocal of the
% regularization level in Algorithm 1, Step 11 and Equation (7). The
% rescaling improves numerical balance; the scalar cone imposes eta>=0.
if regularized
    A(offset+1,zIndex) = 1;
end

% PAPER: Equation (4) maximizes lambda. The regularized dual maximizes
% lambda-alpha*eta. ADMM is written as a minimization, so cobj encodes
% -lambda (unregularized) or -lambda+alpha*eta (regularized).
cobj = zeros(q,1);
cobj(1) = -1;
if regularized
    cobj(zIndex) = 1;
end

end

% -------------------------------------------------------------------------
function projected = localProjectProductCone(value,blocks)
%localProjectProductCone Project onto the PSD cones in Equations (4)/(7).
%
% This eigenvalue-thresholding operation is the Euclidean projection onto a
% PSD cone. It is an ADMM implementation detail, not a separate statistical
% construction in the paper.

projected = zeros(size(value));
offset = 0;
for k = 1:numel(blocks)
    sk = blocks(k);
    idx = offset+(1:sk^2);
    if sk==1
        projected(idx) = max(value(idx),0);
    else
        M = reshape(value(idx),sk,sk);
        M = (M+M')/2;
        [V,D] = eig(M);
        eigenvalues = real(diag(D));
        M = V*diag(max(eigenvalues,0))*V';
        projected(idx) = reshape(real((M+M')/2),[],1);
    end
    offset = offset+sk^2;
end

end

% -------------------------------------------------------------------------
function [Sigma,z] = localUnpackSolution( ...
    x,d,offIndex,zIndex,regularizationAlpha)
%localUnpackSolution Recover the common matrix in Proposition 5, Eq. (4).
%
% The common diagonal is x(1)=lambda. Dividing this matrix by lambda in the
% main routine produces the compatible correlation matrix hat{Q} appearing
% in Equation (5). In the regularized problem z is the auxiliary eta slack.

Sigma = x(1)*eye(d);
for j = 2:d
    for i = 1:j-1
        Sigma(i,j) = x(offIndex(i,j));
        Sigma(j,i) = Sigma(i,j);
    end
end
Sigma = (Sigma+Sigma')/2;
if isempty(zIndex)
    z = 0;
else
    z = x(zIndex)/regularizationAlpha;
end

end

% -------------------------------------------------------------------------
function [R,Z] = localSampleCorrelation(X)
%localSampleCorrelation Compute SampleCorr and standardized observations.
%
% PAPER: Section 2 (Annals p. 2209) denotes these quantities by
% hat{Sigma}_S, hat{mu}_S and hat{sigma}_S^2. The routine is used in
% Algorithm 1,
% Step 2 for the observed data and Step 10 for each bootstrap sample.

nr = size(X,1);
% hat{mu}_S and square roots of the diagonal sample variances.
mu = mean(X,1);
scale = std(X,0,1);
if any(~isfinite(scale)) || any(scale<=0)
    error('FSDA:mdBBtest:ZeroPatternVariance', ...
        ['At least one variable has zero or nonfinite variance within a ' ...
         'retained missingness pattern.']);
end
% Row-vector form of diag(hat{sigma}_S^2)^(-1/2)(X_S-hat{mu}_S).
Z = (X-mu)./scale;
% Pearson sample correlation hat{Sigma}_S=SampleCorr(X_S).
R = (Z'*Z)/(nr-1);
R = real((R+R')/2);
R(1:size(R,1)+1:end) = 1;

end

% -------------------------------------------------------------------------
function [R,clipMagnitude] = localCov2Corr(S)
%localCov2Corr Repair and normalize the compatible matrix in Equation (5).
%
% In exact arithmetic SigmaOpt/lambda is a correlation matrix. The ADMM
% optimizer is recovered from the affine variable x rather than from its
% projection onto the cone, so it satisfies Sigma>=0 only up to the primal
% residual and can carry a small negative eigenvalue. Every use of the
% fitted compatible matrix downstream requires a genuine correlation matrix:
% Algorithm 1, Step 7 takes hat{Q}_S^(1/2), and a matrix that is not
% positive semidefinite has no real square root.
%
% The spectrum is therefore clipped at zero and the result is rescaled to a
% unit diagonal. The order matters: rescaling preserves positive
% semidefiniteness, whereas overwriting the diagonal after clipping would
% destroy it. The operation is a no-op to machine precision when the solve
% is accurate, and clipMagnitude reports how much repair was needed so that
% an inaccurate solve is visible rather than silent.

S = real((S+S')/2);
[V,D] = eig(S);
eigenvalues = real(diag(D));
clipMagnitude = max(0,-min(eigenvalues));
S = V*diag(max(eigenvalues,0))*V';
S = real((S+S')/2);

diagonal = real(diag(S));
if any(diagonal<=0) || any(~isfinite(diagonal))
    error('FSDA:mdBBtest:InvalidCompatibleDiagonal', ...
        'The fitted compatible matrix has a nonpositive diagonal element.');
end
Dinv = diag(1./sqrt(diagonal));
R = Dinv*S*Dinv;
R = real((R+R')/2);
R(1:size(R,1)+1:end) = 1;

end

% -------------------------------------------------------------------------
function root = localMatrixSquareRoot(S,requirePositiveDefinite,tol)
%localMatrixSquareRoot Compute the matrix roots in Algorithm 1, Step 7.
%
% The transformation uses hat{Sigma}_S^(-1/2), which requires positive
% definiteness, and hat{Q}_S^(1/2), for which positive semidefiniteness is
% sufficient. Eigenvalue truncation removes only numerical negative drift.

S = real((S+S')/2);
[V,D] = eig(S);
eigenvalues = real(diag(D));
scale = max(1,max(abs(eigenvalues)));
if min(eigenvalues)<-100*tol*scale
    error('FSDA:mdBBtest:NonPSDMatrix', ...
        'A matrix required by the bootstrap transformation is not PSD.');
end
if requirePositiveDefinite && min(eigenvalues)<=tol^2*scale
    error('FSDA:mdBBtest:NonPDMatrix', ...
        ['A pattern-specific correlation matrix required by the bootstrap ' ...
         'transformation is not positive definite.']);
end
root = V*diag(sqrt(max(eigenvalues,0)))*V';
root = real((root+root')/2);

end

% -------------------------------------------------------------------------
function localPlotResults(out,originalPatternIndex)
%localPlotResults Plot numerical results from Algorithm 1.
%
% These figures are FSDA diagnostics and are not part of Algorithm 1. The
% histogram displays the Step-11 bootstrap statistics against the Step-3
% observed statistic; the second plot visualizes the retained collection S.

if ~isempty(out.Rboot)
    figure('Name','mdBBtest bootstrap distribution','Color','w')
    histogram(out.Rboot)
    hold on
    xline(out.stat,'r','LineWidth',1.5)
    xlabel('Regularized incompatibility index')
    ylabel('Frequency')
    title(sprintf('Observed R = %.4g; bootstrap p-value = %.4g', ...
        out.stat,out.pvalue))
    box on
end

figure('Name','mdBBtest retained patterns','Color','w')
imagesc(double(out.patterns))
clim([0 1])
colormap(gca,[0.85 0.85 0.85; 0.20 0.45 0.75])
colorbar('Ticks',[0.25 0.75], ...
    'TickLabels',{'Missing','Observed'})
xticks(1:out.p)
xticklabels(out.variableNames)
xtickangle(45)
yticks(1:out.npatterns)
yticklabels("P"+string(originalPatternIndex)+ ...
    " (n="+string(out.patternCounts)+")")
xlabel('Variable')
ylabel('Retained missingness pattern')
title('Observed-variable patterns used by mdBBtest')

end
