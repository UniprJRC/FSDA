function out = mdMAARtest(Y, varargin)
%mdMAARtest Diagnostic tests for Missing Always At Random (MAAR).
%
%<a href="matlab: docsearchFS('mdMAARtest')">Link to the help function</a>
%
%  mdMAARtest implements the three diagnostic procedures proposed by
%  Bojinov, Pillai and Rubin (2020) for investigating whether partially
%  observed variables are associated with missingness indicators after
%  conditioning on the fully observed variables.
%
%  The available procedures are:
%
%  1) 'ccm': comparison of conditional means;
%  2) 'dtmm': direct test of a postulated missingness mechanism;
%  3) 'cop': semiparametric Gaussian-copula diagnostic.
%
%  The null hypothesis is that the restrictions implied by MAAR, row
%  exchangeability and conditional independence of the missingness
%  indicators are compatible with the observed data. Rejection indicates
%  that at least one of these assumptions is violated. The procedures are
%  diagnostic tools and do not establish that a nonrejected mechanism is
%  MAAR.
%
%  Required input arguments:
%
%    Y : Input data. Matrix, table or timetable. n x p data matrix possibly
%        containing NaN values. Rows of Y represent observations and columns
%        represent variables.
%        Data Types - single | double | table | timetable
%
%  Optional input arguments:
%
%       method : Diagnostic method. Character or string.
%                Possible values are:
%                'ccm'  = comparison of conditional means;
%                'dtmm' = direct test of the missingness mechanism;
%                'cop'  = Gaussian-copula diagnostic.
%                The default value is 'ccm'.
%                Example - 'method','cop'
%                Data Types - char | string
%
%        alpha : Significance level. Scalar in the interval (0,1).
%                Bonferroni corrections are applied internally to the
%                component tests. The default value is 0.05.
%                Example - 'alpha',0.01
%                Data Types - single | double
%
%      imputed : Completed data sets used by method 'dtmm'. Cell array or
%                n x p x M numeric array. Each completed data set must have
%                the same dimensions as Y and must contain no missing values.
%                If empty, mdEM and mdImputeStochastic are used to generate
%                the imputations. The default value is [].
%                Example - 'imputed',{Yimp1,Yimp2,Yimp3}
%                Data Types - cell | single | double
%
% nimputations : Number of stochastic imputations for method 'dtmm' when
%                option imputed is empty. Positive integer greater than 1.
%                The default value is 10.
%                Example - 'nimputations',20
%                Data Types - single | double
%
%      maxiter : Maximum number of iterations passed to mdEM when imputations
%                must be generated. The default value is 100.
%                Example - 'maxiter',200
%                Data Types - single | double
%
%          tol : Convergence tolerance passed to mdEM. The default value is
%                1e-5.
%                Example - 'tol',1e-8
%                Data Types - single | double
%
%        nsamp : Number of MCMC iterations for method 'cop'. The default
%                value is 2000.
%                Example - 'nsamp',5000
%                Data Types - single | double
%
%        nburn : Number of initial MCMC iterations discarded for method
%                'cop'. The default value is 500.
%                Example - 'nburn',1000
%                Data Types - single | double
%
%         thin : MCMC thinning interval for method 'cop'. The default value
%                is 1.
%                Example - 'thin',5
%                Data Types - single | double
%
% pluginthreshold : Threshold used by the Gaussian-copula sampler. Margins
%                having more than pluginthreshold distinct observed values
%                use fixed normal scores, as in sbgcop. The default value is
%                100.
%                Example - 'pluginthreshold',50
%                Data Types - single | double
%
%
%        ridge : Nonnegative numerical regularization used only when a
%                covariance or information matrix is nearly singular. The
%                default value is 1e-8.
%                Example - 'ridge',1e-6
%                Data Types - single | double
%
%          msg : Display a summary of the results. Boolean. The default value
%                is true.
%                Example - 'msg',false
%                Data Types - logical | single | double
%
%        plots : Produce a diagnostic plot. Boolean. The default value is
%                false.
%                Example - 'plots',true
%                Data Types - logical | single | double
%
%  savesamples : Save the posterior conditional-covariance draws produced by
%                method 'cop'. Boolean. The default value is false.
%                Example - 'savesamples',true
%                Data Types - logical | single | double
%
%  Output:
%
%    out : Structure containing the following fields:
%
%          out.method        = selected diagnostic method.
%          out.reject        = true if at least one component test rejects.
%          out.whichReject   = logical vector associated with the variables
%                              containing missing values. Element k is true
%                              when the missingness indicator of that variable
%                              is implicated by the diagnostic.
%          out.pvalue        = component p-values. For 'ccm' and 'cop' this
%                              is a q x q matrix, where q is the number of
%                              partially observed variables. For 'dtmm' it is
%                              a q x 1 vector.
%          out.stat          = method-specific test statistics.
%          out.alpha         = nominal significance level.
%          out.alphaAdjusted = Bonferroni-adjusted component level.
%          out.missingVariables = indices of variables containing NaNs.
%          out.fullyObservedVariables = indices of fully observed variables.
%          out.variableNames = variable names.
%          out.missingIndicator = n x q matrix; 1 denotes a missing value.
%          out.n             = number of observations.
%          out.p             = number of variables.
%          out.nvarmiss      = number q of partially observed variables.
%          out.interpretation = concise interpretation of the result.
%
%          Additional fields for method 'ccm':
%          out.df1, out.df2  = degrees of freedom of the nested-model F tests.
%          out.rejectPairs   = q x q logical matrix of rejected comparisons.
%
%          Additional fields for method 'dtmm':
%          out.df1, out.df2  = D3 reference-distribution degrees of freedom.
%          out.relativeIncrease = relative increase in variance due to
%                              nonresponse in the pooled likelihood-ratio test.
%          out.deviance      = detailed deviances used in the D3 calculation.
%          out.imputationInfo = information about the completed data sets.
%
%          Additional fields for method 'cop':
%          out.CI            = q x q x 2 array of simultaneous credible limits.
%          out.posteriorMean = posterior means of conditional covariances.
%          out.posteriorMedian = posterior medians.
%          out.rejectPairs   = q x q logical matrix of rejected associations.
%          out.nSaved        = number of retained MCMC draws.
%          out.samples       = posterior draws when savesamples is true.
%
%  More About:
%
% This function implements the three diagnostics from Bojinov, Pillai and
% Rubin:
% 
% ccm: comparison of conditional means using nested Gaussian linear models
%       and Bonferroni correction, following diagMAAR.ccm. 
% dtmm: direct testing of the missingness mechanism through logistic
%       models and multiple imputation. The likelihood-ratio statistics are
%       combined using the Meng– Rubin $D_3$ procedure. ​
% cop: Gaussian-copula diagnostic based on posterior conditional
%      covariances between partially observed variables and missingness
%       indicators.
%
%  Let $Jm$ denote the set of variables that contain missing values and $Jf$ the
%  set of fully observed variables. Under the assumptions considered by
%  Bojinov, Pillai and Rubin, a missingness indicator $R_k$ should not depend on
%  another partially observed variable $Y_j$ after conditioning on $Y_Jf$.
%
%  The comparison-of-conditional-means procedure tests this implication by
%  comparing, for every $j$ different from $k$, the nested linear models
%
%    \[
%     Y_j ~ Y_{Jf}
%    \]
%
%  and
%
%    \[
%     Y_j ~ Y_Jf + R_k + Y_Jf:R_k.
%    \]
%
%  The direct procedure compares logistic models for each $R_k$. The reduced
%  model contains only the fully observed variables, whereas the full model
%  also contains the imputed partially observed variables. The likelihood
%  ratio statistics are combined across imputations using the D3 procedure
%  of Meng and Rubin (1992), as used by the original diagMAAR implementation.
%
%  The Gaussian-copula procedure uses the extended rank-likelihood sampler of
%  Hoff (2007). For every retained posterior correlation matrix, the function
%  computes the covariance between the partially observed outcomes and their
%  missingness indicators conditional on the fully observed variables. A
%  component is rejected when its Bonferroni-adjusted credible interval does
%  not contain zero.
%
%  The original R implementation uses mice for method 'dtmm'. When option
%  imputed is empty, this MATLAB implementation instead uses the joint-normal
%  EM and stochastic-imputation functions already available in FSDA. To use
%  exactly the same completed data sets as another implementation, supply
%  them explicitly through option imputed.
%
%  See also: mdEM, mdImputeStochastic, mdMCARtest, mdLittleTest, mdJJtest
%
%  References:
%
%  Bojinov, I., Pillai, N. S. and Rubin, D. B. (2020), "Diagnosing missing
%  always at random in multivariate data", Biometrika, Vol. 107, pp. 246-253.
%
%  Hoff, P. D. (2007), "Extending the rank likelihood for semiparametric
%  copula estimation", Annals of Applied Statistics, Vol. 1, pp. 265-283.
%
%  Meng, X.-L. and Rubin, D. B. (1992), "Performing likelihood ratio tests
%  with multiply-imputed data sets", Biometrika, Vol. 79, pp. 103-111.
%
%  Copyright 2008-2026.
%  Written by FSDA team
%
%<a href="matlab: docsearchFS('mdMAARtest')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
   %% Example 1: MAAR data and comparison of conditional means.
   rng(1)
   n = 500;
   p = 5;
   Yfull = randn(n,p);
   eta = -1 + Yfull(:,4) - Yfull(:,5);
   prob = 1./(1+exp(-eta));
   Y = Yfull;
   for j = 1:3
       Y(rand(n,1)<prob,j) = NaN;
   end
   out = mdMAARtest(Y,'method','ccm','plots',true);
   disp(out.pvalue)
%}
%
%{
   %% Example 2: Violation of MAAR detected by conditional means.
   rng(2)
   n = 800;
   Yfull = randn(n,5);
   Y = Yfull;
   p1 = 1./(1+exp(-(-1+Yfull(:,4)-Yfull(:,5))));
   R1 = rand(n,1)<p1;
   p2 = 1./(1+exp(-(-1+1.5*Yfull(:,1)+Yfull(:,4)-Yfull(:,5))));
   R2 = rand(n,1)<p2;
   p3 = 1./(1+exp(-(-1+Yfull(:,1)+Yfull(:,2)+Yfull(:,4)-Yfull(:,5))));
   R3 = rand(n,1)<p3;
   Y(R1,1)=NaN; Y(R2,2)=NaN; Y(R3,3)=NaN;
   out = mdMAARtest(Y,'method','ccm','plots',true);
   disp(out.whichReject)
%}
%
%{
   %% Example 3: Direct test with multiple stochastic imputations.
   rng(3)
   n = 500;
   Yfull = randn(n,5);
   Y = Yfull;
   pr = 1./(1+exp(-(-1+Yfull(:,4)-Yfull(:,5))));
   for j=1:3
       Y(rand(n,1)<pr,j)=NaN;
   end
   out = mdMAARtest(Y,'method','dtmm','nimputations',10);
   disp(out.pvalue)
%}

%{
   %% Example 4: Gaussian-copula diagnostic.
   % A small number of iterations is used only to keep the example fast.
   rng(4)
   n = 300;
   Yfull = randn(n,5);
   Y = Yfull;
   pr = 1./(1+exp(-(-1+Yfull(:,4)-Yfull(:,5))));
   for j=1:3
       Y(rand(n,1)<pr,j)=NaN;
   end
   out = mdMAARtest(Y,'method','cop','nsamp',600,'nburn',200, ...
       'plots',true);
   disp(out.posteriorMedian)
%}

%% Beginning of code

if nargin < 1
    error('FSDA:mdMAARtest:TooFewInputs', ...
        'At least one input argument is required.');
end

% Preserve variable names before converting tables to arrays.
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
    error('FSDA:mdMAARtest:WrongInput', ...
        'Input argument Y must be a numeric matrix, table or timetable.');
end

Y = double(Y);
[n,p] = size(Y);
if n < 5 || p < 2
    error('FSDA:mdMAARtest:SmallData', ...
        'At least five observations and two variables are required.');
end
if any(all(isnan(Y),1))
    error('FSDA:mdMAARtest:AllMissingVariable', ...
        'At least one variable contains only missing values.');
end

% Default options.
method = 'ccm';
alpha = 0.05;
imputed = [];
nimputations = 10;
maxiter = 100;
tol = 1e-5;
nsamp = 2000;
nburn = 500;
thin = 1;
pluginthreshold = 100;
seed = [];
ridge = 1e-8;
msg = true;
plots = false;
savesamples = false;

if ~isempty(varargin)
    options = struct('method',method,'alpha',alpha,'imputed',imputed, ...
        'nimputations',nimputations,'maxiter',maxiter,'tol',tol, ...
        'nsamp',nsamp,'nburn',nburn,'thin',thin, ...
        'pluginthreshold',pluginthreshold,'seed',seed,'ridge',ridge, ...
        'msg',msg,'plots',plots,'savesamples',savesamples);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mdMAARtest:WrongInputOpt', ...
            'Number of supplied options is invalid. Values may be missing.');
    end
    if ~isempty(UserOptions)
        aux.chkoptions(options,UserOptions)
    end
    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end

    method = options.method;
    alpha = options.alpha;
    imputed = options.imputed;
    nimputations = options.nimputations;
    maxiter = options.maxiter;
    tol = options.tol;
    nsamp = options.nsamp;
    nburn = options.nburn;
    thin = options.thin;
    pluginthreshold = options.pluginthreshold;
    seed = options.seed;
    ridge = options.ridge;
    msg = options.msg;
    plots = options.plots;
    savesamples = options.savesamples;
end

method = lower(char(method));
if ~ismember(method,{'ccm','dtmm','cop'})
    error('FSDA:mdMAARtest:WrongMethod', ...
        'Option method must be ''ccm'', ''dtmm'' or ''cop''.');
end
if ~isscalar(alpha) || ~isnumeric(alpha) || alpha <= 0 || alpha >= 1
    error('FSDA:mdMAARtest:WrongAlpha', ...
        'Option alpha must be a scalar in the interval (0,1).');
end
if ~isscalar(nimputations) || nimputations <= 1 || ...
        nimputations ~= floor(nimputations)
    error('FSDA:mdMAARtest:WrongNImputations', ...
        'Option nimputations must be an integer greater than 1.');
end
if ~isscalar(maxiter) || maxiter <= 0 || maxiter ~= floor(maxiter)
    error('FSDA:mdMAARtest:WrongMaxiter', ...
        'Option maxiter must be a positive integer.');
end
if ~isscalar(tol) || tol <= 0
    error('FSDA:mdMAARtest:WrongTol', ...
        'Option tol must be a positive scalar.');
end
if ~isscalar(nsamp) || nsamp <= 1 || nsamp ~= floor(nsamp)
    error('FSDA:mdMAARtest:WrongNsamp', ...
        'Option nsamp must be an integer greater than 1.');
end
if ~isscalar(nburn) || nburn < 0 || nburn ~= floor(nburn) || nburn >= nsamp
    error('FSDA:mdMAARtest:WrongNburn', ...
        'Option nburn must be a nonnegative integer smaller than nsamp.');
end
if ~isscalar(thin) || thin < 1 || thin ~= floor(thin)
    error('FSDA:mdMAARtest:WrongThin', ...
        'Option thin must be a positive integer.');
end
if ~isscalar(pluginthreshold) || pluginthreshold < 2 || ...
        pluginthreshold ~= floor(pluginthreshold)
    error('FSDA:mdMAARtest:WrongPluginThreshold', ...
        'Option pluginthreshold must be an integer greater than 1.');
end
if ~isempty(seed) && (~isscalar(seed) || ~isnumeric(seed) || ~isfinite(seed))
    error('FSDA:mdMAARtest:WrongSeed', ...
        'Option seed must be empty or a finite numeric scalar.');
end
if ~isscalar(ridge) || ~isnumeric(ridge) || ridge < 0
    error('FSDA:mdMAARtest:WrongRidge', ...
        'Option ridge must be a nonnegative scalar.');
end
if ~(isscalar(msg) && (islogical(msg) || isnumeric(msg)))
    error('FSDA:mdMAARtest:WrongMsg', ...
        'Option msg must be a logical or numeric scalar.');
end
if ~(isscalar(plots) && (islogical(plots) || isnumeric(plots)))
    error('FSDA:mdMAARtest:WrongPlots', ...
        'Option plots must be a logical or numeric scalar.');
end
if ~(isscalar(savesamples) && (islogical(savesamples) || isnumeric(savesamples)))
    error('FSDA:mdMAARtest:WrongSaveSamples', ...
        'Option savesamples must be a logical or numeric scalar.');
end
msg = logical(msg);
plots = logical(plots);
savesamples = logical(savesamples);


missingMask = isnan(Y);
locmis = find(any(missingMask,1));
locfull = find(~any(missingMask,1));
q = numel(locmis);
if q == 0
    error('FSDA:mdMAARtest:NoMissingValues', ...
        'Input data do not contain missing values.');
end

% Constant variables make all three diagnostics ill-defined.
for j = 1:p
    yj = Y(~isnan(Y(:,j)),j);
    if numel(unique(yj)) < 2
        error('FSDA:mdMAARtest:ConstantVariable', ...
            'Variable %d is constant on its observed values.',j);
    end
end

R = double(missingMask(:,locmis));

% Common output fields.
out = struct;
out.method = method;
out.alpha = alpha;
out.n = n;
out.p = p;
out.nvarmiss = q;
out.variableNames = variableNames;
out.missingVariables = locmis(:);
out.fullyObservedVariables = locfull(:);
out.missingIndicator = R;
out.reject = false;
out.whichReject = false(q,1);
out.pvalue = [];
out.stat = [];
out.alphaAdjusted = alpha;
out.interpretation = '';

switch method
    case 'ccm'
        res = ccmTest(Y,R,locmis,locfull,alpha,ridge);
        out = copyFields(out,res);

    case 'dtmm'
        [imputations,impInfo] = prepareImputations(Y,imputed,nimputations, ...
            maxiter,tol);
        res = dtmmTest(imputations,R,locmis,locfull,alpha,ridge,maxiter,tol);
        out = copyFields(out,res);
        out.imputationInfo = impInfo;

    case 'cop'
        res = copulaTest(Y,R,locmis,locfull,alpha,nsamp,nburn,thin, ...
            pluginthreshold,ridge,savesamples,msg);
        out = copyFields(out,res);
end

out.interpretation = interpretMAAR(out);

if msg
    printMAAR(out)
end
if plots
    plotMAAR(out)
end

end

% -------------------------------------------------------------------------
function out = copyFields(out,res)
fields = fieldnames(res);
for i = 1:numel(fields)
    out.(fields{i}) = res.(fields{i});
end
end

% -------------------------------------------------------------------------
function res = ccmTest(Y,R,locmis,locfull,alpha,ridge)
% Comparison of conditional means, translated from diagMAAR.ccm.
q = numel(locmis);
pvalue = ones(q,q);
stat = zeros(q,q);
df1 = zeros(q,q);
df2 = zeros(q,q);

Xfull = Y(:,locfull);
for jj = 1:q
    response = Y(:,locmis(jj));
    for kk = 1:q
        if kk == jj
            continue
        end
        rgroup = R(:,kk);
        use = ~isnan(response);
        if ~isempty(Xfull)
            use = use & all(~isnan(Xfull),2);
        end
        y = response(use);
        xf = Xfull(use,:);
        rg = rgroup(use);
        Xsmall = [ones(numel(y),1) xf];
        Xlarge = [Xsmall rg xf.*rg];
        test = nestedLinearF(y,Xsmall,Xlarge,ridge);
        pvalue(kk,jj) = test.pvalue;
        stat(kk,jj) = test.stat;
        df1(kk,jj) = test.df1;
        df2(kk,jj) = test.df2;
    end
end

nTests = max(q*(q-1),1);
alphaAdjusted = alpha/nTests;
rejectPairs = pvalue < alphaAdjusted;
rejectPairs(1:q+1:end) = false;
whichReject = any(rejectPairs,2);

if q == 1
    warning('FSDA:mdMAARtest:CCMOneMissingVariable', ...
        ['The CCM procedure requires cross-comparisons between partially ' ...
        'observed variables. With one such variable there is no component test.']);
end

res = struct;
res.pvalue = pvalue;
res.stat = stat;
res.df1 = df1;
res.df2 = df2;
res.alphaAdjusted = alphaAdjusted;
res.rejectPairs = rejectPairs;
res.whichReject = whichReject;
res.reject = any(whichReject);
end

% -------------------------------------------------------------------------
function test = nestedLinearF(y,Xsmall,Xlarge,ridge)
% F test for two nested Gaussian linear models.
ns = size(Xsmall,1);
rs = rank(Xsmall);
rl = rank(Xlarge);
dfdiff = rl-rs;
dfres = ns-rl;

if dfdiff <= 0 || dfres <= 0
    test = struct('stat',NaN,'df1',dfdiff,'df2',dfres,'pvalue',NaN);
    return
end

bs = stableLeastSquares(Xsmall,y,ridge);
bl = stableLeastSquares(Xlarge,y,ridge);
rsdSmall = y-Xsmall*bs;
rsdLarge = y-Xlarge*bl;
sseSmall = sum(rsdSmall.^2);
sseLarge = sum(rsdLarge.^2);
num = max(sseSmall-sseLarge,0)/dfdiff;
den = sseLarge/dfres;
if den <= 0
    F = Inf;
    pval = 0;
else
    F = num/den;
    pval = fSurvival(F,dfdiff,dfres);
end

test = struct('stat',F,'df1',dfdiff,'df2',dfres,'pvalue',pval);
end

% -------------------------------------------------------------------------
function b = stableLeastSquares(X,y,ridge)
% Solve least squares and regularize only in nearly singular cases.
[~,R] = qr(X,0);
if isempty(R) || size(R,1) ~= size(R,2) || rcond(R) < 1e-12
    scale = max(1,trace(X'*X)/size(X,2));
    ridgeUse = max(ridge,eps(scale));
    b = (X'*X+ridgeUse*scale*eye(size(X,2)))\(X'*y);
else
    b = X\y;
end
end

% -------------------------------------------------------------------------
function [imputations,info] = prepareImputations(Y,imputed,M,maxiter,tol)
% Normalize supplied imputations or create FSDA stochastic imputations.
[n,p] = size(Y);
if isempty(imputed)
    if exist('mdEM','file') ~= 2 || exist('mdImputeStochastic','file') ~= 2
        error('FSDA:mdMAARtest:MissingFSDAFunctions', ...
            ['Functions mdEM and mdImputeStochastic must be on the MATLAB ' ...
            'path when option imputed is empty.']);
    end
    emOut = mdEM(Y,'maxiter',maxiter,'tol',tol);
    imputations = cell(M,1);
    for m = 1:M
        imputations{m} = mdImputeStochastic(Y,emOut.loc,emOut.cov);
    end
    info = struct('source','mdEM/mdImputeStochastic', ...
        'nimputations',M,'loc',emOut.loc,'cov',emOut.cov);
else
    if iscell(imputed)
        imputations = imputed(:);
    elseif isnumeric(imputed) && ndims(imputed) <= 3
        if ismatrix(imputed)
            imputations = {double(imputed)};
        else
            imputations = cell(size(imputed,3),1);
            for m = 1:size(imputed,3)
                imputations{m} = double(imputed(:,:,m));
            end
        end
    else
        error('FSDA:mdMAARtest:WrongImputed', ...
            'Option imputed must be a cell array or an n x p x M numeric array.');
    end

    for m = 1:numel(imputations)
        if istable(imputations{m})
            imputations{m} = table2array(imputations{m});
        end
        if ~isnumeric(imputations{m}) || ...
                ~isequal(size(imputations{m}),[n p])
            error('FSDA:mdMAARtest:WrongImputedSize', ...
                'Every completed data set must have the same size as Y.');
        end
        imputations{m} = double(imputations{m});
        if any(isnan(imputations{m}(:)))
            error('FSDA:mdMAARtest:IncompleteImputation', ...
                'A supplied completed data set still contains NaN values.');
        end
        obs = ~isnan(Y);
        tolerance = 100*eps(max(1,max(abs(Y(obs)))));
        if any(abs(imputations{m}(obs)-Y(obs)) > tolerance)
            warning('FSDA:mdMAARtest:ObservedValuesChanged', ...
                'Supplied imputation %d changes at least one observed value.',m);
        end
    end
    if numel(imputations) < 2
        error('FSDA:mdMAARtest:TooFewImputations', ...
            'At least two completed data sets are required for method dtmm.');
    end
    info = struct('source','user supplied','nimputations',numel(imputations), ...
        'loc',[],'cov',[]);
end
end

% -------------------------------------------------------------------------
function res = dtmmTest(imputations,R,locmis,locfull,alpha,ridge,maxiter,tol)
% Direct missingness-mechanism test with Meng-Rubin D3 pooling.
M = numel(imputations);
q = numel(locmis);
pvalue = NaN(q,1);
stat = NaN(q,1);
df1 = q*ones(q,1);
df2 = NaN(q,1);
relativeIncrease = NaN(q,1);
deviance = repmat(struct('fullMLE',[],'smallMLE',[], ...
    'fullPooled',[],'smallPooled',[]),q,1);

for kk = 1:q
    response = R(:,kk);
    bFull = [];
    bSmall = [];
    devFullMLE = zeros(M,1);
    devSmallMLE = zeros(M,1);

    for m = 1:M
        Ym = imputations{m};
        Xsmall = [ones(size(Ym,1),1) Ym(:,locfull)];
        Xfull = [Xsmall Ym(:,locmis)];
        fitSmall = logisticFit(Xsmall,response,maxiter,tol,ridge);
        fitFull = logisticFit(Xfull,response,maxiter,tol,ridge);
        if m == 1
            bSmall = zeros(numel(fitSmall.beta),M);
            bFull = zeros(numel(fitFull.beta),M);
        end
        bSmall(:,m) = fitSmall.beta;
        bFull(:,m) = fitFull.beta;
        devSmallMLE(m) = fitSmall.deviance;
        devFullMLE(m) = fitFull.deviance;
    end

    qbarSmall = mean(bSmall,2);
    qbarFull = mean(bFull,2);
    devSmallPooled = zeros(M,1);
    devFullPooled = zeros(M,1);
    for m = 1:M
        Ym = imputations{m};
        Xsmall = [ones(size(Ym,1),1) Ym(:,locfull)];
        Xfull = [Xsmall Ym(:,locmis)];
        devSmallPooled(m) = logisticDeviance(Xsmall,response,qbarSmall);
        devFullPooled(m) = logisticDeviance(Xfull,response,qbarFull);
    end

    devM = mean(devSmallMLE-devFullMLE);
    devL = mean(devSmallPooled-devFullPooled);
    k = size(bFull,1)-size(bSmall,1);
    r = ((M+1)/(k*(M-1)))*(devM-devL);
    % A slightly negative value can occur because of numerical error.
    r = max(r,0);
    Dm = max(devL,0)/(k*(1+r));
    v = k*(M-1);
    if r <= sqrt(eps)
        w = Inf;
    elseif v > 4
        w = 4+(v-4)*(1+(1-2/v)/r)^2;
    else
        w = v*(1+1/k)*(1+1/r)^2/2;
    end

    pvalue(kk) = fSurvival(Dm,k,w);
    stat(kk) = Dm;
    df1(kk) = k;
    df2(kk) = w;
    relativeIncrease(kk) = r;
    deviance(kk).fullMLE = devFullMLE;
    deviance(kk).smallMLE = devSmallMLE;
    deviance(kk).fullPooled = devFullPooled;
    deviance(kk).smallPooled = devSmallPooled;
end

alphaAdjusted = alpha/q;
whichReject = pvalue < alphaAdjusted;
res = struct;
res.pvalue = pvalue;
res.stat = stat;
res.df1 = df1;
res.df2 = df2;
res.relativeIncrease = relativeIncrease;
res.deviance = deviance;
res.alphaAdjusted = alphaAdjusted;
res.whichReject = whichReject;
res.reject = any(whichReject);
end

% -------------------------------------------------------------------------
function fit = logisticFit(X,y,maxiter,tol,ridge)
% Logistic maximum likelihood by safeguarded iteratively reweighted squares.
y = double(y(:));
X = double(X);
k = size(X,2);
beta = zeros(k,1);
ll = logisticLogLik(X,y,beta);
converged = false;

for iter = 1:maxiter
    eta = X*beta;
    mu = logisticCDF(eta);
    w = max(mu.*(1-mu),1e-10);
    gradient = X'*(y-mu);
    information = X'*(w.*X);
    scale = max(1,trace(information)/max(k,1));
    if rcond(information) < 1e-12
        ridgeUse = max(ridge,eps(scale));
        information = information+ridgeUse*scale*eye(k);
    end
    step = information\gradient;

    stepFactor = 1;
    betaNew = beta+step;
    llNew = logisticLogLik(X,y,betaNew);
    while llNew < ll && stepFactor > 2^-20
        stepFactor = stepFactor/2;
        betaNew = beta+stepFactor*step;
        llNew = logisticLogLik(X,y,betaNew);
    end

    if max(abs(betaNew-beta)) <= tol*(1+max(abs(beta)))
        beta = betaNew;
        ll = llNew;
        converged = true;
        break
    end
    beta = betaNew;
    ll = llNew;
end

if ~converged
    warning('FSDA:mdMAARtest:LogisticNoConvergence', ...
        'A logistic model did not converge in %d iterations.',maxiter);
end
fit = struct('beta',beta,'loglik',ll,'deviance',-2*ll, ...
    'iter',iter,'converged',converged);
end

% -------------------------------------------------------------------------
function dev = logisticDeviance(X,y,beta)
dev = -2*logisticLogLik(X,double(y(:)),beta);
end

% -------------------------------------------------------------------------
function ll = logisticLogLik(X,y,beta)
eta = X*beta;
ll = sum(y.*(-softplus(-eta))+(1-y).*(-softplus(eta)));
end

% -------------------------------------------------------------------------
function y = softplus(x)
y = max(x,0)+log1p(exp(-abs(x)));
end

% -------------------------------------------------------------------------
function p = logisticCDF(x)
p = zeros(size(x));
pos = x >= 0;
p(pos) = 1./(1+exp(-x(pos)));
ex = exp(x(~pos));
p(~pos) = ex./(1+ex);
end

% -------------------------------------------------------------------------
function res = copulaTest(Y,R,locmis,locfull,alpha,nsamp,nburn,thin, ...
    pluginthreshold,ridge,savesamples,msg)
% Gaussian-copula diagnostic based on the extended rank likelihood.
q = numel(locmis);
W = [Y(:,locmis) R Y(:,locfull)];
[n,d] = size(W);
nSaved = floor((nsamp-nburn)/thin);
if nSaved < 20
    warning('FSDA:mdMAARtest:FewCopulaDraws', ...
        'Only %d posterior draws will be retained.',nSaved);
end

% Starting latent values and ordinal levels.
Z = zeros(n,d);
levelIndex = zeros(n,d);
nLevels = zeros(1,d);
plugin = false(1,d);
for j = 1:d
    observed = ~isnan(W(:,j));
    x = W(observed,j);
    [~,~,lev] = unique(x,'sorted');
    levelIndex(observed,j) = lev;
    nLevels(j) = max(lev);
    plugin(j) = numel(unique(x)) > pluginthreshold;
    rmax = rankMaximum(x);
    u = rmax/(numel(x)+1);
    Z(observed,j) = normalInv(u);
    Z(~observed,j) = randn(sum(~observed),1);
end

S = cov(Z,1);
S = stabilizeSPD(S,ridge);
interestSamples = zeros(q,q,nSaved);
saveIndex = 0;
S0 = eye(d);
n0 = d+2;

for iteration = 1:nsamp
    order = randperm(d);
    for jj = 1:d
        j = order(jj);
        other = [1:j-1 j+1:d];
        if isempty(other)
            conditionalMean = zeros(n,1);
            conditionalSD = sqrt(max(S(j,j),eps));
        else
            Soo = stabilizeSPD(S(other,other),ridge);
            regression = S(j,other)/Soo;
            conditionalMean = Z(:,other)*regression';
            conditionalVariance = S(j,j)-regression*S(other,j);
            conditionalSD = sqrt(max(conditionalVariance,eps));
        end

        observed = ~isnan(W(:,j));
        if ~plugin(j)
            for lev = 1:nLevels(j)
                idx = observed & levelIndex(:,j)==lev;
                if ~any(idx)
                    continue
                end
                if lev == 1
                    lower = -Inf;
                else
                    lower = max(Z(observed & levelIndex(:,j)==lev-1,j));
                end
                if lev == nLevels(j)
                    upper = Inf;
                else
                    upper = min(Z(observed & levelIndex(:,j)==lev+1,j));
                end
                Z(idx,j) = truncatedNormal(conditionalMean(idx), ...
                    conditionalSD,lower,upper);
            end
        end
        idxMissing = ~observed;
        if any(idxMissing)
            Z(idxMissing,j) = conditionalMean(idxMissing)+ ...
                conditionalSD*randn(sum(idxMissing),1);
        end
    end

    Psi = n0*S0+Z'*Z;
    invPsi = stableInverse(Psi,ridge);
    precisionDraw = wishartRandom(invPsi,n0+n,ridge);
    S = stableInverse(precisionDraw,ridge);

    if iteration > nburn && mod(iteration-nburn,thin)==0
        saveIndex = saveIndex+1;
        sdS = sqrt(max(diag(S),eps));
        C = S./(sdS*sdS');
        C = (C+C')/2;
        nv = 2*q;
        if nv < d
            C22 = stabilizeSPD(C(nv+1:end,nv+1:end),ridge);
            Vcond = C(1:nv,1:nv)- ...
                C(1:nv,nv+1:end)/C22*C(nv+1:end,1:nv);
        else
            Vcond = C(1:nv,1:nv);
        end
        interestSamples(:,:,saveIndex) = Vcond(1:q,q+1:2*q);
    end

    if msg && nsamp >= 10 && mod(iteration,max(1,floor(nsamp/10)))==0
        fprintf('mdMAARtest copula sampler: %d%% completed.\n', ...
            round(100*iteration/nsamp));
    end
end

interestSamples = interestSamples(:,:,1:saveIndex);
nSaved = saveIndex;
nTests = max(q^2-q,1);
alphaAdjusted = alpha/nTests;
lowerProb = alphaAdjusted/2;
upperProb = 1-lowerProb;

lower = zeros(q,q);
medianValue = zeros(q,q);
upper = zeros(q,q);
posteriorMean = zeros(q,q);
pvalue = ones(q,q);
for j = 1:q
    for k = 1:q
        draws = squeeze(interestSamples(j,k,:));
        qq = quantile(draws,[lowerProb 0.5 upperProb]);
        lower(j,k) = qq(1);
        medianValue(j,k) = qq(2);
        upper(j,k) = qq(3);
        posteriorMean(j,k) = mean(draws);
        pvalue(j,k) = min(1,2*min(mean(draws<=0),mean(draws>=0)));
    end
end

rejectPairs = (lower>0 | upper<0);
if q > 1
    rejectPairs(1:q+1:end) = false;
    pvalue(1:q+1:end) = 1;
end
whichReject = any(rejectPairs,1)';

res = struct;
res.pvalue = pvalue;
res.stat = medianValue;
res.CI = cat(3,lower,upper);
res.posteriorMean = posteriorMean;
res.posteriorMedian = medianValue;
res.alphaAdjusted = alphaAdjusted;
res.rejectPairs = rejectPairs;
res.whichReject = whichReject;
res.reject = any(whichReject);
res.nSaved = nSaved;
res.pluginMarginal = plugin;
if savesamples
    res.samples = interestSamples;
else
    res.samples = [];
end
end

% -------------------------------------------------------------------------
function ranks = rankMaximum(x)
% Ranks with ties assigned their maximum rank.
x = x(:);
[sorted,order] = sort(x);
ranksSorted = zeros(size(x));
first = 1;
while first <= numel(x)
    last = first;
    while last < numel(x) && sorted(last+1)==sorted(first)
        last = last+1;
    end
    ranksSorted(first:last) = last;
    first = last+1;
end
ranks = zeros(size(x));
ranks(order) = ranksSorted;
end

% -------------------------------------------------------------------------
function x = truncatedNormal(mu,sigma,lower,upper)
% Draw N(mu,sigma^2) variates restricted to [lower,upper].
a = (lower-mu)/sigma;
b = (upper-mu)/sigma;
pLower = normalCDF(a);
pUpper = normalCDF(b);
width = pUpper-pLower;
u = pLower+rand(size(mu)).*width;
u = min(max(u,realmin),1-eps);
x = mu+sigma*normalInv(u);

bad = ~isfinite(x) | width <= 10*eps(max(pUpper,pLower));
if any(bad)
    badIndex = find(bad);
    for ii = 1:numel(badIndex)
        i = badIndex(ii);
        accepted = false;
        for attempt = 1:10000
            proposal = mu(i)+sigma*randn;
            if proposal >= lower && proposal <= upper
                x(i) = proposal;
                accepted = true;
                break
            end
        end
        if ~accepted
            lo = lower;
            hi = upper;
            if isinf(lo), lo = mu(i)-8*sigma; end
            if isinf(hi), hi = mu(i)+8*sigma; end
            x(i) = min(max(mu(i),lo+sqrt(eps)),hi-sqrt(eps));
        end
    end
end
end

% -------------------------------------------------------------------------
function p = normalCDF(x)
p = 0.5*erfc(-x/sqrt(2));
end

% -------------------------------------------------------------------------
function x = normalInv(p)
p = min(max(p,realmin),1-eps);
x = -sqrt(2)*erfcinv(2*p);
end

% -------------------------------------------------------------------------
function W = wishartRandom(V,nu,ridge)
% Bartlett decomposition for W ~ Wishart(V,nu).
p = size(V,1);
if nu < p
    error('FSDA:mdMAARtest:WishartDF', ...
        'Wishart degrees of freedom must be at least the matrix dimension.');
end
V = stabilizeSPD(V,ridge);
L = chol(V,'lower');
A = zeros(p,p);
for i = 1:p
    A(i,i) = sqrt(2*randg((nu-i+1)/2));
    if i > 1
        A(i,1:i-1) = randn(1,i-1);
    end
end
W = L*(A*A')*L';
W = (W+W')/2;
end

% -------------------------------------------------------------------------
function A = stabilizeSPD(A,ridge)
A = (A+A')/2;
if isempty(A)
    return
end
base = max(1,trace(abs(A))/size(A,1));
[~,flag] = chol(A);
if flag ~= 0 || rcond(A) < 1e-12
    [V,D] = eig(A,'vector');
    floorValue = max(ridge*base,eps(base));
    D = max(real(D),floorValue);
    A = V*diag(D)*V';
    A = (A+A')/2;
end
end

% -------------------------------------------------------------------------
function Ainv = stableInverse(A,ridge)
A = stabilizeSPD(A,ridge);
Ainv = A\eye(size(A));
Ainv = (Ainv+Ainv')/2;
end

% -------------------------------------------------------------------------
function p = fSurvival(x,df1,df2)
% Upper-tail F probability, including the df2=Inf limiting case.
if isnan(x) || x < 0
    p = NaN;
elseif isinf(x)
    p = 0;
elseif isinf(df2)
    p = gammainc(df1*x/2,df1/2,'upper');
else
    z = df2/(df2+df1*x);
    p = betainc(z,df2/2,df1/2);
end
end

% -------------------------------------------------------------------------
function interpretation = interpretMAAR(out)
if out.reject
    implicated = out.missingVariables(out.whichReject);
    if isempty(implicated)
        interpretation = ['At least one component test rejects the restrictions ' ...
            'implied by MAAR and the accompanying assumptions.'];
    else
        interpretation = sprintf(['The diagnostic rejects the restrictions ' ...
            'implied by MAAR. Missingness indicators associated with variables %s ' ...
            'are implicated.'],mat2str(implicated(:)'));
    end
else
    interpretation = ['The diagnostic does not reject the restrictions implied ' ...
        'by MAAR. This result does not prove that the mechanism is MAAR.'];
end
end

% -------------------------------------------------------------------------
function printMAAR(out)
fprintf('\nMAAR diagnostic: %s\n',upper(out.method));
fprintf('Observations: %d; variables: %d; partially observed variables: %d\n', ...
    out.n,out.p,out.nvarmiss);
fprintf('Nominal alpha: %.4g; component alpha: %.4g\n', ...
    out.alpha,out.alphaAdjusted);
fprintf('Overall rejection: %s\n',char(string(out.reject)));
if any(out.whichReject)
    fprintf('Implicated variable indices: %s\n', ...
        mat2str(out.missingVariables(out.whichReject)'));
else
    fprintf('No missingness indicator was implicated.\n');
end
fprintf('%s\n\n',out.interpretation);
end

% -------------------------------------------------------------------------
function plotMAAR(out)
q = out.nvarmiss;
labels = out.variableNames(out.missingVariables);
figure('Name',['mdMAARtest: ' upper(out.method)],'Color','w');

switch out.method
    case 'ccm'
        values = -log10(max(out.pvalue,realmin));
        values(1:q+1:end) = 0;
        imagesc(values)
        colorbar
        axis tight
        xlabel('Partially observed response variable')
        ylabel('Missingness indicator')
        title('-log_{10} p-values: comparison of conditional means')
        xticks(1:q); yticks(1:q)
        xticklabels(labels); yticklabels(labels)

    case 'dtmm'
        bar(-log10(max(out.pvalue,realmin)))
        hold on
        yline(-log10(out.alphaAdjusted),'--','Adjusted threshold');
        hold off
        xticks(1:q); xticklabels(labels)
        ylabel('-log_{10} p-value')
        xlabel('Missingness indicator')
        title('Direct missingness-mechanism diagnostic')

    case 'cop'
        imagesc(out.posteriorMedian)
        colorbar
        axis tight
        xlabel('Missingness indicator')
        ylabel('Partially observed outcome')
        title('Posterior median conditional covariance')
        xticks(1:q); yticks(1:q)
        xticklabels(labels); yticklabels(labels)
        hold on
        [rr,cc] = find(out.rejectPairs);
        plot(cc,rr,'kx','MarkerSize',12,'LineWidth',2)
        hold off
end
end

%FScategory:MULT-MissingData
