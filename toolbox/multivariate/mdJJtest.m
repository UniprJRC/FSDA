function out = mdJJtest(Y, varargin)
%mdJJtest Jamshidian-Jalal test for MCAR in incomplete continuous data.
%
%<a href="matlab: docsearchFS('mdJJtest')">Link to the help function</a>
%
%  mdJJtest implements the procedure of Jamshidian and Jalal (2010) for
%  testing the Missing Completely At Random (MCAR) assumption. The testing
%  backend follows the implementation of function mcar in the R package
%  mice. The test is intended for continuous variables.
%
%  Required input arguments:
%
%    Y :           Input data. Matrix or table. n x p data matrix possibly
%                  containing missing values (NaN's). Rows of Y represent
%                  observations and columns represent variables.
%                  Data Types - single | double | table
%
%  Optional input arguments:
%
%    imputed :     Completed data sets. Empty value, cell array or 3D array.
%                  If supplied, imputed can be:
%                  (i) an n x p x M numeric array; or
%                  (ii) a cell array containing M n x p numeric matrices or
%                  tables. Observed entries should not be changed. If empty,
%                  mdEM and mdImputeStochastic are used internally.
%                  The default value is [].
%                  Example - 'imputed',Yimp
%                  Data Types - double | single | cell
%
% nimputations :   Number of stochastic imputations. Positive integer.
%                  This option is used only when imputed is empty.
%                  The default value is 5.
%                  Example - 'nimputations',10
%                  Data Types - double | single
%
%       minn :     Minimum pattern size. Positive integer greater than 1.
%                  Missingness patterns whose frequency is less than 
%                  minn are removed. The default value is 6.
%                  Example - 'minn',10
%                  Data Types - double | single
%
%     method :     Test method. Character vector or string.
%                  Possible values are:
%                  'auto'          : Hawkins' test is computed first. If at
%                                    least one Hawkins p-value is smaller
%                                    than alpha, the nonparametric
%                                    Anderson-Darling test is also computed;
%                  'hawkins'       : only Hawkins' test is reported;
%                  'nonparametric' : the Anderson-Darling rank test is
%                                    always reported.
%                  The default value is 'auto'.
%                  Example - 'method','nonparametric'
%                  Data Types - char | string
%
%     nsimul :     Number of Monte Carlo replications used to calibrate the
%                  fourth-order Neyman smooth test in small patterns.
%                  The default value is 10000.
%                  Example - 'nsimul',50000
%                  Data Types - double | single
%
%  usechisq :      Threshold for the asymptotic Neyman reference law.
%                  For a pattern containing at least usechisq observations,
%                  the chi-square approximation is used instead of Monte
%                  Carlo calibration. The default value is 30.
%                  Example - 'usechisq',50
%                  Data Types - double | single
%
%      alpha :     Significance level. Scalar in the interval (0,1).
%                  The default value is 0.05.
%                  Example - 'alpha',0.01
%                  Data Types - double | single
%
%
%    maxiter :     Maximum number of iterations for mdEM. Positive integer.
%                  This option is used only when imputed is empty.
%                  The default value is 200.
%                  Example - 'maxiter',500
%                  Data Types - double | single
%
%        tol :     Convergence tolerance for mdEM. Positive scalar.
%                  This option is used only when imputed is empty.
%                  The default value is 1e-7.
%                  Example - 'tol',1e-8
%                  Data Types - double | single
%
%      ridge :     Numerical regularization. Nonnegative scalar.
%                  A small diagonal ridge is used to stabilize the pooled
%                  covariance matrix in Hawkins' test. The default value is
%                  1e-10.
%                  Example - 'ridge',1e-8
%                  Data Types - double | single
%
%        msg :     Display test results. Boolean.
%                  If true, a concise summary is displayed.
%                  The default value is true.
%                  Example - 'msg',false
%                  Data Types - logical | double
%
%      plots :     Plotting flag. Boolean.
%                  If true, the function displays the p-value histograms
%                  and the retained missingness-pattern matrix.
%                  The default value is false.
%                  Example - 'plots',true
%                  Data Types - logical | double
%
%  Output:
%
%    out :         Structure containing the following fields:
%
%          out.stat          = main test statistic. With multiple
%                              imputations, this is the median statistic;
%          out.pvalue        = main p-value. With multiple imputations, this
%                              is the median p-value;
%          out.df            = degrees of freedom of the main test. This is
%                              empty for the Anderson-Darling test;
%          out.method        = method used;
%          out.alpha         = significance level;
%          out.hawkinsStat   = M x 1 vector containing the combined Hawkins
%                              statistics, one for each imputation;
%          out.hawkinsDF     = degrees of freedom of the combined Hawkins
%                              statistic;
%          out.hawkinsP      = M x 1 vector of Hawkins p-values;
%          out.ADstat        = M x 1 vector of Anderson-Darling statistics,
%                              when computed;
%          out.ADpvalue      = M x 1 vector of Anderson-Darling p-values,
%                              when computed;
%          out.patterns      = g x p logical matrix containing the retained
%                              missingness patterns. A true entry denotes a
%                              missing variable;
%          out.patternCounts = g x 1 vector containing the frequencies of
%                              the retained missingness patterns;
%          out.patternInfo   = table containing pattern frequency and number
%                              of missing variables. Row names identify the
%                              missingness patterns;
%          out.npatterns     = number of retained missingness patterns;
%          out.idxPatterns   = nused x 1 pattern membership vector;
%          out.removedRows   = n x 1 logical vector identifying observations
%                              removed because of sparse patterns;
%          out.removedPatterns = matrix containing removed patterns;
%          out.loc           = EM estimate of location used for the default
%                              stochastic imputations. Empty if imputed is
%                              supplied;
%          out.cov           = EM estimate of covariance used for the default
%                              stochastic imputations. Empty if imputed is
%                              supplied;
%          out.details       = cell array with detailed results for each
%                              completed data set;
%          out.interpretation = concise interpretation of the test.
%
%  More About:
%
%  The procedure first completes the data and then divides each completed
%  data set into groups according to the missingness patterns in the
%  original data. Hawkins' transformation produces, within each pattern,
%  values which should be uniformly distributed when multivariate
%  normality and covariance homogeneity hold. Uniformity is assessed by a
%  fourth-order Neyman smooth test. Pattern-specific p-values are combined
%  by Fisher's method.
%
%  If method='auto' and Hawkins' test is significant, the k-sample
%  Anderson-Darling rank test is applied to distinguish a departure caused
%  by nonnormality from evidence against MCAR. The procedure only addresses
%  the MCAR-versus-MAR framework. It cannot distinguish MCAR from MNAR, and
%  a nonsignificant result does not rule out MNAR.
%
%  When imputed is empty, this function estimates location and covariance
%  by mdEM and generates stochastic completions using
%  mdImputeStochastic. This is a joint multivariate-normal imputation and is
%  not identical to the default chained normal-regression imputation used
%  by mice::mcar. To compare MATLAB and R using exactly the same completed
%  data sets, supply them through option imputed.
%
%  See also: mdEM, mdImputeStochastic, mdLittleTest, mdMCARtest
%
%  References:
%
%  Jamshidian, M. and Jalal, S. (2010), "Tests of Homoscedasticity,
%  Normality, and Missing Completely at Random for Incomplete Multivariate
%  Data", Psychometrika, Vol. 75, pp. 649-674.
%
%  Van Buuren, S. and Groothuis-Oudshoorn, K. (2011), "mice: Multivariate
%  Imputation by Chained Equations in R", Journal of Statistical Software,
%  Vol. 45, pp. 1-67.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdJJtest')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
   %% Example 1: Jamshidian-Jalal test under MCAR.
   rng(10)
   n = 500;
   p = 5;
   Y = randn(n,p);
   Y(rand(n,p)<0.20) = NaN;
   out = mdJJtest(Y);
   disp(out.pvalue)
%}
%
%{
   %% Example 2: Display diagnostic plots.
   rng(20)
   n = 600;
   p = 6;
   Y = randn(n,p);
   Y(rand(n,p)<0.15) = NaN;
   out = mdJJtest(Y,'plots',true);
%}
%
%{
   %% Example 3: Use the nonparametric branch directly.
   rng(30)
   n = 500;
   p = 5;
   Y = trnd(5,n,p);
   Y(rand(n,p)<0.20) = NaN;
   out = mdJJtest(Y,'method','nonparametric');
   disp(out.ADpvalue)
%}
%
%{
   %% Example 4: Supply completed data sets.
   rng(40)
   n = 300;
   p = 4;
   Y = randn(n,p);
   Y(rand(n,p)<0.10) = NaN;
   outEM = mdEM(Y);
   Yimp = cell(3,1);
   for j = 1:3
       Yimp{j} = mdImputeStochastic(Y,outEM.loc,outEM.cov);
   end
   out = mdJJtest(Y,'imputed',Yimp,'method','hawkins');
%}

%% Beginning of code

% Input parameters checking
if nargin < 1
    error('FSDA:mdJJtest:TooFewInputs', ...
        'At least one input argument is required.');
end

% Store variable names before transforming a table into a numeric matrix.
if istable(Y)
    Y = Y{:,:};
end

if ~ismatrix(Y) || ~isnumeric(Y)
    error('FSDA:mdJJtest:WrongInput', ...
        'Input argument Y must be a numeric matrix or a numeric table.');
end

Y = double(Y);
[nOriginal,p] = size(Y);

if nOriginal < 4 || p < 2
    error('FSDA:mdJJtest:TooSmall', ...
        'At least four observations and two variables are required.');
end

if any(all(isnan(Y),1))
    error('FSDA:mdJJtest:AllMissingColumn', ...
        'At least one variable contains only missing values.');
end

if any(all(isnan(Y),2))
    warning('FSDA:mdJJtest:AllMissingRow', ...
        ['Some rows contain only missing values. The test may be unreliable ' ...
         'for these observations.']);
end

if any(columnRangeIgnoringNaN(Y)==0)
    error('FSDA:mdJJtest:ConstantColumn', ...
        'At least one variable is constant on its observed values.');
end

% Default options
imputed = [];
nimputations = 5;
minn = 6;
method = 'auto';
nsimul = 10000;
usechisq = 30;
alpha = 0.05;
maxiter = 200;
tol = 1e-7;
ridge = 1e-10;
msg = true;
plots = false;

% Optional arguments
if ~isempty(varargin)
    options = struct('imputed',imputed,'nimputations',nimputations, ...
        'minn',minn,'method',method,'nsimul',nsimul, ...
        'usechisq',usechisq,'alpha',alpha, ...
        'maxiter',maxiter,'tol',tol,'ridge',ridge, ...
        'msg',msg,'plots',plots);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));

    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:mdJJtest:WrongInputOpt', ...
                ['Number of supplied options is invalid. Probably values ' ...
                 'for some parameters are missing.']);
        end
        aux.chkoptions(options,UserOptions)
    end

    for j = 1:2:length(varargin)
        options.(varargin{j}) = varargin{j+1};
    end

    imputed = options.imputed;
    nimputations = options.nimputations;
    minn = options.minn;
    method = options.method;
    nsimul = options.nsimul;
    usechisq = options.usechisq;
    alpha = options.alpha;
    maxiter = options.maxiter;
    tol = options.tol;
    ridge = options.ridge;
    msg = options.msg;
    plots = options.plots;
end

% Basic checks on options
if ~isempty(imputed) && ~(iscell(imputed) || isnumeric(imputed))
    error('FSDA:mdJJtest:WrongImputed', ...
        'Option ''imputed'' must be empty, a cell array or a numeric array.');
end

if ~isscalar(nimputations) || nimputations < 1 || ...
        nimputations ~= floor(nimputations)
    error('FSDA:mdJJtest:WrongNimputations', ...
        'Option ''nimputations'' must be a positive integer.');
end

if ~isscalar(minn) || minn <= 1 || minn ~= floor(minn)
    error('FSDA:mdJJtest:WrongMinn', ...
        'Option ''minn'' must be an integer greater than 1.');
end

if ~(ischar(method) || (isstring(method) && isscalar(method)))
    error('FSDA:mdJJtest:WrongMethod', ...
        'Option ''method'' must be a character vector or a string scalar.');
end
method = lower(char(method));
if ~ismember(method,{'auto','hawkins','nonparametric'})
    error('FSDA:mdJJtest:WrongMethod', ...
        ['Option ''method'' must be ''auto'', ''hawkins'' or ' ...
         '''nonparametric''.']);
end

if ~isscalar(nsimul) || nsimul < 100 || nsimul ~= floor(nsimul)
    error('FSDA:mdJJtest:WrongNsimul', ...
        'Option ''nsimul'' must be an integer greater than or equal to 100.');
end

if ~isscalar(usechisq) || usechisq < 2 || usechisq ~= floor(usechisq)
    error('FSDA:mdJJtest:WrongUsechisq', ...
        'Option ''usechisq'' must be an integer greater than or equal to 2.');
end

if ~isscalar(alpha) || ~isnumeric(alpha) || alpha <= 0 || alpha >= 1
    error('FSDA:mdJJtest:WrongAlpha', ...
        'Option ''alpha'' must be a scalar in the interval (0,1).');
end

if ~isscalar(maxiter) || maxiter < 1 || maxiter ~= floor(maxiter)
    error('FSDA:mdJJtest:WrongMaxiter', ...
        'Option ''maxiter'' must be a positive integer.');
end

if ~isscalar(tol) || ~isnumeric(tol) || tol <= 0
    error('FSDA:mdJJtest:WrongTol', ...
        'Option ''tol'' must be a positive scalar.');
end

if ~isscalar(ridge) || ~isnumeric(ridge) || ridge < 0
    error('FSDA:mdJJtest:WrongRidge', ...
        'Option ''ridge'' must be a nonnegative scalar.');
end

if ~(islogical(msg) || (isnumeric(msg) && isscalar(msg)))
    error('FSDA:mdJJtest:WrongMsg', ...
        'Option ''msg'' must be a logical or numeric scalar.');
end
if ~(islogical(plots) || (isnumeric(plots) && isscalar(plots)))
    error('FSDA:mdJJtest:WrongPlots', ...
        'Option ''plots'' must be a logical or numeric scalar.');
end
msg = logical(msg);
plots = logical(plots);


% Identify the original missingness patterns.
missingMaskOriginal = isnan(Y);
[patternsOriginal,idxPatternsOriginal,countsOriginal] = ...
    patternGroups(missingMaskOriginal);

% Construct or validate completed data sets.
loc = [];
covmat = [];
if isempty(imputed)
    outEM = mdEM(Y,'maxiter',maxiter,'tol',tol, ...
        'Patterns',patternsOriginal,'idxPatterns',idxPatternsOriginal);
    loc = outEM.loc;
    covmat = outEM.cov;
    imputations = cell(nimputations,1);
    for j = 1:nimputations
        imputations{j} = mdImputeStochastic(Y,loc,covmat);
    end
else
    imputations = normalizeImputedInput(imputed,Y);
end

% Remove sparse patterns. The current mice source removes groups with
% frequency less than or equal to minn.
removePattern = countsOriginal <= minn;
removedRows = removePattern(idxPatternsOriginal);
removedPatterns = patternsOriginal(removePattern,:);
removedPatternCounts = countsOriginal(removePattern);

if all(removedRows)
    error('FSDA:mdJJtest:NoRows', ...
        ['All observations belong to missingness patterns with frequency ' ...
         'less than or equal to minn. Lower minn with caution.']);
end

Yused = Y(~removedRows,:);
imputations = cellfun(@(Z) Z(~removedRows,:),imputations, ...
    'UniformOutput',false);

[patterns,idxPatterns,patternCounts] = patternGroups(isnan(Yused));
g = size(patterns,1);
nused = size(Yused,1);

if g < 2
    error('FSDA:mdJJtest:OnePattern', ...
        'At least two retained missingness patterns are required.');
end

if any(patternCounts < 2)
    error('FSDA:mdJJtest:TinyPattern', ...
        'At least two observations are required in every retained pattern.');
end

if nused-g-p <= 0
    error('FSDA:mdJJtest:InsufficientDF', ...
        'Hawkins'' transformation requires nused-g-p to be positive.');
end

% Compute Hawkins' quantities for every completed data set.
M = numel(imputations);
details = cell(M,1);
for j = 1:M
    details{j} = hawkinsCore(imputations{j},idxPatterns,ridge);
end

hawkinsStat = [];
hawkinsP = [];
hawkinsDF = 2*g;

if strcmp(method,'auto') || strcmp(method,'hawkins')
    hawkinsStat = zeros(M,1);
    hawkinsP = zeros(M,1);

    for j = 1:M
        groupP = zeros(g,1);
        for r = 1:g
            groupP(r) = neymanPValue(details{j}.A{r},nsimul,usechisq);
            if groupP(r)==0
                groupP(r) = 1/nsimul;
            end
        end

        hawkinsStat(j) = -2*sum(log(groupP));
        hawkinsP(j) = chi2cdf(hawkinsStat(j),hawkinsDF,"upper");
        details{j}.NeymanP = groupP;
    end
end

% Compute the k-sample Anderson-Darling test when required.
ADstat = [];
ADpvalue = [];
runAD = strcmp(method,'nonparametric') || ...
    (strcmp(method,'auto') && any(hawkinsP<alpha));

if runAD
    ADstat = zeros(M,1);
    ADpvalue = zeros(M,1);
    for j = 1:M
        ADout = andersonDarlingCore(details{j}.Fij);
        ADstat(j) = ADout.Statistic;
        ADpvalue(j) = ADout.PValue;
        details{j}.AndersonDarling = ADout;
    end
end

% Select the main statistic and p-value.
if strcmp(method,'hawkins') || (strcmp(method,'auto') && isempty(ADpvalue))
    stat = median(hawkinsStat);
    pvalue = median(hawkinsP);
    df = hawkinsDF;
else
    stat = median(ADstat);
    pvalue = median(ADpvalue);
    df = [];
end

% Create the pattern-information table.
patternNames = cell(g,1);
for r = 1:g
    patternNames{r} = sprintf('%d',patterns(r,:));
end
patternInfo = table(patternCounts,sum(patterns,2), ...
    'VariableNames',{'Frequency','NMissing'}, ...
    'RowNames',patternNames);

% Store results.
out = struct;
out.stat = stat;
out.pvalue = pvalue;
out.df = df;
out.method = method;
out.alpha = alpha;
out.minn = minn;
out.nsimul = nsimul;
out.usechisq = usechisq;
out.nimputations = M;
out.n = nOriginal;
out.nused = nused;
out.p = p;
out.patterns = patterns;
out.patternCounts = patternCounts;
out.patternInfo = patternInfo;
out.npatterns = g;
out.idxPatterns = idxPatterns;
out.removedRows = removedRows;
out.removedPatterns = removedPatterns;
out.removedPatternCounts = removedPatternCounts;
out.hawkinsStat = hawkinsStat;
out.hawkinsDF = hawkinsDF;
out.hawkinsP = hawkinsP;
out.ADstat = ADstat;
out.ADpvalue = ADpvalue;
out.loc = loc;
out.cov = covmat;
out.details = details;
out.interpretation = interpretResult(method,hawkinsP,ADpvalue,alpha);

if msg
    printSummary(out)
end

if plots
    plotResults(out)
end

end

% -------------------------------------------------------------------------
function imputations = normalizeImputedInput(imputed,Y)
%normalizeImputedInput Convert supplied completed data sets to a cell array.

[n,p] = size(Y);

if iscell(imputed)
    imputations = imputed(:);
elseif isnumeric(imputed) && ndims(imputed)<=3
    if ismatrix(imputed)
        imputations = {double(imputed)};
    else
        M = size(imputed,3);
        imputations = cell(M,1);
        for j = 1:M
            imputations{j} = double(imputed(:,:,j));
        end
    end
else
    error('FSDA:mdJJtest:WrongImputed', ...
        ['Option ''imputed'' must be a cell array or an n x p x M ' ...
         'numeric array.']);
end

observed = ~isnan(Y);
observedValues = Y(observed);
if isempty(observedValues)
    toleranceObserved = 0;
else
    toleranceObserved = 100*eps(max(1,max(abs(observedValues))));
end

for j = 1:numel(imputations)
    if istable(imputations{j})
        imputations{j} = imputations{j}{:,:};
    end

    if ~isnumeric(imputations{j}) || ...
            ~isequal(size(imputations{j}),[n p])
        error('FSDA:mdJJtest:WrongImputedSize', ...
            'Every completed data set must have the same size as Y.');
    end

    imputations{j} = double(imputations{j});

    if any(isnan(imputations{j}(:)))
        error('FSDA:mdJJtest:IncompleteImputation', ...
            'A supplied completed data set still contains NaN values.');
    end

    if any(abs(imputations{j}(observed)-Y(observed))>toleranceObserved)
        warning('FSDA:mdJJtest:ObservedChanged', ...
            'Supplied imputation %d changes at least one observed value.',j);
    end
end
end

% -------------------------------------------------------------------------
function [patterns,idxPatterns,counts] = patternGroups(missingMask)
%patternGroups Group rows having the same missingness pattern.

[patterns,~,idxPatterns] = unique(missingMask,'rows','stable');
counts = accumarray(idxPatterns,1,[size(patterns,1) 1]);
end

% -------------------------------------------------------------------------
function out = hawkinsCore(Y,idxPatterns,ridge)
%hawkinsCore Compute Hawkins' transformed quantities.

[n,p] = size(Y);
groups = unique(idxPatterns,'stable');
g = length(groups);

pooledSS = zeros(p,p);
centeredGroups = cell(g,1);
groupSizes = zeros(g,1);

for r = 1:g
    Yr = Y(idxPatterns==groups(r),:);
    groupSizes(r) = size(Yr,1);
    Yr = Yr-mean(Yr,1);
    centeredGroups{r} = Yr;
    pooledSS = pooledSS+Yr'*Yr;
end

pooledCov = pooledSS/(n-g);
pooledCov = stabilizeSPD(pooledCov,ridge,'pooled covariance');
invPooledCov = pooledCov\eye(p);
invPooledCov = (invPooledCov+invPooledCov')/2;

Fij = cell(g,1);
A = cell(g,1);
df2 = n-g-p;

for r = 1:g
    Yr = centeredGroups{r};
    mahal = sum((Yr*invPooledCov).*Yr,2);
    h = groupSizes(r)*mahal;
    den = p*((groupSizes(r)-1)*(n-g)-h);

    if any(den<=0)
        error('FSDA:mdJJtest:HawkinsDenominator', ...
            ['A nonpositive denominator occurred in Hawkins'' ' ...
             'transformation. Check pattern sizes, collinearity or the ' ...
             'supplied imputations.']);
    end

    Fij{r} = ((n-g-p)*h)./den;
    A{r} = fSurvival(Fij{r},p,df2);
end

out = struct;
out.Fij = Fij;
out.A = A;
out.groupSizes = groupSizes;
out.pooledCov = pooledCov;
end

% -------------------------------------------------------------------------
function pvalue = neymanPValue(u,nsimul,usechisq)
%neymanPValue Fourth-order Neyman smooth test for uniformity.

u = u(:);
n = length(u);
Phi = shiftedLegendre(u,4);
stat = sum(sum(Phi,1).^2)/n;

if n<usechisq
    U = rand(n,nsimul);
    z = 2*U-1;
    P1 = z;
    P2 = (3*z.^2-1)/2;
    P3 = (5*z.*P2-2*P1)/3;
    P4 = (7*z.*P3-3*P2)/4;
    sumPol = [sqrt(3)*sum(P1,1); sqrt(5)*sum(P2,1); ...
        sqrt(7)*sum(P3,1); 3*sum(P4,1)];
    simulated = sum(sumPol.^2,1)/n;
    pvalue = mean(simulated>stat);
else
    pvalue = chi2cdf(stat,4,"upper");
end
end

% -------------------------------------------------------------------------
function Phi = shiftedLegendre(u,order)
%shiftedLegendre Orthonormal shifted Legendre polynomials.

u = u(:);
z = 2*u-1;
Phi = zeros(length(u),order);
Pprev = ones(size(z));
Pcurr = z;
Phi(:,1) = sqrt(3)*Pcurr;

for degree = 2:order
    Pnext = ((2*degree-1)*z.*Pcurr-(degree-1)*Pprev)/degree;
    Phi(:,degree) = sqrt(2*degree+1)*Pnext;
    Pprev = Pcurr;
    Pcurr = Pnext;
end
end

% -------------------------------------------------------------------------
function out = andersonDarlingCore(Fij)
%andersonDarlingCore k-sample Anderson-Darling test 

allF = vertcat(Fij{:});
groupSizes = cellfun(@numel,Fij);
k = length(groupSizes);
n = length(allF);

if k<2
    error('FSDA:mdJJtest:ADGroups', ...
        'At least two groups are required for the Anderson-Darling test.');
end

if n<4
    error('FSDA:mdJJtest:ADSample', ...
        'At least four observations are required for the Anderson-Darling test.');
end

xsort = sort(allF);
xsort(end) = [];
isNew = [true; diff(xsort)~=0];
zj = xsort(isNew);
first = find(isNew);
last = [first(2:end)-1; length(xsort)];
hj = last-first+1;
hn = cumsum(hj);

ADgroupRaw = zeros(k,1);
for r = 1:k
    fr = Fij{r}(:);
    counts = zeros(length(zj),1);
    for j = 1:length(zj)
        counts(j) = sum(fr==zj(j));
    end
    mij = cumsum(counts);
    num = (n*mij-groupSizes(r)*hn).^2;
    den = hn.*(n-hn);
    ADgroupRaw(r) = sum(hj.*(num./den))/groupSizes(r);
end

ADstat = sum(ADgroupRaw)/n;
ADgroup = ADgroupRaw/n;
J = sum(1./groupSizes);
harmonic = cumsum(1./(1:(n-1)));
H = harmonic(end);
idx = 1:(n-2);
G = sum((H-harmonic(idx))./(n-idx));

a = (4*G-6)*(k-1)+(10-6*G)*J;
b = (2*G-4)*k^2+8*H*k+(2*G-14*H-4)*J-8*H+4*G-6;
c = (6*H+2*G-2)*k^2+(4*H-4*G+6)*k+(2*H-6)*J+4*H;
d = (2*H+6)*k^2-4*H*k;
variance = max(((a*n^3)+(b*n^2)+(c*n)+d)/ ...
    ((n-1)*(n-2)*(n-3)),0);

if variance==0
    standardized = sign(ADstat-(k-1))*Inf;
else
    standardized = (ADstat-(k-1))/sqrt(variance);
end

b0 = [0.675 1.281 1.645 1.960 2.326];
b1 = [-0.245 0.250 0.678 1.149 1.822];
b2 = [-0.105 -0.305 -0.362 -0.391 -0.396];
c0 = [1.09861228866811 2.19722457733622 2.94443897916644 ...
    3.66356164612965 4.59511985013459];
qnt = b0+b1/sqrt(k-1)+b2/(k-1);

if standardized<=qnt(3)
    take = 2:5;
else
    take = 1:4;
end

yy = spline(qnt(take),c0(take),standardized);
if yy>=0
    ey = exp(-yy);
    pvalue = ey/(1+ey);
else
    pvalue = 1/(1+exp(yy));
end

out = struct;
out.PValue = pvalue;
out.Statistic = ADstat;
out.GroupContributions = ADgroup;
out.Standardized = standardized;
out.Variance = variance;
end

% -------------------------------------------------------------------------
function A = stabilizeSPD(A,ridge,label)
%stabilizeSPD Symmetrize and regularize a covariance matrix if necessary.

A = (A+A')/2;
if isempty(A)
    return
end

base = max(1,trace(abs(A))/size(A,1));
if ridge>0
    A = A+ridge*base*eye(size(A));
end

[~,flag] = chol(A);
if flag~=0
    [V,D] = eig(A);
    d = real(diag(D));
    floorValue = max(ridge*base,eps(base));
    d = max(d,floorValue);
    A = V*diag(d)*V';
    A = (A+A')/2;
    warning('FSDA:mdJJtest:Regularized', ...
        'The %s was projected onto the positive-definite cone.',label);
end
end

% -------------------------------------------------------------------------
function pvalue = fSurvival(x,df1,df2)
%fSurvival Upper-tail probability of the F distribution.

x = max(x,0);
z = df2./(df2+df1*x);
pvalue = betainc(z,df2/2,df1/2);
pvalue(isinf(x)) = 0;
end


% -------------------------------------------------------------------------
function txt = interpretResult(method,hawkinsP,ADpvalue,alpha)
%interpretResult Produce a concise interpretation of the selected test.

if strcmp(method,'hawkins')
    if median(hawkinsP)<alpha
        txt = ['Hawkins'' test is significant. Under multivariate ' ...
            'normality this is evidence against MCAR.'];
    else
        txt = ['Hawkins'' test is not significant: there is no evidence ' ...
            'against multivariate normality or MCAR.'];
    end
elseif strcmp(method,'nonparametric')
    if median(ADpvalue)<alpha
        txt = ['The Anderson-Darling rank test is significant: reject ' ...
            'MCAR within the MCAR-versus-MAR framework.'];
    else
        txt = ['The Anderson-Darling rank test is not significant: ' ...
            'there is no evidence against MCAR.'];
    end
else
    if isempty(ADpvalue)
        txt = ['Hawkins'' test is not significant: there is no evidence ' ...
            'against multivariate normality or MCAR.'];
    elseif median(ADpvalue)<alpha
        txt = ['The Hawkins and Anderson-Darling results indicate a ' ...
            'departure from MCAR within the MCAR-versus-MAR framework.'];
    else
        txt = ['Hawkins'' test is significant but the nonparametric ' ...
            'Anderson-Darling test is not. This suggests nonnormality, ' ...
            'with no evidence against MCAR.'];
    end
end
end

% -------------------------------------------------------------------------
function printSummary(out)
%printSummary Display a concise summary in the Command Window.

fprintf('\nJamshidian-Jalal MCAR test\n')
fprintf('Retained patterns: %d\n',size(out.patterns,1))
fprintf('Cases used: %d of %d\n',out.nused,out.n)

if any(out.removedRows)
    fprintf('Cases removed because of sparse patterns: %d\n', ...
        sum(out.removedRows))
end

if ~isempty(out.hawkinsP)
    fprintf(['Hawkins: median chi-square = %.6g, df = %d, ' ...
        'median p-value = %.6g\n'],median(out.hawkinsStat), ...
        out.hawkinsDF,median(out.hawkinsP))

    if any(out.hawkinsP<out.alpha) && median(out.hawkinsP)>=out.alpha
        fprintf('Warning: at least one Hawkins p-value is below alpha.\n')
    end
end

if ~isempty(out.ADpvalue)
    fprintf(['Anderson-Darling: median statistic = %.6g, ' ...
        'median p-value = %.6g\n'],median(out.ADstat), ...
        median(out.ADpvalue))

    if any(out.ADpvalue<out.alpha) && median(out.ADpvalue)>=out.alpha
        fprintf(['Warning: at least one Anderson-Darling p-value is ' ...
            'below alpha.\n'])
    end
end

fprintf('%s\n',out.interpretation)
fprintf(['Caution: the procedure does not distinguish MCAR from MNAR; ' ...
    'a nonsignificant result does not rule out MNAR.\n\n'])
end

% -------------------------------------------------------------------------
function plotResults(out)
%plotResults Display p-value histograms and retained patterns.

ntests = (~isempty(out.hawkinsP))+(~isempty(out.ADpvalue));

if ntests>0
    figure('Name','mdJJtest p-values','Color','w')
    tiledlayout(ntests,1,'TileSpacing','compact','Padding','compact')

    if ~isempty(out.hawkinsP)
        nexttile
        histogram(out.hawkinsP,'BinLimits',[0 1])
        xline(out.alpha,'--')
        xlabel('Hawkins p-values')
        ylabel('Frequency')
        title(sprintf('%.1f%% below alpha', ...
            100*mean(out.hawkinsP<out.alpha)))
    end

    if ~isempty(out.ADpvalue)
        nexttile
        histogram(out.ADpvalue,'BinLimits',[0 1])
        xline(out.alpha,'--')
        xlabel('Anderson-Darling p-values')
        ylabel('Frequency')
        title(sprintf('%.1f%% below alpha', ...
            100*mean(out.ADpvalue<out.alpha)))
    end
end

figure('Name','Retained missingness patterns','Color','w')
imagesc(~out.patterns)
colormap([0.85 0.25 0.20; 0.20 0.50 0.80])
caxis([0 1])
cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = {'Missing','Observed'};
xlabel('Variable')
ylabel('Missingness pattern')
xticks(1:out.p)
varnames=cellstr("Y"+(1:out.p));
xticklabels(varnames)
yticks(1:size(out.patterns,1))
labels = arrayfun(@(j) sprintf('%d (n=%d)',j, ...
    out.patternInfo.Frequency(j)),1:size(out.patterns,1), ...
    'UniformOutput',false);
yticklabels(labels)
title('Patterns retained for the MCAR test')
end

% -------------------------------------------------------------------------
function ranges = columnRangeIgnoringNaN(Y)
%columnRangeIgnoringNaN Range of each variable over its observed values.

ranges = zeros(1,size(Y,2));
for j = 1:size(Y,2)
    values = Y(~isnan(Y(:,j)),j);
    ranges(j) = max(values)-min(values);
end
end

%FScategory:MULT-MissingData
