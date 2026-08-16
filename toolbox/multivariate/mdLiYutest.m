function out = mdLiYutest(Y, varargin)
%mdLiYutest Li-Yu nonparametric test for Missing Completely At Random (MCAR).
%
%<a href="matlab: docsearchFS('mdLiYutest')">Link to the help function</a>
%
%  mdLiYutest implements the nonparametric MCAR test proposed by Li and
%  Yu (2014). The procedure compares the distributions of the variables
%  jointly observed by pairs of missingness-pattern groups using the
%  energy-distance dissimilarity of Rizzo and Szekely (2010). The null
%  distribution of the resulting Q statistic is approximated by bootstrap
%  resampling from the complete cases while preserving the observed
%  missingness patterns and their frequencies.
%
%  A rejection implies that the observed-data distributions differ across
%  at least two comparable missingness-pattern groups and therefore that
%  MCAR cannot hold. Failure to reject does not prove MCAR, because
%  distributional differences may occur only in variables that are missing
%  within a pattern and are therefore not directly testable from the
%  observed data.
%
%  Required input arguments:
%
%    Y :           Input data. Matrix, table or timetable. n x p data matrix
%                  possibly containing missing values (NaN's). Rows of Y
%                  represent observations and columns represent variables.
%                  At least one complete case is required because the
%                  bootstrap distribution is generated from the empirical
%                  distribution of the complete cases.
%                  Data Types - single | double | table | timetable
%
%  Optional input arguments:
%
%      alpha :     Significance level. Scalar in the interval (0,1).
%                  The default value is 0.05.
%                  Example - 'alpha',0.01
%                  Data Types - single | double
%
%      nboot :     Number of bootstrap replications. Positive integer.
%                  In each replication n complete observations are sampled
%                  with replacement from the complete-case group and the
%                  original missingness-pattern frequencies are imposed on
%                  the resampled complete data. Li and Yu (2014) used 500
%                  bootstrap replications in their main simulation studies.
%                  The default value is 500.
%                  Example - 'nboot',1000
%                  Data Types - single | double
%
% standardize :     Standardize variables before computing distances. Boolean.
%                  If true, each variable is centered by the median and
%                  divided by the normalized median absolute deviation
%                  (MAD) computed using the complete cases only. The same
%                  location and scale are applied to all missingness-pattern
%                  groups and to the bootstrap resamples. The normalized MAD
%                  is 1.482602218505602 times the raw MAD. The default value
%                  is false, which reproduces the Li-Yu statistic on the
%                  original measurement scales.
%                  Example - 'standardize',true
%                  Data Types - logical | single | double
%
%        msg :     Display a summary of the test. Boolean.
%                  If true, the Q statistic, bootstrap critical value,
%                  bootstrap p-value and conclusion are displayed. The
%                  default value is true.
%                  Example - 'msg',false
%                  Data Types - logical | single | double
%
%      plots :     Produce a graphical summary. Boolean.
%                  If true, a histogram of the bootstrap null distribution
%                  of Q is shown together with vertical lines for the
%                  observed Q statistic and the bootstrap critical value.
%                  The default value is false.
%                  Example - 'plots',true
%                  Data Types - logical | single | double
%
%  Output:
%
%    out :         Structure containing the following fields:
%
%          out.stat          = observed Li-Yu Q test statistic.
%          out.pvalue        = empirical upper-tail bootstrap probability,
%                              computed as the proportion of bootstrap Q
%                              statistics not smaller than out.stat.
%          out.reject        = true if out.stat is larger than the empirical
%                              upper alpha critical value.
%          out.criticalValue = empirical upper alpha critical value obtained
%                              from the bootstrap distribution.
%          out.alpha         = nominal significance level.
%          out.nboot         = number of bootstrap replications.
%          out.standardize   = logical flag indicating whether robust
%                              standardization was applied.
%          out.center        = 1 x p vector containing the complete-case
%                              medians used for standardization. It is empty
%                              when standardize is false.
%          out.scale         = 1 x p vector containing the complete-case
%                              normalized MADs used for standardization. It
%                              is empty when standardize is false.
%          out.bootstrapStat = nboot x 1 vector containing the bootstrap Q
%                              statistics.
%          out.B             = between-pattern dissimilarity measure entering
%                              the numerator of Q.
%          out.W             = within-pattern dissimilarity measure entering
%                              the denominator of Q.
%          out.n             = number of observations.
%          out.p             = number of variables.
%          out.npatterns     = number of distinct missingness patterns.
%          out.patterns      = npatterns x p logical matrix containing the
%                              distinct observed-data patterns. A true entry
%                              means that the corresponding variable is
%                              observed. The complete-case pattern is placed
%                              in the first row.
%          out.patternCounts = npatterns x 1 vector containing the frequency
%                              of each missingness pattern.
%          out.patternMembership = n x 1 vector identifying the reordered
%                              missingness pattern of each observation.
%          out.patternInfo   = table with one row for each missingness
%                              pattern and columns nPattern, pObserved,
%                              pMissing and WContribution.
%          out.pairInfo      = table with one row for each pair of patterns
%                              sharing at least one observed variable. The
%                              table contains the two pattern indices, their
%                              sample sizes, the number and names of common
%                              observed variables, the energy dissimilarity
%                              and the contribution to B.
%          out.completePattern = index of the complete-case pattern in the
%                              returned ordering. It is always equal to 1.
%          out.nComplete     = number of complete cases.
%          out.completeFraction = fraction of complete cases in the data.
%          out.variableNames = variable names used in the output tables.
%          out.interpretation = concise interpretation of the result.
%
%  More About:
%
%  Let the incomplete data be partitioned into s distinct missingness-pattern
%  groups. Group i contains $n_i$ observations, and $o_i$ denotes the set of
%  variables observed in that group. For two groups $i$ and $j$, define
%
%  \[
%      o_{ij}=o_i\cap o_j.
%  \]
%
%  Only pairs for which $o_{ij}$ is nonempty can be compared directly from
%  the observed data. If $X$ and $Z$ are two random samples with $n_1$ and
%  $n_2$ observations, respectively, define
%
%  \[
%      g(X,Z)
%      =
%      \frac{1}{n_1n_2}
%      \sum_{a=1}^{n_1}\sum_{b=1}^{n_2}
%      \lVert x_a-z_b\rVert,
%  \]
%
%  and the sample energy dissimilarity
%
%  \[
%      d(X,Z)=2g(X,Z)-g(X,X)-g(Z,Z).
%  \]
%
%  For every comparable pair of missingness-pattern groups, mdLiYutest
%  computes this dissimilarity using only the variables in $o_{ij}$. The
%  overall between-pattern measure is
%
%  \[
%      B
%      =
%      \sum_{1\leq i<j\leq s,\;o_{ij}\neq\varnothing}
%      \frac{n_i n_j}{2n}
%      d(Y_{i,o_{ij}},Y_{j,o_{ij}}).
%  \]
%
%  The within-pattern measure is
%
%  \[
%      W
%      =
%      \frac{1}{2}
%      \sum_{i=1}^{s}
%      n_i g(Y_{i,o_i},Y_{i,o_i}).
%  \]
%
%  The Li-Yu test statistic is
%
%  \[
%      Q
%      =
%      \frac{B/(s-1)}{W/(n-s)}.
%  \]
%
%  Large values of $Q$ indicate that between-pattern distributional
%  differences are large relative to within-pattern dispersion. Li and Yu
%  (2014) show that the test is consistent against alternatives that induce
%  distributional differences in the observed variables and have finite
%  second moments.
%
%  Since the energy dissimilarity is based on Euclidean distances, the
%  statistic depends on the measurement scales of the variables. Li and Yu
%  (2014) do not prescribe a standardization step. Therefore, the default
%  standardize=false computes the statistic on the original scales. If
%  standardize=true, mdLiYutest uses the complete cases to compute, for each
%  variable j, the median $m_j$ and normalized MAD $s_j$,
%
%  \[
%      s_j
%      =
%      1.4826
%      \operatorname{median}_{i\in\mathcal C}
%      |Y_{ij}-m_j|,
%  \]
%
%  where $\mathcal C$ is the set of complete cases. All observations are
%  then transformed using the same robust location and scale,
%
%  \[
%      Y_{ij}^{*}=\frac{Y_{ij}-m_j}{s_j}.
%  \]
%
%  The transformation is common to all missingness-pattern groups; scales
%  are never estimated separately within a pattern because differences in
%  dispersion across patterns are part of the alternatives that the test is
%  intended to detect. The bootstrap resamples are generated from the same
%  standardized complete-case empirical distribution.
%
%  The null distribution of $Q$ is calibrated by the bootstrap procedure in
%  Section 3 of Li and Yu (2014). Let $Y_1$ denote the complete-case group.
%  In each bootstrap replication, $n$ rows are sampled with replacement from
%  $Y_1$. The first $n_1$ sampled rows are assigned the complete-case pattern,
%  the next $n_2$ rows are assigned the second pattern, and so on. Variables
%  missing in each original pattern are then deleted from the corresponding
%  bootstrap rows. Consequently, every bootstrap sample has exactly the
%  same missingness-pattern frequencies as the observed data.
%
%  The empirical upper alpha quantile of the bootstrap $Q$ statistics is used
%  as the critical value. The paper bases the formal decision on whether the
%  observed $Q$ exceeds this critical value. out.pvalue is supplied as a
%  convenient empirical upper-tail summary of the same bootstrap
%  distribution.
%
%  The bootstrap calibration requires a reasonable number of complete
%  cases. Li and Yu (2014) found in their simulations that very small
%  complete-case samples can make the test conservative. For $p=4$, about 50
%  complete cases performed well in the settings they studied, but this is
%  not a universal threshold and larger dimensions can require more complete
%  cases.
%
%  The test evaluates the null hypothesis that all distributions observable
%  in pairwise common variable subsets are equal across missingness-pattern
%  groups. MCAR implies this null hypothesis, so rejection rules out MCAR.
%  Nonrejection does not by itself establish MCAR because differences may be
%  confined to unobserved components. Under an additional MAR assumption,
%  Li and Yu discuss nonrejection as compatible with MCAR.
%
%  See also: mdLittletest, mdJJtest, mdMAARtest, mdpattern
%
%  References:
%
%  Li, J. and Yu, Y. (2014), "A Nonparametric Test of Missing Completely at
%  Random for Incomplete Multivariate Data", Psychometrika,
%  doi: 10.1007/s11336-014-9410-4.
%
%  Rizzo, M. L. and Szekely, G. J. (2010), "DISCO analysis: A nonparametric
%  extension of analysis of variance", The Annals of Applied Statistics,
%  Vol. 4, pp. 1034-1055.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdLiYutest')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% Example 1: Li-Yu MCAR test with default options.
    % Generate Gaussian data and impose missingness independently of the
    % data. The generated mechanism is MCAR.
    rng(1)
    n = 250;
    p = 4;
    Y = randn(n,p);
    miss = rand(n,p) < 0.10;
    miss(1:80,:) = false;
    Y(miss) = NaN;

    out = mdLiYutest(Y);
    disp(out.stat)
    disp(out.pvalue)
%}
%
%{
    %% Example 2: A non-MCAR alternative driven by an observed variable.
    % Variable 1 is always observed. Missingness in variable 2 is made much
    % more likely for small values of variable 1, creating different
    % observed-data distributions across missingness-pattern groups.
    rng(2)
    n = 400;
    Y = randn(n,4);
    probMiss = 1./(1+exp(1.5+2*Y(:,1)));
    miss2 = rand(n,1) < probMiss;
    Y(miss2,2) = NaN;

    out = mdLiYutest(Y,'nboot',500,'plots',true);
    disp(out.reject)
    disp(out.pairInfo)
%}
%
%{
    %% Example 3: Preserve variable names supplied through a table.
    rng(3)
    Y = array2table(randn(300,5), ...
        'VariableNames',{'Income','Age','Score1','Score2','Score3'});
    M = false(height(Y),width(Y));
    M(1:90,:) = false;
    M(91:end,3:5) = rand(height(Y)-90,3) < 0.12;
    Yarray = table2array(Y);
    Yarray(M) = NaN;
    Y = array2table(Yarray,'VariableNames',Y.Properties.VariableNames);

    out = mdLiYutest(Y,'nboot',250,'msg',false);
    disp(out.patternInfo)
%}
%
%{
    %% Example 4: Robust standardization based on complete cases.
    % Use a common normalized MAD for each variable, estimated using only
    % complete observations. This is useful when variables have very
    % different measurement scales.
    rng(4)
    n = 300;
    Y = [randn(n,1), 100*randn(n,1), 0.01*randn(n,1)];
    miss = rand(n,3) < 0.10;
    miss(1:100,:) = false;
    Y(miss) = NaN;

    out = mdLiYutest(Y,'standardize',true,'msg',false);
    disp(out.scale)
    disp(out.pvalue)
%}

%% Beginning of code

% Input parameters checking.
if nargin < 1
    error('FSDA:mdLiYutest:TooFewInputs', ...
        'At least one input argument is required.');
end

% Preserve variable names before converting tables to numeric arrays.
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
    error('FSDA:mdLiYutest:WrongInput', ...
        'Input argument Y must be a numeric matrix, table or timetable.');
end

Y = double(Y);
[n,p] = size(Y);

if n < 3 || p < 1
    error('FSDA:mdLiYutest:SmallData', ...
        'At least three observations and one variable are required.');
end

observedValues = Y(~isnan(Y));
if any(~isfinite(observedValues))
    error('FSDA:mdLiYutest:NonFiniteData', ...
        'Observed entries of Y must be finite. Missing values must be NaN.');
end

% Default options.
alpha = 0.05;
nboot = 500;
standardize = false;
msg = true;
plots = false;

% Optional arguments.
if ~isempty(varargin)
    options = struct('alpha',alpha,'nboot',nboot, ...
        'standardize',standardize,'msg',msg,'plots',plots);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));

    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mdLiYutest:WrongInputOpt', ...
            ['Number of supplied options is invalid. ' ...
            'Values may be missing.']);
    end

    if ~isempty(UserOptions)
        aux.chkoptions(options,UserOptions)
    end

    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end

    alpha = options.alpha;
    nboot = options.nboot;
    standardize = options.standardize;
    msg = options.msg;
    plots = options.plots;
end

% Validate options also when defaults are used.
if ~isscalar(alpha) || ~isnumeric(alpha) || alpha <= 0 || alpha >= 1
    error('FSDA:mdLiYutest:WrongAlpha', ...
        'Option alpha must be a numeric scalar in the interval (0,1).');
end

if ~isscalar(nboot) || ~isnumeric(nboot) || nboot < 1 || ...
        nboot ~= floor(nboot)
    error('FSDA:mdLiYutest:WrongNboot', ...
        'Option nboot must be a positive integer.');
end
if nboot < 20
    warning('FSDA:mdLiYutest:FewBootstrapReplications', ...
        ['Only %d bootstrap replications were requested. The empirical ' ...
        'upper-tail calibration can be very coarse.'],nboot);
end

if ~(isscalar(standardize) && (islogical(standardize) || isnumeric(standardize)))
    error('FSDA:mdLiYutest:WrongStandardize', ...
        'Option standardize must be a logical or numeric scalar.');
end

if isnumeric(standardize) && ~ismember(standardize,[0 1])
    error('FSDA:mdLiYutest:WrongStandardize', ...
        'Numeric option standardize must be equal to 0 or 1.');
end

if ~(isscalar(msg) && (islogical(msg) || isnumeric(msg)))
    error('FSDA:mdLiYutest:WrongMsg', ...
        'Option msg must be a logical or numeric scalar.');
end

if ~(isscalar(plots) && (islogical(plots) || isnumeric(plots)))
    error('FSDA:mdLiYutest:WrongPlots', ...
        'Option plots must be a logical or numeric scalar.');
end

standardize = logical(standardize);
msg = logical(msg);
plots = logical(plots);

% Identify distinct observed-data patterns. A true entry means observed.
obsMatrix = ~isnan(Y);
[Patterns0,~,idx0] = unique(obsMatrix,'rows','stable');
counts0 = accumarray(idx0,1,[size(Patterns0,1) 1]);

% The Li-Yu bootstrap requires a complete-case pattern.
complete0 = find(all(Patterns0,2),1);
if isempty(complete0)
    error('FSDA:mdLiYutest:NoCompleteCases', ...
        ['The Li-Yu bootstrap requires at least one complete case. ' ...
        'No row of Y is fully observed.']);
end

% Put the complete-case pattern first, as in the notation of Li and Yu.
other = setdiff((1:size(Patterns0,1))',complete0,'stable');
order = [complete0; other];
Patterns = Patterns0(order,:);
patternCounts = counts0(order);

% Map old pattern numbers to the reordered pattern numbers.
oldToNew = zeros(size(Patterns0,1),1);
oldToNew(order) = 1:numel(order);
idxPatterns = oldToNew(idx0);

s = size(Patterns,1);
if s < 2
    error('FSDA:mdLiYutest:NoMissingValues', ...
        ['Only the complete-case pattern is present. The Li-Yu MCAR test ' ...
        'requires at least two distinct missingness patterns.']);
end

if n <= s
    error('FSDA:mdLiYutest:NoWithinPatternDf', ...
        ['The statistic W/(n-s) is undefined because n is not larger ' ...
        'than the number of missingness patterns.']);
end

nComplete = patternCounts(1);
if nComplete < 50
    warning('FSDA:mdLiYutest:FewCompleteCases', ...
        ['Only %d complete cases are available. Li and Yu (2014) show ' ...
        'that a small complete-case sample can make the bootstrap test ' ...
        'conservative. Their p=4 simulations performed well with about ' ...
        '50 complete cases, but the required size depends on dimension ' ...
        'and distribution.'],nComplete);
end

% Optionally standardize all variables using a single robust location and
% scale estimated from the complete cases. The same transformation is used
% for every missingness pattern and is therefore also inherited by the
% bootstrap resamples generated below.
if standardize
    completeRows = idxPatterns == 1;
    YcompleteOriginal = Y(completeRows,:);

    center = median(YcompleteOriginal,1);
    rawMAD = median(abs(YcompleteOriginal-center),1);
    scaleMAD = 1.482602218505602*rawMAD;

    if any(~isfinite(scaleMAD) | scaleMAD <= 0)
        badVariables = find(~isfinite(scaleMAD) | scaleMAD <= 0);
        error('FSDA:mdLiYutest:ZeroCompleteCaseMAD', ...
            ['The normalized MAD computed from the complete cases is zero ' ...
            'or undefined for variable(s) %s. Robust standardization ' ...
            'cannot be performed.'],strjoin(string(badVariables),', '));
    end

    Y = (Y-center)./scaleMAD;
else
    center = [];
    scaleMAD = [];
end

% Store the observations belonging to each pattern.
groups = cell(s,1);
for r = 1:s
    groups{r} = Y(idxPatterns==r,:);
end

% Compute the observed Q statistic and diagnostic contributions.
[Qobs,Bobs,Wobs,WContribution,pairDetail] = ...
    liYuStatistic(groups,Patterns,patternCounts,variableNames);

if ~isfinite(Qobs)
    error('FSDA:mdLiYutest:UndefinedStatistic', ...
        ['The Li-Yu Q statistic is not finite. This can occur when there ' ...
        'is no within-pattern dispersion or no informative comparison.']);
end

% Complete cases supply the empirical distribution used by the bootstrap.
Ycomplete = groups{1};

% Generate the bootstrap null distribution while preserving the original
% missingness-pattern frequencies exactly.
bootstrapStat = NaN(nboot,1);
for b = 1:nboot
    draw = randi(nComplete,n,1);
    YstarComplete = Ycomplete(draw,:);

    groupsStar = cell(s,1);
    first = 1;
    for r = 1:s
        last = first + patternCounts(r) - 1;
        Yr = YstarComplete(first:last,:);
        Yr(:,~Patterns(r,:)) = NaN;
        groupsStar{r} = Yr;
        first = last + 1;
    end

    bootstrapStat(b) = liYuStatistic(groupsStar,Patterns,patternCounts, ...
        variableNames);
end

if any(~isfinite(bootstrapStat))
    error('FSDA:mdLiYutest:UndefinedBootstrapStatistic', ...
        ['At least one bootstrap Q statistic is not finite. Check whether ' ...
        'the complete cases have sufficient within-sample variability.']);
end

% Empirical upper alpha quantile. This order-statistic definition uses the
% smallest observed bootstrap value whose empirical CDF is at least 1-alpha.
bootstrapSorted = sort(bootstrapStat);
criticalIndex = ceil((1-alpha)*nboot);
criticalIndex = max(1,min(nboot,criticalIndex));
criticalValue = bootstrapSorted(criticalIndex);

% The formal paper decision uses Q > criticalValue.
reject = Qobs > criticalValue;

% Convenient empirical upper-tail bootstrap probability. Li and Yu formulate
% the test through the critical value rather than a bootstrap p-value.
pvalue = mean(bootstrapStat >= Qobs);

% Pattern-level diagnostic table.
pObserved = sum(Patterns,2);
pMissing = p-pObserved;
patternNames = join(string(double(Patterns)),"",2);
patternInfo = table(patternCounts,pObserved,pMissing,WContribution, ...
    'VariableNames',{'nPattern','pObserved','pMissing','WContribution'}, ...
    'RowNames',cellstr(patternNames));

% Pairwise diagnostic table.
if isempty(pairDetail)
    pairInfo = table('Size',[0 8], ...
        'VariableTypes',{'double','double','double','double','double', ...
        'string','double','double'}, ...
        'VariableNames',{'Pattern1','Pattern2','n1','n2','pCommon', ...
        'CommonVariables','Distance','BContribution'});
else
    pairInfo = struct2table(pairDetail);
end

if reject
    interpretation = [ ...
        'The Li-Yu test rejects the observed-data homogeneity restrictions ' ...
        'implied by MCAR. Therefore the missingness mechanism is not MCAR.'];
else
    interpretation = [ ...
        'The Li-Yu test does not reject the observed-data homogeneity ' ...
        'restrictions implied by MCAR. This does not prove MCAR because ' ...
        'differences confined to unobserved components cannot be tested.'];
end

% Store output.
out = struct;
out.stat = Qobs;
out.pvalue = pvalue;
out.reject = reject;
out.criticalValue = criticalValue;
out.alpha = alpha;
out.nboot = nboot;
out.standardize = standardize;
out.center = center;
out.scale = scaleMAD;
out.bootstrapStat = bootstrapStat;
out.B = Bobs;
out.W = Wobs;
out.n = n;
out.p = p;
out.npatterns = s;
out.patterns = Patterns;
out.patternCounts = patternCounts;
out.patternMembership = idxPatterns;
out.patternInfo = patternInfo;
out.pairInfo = pairInfo;
out.completePattern = 1;
out.nComplete = nComplete;
out.completeFraction = nComplete/n;
out.variableNames = variableNames;
out.interpretation = interpretation;

if msg
    printLiYu(out)
end

if plots
    plotLiYu(out)
end

end

% -------------------------------------------------------------------------
function [Q,B,W,WContribution,pairDetail] = ...
    liYuStatistic(groups,Patterns,patternCounts,variableNames)
%liYuStatistic computes the Li-Yu B, W and Q statistics.
%
% Inputs:
%   groups        : cell array containing one n_i-by-p matrix per pattern.
%   Patterns      : s-by-p logical matrix; true means observed.
%   patternCounts : s-by-1 vector containing n_i.
%   variableNames : p-by-1 or 1-by-p string array of variable names.
%
% Outputs:
%   Q             : Li-Yu test statistic.
%   B             : between-pattern dissimilarity measure.
%   W             : within-pattern dissimilarity measure.
%   WContribution : s-by-1 vector of pattern contributions to W.
%   pairDetail    : structure array containing comparable-pair details.

s = size(Patterns,1);
n = sum(patternCounts);

wantWContribution = nargout >= 4;
wantPairDetail = nargout >= 5;

if wantWContribution
    WContribution = zeros(s,1);
else
    WContribution = [];
end
W = 0;

% W = 1/2 sum_i n_i g(Y_i,oi,Y_i,oi).
for i = 1:s
    obs = Patterns(i,:);
    ni = patternCounts(i);

    if any(obs) && ni > 1
        Xi = groups{i}(:,obs);
        contributionW = sum(pdist(Xi,'euclidean'))/ni;
    else
        contributionW = 0;
    end

    WContribution(i) = contributionW;
    W = W + contributionW;
end

% B = sum_{i<j,oij nonempty} n_i n_j/(2n) d_ij.
B = 0;
if wantPairDetail
    pairDetail = struct('Pattern1',{},'Pattern2',{},'n1',{},'n2',{}, ...
        'pCommon',{},'CommonVariables',{},'Distance',{},'BContribution',{});
else
    pairDetail = [];
end

k = 0;
nComparablePairs = 0;
for i = 1:s-1
    for j = i+1:s
        common = Patterns(i,:) & Patterns(j,:);
        if ~any(common)
            continue
        end

        Xi = groups{i}(:,common);
        Xj = groups{j}(:,common);

        gij = averageEuclideanDistance(Xi,Xj);
        gii = averageEuclideanDistance(Xi,Xi);
        gjj = averageEuclideanDistance(Xj,Xj);

        dij = 2*gij-gii-gjj;

        % The empirical energy distance is nonnegative in exact arithmetic.
        % Remove only tiny negative values caused by floating-point roundoff.
        scale = max([1,abs(2*gij),abs(gii),abs(gjj)]);
        if dij < 0 && abs(dij) <= 100*eps(scale)
            dij = 0;
        end

        contribution = patternCounts(i)*patternCounts(j)/(2*n)*dij;
        B = B + contribution;
        nComparablePairs = nComparablePairs + 1;

        if wantPairDetail
            k = k+1;
            pairDetail(k).Pattern1 = i;
            pairDetail(k).Pattern2 = j;
            pairDetail(k).n1 = patternCounts(i);
            pairDetail(k).n2 = patternCounts(j);
            pairDetail(k).pCommon = sum(common);
            pairDetail(k).CommonVariables = ...
                strjoin(variableNames(common),", ");
            pairDetail(k).Distance = dij;
            pairDetail(k).BContribution = contribution;
        end
    end
end

if nComparablePairs == 0
    Q = NaN;
    return
end

if W <= 0
    if B > 0
        Q = Inf;
    else
        Q = NaN;
    end
    return
end

Q = (B/(s-1))/(W/(n-s));

end

% -------------------------------------------------------------------------
function g = averageEuclideanDistance(X,Z)
%averageEuclideanDistance computes the mean of all pairwise Euclidean distances.
%
% Inputs:
%   X : n1-by-d numeric matrix.
%   Z : n2-by-d numeric matrix.
%
% Output:
%   g : (n1*n2)^(-1) times the sum of ||X_i-Z_j|| over all pairs.
%
% The computation is performed in blocks to avoid allocating a full
% n1-by-n2-by-d array. Diagonal zero distances are included when X and Z
% are the same sample, exactly as in the V-statistic definition used by
% Li and Yu (2014).

n1 = size(X,1);
n2 = size(Z,1);

if size(X,2) ~= size(Z,2)
    error('FSDA:mdLiYutest:DistanceDimensionMismatch', ...
        'The two samples used in an energy distance must have equal dimension.');
end

if size(X,2) == 0
    g = 0;
    return
end

% A 512-by-512 distance block requires about 2 MB in double precision.
blockSize = 512;
sumDistance = 0;

for i1 = 1:blockSize:n1
    i2 = min(i1+blockSize-1,n1);
    Xi = X(i1:i2,:);

    for j1 = 1:blockSize:n2
        j2 = min(j1+blockSize-1,n2);
        Zj = Z(j1:j2,:);

        % Accumulate squared coordinate differences directly. This is more
        % stable than the identity ||x-z||^2=||x||^2+||z||^2-2*x'*z
        % when observations contain large common offsets.
        squaredDistance = zeros(size(Xi,1),size(Zj,1));
        for h = 1:size(Xi,2)
            delta = bsxfun(@minus,Xi(:,h),Zj(:,h)');
            squaredDistance = squaredDistance + delta.^2;
        end

        distanceBlock = sqrt(squaredDistance);
        sumDistance = sumDistance + sum(distanceBlock(:));
    end
end

g = sumDistance/(n1*n2);

end

% -------------------------------------------------------------------------
function printLiYu(out)
%printLiYu displays a concise summary of the Li-Yu test.

fprintf('\nLi-Yu nonparametric MCAR test\n');
fprintf('Observations: %d; variables: %d; missingness patterns: %d\n', ...
    out.n,out.p,out.npatterns);
fprintf('Complete cases: %d (%.1f%%)\n', ...
    out.nComplete,100*out.completeFraction);
if out.standardize
    fprintf('Standardization: complete-case median / normalized MAD\n');
else
    fprintf('Standardization: none\n');
end
fprintf('Q statistic: %.6g\n',out.stat);
fprintf('Bootstrap critical value (alpha = %.4g): %.6g\n', ...
    out.alpha,out.criticalValue);
fprintf('Empirical bootstrap p-value: %.6g\n',out.pvalue);
fprintf('Reject MCAR restrictions: %s\n',char(string(out.reject)));
fprintf('%s\n\n',out.interpretation);

end

% -------------------------------------------------------------------------
function plotLiYu(out)
%plotLiYu plots the bootstrap null distribution and observed statistic.

figure('Name','mdLiYutest: bootstrap null distribution','Color','w');
histogram(out.bootstrapStat)
hold on
xline(out.criticalValue,'--', ...
    sprintf('critical value = %.4g',out.criticalValue), ...
    'LabelVerticalAlignment','middle');
xline(out.stat,'-', ...
    sprintf('observed Q = %.4g',out.stat), ...
    'LabelVerticalAlignment','bottom');

xlabel('Bootstrap Q statistic')
ylabel('Frequency')
title({ ...
    'Li-Yu nonparametric MCAR test', ...
    sprintf('B = %d bootstrap replications; empirical p-value = %.4g', ...
    out.nboot,out.pvalue)})
grid on
box on
hold off

end

%FScategory:MULT-MissingData
