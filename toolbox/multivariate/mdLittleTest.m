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
%    Y :           Input data. Array or table. n x p data matrix; n observations and
%                  p variables possibly containing missing values (NaN's).
%                  Rows of Y represent observations, and columns represent
%                  variables.
%                  Data Types - single | double
%
%  Optional input arguments:
%
%    emOut :       Structure containing EM estimates.
%                  If supplied, it must contain fields:
%                  emOut.loc   : p x 1 estimated mean vector
%                  emOut.cov: p x p estimated covariance matrix
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
%    plots :       Produce a graphical summary of the missingness patterns.
%                  Boolean. If true, the patterns having the largest
%                  percentage contributions to Little's statistic are
%                  displayed in decreasing order. The percentage for
%                  pattern $r$ is $100T_r/T$, where $T_r$ is its individual
%                  contribution and $T$ is Little's overall statistic.
%                  At most the ten most important informative patterns are
%                  shown. The default is false.
%                  Example - 'plots',true
%                  Data Types - logical | single | double
%
%    tol :         Convergence tolerance for mdEM.
%                  Used only if option 'emOut' is empty.
%                  Default is 1e-7.
%                  Example - 'tol',1e-6
%                  Data Types - single | double
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
%      out.loc          = Global MLE of the mean vector (from EM).
%      out.cov          = Global MLE of the covariance matrix (from EM).
%      out.npatterns   = Number of distinct missingness patterns.
%      out.patterns    = R x p logical matrix. Each row is a distinct
%                        missingness pattern; true means observed entry.
%                        R is the number of distinct missingness patterns.
%      out.patternInfo = Table with one row for each distinct missingness
%                        pattern. Row names encode the observed variables:
%                        for example, with five variables, row name 11000
%                        indicates that only variables 1 and 2 are observed.
%                        The columns are:
%                        nPattern    = number of observations in the pattern;
%                        pObserved   = number of observed variables;
%                        Contribution= contribution to Little's statistic;
%                        CondSigma   = condition number of the covariance
%                                      submatrix used for that pattern.
%
%
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
%      T
%      =
%      \sum_{r=1}^{R}
%      n_r
%      (\bar{y}_r-\widehat{\mu}_r)^{\mathsf T}
%      \widehat{\Sigma}_r^{-1}
%      (\bar{y}_r-\widehat{\mu}_r).
%  \]
%
%  Under the null hypothesis of MCAR, T is asymptotically distributed as a
%  chi-square random variable with degrees of freedom
%
%  \[
%      \mathrm{df}
%      =
%      \sum_{r=1}^{R}p_r-p.
%  \]
%
%  where $p$ is the total number of variables.
%
%  If option plots is true, patterns are ranked according to their individual
%  contributions to Little's statistic. Thus, a pattern is regarded as more
%  important when
%
%  \[
%      T_r
%      =
%      n_r
%      (\bar{y}_r-\widehat{\mu}_r)^{\mathsf T}
%      \widehat{\Sigma}_r^{-1}
%      (\bar{y}_r-\widehat{\mu}_r)
%  \]
%
%  is large. The plot therefore identifies the missingness patterns that
%  contribute most strongly to the overall evidence against MCAR.
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

%{
    %% Example 5: Graphical identification of influential patterns.
    rng(6);
    Y = randn(300,5);

    Y(1:50,1) = NaN;
    Y(51:100,2) = NaN;
    Y(101:140,[1 3]) = NaN;
    Y(141:180,[2 4]) = NaN;
    Y(181:210,[3 5]) = NaN;

    out = mdLittleTest(Y,'plots',true);
%}

%% Beginning of code

% Input parameters checking
if nargin < 1
    error('FSDA:mdLittleTest:TooFewInputs', ...
        'At least one input argument is required.');
end

p = size(Y,2);

if istable(Y)
    Y=Y{:,:};
end

if ~ismatrix(Y) || ~isnumeric(Y)
    error('FSDA:mdLittleTest:WrongInput', ...
        'Input argument Y must be a numeric matrix.');
end


% Default options
maxiter   = 200;
tol       = 1e-7;
emOut=[];
plots=false;

% Optional arguments
if ~isempty(varargin)

options = struct('maxiter',maxiter,'tol',tol, ...
    'emOut',emOut,'plots',plots);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:mdLittleTest:WrongInputOpt', ...
                ['Number of supplied options is invalid. ' ...
                'Values may be missing.']);
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
    emOut=options.emOut;
    plots = options.plots;


    % Basic checks on options
    if ~isscalar(maxiter) || maxiter <= 0
        error('FSDA:mdLittleTest:WrongMaxiter', ...
            'Option ''maxiter'' must be a positive scalar.');
    end

    if ~isscalar(tol) || tol <= 0
        error('FSDA:mdLittleTest:WrongTol', ...
            'Option ''tol'' must be a positive scalar.');
    end

    if ~(isscalar(plots) && (islogical(plots) || isnumeric(plots)))
        error('FSDA:mdLittleTest:WrongPlots', ...
            'Option ''plots'' must be a logical or numeric scalar.');
    end
    plots = logical(plots);

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
    % nr = number of units with pattern r
    nr = sum(idxr);
    obsr = Patterns(r,:);
    % number of observed variables in pattern r
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
    pvalue = chi2cdf(stat,df,'upper');
end


rowNames = string(double(Patterns));              % convert to string matrix
rowNames = join(rowNames,"",2);          % join columns with space

patternInfo = table( nrVec, prVec, termVec, condVec, ...
    'VariableNames', {'nPattern','pObserved','Contribution','CondSigma'}, ...
    'RowNames',rowNames);

% Produce a graphical summary of the patterns that contribute most to
% Little's statistic.
if plots

% Produce a graphical summary of the patterns that contribute most to
% Little's statistic.
if plots

    % Completely missing patterns have no contribution to Little's
    % statistic and are therefore excluded from the graphical ranking.
    informative = patternInfo.pObserved > 0;

    contribution = patternInfo.Contribution(informative);
    nPattern = patternInfo.nPattern(informative);
    pObserved = patternInfo.pObserved(informative);
    patternNames = string(patternInfo.Properties.RowNames(informative));

    if ~isempty(contribution)

        % Express each contribution as a percentage of Little's statistic.
        if stat > 0
            percentageContribution = 100*contribution/stat;
        else
            percentageContribution = zeros(size(contribution));
        end

        % Rank patterns according to their contribution to Little's
        % statistic. Display at most the ten largest contributions.
        [percentageSorted,ord] = sort(percentageContribution,'descend');

        nToPlot = min(10,numel(ord));
        ord = ord(1:nToPlot);

        percentagePlot = percentageSorted(1:nToPlot);
        patternPlot = patternNames(ord);
        nPatternPlot = nPattern(ord);
        pObservedPlot = pObserved(ord);

        figure('Name','mdLittleTest: pattern contributions','Color','w');

        barh(percentagePlot)

        ax = gca;
        ax.YTick = 1:nToPlot;
        ax.YTickLabel = patternPlot;
        ax.YDir = 'reverse';

        xlabel('Percentage contribution to Little''s statistic')
        ylabel('Missingness pattern')

        title({ ...
            'Most important missingness patterns', ...
            sprintf('Little''s statistic = %.4g; p-value = %.4g', ...
            stat,pvalue)})

        grid on
        box on

        % Add the percentage contribution, number of observations and
        % number of observed variables to each displayed pattern.
        xMaxPlot = max(percentagePlot);
        if xMaxPlot <= 0
            xMaxPlot = 1;
        end

        dx = 0.015*xMaxPlot;

        patternDescription = compose( ...
            '%.1f%%   n = %d, p_{obs} = %d', ...
            percentagePlot,nPatternPlot,pObservedPlot);

        text(percentagePlot+dx,1:nToPlot,patternDescription, ...
            'VerticalAlignment','middle', ...
            'HorizontalAlignment','left')

        xlim([0 max(percentagePlot)+0.35*xMaxPlot])
    end
end
end

% Store output
out = struct;
out.stat        = stat;
out.df          = df;
out.pvalue      = pvalue;
out.loc          = muhat;
out.cov       = Sigmahat;
out.patterns    = Patterns;
out.npatterns   = nPatterns;
out.patternInfo = patternInfo;

end