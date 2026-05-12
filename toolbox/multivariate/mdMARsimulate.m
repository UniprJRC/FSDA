function [Ymar,out] = mdMARsimulate(Y,varargin)
%mdMARsimulate generates missing values under a MAR logistic mechanism.
%
%<a href="matlab: docsearchFS('mdMARsimulate')">Link to the help function</a>
%
%  This function introduces missing values in a data matrix using a Missing At
%  Random (MAR) mechanism. The probability that an entry is missing is modeled
%  through a logistic regression whose covariates are observed columns of X.
%
%  The logistic link is evaluated using the Statistics and Machine Learning
%  Toolbox probability distribution object created by
%
%      pdLogistic = makedist('Logistic','mu',0,'sigma',1);
%
%  Therefore, if $\eta_{ij} = \alpha_j + x_i' \beta_j$, the missingness probability
%  is computed as
%
%      $Pr(M_{ij}=1 | x_i)$ = cdf(pdLogistic,eta_ij),
%
%  which is equal to
%
%  \[
%      1/(1+exp(-\eta_{ij})).
%   \]
%
%  If obsCols=1 and missCols=2:p, the mechanism is
%
%   \[ 
%      Pr(M_{ij}=1 | X_{i1}) = cdf(pdLogistic, \alpha_j + \beta_j*X_{i1}), \qquad
%      j = 2, ..., p,
%    \]
%
%  where $M_{ij}=1$ denotes a missing entry. The intercept alpha_j is chosen so
%  that the expected missingness proportion in the corresponding missing
%  column is equal to missRate.
%
%  Required input arguments:
%
%    Y: Input data matrix. Matrix. n-by-p numeric matrix. Rows are
%        observations and columns are variables. Missing values already
%        present in Y are preserved. The columns used to drive the MAR
%        mechanism, specified by option obsCols, must contain finite values.
%
%  Optional input arguments:
%
%
%    obsCols : Observed columns driving the MAR mechanism. Vector.
%        Vector containing the indices of the fully observed variables used
%        in the logistic missingness model. The default value is 1.
%        Example - 'obsCols',[1 3]
%        Data Types - double
%
%    missCols : Columns in which missing values are generated. Vector.
%        Vector containing the column indices in which additional missing
%        values are introduced. The default value is 2:p.
%        Example - 'missCols',2:5
%        Data Types - double
%
%    missRate: Desired missingness proportion in missCols. Scalar or vector.
%        Scalar in the interval [0,1], or a row/column vector with one value
%        for each column specified in missCols. If missRate is scalar, the same
%        target missingness proportion is used for all columns in missCols.
%        The default value is 0.3.
%        Example - 'missRate',0.2
%        Data Types - double
%
%    beta : Logistic regression coefficients. Scalar, vector or matrix.
%        If beta is scalar, the same coefficient is used for all variables in
%        obsCols and missCols. If beta is a vector with length equal to
%        numel(obsCols), the same linear predictor is used for all columns in
%        missCols. If beta is a vector with length equal to numel(missCols)
%        and numel(obsCols)=1, a different coefficient is used for each
%        missing column. If beta is a matrix of size
%        numel(obsCols)-by-numel(missCols), each missing column has its own
%        logistic slope vector. The default value is 1.5.
%        Example - 'beta',1.5
%        Data Types - double
%
%    alphaInterval : Initial interval for alpha search. Vector.
%        Two-element vector used as the starting bracketing interval for
%        fzero. If the root is not bracketed, the interval is expanded
%        automatically. The default value is [-10 10].
%        Example - 'alphaInterval',[-20 20]
%        Data Types - double
%
%    plots : Plot missingness patterns. Boolean.
%        If plots is equal to 1 (true), a bar plot of the missingness patterns is
%        produced. If plots is equal to 0 (false), no plot is produced. The default
%        value is 0 (false).
%        Example - 'plots',true
%        Data Types - double
%
%    msg : Level of output to display. Boolean.
%        If msg is true, a compact summary of the target and obtained
%        missingness proportions is printed on the screen. The default value
%        is false.
%        Example - 'msg',true
%        Data Types - logical
%
%  Output:
%
%    Ymar : Data matrix with MAR missing values. Matrix.
%        n-by-p matrix equal to Y except for the additional NaN values
%        generated in the columns specified by missCols.
%
%    out : Structure which contains the following compact diagnostic fields:
%
%        out.alpha = 1-by-numel(missCols) vector containing the intercepts
%            used in the logistic missingness models.
%        out.beta = numel(obsCols)-by-numel(missCols) matrix containing the
%            logistic slope coefficients used for each missing column.
%        out.missRateTarget = 1-by-numel(missCols) vector containing the target
%            missingness proportions.
%        out.missRateGeneratedMissCols = scalar containing the proportion of
%            generated Bernoulli missing indicators in missCols.
%        out.missRateNewMissCols = scalar containing the proportion of newly
%            generated NaN values in missCols, excluding cells that were
%            already NaN in Y.
%        out.missRateObtainedMissCols = scalar containing the final proportion
%            of NaN values in missCols.
%        out.missRateObtainedAll = scalar containing the final proportion of
%            NaN values in the whole matrix Ymar.
%        out.patternTable = table summarizing missingness patterns, counts
%            and proportions. In the pattern strings, 0 means observed and 1
%            means missing.
%        out.obsCols = columns used to drive the MAR mechanism.
%        out.missCols = columns in which missing values were generated.
%        out.class = 'mdMARsimulate'.
%
%        Large diagnostic objects such as the full probability matrix, the
%        missingness masks and the LogisticDistribution object are not stored
%        in out, in order to keep the output structure compact.
%
%  See also: mdpattern, mdPartialMD, mdEM, mdTEM, mdLittleTest
%
%  References:
%
%  Little, R.J.A. and Rubin, D.B. (2019), "Statistical Analysis with Missing
%  Data", 3rd edition, Wiley.
%
%  Rubin, D.B. (1976), Inference and missing data, "Biometrika", Vol. 63,
%  pp. 581-592.
%
%  Copyright 2008-2026.
%  Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mdMARsimulate')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%  Examples:

%{
  %% MAR missingness driven by the first variable.
  rng(1708)
  n = 1000;
  p = 4;
  Y= randn(n,p);
  [Ymar,out] = mdMARsimulate(Y,'missRate',0.3,'beta',1.5, ...
      'obsCols',1,'missCols',2:p,'plots',1);
  mdpattern(Ymar)
  disp(out.patternTable)
%}

%{
  %% Different logistic slopes for different missing columns.
  rng(1708)
  Y= randn(1000,4);
  [Ymar,out] = mdMARsimulate(Y,'missRate',0.3,'obsCols',1, ...
      'missCols',2:4,'beta',[0.5 1.5 3],'msg',true);
  disp(out.alpha)
%}

%{
  %% MAR mechanism driven by two observed variables.
  rng(1708)
  Y= randn(1000,5);
  B = [1.2 0.5 1.0; -0.7 1.5 0.2];
  [Ymar,out] = mdMARsimulate(Y,'missRate',[0.2 0.3 0.4], ...
      'obsCols',[1 2],'missCols',3:5,'beta',B);
  disp(out.patternTable)
%}

%%  Beginning of code.

% Input parameters checking
if nargin < 1
    error('FSDA:mdMARsimulate:MissingInput','Input matrix Yis required.');
end

if ~isnumeric(Y) || ~ismatrix(Y)
    error('FSDA:mdMARsimulate:WrongInput','Y must be a numeric matrix.');
end

[n,p] = size(Y);

if p < 2
    error('FSDA:mdMARsimulate:WrongInput','Y must contain at least two columns.');
end

%% User options
options = struct;
options.missRate = 0.3;
options.obsCols = 1;
options.missCols = 2:p;
options.beta = 1.5;
options.alphaInterval = [-10 10];
options.plots = 0;
options.msg = false;

if ~isempty(varargin)

    [varargin{:}] = convertStringsToChars(varargin{:});

    UserOptions = varargin(1:2:length(varargin));

    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mdMARsimulate:WrongInputOpt', ...
            'Number of supplied options is invalid. Probably values for some parameters are missing.');
    end

    % Check if user options are valid options.
    aux.chkoptions(options,UserOptions)

    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
end

    missRate = options.missRate;
    obsCols = options.obsCols(:)';
    missCols = options.missCols(:)';
    beta = options.beta;
    alphaInterval = options.alphaInterval;
    plots = options.plots;
    msg = options.msg;


%% Further checks
if isempty(obsCols) || any(obsCols < 1) || any(obsCols > p) || any(obsCols ~= floor(obsCols))
    error('FSDA:mdMARsimulate:WrongInputOpt','obsCols must contain valid column indices of X.');
end

if isempty(missCols) || any(missCols < 1) || any(missCols > p) || any(missCols ~= floor(missCols))
    error('FSDA:mdMARsimulate:WrongInputOpt','missCols must contain valid column indices of X.');
end

if any(ismember(missCols,obsCols))
    error('FSDA:mdMARsimulate:WrongInputOpt','obsCols and missCols must be disjoint.');
end

finiteObs = isfinite(Y(:,obsCols));
if any(~finiteObs(:))
    error('FSDA:mdMARsimulate:WrongInput','Columns specified in obsCols must contain finite observed values.');
end

if ~isnumeric(missRate) || isempty(missRate) || any(missRate(:) < 0) || any(missRate(:) > 1)
    error('FSDA:mdMARsimulate:WrongInputOpt','missRate must contain values in the interval [0,1].');
end

nmiss = numel(missCols);
nobs = numel(obsCols);

if isscalar(missRate)
    missRate = repmat(missRate,1,nmiss);
else
    missRate = missRate(:)';
    if numel(missRate) ~= nmiss
        error('FSDA:mdMARsimulate:WrongInputOpt', ...
            'If missRate is not scalar, it must have one value for each column in missCols.');
    end
end

if ~isnumeric(beta) || isempty(beta) || any(~isfinite(beta(:)))
    error('FSDA:mdMARsimulate:WrongInputOpt','beta must be a finite numeric scalar, vector.');
end



if isscalar(beta)
    beta = repmat(beta,nobs,nmiss);
elseif isvector(beta) && numel(beta) == nobs
    beta = repmat(beta(:),1,nmiss);
elseif isvector(beta) && nobs == 1 && numel(beta) == nmiss
    beta = reshape(beta,1,nmiss);
elseif isequal(size(beta),[nobs nmiss])
    % beta already has the required size
else
    error('FSDA:mdMARsimulate:WrongInputOpt', ...
        ['beta must be scalar, a vector of length numel(obsCols), ', ...
        'a vector of length numel(missCols) when numel(obsCols)=1, ', ...
        'or a numel(obsCols)-by-numel(missCols) matrix.']);
end

if ~isnumeric(alphaInterval) || numel(alphaInterval) ~= 2 || alphaInterval(1) >= alphaInterval(2)
    error('FSDA:mdMARsimulate:WrongInputOpt', ...
        'alphaInterval must be a two-element increasing numeric vector.');
end

if ~(isscalar(plots) && (islogical(plots) || isnumeric(plots)) && any(double(plots) == [0 1]))
    error('FSDA:mdMARsimulate:WrongInputOpt','plots must be a logical or numeric scalar equal to 0 or 1.');
end

if ~(isscalar(msg) && (islogical(msg) || isnumeric(msg)))
    error('FSDA:mdMARsimulate:WrongInputOpt','msg must be a logical or numeric scalar.');
end
msg = logical(msg);


%% Create the standard logistic distribution object
% The CDF of this object is the inverse-logit link:
% cdf(pdLogistic,z) = 1/(1+exp(-z)).
pdLogistic = makedist('Logistic','mu',0,'sigma',1);

%% Compute alpha and missingness probabilities
Yobs = Y(:,obsCols);
alpha = zeros(1,nmiss);
probs = zeros(n,nmiss);

for j = 1:nmiss
    etaNoIntercept = Yobs * beta(:,j);
    targetj = missRate(j);

    if targetj == 0
        alpha(j) = -Inf;
        probs(:,j) = zeros(n,1);
    elseif targetj == 1
        alpha(j) = Inf;
        probs(:,j) = ones(n,1);
    else
        fun = @(a) mean(cdf(pdLogistic,a + etaNoIntercept)) - targetj;
        bracket = alphaInterval(:)';
        f1 = fun(bracket(1));
        f2 = fun(bracket(2));

        % Expand the bracket if necessary.
        kexp = 0;
        while f1*f2 > 0 && kexp < 50
            center = mean(bracket);
            halfWidth = diff(bracket);
            bracket = [center-halfWidth center+halfWidth];
            f1 = fun(bracket(1));
            f2 = fun(bracket(2));
            kexp = kexp + 1;
        end

        if f1*f2 > 0
            error('FSDA:mdMARsimulate:NoRoot', ...
                'Unable to bracket alpha for missing column %d. Try a wider alphaInterval.',missCols(j));
        end

        alpha(j) = fzero(fun,bracket);
        probs(:,j) = cdf(pdLogistic,alpha(j) + etaNoIntercept);
    end
end

%% Generate missingness indicators and apply them to Y
Ymar = Y;

% Each generated indicator is Bernoulli with probability probs(i,j).
% This is equivalent to rand(n,nmiss) < probs, but uses binornd from the
% Statistics and Machine Learning Toolbox and is closer to R's rbinom.
newMissSmall = logical(binornd(1,probs));

Xpart = Ymar(:,missCols);
oldMissSmall = isnan(Xpart);
Xpart(newMissSmall) = NaN;
Ymar(:,missCols) = Xpart;

newMissEffectiveSmall = newMissSmall & ~oldMissSmall;
missMask = isnan(Ymar);

%% Missingness patterns
[patterns,~,groupID] = unique(missMask,'rows');
counts = accumarray(groupID,1);
[countsSorted,ord] = sort(counts,'descend');
patternsSorted = patterns(ord,:);
patternStrings = cell(size(patternsSorted,1),1);

for i = 1:size(patternsSorted,1)
    patternStrings{i} = sprintf('%d',patternsSorted(i,:));
end

patternTable = table(patternStrings,countsSorted,countsSorted/n, ...
    'VariableNames',{'Pattern_0obs_1mis','Count','Proportion'});

%% Store compact output structure
out = struct;
out.alpha = alpha;
out.beta = beta;
out.missRateTarget = missRate;
out.missRateGeneratedMissCols = mean(newMissSmall(:));
out.missRateNewMissCols = mean(newMissEffectiveSmall(:));
tmpObtainedMissCols = isnan(Ymar(:,missCols));
out.missRateObtainedMissCols = mean(tmpObtainedMissCols(:));
out.missRateObtainedAll = mean(isnan(Ymar(:)));
out.patternTable = patternTable;
out.obsCols = obsCols;
out.missCols = missCols;
out.class = 'mdMARsimulate';

%% Display summary if requested
if msg
    fprintf('Target NA proportion in missCols:     %s\n',num2str(missRate));
    fprintf('Generated NA proportion in missCols:  %.4f\n',out.missRateGeneratedMissCols);
    fprintf('New NA proportion in missCols:        %.4f\n',out.missRateNewMissCols);
    fprintf('Obtained NA proportion in missCols:   %.4f\n',out.missRateObtainedMissCols);
    fprintf('Obtained NA proportion in all Ymar:   %.4f\n',out.missRateObtainedAll);
end

%% Plot missingness patterns if requested
if plots == 1
    figure
    bar(countsSorted)
    set(gca,'XTick',1:length(countsSorted),'XTickLabel',patternStrings)
    xlabel('Missingness pattern: 0 = observed, 1 = missing')
    ylabel('Frequency')
    title('Missingness patterns')
    grid on
end

end


%FScategory:MULT-MissingData
