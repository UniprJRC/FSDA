function [Xmar,out] = mdMARsimulate(X,varargin)
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
%  Therefore, if eta_ij = alpha_j + x_i'beta_j, the missingness probability
%  is computed as
%
%      Pr(M_ij=1 | x_i) = cdf(pdLogistic,eta_ij),
%
%  which is equal to
%
%      1/(1+exp(-eta_ij)).
%
%  If obsCols=1 and missCols=2:p, the mechanism is
%
%      Pr(M_ij=1 | X_i1) = cdf(pdLogistic,alpha_j + beta_j*X_i1),
%      j = 2, ..., p,
%
%  where M_ij=1 denotes a missing entry. The intercept alpha_j is chosen so
%  that the expected missingness proportion in the corresponding missing
%  column is equal to naProp.
%
%  Required input arguments:
%
%    X : Input data matrix. Matrix. n-by-p numeric matrix. Rows are
%        observations and columns are variables. Missing values already
%        present in X are preserved. The columns used to drive the MAR
%        mechanism, specified by option obsCols, must contain finite values.
%
%  Optional input arguments:
%
%    naProp : Desired missingness proportion. Scalar or vector.
%        Scalar in the interval [0,1], or a row/column vector with one value
%        for each column specified in missCols. If naProp is scalar, the same
%        target missingness proportion is used for all columns in missCols.
%        The default value is 0.3.
%        Example - 'naProp',0.2
%        Data Types - double
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
%    seed : Random seed. Empty value or scalar.
%        If seed is not empty, rng(seed) is called before generating the
%        missingness indicators. If seed is empty, the current random stream
%        is used. The default value is empty.
%        Example - 'seed',1708
%        Data Types - double
%
%    plots : Plot missingness patterns. Scalar.
%        If plots is equal to 1, a bar plot of the missingness patterns is
%        produced. If plots is equal to 0, no plot is produced. The default
%        value is 0.
%        Example - 'plots',1
%        Data Types - double
%
%    msg : Level of output to display. Boolean.
%        If msg is true, a compact summary of the target and obtained
%        missingness proportions is printed on the screen. The default value
%        is true.
%        Example - 'msg',false
%        Data Types - logical
%
%  Output:
%
%    Xmar : Data matrix with MAR missing values. Matrix.
%        n-by-p matrix equal to X, except for the additional NaN values
%        generated in the columns specified by missCols.
%
%    out : Structure which contains the following fields:
%
%        out.alpha = 1-by-numel(missCols) vector containing the intercepts
%            used in the logistic missingness models.
%        out.beta = numel(obsCols)-by-numel(missCols) matrix containing the
%            logistic slope coefficients used for each missing column.
%        out.probs = n-by-numel(missCols) matrix containing the missingness
%            probabilities used to generate the Bernoulli missingness
%            indicators. The probabilities are computed as
%            cdf(out.pdLogistic,alpha_j + X(:,obsCols)*beta_j).
%        out.pdLogistic = standard logistic distribution object created by
%            makedist('Logistic','mu',0,'sigma',1).
%        out.newMissMask = n-by-p logical matrix. True entries identify the
%            newly generated missing cells.
%        out.missMask = n-by-p logical matrix. True entries identify all NaN
%            values in Xmar.
%        out.naPropTarget = 1-by-numel(missCols) vector containing the target
%            missingness proportions.
%        out.naPropGeneratedMissCols = scalar containing the proportion of
%            generated missing indicators in missCols.
%        out.naPropNewMissCols = scalar containing the proportion of newly
%            generated NaN values in missCols, excluding cells that were
%            already NaN in X.
%        out.naPropObtainedMissCols = scalar containing the final proportion
%            of NaN values in missCols.
%        out.naPropObtainedAll = scalar containing the final proportion of
%            NaN values in the whole matrix Xmar.
%        out.patterns = matrix containing the observed missingness patterns.
%            0 means observed and 1 means missing.
%        out.counts = frequency of each missingness pattern.
%        out.patternTable = table summarizing patterns, counts and
%            proportions.
%        out.obsCols = columns used to drive the MAR mechanism.
%        out.missCols = columns in which missing values were generated.
%        out.class = 'mdMARsimulate'.
%
%  See also: makedist, cdf, binornd, fzero, rng, isnan
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
%
%{
   %% MAR missingness driven by the first variable.
   rng(1708)
   n = 1000;
   p = 4;
   X = randn(n,p);
   [Xmar,out] = mdMARsimulate(X,'naProp',0.3,'beta',1.5, ...
       'obsCols',1,'missCols',2:p,'plots',1);
   disp(out.patternTable)
%}
%
%{
   %% Different logistic slopes for different missing columns.
   rng(1708)
   X = randn(1000,4);
   [Xmar,out] = mdMARsimulate(X,'naProp',0.3,'obsCols',1, ...
       'missCols',2:4,'beta',[0.5 1.5 3],'msg',true);
   disp(out.alpha)
%}
%
%{
   %% MAR mechanism driven by two observed variables.
   rng(1708)
   X = randn(1000,5);
   B = [1.2 0.5 1.0; -0.7 1.5 0.2];
   [Xmar,out] = mdMARsimulate(X,'naProp',[0.2 0.3 0.4], ...
       'obsCols',[1 2],'missCols',3:5,'beta',B);
   disp(out.patternTable)
%}
%
%  Beginning of code.

%% Input parameters checking
if nargin < 1
    error('FSDA:mdMARsimulate:MissingInput','Input matrix X is required.');
end

if ~isnumeric(X) || ~ismatrix(X)
    error('FSDA:mdMARsimulate:WrongInput','X must be a numeric matrix.');
end

[n,p] = size(X);

if p < 2
    error('FSDA:mdMARsimulate:WrongInput','X must contain at least two columns.');
end


%% User options
options = struct;
options.naProp = 0.3;
options.obsCols = 1;
options.missCols = 2:p;
options.beta = 1.5;
options.alphaInterval = [-10 10];
options.seed = [];
options.plots = 0;
options.msg = true;

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

naProp = options.naProp;
obsCols = options.obsCols(:)';
missCols = options.missCols(:)';
beta = options.beta;
alphaInterval = options.alphaInterval;
seed = options.seed;
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

finiteObs = isfinite(X(:,obsCols));
if any(~finiteObs(:))
    error('FSDA:mdMARsimulate:WrongInput','Columns specified in obsCols must contain finite observed values.');
end

if ~isnumeric(naProp) || isempty(naProp) || any(naProp(:) < 0) || any(naProp(:) > 1)
    error('FSDA:mdMARsimulate:WrongInputOpt','naProp must contain values in the interval [0,1].');
end

nmiss = numel(missCols);
nobs = numel(obsCols);

if isscalar(naProp)
    naProp = repmat(naProp,1,nmiss);
else
    naProp = naProp(:)';
    if numel(naProp) ~= nmiss
        error('FSDA:mdMARsimulate:WrongInputOpt', ...
            'If naProp is not scalar, it must have one value for each column in missCols.');
    end
end

if ~isnumeric(beta) || isempty(beta)
    error('FSDA:mdMARsimulate:WrongInputOpt','beta must be a numeric scalar, vector or matrix.');
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

if ~(isscalar(plots) && isnumeric(plots) && any(plots == [0 1]))
    error('FSDA:mdMARsimulate:WrongInputOpt','plots must be a numeric scalar equal to 0 or 1.');
end

if ~(isscalar(msg) && (islogical(msg) || isnumeric(msg)))
    error('FSDA:mdMARsimulate:WrongInputOpt','msg must be a logical or numeric scalar.');
end
msg = logical(msg);

%% Initialize random generator if requested
if ~isempty(seed)
    if ~(isnumeric(seed) && isscalar(seed) && isfinite(seed))
        error('FSDA:mdMARsimulate:WrongInputOpt','seed must be empty or a finite numeric scalar.');
    end
    rng(seed)
end

%% Create the standard logistic distribution object
% The CDF of this object is the inverse-logit link:
% cdf(pdLogistic,z) = 1/(1+exp(-z)).
pdLogistic = makedist('Logistic','mu',0,'sigma',1);

%% Compute alpha and missingness probabilities
Xobs = X(:,obsCols);
alpha = zeros(1,nmiss);
probs = zeros(n,nmiss);

for j = 1:nmiss
    etaNoIntercept = Xobs * beta(:,j);
    targetj = naProp(j);

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

%% Generate missingness indicators and apply them to X
Xmar = X;

% Each generated indicator is Bernoulli with probability probs(i,j).
% This is equivalent to rand(n,nmiss) < probs, but uses binornd from the
% Statistics and Machine Learning Toolbox and is closer to R's rbinom.
newMissSmall = logical(binornd(1,probs));

Xpart = Xmar(:,missCols);
oldMissSmall = isnan(Xpart);
Xpart(newMissSmall) = NaN;
Xmar(:,missCols) = Xpart;

newMissMask = false(n,p);
newMissMask(:,missCols) = newMissSmall & ~oldMissSmall;
missMask = isnan(Xmar);

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

%% Store output structure
% TODO: reduce the output structure size
out = struct;
out.alpha = alpha;
out.beta = beta;
out.probs = probs;
out.pdLogistic = pdLogistic;
out.newMissMask = newMissMask;
out.missMask = missMask;
out.naPropTarget = naProp;
out.naPropGeneratedMissCols = mean(newMissSmall(:));
tmpNewMissCols = newMissMask(:,missCols);
out.naPropNewMissCols = mean(tmpNewMissCols(:));
tmpMissCols = isnan(Xmar(:,missCols));
out.naPropObtainedMissCols = mean(tmpMissCols(:));
out.naPropObtainedAll = mean(isnan(Xmar(:)));
out.patterns = patternsSorted;
out.counts = countsSorted;
out.patternTable = patternTable;
out.obsCols = obsCols;
out.missCols = missCols;
out.class = 'mdMARsimulate';

%% Display summary if requested
if msg
    fprintf('Target NA proportion in missCols:     %s\n',num2str(naProp));
    fprintf('Generated NA proportion in missCols:  %.4f\n',out.naPropGeneratedMissCols);
    fprintf('New NA proportion in missCols:        %.4f\n',out.naPropNewMissCols);
    fprintf('Obtained NA proportion in missCols:   %.4f\n',out.naPropObtainedMissCols);
    fprintf('Obtained NA proportion in all Xmar:   %.4f\n',out.naPropObtainedAll);
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