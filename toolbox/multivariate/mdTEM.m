function out = mdTEM(Y, varargin)
% mdTEM  EM algorithm with trimming (TEM) for data with missing values.
%
%
%<a href="matlab: docsearchFS('mdTEM')">Link to the help function</a>
%
% The algorithm:
%  - At each iteration compute adjusted partial Mahalanobis distances;
%  - Rank them and set weights w_i = 1 for the lowest n*(1-alpha) rows,
%       else 0;
%  - Run E-step and M-step using these weights;
%  - Apply a consistency correction for the truncation induced by
%       trimming (see input option consistencyfactor);
%  - Repeat until convergence or maxiter.
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
%   method : method used to rescale the distances. String scalar or char vector.
%            Possible values are.
%
%            'pri'      = principled EM rescaling (default),
%                         d2_partial + (p - pobs).
%
%            'expScale' = expectation scaling,
%                         d2_partial * (p / pobs).
%
%            'zMap'     = standardization mapping,
%                         p + sqrt(2*p) * ((d2_partial - pobs) ./ sqrt(2*pobs)).
%
%            'detMap'   = determinant-based rescaling,
%                         d2_partial * (p / pobs) * (g_full / g_obs).
%
%            'chiMap'   = chi-square quantile mapping. Use the cdf and
%                         inverse of the cdf of Chi2 distribution.
%
%            'betaMap'  = Beta quantile mapping. Use the cdf and
%                         inverse of the cdf of Beta distribution.
%            'impMD'    = MD on EM-imputed data.
%            Example - 'method','chiMap'
%            Data Types - string scalar | char vector
%
%   consistencyfactor : treatment of the truncation bias of the scatter
%            estimate. Character vector or string scalar. Possible values:
%
%            'global'   = (default) single scalar factor
%                         k = (n/h)*F_{chi2_{p+2}}(a), a = chi2inv(h/n,p),
%                         applied to the whole scatter matrix. This is the
%                         complete-data Tallis factor evaluated at the full
%                         dimension p and at the global retained fraction
%                         h/n. It is exact only for the rows without
%                         missing entries.
%
%            'pattern'  = exact pattern-wise correction. Trimming on any
%                         strictly increasing adjustment induces, within
%                         each missingness pattern g, an exact radial
%                         truncation of the corresponding Gaussian
%                         marginal at the threshold
%                             a_g = phi_{p_g}^{-1}(c),
%                         where c is the h-th smallest adjusted distance
%                         and phi is the adjustment. Tallis's theorem then
%                         applies exactly with dimension p_g, giving
%                             gamma_g = F_{chi2_{p_g}}(a_g),
%                             k_g     = F_{chi2_{p_g+2}}(a_g)/gamma_g.
%                         The correction is applied to the data-driven part
%                         of the expected second moment only; the
%                         conditional-variance term Sigma_{m|o} is a
%                         model-based quantity and is left uncorrected.
%                         The location estimate needs no correction.
%
%            'weighted' = single scalar, but computed as the
%                         information-weighted average of the exact
%                         pattern-wise factors,
%                             kbar = sum_g h_g p_g k_g / sum_g h_g p_g,
%                         which reduces to 'global' when there are no
%                         missing values.
%
%            'none'     = no correction.
%
%            Example - 'consistencyfactor','pattern'
%            Data Types - char | string
%
%   condmeanimp :  Also give the matrix of conditional mean imputed values. Boolean.
%                 if true structure out also contains the matrix of imputed values called Yimp.
%                 The default value of condmeanimp is false.
%                 Example - 'condmeanimp',true
%                 Data Types - logical
%
%   stochimp :     Also give the matrix of stochastic imputed values. Boolean.
%                 if true structure out also contains the matrix of imputed values called stochYimp.
%                 The default value of stochimp is false.
%                 Example - 'stochimp',true
%                 Data Types - logical
%
%   storeobj :    Compute value of the objective function in each iteration.
%                 Boolean. If true structure out also contains the trimmed
%                 sum of the smallest adjusted distances in each iteration.
%                 The default value of storeobj is true.
%                 Example - 'storeobj',false
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
%              out.weights = n x 1 vector of final 0/1 trimming weights.
%              out.Yimp = empty value of matrix Y with imputed values
%                   (depending on input option condmeanimp)
%              out.stochYimp = empty value of matrix Y with imputed values
%                   (only if input option stochimp is true)
%              out.obj  = empty value or value of the objective function (trimmed sum of
%                   smallest MD) in each iteration
%                   (only if input option storeobj is true)
%              out.kfactor = scalar consistency factor actually applied to
%                   the whole scatter matrix ('global', 'weighted' and
%                   'none'), or the information-weighted average kbar of
%                   the pattern-wise factors ('pattern').
%              out.kinfo = table of the pattern-wise quantities at the last
%                   iteration, with variables pobs (p_g), nkept (h_g),
%                   athr (a_g), gammag (gamma_g) and kg (k_g). Rows with
%                   athr <= 0 correspond to patterns that the adjustment
%                   excludes entirely. Empty when consistencyfactor is
%                   'global' or 'none'.
%
%
% More About:
%
%   The pattern-wise correction is exact under the Gaussian model with
%   known parameters and MCAR missingness. With plug-in estimates it is a
%   plug-in quantity, exactly as in the complete-data case. Under MAR the
%   partial distances are not, in general, chi-squared given the pattern,
%   and both gamma_g and k_g become approximations.
%
%   Patterns for which a_g <= 0 cannot contribute any retained unit. This
%   occurs for the additive adjustment 'pri', for which the adjusted
%   distance is bounded below by p - p_g, when trimming is severe and p_g
%   is small. Such patterns are reported in out.kinfo.
%
%   Patterns with very few retained units give unstable k_g. Their factor
%   is shrunk towards kbar; see the local function localPatternFactors.
%
%
% See also: mdEM, mdImputeCondMean.m, mdPartialMD.m, mdPartialMD2full
%
% References:
%
% Little, R. J. A., & Rubin, D. B. (2019). Statistical Analysis with
% Missing Data (3rd ed.). Hoboken, NJ: John Wiley & Sons.
% Tallis, G. M. (1963). Elliptical and radial truncation in normal
% samples. Annals of Mathematical Statistics, 34, 940-944.
% van Buuren, S. (2018). Flexible Imputation of Missing Data (2nd ed.).
% Boca Raton, FL: Chapman & Hall/CRC (Taylor & Francis Group).
% Templ, M. (2023). Visualization and Imputation of Missing Values: With
% Applications in R. Cham, Switzerland: Springer Nature.
%
%
% Copyright 2008-2026.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('mdTEM')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Call to mdTEM with all the default options.
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
    out=mdTEM(Y);
    % Show true means and inputed means
    scatter(out.loc,muTrue)
    refline(1)
    xlabel('Imputed means')
    ylabel('True means')
%}

%{
    %% Example of use of option condmeanimp.
    % number of variables
    p = 15;
    % number of observations
    n = 1000;
    % target pairwise correlation (0<rho<1)
    rho = 0.9;
    % Covariance matrix (unit variances)
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);      % upper-triangular such that Sigma = R'*R
    % Generate samples ~ N(0,Sigma)
    Yfull = randn(n,p) * R;   % Strong positive correlation between the vars
    missRate = 0.25;     % MCAR missing probability per entry
    missMask = rand(n,p) < missRate;
    Y=Yfull;
    Y(missMask) = NaN;
    % md with missing imputation
    out=mdTEM(Y,'condmeanimp',true);
    % Mahalanobis distances using original matrix
    d2Ori=mahalFS(Yfull,mean(Yfull),cov(Yfull));
    % Calculate the Mahalanobis distance for the imputed data
    d2Imp = mahalFS(out.Yimp, mean(out.Yimp), cov(out.Yimp));
    % Compare original with distances for the imputed data
    scatter(d2Ori,d2Imp)
    xlabel('Original Mahalanobis Distances');
    ylabel('Imputed Mahalanobis Distances');
    grid on
%}

%{
    %% Exact pattern-wise consistency correction.
    % The single global factor is exact only for the complete rows. The
    % pattern-wise correction uses the exact factor for every missingness
    % pattern and removes most of the residual downward bias of the
    % scatter estimate.
    rng(7)
    p = 5;  n = 2000;  rho = 0.5;
    SigmaTrue = (1-rho)*eye(p) + rho*ones(p);
    Yfull = randn(n,p)*chol(SigmaTrue);
    Y = Yfull;
    Y(rand(n,p) < 0.35) = NaN;
    % rows left with no observed entry are put back
    allmiss = all(isnan(Y),2);
    Y(allmiss,1) = Yfull(allmiss,1);

    outG = mdTEM(Y,'alpha',0.5,'consistencyfactor','global');
    outP = mdTEM(Y,'alpha',0.5,'consistencyfactor','pattern');

    errG = max(max(abs(outG.cov-SigmaTrue)));
    errP = max(max(abs(outP.cov-SigmaTrue)));
    disp(['max |Sigma-hat - Sigma|, global  factor : ' num2str(errG)])
    disp(['max |Sigma-hat - Sigma|, pattern factors: ' num2str(errP)])
    % pattern-wise quantities of the last iteration
    disp(outP.kinfo)
%}

%{
    % Diagnostic: are incomplete rows retained at the right rate?
    % Under the null the share of incomplete rows in the retained subset
    % should match their share in the whole sample.
    rng(1)
    p = 7;  n = 10000;
    Y = randn(n,p);
    Y(rand(n,p) < 0.30) = NaN;
    allmiss = all(isnan(Y),2);  Y(allmiss,1) = randn(sum(allmiss),1);
    incomplete = any(isnan(Y),2);
    for cf = ["global" "pattern"]
        o = mdTEM(Y,'alpha',0.5,'consistencyfactor',cf);
        fprintf('%-8s share incomplete kept %.3f  (sample %.3f)\n', ...
            cf, mean(incomplete(o.weights==1)), mean(incomplete));
    end
%}

%% Beginning of code
alpha=0.5;
mus=[];
sigs=[];
maxiter=100;
tol=1e-5;
tol_sigma=true;
condmeanimp=false;
stochimp=false;
storeobj=true;
method='pri';
consistencyfactor='global';

if nargin>1
    options=struct('storeobj',storeobj,'alpha',alpha,'mus',mus,'sigs',sigs,'maxiter',maxiter,'tol',tol, ...
        'tol_sigma',tol_sigma,'condmeanimp',condmeanimp,'stochimp',stochimp,'method',method, ...
        'consistencyfactor',consistencyfactor);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:mdTEM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    condmeanimp=options.condmeanimp;
    stochimp=options.stochimp;
    method=string(options.method);
    storeobj=options.storeobj;
    consistencyfactor=char(string(options.consistencyfactor));
end

validCF={'global','pattern','weighted','none'};
if ~any(strcmpi(consistencyfactor,validCF))
    error('FSDA:mdTEM:WrongInputOpt', ...
        ['Option ''consistencyfactor'' must be one of ''global'', ', ...
         '''pattern'', ''weighted'' or ''none''.']);
end
consistencyfactor=lower(consistencyfactor);

if storeobj==true
    obj=zeros(maxiter,1);
else
    obj=[];
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

% number to keep:
keep_count = max(0, floor(n * (1 - alpha)));

% missingness mask and patterns are intrinsic to the data: computed once
nanY = isnan(Y);

w = zeros(n,1);
kinfo = [];
kfactor = 1;

while (dif > tol) && (iter < maxiter)
    iter = iter + 1;
    mus_old = mus;
    sigs_old = sigs;


    if method=="impMD"
        Yimp=mdImputeCondMean(Y, mus, sigs);
        % In this case compute Mahalanobis distances on imputed data
        d2_adj=mahalFS(Yimp,mus',sigs);
        poss = sum(~nanY,2);

    elseif method=="detMap"
        [d2, poss] = mdPartialMD(Y, mus, sigs);
        d2_adj = mdPartialMD2full(d2, p, poss,'method',method,'Y',Y,'Sigma',sigs);

    else
        % Trimming step: compute adjusted partial Mahalanobis distances
        [d2, poss] = mdPartialMD(Y, mus, sigs);
        d2_adj = mdPartialMD2full(d2, p, poss,'method',method);
    end

    % rank and select the smallest n*(1-alpha)
    % We treat NaN distances as large (so they're trimmed)
    nan_mask = isnan(d2_adj);

    % find indices of smallest distances
    % create sorted index from available (non-NaN) adj distances
    [~, idx_sorted] = sort(d2_adj, 'ascend', 'MissingPlacement', 'last');
    keep_idx = idx_sorted(1:min(keep_count, sum(~nan_mask)));


    w = zeros(n,1);
    w(keep_idx) = 1;
    mm = sum(w);

    % Common threshold implied by the concentration step: the largest
    % retained adjusted distance. Every retention event is equivalent to
    % d2_partial <= phi_{p_g}^{-1}(c) within each pattern.
    cthr = max(d2_adj(keep_idx));

    switch consistencyfactor

        case 'pattern'
            % Exact pattern-wise correction. The corrected accumulation
            % replaces the standard E- and M-steps: the data-driven part
            % of each pattern's expected second moment is divided by k_g,
            % the conditional-variance part is not, and the location
            % update is left uncorrected.
            kinfo = localPatternFactors(nanY, w, poss, cthr, p, n, ...
                sigs, method);
            [mus, sigs, kfactor] = localCorrectedStep(Y, nanY, w, mus, ...
                sigs, kinfo);

        otherwise
            % Standard E-step and M-step, then a scalar correction
            [T1, T2] = aux.NAcompute_expected_stats(Y, mus, sigs, w);
            [mus, sigs] = aux.NAmaximization_step(T1, T2, w);

            switch consistencyfactor
                case 'global'
                    a = chi2inv(mm/n,p);
                    kfactor = (n./mm).*(chi2cdf(a,p+2));
                    kinfo = [];
                case 'weighted'
                    kinfo = localPatternFactors(nanY, w, poss, cthr, p, ...
                        n, sigs, method);
                    kfactor = localWeightedFactor(kinfo);
                case 'none'
                    kfactor = 1;
                    kinfo = [];
            end
            if isfinite(kfactor) && kfactor > 0
                sigs = sigs/kfactor;
            end
    end

    if storeobj==true
        obj(iter)=sum(d2_adj(keep_idx))/kfactor;
    end

    % convergence check
    mu_diff = max(abs(mus(:) - mus_old(:)));
    sigma_diff = max(abs(sigs(:) - sigs_old(:)));
    if tol_sigma
        dif = max(mu_diff, sigma_diff);
    else
        dif = mu_diff;
    end
end

%% EM imputation of missing values (conditional means or stochastic imputation)
if condmeanimp ==true
    Yimp = mdImputeCondMean(Y, mus, sigs);
else
    Yimp=[];
end

if stochimp == true
    stochYimp = mdImputeStochastic(Y, mus, sigs);
else
    stochYimp=[];
end

if storeobj==true
    obj=obj(1:iter);
end

out.loc = mus;
out.cov = sigs;
out.iter = iter;
out.weights = w;
out.Yimp=Yimp;
out.stochYimp=stochYimp;
out.obj=obj;
out.kfactor = kfactor;
out.kinfo = kinfo;

end

%% ------------------------------------------------------------------
function kinfo = localPatternFactors(nanY, w, poss, cthr, p, n, Sigma, method)
% Exact Tallis factors, one per missingness pattern.
%
% Trimming retains the units with adjusted distance <= cthr. Within
% pattern g this is equivalent to d2_partial <= a_g = phi_{p_g}^{-1}(cthr),
% a radial truncation of the Gaussian marginal in dimension p_g. Hence
%    gamma_g = F_{chi2_{p_g}}(a_g),
%    k_g     = F_{chi2_{p_g+2}}(a_g)/gamma_g.
% Patterns with a_g <= 0 cannot retain any unit; their factor is set to
% NaN and later replaced by kbar. Patterns with too few retained units
% give unstable k_g and are shrunk towards kbar.

keep = w > 0;
patt = unique(nanY(keep,:),'rows');
G = size(patt,1);

pobs   = zeros(G,1);
nkept  = zeros(G,1);
athr   = zeros(G,1);
gammag = zeros(G,1);
kg     = nan(G,1);

for g = 1:G
    pg = patt(g,:);
    rows = keep & all(nanY == pg,2);
    obs = ~pg;
    pobs(g)  = sum(obs);
    nkept(g) = sum(rows);

    if pobs(g) == 0
        athr(g) = 0; gammag(g) = 0; kg(g) = NaN;
        continue
    end

    athr(g) = localInvAdjust(cthr, pobs(g), p, n, Sigma, obs, method);

    if athr(g) <= 0 || ~isfinite(athr(g))
        gammag(g) = 0; kg(g) = NaN;
    else
        gammag(g) = chi2cdf(athr(g), pobs(g));
        if gammag(g) <= 0
            kg(g) = NaN;
        else
            kg(g) = chi2cdf(athr(g), pobs(g)+2)/gammag(g);
        end
    end
end

kinfo = table(pobs, nkept, athr, gammag, kg);

% Stabilization: patterns with too few retained units, or with an
% undefined factor, borrow the information-weighted average.
kbar = localWeightedFactor(kinfo);
bad = ~isfinite(kinfo.kg) | kinfo.kg <= 0 | kinfo.kg > 1 | ...
      kinfo.nkept < max(2, kinfo.pobs);
kinfo.kg(bad) = kbar;
end

%% ------------------------------------------------------------------
function kbar = localWeightedFactor(kinfo)
% Information-weighted average of the exact pattern-wise factors,
%    kbar = sum_g h_g p_g k_g / sum_g h_g p_g,
% computed over the patterns for which k_g is well defined.
ok = isfinite(kinfo.kg) & kinfo.kg > 0 & kinfo.kg <= 1 & kinfo.nkept > 0;
if ~any(ok)
    kbar = 1;
    return
end
wgt = kinfo.nkept(ok).*kinfo.pobs(ok);
kbar = sum(wgt.*kinfo.kg(ok))/sum(wgt);
if ~isfinite(kbar) || kbar <= 0
    kbar = 1;
end
end

%% ------------------------------------------------------------------
function a = localInvAdjust(c, pg, p, n, Sigma, obs, method)
% Inverse of the adjustment for a pattern with pg observed variables:
% the value a of the partial squared distance such that phi_{pg}(a) = c.

method = string(method);

switch method
    case "pri"
        a = c - (p - pg);

    case "expScale"
        a = c*pg/p;

    case "zMap"
        % c = p + sqrt(2p)*(a-pg)/sqrt(2*pg)
        a = pg + sqrt(pg/p)*(c - p);

    case "detMap"
        % c = a*(p/pg)*(|Sigma|^{1/p} / |Sigma_oo|^{1/pg})
        gfull = exp(localLogdetSPD(Sigma)/p);
        gobs  = exp(localLogdetSPD(Sigma(obs,obs))/pg);
        a = c*(pg/p)*(gobs/gfull);

    case "chiMap"
        u = chi2cdf(c, p);
        u = min(max(u, eps), 1-eps);
        a = chi2inv(u, pg);

    case "betaMap"
        cn = (n-1)^2/n;
        if n <= p+1 || n <= pg+1
            a = c;  % not defined; fall back on the identity
            return
        end
        u = min(max(c/cn, 0), 1-eps);
        al = betacdf(u, p/2, (n-p-1)/2);
        al = min(max(al, eps), 1-eps);
        a = cn*betainv(al, pg/2, (n-pg-1)/2);

    case "impMD"
        % With the same fitted parameters used for imputation and for the
        % distance, the imputed-data distance equals the partial distance,
        % so the adjustment is the identity.
        a = c;

    otherwise
        a = c;
end
end

%% ------------------------------------------------------------------
function [musNew, sigsNew, kbar] = localCorrectedStep(Y, nanY, w, mus, ...
    sigs, kinfo)
% Corrected E- and M-step, grouped by missingness pattern.
%
% For pattern g with observed set o and missing set m, let
%    A_g = Sigma_{mo} Sigma_{oo}^{-1},
%    C_g = Sigma_{mm} - A_g Sigma_{om},
%    Z_i = y_{i,o} - mu_o.
% The expected second moment of the completed row about mu is
%    B_g Z_i Z_i' B_g' + Q_g C_g Q_g',   B_g = P_g + Q_g A_g,
% of which only the first term is affected by truncation. Dividing that
% term by k_g and leaving C_g untouched restores consistency. The blocks
% are accumulated directly, without forming B_g.

p = size(Y,2);
mus = mus(:);
keep = w > 0;

S  = zeros(p,p);   % corrected second moment about the current mus
m1 = zeros(p,1);   % mean deviation from the current mus
h  = 0;

patt = unique(nanY(keep,:),'rows');

for g = 1:size(patt,1)
    pg = patt(g,:);
    rows = keep & all(nanY == pg,2);
    hg = sum(rows);
    if hg == 0, continue, end

    o = find(~pg);
    m = find(pg);
    if isempty(o), continue, end

    % pattern factor
    idx = find(kinfo.pobs == numel(o) & kinfo.nkept == hg, 1);
    if isempty(idx)
        kg = 1;
    else
        kg = kinfo.kg(idx);
    end
    if ~isfinite(kg) || kg <= 0, kg = 1; end

    Z   = Y(rows,o) - mus(o)';        % hg x pg
    Szz = (Z'*Z)/kg;                  % corrected data-driven part
    sZ  = sum(Z,1)';                  % pg x 1, NOT corrected

    S(o,o) = S(o,o) + Szz;
    m1(o)  = m1(o)  + sZ;

    if ~isempty(m)
        Soo = sigs(o,o);
        Ag  = sigs(m,o)/Soo;                       % (p-pg) x pg
        Cg  = sigs(m,m) - Ag*sigs(o,m);            % conditional variance
        Cg  = (Cg + Cg')/2;

        SzzAg = Szz*Ag';
        S(o,m) = S(o,m) + SzzAg;
        S(m,o) = S(m,o) + SzzAg';
        S(m,m) = S(m,m) + Ag*SzzAg + hg*Cg;
        m1(m)  = m1(m)  + Ag*sZ;
    end

    h = h + hg;
end

if h == 0
    musNew = mus; sigsNew = sigs; kbar = 1; return
end

m1 = m1/h;
S  = S/h;

musNew  = mus + m1;
sigsNew = S - (m1*m1');
sigsNew = (sigsNew + sigsNew')/2;

kbar = localWeightedFactor(kinfo);
end

%% ------------------------------------------------------------------
function val = localLogdetSPD(S)
S = (S + S')/2;
[R,flag] = chol(S);
if flag == 0
    val = 2*sum(log(diag(R)));
else
    val = log(max(det(S), realmin));
end
end

%FScategory:MULT-MissingData
