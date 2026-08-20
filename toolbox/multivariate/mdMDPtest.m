function out = mdMDPtest(Y, varargin)
%mdMDPtest Mahalanobis Distance Perturbation test for MCAR
%
%<a href="matlab: docsearchFS('mdMDPtest')">Link to the help function</a>
%
%  mdMDPtest implements the Mahalanobis Distance Perturbation (MDP) test
%  for assessing the Missing Completely At Random (MCAR) assumption. The
%  test is based on the change in Mahalanobis distances for the units
%  without missing values when location and scatter are estimated:
%
%    1) using only the complete rows;
%    2) using all rows through EM/TEM in the presence of missing values.
%
%  The bootstrap null hypothesis is that the observed perturbation is
%  compatible with MCAR. The null distribution is generated from a Gaussian
%  model fitted on the complete rows and then the observed missingness mask
%  is imposed on the generated data.
%
%
%  Required input arguments:
%
%    Y : Input data. Matrix. n x p data matrix possibly containing NaNs.
%        Rows of Y represent observations and columns represent variables.
%        Data Types - double
%
%
%  Optional input arguments:
%
%    alpha : Trimming/robustness level. Scalar in the interval [0,0.5].
%            The default value is 0.
%
%            If alpha=0, classical complete-case estimates are used and
%            mdEM is used for the incomplete-data fit.
%
%            If alpha>0 and robust='FS', the Forward Search on the
%            complete cases starts monitoring at
%
%                h = floor(nComplete*(1-alpha)).
%
%            where nComplete is the number of rows with complete
%            observations.
%
%            The final complete-case location and scatter are not based
%            necessarily on exactly h observations. FSM applies its
%            signal-detection and stopping rules and the final estimates
%            are computed from the observations not declared as outliers.
%
%            If alpha>0 and robust='MCD', alpha is used as the MCD
%            breakdown-point option. If robust='MM', alpha is used as the
%            breakdown point of the preliminary S estimator.
%
%            For all robust estimators, outliers are declared using a
%            Bonferronized simultaneous confidence level of 0.99. For FSM
%            two signal/stopping rules are available. With bonflev=0.99,
%            the direct 99 percent Bonferroni bound is used. With
%            bonflev=[], the standard FSM signal-detection, validation and
%            envelope-resuperimposition rules are used. For MCD and MM,
%            the individual confidence level is 1-0.01/nComplete.
%            Example - 'alpha', 0.1
%            Data Types - double inside [0 0.5]
%
%
%   method : Rescaling method used inside mdTEM and mdPartialMD2full.
%            The default value is 'pri'.
%            Example - 'method','betaMap'
%            Data Types - char | string
%
% consistencyfactor : TEM truncation-bias correction used when alpha>0.
%            Character vector or string scalar. Possible values are
%            'pattern', 'adaptive', 'global', 'weighted' or 'none'. The
%            default is 'pattern', which preserves the previous mdMDPtest
%            behavior and uses Gaussian pattern-wise Tallis factors.
%
%            With 'adaptive', mdTEM estimates the pattern-wise radial
%            consistency factors from the complete cases, using the
%            complete-case location and scatter fitted by the selected
%            robust estimator. Thus the all-data TEM scatter is calibrated
%            to the same elliptical scatter target as the complete-case
%            estimator without specifying the radial distribution. The
%            reference sample and reference fit are reconstructed inside
%            every bootstrap replicate.
%
%            The options 'global', 'weighted' and 'none' have the same
%            meaning as in mdTEM. This option is ignored when alpha=0.
%            The current coupledtrim sensitivity implementation supports
%            only consistencyfactor='pattern'.
%            Example - 'consistencyfactor','adaptive'
%            Data Types - char | string
%
% adaptivepool : Pool adaptive projected reference radii across patterns
%            having the same observed dimension. Logical scalar. The
%            default is false. This option is used only when alpha>0 and
%            consistencyfactor='adaptive'.
%            Example - 'adaptivepool',true
%            Data Types - logical
%
% adaptiveminref : Minimum number of complete-case reference radii inside
%            a pattern cutoff required by the adaptive correction. Positive
%            integer. The default is 20. This option is used only when
%            alpha>0 and consistencyfactor='adaptive'.
%            Example - 'adaptiveminref',10
%            Data Types - double
%
%   nsimul : Number of bootstrap simulations. Scalar integer.
%            The default value is 499.
%            Example - 'nsimul',999
%            Data Types - double
%
%  conflev : Confidence level used to compute bootstrap confidence intervals.
%            Scalar in the interval (0,1). The default value is 0.95.
%            Example - 'conflev',0.99
%            Data Types - double
%       
%    eps0 : Positive numerical stabilizer used in the log-distance ratio.
%           The statistic is
%
%             log((d2_all + eps0) ./ (d2_cc + eps0)).
%
%           The default value is 1e-12.
%           Example - 'eps0',1e-10
%           Data Types - double
%
% aggregation : Mean aggregation rule. Character or string. Possible values are
%           'ordinary', 'inlier' or 'fixed'. The default is 'ordinary'.
%
%           'ordinary' computes the two mean statistics using all complete
%           cases (MDP-U).
%
%           'inlier' computes the two mean statistics using only complete
%           cases not declared as outliers by the selected robust
%           complete-case estimator (MDP-I). This option requires alpha>0.
%
%           'fixed' computes the two mean statistics using only complete
%           cases whose complete-case squared Mahalanobis distance is not
%           larger than chi2inv(filterlev,p) (MDP-F).
%
%           The median statistics are not affected by this option. With
%           coupledtrim=false (default), the aggregation rule acts only on
%           the final two means and does not alter the EM/TEM fit. The
%           selected mean aggregation rule is reconstructed independently
%           in every bootstrap sample.
%           Example - 'aggregation','inlier'
%           Data Types - char | string
%
% filterlev : Probability level defining the fixed robust-distance gate for
%           aggregation='fixed'. Scalar in the interval (0,1). The cutoff
%           is c=chi2inv(filterlev,p). The default value is 0.99.
%           This option is ignored for the other aggregation rules.
%           Example - 'filterlev',0.975
%           Data Types - double
%
% coupledtrim : Couple complete-case outlier decisions to TEM. Logical.
%           The default is false. This option is available only when
%           alpha>0 and aggregation='inlier'. When true, complete cases
%           declared as outliers by the selected robust complete-case
%           estimator are forced to have zero weight in every TEM
%           concentration step. TEM still applies its own trimming rule to
%           the full data; the coupled weight is the product of the TEM
%           trimming indicator and the fixed complete-case eligibility
%           indicator. This mode is intended for sensitivity analysis and
%           is not covered by the analytical reference law returned for
%           the uncoupled procedure.
%           Example - 'coupledtrim',true
%           Data Types - logical
%
%  robust  : Robust estimator for the complete cases. Character, string
%            or structure. This option is used only when alpha is strictly
%            positive. If robust is a character vector or string scalar,
%            possible values are 'FS', 'MCD' or 'MM'. The default is 'FS'.
%            If robust is a structure, the following fields can be supplied:
%              robust.class = robust estimator. Possible values are 'FS',
%                             'MCD' or 'MM'. This field is required.
%              robust.eff   = nominal efficiency of the MM estimator. This
%                             field is used only if robust.class='MM'. The
%                             default value is 0.95.
%              robust.bonflev = Forward Search signal/stopping rule. This
%                             field is used only if robust.class='FS'.
%                             Use [] (default) for the standard FSM
%                             signal-detection, validation and envelope-
%                             resuperimposition rules, or 0.99 for the
%                             direct 99 percent Bonferroni-bound rule.
%            The robustness level is controlled by alpha. In the FS case,
%            floor(nComplete*(1-alpha)) is used as the FSM initial
%            monitoring step; in the MCD case alpha is passed as bdp; in
%            the MM case alpha is passed as Sbdp. For all three estimators,
%            the confidence level used to declare outliers is fixed at the
%            Bonferronized simultaneous level 0.99. In the FS case the user
%            can choose between the direct Bonferroni rule and the standard
%            envelope-resuperimposition rule through robust.bonflev.
%            Example - 'robust','MCD' | r=struct; r.class='FS'; r.bonflev=[]; 'robust',r
%            Data Types - char | string | struct
%
%      tol : Convergence tolerance passed to mdTEM. Scalar.
%            The default value is 1e-10.
%            Example - 'tol',1e-8
%            Data Types - double
%
%    plots : Flag to produce the output plot.
%            The default value is false.
%            Example - 'plots',true
%            Data Types - logical
%
%
%  Output:
%
%    out : Structure containing the following fields:
%
%          out.pvalue      = 1 x 4 vector containing the bootstrap p-values
%                            for the four statistics.
%          out.Tobs        = 1 x 4 vector containing the observed values of
%                            the four statistics.
%          out.Tboot       = nsimul x 4 matrix containing the bootstrap
%                            values of the four statistics.
%          out.alpha       = Value of input option alpha.
%          out.method      = Value of input option method.
%          out.consistencyfactor = TEM consistency-factor option actually
%                            requested. It is relevant only when alpha>0.
%          out.adaptivepool = Logical adaptive pooling option.
%          out.adaptiveminref = Minimum adaptive reference count.
%          out.TEMkfactor  = Scalar summary of the final TEM consistency
%                            correction. NaN for alpha=0.
%          out.TEMkinfo    = Final pattern-wise TEM consistency-factor
%                            diagnostics. Empty for alpha=0 or whenever the
%                            selected mdTEM correction does not return them.
%          out.aggregation = Mean aggregation rule actually used:
%                            'ordinary', 'inlier' or 'fixed'.
%          out.filterlev   = Value of input option filterlev.
%          out.filterCutoff = Fixed chi-square cutoff used by MDP-F. It is
%                            NaN unless aggregation='fixed'.
%          out.coupledtrim = Logical flag indicating whether estimator-
%                            declared complete-case outliers are forced to
%                            zero weight inside TEM.
%          out.coupledRows = Original row numbers forced to zero TEM weight.
%          out.nCoupledCC = Number of observed complete cases forced to zero
%                            TEM weight.
%          out.nCoupledBoot = Number of complete cases forced to zero TEM
%                            weight in each retained bootstrap sample.
%          out.TEMweights = Final TEM weights used for the observed data.
%          out.TEMinternalWeights = Final unconstrained TEM concentration
%                            indicators before multiplication by the fixed
%                            complete-case eligibility mask. For uncoupled
%                            fits this equals out.TEMweights.
%          out.meanWeightsCC = Logical vector of length out.nComplete. True
%                            entries identify complete cases entering the
%                            two mean statistics.
%          out.meanRows    = Original row numbers entering the two mean
%                            statistics.
%          out.nMeanCC     = Number of complete cases entering the two mean
%                            statistics.
%          out.nMeanBoot   = Number of complete cases entering the two mean
%                            statistics in each retained bootstrap sample.
%          out.nComplete   = Number of complete rows.
%          out.completeIdx = Logical index of complete rows.
%          out.locCC       = Complete-case estimate of location.
%          out.covCC       = Complete-case estimate of scatter.
%          out.robust      = Value of input option robust.
%          out.robustClass = Robust estimator actually used for complete
%                            cases: 'FS', 'MCD' or 'MM'.
%          out.robustEff   = MM nominal efficiency. This field is relevant
%                            only when out.robustClass='MM'.
%          out.robustBonflev = Value actually used for the FS
%                            signal/stopping rule. It is 0.99 for the
%                            direct Bonferroni rule and [] for the standard
%                            FSM envelope-resuperimposition rule. If fewer
%                            than four monitoring steps are available, []
%                            is automatically changed to 0.99.
%          out.FSrule      = Forward Search signal/stopping rule actually
%                            used: 'bonferroni' or 'resuperimposition'. It
%                            is empty if out.robustClass is not 'FS'.
%          out.outlierConflev = Simultaneous Bonferronized confidence level
%                            used to declare complete-case outliers (0.99).
%          out.outliersCC  = Complete-case units declared as outliers by
%                            the selected robust estimator.
%          out.nOutliersCC = Number of units in out.outliersCC.
%          out.d2_cc       = Mahalanobis distances computed from complete
%                            rows only.
%          out.d2_all      = Mahalanobis distances for the same complete
%                            rows when parameters are estimated from all
%                            the data through EM/TEM.
%          out.ciBoot      = 2 x 4 matrix containing the bootstrap
%                            confidence intervals for the four statistics.
%          out.loc          = Estimated location from EM/TEM fit on all data.
%          out.cov          = Estimated scatter from EM/TEM fit on all data.
%          out.eps0        = Numerical stabilizer used in the log-distance
%                            ratio.
%          out.asympt      = Structure containing the analytical Gaussian
%                            benchmark for the two selected mean statistics.
%                            For MDP-U the ordinary first-order coefficients
%                            are used. MDP-I inherits the same coefficients
%                            under the condition that the complete-case rule
%                            excludes only O_p(1) clean observations. For
%                            MDP-F the base estimator-discrepancy variance is
%                            multiplied by kc^2 for the mean distance
%                            difference and by lambda^2 for the stabilized
%                            mean log-ratio, where
%                              kc=F_{p+2}(c)/F_p(c)
%                            and lambda is the corresponding truncated radial
%                            coefficient. Fields aggregation, filterCutoff,
%                            gammaFilter, kc, lambda, baseSigmaD2 and
%                            baseKappa document the applied scaling. For
%                            alpha=0, the classical closed-form covariance
%                            information formula is used. For alpha>0,
%                            method='pri' and consistencyfactor='pattern',
%                            the Gaussian pattern-wise TEM influence function
%                            and effective Jacobian are evaluated analytically.
%                            Analytical calibration for consistencyfactor=
%                            'adaptive' is deliberately not returned in this
%                            implementation; its empirical-factor influence
%                            contribution will be added separately. The scalar
%                            influence of the robust
%                            complete-case scatter estimator is evaluated by
%                            a delete-one jackknife and combined with the TEM
%                            contribution in the end-to-end sandwich variance.
%                            The fields TDmean and TLmean contain the
%                            asymptotic variance of sqrt(n)T (sigma2), the
%                            standard error of T (se), the studentized value
%                            (z), and the two-sided asymptotic p-value. For
%                            alpha>0, out.asympt.TEM contains diagnostics for
%                            the analytical TEM contribution and
%                            out.asympt.completeCase documents the numerical
%                            influence evaluation. The adaptive FS branch is
%                            flagged as experimental because its stopping rule
%                            is outside the fixed-fraction FS asymptotic theorem.
%
%
%  More About:
%
%  Let d2_cc denote the squared Mahalanobis distances computed on the
%  complete rows using the complete-case estimates, and let d2_all denote
%  the distances for the same rows when location and scatter are estimated
%  from all the data using EM/TEM. The function monitors the following four
%  statistics:
%
%    1) median( log((d2_all + eps0) ./ (d2_cc + eps0)) );
%    2) selected mean of log((d2_all + eps0) ./ (d2_cc + eps0));
%    3) median( d2_all - d2_cc );
%    4) selected mean of d2_all - d2_cc.
%
%  The selected means use all complete cases for aggregation='ordinary',
%  the estimator-defined complete-case inliers for aggregation='inlier',
%  and the fixed robust-distance gate for aggregation='fixed'. The same
%  rule is re-estimated in every bootstrap sample. With coupledtrim=true,
%  estimator-declared complete-case outliers are additionally forced to zero
%  weight in every TEM concentration step.
%
%  When consistencyfactor='adaptive', the complete-case observations are
%  also the reference sample for the empirical radial correction. The
%  reference location and scatter are the same complete-case estimates used
%  to compute d2_cc. Consequently the adaptive TEM correction targets the
%  complete-case scatter functional under an elliptical MCAR model. The
%  mask-preserving parametric bootstrap remains Gaussian-model based; it
%  reruns both the robust complete-case fit and adaptive radial correction in
%  every bootstrap sample.
%
%  Small p-values indicate that the change in distances is larger than what
%  is expected under the MCAR bootstrap model.
%
%  The analytical Gaussian benchmark follows the selected mean aggregation.
%  MDP-I uses the same first-order variance as MDP-U under O_p(1) clean
%  exclusions. MDP-F multiplies the base mean-difference variance by kc^2
%  and the base stabilized-log-ratio variance by the appropriate truncated
%  radial coefficient lambda^2. When coupledtrim=true, the analytical
%  benchmark is deliberately disabled because the external eligibility mask
%  changes the TEM estimating equations; bootstrap calibration is used. The
%  bootstrap p-values in out.pvalue remain the primary finite-sample
%  calibration.
%
% See also: mdEM, mdImputeCondMean.m, mdPartialMD.m, mdPartialMD2full
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
% Copyright 2008-2026.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mdMDPtest')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example 1: Basic call with default options.
    % Load data with missing values and run the test with default settings.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X);
    % Display observed statistics and p-values
    disp(out.Tobs)
    disp(out.pvalue)
%}

%{
    %% Example 2: Test with trimming.
    % Run the bootstrap test using TEM with trimming level alpha=0.25.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'nsimul',199);
    % Display p-values
    disp(out.pvalue)
%}

%{
    %% Example 3: Simulated data under MCAR.
    % Generate Gaussian data with MCAR missingness and apply the test.
    rng(1)
    n = 300;
    p = 5;
    rho = 0.5;
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    mu = zeros(1,p);
    
    Yfull = mvnrnd(mu,Sigma,n);
    
    missRate = 0.10;
    missMask = rand(n,p) < missRate;
    Y = Yfull;
    Y(missMask) = NaN;
    % Show also the output plot
    out = mdMDPtest(Y,'nsimul',199,'plots',true);
    
    disp('Observed statistics:')
    disp(out.Tobs)
    disp('Bootstrap p-values:')
    disp(out.pvalue)
%}


%{
    % Example 4: Comparison of several trimming levels.
    rng(1)
    n = 300;
    p = 5;
    rho = 0.5;
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    mu = zeros(1,p);
    Yfull = mvnrnd(mu,Sigma,n);
    missRate = 0.10;
    missMask = rand(n,p) < missRate;
    Y = Yfull;
    Y(missMask) = NaN;    
    Alpha = [0 0.10 0.25 0.50]';
    pval = zeros(length(Alpha),4);
    for i=1:length(Alpha)
        out = mdMDPtest(Y,'alpha',Alpha(i),'nsimul',199);
        pval(i,:) = out.pvalue;
    end
    pvalTable = array2table(pval, ...
        'VariableNames',{'medLogRatio','meanLogRatio','medDiff','meanDiff'}, ...
        'RowNames',string(Alpha));
    disp(pvalTable)
%}

%{
    % Example 5: Different rescaling method.
    % Use method betaMap instead of the default pri.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'method','betaMap','nsimul',199);
    disp(out.pvalue)
%}

%{
    %% Example 6: MM estimator for the complete cases.
    % Use an MM estimator with 95 percent nominal efficiency. The
    % preliminary S estimator has breakdown point alpha=0.25. Outliers are
    % declared using the fixed Bonferronized simultaneous level 0.99.
    load cows2026
    X = cows2026{:,:};
    r = struct;
    r.class = 'MM';
    r.eff = 0.95;
    out = mdMDPtest(X,'alpha',0.25,'robust',r,'nsimul',199);
    disp(out.pvalue)
%}

%{
    %% Example 7: Standard FSM envelope-resuperimposition rule.
    % Use the standard FSM signal-detection and validation procedure rather
    % than the direct Bonferroni-bound stopping rule.
    load cows2026
    X = cows2026{:,:};
    r = struct;
    r.class = 'FS';
    r.bonflev = [];
    out = mdMDPtest(X,'alpha',0.25,'robust',r,'nsimul',199);
    disp(out.FSrule)
    disp(out.pvalue)
%}

%{
    %% Example 8: Estimator-defined inlier mean (MDP-I).
    % The medians use all complete cases. The two means use only complete
    % cases not declared as outliers by the MCD complete-case fit.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'robust','MCD', ...
        'aggregation','inlier','nsimul',199);
    disp(out.nMeanCC)
    disp(out.pvalue)
%}

%{
    %% Example 9: Fixed robust-distance filtered mean (MDP-F).
    % Use the chi-square 0.99 cutoff on the complete-case robust distances.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'robust','MCD', ...
        'aggregation','fixed','filterlev',0.99,'nsimul',199);
    disp(out.filterCutoff)
    disp(out.nMeanCC)
    disp(out.pvalue)
    disp([out.asympt.kc out.asympt.lambda])
    disp([out.asympt.TDmean.sigma2 out.asympt.TLmean.sigma2])
%}

%{
    %% Example 10: Coupled complete-case outlier trimming in TEM.
    % This sensitivity mode uses the MDP-I aggregation rule and additionally
    % forces estimator-declared complete-case outliers to zero TEM weight.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'robust','MCD', ...
        'aggregation','inlier','coupledtrim',true,'nsimul',199);
    disp(out.nCoupledCC)
    disp(out.coupledRows)
    disp(out.pvalue)
    disp(out.asympt.reason)
%}

%{
    %% Example 11: Distribution-adaptive TEM correction.
    % Use the MCD complete-case fit as the reference geometry for empirical
    % pattern-wise radial consistency factors. The same pipeline is rerun
    % independently inside every bootstrap sample.
    load cows2026
    X = cows2026{:,:};
    out = mdMDPtest(X,'alpha',0.25,'robust','MCD', ...
        'consistencyfactor','adaptive','nsimul',199);
    disp(out.TEMkinfo)
    disp(out.pvalue)
    disp(out.asympt.reason)
%}

%% Beginning of code

if ~ismatrix(Y) || ~isnumeric(Y)
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Input argument Y must be a numeric matrix.');
end

% Default options
options = struct;
options.alpha   = 0;
options.method  = 'pri';
options.nsimul  = 499;
options.conflev = 0.95;
options.tol     = 1e-10;
options.plots   = false;
options.robust  ="FS";
options.eps0    = 1e-12;
options.aggregation = 'ordinary';
options.filterlev = 0.99;
options.coupledtrim = false;
options.consistencyfactor = 'pattern';
options.adaptivepool = false;
options.adaptiveminref = 20;

% Check supplied options
if ~isempty(varargin)
    if mod(length(varargin),2) ~= 0
        error('FSDA:mdMDPtest:WrongInputOpt', ...
            'Optional arguments must be supplied in name/value pairs.');
    end

    UserOptions = varargin(1:2:end);
    if ~isempty(UserOptions)
        aux.chkoptions(options,UserOptions);
    end

    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
end

alpha   = options.alpha;
method  = char(options.method);
nsimul  = options.nsimul;
conflev = options.conflev;
tol     = options.tol;
plots   = options.plots;
robust  = options.robust;
eps0 = options.eps0;
aggregation = local_parse_aggregation(options.aggregation);
filterlev = options.filterlev;
coupledtrim = options.coupledtrim;
consistencyfactor = local_parse_consistencyfactor(options.consistencyfactor);
adaptivepool = options.adaptivepool;
adaptiveminref = options.adaptiveminref;

if ~isscalar(alpha) || ~isnumeric(alpha) || alpha < 0 || alpha > 0.5
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option alpha must be a scalar in the interval [0,0.5].');
end

if ~isscalar(nsimul) || nsimul <= 0 || nsimul ~= floor(nsimul)
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option nsimul must be a positive integer.');
end

if ~isscalar(conflev) || conflev <= 0 || conflev >= 1
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option conflev must be a scalar in the interval (0,1).');
end

if ~isscalar(tol) || tol <= 0
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option tol must be a positive scalar.');
end

if ~isscalar(eps0) || ~isnumeric(eps0) || ~isfinite(eps0) || eps0 <= 0
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option eps0 must be a finite positive scalar.');
end

if ~isscalar(filterlev) || ~isnumeric(filterlev) || ...
        ~isfinite(filterlev) || filterlev <= 0 || filterlev >= 1
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option filterlev must be a finite scalar in the interval (0,1).');
end

if ~(islogical(coupledtrim) && isscalar(coupledtrim))
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option coupledtrim must be a logical scalar.');
end

if ~(islogical(adaptivepool) && isscalar(adaptivepool))
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option adaptivepool must be a logical scalar.');
end

if ~isscalar(adaptiveminref) || ~isnumeric(adaptiveminref) || ...
        ~isfinite(adaptiveminref) || adaptiveminref < 1 || ...
        adaptiveminref ~= floor(adaptiveminref)
    error('FSDA:mdMDPtest:WrongInputOpt', ...
        'Option adaptiveminref must be a positive integer.');
end

if coupledtrim && ~strcmp(consistencyfactor,'pattern')
    error('FSDA:mdMDPtest:CoupledConsistencyFactor', ...
        ['The current coupledtrim sensitivity implementation supports only ' ...
         'consistencyfactor=''pattern''.']);
end

if coupledtrim && alpha == 0
    error('FSDA:mdMDPtest:CoupledNeedsRobust', ...
        'Option coupledtrim=true requires alpha>0.');
end

if coupledtrim && ~strcmp(aggregation,'inlier')
    error('FSDA:mdMDPtest:CoupledNeedsInlier', ...
        ['Option coupledtrim=true is currently defined only for ' ...
         'aggregation=''inlier''.']);
end

if strcmp(aggregation,'inlier') && alpha == 0
    error('FSDA:mdMDPtest:InlierNeedsRobust', ...
        ['Aggregation ''inlier'' requires alpha>0 because the complete-case ' ...
         'outlier set must be supplied by a robust complete-case fit.']);
end

% Parse the complete-case robust estimator. A simple character/string
% specification is retained for convenience; a structure permits
% estimator-specific tuning such as the MM efficiency.
[robustClass,robustEff,robustBonflev] = local_parse_robust(robust);

[n,p] = size(Y);
maskMiss = isnan(Y);
completeIdx = all(~maskMiss,2);
nComplete = sum(completeIdx);

if nComplete < p + 2
    error('FSDA:mdMDPtest:TooFewCompleteRows', ...
        ['Too few complete rows to compute the reference complete-case ' ...
         'covariance matrix.']);
end

if alpha > 0 && strcmpi(robustClass,'FS')

    h = floor(nComplete*(1-alpha));

    if h < p+1
        error('FSDA:mdMDPtest:TooSmallFSInit', ...
            ['The value of alpha produces an FSM initial monitoring step ' ...
             'smaller than p+1. Decrease alpha or increase the number ' ...
             'of complete observations.']);
    end

    nMon = nComplete-h;

    if nMon < 3
        error('FSDA:mdMDPtest:TooLargeFSInit', ...
            ['The value of alpha leaves fewer than three Forward Search ' ...
             'monitoring steps. Increase alpha.']);

    elseif nMon < 4 && isempty(robustBonflev)

        robustBonflev = 0.99;

        warning('FSDA:mdMDPtest:FSBonferroniFallback', ...
            ['Fewer than four Forward Search monitoring steps are ' ...
             'available. The direct 99 percent Bonferroni stopping ' ...
             'rule is used instead of the standard FSM ' ...
             'envelope-resuperimposition rule.']);
    end
end


% Extraction of complete cases
Ycc = Y(completeIdx,:);

% Robust mu and Sigma based on complete cases
[muCC,SigCC,outRob] = local_complete_case_fit(Ycc,alpha, ...
    robustClass,robustEff,robustBonflev);

% Distances based on complete rows only
d2_cc = mahalFS(Ycc, muCC, SigCC);

% Construct the complete-case eligibility mask used by the optional
% coupled TEM sensitivity mode. With coupledtrim=false no row is forced out.
forcedZeroTEM = false(n,1);
if coupledtrim
    ccOutlierMask = local_cc_outlier_mask(outRob,nComplete);
    completeRows = find(completeIdx);
    forcedZeroTEM(completeRows(ccOutlierMask)) = true;
end

% Distances based on EM/TEM fit using all rows. The fitted object is kept
% because the alpha>0 analytical benchmark needs the final TEM weights and
% pattern-wise trimming information. In coupled mode the fixed eligibility
% mask is multiplied by TEM's own concentration weights at every iteration.
[d2_all_cc, muHat, SigHat, outFit] = local_fit_and_get_complete_distances( ...
    Y, completeIdx, alpha, method, tol, forcedZeroTEM, ...
    consistencyfactor, Ycc, muCC, SigCC, adaptivepool, adaptiveminref);

% Construct the observed-data aggregation rule for the two mean statistics.
% The median statistics always use all complete cases.
[meanWeightsCC,aggInfo] = local_aggregation_weights(d2_cc,outRob,alpha, ...
    aggregation,filterlev,p);

if aggInfo.nSelected == 0
    error('FSDA:mdMDPtest:NoAggregationRows', ...
        'The selected mean aggregation rule retains no complete cases.');
end

% Observed statistics
Tobs = local_statistic(d2_cc,d2_all_cc,eps0,meanWeightsCC);

% Analytical Gaussian benchmark for the selected mean statistics. First
% compute the base estimator-discrepancy variance and then apply the scalar
% first-order coefficient implied by the selected aggregation rule.
if coupledtrim
    asympt = local_unavailable_asymptotic(nComplete/n);
    asympt.reason = ['Analytical calibration is not supplied when ' ...
        'coupledtrim=true because the external complete-case eligibility ' ...
        'constraint changes the TEM estimating equations. Use the ' ...
        'mask-preserving bootstrap for this sensitivity mode.'];
    asympt.mode = 'coupled TEM sensitivity';
    asympt.theoryStatus = ['Bootstrap calibration only for coupled TEM. ' ...
        'The uncoupled analytical reference law is not applied.'];
    asympt.aggregation = aggregation;
else
    if alpha == 0
        asympt = local_classical_asymptotic(maskMiss, SigHat, nComplete, ...
            Tobs, eps0);
    elseif strcmp(consistencyfactor,'pattern')
        asympt = local_tem_asymptotic(Y, maskMiss, completeIdx, outFit, ...
            alpha, method, robustClass, robustEff, robustBonflev, ...
            Tobs, eps0);
    else
        asympt = local_unavailable_asymptotic(nComplete/n);
        if strcmp(consistencyfactor,'adaptive')
            asympt.reason = ['The adaptive TEM fit is available, but its ' ...
                'analytical sandwich calibration is not yet implemented in ' ...
                'mdMDPtest. The current adaptive result is calibrated by ' ...
                'the mask-preserving bootstrap.'];
            asympt.mode = 'adaptive TEM bootstrap calibration';
            asympt.theoryStatus = ['Adaptive radial consistency factors are ' ...
                'used in the TEM fit. Their first-order influence contribution ' ...
                'is deliberately not replaced by the Gaussian pattern formula.'];
        else
            asympt.reason = ['For alpha>0 the analytical TEM benchmark is ' ...
                'currently implemented only for consistencyfactor=''pattern''.'];
            asympt.mode = 'TEM bootstrap calibration';
            asympt.theoryStatus = 'Bootstrap calibration for the selected TEM correction.';
        end
    end

    asympt = local_apply_aggregation_asymptotic(asympt,aggregation,aggInfo, ...
        p,eps0,n,Tobs);
end

% Bootstrap under MCAR
Tboot = NaN(nsimul,4);
nMeanBoot = NaN(nsimul,1);
nCoupledBoot = NaN(nsimul,1);

% Use robust complete-case fit to generate bootstrap samples
SigGen = local_make_spd(SigCC);
R = chol(SigGen,'upper');

for j = 1:nsimul

    % Generate full data from Gaussian model
    YfullStar = randn(n,p) * R + muCC(:)';

    % Impose the observed missingness pattern
    Ystar = YfullStar;
    Ystar(maskMiss) = NaN;

    % Complete-case reference distances in bootstrap world
    YccStar = YfullStar(completeIdx,:);

    [muCCStar,SigCCStar,outRobStar] = local_complete_case_fit(YccStar,alpha, ...
        robustClass,robustEff,robustBonflev);

    d2_cc_star = mahalFS(YccStar, muCCStar, SigCCStar);

    % Reconstruct the coupled complete-case eligibility mask independently
    % in every bootstrap sample.
    forcedZeroStar = false(n,1);
    if coupledtrim
        ccOutlierMaskStar = local_cc_outlier_mask(outRobStar,nComplete);
        completeRows = find(completeIdx);
        forcedZeroStar(completeRows(ccOutlierMaskStar)) = true;
    end
    nCoupledBoot(j) = sum(forcedZeroStar);

    % EM/TEM distances for the same complete rows.
    d2_all_cc_star = local_fit_and_get_complete_distances( ...
        Ystar, completeIdx, alpha, method, tol, forcedZeroStar, ...
        consistencyfactor, YccStar, muCCStar, SigCCStar, ...
        adaptivepool, adaptiveminref);

    % Reconstruct the selected aggregation rule inside this bootstrap sample.
    [meanWeightsStar,aggInfoStar] = local_aggregation_weights( ...
        d2_cc_star,outRobStar,alpha,aggregation,filterlev,p);
    nMeanBoot(j) = aggInfoStar.nSelected;

    % Store the four statistics. If the selected mean rule retains no rows,
    % the two mean statistics are NaN and this bootstrap sample is discarded.
    Tboot(j,:) = local_statistic(d2_cc_star,d2_all_cc_star,eps0, ...
        meanWeightsStar);
end

% Remove bootstrap samples containing NaNs
validBoot = all(~isnan(Tboot),2);
Tboot = Tboot(validBoot,:);
nMeanBoot = nMeanBoot(validBoot);
nCoupledBoot = nCoupledBoot(validBoot);

if isempty(Tboot)
    error('FSDA:mdMDPtest:NoValidBootstrap', ...
        'All bootstrap replicates failed.');
end

% Bootstrap p-values
pvalue = (1 + sum(abs(Tboot) >= abs(Tobs),1)) / (size(Tboot,1) + 1);

% Bootstrap confidence intervals
alphaCI = 1 - conflev;
ciBoot = quantile(Tboot,[alphaCI/2 1-alphaCI/2],1);

% Store output
out = struct;
out.pvalue      = pvalue;
out.Tobs        = Tobs;
out.Tboot       = Tboot;
out.alpha       = alpha;
out.method      = method;
out.consistencyfactor = consistencyfactor;
out.adaptivepool = adaptivepool;
out.adaptiveminref = adaptiveminref;
if isfield(outFit,'kfactor')
    out.TEMkfactor = outFit.kfactor;
else
    out.TEMkfactor = NaN;
end
if isfield(outFit,'kinfo')
    out.TEMkinfo = outFit.kinfo;
else
    out.TEMkinfo = [];
end
out.aggregation = aggregation;
out.filterlev = filterlev;
out.filterCutoff = aggInfo.cutoff;
out.coupledtrim = coupledtrim;
out.coupledRows = find(forcedZeroTEM);
out.nCoupledCC = sum(forcedZeroTEM);
out.nCoupledBoot = nCoupledBoot;
if isfield(outFit,'weights')
    out.TEMweights = outFit.weights(:);
    if isfield(outFit,'internalWeights')
        out.TEMinternalWeights = outFit.internalWeights(:);
    else
        out.TEMinternalWeights = outFit.weights(:);
    end
else
    % Untrimmed EM uses every row; expose this as unit TEM weights so that
    % the output fields have a consistent length for alpha=0 as well.
    out.TEMweights = true(n,1);
    out.TEMinternalWeights = true(n,1);
end
out.meanWeightsCC = meanWeightsCC;
completeRows = find(completeIdx);
out.meanRows = completeRows(meanWeightsCC);
out.nMeanCC = aggInfo.nSelected;
out.nMeanBoot = nMeanBoot;
out.nComplete   = nComplete;
out.completeIdx = completeIdx;

out.locCC = muCC;
out.covCC = SigCC;
out.robust = robust;
out.robustClass = robustClass;
out.robustEff = robustEff;
out.robustBonflev = robustBonflev;
out.outlierConflev = 0.99;

if alpha > 0 && strcmpi(robustClass,'FS')
    if isempty(robustBonflev)
        out.FSrule = 'resuperimposition';
    else
        out.FSrule = 'bonferroni';
    end
else
    out.FSrule = '';
end

if alpha > 0
    out.outliersCC = outRob.outliers;
    out.nOutliersCC = sum(local_cc_outlier_mask(outRob,nComplete));
else
    out.outliersCC = [];
    out.nOutliersCC = 0;
end

out.d2_cc       = d2_cc;
out.d2_all      = d2_all_cc;

out.ciBoot      = ciBoot;
out.loc       = muHat;
out.cov      = SigHat;
out.eps0 = eps0;
out.asympt = asympt;

% Optional plots
if plots
    statNames = {'median log-ratio', ...
                 'mean log-ratio', ...
                 'median difference', ...
                 'mean difference'};

    figure;
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    % Left part: scatter plot of the two distances
    nexttile([2 1])
    scatter(d2_cc,d2_all_cc,'o')
    refline(1,0)
    xlabel('Complete-case distance')
    ylabel('Distance from all-data EM/TEM')
    title(['\alpha=' num2str(alpha) ', aggregation=' aggregation])
    box on

    % Right part: bootstrap distributions of the four statistics
    for j = 1:4
        nexttile
        histogram(Tboot(:,j))
        hold on
        xline(Tobs(j),'r','LineWidth',1.5)
        title(['p=' num2str(pvalue(j),4)])
        xlabel(statNames{j})
        box on
    end
end

end

% -------------------------------------------------------------------------

function [muCC,SigCC,outRob] = ...
    local_complete_case_fit(Ycc,alpha,robustClass,robustEff,robustBonflev)
%local_complete_case_fit estimates location and scatter from complete cases.
%
% For all robust estimators, outliers are declared at simultaneous
% Bonferronized confidence level 0.99. For FSM, robustBonflev=0.99 uses
% the direct Bonferroni bound, whereas robustBonflev=[] uses the standard
% signal-detection, validation and envelope-resuperimposition rules. MCD
% and MM use the corresponding individual level 1-0.01/ncc.

ncc = size(Ycc,1);
conflevOut = 1-0.01/ncc;

if alpha == 0

    muCC = mean(Ycc,1);
    SigCC = cov(Ycc);
    outRob = [];

elseif strcmpi(robustClass,'FS')

    h = floor(ncc*(1-alpha));

    if isempty(robustBonflev)
        % Standard FSM signal detection and envelope resuperimposition.
        outRob = FSM(Ycc, ...
            'init',h, ...
            'plots',0, ...
            'msg',false);
    else
        % Direct Bonferroni-bound stopping rule.
        outRob = FSM(Ycc, ...
            'init',h, ...
            'bonflev',robustBonflev, ...
            'plots',0, ...
            'msg',false);
    end

    muCC = outRob.loc;
    SigCC = outRob.cov;

elseif strcmpi(robustClass,'MCD')

    outRob = mcd(Ycc, ...
        'bdp',alpha, ...
        'conflev',conflevOut, ...
        'plots',0, ...
        'msg',0);

    muCC = outRob.loc;
    SigCC = outRob.cov;

elseif strcmpi(robustClass,'MM')

    outRob = MMmult(Ycc, ...
        'Sbdp',alpha, ...
        'eff',robustEff, ...
        'conflev',conflevOut, ...
        'plots',0, ...
        'Smsg',0);

    muCC = outRob.loc;
    SigCC = outRob.cov;

else
    error('FSDA:mdMDPtest:WrongRobust', ...
        'Robust estimator must be ''FS'', ''MCD'' or ''MM''.');
end

end

% -------------------------------------------------------------------------
function [robustClass,robustEff,robustBonflev] = local_parse_robust(robust)
%local_parse_robust parses the robust option.

robustEff = 0.95;
robustBonflev = []; % default is envelope resuperimposition

if isstruct(robust)
    if ~isscalar(robust)
        error('FSDA:mdMDPtest:WrongRobust', ...
            'Option robust must be a scalar structure.');
    end

    allowedFields = {'class','eff','bonflev'};
    suppliedFields = fieldnames(robust);
    wrongFields = setdiff(suppliedFields,allowedFields);
    if ~isempty(wrongFields)
        error('FSDA:mdMDPtest:WrongRobustField', ...
            'Unknown field in option robust: %s.',wrongFields{1});
    end

    if ~isfield(robust,'class') || isempty(robust.class)
        error('FSDA:mdMDPtest:MissingRobustClass', ...
            ['When robust is a structure, field robust.class must be ' ...
            'specified as ''FS'', ''MCD'' or ''MM''.']);
    end
    robustClass = robust.class;

    if isfield(robust,'eff') && ~isempty(robust.eff)
        robustEff = robust.eff;
    end

    if isfield(robust,'bonflev')
        robustBonflev = robust.bonflev;
    end
else
    robustClass = robust;
end

if isstring(robustClass)
    if ~isscalar(robustClass)
        error('FSDA:mdMDPtest:WrongRobust', ...
            'Option robust.class must be a character vector or string scalar.');
    end
    robustClass = char(robustClass);
end

if ~ischar(robustClass) || size(robustClass,1) ~= 1
    error('FSDA:mdMDPtest:WrongRobust', ...
        'Option robust must specify ''FS'', ''MCD'' or ''MM''.');
end

robustClass = upper(strtrim(robustClass));
if ~any(strcmp(robustClass,{'FS','MCD','MM'}))
    error('FSDA:mdMDPtest:WrongRobust', ...
        'Option robust must specify ''FS'', ''MCD'' or ''MM''.');
end

if ~isscalar(robustEff) || ~isnumeric(robustEff) || ...
        ~isfinite(robustEff) || robustEff < 0.5 || robustEff > 0.99
    error('FSDA:mdMDPtest:WrongRobustEff', ...
        'Field robust.eff must be a scalar in the interval [0.5,0.99].');
end

if strcmp(robustClass,'FS')
    if ~(isempty(robustBonflev) || ...
            (isnumeric(robustBonflev) && isscalar(robustBonflev) && ...
             isfinite(robustBonflev) && robustBonflev == 0.99))
        error('FSDA:mdMDPtest:WrongRobustBonflev', ...
            ['Field robust.bonflev must be either 0.99 or empty. ' ...
             'Use 0.99 for the direct Bonferroni rule and [] for the ' ...
             'standard FSM envelope-resuperimposition rule.']);
    end
else
    if isstruct(robust) && isfield(robust,'bonflev')
        error('FSDA:mdMDPtest:WrongRobustBonflev', ...
            'Field robust.bonflev can be supplied only when robust.class=''FS''.');
    end
    robustBonflev = [];
end

end


function [d2_all_cc, muHat, SigHat, outFit] = local_fit_and_get_complete_distances( ...
    Y, completeIdx, alpha, method, tol, forcedZero, consistencyfactor, ...
    Yref, muref, sigmaref, adaptivepool, adaptiveminref)
%local_fit_and_get_complete_distances fits EM/TEM and returns complete-row distances.
%
% forcedZero is a logical n-vector. When alpha>0 and at least one entry is
% true, TEM's own concentration indicator is multiplied by ~forcedZero at
% every iteration. This implements the coupled sensitivity rule without
% changing the default mdTEM fit. For the adaptive correction, Yref, muref
% and sigmaref are the complete-case reference sample and fitted geometry.

p = size(Y,2);
n = size(Y,1);
if isempty(forcedZero)
    forcedZero = false(n,1);
else
    forcedZero = logical(forcedZero(:));
    if numel(forcedZero) ~= n
        error('FSDA:mdMDPtest:WrongCoupledMask', ...
            'The coupled TEM mask must have one entry for every row of Y.');
    end
end

if alpha == 0
    if any(forcedZero)
        error('FSDA:mdMDPtest:CoupledNeedsRobust', ...
            'Coupled TEM is available only when alpha>0.');
    end
    outFit = mdEM(Y);
elseif any(forcedZero)
    outFit = local_coupled_tem(Y,forcedZero,alpha,method,tol);
else
    temArgs = {'method',method,'alpha',alpha,'tol',tol, ...
        'consistencyfactor',consistencyfactor};
    if strcmp(consistencyfactor,'adaptive')
        temArgs = [temArgs, {'Yref',Yref,'muref',muref(:), ...
            'sigmaref',sigmaref,'adaptivepool',adaptivepool, ...
            'adaptiveminref',adaptiveminref}];
    end
    outFit = mdTEM(Y,temArgs{:});
    outFit.internalWeights = outFit.weights;
    outFit.forcedZero = forcedZero;
end

muHat = outFit.loc;
SigHat = outFit.cov;

[d2_part, poss] = mdPartialMD(Y, muHat, SigHat);
d2_full = mdPartialMD2full(d2_part, p, poss, 'method', method);
d2_all_cc = d2_full(completeIdx);

end

% -------------------------------------------------------------------------
function consistencyfactor = local_parse_consistencyfactor(consistencyfactor)
%local_parse_consistencyfactor parses the mdTEM consistency-factor option.

if isstring(consistencyfactor)
    if ~isscalar(consistencyfactor)
        error('FSDA:mdMDPtest:WrongConsistencyFactor', ...
            ['Option consistencyfactor must be a character vector or ' ...
             'string scalar.']);
    end
    consistencyfactor = char(consistencyfactor);
end

if ~ischar(consistencyfactor) || size(consistencyfactor,1) ~= 1
    error('FSDA:mdMDPtest:WrongConsistencyFactor', ...
        ['Option consistencyfactor must be ''pattern'', ''adaptive'', ' ...
         '''global'', ''weighted'' or ''none''.']);
end

consistencyfactor = lower(strtrim(consistencyfactor));
valid = {'pattern','adaptive','global','weighted','none'};
if ~any(strcmp(consistencyfactor,valid))
    error('FSDA:mdMDPtest:WrongConsistencyFactor', ...
        ['Option consistencyfactor must be ''pattern'', ''adaptive'', ' ...
         '''global'', ''weighted'' or ''none''.']);
end
end

% -------------------------------------------------------------------------
function aggregation = local_parse_aggregation(aggregation)
%local_parse_aggregation parses the mean aggregation option.

if isstring(aggregation)
    if ~isscalar(aggregation)
        error('FSDA:mdMDPtest:WrongAggregation', ...
            'Option aggregation must be a character vector or string scalar.');
    end
    aggregation = char(aggregation);
end


if ~ischar(aggregation) || size(aggregation,1) ~= 1
    error('FSDA:mdMDPtest:WrongAggregation', ...
        'Option aggregation must be ''ordinary'', ''inlier'' or ''fixed''.');
end


aggregation = lower(strtrim(aggregation));
if ~any(strcmp(aggregation,{'ordinary','inlier','fixed'}))
    error('FSDA:mdMDPtest:WrongAggregation', ...
        'Option aggregation must be ''ordinary'', ''inlier'' or ''fixed''.');
end
end

% -------------------------------------------------------------------------
function [A,info] = local_aggregation_weights(d2_cc,outRob,alpha, ...
    aggregation,filterlev,p)
%local_aggregation_weights builds the complete-case mean aggregation rule.
%
% A is a logical vector indexed within the complete-case sample. The median
% statistics do not use A. For MDP-I, outlier labels are reconstructed from
% the complete-case robust fit. For MDP-F, the cutoff is fixed at the
% chi-square filterlev quantile and is applied to the complete-case distances.

ncc = numel(d2_cc);
A = true(ncc,1);
cutoff = NaN;

switch aggregation
    case 'ordinary'
        % MDP-U: all complete cases enter the two mean statistics.

    case 'inlier'
        % MDP-I: remove complete cases declared as outliers by the selected
        % robust complete-case estimator.
        if alpha == 0 || isempty(outRob) || ~isfield(outRob,'outliers')
            error('FSDA:mdMDPtest:MissingOutlierSet', ...
                ['Aggregation ''inlier'' requires a robust complete-case ' ...
                 'fit returning the field outliers.']);
        end

        A = ~local_cc_outlier_mask(outRob,ncc);

    case 'fixed'
        % MDP-F: retain only complete cases inside the fixed robust-distance
        % region. The same fixed cutoff is used in every bootstrap sample.
        cutoff = chi2inv(filterlev,p);
        A = isfinite(d2_cc(:)) & d2_cc(:) <= cutoff;
end

info = struct;
info.nSelected = sum(A);
info.fractionSelected = info.nSelected/ncc;
info.cutoff = cutoff;
end

% -------------------------------------------------------------------------
function outlierMask = local_cc_outlier_mask(outRob,ncc)
%local_cc_outlier_mask converts robust complete-case outliers to a logical mask.

outlierMask = false(ncc,1);
if isempty(outRob) || ~isstruct(outRob) || ~isfield(outRob,'outliers')
    return
end

outliers = outRob.outliers;
if islogical(outliers) && numel(outliers) == ncc
    outlierMask = outliers(:);
elseif ~isempty(outliers)
    outliers = outliers(:);
    if isnumeric(outliers)
        outliers = outliers(isfinite(outliers) & outliers == floor(outliers) & ...
            outliers >= 1 & outliers <= ncc);
        outlierMask(unique(outliers)) = true;
    end
end
end

% -------------------------------------------------------------------------
function T = local_statistic(d2_cc, d2_all, eps0, meanWeights)
%local_statistic computes the four MDP test statistics.
%
% The median statistics use all complete cases. The two mean statistics use
% the logical selection vector meanWeights.

d2_cc = d2_cc(:);
d2_all = d2_all(:);
meanWeights = logical(meanWeights(:));

if numel(d2_cc) ~= numel(d2_all) || numel(meanWeights) ~= numel(d2_cc)
    error('FSDA:mdMDPtest:WrongAggregationWeights', ...
        'Distance vectors and aggregation weights must have the same length.');
end

rat = log((d2_all + eps0) ./ (d2_cc + eps0));
dif = d2_all - d2_cc;

T1 = median(rat);
T3 = median(dif);

if any(meanWeights)
    T2 = mean(rat(meanWeights));
    T4 = mean(dif(meanWeights));
else
    T2 = NaN;
    T4 = NaN;
end

T = [T1 T2 T3 T4];
end

% -------------------------------------------------------------------------
function asympt = local_apply_aggregation_asymptotic(asympt,aggregation, ...
    aggInfo,p,eps0,n,Tobs)
%local_apply_aggregation_asymptotic applies MDP-U/I/F first-order scaling.
%
% The estimator-discrepancy calculation in local_classical_asymptotic or
% local_tem_asymptotic supplies the base variance sigma_D^2 associated with
%   -tr{Sigma^{-1}(SigmaHat_A-SigmaHat_C)}.
% The selected mean aggregation changes only its scalar first-order
% coefficient:
%   MDP-U: TD coefficient 1,       TL coefficient kappa;
%   MDP-I: TD coefficient 1,       TL coefficient kappa;
%   MDP-F: TD coefficient k_c,     TL coefficient lambda.
%
% For MDP-I the equality with MDP-U is conditional on the clean
% complete-case rule excluding only O_p(1) observations and on the maximal
% omitted contrasts being o_p(sqrt(n)). For MDP-F, c is fixed as n grows.

asympt.aggregation = aggregation;
asympt.filterCutoff = NaN;
asympt.kc = NaN;
asympt.lambda = NaN;
asympt.gammaFilter = NaN;

% If the base estimator-discrepancy benchmark is unavailable, retain its
% reason and only document the requested aggregation.
if ~isfield(asympt,'available') || ~asympt.available
    if strcmp(aggregation,'fixed')
        asympt.filterCutoff = aggInfo.cutoff;
    end
    return
end

baseSigmaD2 = asympt.TDmean.sigma2;
baseDegenerate = asympt.degenerate;
if ~isfinite(baseSigmaD2) || baseSigmaD2 < 0
    asympt.available = false;
    asympt.reason = 'The base analytical variance is not finite.';
    return
end

% Store the unscaled estimator-discrepancy variance for diagnostics.
asympt.baseSigmaD2 = baseSigmaD2;
if isfield(asympt,'TLmean') && isfield(asympt.TLmean,'kappa')
    baseKappa = asympt.TLmean.kappa;
else
    baseKappa = NaN;
end
asympt.baseKappa = baseKappa;

switch aggregation
    case 'ordinary'
        coefD = 1;
        coefL = baseKappa;
        asympt.kc = 1;
        asympt.lambda = baseKappa;
        asympt.gammaFilter = 1;

    case 'inlier'
        coefD = 1;
        coefL = baseKappa;
        asympt.kc = 1;
        asympt.lambda = baseKappa;
        asympt.gammaFilter = 1;
        note = ['Aggregation=''inlier'' uses the MDP-U first-order law ' ...
            'under the condition that the number of falsely excluded clean ' ...
            'complete cases is O_p(1) and the maximal omitted contrasts are ' ...
            'o_p(sqrt(n)).'];
        if isfield(asympt,'theoryStatus') && ~isempty(asympt.theoryStatus)
            asympt.theoryStatus = [asympt.theoryStatus ' ' note];
        else
            asympt.theoryStatus = note;
        end

    case 'fixed'
        c = aggInfo.cutoff;
        gamma = chi2cdf(c,p);
        if ~isfinite(c) || c <= 0 || ~isfinite(gamma) || gamma <= 0
            asympt.available = false;
            asympt.reason = ['Unable to evaluate the fixed-gate analytical ' ...
                'coefficients.'];
            return
        end

        kc = chi2cdf(c,p+2)/gamma;
        lambda = integral(@(q) local_kappa_integrand(q,p,eps0),0,c, ...
            'RelTol',1e-10,'AbsTol',1e-12)/(p*gamma);

        if ~isfinite(kc) || kc <= 0 || ~isfinite(lambda) || lambda <= 0
            asympt.available = false;
            asympt.reason = ['The fixed-gate analytical coefficients are ' ...
                'not finite and positive.'];
            return
        end

        coefD = kc;
        coefL = lambda;
        asympt.filterCutoff = c;
        asympt.kc = kc;
        asympt.lambda = lambda;
        asympt.gammaFilter = gamma;

        note = ['Aggregation=''fixed'' applies the fixed-gate first-order ' ...
            'Tallis scaling to the base estimator-discrepancy variance.'];
        if isfield(asympt,'theoryStatus') && ~isempty(asympt.theoryStatus)
            asympt.theoryStatus = [asympt.theoryStatus ' ' note];
        else
            asympt.theoryStatus = note;
        end

    otherwise
        error('FSDA:mdMDPtest:WrongAggregation', ...
            'Unknown aggregation rule in analytical scaling.');
end

if ~isfinite(coefL) || coefL <= 0
    asympt.available = false;
    asympt.reason = 'The analytical log-ratio coefficient is not finite.';
    return
end

% Apply the selected scalar coefficient to the common estimator-discrepancy
% influence law.
asympt.TDmean.coefficient = coefD;
asympt.TDmean.sigma2 = coefD^2*baseSigmaD2;
asympt.TDmean.se = sqrt(asympt.TDmean.sigma2/n);

asympt.TLmean.coefficient = coefL;
asympt.TLmean.sigma2 = coefL^2*baseSigmaD2;
asympt.TLmean.se = sqrt(asympt.TLmean.sigma2/n);

% kappa remains the unfiltered stabilized-log coefficient. For the fixed
% gate, lambda is reported separately and is the coefficient actually used.
if ~isfield(asympt.TLmean,'kappa')
    asympt.TLmean.kappa = baseKappa;
end
asympt.TLmean.lambda = coefL;

% TEM variance components, when available, inherit the same TD multiplier.
if isfield(asympt.TDmean,'components') && isstruct(asympt.TDmean.components)
    fn = fieldnames(asympt.TDmean.components);
    for j = 1:numel(fn)
        val = asympt.TDmean.components.(fn{j});
        if isnumeric(val) && isscalar(val) && isfinite(val)
            asympt.TDmean.components.(fn{j}) = coefD^2*val;
        end
    end
end

% Preserve the base TEM and complete-case diagnostics for backward
% compatibility, and add their contribution on the selected TD scale.
if isfield(asympt,'TEM') && isstruct(asympt.TEM) && ...
        isfield(asympt.TEM,'sigma2') && isfinite(asympt.TEM.sigma2)
    asympt.TEM.sigma2Selected = coefD^2*asympt.TEM.sigma2;
end
if isfield(asympt,'completeCase') && isstruct(asympt.completeCase)
    if isfield(asympt.completeCase,'sigma2') && ...
            isfinite(asympt.completeCase.sigma2)
        asympt.completeCase.sigma2Selected = ...
            coefD^2*asympt.completeCase.sigma2;
    end
    if isfield(asympt.completeCase,'crossTEMCC') && ...
            isfinite(asympt.completeCase.crossTEMCC)
        asympt.completeCase.crossTEMCCSelected = ...
            coefD^2*asympt.completeCase.crossTEMCC;
    end
end

% Degeneracy is unchanged by a strictly positive scalar multiplier. Keep
% the numerical decision made by the base analytical routine.
asympt.degenerate = baseDegenerate;

if asympt.degenerate
    asympt.TDmean.z = NaN;
    asympt.TDmean.pvalue = NaN;
    asympt.TLmean.z = NaN;
    asympt.TLmean.pvalue = NaN;
else
    asympt.TDmean.z = Tobs(4)/asympt.TDmean.se;
    asympt.TDmean.pvalue = erfc(abs(asympt.TDmean.z)/sqrt(2));
    asympt.TLmean.z = Tobs(2)/asympt.TLmean.se;
    asympt.TLmean.pvalue = erfc(abs(asympt.TLmean.z)/sqrt(2));
end
end

% -------------------------------------------------------------------------
function asympt = local_tem_asymptotic(Y, maskMiss, completeIdx, outFit, ...
    alpha, method, robustClass, robustEff, robustBonflev, Tobs, eps0)
%local_tem_asymptotic Analytical alpha>0 TEM sandwich benchmark.
%
% The implementation follows the scalar influence-function representation
%
%   J_eff' * b = aSigma,
%   -aSigma' * psi_TEM(W) = b' * xi(W),
%
% where xi(W)=vech{w_R(U_R-Sigma)}. The complete-case contribution to the
% MDP statistic is C*zeta_C/q. The robust complete-case scalar influence
% zeta_C is evaluated by a delete-one jackknife. Consequently the TEM part
% of the calculation is analytical, while the estimator-specific robust
% complete-case influence is evaluated numerically.
%
% The current Jacobian implementation is restricted to method='pri'. This
% is the parameter-free additive distance adjustment used in the theoretical
% derivation. In particular, detMap would require additional derivatives
% because its adjustment depends on Sigma.

[n,p] = size(Y);
qhat = sum(completeIdx)/n;
s = p*(p+1)/2;

asympt = local_unavailable_asymptotic(qhat);
asympt.reason = '';
asympt.VA = NaN;
asympt.mode = 'TEM influence-function sandwich';
asympt.theoryStatus = 'available';
asympt.TEM = struct('available',false,'sigma2',NaN,'threshold',NaN, ...
    'Jrcond',NaN,'nActivePatterns',0,'nPatterns',0);
asympt.completeCase = struct('influenceMethod','delete-one jackknife', ...
    'sigma2',NaN,'crossTEMCC',NaN,'nComplete',sum(completeIdx));

if ~strcmpi(method,'pri')
    asympt.reason = ['For alpha>0 the analytical TEM Jacobian is currently ' ...
        'implemented only for method=''pri''.'];
    return
end

if isempty(outFit) || ~isfield(outFit,'weights') || ...
        ~isfield(outFit,'loc') || ~isfield(outFit,'cov')
    asympt.reason = 'The fitted TEM object does not contain the required fields.';
    return
end

muHat = outFit.loc(:);
SigmaHat = (outFit.cov + outFit.cov')/2;
if any(~isfinite(SigmaHat(:))) || rcond(SigmaHat) <= 1e-12
    asympt.reason = 'The fitted TEM scatter matrix is numerically singular.';
    return
end

% Duplication matrix and derivative of log|Sigma|.
Dp = local_duplication_matrix(p);
K = SigmaHat\eye(p);
K = (K+K')/2;
aSigma = Dp' * K(:);

% Reconstruct the common adjusted-distance threshold used by TEM. For the
% pri adjustment, c = a_g + p-p_g. The kinfo values are from the last TEM
% concentration step and are therefore preferable to recomputing the
% threshold from rounded final distances. A final-distance fallback is kept
% for safety.
cthr = NaN;
if isfield(outFit,'kinfo') && istable(outFit.kinfo) && ~isempty(outFit.kinfo) && ...
        all(ismember({'pobs','athr'},outFit.kinfo.Properties.VariableNames))
    okc = isfinite(outFit.kinfo.athr) & outFit.kinfo.athr > 0;
    if any(okc)
        cvals = outFit.kinfo.athr(okc) + p - outFit.kinfo.pobs(okc);
        cthr = median(cvals(isfinite(cvals)));
    end
end

[d2part,poss] = mdPartialMD(Y,muHat,SigmaHat);
d2adj = mdPartialMD2full(d2part,p,poss,'method','pri');
w = outFit.weights(:) > 0;

if ~isfinite(cthr)
    if ~any(w & isfinite(d2adj))
        asympt.reason = 'Unable to reconstruct the final TEM trimming threshold.';
        return
    end
    cthr = max(d2adj(w & isfinite(d2adj)));
end

% Missingness patterns and empirical pattern probabilities.
[patt,~,ic] = unique(maskMiss,'rows');
G = size(patt,1);
asympt.TEM.nPatterns = G;

J = zeros(s,s);
active = false(G,1);
kgVec = NaN(G,1);
Lcell = cell(G,1);
Vcell = cell(G,1);
obsCell = cell(G,1);

% Build the effective scatter Jacobian pattern by pattern. For pri,
% a_g(c)=c-(p-p_g) and therefore dot(a_g)=1.
for g = 1:G
    obs = find(~patt(g,:));
    mis = find(patt(g,:));
    pg = numel(obs);
    pig = sum(ic==g)/n;
    obsCell{g} = obs;

    if pg == 0
        continue
    end

    ag = cthr - (p-pg);
    if ~isfinite(ag) || ag <= 0
        continue
    end

    gammag = chi2cdf(ag,pg);
    fg = chi2pdf(ag,pg);
    kg = chi2cdf(ag,pg+2)/gammag;

    if ~isfinite(kg) || kg <= 0 || ~isfinite(fg)
        continue
    end

    SigmaG = SigmaHat(obs,obs);
    SigmaG = (SigmaG+SigmaG')/2;
    if rcond(SigmaG) <= 1e-12
        asympt.reason = ['A pattern-specific covariance block is numerically ' ...
            'singular; TEM analytical benchmark not computed.'];
        return
    end
    BG = SigmaG\eye(pg);
    BG = (BG+BG')/2;

    Lg = zeros(p,pg);
    Lg(obs,:) = eye(pg);
    Vg = zeros(p,p);

    if ~isempty(mis)
        Ag = SigmaHat(mis,obs)/SigmaG;
        Cg = SigmaHat(mis,mis)-Ag*SigmaHat(obs,mis);
        Cg = (Cg+Cg')/2;
        Lg(mis,:) = Ag;
        Vg(mis,mis) = Cg;
    end

    ug = -gammag + 2*ag^2*fg/(pg*(pg+2)*kg);
    vg = fg*(ag^2/(pg*(pg+2)*kg)-ag/pg);
    LSigmaL = Lg*SigmaG*Lg';

    for h = 1:s
        H = reshape(Dp(:,h),p,p);
        Hg = H(obs,obs);
        tg = trace(BG*Hg);
        JH = ug*(Lg*Hg*Lg') + vg*tg*LSigmaL;
        J(:,h) = J(:,h) + pig*local_vech(JH);
    end

    active(g) = true;
    kgVec(g) = kg;
    Lcell{g} = Lg;
    Vcell{g} = Vg;
end

asympt.TEM.threshold = cthr;
asympt.TEM.nActivePatterns = sum(active);

if ~any(active)
    asympt.reason = 'No active missingness pattern is available at the TEM threshold.';
    return
end

Jrcond = rcond(J);
asympt.TEM.Jrcond = Jrcond;
if ~isfinite(Jrcond) || Jrcond <= 1e-12
    asympt.reason = ['The effective TEM scatter Jacobian is numerically ' ...
        'singular; analytical benchmark not computed.'];
    return
end

% b gives the scalar projection needed by TD,mean without constructing the
% full TEM covariance matrix.
b = J'\aSigma;

% Empirical values of xi_i=vech{w_i(U_i-Sigma)}. The theoretical exact
% pattern factor is used here, rather than the finite-sample shrinkage of
% unstable factors that mdTEM may apply internally.
Xi = zeros(n,s);
for i = 1:n
    if ~w(i)
        continue
    end
    g = ic(i);
    if ~active(g)
        continue
    end

    obs = obsCell{g};
    zi = Y(i,obs)' - muHat(obs);
    xhat = Lcell{g}*zi;
    Ui = (xhat*xhat')/kgVec(g) + Vcell{g};
    Xi(i,:) = local_vech(Ui-SigmaHat)';
end

% b'xi is the TEM contribution to the influence function of TD,mean,
% because TD,mean has the opposite sign of aSigma'psi_TEM.
temContribution = Xi*b;
temContribution = temContribution-mean(temContribution);
sigmaTEM2 = mean(temContribution.^2);
asympt.TEM.available = true;
asympt.TEM.sigma2 = sigmaTEM2;

% Evaluate the scalar influence of the robust complete-case scatter by a
% delete-one jackknife. The same random-number state is reset before every
% robust refit so that stochastic subsampling in MCD/MM does not inject
% avoidable algorithmic noise into the jackknife differences.
Ycc = Y(completeIdx,:);
[zetaC,ccReason] = local_cc_scalar_jackknife(Ycc,alpha,robustClass, ...
    robustEff,robustBonflev,aSigma);

if ~isempty(ccReason)
    asympt.reason = ccReason;
    return
end

ccContribution = zeros(n,1);
ccContribution(completeIdx) = zetaC/qhat;
ccContribution = ccContribution-mean(ccContribution);

% End-to-end influence of TD,mean:
%   zeta_D = b'xi + C*zeta_C/q.
zetaD = temContribution + ccContribution;
zetaD = zetaD-mean(zetaD);

sigmaCC2 = mean(ccContribution.^2);
crossTEMCC = mean(temContribution.*ccContribution);
sigmaD2 = mean(zetaD.^2);

asympt.completeCase.sigma2 = sigmaCC2;
asympt.completeCase.crossTEMCC = crossTEMCC;
asympt.completeCase.nComplete = size(Ycc,1);

if strcmpi(robustClass,'FS')
    asympt.theoryStatus = ['experimental for adaptive FS: the TEM influence ' ...
        'calculation is analytical, but the current FS signal/stopping rule ' ...
        'is outside the fixed-fraction FS asymptotic theorem.'];
else
    asympt.theoryStatus = ['TEM influence calculation is analytical; the ' ...
        'robust complete-case scalar influence is evaluated by delete-one ' ...
        'jackknife.'];
end

if ~isfinite(sigmaD2) || sigmaD2 < 0
    asympt.reason = 'The end-to-end TEM sandwich variance is not finite.';
    return
end

% Stabilized log-ratio coefficient.
kappa = integral(@(q) local_kappa_integrand(q,p,eps0),0,Inf, ...
    'RelTol',1e-10,'AbsTol',1e-12)/p;

tolVar = 1e-10*max(1,sigmaTEM2+sigmaCC2);
asympt.available = true;
asympt.reason = '';
asympt.qhat = qhat;
asympt.degenerate = sigmaD2 <= tolVar;

asympt.TDmean.sigma2 = sigmaD2;
asympt.TDmean.se = sqrt(sigmaD2/n);
asympt.TDmean.components = struct('TEM',sigmaTEM2, ...
    'completeCase',sigmaCC2,'twiceCross',2*crossTEMCC);

asympt.TLmean.kappa = kappa;
asympt.TLmean.sigma2 = kappa^2*sigmaD2;
asympt.TLmean.se = sqrt(asympt.TLmean.sigma2/n);

if asympt.degenerate
    asympt.TDmean.z = NaN;
    asympt.TDmean.pvalue = NaN;
    asympt.TLmean.z = NaN;
    asympt.TLmean.pvalue = NaN;
else
    asympt.TDmean.z = Tobs(4)/asympt.TDmean.se;
    asympt.TDmean.pvalue = erfc(abs(asympt.TDmean.z)/sqrt(2));
    asympt.TLmean.z = Tobs(2)/asympt.TLmean.se;
    asympt.TLmean.pvalue = erfc(abs(asympt.TLmean.z)/sqrt(2));
end
end

% -------------------------------------------------------------------------
function [zetaC,reason] = local_cc_scalar_jackknife(Ycc,alpha,robustClass, ...
    robustEff,robustBonflev,aSigma)
%local_cc_scalar_jackknife Delete-one estimate of the robust scatter IF.
%
% If theta_hat is asymptotically linear on ncc complete observations, then
% (ncc-1)*(theta_hat-theta_hat_{(-i)}) estimates the centered influence.
% After centering, the common full-sample estimate cancels, so it is enough
% to center the scalar leave-one-out estimates themselves.

ncc = size(Ycc,1);
zetaC = [];
reason = '';

if ncc < 5
    reason = 'Too few complete observations for the robust delete-one influence.';
    return
end

looScalar = NaN(ncc,1);
rngState = rng;
cleanupObj = onCleanup(@() rng(rngState)); %#ok<NASGU>

for i = 1:ncc
    keep = true(ncc,1);
    keep(i) = false;
    try
        % Reset the state to pair stochastic robust fits across deletions.
        rng(rngState)
        [~,SigMinus] = local_complete_case_fit(Ycc(keep,:),alpha, ...
            robustClass,robustEff,robustBonflev);
        if any(~isfinite(SigMinus(:)))
            reason = sprintf(['Nonfinite robust scatter in delete-one fit ' ...
                'for complete observation %d.'],i);
            zetaC = [];
            return
        end
        looScalar(i) = aSigma' * local_vech((SigMinus+SigMinus')/2);
    catch ME
        reason = sprintf(['Robust delete-one fit failed for complete ' ...
            'observation %d: %s'],i,ME.message);
        zetaC = [];
        return
    end
end

if any(~isfinite(looScalar))
    reason = 'The robust delete-one influence contains nonfinite values.';
    return
end

zetaC = -(ncc-1)*(looScalar-mean(looScalar));
end

% -------------------------------------------------------------------------
function v = local_vech(S)
%local_vech Lower-triangular half-vectorization, column by column.

p = size(S,1);
v = zeros(p*(p+1)/2,1);
k = 0;
for j = 1:p
    for i = j:p
        k = k+1;
        v(k) = S(i,j);
    end
end
end

% -------------------------------------------------------------------------
function asympt = local_classical_asymptotic(maskMiss, SigmaHat, ...
    nComplete, Tobs, eps0)
%local_classical_asymptotic Analytical benchmark for alpha=0.
%
% Under Gaussian MCAR,
%   sqrt(n)*TDmean -> N(0,sigmaD2),
% where
%   sigmaD2 = 2*p/qhat - aSigma'*(IA\aSigma).
% The mean stabilized log-ratio has asymptotic variance
%   sigmaL2 = kappa^2*sigmaD2,
% with
%   kappa = E{Q/(Q+eps0)}/p, Q~chi2_p.

[n,p] = size(maskMiss);
qhat = nComplete/n;
s = p*(p+1)/2;

asympt = local_unavailable_asymptotic(qhat);
asympt.reason = '';
asympt.VA = NaN;

% Duplication matrix and derivative of log|Sigma|.
Dp = local_duplication_matrix(p);
SigmaHat = (SigmaHat + SigmaHat')/2;
K = SigmaHat\eye(p);
K = (K + K')/2;

% a_\Sigma=D_p^{\mathsf T}\operatorname{vec}(\Sigma^{-1}),
% Equation (25)
aSigma = Dp' * K(:);

% Observed-data covariance information, averaged over the empirical
% missingness-pattern distribution.
IA = zeros(s,s);
[patt,~,ic] = unique(maskMiss,'rows');
idxFull = reshape(1:p*p,p,p);

for g = 1:size(patt,1)
    obs = ~patt(g,:);
    pg = sum(obs);
    if pg == 0
        continue
    end

    ng = sum(ic == g);
    pig = ng/n;

    SigmaG = SigmaHat(obs,obs);
    SigmaG = (SigmaG + SigmaG')/2;
    BG = SigmaG\eye(pg);
    BG = (BG + BG')/2;

    % Dg maps vech(Sigma) to vec(Sigma_g).
    idxG = idxFull(obs,obs);
    Dg = Dp(idxG(:),:);
    Ir = 0.5 * Dg' * kron(BG,BG) * Dg;
    IA = IA + pig*Ir;
end

IA = (IA + IA')/2;

if rcond(IA) <= 1e-12
    asympt.reason = ['Observed-data covariance information is numerically ' ...
        'singular; analytical benchmark not computed.'];
    return
end

v = IA\aSigma;
VA = aSigma' * v;
sigmaD2 = 2*p/qhat - VA;

% The theoretical variance is nonnegative. Remove only tiny negative values
% caused by floating-point roundoff.
tolVar = 1e-10 * max(1,2*p/qhat);
if sigmaD2 < 0 && sigmaD2 >= -tolVar
    sigmaD2 = 0;
elseif sigmaD2 < -tolVar
    asympt.reason = ['Estimated analytical variance is negative beyond ' ...
        'numerical tolerance; analytical benchmark not computed.'];
    asympt.VA = VA;
    return
end

% Stabilized log-ratio coefficient. The local integrand handles the
% endpoints explicitly, avoiding the indeterminate product 0*Inf for p<2.
kappa = integral(@(q) local_kappa_integrand(q,p,eps0), 0, Inf, ...
    'RelTol',1e-10,'AbsTol',1e-12) / p;

asympt.available = true;
asympt.reason = '';
asympt.qhat = qhat;
asympt.VA = VA;
asympt.degenerate = sigmaD2 <= tolVar;

asympt.TDmean.sigma2 = sigmaD2;
asympt.TDmean.se = sqrt(sigmaD2/n);

asympt.TLmean.kappa = kappa;
asympt.TLmean.sigma2 = kappa^2*sigmaD2;
asympt.TLmean.se = sqrt(asympt.TLmean.sigma2/n);

if asympt.degenerate
    asympt.TDmean.z = NaN;
    asympt.TDmean.pvalue = NaN;
    asympt.TLmean.z = NaN;
    asympt.TLmean.pvalue = NaN;
else
    % Tobs(4) is mean difference; Tobs(2) is mean log-ratio.
    asympt.TDmean.z = Tobs(4)/asympt.TDmean.se;
    asympt.TDmean.pvalue = erfc(abs(asympt.TDmean.z)/sqrt(2));

    asympt.TLmean.z = Tobs(2)/asympt.TLmean.se;
    asympt.TLmean.pvalue = erfc(abs(asympt.TLmean.z)/sqrt(2));
end
end

% -------------------------------------------------------------------------
function y = local_kappa_integrand(q,p,eps0)
%local_kappa_integrand Integrand for E{Q/(Q+eps0)}, Q~chi2_p.

ratio = q./(q+eps0);
ratio(q == 0) = 0;
ratio(isinf(q)) = 1;
y = ratio.*chi2pdf(q,p);
y(~isfinite(y)) = 0;
end

% -------------------------------------------------------------------------
function asympt = local_unavailable_asymptotic(qhat)
%local_unavailable_asymptotic Initialize analytical benchmark output.

asympt = struct;
asympt.available = false;
asympt.reason = 'Analytical benchmark is not available for this configuration.';
asympt.qhat = qhat;
asympt.VA = NaN;
asympt.degenerate = false;
asympt.aggregation = '';
asympt.filterCutoff = NaN;
asympt.kc = NaN;
asympt.lambda = NaN;
asympt.gammaFilter = NaN;
asympt.baseSigmaD2 = NaN;
asympt.baseKappa = NaN;
asympt.TDmean = struct('coefficient',NaN,'sigma2',NaN,'se',NaN, ...
    'z',NaN,'pvalue',NaN);
asympt.TLmean = struct('kappa',NaN,'lambda',NaN,'coefficient',NaN, ...
    'sigma2',NaN,'se',NaN,'z',NaN,'pvalue',NaN);
end

% -------------------------------------------------------------------------
function D = local_duplication_matrix(p)
%local_duplication_matrix Duplication matrix for column-wise vech.
%
% D satisfies vec(H)=D*vech(H) for every symmetric p-by-p matrix H, where
% vech stacks the lower triangular part column by column.

s = p*(p+1)/2;
D = zeros(p*p,s);
k = 0;
for j = 1:p
    for i = j:p
        k = k + 1;
        E = zeros(p,p);
        E(i,j) = 1;
        if i ~= j
            E(j,i) = 1;
        end
        D(:,k) = E(:);
    end
end
end

% -------------------------------------------------------------------------
function Sspd = local_make_spd(S)
% Make covariance matrix symmetric positive definite if needed.

S = (S + S')/2;
[~,flag] = chol(S);

if flag == 0
    Sspd = S;
    return
end

lam = 1e-8 * trace(S) / size(S,1);
if lam <= 0 || ~isfinite(lam)
    lam = 1e-8;
end

I = eye(size(S));
for k = 1:8
    Stmp = S + lam*I;
    [~,flag] = chol(Stmp);
    if flag == 0
        Sspd = Stmp;
        return
    end
    lam = 10*lam;
end

error('FSDA:mdMDPtest:NonSPD', ...
    'Unable to regularize covariance matrix to positive definiteness.');
end


% -------------------------------------------------------------------------
function out = local_coupled_tem(Y,forcedZero,alpha,method,tol)
%local_coupled_tem Pattern-corrected TEM with fixed forced-zero rows.
%
% This local sensitivity implementation mirrors the pattern-corrected mdTEM
% concentration step. At each iteration it first constructs TEM's ordinary
% indicator wT from the h=floor(n*(1-alpha)) smallest adjusted distances and
% then applies the fixed eligibility constraint
%
%   w = wT .* (~forcedZero).
%
% Consequently a complete case rejected by the robust complete-case fit can
% never contribute to the TEM E/M step. The internal threshold is still the
% h-th ordinary TEM adjusted-distance threshold. Because the extra exclusion
% is not generated solely by that radial threshold, the usual pattern-wise
% Tallis correction is only a working correction in this sensitivity mode;
% no analytical reference law is attached to the resulting MDP statistic.

[n,p] = size(Y);
forcedZero = logical(forcedZero(:));
if numel(forcedZero) ~= n
    error('FSDA:mdMDPtest:WrongCoupledMask', ...
        'The coupled TEM mask must have one entry for every row of Y.');
end

maxiter = 100;
tol_sigma = true;
method = string(method);

% Same initialization used by mdTEM.
mus = mean(Y,1,"omitmissing")';
X0 = Y;
for j = 1:p
    miss = isnan(X0(:,j));
    X0(miss,j) = mus(j);
end
sigs = cov(X0,1);

keep_count = max(0,floor(n*(1-alpha)));
nanY = isnan(Y);
w = zeros(n,1);
wT = zeros(n,1);
kinfo = [];
kfactor = 1;
dif = Inf;
iter = 0;

while dif > tol && iter < maxiter
    iter = iter + 1;
    mus_old = mus;
    sigs_old = sigs;

    if method == "impMD"
        Yimp = mdImputeCondMean(Y,mus,sigs);
        d2_adj = mahalFS(Yimp,mus',sigs);
        poss = sum(~nanY,2);
    elseif method == "detMap"
        [d2,poss] = mdPartialMD(Y,mus,sigs);
        d2_adj = mdPartialMD2full(d2,p,poss,'method',method, ...
            'Y',Y,'Sigma',sigs);
    else
        [d2,poss] = mdPartialMD(Y,mus,sigs);
        d2_adj = mdPartialMD2full(d2,p,poss,'method',method);
    end

    nanMask = isnan(d2_adj);
    [~,idxSorted] = sort(d2_adj,'ascend','MissingPlacement','last');
    keepIdx = idxSorted(1:min(keep_count,sum(~nanMask)));
    if isempty(keepIdx)
        error('FSDA:mdMDPtest:NoTEMRows', ...
            'No finite row is available for the coupled TEM concentration step.');
    end

    % Ordinary TEM concentration indicator and coupled product weight.
    wT = zeros(n,1);
    wT(keepIdx) = 1;
    w = wT;
    w(forcedZero) = 0;
    if ~any(w)
        error('FSDA:mdMDPtest:NoTEMRows', ...
            'The coupled TEM constraint removes all retained rows.');
    end

    % Keep the ordinary TEM threshold. The additional fixed exclusion is
    % represented only through w, exactly as in the product definition.
    cthr = max(d2_adj(keepIdx));

    kinfo = local_coupled_pattern_factors(nanY,w,poss,cthr,p,n,sigs,method);
    [mus,sigs,kfactor] = local_coupled_corrected_step(Y,nanY,w,mus, ...
        sigs,kinfo);

    muDiff = max(abs(mus(:)-mus_old(:)));
    sigmaDiff = max(abs(sigs(:)-sigs_old(:)));
    if tol_sigma
        dif = max(muDiff,sigmaDiff);
    else
        dif = muDiff;
    end
end

out = struct;
out.loc = mus;
out.cov = sigs;
out.iter = iter;
out.weights = w;
out.internalWeights = wT;
out.forcedZero = forcedZero;
out.kfactor = kfactor;
out.kinfo = kinfo;
end

% -------------------------------------------------------------------------
function kinfo = local_coupled_pattern_factors(nanY,w,poss,cthr,p,n,Sigma,method)
%local_coupled_pattern_factors Pattern-wise Tallis factors for coupled TEM.

keep = w > 0;
patt = unique(nanY(keep,:),'rows');
G = size(patt,1);
pobs = zeros(G,1);
nkept = zeros(G,1);
athr = zeros(G,1);
gammag = zeros(G,1);
kg = nan(G,1);

for g = 1:G
    pg = patt(g,:);
    rows = keep & all(nanY == pg,2);
    obs = ~pg;
    pobs(g) = sum(obs);
    nkept(g) = sum(rows);

    if pobs(g) == 0
        athr(g) = 0;
        gammag(g) = 0;
        kg(g) = NaN;
        continue
    end

    athr(g) = local_coupled_inv_adjust(cthr,pobs(g),p,n,Sigma,obs,method);
    if athr(g) <= 0 || ~isfinite(athr(g))
        gammag(g) = 0;
        kg(g) = NaN;
    else
        gammag(g) = chi2cdf(athr(g),pobs(g));
        if gammag(g) <= 0
            kg(g) = NaN;
        else
            kg(g) = chi2cdf(athr(g),pobs(g)+2)/gammag(g);
        end
    end
end

kinfo = table(pobs,nkept,athr,gammag,kg);
kbar = local_coupled_weighted_factor(kinfo);
bad = ~isfinite(kinfo.kg) | kinfo.kg <= 0 | kinfo.kg > 1 | ...
    kinfo.nkept < max(2,kinfo.pobs);
kinfo.kg(bad) = kbar;
end

% -------------------------------------------------------------------------
function kbar = local_coupled_weighted_factor(kinfo)
%local_coupled_weighted_factor Information-weighted pattern factor.

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

% -------------------------------------------------------------------------
function a = local_coupled_inv_adjust(c,pg,p,n,Sigma,obs,method)
%local_coupled_inv_adjust Inverse adjusted-distance map for one pattern.

method = string(method);
switch method
    case "pri"
        a = c-(p-pg);
    case "expScale"
        a = c*pg/p;
    case "zMap"
        a = pg+sqrt(pg/p)*(c-p);
    case "detMap"
        gfull = exp(local_coupled_logdet_spd(Sigma)/p);
        gobs = exp(local_coupled_logdet_spd(Sigma(obs,obs))/pg);
        a = c*(pg/p)*(gobs/gfull);
    case "chiMap"
        u = chi2cdf(c,p);
        u = min(max(u,eps),1-eps);
        a = chi2inv(u,pg);
    case "betaMap"
        cn = (n-1)^2/n;
        if n <= p+1 || n <= pg+1
            a = c;
            return
        end
        u = min(max(c/cn,0),1-eps);
        al = betacdf(u,p/2,(n-p-1)/2);
        al = min(max(al,eps),1-eps);
        a = cn*betainv(al,pg/2,(n-pg-1)/2);
    case "impMD"
        a = c;
    otherwise
        a = c;
end
end

% -------------------------------------------------------------------------
function [musNew,sigsNew,kbar] = local_coupled_corrected_step(Y,nanY,w, ...
    mus,sigs,kinfo)
%local_coupled_corrected_step Pattern-wise corrected E/M update.

p = size(Y,2);
mus = mus(:);
keep = w > 0;
S = zeros(p,p);
m1 = zeros(p,1);
h = 0;
patt = unique(nanY(keep,:),'rows');

for g = 1:size(patt,1)
    pg = patt(g,:);
    rows = keep & all(nanY == pg,2);
    hg = sum(rows);
    if hg == 0
        continue
    end

    o = find(~pg);
    m = find(pg);
    if isempty(o)
        continue
    end

    idx = find(kinfo.pobs == numel(o) & kinfo.nkept == hg,1);
    if isempty(idx)
        kg = 1;
    else
        kg = kinfo.kg(idx);
    end
    if ~isfinite(kg) || kg <= 0
        kg = 1;
    end

    Z = Y(rows,o)-mus(o)';
    Szz = (Z'*Z)/kg;
    sZ = sum(Z,1)';

    S(o,o) = S(o,o)+Szz;
    m1(o) = m1(o)+sZ;

    if ~isempty(m)
        Soo = sigs(o,o);
        Ag = sigs(m,o)/Soo;
        Cg = sigs(m,m)-Ag*sigs(o,m);
        Cg = (Cg+Cg')/2;

        SzzAg = Szz*Ag';
        S(o,m) = S(o,m)+SzzAg;
        S(m,o) = S(m,o)+SzzAg';
        S(m,m) = S(m,m)+Ag*SzzAg+hg*Cg;
        m1(m) = m1(m)+Ag*sZ;
    end
    h = h+hg;
end

if h == 0
    musNew = mus;
    sigsNew = sigs;
    kbar = 1;
    return
end

m1 = m1/h;
S = S/h;
musNew = mus+m1;
sigsNew = S-(m1*m1');
sigsNew = (sigsNew+sigsNew')/2;
kbar = local_coupled_weighted_factor(kinfo);
end

% -------------------------------------------------------------------------
function val = local_coupled_logdet_spd(S)
%local_coupled_logdet_spd Stable log determinant for a positive scatter matrix.

S = (S+S')/2;
[R,flag] = chol(S);
if flag == 0
    val = 2*sum(log(diag(R)));
else
    val = log(max(det(S),realmin));
end
end

%FScategory:MULT-MissingData