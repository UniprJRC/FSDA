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
%                            benchmark for the two mean statistics. For
%                            alpha=0, out.asympt.available is true when the
%                            analytical variance is numerically available.
%                            The structure contains qhat, VA, and the fields
%                            TDmean and TLmean. Each statistic contains the
%                            asymptotic variance of sqrt(n)T (sigma2), the
%                            standard error of T (se), the studentized value
%                            (z), and the two-sided asymptotic p-value. The
%                            TLmean structure also contains kappa. For
%                            alpha>0, out.asympt.available is false because
%                            the robust TEM analytical variance is not yet
%                            implemented in this routine.
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
%    2) mean  ( log((d2_all + eps0) ./ (d2_cc + eps0)) );
%    3) median( d2_all - d2_cc );
%    4) mean  ( d2_all - d2_cc );
%
%  Small p-values indicate that the change in distances is larger than what
%  is expected under the MCAR bootstrap model.
%
%  For alpha=0 the function also reports an analytical Gaussian benchmark
%  for the mean difference and mean log-ratio. The bootstrap p-values in
%  out.pvalue remain the primary finite-sample calibration.
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

% Distances based on EM/TEM fit using all rows
[d2_all_cc, muHat, SigHat] = local_fit_and_get_complete_distances( ...
    Y, completeIdx, alpha, method, tol);

% Observed statistics
Tobs = local_statistic(d2_cc, d2_all_cc,eps0);

% Analytical Gaussian benchmark for the two mean statistics. The
% bootstrap remains the primary finite-sample calibration.
if alpha == 0
    asympt = local_classical_asymptotic(maskMiss, SigHat, nComplete, ...
        Tobs, eps0);
else
    asympt = local_unavailable_asymptotic(nComplete/n);
end

% Bootstrap under MCAR
Tboot = NaN(nsimul,4);

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

    [muCCStar,SigCCStar] = local_complete_case_fit(YccStar,alpha, ...
        robustClass,robustEff,robustBonflev);

 
    d2_cc_star = mahalFS(YccStar, muCCStar, SigCCStar);

    % EM/TEM distances for the same complete rows
    d2_all_cc_star = local_fit_and_get_complete_distances( ...
        Ystar, completeIdx, alpha, method, tol);

    % Store the four statistics
    Tboot(j,:) = local_statistic(d2_cc_star, d2_all_cc_star, eps0);
end

% Remove bootstrap samples containing NaNs
Tboot = Tboot(all(~isnan(Tboot),2),:);

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
    out.nOutliersCC = numel(outRob.outliers);
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
    title(['\alpha=' num2str(alpha)])
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


function [d2_all_cc, muHat, SigHat] = local_fit_and_get_complete_distances( ...
    Y, completeIdx, alpha, method, tol)
% Compute distances for complete rows after fitting EM/TEM on all data.

p = size(Y,2);

if alpha == 0
    outFit = mdEM(Y);
else
    outFit = mdTEM(Y,'method',method,'alpha',alpha,'tol',tol,'consistencyfactor','pattern');
end

muHat = outFit.loc;
SigHat = outFit.cov;

[d2_part, poss] = mdPartialMD(Y, muHat, SigHat);
d2_full = mdPartialMD2full(d2_part, p, poss, 'method', method);
d2_all_cc = d2_full(completeIdx);

end

% -------------------------------------------------------------------------
function T = local_statistic(d2_cc, d2_all, eps0)
%local_statistic Compute the four MDP test statistics.

rat = log((d2_all + eps0) ./ (d2_cc + eps0));
dif = d2_all - d2_cc;

T1 = median(rat);
T2 = mean(rat);
T3 = median(dif);
T4 = mean(dif);

T = [T1 T2 T3 T4];
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
asympt.reason = ['Analytical benchmark currently implemented only for ' ...
    'the classical alpha=0 case.'];
asympt.qhat = qhat;
asympt.VA = NaN;
asympt.degenerate = false;
asympt.TDmean = struct('sigma2',NaN,'se',NaN,'z',NaN,'pvalue',NaN);
asympt.TLmean = struct('kappa',NaN,'sigma2',NaN,'se',NaN, ...
    'z',NaN,'pvalue',NaN);
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

%FScategory:MULT-MissingData