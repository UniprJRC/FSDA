function out = benford(x, varargin)
%Benford-law analysis for leading digits, first-two digits, and mantissa
%
%<a href="matlab: docsearchFS('benford')">Link to the help page for this function</a>
%
%   benford analyzes a positive numeric sample according to Benford's law.
%   The function extracts the first digit, second digit, first two digits,
%   and mantissa, compares the observed frequencies with the Benford
%   theoretical frequencies, and returns chi-square, two-stage, and
%   mantissa-based test results. It can also produce comparison plots for
%   the selected digit groups. Before the analysis, invalid observations
%   are removed and the digit extraction is checked for consistency using
%   independent methods.
%
% Required input arguments:
%
%       x:  Numeric vector of observations. Positive values are used in
%           the analysis.
%           Data Types - single|double
%
% Optional input arguments:
%
%    Digits:  Numeric vector specifying the digit groups to analyze.
%             Allowed values are 1 (first digit), 2 (second digit) and
%             12 (first two digits). Default is [1 2 12].
%             Example - 'Digits',[1 12]
%             Data Types - double
%
%     Plots:  Logical or numeric scalar, or text, specifying the requested
%             plots. The plots are drawn by function benfordplot at the end
%             of the analysis, so those which accompany a test need that
%             test in option Tests. True (or 1) stands for "all", false (or
%             0) for no plot. Valid labels (not case sensitive) are:
%               - "first", "second", "firstTwo": observed against expected
%                 counts of the digit group (need the group in Digits);
%               - "deviation": observed minus expected proportions with
%                 simultaneous bands and the MAD conformity category
%                 (Nigrini, 2012);
%               - "heatmap": Pearson residuals of the 90 first-two-digit
%                 cells (needs 12 in Digits);
%               - "magnitude": histogram of log10(x) and counts per order
%                 of magnitude;
%               - "mantissa": mantissa histogram and ECDF difference with
%                 the KS band and the Kuiper strip;
%               - "duplication": duplicate excess by value type;
%               - "twostage", "savvy", "lacunary": plots of the
%                 corresponding test;
%               - "all": all the plots above whose ingredients are
%                 available.
%             See benfordplot for the description of each plot. Default is
%             ["first","second"]. To redraw plots, or to change the level
%             of their bands, call benfordplot with the output of benford.
%             Example - 'Plots',["deviation","heatmap"]
%             Data Types - logical|double|string|char|cell
%
%     Tests:  Text specifying the requested tests. Implemented tests are
%             "chi2", "twostage", "mantissa", "savvy" and "lacunary".
%             Specifically:
%               - "chi2"     : Chi-square test for the selected digit groups
%               - "twostage" : Two-stage Benford test based on MAD and chi2
%               - "mantissa" : Mantissa-based tests for uniformity (arc,
%                              mean, variance, joint, Kolmogorov-Smirnov,
%                              Watson U2, d* and Kuiper)
%               - "savvy"    : Tests for Benford-savvy (adversarial)
%                              manipulation, Barabesi et al. (2026): KS2
%                              and KU2 on the fractional part of the
%                              significand, Q_Delta = Q2 - Q1 and the
%                              combined tests G_KS and G_KU (min p-value).
%                              Q1 (Pearson) and Q2 (Hotelling) are also
%                              reported.
%               - "lacunary" : Single-observation lacunary harmonic test of
%                              Ferrari (2026), applied to the base-10
%                              mantissa of EACH observation truncated to
%                              m = LacunaryM digits. It returns one
%                              statistic and one chi2(18) p-value per
%                              observation, and flags observations whose
%                              multiscale digit structure is atypical.
%             Default is ["chi2","twostage"].
%             Example - 'Tests',["chi2","mantissa","savvy"]
%             Data Types - string|char|cell
%
%     Alpha:  Significance level for the chi-square tests. Scalar in (0,1).
%             Default value is 0.05.
%             Example - 'Alpha',0.01
%             Data Types - single|double
%
%         B:  Number of Monte Carlo replications for the mantissa-based
%             and Benford-savvy tests. Positive integer or empty.
%             If empty (default) it is chosen automatically from options
%             MCHalfWidth, MCCL and MCPRef.
%             Example - 'B',10000
%             Data Types - single|double
%
% LacunaryM:  Number of mantissa digits m used by the lacunary test.
%             Integer between 1 and 15. Default is 10. The mantissa of a
%             double is reliable to about 13 digits only, so values of m
%             above 12 are not recommended. The chi2 calibration of the
%             test is asymptotic in m (see Ferrari, 2026).
%             Example - 'LacunaryM',12
%             Data Types - double
%
%      PVTS:  Percentile used for the two-stage critical values. Scalar
%             in (0,100). Default is 5.
%             Example - 'PVTS',1
%             Data Types - single|double
%
%   Display:  Logical scalar controlling printed output. If true, a
%             summary of each test result is printed to the command window.
%             Default is true.
%             Example - 'Display',false
%             Data Types - logical
%
% Output:
%
%     out:  Struct containing the results of the Benford analysis.
%           Fields are created only for the selected digit groups and
%           requested tests. The fields are:
%
%             out.data      = observations used in the analysis.
%             out.n         = sample size.
%             out.options   = options used (derived internal flags
%                             included).
%             out.counts    = structure with fields first, second and
%                             firstTwo: observed counts of the digit
%                             groups in Digits (columns 1:9, 0:9, 10:99).
%             out.expected  = structure with the same fields: expected
%                             counts under the Benford law.
%             out.chi2      = (Tests "chi2") structure with fields first,
%                             second and firstTwo, each with fields stat,
%                             df, pValue and reject.
%             out.twoStage  = (Tests "twostage") structure with fields
%                             labels, cvMad, cvChi2, obsMad, obsChi2,
%                             rejectMad, rejectChi2, pValueMad,
%                             pValueChi2 (rows ordered as labels), and
%                             madSim and chi2Sim (B-by-3 null samples
%                             used by the plot).
%             out.mantissa  = n-by-1 vector of the mantissae
%                             frac(log10(x)).
%             out.mantissaTests = (Tests "mantissa") structure with fields
%                             arc, mean, variance, joint, ks, u2, dstar
%                             and kuiper, each with fields stat, pAsym
%                             (asymptotic p-value) and pMC (Monte Carlo
%                             p-value).
%             out.savvy     = (Tests "savvy") structure with fields n, B,
%                             Q1, Q2, KS2, KU2, QDelta, GKS and GKU, each
%                             with fields stat, pAsym and pMC (pAsym is
%                             NaN for GKS and GKU, whose statistic is a
%                             p-value and which reject for small values),
%                             and nullStats (B-by-4 null sample with
%                             columns KS2, KU2, Q1, Q2).
%             out.lacunary  = (Tests "lacunary") structure with fields m,
%                             base, df, alpha, T and pValue (one entry per
%                             observation), reject, idxReject, nReject,
%                             nExpected and pBinom (rough dataset-level
%                             summary, not part of the original paper).
%             out.duplication = (only if requested in Plots) duplication
%                             diagnostics.
%             out.extractionCheck, out.countsCheck = consistency checks of
%                             the digit extraction (option CheckExtraction).
%             out.figures   = handles of the figures drawn by benfordplot.
%
% See also: benfordplot
%
% References:
%
% The references are grouped by the test (or option) they refer to.
%
% Benford's law, general background
%
%    Benford, F. (1938). The law of anomalous numbers. Proceedings of the
%    American Philosophical Society, 78(4), 551-572.
%
%    Berger, A. and Hill, T.P. (2015). An Introduction to Benford's Law.
%    Princeton University Press, Princeton, NJ.
%
%    Nigrini, M.J. (2012). Benford's Law: Applications for Forensic
%    Accounting, Auditing, and Fraud Detection. Wiley, Hoboken, NJ.
%
%    Barabesi, L., Cerioli, A. and Di Marzio, M. (2023). Statistical models
%    and the Benford hypothesis: a unified framework. TEST, 32, 1479-1507.
%    doi:10.1007/s11749-023-00881-y
%
% Tests = "chi2" (chi-square test on first, second and first-two digits)
%
%    Pearson, K. (1900). On the criterion that a given system of
%    deviations from the probable in the case of a correlated system of
%    variables is such that it can be reasonably supposed to have arisen
%    from random sampling. Philosophical Magazine, Series 5, 50(302),
%    157-175.
%
%    Nigrini, M.J. (2012), see above (use of the chi-square statistic on
%    the digit counts, and of the MAD conformity ranges).
%
% Tests = "twostage" (two-stage test based on MAD and chi2)
%
%    Barabesi, L., Cerasa, A., Cerioli, A. and Perrotta, D. (2018).
%    Goodness-of-fit testing for the Newcomb-Benford law with application
%    to the detection of customs fraud. Journal of Business and Economic
%    Statistics, 36(2), 346-358. doi:10.1080/07350015.2016.1172014
%
%    Cerioli, A., Barabesi, L., Cerasa, A., Menegatti, M. and Perrotta, D.
%    (2019). Newcomb-Benford law and the detection of frauds in
%    international trade. Proceedings of the National Academy of Sciences
%    of the USA, 116(1), 106-115. doi:10.1073/pnas.1806617115
%
% Tests = "mantissa" (tests on the mantissa m = frac(log10(x)), uniform on
% [0,1) under the Benford hypothesis)
%
%   arc test (Rayleigh test of uniformity on the circle):
%    Alexander, J.C. (2009). Remarks on the use of Benford's Law. Working
%    Paper, Case Western Reserve University. SSRN 1505147.
%    doi:10.2139/ssrn.1505147
%
%    Rayleigh, Lord (1919). On the problem of random vibrations, and of
%    random flights in one, two, or three dimensions. Philosophical
%    Magazine, Series 6, 37(220), 321-347.
%
%    Mardia, K.V. and Jupp, P.E. (2000). Directional Statistics. Wiley,
%    Chichester.
%
%   mean, variance and joint (mean and variance) tests: normal tests on
%   the sample mean and variance of the mantissa (null values 1/2 and
%   1/12), in the spirit of the digit-based conformity tests of
%    Cerqueti, R. and Lupi, C. (2021). Some new tests of conformity with
%    Benford's law. Stats, 4(3), 745-761. doi:10.3390/stats4030044
%
%   Kolmogorov-Smirnov test:
%    Kolmogorov, A.N. (1933). Sulla determinazione empirica di una legge
%    di distribuzione. Giornale dell'Istituto Italiano degli Attuari, 4,
%    83-91.
%
%   Watson U2 test:
%    Watson, G.S. (1961). Goodness-of-fit tests on a circle. Biometrika,
%    48(1/2), 109-114.
%
%   Kuiper test:
%    Kuiper, N.H. (1960). Tests concerning random points on a circle.
%    Proceedings of the Koninklijke Nederlandse Akademie van
%    Wetenschappen, Series A, 63, 38-47.
%
%   d* test (Euclidean distance on the first-digit distribution,
%   normalised by its maximum):
%    Tam Cho, W.K. and Gaines, B.J. (2007). Breaking the (Benford) law:
%    statistical fraud detection in campaign finance. The American
%    Statistician, 61(3), 218-223. doi:10.1198/000313007X223496
%
%    Morrow, J. (2014). Benford's law, families of distributions and a
%    test basis. CEP Discussion Paper No. 1291, Centre for Economic
%    Performance, London School of Economics and Political Science.
%
%   Finite-sample corrections used for the asymptotic p-values of the KS,
%   Kuiper and Watson U2 tests:
%    Stephens, M.A. (1970). Use of the Kolmogorov-Smirnov, Cramer-von
%    Mises and related statistics without extensive tables. Journal of
%    the Royal Statistical Society, Series B, 32(1), 115-122.
%
% Tests = "savvy" (KS2, KU2, Q1, Q2, Q_Delta, G_KS, G_KU)
%
%    Barabesi, L., Cerioli, A., Cerasa, A. and Perrotta, D. (2026).
%    Statistical Detection of Adversarial Compliance With Benford's Law.
%    Statistical Analysis and Data Mining, 19(3), e70081.
%    doi:10.1002/sam.70081 (preprint: Robust inference under Benford's
%    law, arXiv:2507.08650, 2025)
%
%    Barabesi, L., Cerasa, A., Cerioli, A. and Perrotta, D. (2022). On
%    characterizations and tests of Benford's law. Journal of the American
%    Statistical Association, 117(540), 1887-1903.
%    doi:10.1080/01621459.2021.1891927
%    (sum-invariance characterization and Hotelling statistic Q2)
%
%    Kuiper (1960), Kolmogorov (1933) and Stephens (1970), see above (KS2
%    and KU2 are the KS and Kuiper statistics on the fractional part of
%    the significand).
%
% Tests = "lacunary" (single-observation lacunary harmonic test) and
% option LacunaryM
%
%    Ferrari, D. (2026). Single-Observation Uniformity Testing under
%    Increasing Precision via Lacunary Harmonics. arXiv:2604.27772
%    [stat.ME], submitted 30 April 2026.
%
%    Salem, R. and Zygmund, A. (1947). On lacunary trigonometric series.
%    Proceedings of the National Academy of Sciences of the USA, 33(11),
%    333-338. (lacunary central limit theorem underlying the chi2 null
%    distribution)
%
% Digit extraction, duplication diagnostics and plots implement no
% specific published procedure.
%
% Internal dependencies:
%
%   Digit extraction:
%       ExtractDigits
%       ExtractDigitsAC
%       ExtractDigitsStr
%       ExtractMantissa
%
%   Statistical distributions and tests:
%       benfordProbabilities
%       Chi2onDigitDistrib
%       benfordPatternSupport
%       MarginalDistrib
%       SubHypotheses
%       TwoStageTestCV
%       mantissa_tests
%       mc_pvalue
%
%   Kuiper and Benford-savvy tests, lacunary test:
%       kolmogorov_pvalue
%       kuiper_pvalue
%       chi2_upper
%       upper_tail_p
%       savvy_cdf
%       savvy_null
%       savvy_tests
%       lacunary_tests
%
%   Auxiliary tools:
%       normalizeTextList
%       countDigits
%       duplication_tests
%
%   Plotting (separate file):
%       benfordplot
%
%   External non-basic MATLAB dependencies:
%       mnrnd
%       combnk
%       chi2cdf
%       prctile
%       uniquetol
%
% Acknowledgements:
%
% Copyright 2026.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('benford')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% A clean Benford-like sample, with the default tests and plots.
    % The amounts are log-uniform over three orders of magnitude, hence
    % exactly Benford. By default the chi-square and two-stage tests are
    % run and the first and second digit plots are drawn. None of the
    % tests should reject.
    rng(1)
    V = 10.^(3*rand(2000,1));
    out = benford(V);
%}
%
%{
    % Selecting digit groups, significance level and plots.
    % All three digit groups are analysed with the chi-square test at
    % level 0.01, and three bar charts are drawn.
    rng(1)
    V = 10.^(3*rand(2000,1));
    out = benford(V,'Digits',[1 2 12],'Tests',"chi2",'Alpha',0.01, ...
        'Plots',["first","second","firstTwo"]);
%}
%
%{
    %% Running silently and using the output structure.
    % Display is suppressed and plots are disabled; results are accessed
    % directly from the output struct.
    rng(1)
    V = 10.^(3*rand(2000,1));
    out = benford(V,'Tests',"chi2",'Plots',false,'Display',false);
    out.chi2.firstTwo            % statistic, df, p-value and decision
    out.counts.first             % observed first-digit counts
    out.expected.first           % expected first-digit counts
%}
%
%{
    % A sample contaminated with round amounts.
    % 5% of the amounts are replaced by 50, 500 or 5000. The chi-square and
    % two-stage tests reject. The deviation plot shows which digits are
    % responsible and whether the departure is relevant in practice (MAD
    % conformity category), the heatmap pinpoints the single hot cell (first
    % digit 5, second digit 0), the two-stage plot shows the observed
    % statistics against the null distributions and the critical values, and
    % the duplication plot shows the repeated values.
    rng(1)
    n = 2000;
    V = 10.^(3*rand(n,1));
    idx = randperm(n,100);
    V(idx) = 50*10.^randi([0 2],100,1);
    out = benford(V,'Tests',["chi2","twostage"], ...
        'Plots',["deviation","heatmap","twostage","duplication"]);
%}
%
%{
    % Is the Benford hypothesis plausible at all? Order of magnitude.
    % A sample spanning less than one decade cannot follow Benford's law:
    % all the values start with 1. The magnitude plot shows the narrow range.
    rng(1)
    V = 100 + 15*rand(1000,1);
    out = benford(V,'Tests',"chi2",'Plots',["magnitude","first"]);
%}
%
%{
    %% Mantissa-based tests.
    % Include arc, mean, variance, joint, KS, U2, d*, Kuiper.
    % Benford sample against a sample uniform in (1,1000), whose mantissa
    % is not uniform. The mantissa plot shows the ECDF difference with the
    % KS band and the Kuiper strip: the curve of the uniform sample leaves
    % both. The p-values (asymptotic and Monte Carlo) are in the output.
    rng(1)
    Vb = 10.^(3*rand(1000,1));
    Vu = 1 + 999*rand(1000,1);
    outB = benford(Vb,'Tests',"mantissa",'Plots',"mantissa");
    outU = benford(Vu,'Tests',"mantissa",'Plots',"mantissa");
    outU.mantissaTests.kuiper
%}
%
%{
    % Controlling the Monte Carlo accuracy.
    % B: replications of the mantissa-based and savvy tests (if empty it is
    % chosen from MCHalfWidth, MCCL and MCPRef). BTwoStage: replications of
    % the two-stage test. PVTS: percentage used for its critical values.
    rng(1)
    V = 10.^(3*rand(1000,1));
    out = benford(V,'Tests',["twostage","mantissa"],'B',10000, ...
        'BTwoStage',50000,'PVTS',1,'Plots',false);
%}
%
%{
    %% Benford-savvy tests: only the first digit mimics Benford's law.
    % The first digit is Benford but the fractional part of the significand
    % is uniform (manipulated-Benford scheme). The first-digit Pearson
    % statistic Q1 is within its null range, whereas Q_Delta, KS2, KU2 and
    % the combined tests G_KS and G_KU typically reject. The savvy plot
    % shows the observed point against the Monte Carlo null sample.
    rng(1)
    n = 1000;
    V = floor(10.^rand(n,1)) + rand(n,1);
    out = benford(V,'Tests',["chi2","savvy"],'Plots',"savvy");
%}
%
%{
    %% Lacunary single-observation test: flagging individual observations.
    % 30 amounts are exact powers of ten, whose mantissa is exactly zero
    % (all the digits are 0): they are flagged with a tiny p-value. About
    % alpha*n further observations are flagged by chance, as the QQ plot of
    % the p-values shows (points above the band are not explained by
    % chance). LacunaryM is the number of mantissa digits used.
    rng(1)
    V = 10.^(3*rand(1000,1));
    idx = randperm(1000,30);
    V(idx) = 10.^randi([1 3],30,1);
    out = benford(V,'Tests',"lacunary",'LacunaryM',10,'Plots',"lacunary");
    all(ismember(idx,out.lacunary.idxReject))   % true: all are flagged
    V(out.lacunary.idxReject(1:5))
%}
%
%{
    % All the tests and all the plots on a clean sample.
    % The duplication diagnostics are computed because "all" includes them.
    rng(1)
    V = 10.^(3*rand(1000,1));
    out = benford(V,'Tests',["chi2","twostage","mantissa","savvy","lacunary"], ...
        'Plots',"all");
%}
%
%{
    % Redrawing selected plots from a stored output.
    % The plots can be drawn again, for example with another level for the
    % bands, with benfordplot, without repeating the tests.
    rng(1)
    V = 10.^(3*rand(1000,1));
    out = benford(V,'Tests',["chi2","twostage","mantissa"], ...
        'Plots',false,'Display',false);
    h = benfordplot(out,'plots',["deviation","mantissa"],'alpha',0.01, ...
        'docked',false);
%}
%
%{
    % Cross-checking the digit extraction with independent methods.
    % Three independent extraction routines are compared for consistency.
    rng(1)
    V = 10.^(3*rand(1000,1));
    out = benford(V,'CheckExtraction',true,'Plots',false);
    out.extractionCheck.isConsistent
%}


%% Beginning of code

% Options parsing: default options
options = struct( ...
    'Digits', [1 2 12], ...
    'Plots', ["first", "second"], ...
    'Tests', ["chi2", "twostage"], ...
    'Alpha', 0.05, ...
    'B', [], ... 
    'BTwoStage', [], ...         
    'BTwoStageTailN', 5000, ...  
    'MCHalfWidth', 0.005, ...
    'MCCL', 0.95, ...
    'MCPRef', 0.05, ...
    'PVTS', 5, ...
    'LacunaryM', 10, ...        
    'Display', true, ...
    'Overview', false, ...
    'CheckExtraction', false);


% Options parsing: varargin handling.
if ~isempty(varargin)
    [varargin{:}] = convertStringsToChars(varargin{:});
end

UserOptions = varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    if numel(varargin) ~= 2 * numel(UserOptions)
        error('FSDA:benford:WrongInputOpt', ...
            'Number of supplied options is invalid. Probably values for some parameters are missing.');
    end

    inpchk = isfield(options, UserOptions);
    WrongOptions = UserOptions(inpchk == 0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:benford:NonExistInputOpt', ...
            'In total %d non-existent user options found.', numel(WrongOptions));
    end
end

for i = 1:2:numel(varargin)
    options.(varargin{i}) = varargin{i+1};
end

% Validate scalar numeric options.
if ~isscalar(options.Alpha) || ~isfinite(options.Alpha) || options.Alpha <= 0 || options.Alpha >= 1
    error('FSDA:benford:WrongAlpha', 'Alpha must be a scalar between 0 and 1.');
end

if ~isscalar(options.LacunaryM) || ~isfinite(options.LacunaryM) || ...
        options.LacunaryM ~= round(options.LacunaryM) || ...
        options.LacunaryM < 1 || options.LacunaryM > 15
    error('FSDA:benford:WrongLacunaryM', ...
        'LacunaryM must be an integer between 1 and 15.');
end

if ~(isempty(options.B) || (isscalar(options.B) && isfinite(options.B) && options.B > 0))
    error('FSDA:benford:WrongB', 'B must be a positive scalar or empty.');
end

if ~(isempty(options.BTwoStage) || (isscalar(options.BTwoStage) && isfinite(options.BTwoStage) && options.BTwoStage > 0))
    error('FSDA:benford:WrongBTwoStage', 'BTwoStage must be a positive scalar or empty.');
end

if ~isscalar(options.BTwoStageTailN) || ~isfinite(options.BTwoStageTailN) || options.BTwoStageTailN <= 0
    error('FSDA:benford:WrongBTwoStageTailN', ...
        'BTwoStageTailN must be a positive scalar.');
end

if ~isscalar(options.MCHalfWidth) || ~isfinite(options.MCHalfWidth) || ...
        options.MCHalfWidth <= 0 || options.MCHalfWidth >= 0.5
    error('FSDA:benford:WrongMCHalfWidth', ...
        'MCHalfWidth must be a scalar in (0, 0.5).');
end

if ~isscalar(options.MCCL) || ~isfinite(options.MCCL) || options.MCCL <= 0 || options.MCCL >= 1
    error('FSDA:benford:WrongMCCL', ...
        'MCCL must be a scalar in (0, 1).');
end

if ~isscalar(options.MCPRef) || ~isfinite(options.MCPRef) || options.MCPRef < 0 || options.MCPRef > 1
    error('FSDA:benford:WrongMCPRef', ...
        'MCPRef must be a scalar in [0, 1].');
end

if ~isscalar(options.PVTS) || ~isfinite(options.PVTS) || options.PVTS <= 0 || options.PVTS >= 100
    error('FSDA:benford:WrongPVTS', 'PVTS must be a scalar in (0, 100).');
end

if isnumeric(options.Display) && isscalar(options.Display) && ismember(options.Display, [0 1])
    options.Display = logical(options.Display);
end

if ~isscalar(options.Display) || ~islogical(options.Display)
    error('FSDA:benford:WrongDisplay', 'Display must be a logical scalar.');
end

if ~isscalar(options.Overview) || ~islogical(options.Overview)
    error('FSDA:benford:WrongOverview', 'Overview must be a logical scalar.');
end

if ~isscalar(options.CheckExtraction) || ~islogical(options.CheckExtraction)
    error('FSDA:benford:WrongCheckExtraction', ...
        'CheckExtraction must be a logical scalar.');
end


% Validate Digits.
allowedDigits = [1 2 12];
if ~isnumeric(options.Digits) || ~isvector(options.Digits) || isempty(options.Digits)
    error('FSDA:benford:WrongDigits', ...
        'Digits must be a numeric vector containing values from [1 2 12].');
end
requestedDigits = unique(options.Digits(:).', 'stable');
if any(~ismember(requestedDigits, allowedDigits))
    error('FSDA:benford:WrongDigits', ...
        'Digits must contain only values from [1 2 12].');
end

useFirst    = ismember(1, requestedDigits);
useSecond   = ismember(2, requestedDigits);
useFirstTwo = ismember(12, requestedDigits);

% Validate text options.
allowedPlots = ["first", "second", "firsttwo", "deviation", "heatmap", ...
    "magnitude", "mantissa", "duplication", "twostage", "savvy", "lacunary", "all"];
allowedTests = ["chi2", "twostage", "mantissa", "savvy", "lacunary"];

% A numeric 0/1 is accepted as a logical value, for Plots and Display
if isnumeric(options.Plots) && isscalar(options.Plots) && ismember(options.Plots, [0 1])
    options.Plots = logical(options.Plots);
end

if islogical(options.Plots) && isscalar(options.Plots)
    opts.doPlots = options.Plots;
    if opts.doPlots
        opts.plots = "all";
    else
        opts.plots = string.empty(1, 0);
    end
else
    opts.plots = normalizeTextList(options.Plots);
    if any(~ismember(opts.plots, allowedPlots))
        error('FSDA:benford:WrongPlots', ...
            'Plots must contain only: %s.', strjoin(cellstr(allowedPlots), ', '));
    end
    opts.doPlots = ~isempty(opts.plots);
end

opts.tests = normalizeTextList(options.Tests);
if any(~ismember(opts.tests, allowedTests))
    error('FSDA:benford:WrongTests', ...
        'Tests must contain only "chi2", "twostage", "mantissa", "savvy", and "lacunary".');
end

% Derived internal options.
opts.alpha    = options.Alpha;

if isempty(options.B)
    opts.B = chooseMonteCarloB(options.MCHalfWidth, options.MCCL, options.MCPRef);
else
    opts.B = round(options.B);
end
if isempty(options.BTwoStage)
    opts.BTwoStage = chooseTwoStageB(options.PVTS, options.BTwoStageTailN);
else
    opts.BTwoStage = round(options.BTwoStage);
end
opts.BTwoStageTailN = options.BTwoStageTailN;

opts.PVTS     = options.PVTS;
opts.Display  = options.Display;
opts.Overview = options.Overview;
opts.Digits   = requestedDigits;

% Derived internal flags.
opts.doChi2         = any(opts.tests == "chi2");
opts.doTwoStage     = any(opts.tests == "twostage");
opts.doMantissa     = any(opts.tests == "mantissa");
opts.doSavvy        = any(opts.tests == "savvy");
opts.doLacunary     = any(opts.tests == "lacunary");
opts.LacunaryM      = options.LacunaryM;
opts.doDuplication  = opts.doPlots && any(ismember(opts.plots, ["duplication", "all"]));

% Extra digit extraction methods
opts.CheckExtraction = options.CheckExtraction;


%% Initialize the output structure.
out = struct();
out.data = x;
out.options = opts;
out.digits = struct();
out.counts = struct();
out.expected = struct();
out.figures = struct();
out.chi2 = struct();
out.twoStage = struct();
out.extractionCheck = struct( ...
    'isConsistent', true, ...
    'nMismatch', 0, ...
    'idxMismatch', [], ...
    'xMismatch', [], ...
    'table', table());

%% Digits extraction

% Sample size after cleaning.
nObs = numel(x);
out.n = nObs;

% Extract the digits with uo to three independent methods.
if opts.doSavvy
    [D1, D2, D12, Sig] = ExtractDigits(x);    % significand needed by the savvy tests
else
    [D1, D2, D12] = ExtractDigits(x);
    Sig = [];
end

if opts.CheckExtraction
    [D1ac, D2ac, D12ac] = ExtractDigitsAC(x);
    [D1s, D2s, D12s]    = ExtractDigitsStr(x);

    digitsMain = [D1(:), D2(:), D12(:)];
    digitsAC   = [D1ac(:), D2ac(:), D12ac(:)];
    digitsStr  = [D1s(:), D2s(:), D12s(:)];

    mismatchRows = any(digitsMain ~= digitsAC, 2) | any(digitsMain ~= digitsStr, 2);

    out.extractionCheck = struct( ...
        'isConsistent', ~any(mismatchRows), ...
        'nMismatch', nnz(mismatchRows), ...
        'idxMismatch', find(mismatchRows), ...
        'xMismatch', x(mismatchRows), ...
        'table', table());

    if any(mismatchRows)
        out.extractionCheck.table = table( ...
            out.extractionCheck.idxMismatch, ...
            out.extractionCheck.xMismatch, ...
            digitsMain(mismatchRows, 1), digitsMain(mismatchRows, 2), digitsMain(mismatchRows, 3), ...
            digitsAC(mismatchRows, 1), digitsAC(mismatchRows, 2), digitsAC(mismatchRows, 3), ...
            digitsStr(mismatchRows, 1), digitsStr(mismatchRows, 2), digitsStr(mismatchRows, 3), ...
            'VariableNames', { ...
                'idx', 'x', ...
                'D1', 'D2', 'D12', ...
                'D1ac', 'D2ac', 'D12ac', ...
                'D1str', 'D2str', 'D12str'});
    end
end

% ---------------------------------------------------------------------
% Observed and theoretical digit distributions
% ---------------------------------------------------------------------

% Main counting path. Compute the Benford expected frequencies.
if useFirst
    out.counts.first = histcounts(D1, 0.5:1:9.5).';
    out.expected.first = benfordProbabilities(1) * nObs;
end
if useSecond
    out.counts.second = histcounts(D2, -0.5:1:9.5).';
    out.expected.second = benfordProbabilities(2) * nObs;
end
if useFirstTwo
    out.counts.firstTwo = histcounts(D12, 9.5:1:99.5).';
    out.expected.firstTwo = benfordProbabilities(12) * nObs;
end

% Counting of potential mismatches n digit extraction
out.countsCheck = struct( ...
    'isConsistent', true, ...
    'firstMismatch', false, ...
    'secondMismatch', false, ...
    'firstTwoMismatch', false);

if opts.CheckExtraction
    cFirstAC    = countDigits(D1ac, 1:9);
    cSecondAC   = countDigits(D2ac, 0:9);
    cFirstTwoAC = countDigits(D12ac, 10:99);

    mismatchFirst    = useFirst    && any(out.counts.first    ~= cFirstAC);
    mismatchSecond   = useSecond   && any(out.counts.second   ~= cSecondAC);
    mismatchFirstTwo = useFirstTwo && any(out.counts.firstTwo ~= cFirstTwoAC);

    out.countsCheck.isConsistent    = ~(mismatchFirst || mismatchSecond || mismatchFirstTwo);
    out.countsCheck.firstMismatch   = mismatchFirst;
    out.countsCheck.secondMismatch  = mismatchSecond;
    out.countsCheck.firstTwoMismatch = mismatchFirstTwo;
end


% Mantissa extraction.
out.mantissa      = ExtractMantissa(x);
out.mantissaTests = struct();

if opts.doMantissa
    out.mantissaTests = mantissa_tests(x, opts.B);

    if opts.Display
        fprintf('\nMantissa tests\n');
        fprintf('%-12s %12s %12s %12s\n', 'Test', 'Statistic', 'p_asym', 'p_mc');

        mantissaFields = {'arc', 'mean', 'variance', 'joint', 'ks', 'u2', 'dstar', 'kuiper'};
        mantissaLabels = {'Arc', 'Mean', 'Variance', 'Joint', 'KS', 'U2', 'd*', 'Kuiper'};

        for k = 1:numel(mantissaFields)
            S = out.mantissaTests.(mantissaFields{k});
            fprintf('%-12s %12.4f %12.4g %12.4g\n', ...
                mantissaLabels{k}, S.stat, S.pAsym, S.pMC);
        end
    end
end


% ---------------------------------------------------------------------
% Benford-savvy tests (Barabesi, Cerioli, Cerasa and Perrotta, 2026)
% ---------------------------------------------------------------------

out.savvy = struct();

if opts.doSavvy
    out.savvy = savvy_tests(D1, Sig, opts.B);

    if opts.Display
        fprintf('\nBenford-savvy tests (n = %d, B = %d)\n', out.savvy.n, out.savvy.B);
        fprintf('%-12s %12s %12s %12s\n', 'Test', 'Statistic', 'p_asym', 'p_mc');

        savvyFields = {'Q1', 'Q2', 'KS2', 'KU2', 'QDelta', 'GKS', 'GKU'};
        savvyLabels = {'Q1', 'Q2', 'KS2', 'KU2', 'Q_Delta', 'G_KS', 'G_KU'};

        for k = 1:numel(savvyFields)
            S = out.savvy.(savvyFields{k});
            fprintf('%-12s %12.4f %12.4g %12.4g\n', ...
                savvyLabels{k}, S.stat, S.pAsym, S.pMC);
        end
    end
end

% ---------------------------------------------------------------------
% Lacunary single-observation test (Ferrari, 2026)
% ---------------------------------------------------------------------

out.lacunary = struct();

if opts.doLacunary
    out.lacunary = lacunary_tests(out.mantissa, opts.LacunaryM, opts.alpha);

    if opts.Display
        L = out.lacunary;
        fprintf('\nLacunary single-observation test (m = %d digits, chi2 with %d df)\n', ...
            L.m, L.df);
        fprintf('Observations flagged at alpha = %g: %d of %d (about %.1f expected under H0)\n', ...
            L.alpha, L.nReject, numel(L.pValue), L.nExpected);
        fprintf('Binomial p-value for the number of flagged observations: %.4g\n', L.pBinom);
    end
end

% ---------------------------------------------------------------------
% Plots: drawn at the end of the function, by benfordplot, once all the
% requested tests have been run
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% Chi-square tests
% ---------------------------------------------------------------------

% Compute the chi-square test only for the requested hypothesis set.
if opts.doChi2
    if useFirst
        % chi2
        out.chi2.first.stat   = Chi2onDigitDistrib(D1, 1);
        % Degrees of freedom.
        out.chi2.first.df     = 8;
        % Convert the chi-square statistic into p-values.
        out.chi2.first.pValue = 1 - chi2cdf(out.chi2.first.stat, out.chi2.first.df);
        % Apply the significance-level decision rule.
        out.chi2.first.reject = out.chi2.first.pValue < opts.alpha;
        % Optional command-window summary.
        if opts.Display
            fprintf('Chi2 first digit  = %.4f, p = %.4g, reject = %d\n', ...
                out.chi2.first.stat, out.chi2.first.pValue, out.chi2.first.reject);
        end
    end

    if useSecond
        out.chi2.second.stat   = Chi2onDigitDistrib(D2, 2);
        out.chi2.second.df   = 9;
        out.chi2.second.pValue   = 1 - chi2cdf(out.chi2.second.stat, out.chi2.second.df);
        out.chi2.second.reject   = out.chi2.second.pValue < opts.alpha;
        if opts.Display
            fprintf('Chi2 second digit = %.4f, p = %.4g, reject = %d\n', ...
                out.chi2.second.stat, out.chi2.second.pValue, out.chi2.second.reject);
        end
    end

    if useFirstTwo
        out.chi2.firstTwo.stat = Chi2onDigitDistrib(D12, 12);
        out.chi2.firstTwo.df = 89;
        out.chi2.firstTwo.pValue = 1 - chi2cdf(out.chi2.firstTwo.stat, out.chi2.firstTwo.df);
        out.chi2.firstTwo.reject = out.chi2.firstTwo.pValue < opts.alpha;
        if opts.Display
            fprintf('Chi2 first 2 digs = %.4f, p = %.4g, reject = %d\n', ...
                out.chi2.firstTwo.stat, out.chi2.firstTwo.pValue, out.chi2.firstTwo.reject);
        end
    end
end

% ---------------------------------------------------------------------
% Two-stage tests
% ---------------------------------------------------------------------

% Compute critical values and observed statistics for the two-stage test.
if opts.doTwoStage
    if ~all([useFirst, useSecond, useFirstTwo])
        error('FSDA:benford:DigitsRequiredForTwostage', ...
            'Two-stage test requires Digits = [1 2 12].');
    end
    nD = 2;
    [labels, cv_mad, cv_chi2, madSim, chi2Sim] = ...
        TwoStageTestCV(nObs, opts.BTwoStage, nD, opts.PVTS);

    % Store observed test statistics for the three nested hypotheses.
    obs_mad = zeros(3, 1);
    obs_chi2 = zeros(3, 1);

    % Mean absolute deviation statistics.
    obs_mad(1) = sum(abs(out.counts.first - out.expected.first)) / ...
        (nObs * numel(out.expected.first));
    obs_mad(2) = sum(abs(out.counts.second - out.expected.second)) / ...
        (nObs * numel(out.expected.second));
    obs_mad(3) = sum(abs(out.counts.firstTwo - out.expected.firstTwo)) / ...
        (nObs * numel(out.expected.firstTwo));

    % Chi-square-type statistics.
    obs_chi2(1) = sum(((out.counts.first - out.expected.first) .^ 2) ./ out.expected.first);
    obs_chi2(2) = sum(((out.counts.second - out.expected.second) .^ 2) ./ out.expected.second);
    obs_chi2(3) = sum(((out.counts.firstTwo - out.expected.firstTwo) .^ 2) ./ out.expected.firstTwo);

    % % Empirical Monte Carlo p-values.
    %{ 
        % Approccio di Andrea
        outAC = struct;
        outAC.twoStage.pValueMad = zeros(3, 1);
        outAC.twoStage.pValueChi2 = zeros(3, 1);
    
        for j = 1:3
            outAC.twoStage.pValueMad(j)  = (sum(madSim(j, :)  >= obs_mad(j))  + 1) / (size(madSim, 2)  + 1);
            outAC.twoStage.pValueChi2(j) = (sum(chi2Sim(j, :) >= obs_chi2(j)) + 1) / (size(chi2Sim, 2) + 1);
        end
    %}

    
    % Andrea's version optimised for speed
    % (the line "out = struct;" that used to be here discarded every field
    % computed before the two-stage test, and has been removed)
    out.twoStage.pValueMad  = mc_pvalue(madSim.', obs_mad);
    out.twoStage.pValueChi2 = mc_pvalue(chi2Sim.', obs_chi2);
    
    out.twoStage.labels = labels;
    out.twoStage.madSim  = madSim.';     % B-by-3 simulated MAD under H0 (used by benfordplot)
    out.twoStage.chi2Sim = chi2Sim.';    % B-by-3 simulated chi2 statistics under H0
    out.twoStage.cvMad = cv_mad(:);
    out.twoStage.cvChi2 = cv_chi2(:);
    out.twoStage.obsMad = obs_mad;
    out.twoStage.obsChi2 = obs_chi2;
    out.twoStage.rejectMad = obs_mad >= cv_mad(:);
    out.twoStage.rejectChi2 = obs_chi2 >= cv_chi2(:);

    % IMPORTANT REMARK: TO CHECK !! 
    % The critical values for H1 and H2 are conditional on H12 being
    % rejected. rejectMad/rejectChi2 for H1 and H2 were reported even when
    % H12 is not rejected. Also, pValueMad(1:2) were unconditional Monte
    % Carlo p-values, while the critical values were conditional, so the
    % two are not coherent. We should report reject(1:2) = reject(3) & (obs
    % >= cv) and compute the stage-2 p-values on the same conditioning
    % subset.

    %{
        outDP.twoStage.labels   = labels;
        outDP.twoStage.cvMad    = cv_mad(:);
        outDP.twoStage.cvChi2   = cv_chi2(:);
        outDP.twoStage.obsMad   = obs_mad;
        outDP.twoStage.obsChi2  = obs_chi2;
        %
        outDP.twoStage.rejectMad  = false(3, 1);
        outDP.twoStage.rejectChi2 = false(3, 1);
        %
        outDP.twoStage.rejectMad(3)    = obs_mad(3)  >= cv_mad(3);
        outDP.twoStage.rejectChi2(3)   = obs_chi2(3) >= cv_chi2(3);
        %
        outDP.twoStage.rejectMad(1:2)  = out.twoStage.rejectMad(3)  & (obs_mad(1:2)  >= cv_mad(1:2));
        outDP.twoStage.rejectChi2(1:2) = out.twoStage.rejectChi2(3) & (obs_chi2(1:2) >= cv_chi2(1:2));
        %
        outDP.twoStage.pValueMad  = nan(3, 1);
        outDP.twoStage.pValueChi2 = nan(3, 1);
        %
        outDP.twoStage.pValueMad(3)  = mc_pvalue(madSim(3, :).', obs_mad(3));
        outDP.twoStage.pValueChi2(3) = mc_pvalue(chi2Sim(3, :).', obs_chi2(3));
        
        maskMad  = madSim(3, :)  >= cv_mad(3);
        maskChi2 = chi2Sim(3, :) >= cv_chi2(3);
        
        if any(maskMad)
            outDP.twoStage.pValueMad(1) = mc_pvalue(madSim(1, maskMad).', obs_mad(1));
            outDP.twoStage.pValueMad(2) = mc_pvalue(madSim(2, maskMad).', obs_mad(2));
        end
        
        if any(maskChi2)
            outDP.twoStage.pValueChi2(1) = mc_pvalue(chi2Sim(1, maskChi2).', obs_chi2(1));
            outDP.twoStage.pValueChi2(2) = mc_pvalue(chi2Sim(2, maskChi2).', obs_chi2(2));
        end
    %}


    % Optional command-window summary.
    if opts.Display
        fprintf('\nTwo-stage test (MAD)\n');
        for k = 1:numel(labels)
            fprintf('H%s: stat = %.4f, cv = %.4f, p = %.4g, reject = %d\n', ...
                labels{k}, obs_mad(k), cv_mad(k), out.twoStage.pValueMad(k), ...
                out.twoStage.rejectMad(k));
        end

        fprintf('\nTwo-stage test (chi2)\n');
        for k = 1:numel(labels)
            fprintf('H%s: stat = %.4f, cv = %.4f, p = %.4g, reject = %d\n', ...
                labels{k}, obs_chi2(k), cv_chi2(k), out.twoStage.pValueChi2(k), ...
                out.twoStage.rejectChi2(k));
        end
    end

end

% ---------------------------------------------------------------------
% Plots
% ---------------------------------------------------------------------

if opts.doPlots
    % The duplication diagnostics are computed only if they are going to be plotted
    if opts.doDuplication
        out.duplication = duplication_tests(x);
    end

    % Digit-group plots of groups which were not analysed are dropped silently
    plotsToDraw = opts.plots;
    if ~useFirst
        plotsToDraw(plotsToDraw == "first") = [];
    end
    if ~useSecond
        plotsToDraw(plotsToDraw == "second") = [];
    end
    if ~useFirstTwo
        plotsToDraw(plotsToDraw == "firsttwo") = [];
    end

    if ~isempty(plotsToDraw)
        out.figures = benfordplot(out, 'plots', plotsToDraw, 'alpha', opts.alpha);
    end
end

end

% -------------------------------------------------------------------

%% Internal Dependencies: digit extraction 

function [D1, D2, D12, S] = ExtractDigits(V)
%
%   [D1,D2,D12,S] = ExtractDigits(V) returns:
%       D1   first digit of each observation
%       D2   second digit of each observation
%       D12  first two digits of each observation
%       S    (optional) significand in [1,10), consistent with D1, D2, D12
%
%   The input V is assumed to contain positive numeric values. Each value
%   is scaled so that its first two significant digits lie in the interval
%   [10,99]. The first digit is then obtained from D12 and the second digit
%   from the remainder.
%
%   The input V must contain positive finite values.

sz = size(V);
V = V(:);

if any(~isfinite(V)) || any(V <= 0)
    error('FSDA:benford:InvalidData', ...
        'ExtractDigits requires positive finite values.');
end

% Scale each value so that its first two significant digits lie in [10,99].
expo   = floor(log10(V));
scaled = V ./ 10.^expo;

scaled(scaled < 1)   = scaled(scaled < 1)  * 10;
scaled(scaled >= 10) = scaled(scaled >= 10) / 10;
tol = 1e-12;
if isa(V,'single'), tol = 1e-6; end
D12 = uint8(min(floor(10 .* scaled .* (1 + tol)), 99));

% Split the first two digits into first and second digit.
% D1 = floor(D12 ./ 10);
% D2 = D12 - 10 .* D1;
D12d = double(D12); 
D1   = floor(D12d./10); 
D2   = D12d - 10*D1;

D1  = reshape(D1, sz);
D2  = reshape(D2, sz);
D12 = reshape(D12, sz);

if any(D12(:) < 10 | D12(:) > 99) || any(D1(:) < 1 | D1(:) > 9) || ...
   any(D2(:) < 0 | D2(:) > 9) || any(10*double(D1(:)) + double(D2(:)) ~= double(D12(:)))
    error('FSDA:benford:InvalidDigits', ...
        'ExtractDigits produced out-of-range or inconsistent digits.');
end

% Significand, consistent with the digits above (same tolerance, and kept
% strictly below 10 as D12 is clamped at 99).
if nargout > 3
    S = scaled .* (1 + tol);
    S = min(S, 10 * (1 - eps(class(S))));
    S = reshape(S, sz);
end



end


function [D1, D2, D12] = ExtractDigitsAC(V)
%ExtractDigits Extract the first digit, second digit, and first two digits
%
%   [D1,D2,D12] = ExtractDigits(V) returns:
%       D1   first digit of each observation
%       D2   second digit of each observation
%       D12  first two digits of each observation
%
%   The input V is assumed to contain positive numeric values. Each value
%   is scaled so that its first two significant digits lie in the interval
%   [10,99]. The first digit is then obtained from D12 and the second digit
%   from the remainder.

[n, nsim] = size(V);
D12 = zeros(n, nsim);

for i = 1:n
    for j = 1:nsim
        V_ij = V(i, j);

        % Scale the value until its first significant digits are visible.
        while V_ij < 10
            V_ij = V_ij * 10;
        end

        % Reduce the value until only the first two digits remain.
        while V_ij >= 100
            V_ij = V_ij / 10;
        end

        % Store the first two digits.
        D12(i, j) = floor(V_ij);
    end
end

% Split the first two digits into first and second digit.
D1 = floor(D12 / 10);
D2 = D12 - 10 * floor(D12 / 10);
end

% -------------------------------------------------------------------------

function [D1, D2, D12] = ExtractDigitsStr(V)
%ExtractDigitsStr Extract digits using scientific-notation strings
%
%   [D1,D2,D12] = ExtractDigitsStr(V) returns the first digit, second
%   digit, and first two digits using a string-based extraction.

sz = size(V);
V = V(:);

if any(~isfinite(V)) || any(V <= 0)
    error('FSDA:benford:InvalidData', ...
        'ExtractDigitsStr requires positive finite values.');
end

if isa(V, 'single')
    fmt = '%.5e';
else
    fmt = '%.14e';
end

s = compose(fmt, V);
mant = extractBefore(s, "e");   % e.g. "1.234567890123456e+03"

D1 = str2double(extractBetween(mant, 1, 1));
D2 = str2double(extractBetween(mant, 3, 3));
D12 = 10 .* D1 + D2;

D1  = reshape(D1, sz);
D2  = reshape(D2, sz);
D12 = reshape(D12, sz);
end


% -------------------------------------------------------------------------

function m = ExtractMantissa(V)
%ExtractMantissa Compute the base-10 mantissa in [0,1)

sz = size(V);
V = V(:);
m = mod(log10(V), 1);
m = reshape(m, sz);
end


%% Internal Dependencies: statistical distributions and tests

function fben = benfordProbabilities(nd)
%benfordProbabilities Benford probabilities for standard digit groups
%
%   fben = benfordProbabilities(nd) returns the Benford probability vector
%   for a standard digit group:
%       nd = 1   first digit (1 to 9)
%       nd = 2   second digit (0 to 9)
%       nd = 12  first two digits (10 to 99)
%
%   This function is used to build the expected frequencies for the main
%   digit-distribution tests.

log10base = log(10);

if nd == 1
    % First digit probabilities: P(D1 = d) = log10(1 + 1/d), d = 1,...,9
    d = (1:9)';
    fben = log(1 + 1 ./ d) ./ log10base;

elseif nd == 2
    % Second digit probabilities are obtained by aggregating over 10..99.
    d = 10:99;
    fben12 = log(1 + 1 ./ d) ./ log10base;
    fben12 = reshape(fben12, [10 9]);
    fben = sum(fben12, 2);

elseif nd == 12
    % First-two-digit probabilities: P(D12 = d) = log10(1 + 1/d)
    d = (10:99)';
    fben = log(1 + 1 ./ d) ./ log10base;

else
    % Unsupported digit level.
    fben = [];
end
end

% -------------------------------------------------------------------------

function [D, Theo] = benfordPatternSupport(nD)
%benfordPatternSupport Theoretical support for an nD-digit Benford pattern
%
%   [D,Theo] = benfordPatternSupport(nD) builds the full set of digit patterns
%   of length nD and their corresponding Benford probabilities.
%
%   D contains all possible digit combinations in lexicographic order.
%   Theo contains the Benford probability associated with each pattern.
%
%   This is an internal helper used by the two-stage test machinery to
%   work with the full pattern space and its marginal sub-hypotheses.

first = 10^(nD - 1);
last = 10^nD - 1;
nn = last - first + 1;

D = zeros(nn, nD);
Theo = zeros(nn, 1);

for i = first:last
    den = 0;
    itext = num2str(i);

    % Store each digit of the pattern and rebuild its numeric value.
    for j = 1:nD
        d_ij = str2num(itext(j)); %#ok<ST2NM>
        D(i - first + 1, j) = d_ij;
        den = den + 10^(nD - j) * d_ij;
    end

    % Benford probability for the pattern.
    Theo(i - first + 1) = log(1 + 1 / den) / log(10);
end
end

% -------------------------------------------------------------------------

function chi2 = Chi2onDigitDistrib(D, nd)
%Chi2onDigitDistrib Chi-square statistic for digit frequencies
%
%   chi2 = Chi2onDigitDistrib(D,nd) computes the chi-square statistic that
%   compares the observed digit distribution in D with the Benford
%   theoretical distribution.
%
%   Input:
%       D   Matrix of observed digits, one column per sample.
%       nd  Either:
%             - a scalar: 1, 2, or 12, indicating the digit group
%             - a vector: the theoretical Benford probabilities
%
%   Output:
%       chi2  Chi-square statistic for each column of D.
%
%   The function first builds the observed frequency table and then
%   compares it with the expected Benford frequencies.

[n, nsim] = size(D);

if isscalar(nd)
    fben = benfordProbabilities(nd);

    switch nd
        case 1
            idx = D;                 % digits 1..9 -> categories 1..9
            nCat = 9;
        case 2
            idx = D + 1;             % digits 0..9 -> categories 1..10
            nCat = 10;
        case 12
            idx = D - 9;             % digits 10..99 -> categories 1..90
            nCat = 90;
        otherwise
            error('FSDA:benford:WrongDigits', ...
                'nd must be 1, 2, or 12.');
    end

else
    fben = nd(:);

    switch numel(fben)
        case 9
            idx = D;
            nCat = 9;
        case 10
            idx = D + 1;
            nCat = 10;
        case 90
            idx = D - 9;
            nCat = 90;
        otherwise
            error('FSDA:benford:WrongDigits', ...
                'Probability vector must have length 9, 10, or 90.');
    end
end

% Count observations per category and simulation.
colIdx = repmat(1:nsim, n, 1);
f = accumarray([idx(:), colIdx(:)], 1, [nCat, nsim], @sum, 0);

% Expected frequencies.
fben = n * fben(:);
fben = repmat(fben, 1, nsim);

% Chi-square statistic for each simulation column.
chi2 = sum(((f - fben) .^ 2) ./ fben, 1);
end

% function chi2 = Chi2onDigitDistrib(D, nd)
% %Chi2onDigitDistrib Chi-square statistic for digit frequencies
% %
% %   chi2 = Chi2onDigitDistrib(D,nd) computes the chi-square statistic that
% %   compares the observed digit distribution in D with the Benford
% %   theoretical distribution.
% %
% %   Input:
% %       D   Matrix of observed digits, one column per sample.
% %       nd  Either:
% %             - a scalar: 1, 2, or 12, indicating the digit group
% %             - a vector: the theoretical Benford probabilities
% %
% %   Output:
% %       chi2  Chi-square statistic for each column of D.
% %
% %   The function first builds the observed frequency table and then
% %   compares it with the expected Benford frequencies.
% 
% [n, nsim] = size(D);
% 
% if isscalar(nd)
%     % Build the observed frequencies from the digit values.
%     if nd == 1
%         % First digit: categories 1,...,9
%         f = zeros(9, nsim);
%         for d = 1:9
%             f(d, :) = sum(D == d);
%         end
% 
%     elseif nd == 2
%         % Second digit: categories 0,...,9
%         f = zeros(10, nsim);
%         for d = 0:9
%             f(d + 1, :) = sum(D == d);
%         end
% 
%     elseif nd == 12
%         % First two digits: categories 10,...,99
%         f = zeros(90, nsim);
%         for d = 10:99
%             f(d - 9, :) = sum(D == d);
%         end
%     end
% 
%     % Theoretical Benford probabilities for the selected digit group.
%     fben = benfordProbabilities(nd);
% 
% else
%     % If nd is already a probability vector, use it directly.
%     if length(nd) == 9
%         f = zeros(9, nsim);
%         for d = 1:9
%             f(d, :) = sum(D == d);
%         end
%     elseif length(nd) == 10
%         f = zeros(10, nsim);
%         for d = 0:9
%             f(d + 1, :) = sum(D == d);
%         end
%     elseif length(nd) == 90
%         f = zeros(90, nsim);
%         for d = 10:99
%             f(d - 9, :) = sum(D == d);
%         end
%     end
% 
%     fben = nd;
% end
% 
% % Convert probabilities into expected frequencies.
% fben = n * fben * ones(1, nsim);
% 
% % Chi-square contribution summed over categories.
% chi2 = sum(((f - fben) .^ 2) ./ fben);
% end


% -------------------------------------------------------------------------

function Marg = MarginalDistrib(D, Freq, cols)
%MarginalDistrib Marginal frequencies for a subset of digit positions
%
%   Marg = MarginalDistrib(D,Freq,cols) aggregates the frequency matrix
%   Freq over the subset of digit positions specified by cols.
%
%   Inputs:
%       D     Matrix of digit patterns.
%       Freq  Frequency matrix, one column per simulation.
%       cols  Logical vector selecting the digit positions to keep.
%
%   Output:
%       Marg  Marginal frequencies for the selected digit subset.
%
%   If all positions are selected, the function returns Freq unchanged.

% If all digits are selected, no marginalization is needed.
if sum(cols) == size(D, 2)
    Marg = Freq;
else
    % Keep only the selected digit positions.
    d = D(:, cols == 1);
    [n, k] = size(d);

    % Convert each digit pattern into a string label for grouping.
    dig = cell(n, 1);
    for i = 1:n
        dig_i = '';
        for j = 1:k
            dig_i = [dig_i num2str(d(i, j))]; %#ok<AGROW>
        end
        dig{i} = dig_i;
    end

    % Identify unique patterns and sum the corresponding frequencies.
    Udig = unique(dig);
    Marg = zeros(length(Udig), size(Freq, 2));
    for i = 1:length(Udig)
        rows = strcmp(dig, Udig{i});
        Marg(i, :) = sum(Freq(rows == 1, :));
    end
end
end

% -------------------------------------------------------------------------

function [labels, cols] = SubHypotheses(nD)
%SubHypotheses Build all non-empty subsets of digit positions
%
%   [labels,cols] = SubHypotheses(nD) returns the set of all non-empty
%   digit-position subsets for a pattern of length nD.
%
%   Outputs:
%       labels  Cell array with the subset labels, for example '1', '2',
%               '12', etc.
%       cols    Binary matrix whose rows identify the selected digit
%               positions for each subset.
%
%   The number of subsets is 2^nD - 1.

n = 2^nD - 1;
labels = cell(n, 1);
cols = zeros(n, nD);
row = 1;

% Enumerate all subsets of size 1 to nD.
for i = 1:nD
    c = combnk(1:nD, i);
    for j = 1:size(c, 1)
        cols(row, c(j, :)) = 1;

        % Build the label by concatenating selected digit indices.
        lab = '';
        for k = 1:i
            lab = [lab num2str(c(j, k))]; %#ok<AGROW>
        end
        labels{row} = lab;
        row = row + 1;
    end
end
end

% -------------------------------------------------------------------------

function [labels, cv_mad, cv_chi2, mad, chi2] = TwoStageTestCV(n, B, nD, pv)
%TwoStageTestCV Monte Carlo critical values for the two-stage Benford test
%
%   Computes critical values for the two-stage Benford test by Monte Carlo
%   simulation.
%
%   Inputs:
%       n   Sample size.
%       B   Number of Monte Carlo replications.
%       nD  Number of digits in the Benford pattern.
%       pv  Percentile used to derive the critical values.
%
%   Outputs:
%       labels   Labels of the sub-hypotheses.
%       cv_mad   Critical values for the MAD statistic.
%       cv_chi2  Critical values for the chi-square statistic.
%       mad      Simulated MAD statistics under the null hypothesis.
%       chi2     Simulated chi-square-type statistics under the null
%                hypothesis.
%
%   The test is based on several nested sub-hypotheses built from the
%   selected digit positions. For each replication, the function simulates
%   multinomial counts under the Benford model and computes the test
%   statistics for each sub-hypothesis.
%

% Build the full theoretical digit space and Benford probabilities.
[D, Theo] = benfordPatternSupport(nD);

% Build the nested sub-hypotheses.
[labels, cols] = SubHypotheses(nD);
nSub = numel(labels);

% Precompute the theoretical marginal for each sub-hypothesis once.
MargThCell = cell(nSub, 1);
nMarg = zeros(nSub, 1);
for i = 1:nSub
    MargThCell{i} = MarginalDistrib(D, n * Theo, cols(i, :));
    nMarg(i) = size(MargThCell{i}, 1);
end

% Preallocate outputs.
cv_mad = zeros(3, 1);
cv_chi2 = zeros(3, 1);
mad = zeros(nSub, B);
chi2 = zeros(nSub, B);

if numel(Theo) < 9000
    % Simulate all replications at once.
    RndF = mnrnd(n, Theo, B).';

    for i = 1:nSub
        c_i = cols(i, :);
        MargTh = MargThCell{i};
        MargSim = MarginalDistrib(D, RndF, c_i);

        MargTh = repmat(MargTh, 1, B);
        AbsDiff = abs(MargSim - MargTh);

        mad(i, :) = sum(AbsDiff, 1) / (n * nMarg(i));
        chi2(i, :) = sum((AbsDiff .^ 2) ./ MargTh, 1);
    end
else
    % Simulate in blocks to limit memory usage.
    blockSize = 100;
    nBlocks = ceil(B / blockSize);
    startIdx = 1;

    for block = 1:nBlocks
        thisB = min(blockSize, B - startIdx + 1);
        RndF = mnrnd(n, Theo, thisB).';

        for i = 1:nSub
            c_i = cols(i, :);
            MargTh = MargThCell{i};
            MargSim = MarginalDistrib(D, RndF, c_i);

            MargTh = repmat(MargTh, 1, thisB);
            AbsDiff = abs(MargSim - MargTh);

            idx = startIdx:(startIdx + thisB - 1);
            mad(i, idx) = sum(AbsDiff, 1) / (n * nMarg(i));
            chi2(i, idx) = sum((AbsDiff .^ 2) ./ MargTh, 1);
        end

        startIdx = startIdx + thisB;
    end
end

% Sequential critical values for the three nested hypotheses.
pct = 100 - pv;

cv_mad(3) = prctile(mad(3, :), pct);
maskMad = mad(3, :) >= cv_mad(3);
if any(maskMad)
    cv_mad(1) = prctile(mad(1, maskMad), pct);
    cv_mad(2) = prctile(mad(2, maskMad), pct);
else
    cv_mad(1) = prctile(mad(1, :), pct);
    cv_mad(2) = prctile(mad(2, :), pct);
end

cv_chi2(3) = prctile(chi2(3, :), pct);
maskChi2 = chi2(3, :) >= cv_chi2(3);
if any(maskChi2)
    cv_chi2(1) = prctile(chi2(1, maskChi2), pct);
    cv_chi2(2) = prctile(chi2(2, maskChi2), pct);
else
    cv_chi2(1) = prctile(chi2(1, :), pct);
    cv_chi2(2) = prctile(chi2(2, :), pct);
end
end


% function [labels, cv_mad, cv_chi2, mad, chi2] = TwoStageTestCV(n, B, nD, pv)
% %TwoStageTestCV Monte Carlo critical values for the two-stage Benford test
% %
% %   [labels,cv_mad,cv_chi2,mad,chi2] = TwoStageTestCV(n,B,nD,pv) computes
% %   critical values for the two-stage Benford test by Monte Carlo simulation.
% %
% %   Inputs:
% %       n   Sample size.
% %       B   Number of Monte Carlo replications.
% %       nD  Number of digits in the Benford pattern.
% %       pv  Percentile used to derive the critical values.
% %
% %   Outputs:
% %       labels   Labels of the sub-hypotheses.
% %       cv_mad   Critical values for the MAD statistic.
% %       cv_chi2  Critical values for the chi-square statistic.
% %       mad      Simulated MAD statistics under the null hypothesis.
% %       chi2     Simulated chi-square-type statistics under the null
% %                hypothesis.
% %
% %   The test is based on several nested sub-hypotheses built from the
% %   selected digit positions. For each replication, the function simulates
% %   multinomial counts under the Benford model and computes the test
% %   statistics for each sub-hypothesis.
% 
% 
% % Preallocate outputs for the three sub-hypotheses:
% % first digit, second digit, and first two digits.
% cv_mad = zeros(3, 1);
% cv_chi2 = zeros(3, 1);
% 
% % Build the full theoretical digit space and Benford probabilities.
% [D, Theo] = benfordPatternSupport(nD);
% 
% % Define the sub-hypotheses and the corresponding digit masks.
% [labels, cols] = SubHypotheses(nD);
% 
% % If the theoretical support is small, simulate all replications at once.
% if length(Theo) < 9000
%     RndF = mnrnd(n, Theo, B)';
% 
%     chi2 = zeros(length(labels), B);
%     mad = zeros(length(labels), B);
% 
%     for i = 1:length(labels)
%         c_i = cols(i, :);
% 
%         % Theoretical and simulated marginal distributions for this subset.
%         MargTh = MarginalDistrib(D, n * Theo, c_i);
%         MargSim = MarginalDistrib(D, RndF, c_i);
% 
%         % Replicate the theoretical marginal across all simulations.
%         MargTh = MargTh * ones(1, B);
% 
%         % Absolute differences between observed and theoretical margins.
%         AbsDiff = abs(MargSim - MargTh);
% 
%         % Mean absolute deviation statistic.
%         mad(i, :) = sum(AbsDiff) / (n * size(MargTh, 1));
% 
%         % Chi-square-type statistic.
%         chi2(i, :) = sum((AbsDiff .^ 2) ./ MargTh);
%     end
% else
%     % For larger supports, simulate in blocks to reduce memory usage.
%     blockSize = 100;
%     nBlocks = ceil(B / blockSize);
% 
%     mad = zeros(length(labels), B);
%     chi2 = zeros(length(labels), B);
% 
%     startIdx = 1;
%     for block = 1:nBlocks
%         t1 = tic;
%         thisB = min(blockSize, B - startIdx + 1);
%         RndF = mnrnd(n, Theo, thisB)';
% 
%         for i = 1:length(labels)
%             c_i = cols(i, :);
% 
%             MargTh = MarginalDistrib(D, n * Theo, c_i);
%             MargSim = MarginalDistrib(D, RndF, c_i);
%             MargTh = MargTh * ones(1, thisB);
% 
%             AbsDiff = abs(MargSim - MargTh);
% 
%             mad_b = sum(AbsDiff) / (n * size(MargTh, 1));
%             chi2_b = sum((AbsDiff .^ 2) ./ MargTh);
% 
%             idx = startIdx:(startIdx + thisB - 1);
%             mad(i, idx) = mad_b;
%             chi2(i, idx) = chi2_b;
%         end
% 
%         startIdx = startIdx + thisB;
% 
%         % Optional progress output for long simulations.
%         t2 = toc(t1);
%         fprintf('Block %d/%d completed in %.2f s\n', block, nBlocks, t2);
%     end
% end
% 
% % The two-stage rule works sequentially:
% % 1) compute the threshold for the most general hypothesis,
% % 2) condition the thresholds for the lower-level hypotheses on it.
% cv_mad(3) = prctile(mad(3, :), 100 - pv);
% cols_mad = (mad(3, :) >= cv_mad(3));
% cv_mad(1) = prctile(mad(1, cols_mad), 100 - pv);
% cv_mad(2) = prctile(mad(2, cols_mad), 100 - pv);
% 
% cv_chi2(3) = prctile(chi2(3, :), 100 - pv);
% cols_chi2 = (chi2(3, :) >= cv_chi2(3));
% cv_chi2(1) = prctile(chi2(1, cols_chi2), 100 - pv);
% cv_chi2(2) = prctile(chi2(2, cols_chi2), 100 - pv);
% end

% -------------------------------------------------------------------------

function S = mantissa_tests(x, B, seed)
%MANTISSA_TESTS Mantissa-based Benford tests with Monte Carlo p-values
%
%   S = mantissa_tests(x, B, seed) computes arc, mean, variance, joint,
%   KS, Watson U2, and d* tests on the mantissa m = frac(log10(x)).
%
%   Output structure S contains one field per test, each with:
%       .stat   observed statistic
%       .pAsym  asymptotic/reference p-value
%       .pMC    Monte Carlo p-value

x = x(:);
x = x(isfinite(x) & x > 0);
N = numel(x);

if nargin >= 3 && ~isempty(seed)
    oldRng = rng;
    rng(seed);
    c = onCleanup(@() rng(oldRng));
end

m = get_mantissa(x);
obs = stats_from_m(m);

nullS = zeros(B, 8);
for r = 1:B
    nullS(r, :) = ext_vector(stats_from_m(rand(N, 1)));
end

obsV = ext_vector(obs);
%pMC = (1 + sum(nullS >= obsV, 1)) / (B + 1);
pMC = mc_pvalue(nullS, obsV);


S = struct();
S.N = N;

S.arc.stat = obs.arc;
S.arc.pAsym = exp(-obs.arc / 2);
S.arc.pMC = pMC(1);

S.mean.stat = obs.zMean;
S.mean.pAsym = erfc(abs(obs.zMean) / sqrt(2));
S.mean.pMC = pMC(2);

S.variance.stat = obs.zVar;
S.variance.pAsym = erfc(abs(obs.zVar) / sqrt(2));
S.variance.pMC = pMC(3);

S.joint.stat = obs.joint;
S.joint.pAsym = exp(-obs.joint / 2);
S.joint.pMC = pMC(4);

S.ks.stat = obs.D;
S.ks.pAsym = kolmogorov_pvalue(obs.D, N);
S.ks.pMC = pMC(5);

S.u2.stat = obs.U2;
S.u2.pAsym = watson_p(obs.U2, N);
S.u2.pMC = pMC(6);

S.dstar.stat = obs.dstar;
S.dstar.pAsym = NaN;
S.dstar.pMC = pMC(7);

% Kuiper test: scale invariant counterpart of the KS test (V = D+ + D-)
S.kuiper.stat = obs.V;
S.kuiper.pAsym = kuiper_pvalue(obs.V, N);
S.kuiper.pMC = pMC(8);

S.arcAngle = obs.arcAngle;

    function m = get_mantissa(x)
        m = mod(log10(x), 1);
        m(m > 1 - 1e-12) = 0;
    end

    function s = stats_from_m(m)
        Nloc = numel(m);

        th = 2 * pi * m;
        c = mean(cos(th));
        sn = mean(sin(th));
        s.arc = 2 * Nloc * (c^2 + sn^2);
        s.arcAngle = mod(atan2(sn, c), 2 * pi) / (2 * pi);

        mbar = mean(m);
        v = var(m, 1);
        s.zMean = sqrt(12 * Nloc) * (mbar - 0.5);
        s.zVar = sqrt(180 * Nloc) * (v - 1 / 12);
        s.joint = s.zMean^2 + s.zVar^2;

        u = sort(m);
        i = (1:Nloc)';
        Dp = max(i / Nloc - u);
        Dm = max(u - (i - 1) / Nloc);
        s.D = max(Dp, Dm);
        s.V = Dp + Dm;

        W2 = 1 / (12 * Nloc) + sum((u - (2 * i - 1) / (2 * Nloc)).^2);
        s.U2 = W2 - Nloc * (mbar - 0.5)^2;

        b = log10(1 + 1 ./ (1:9)');
        d1 = min(max(floor(10.^m), 1), 9);
        p = accumarray(d1, 1, [9 1]) / Nloc;
        s.dstar = dstar_from_p(p, b);
    end

    function ds = dstar_from_p(p, b)
        p = p(:);
        b = b(:);
        dmax = sqrt(1 - 2 * min(b) + sum(b.^2));
        ds = norm(p - b) / dmax;
    end

    function v = ext_vector(s)
        v = [s.arc, abs(s.zMean), abs(s.zVar), s.joint, s.D, s.U2, s.dstar, s.V];
    end

    function p = watson_p(U2, Nloc)
        U2s = (U2 - 0.1 / Nloc + 0.1 / Nloc^2) * (1 + 0.8 / Nloc);
        k = 1:100;
        p = 2 * sum((-1).^(k - 1) .* exp(-2 * k.^2 * pi^2 * U2s));
        p = min(max(p, 0), 1);
    end
end

%% Internal dependencies: plotting and auxiliary tools

function txt = normalizeTextList(x)
if ischar(x) || isstring(x)
    txt = string(x);
elseif iscellstr(x)
    txt = string(x);
else
    txt = string.empty(1, 0);
    return
end

txt = lower(strtrim(txt(:).'));
txt = txt(txt ~= "");
end

% -------------------------------------------------------------------------
% function counts = countDigits(digits, categories)
%     counts = zeros(numel(categories), 1);
%     for i = 1:numel(categories)
%         counts(i) = sum(digits == categories(i));
%     end
% end
% the above is slower. Vectorized with accumarray.
function counts = countDigits(digits, categories)
[~, idx] = ismember(digits(:), categories(:));
idx = idx(idx > 0);
counts = accumarray(idx, 1, [numel(categories) 1], @sum, 0);
end


% -------------------------------------------------------------------------

function pMC = mc_pvalue(nullStat, obsStat)
%MCPVALUE Monte Carlo p-value from simulated null statistics
%
%   pMC = mc_pvalue(nullStat, obsStat) computes
%   (1 + number of simulated statistics >= observed statistic) / (B + 1).
%
%   nullStat is B-by-K, obsStat is 1-by-K or K-by-1.

obsStat = obsStat(:).';
if size(nullStat, 2) ~= numel(obsStat)
    nullStat = nullStat.';
end

pMC = (1 + sum(nullStat >= obsStat, 1)) / (size(nullStat, 1) + 1);
pMC = pMC(:);
end


% -------------------------------------------------------------------------

function out = duplication_tests(x, varargin)
%DUPLICATION_TESTS Duplication diagnostics for Benford analysis
%
%   out = duplication_tests(x) counts exact duplicates.
%   out = duplication_tests(x, 'Tol', tol) uses an absolute tolerance.
%   out = duplication_tests(x, 'RelTol', relTol) uses a relative tolerance.
%
%   The output splits duplicate excess counts into four buckets:
%       dupNegMult10  duplicated negative values that are multiples of 10
%       dupNegOther   duplicated negative values that are not multiples of 10
%       dupPosOther   duplicated nonnegative values that are not multiples of 10
%       dupPosMult10  duplicated nonnegative values that are multiples of 10
%
%   Values that are multiples of 10 are useful to track because they often
%   reveal rounding, coarse reporting, truncation, or unit/measurement
%   limits. A concentration of duplicated values ending in 0 can therefore
%   signal artificial or low-resolution data and may also distort Benford
%   digit patterns.

p = inputParser;
p.addParameter('Tol', [], @(v) isempty(v) || (isscalar(v) && isfinite(v) && v >= 0));
p.addParameter('RelTol', [], @(v) isempty(v) || (isscalar(v) && isfinite(v) && v >= 0));
p.parse(varargin{:});

absTol = p.Results.Tol;
relTol = p.Results.RelTol;

x = x(:);
x = x(isfinite(x));

out = struct( ...
    'dupNegMult10', 0, ...
    'dupNegOther', 0, ...
    'dupPosOther', 0, ...
    'dupPosMult10', 0, ...
    'z', zeros(4, 1), ...
    'nObs', 0, ...
    'nDup', 0);

if isempty(x)
    return
end

% Identify repeated values.
% If no tolerance is requested, use exact matching.
% Otherwise, group values that are close enough according to Tol/RelTol.
if isempty(absTol) && isempty(relTol)
    [u, ~, ic] = unique(x, 'stable');
else
    if isempty(absTol)
        absTol = 0;
    end
    if isempty(relTol)
        relTol = 0;
    end

    scale = max(abs(x));
    if scale == 0
        scale = 1;
    end

    effTol = max(absTol, relTol * scale);
    [u, ~, ic] = uniquetol(x, effTol, 'DataScale', 1, 'ByRows', false);
end

% Count how many times each unique value appears.
freq = accumarray(ic, 1);
dupExcess = max(freq - 1, 0);   % extra copies beyond the first occurrence

% Classify each unique value by sign and by whether it is a multiple of 10.
% The multiple-of-10 bucket is intentionally tracked because many repeated
% values ending in 0 often reflect rounding or coarse measurement/reporting.
isNeg = u < 0;
isMult10 = (u ~= 0) & (mod(abs(u), 10) == 0);

% Count duplicate excess in each bucket.
out.dupNegMult10 = sum(dupExcess(isNeg & isMult10));
out.dupNegOther  = sum(dupExcess(isNeg & ~isMult10));
out.dupPosOther  = sum(dupExcess(~isNeg & ~isMult10));
out.dupPosMult10 = sum(dupExcess(~isNeg & isMult10));

% Simple standardized scores for the four buckets.
counts = [out.dupNegMult10, out.dupNegOther, out.dupPosOther, out.dupPosMult10];
mu = mean(counts);
sg = std(counts);

if sg > 0
    out.z = ((counts - mu) / sg).';
else
    out.z = zeros(4, 1);
end

out.nObs = numel(x);
out.nDup = sum(dupExcess);
end

% ----------------------------------------------------------------

function B = chooseMonteCarloB(halfWidth, confLevel, pRef)
%chooseMonteCarloB Choose MC replications from target CI half-width
%
%   B = chooseMonteCarloB(halfWidth,confLevel,pRef) returns the minimum
%   number of replications needed so that the approximate two-sided Monte
%   Carlo CI for a proportion has half-width at most halfWidth.

z = sqrt(2) * erfinv(confLevel);   % two-sided normal quantile
pRef = min(max(pRef, eps), 1 - eps);

B = ceil((z^2 * pRef * (1 - pRef)) / (halfWidth^2));
B = max(B, 1);
end

% ----------------------------------------------------------------

function B = chooseTwoStageB(pv, nTailTarget)
%chooseTwoStageB Choose MC size for conditional two-stage critical values
%
%   B = chooseTwoStageB(pv,nTailTarget) returns the total number of Monte
%   Carlo replications needed so that the expected number of replications
%   in the conditional tail is at least nTailTarget.

pTail = pv / 100;
B = ceil(nTailTarget / pTail);
B = max(B, 1);
end


% ----------------------------------------------------------------
%% Internal dependencies: Kuiper and Benford-savvy tests

function p = kolmogorov_pvalue(D, n)
%kolmogorov_pvalue Asymptotic p-value of the KS statistic (Stephens' correction)

lam = (sqrt(n) + 0.12 + 0.11 / sqrt(n)) * D;
if lam < 0.3
    p = 1;      % the alternating series is useless here, and p is ~1
    return
end
k = 1:100;
p = 2 * sum((-1).^(k - 1) .* exp(-2 * k.^2 * lam^2));
p = min(max(p, 0), 1);
end

% ----------------------------------------------------------------

function p = kuiper_pvalue(V, n)
%kuiper_pvalue Asymptotic p-value of the Kuiper statistic V = D+ + D-
%
%   Uses Kuiper's limiting distribution with the finite-sample correction
%   of Stephens (1970): lambda = (sqrt(n) + 0.155 + 0.24/sqrt(n)) * V.

lam = (sqrt(n) + 0.155 + 0.24 / sqrt(n)) * V;
if lam < 0.4
    p = 1;
    return
end
k = 1:100;
p = 2 * sum((4 * k.^2 * lam^2 - 1) .* exp(-2 * k.^2 * lam^2));
p = min(max(p, 0), 1);
end

% ----------------------------------------------------------------

function p = chi2_upper(q, df)
%chi2_upper Upper tail probability of a chi2(df) variable (base MATLAB)
p = gammainc(q / 2, df / 2, 'upper');
end

% ----------------------------------------------------------------

function p = upper_tail_p(T)
%upper_tail_p Rank-based p-values p(i) = #{j: T(j) >= T(i)} / numel(T)
%
%   Ties are handled by giving every tied value the same p-value.

T = T(:);
N = numel(T);
[Ts, ord] = sort(T);
isNew = [true; diff(Ts) > 0];
firstPos = find(isNew);
pos = firstPos(cumsum(isNew));       % first sorted position of each tie group
p = zeros(N, 1);
p(ord) = (N - pos + 1) / N;
end

% ----------------------------------------------------------------

function F = savvy_cdf(u)
%savvy_cdf Null cdf of the fractional part of the significand, u in [0,1)
%
%   F(u) = sum_{d=1}^{9} log10((d+u)/d) = log10( Gamma(u+10) / (Gamma(u+1) 9!) )
%   (Proposition 1 in Barabesi et al., 2026).

F = (gammaln(u + 10) - gammaln(u + 1) - gammaln(10)) ./ log(10);
end

% ----------------------------------------------------------------

function T = savvy_null(n, B)
%savvy_null Null replicates of [KS2 KU2 Q1 Q2] under the Benford hypothesis
%
%   T is B-by-4. Samples are generated as S = 10^U, U ~ Uniform(0,1), and
%   processed in blocks of columns to limit memory usage.

C  = log10(exp(1));
dd = (1:9)';
pd = log10(1 + 1 ./ dd);
Sigma2 = C * diag(dd + 0.5) - C^2 * ones(9);   % var(Z2) under H0

iUp = (1:n)' / n;
iLo = ((1:n)' - 1) / n;

blk = max(1, floor(1e6 / n));
T = zeros(B, 4);
done = 0;
while done < B
    b  = min(blk, B - done);
    Sg = 10 .^ rand(n, b);
    Dg = floor(Sg);
    Fr = Sg - Dg;                               % fractional parts

    % KS2 and KU2 on the fractional part
    Fv = savvy_cdf(sort(Fr, 1));
    A  = max(iUp - Fv, [], 1);
    Bm = max(Fv - iLo, [], 1);

    % First-digit counts (Q1) and digit-wise sums of significands (Q2)
    cnt = zeros(9, b);
    sm  = zeros(9, b);
    for d = 1:9
        mask = (Dg == d);
        cnt(d, :) = sum(mask, 1);
        sm(d, :)  = sum(Sg .* mask, 1);
    end
    Q1 = sum((cnt - n * pd).^2 ./ (n * pd), 1);
    Y  = sm / n - C;
    Q2 = n * sum(Y .* (Sigma2 \ Y), 1);

    T(done + (1:b), :) = [max(A, Bm).', (A + Bm).', Q1.', Q2.'];
    done = done + b;
end
end

% ----------------------------------------------------------------

function S = savvy_tests(D1, Sig, B)
%savvy_tests Tests of Benford-savvy (manipulated-Benford) behaviour
%
%   S = savvy_tests(D1, Sig, B) computes the statistics of Barabesi,
%   Cerioli, Cerasa and Perrotta (2026):
%       KS2, KU2  Kolmogorov-Smirnov and Kuiper statistics on the fractional
%                 part of the significand
%       Q1        first-digit Pearson statistic
%       Q2        first-digit Hotelling statistic (sum-invariance)
%       QDelta    Q2 - Q1
%       GKS, GKU  min p-value combinations of (KS2,QDelta) and (KU2,QDelta)
%
%   Inputs:
%       D1   first digits (1..9)
%       Sig  significands in [1,10)
%       B    number of Monte Carlo replications
%
%   Each field of S has .stat, .pAsym and .pMC. All tests reject for large
%   values of the statistic, except GKS and GKU, whose statistic is a
%   p-value and which reject for small values (pMC is their p-value).
%   Monte Carlo p-values treat the observed sample as one more replicate,
%   so they are (1 + #{null >= obs})/(B + 1). The p-values entering GKS and
%   GKU are rank-based and computed on the same pool of B+1 samples.
%
%   Asymptotic p-values: Kolmogorov and Kuiper limits for KS2 and KU2 (their
%   exact null laws are distribution free, because F(<S>) is uniform under
%   H0), chi2(8), chi2(9) and chi2(1) upper tails for Q1, Q2 and QDelta (the
%   latter is a right-tail approximation only). Not available for G.
%
%   The null distribution assumes continuous data. Data with few
%   significant digits (e.g. rounded to 2 digits) will lead to rejection;
%   see the truncated-Benford simulation in the Supplement of the paper.

D1  = double(D1(:));
Sig = double(Sig(:));
n   = numel(D1);

if n * B > 2e9
    warning('FSDA:benford:SlowMC', ...
        ['The Monte Carlo for the savvy tests has n*B = %.2g and may be slow. ' ...
        'Consider a smaller B or the asymptotic p-values.'], n * B);
end

C  = log10(exp(1));
dd = (1:9)';
pd = log10(1 + 1 ./ dd);

% Observed statistics --------------------------------------------------
Fr  = min(max(Sig - D1, 0), 1 - eps);          % fractional part of the significand
Fv  = savvy_cdf(sort(Fr));
iv  = (1:n)';
A   = max(iv / n - Fv);
Bm  = max(Fv - (iv - 1) / n);
KS2 = max(A, Bm);
KU2 = A + Bm;

cnt = accumarray(D1, 1, [9 1]);
sm  = accumarray(D1, Sig, [9 1]);
Q1  = sum((cnt - n * pd).^2 ./ (n * pd));
Sigma2 = C * diag(dd + 0.5) - C^2 * ones(9);
y   = sm / n - C;
Q2  = n * (y' * (Sigma2 \ y));
QD  = Q2 - Q1;

% Monte Carlo pool (B null samples + observed sample in the last row) ---
T    = [savvy_null(n, B); KS2, KU2, Q1, Q2];
Tall = [T(:, 3), T(:, 4), T(:, 1), T(:, 2), T(:, 4) - T(:, 3)];  % Q1 Q2 KS2 KU2 QDelta

P = zeros(size(Tall));
for j = 1:size(Tall, 2)
    P(:, j) = upper_tail_p(Tall(:, j));
end

GKSall = min(P(:, 3), P(:, 5));
GKUall = min(P(:, 4), P(:, 5));
N = B + 1;
pGKS = sum(GKSall <= GKSall(end)) / N;
pGKU = sum(GKUall <= GKUall(end)) / N;

% Output ---------------------------------------------------------------
S = struct();
S.n = n;
S.B = B;
S.nullStats = T(1:B, :);   % null sample, columns KS2, KU2, Q1, Q2 (used by benfordplot)

S.Q1.stat      = Q1;   S.Q1.pAsym     = chi2_upper(Q1, 8);              S.Q1.pMC     = P(end, 1);
S.Q2.stat      = Q2;   S.Q2.pAsym     = chi2_upper(Q2, 9);              S.Q2.pMC     = P(end, 2);
S.KS2.stat     = KS2;  S.KS2.pAsym    = kolmogorov_pvalue(KS2, n);      S.KS2.pMC    = P(end, 3);
S.KU2.stat     = KU2;  S.KU2.pAsym    = kuiper_pvalue(KU2, n);          S.KU2.pMC    = P(end, 4);
S.QDelta.stat  = QD;   S.QDelta.pAsym = chi2_upper(max(QD, 0), 1);      S.QDelta.pMC = P(end, 5);
S.GKS.stat     = GKSall(end);  S.GKS.pAsym = NaN;  S.GKS.pMC = pGKS;
S.GKU.stat     = GKUall(end);  S.GKU.pAsym = NaN;  S.GKU.pMC = pGKU;
end


% ----------------------------------------------------------------

function L = lacunary_tests(U, m, alpha)
%lacunary_tests Single-observation lacunary harmonic test (Ferrari, 2026)
%
%   L = lacunary_tests(U, m, alpha) applies, to each element of the vector
%   of mantissae U in [0,1), the test of uniformity of Ferrari (2026),
%   based on the base-b (b = 10) expansion U_m = 0.D_1...D_m truncated to
%   m digits. With the lacunary frequencies k_{j,t} = j*b^(t-1),
%
%       Z_{j,t} = exp(2*pi*i*k_{j,t}*U_m),   j = 1..b-1,  t = 1..m
%       Zbar_j  = (1/m) * sum_t Z_{j,t}
%       T       = 2*m * sum_j |Zbar_j|^2   ~  chi2 with 2*(b-1) df under H0
%
%   (Theorem 1 of the paper; the limit is for increasing precision m).
%
%   Numerical note: k_{j,t}*U_m is never formed directly (it reaches about
%   9e9 for m = 10 and would lose precision in the argument of cos/sin).
%   With N = U_m*b^m (an integer), k_{j,t}*U_m = j*N/b^(m-t+1), so only the
%   fractional part of j*R_t/b^L is needed, where R_t = mod(N,b^L) holds
%   the digits D_t..D_m and L = m-t+1. This is exact in double precision
%   for m <= 15.
%
%   Output fields: m, base, df, alpha, T and pValue (n-by-1, one per
%   observation), reject (pValue < alpha), idxReject, nReject, nExpected
%   (alpha*n) and pBinom. pBinom is NOT part of the paper: it is the
%   upper-tail Binomial(n,alpha) probability of at least nReject flagged
%   observations, a rough dataset-level summary that assumes independent
%   observations and the nominal size of the single-observation test (in
%   simulation the size is slightly conservative, about 4.4% at 5% for
%   m = 10).

b  = 10;
df = 2 * (b - 1);
U  = U(:);
n  = numel(U);

% A mantissa within 1e-12 of 1 is a power of ten up to round-off
U(U >= 1 - 1e-12) = 0;

bm = b^m;
j  = 1:(b - 1);
T  = zeros(n, 1);

chunk = 2e5;                                   % limits memory (chunk-by-9 arrays)
for s = 1:chunk:n
    idx = s:min(s + chunk - 1, n);
    N   = mod(floor(U(idx) * bm), bm);         % integer N = U_m * b^m
    sC  = zeros(numel(idx), b - 1);
    sS  = zeros(numel(idx), b - 1);
    for t = 1:m
        bL = b^(m - t + 1);
        R  = mod(N, bL);                       % digits D_t, ..., D_m
        ph = mod(R * j, bL) / bL;              % frac(k_{j,t} * U_m), exact
        sC = sC + cos(2 * pi * ph);
        sS = sS + sin(2 * pi * ph);
    end
    T(idx) = (2 / m) * sum(sC.^2 + sS.^2, 2);  % = 2*m*sum_j |Zbar_j|^2
end

pValue = chi2_upper(T, df);
reject = pValue < alpha;
nRej   = sum(reject);

if nRej >= 1
    pBinom = betainc(alpha, nRej, n - nRej + 1);   % P(Bin(n,alpha) >= nRej)
else
    pBinom = 1;
end

L = struct();
L.m         = m;
L.base      = b;
L.df        = df;
L.alpha     = alpha;
L.T         = T;
L.pValue    = pValue;
L.reject    = reject;
L.idxReject = find(reject);
L.nReject   = nRej;
L.nExpected = alpha * n;
L.pBinom    = pBinom;
end


%FScategory:UTISTAT
