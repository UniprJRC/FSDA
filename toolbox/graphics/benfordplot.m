function h = benfordplot(out, varargin)
%benfordplot draws the graphical diagnostics of a Benford analysis
%
%<a href="matlab: docsearchFS('benfordplot')">Link to the help page for this function</a>
%
%   benfordplot takes the structure produced by function benford and draws
%   the requested diagnostic plots. Computation and plotting are
%   separated: no test is recomputed here, so the plots can be redrawn
%   (e.g. with a different significance level) without repeating the
%   Monte Carlo simulations stored in out.
%
%   The plots are of two kinds: those that compare the observed digit
%   frequencies with the Benford law (labels "first", "second",
%   "firsttwo", "deviation" and "heatmap") and those that accompany a
%   specific test of function benford (labels "mantissa", "twostage",
%   "savvy", "lacunary", "duplication" and "magnitude"). A plot whose
%   ingredients are not in out (for example "twostage" when the two-stage
%   test has not been run) is skipped, with a warning if it was
%   explicitly requested.
%
% Required input arguments:
%
%        out:  Structure returned by function benford. It must contain
%              at least the fields counts and mantissa produced by benford.
%              Data Types - struct
%
% Optional input arguments:
%
%      plots:  Text specifying the requested plots. Valid labels are:
%                "first", "second", "firsttwo" - observed against
%                       expected counts of the digit group, with the
%                       cumulative frequencies below;
%                "deviation" - observed minus expected proportion for each
%                       available digit group, with simultaneous
%                       (Bonferroni) acceptance bands under the Benford
%                       hypothesis. Cells outside the band are in red. The
%                       title reports the Mean Absolute Deviation (MAD)
%                       and the conformity category of Nigrini (2012). For
%                       large samples every cell tends to fall outside the
%                       band, and the MAD category then tells whether the
%                       deviation is practically relevant;
%                "heatmap" - Pearson residuals of the 90 first-two-digit
%                       cells arranged by first digit (rows) and second
%                       digit (columns). Cells flagged by the simultaneous
%                       band are labelled in bold. Round-number spikes
%                       and empty cells stand out immediately;
%                "magnitude" - histogram of log10(x) and counts per order
%                       of magnitude. A sample that spans several orders
%                       of magnitude is a prerequisite for the Benford
%                       hypothesis to be plausible;
%                "mantissa" - histogram of the mantissa and plot of the
%                       empirical cdf minus the uniform cdf, with the
%                       Kolmogorov-Smirnov band and the Kuiper strip. The
%                       curve leaves the KS band when KS rejects, and does
%                       not fit within the strip (of height equal to the
%                       Kuiper critical value) when Kuiper rejects;
%                "twostage" - for MAD and chi2, null distribution of the
%                       statistic of the two-stage test for H12 and, for
%                       H1 and H2, null distribution conditional on H12
%                       being rejected, with the critical value and the
%                       observed value;
%                "savvy" - Monte Carlo null sample of (Q1,Q_Delta),
%                       (KS2,Q_Delta) and (KU2,Q_Delta) with the observed
%                       point. Shows how a sample can pass the first-digit
%                       test Q1 and still fail on Q_Delta;
%                "lacunary" - QQ plot of -log10 of the p-values of the
%                       lacunary test against uniform expectation, with
%                       pointwise bands, and histogram of the p-values.
%                       If the points lie above the band, more
%                       observations are flagged than chance explains;
%                "duplication" - duplicate excess by value type;
%                "all" - every plot above for which out contains the
%                       necessary results.
%              Labels are not case sensitive. A logical or numeric scalar
%              is also accepted (true or 1 stands for "all"). Default is
%              "all".
%              Example - 'plots',["deviation","heatmap"]
%              Data Types - logical|string|char|cell
%
%      alpha:  Significance level for the plot bands and critical values.
%              Positive scalar, at most one half. Used for the simultaneous
%              bands of
%              "deviation" and "heatmap", for the critical values of
%              "mantissa", "twostage" and "savvy" and for the bands of
%              "lacunary". The critical values of the two-stage test were
%              computed by benford with its own level (option PVTS of
%              benford) and are not recomputed here. If empty (default),
%              the level stored in out.options.alpha is used, or 0.05 if
%              not available.
%              Example - 'alpha',0.01
%              Data Types - double
%
%      ngrid:  Positive integer. Number of grid points used by the plots
%              which would otherwise draw one marker per observation
%              ("mantissa" and "lacunary"). Default is 1000.
%              Example - 'ngrid',500
%              Data Types - double
%
%     docked:  Logical scalar. If true (default) the figures are docked.
%              If false each figure opens in its own window.
%              Example - 'docked',false
%              Data Types - logical
%
% Output:
%
% Output:
%
%          h:  Structure containing the handles of the figures, one field
%              per drawn plot, named after the label.
%              h.first      = first-digit plot: observed and expected counts,
%                              with cumulative observed and expected
%                              frequencies below.
%              h.second     = second-digit plot: observed and expected
%                              counts, with cumulative observed and
%                              expected frequencies below.
%              h.firstTwo   = first-two-digit plot: observed and expected
%                              counts, with cumulative observed and
%                              expected frequencies below.
%              h.deviation  = observed minus expected proportions for the
%                              available digit groups, with simultaneous
%                              Bonferroni bands under the Benford null.
%                              Bars outside the band are highlighted and
%                              the title reports the MAD and the
%                              conformity category.
%              h.heatmap    = first-two-digit Pearson residual heatmap,
%                              arranged by first digit (rows) and second
%                              digit (columns); cells outside the
%                              simultaneous cutoff are shown in bold.
%              h.magnitude  = magnitude overview: histogram of log10(x)
%                              and counts per order of magnitude.
%              h.mantissa   = mantissa diagnostics: histogram of the
%                              mantissa and ECDF minus uniform cdf, with
%                              KS band and Kuiper strip.
%              h.duplication = duplicate-excess diagnostics by value type
%                              and corresponding standardized z-scores.
%              h.twostage   = two-stage test diagnostics: null
%                              distributions, critical values, and
%                              observed statistics for MAD and chi2.
%              h.savvy      = Benford-savvy Monte Carlo diagnostics:
%                              null sample of Q1/QDelta, KS2/QDelta and
%                              KU2/QDelta with the observed point.
%              h.lacunary   = lacunary-test diagnostics: QQ plot of
%                              p-values against uniform expectation, with
%                              pointwise bands, and a histogram of the
%                              p-values.
%              Fields are present only for plots that have been drawn.
%
% See also: benford
%
% References:
%
%    Nigrini, M.J. (2012). Benford's Law: Applications for Forensic
%    Accounting, Auditing, and Fraud Detection. Wiley, Hoboken, NJ.
%    (MAD conformity ranges used in "deviation")
%
%    Goodman, L.A. (1965). On simultaneous confidence intervals for
%    multinomial proportions. Technometrics, 7(2), 247-254.
%    (simultaneous bands for multinomial proportions; the Bonferroni
%    version based on the normal approximation is used here)
%
%    Kuiper, N.H. (1960). Tests concerning random points on a circle.
%    Proceedings of the Koninklijke Nederlandse Akademie van
%    Wetenschappen, Series A, 63, 38-47.
%
%    Stephens, M.A. (1970). Use of the Kolmogorov-Smirnov, Cramer-von
%    Mises and related statistics without extensive tables. Journal of
%    the Royal Statistical Society, Series B, 32(1), 115-122.
%    (finite-sample critical values for the KS and Kuiper bands)
%
%    For the tests behind the "twostage", "mantissa", "savvy" and
%    "lacunary" plots see the references of function benford.
%
% Copyright 2026.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('benfordplot')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% Deviation plots with simultaneous bands and MAD conformity category.
    % Clean sample: the bars stay within the bands and the MAD is small.
    rng(1)
    V = 10.^(3*rand(2000,1));
    out = benford(V,'Plots',false,'Display',false);
    h = benfordplot(out,'plots',"deviation");
%}
%
%{
    %% Heatmap of the first-two-digit Pearson residuals.
    % 5% of the amounts have been replaced by the round values 50, 500 and
    % 5000: a single hot cell appears in the row of first digit 5, column
    % of second digit 0.
    rng(1)
    n = 2000;
    V = 10.^(3*rand(n,1));
    idx = randperm(n,100);
    V(idx) = 50*10.^randi([0 2],100,1);
    out = benford(V,'Plots',false,'Display',false);
    h = benfordplot(out,'plots',["deviation","heatmap"]);
%}
%
%{
    % Plots of the mantissa and of the order of magnitude.
    % A narrow-range sample: little more than one tenth of a decade.
    rng(1)
    V = 100 + 15*rand(1000,1);
    out = benford(V,'Tests',"mantissa",'Plots',false,'Display',false);
    h = benfordplot(out,'plots',["magnitude","mantissa"]);
%}
%
%{
    % Redrawing with a different level without repeating the simulations.
    % The two-stage and mantissa plots are redrawn at level 0.01.
    rng(1)
    V = 10.^randn(1000,1);
    out = benford(V,'Tests',["chi2","twostage","mantissa"],'Plots',false);
    h = benfordplot(out,'plots',["twostage","mantissa"],'alpha',0.01);
%}
%
%{
    % Plots of all the tests (running time depends on the number of replications).
    % The "duplication" plot is drawn only if benford has been called with
    % that label in option Plots, because it needs out.duplication.
    rng(1)
    V = 10.^randn(1000,1);
    out = benford(V,'Tests',["chi2","twostage","mantissa","savvy","lacunary"], ...
        'Plots',false,'Display',false);
    h = benfordplot(out,'plots',"all");
%}

%% Beginning of code

% Check the input structure
if nargin < 1 || ~isstruct(out) || ~isfield(out, 'counts') || ~isfield(out, 'mantissa')
    error('FSDA:benfordplot:WrongInput', ...
        'The first argument must be the structure returned by benford.');
end

% Options: defaults, user options (not case sensitive)
options = struct('plots', "all", 'alpha', [], 'ngrid', 1000, 'docked', true);

UserOptions = varargin(1:2:length(varargin));
if ~isempty(UserOptions) && (length(varargin) ~= 2 * length(UserOptions))
    error('FSDA:benfordplot:WrongInputOpt', ...
        'Number of supplied options is invalid. Probably values for some parameters are missing.');
end

optNames = fieldnames(options);
for i = 1:2:numel(varargin)
    name = char(string(varargin{i}));
    k = find(strcmpi(name, optNames), 1);
    if isempty(k)
        error('FSDA:benfordplot:WrongInputOpt', ...
            'Non existent user option found -> %s', name);
    end
    options.(optNames{k}) = varargin{i + 1};
end

% Significance level
alpha = options.alpha;
if isempty(alpha)
    alpha = 0.05;
    if isfield(out, 'options') && isstruct(out.options) && ...
            isfield(out.options, 'alpha') && ~isempty(out.options.alpha)
        alpha = out.options.alpha;
    end
end
if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || alpha <= 0 || alpha > 0.5
    error('FSDA:benfordplot:WrongAlpha', 'alpha must be a scalar in (0, 0.5].');
end

% Grid and window style
ngrid = options.ngrid;
if ~isscalar(ngrid) || ~isfinite(ngrid) || ngrid < 10
    error('FSDA:benfordplot:WrongNgrid', 'ngrid must be a scalar not smaller than 10.');
end
ngrid = round(ngrid);
docked = logical(options.docked);

% Requested plots
allowed = ["first", "second", "firsttwo", "deviation", "heatmap", ...
    "magnitude", "mantissa", "duplication", "twostage", "savvy", ...
    "lacunary", "all"];

pl = options.plots;
if islogical(pl) || (isnumeric(pl) && isscalar(pl))
    if pl
        pl = "all";
    else
        pl = string.empty(1, 0);
    end
else
    pl = lower(strtrim(string(pl(:).')));
    pl = pl(pl ~= "");
end

if any(~ismember(pl, allowed))
    error('FSDA:benfordplot:WrongPlots', ...
        'plots must contain only: %s.', strjoin(cellstr(allowed), ', '));
end

explicit = ~any(pl == "all");
if any(pl == "all")
    pl = allowed(1:end - 1);
end

h = struct();

%% Digit comparison plots
for lab = ["first", "second", "firsttwo"]
    if any(pl == lab)
        [fld, xv, ttl, tickLabels] = digitSpec(char(lab));
        if isfield(out.counts, fld) && isfield(out, 'expected') && isfield(out.expected, fld)
            h.(fld) = plotDigitComparison(xv, out.counts.(fld), ...
                out.expected.(fld), ttl, tickLabels, docked);
        else
            skipPlot(explicit, char(lab), ...
                'the counts of this digit group are not in out (see option Digits of benford)');
        end
    end
end

%% Deviation plot with simultaneous bands and MAD category
if any(pl == "deviation")
    groups = {'first', 'second', 'firstTwo'};
    ok = isfield(out.counts, groups);
    if isfield(out, 'expected')
        ok = ok & isfield(out.expected, groups);
    else
        ok(:) = false;
    end
    groups = groups(ok);
    if ~isempty(groups)
        h.deviation = plotDeviation(out, groups, alpha, docked);
    else
        skipPlot(explicit, 'deviation', 'no digit counts in out');
    end
end

%% Heatmap of the first-two-digit Pearson residuals
if any(pl == "heatmap")
    if isfield(out.counts, 'firstTwo') && isfield(out, 'expected') && ...
            isfield(out.expected, 'firstTwo')
        h.heatmap = plotHeatmap(out, alpha, docked);
    else
        skipPlot(explicit, 'heatmap', ...
            'the first-two-digit counts are not in out (Digits must include 12)');
    end
end

%% Magnitude overview
if any(pl == "magnitude")
    if isfield(out, 'data') && ~isempty(out.data)
        h.magnitude = plotMagnitude(out.data, docked);
    else
        skipPlot(explicit, 'magnitude', 'out.data is not available');
    end
end

%% Mantissa: histogram and ECDF difference with KS band and Kuiper strip
if any(pl == "mantissa")
    h.mantissa = plotMantissa(out.mantissa, alpha, ngrid, docked);
end

%% Duplication diagnostics
if any(pl == "duplication")
    if isfield(out, 'duplication') && isstruct(out.duplication) && ...
            isfield(out.duplication, 'dupNegMult10')
        D = out.duplication;
        h.duplication = plotDuplication(D.dupNegMult10, D.dupNegOther, ...
            D.dupPosOther, D.dupPosMult10, D.z, docked);
    else
        skipPlot(explicit, 'duplication', 'out.duplication is not available');
    end
end

%% Two-stage test
if any(pl == "twostage")
    if isfield(out, 'twoStage') && isstruct(out.twoStage) && ...
            isfield(out.twoStage, 'madSim')
        h.twostage = plotTwoStage(out.twoStage, docked);
    else
        skipPlot(explicit, 'twostage', 'the two-stage test has not been run');
    end
end

%% Benford-savvy tests
if any(pl == "savvy")
    if isfield(out, 'savvy') && isstruct(out.savvy) && isfield(out.savvy, 'nullStats')
        h.savvy = plotSavvy(out.savvy, alpha, docked);
    else
        skipPlot(explicit, 'savvy', 'the savvy tests have not been run');
    end
end

%% Lacunary test
if any(pl == "lacunary")
    if isfield(out, 'lacunary') && isstruct(out.lacunary) && isfield(out.lacunary, 'pValue')
        h.lacunary = plotLacunary(out.lacunary, ngrid, docked);
    else
        skipPlot(explicit, 'lacunary', 'the lacunary test has not been run');
    end
end

end

% -------------------------------------------------------------------------
%% Internal dependencies: helpers

function skipPlot(explicit, label, why)
%skipPlot Warn that a requested plot cannot be drawn
if explicit
    warning('FSDA:benfordplot:Skipped', 'Plot "%s" skipped: %s.', label, why);
end
end

% -------------------------------------------------------------------------

function fig = newFigure(name, docked)
%newFigure Create a figure, docked or not
if docked
    fig = figure('Name', name, 'NumberTitle', 'off', 'WindowStyle', 'docked');
else
    fig = figure('Name', name, 'NumberTitle', 'off');
end
end

% -------------------------------------------------------------------------

function [fld, xv, ttl, tickLabels] = digitSpec(lab)
%digitSpec Field of out.counts, support, title and tick labels of a digit group
switch lab
    case 'first'
        fld = 'first';
        xv = 1:9;
        ttl = 'First digit distribution';
        tickLabels = cellstr(string(1:9));
    case 'second'
        fld = 'second';
        xv = 0:9;
        ttl = 'Second digit distribution';
        tickLabels = cellstr(string(0:9));
    otherwise
        fld = 'firstTwo';
        xv = 10:99;
        ttl = 'First two digits distribution';
        tickLabels = cellstr(string(10:99));
end
end

% -------------------------------------------------------------------------

function [xv, name] = groupAxis(f)
%groupAxis Support (column vector) and name of a digit group
switch f
    case 'first'
        xv = (1:9).';
        name = 'First digit';
    case 'second'
        xv = (0:9).';
        name = 'Second digit';
    otherwise
        xv = (10:99).';
        name = 'First two digits';
end
end

% -------------------------------------------------------------------------

function catg = madCategory(f, madv)
%madCategory Conformity category of the MAD according to Nigrini (2012)
switch f
    case 'first'
        thr = [0.006 0.012 0.015];
    case 'second'
        thr = [0.008 0.010 0.012];
    otherwise
        thr = [0.0012 0.0018 0.0022];
end
names = {'close conformity', 'acceptable conformity', ...
    'marginally acceptable conformity', 'nonconformity'};
catg = names{1 + sum(madv > thr)};
end

% -------------------------------------------------------------------------

function q = sortedQuantile(v, p)
%sortedQuantile Empirical quantile without toolbox functions
s = sort(v(:));
q = s(min(numel(s), max(1, ceil(p * numel(s)))));
end

% -------------------------------------------------------------------------

function q = kolmogorovQ(lam)
%kolmogorovQ Upper tail of the limiting Kolmogorov distribution
if lam < 0.3
    q = 1;
    return
end
k = 1:100;
q = 2 * sum((-1).^(k - 1) .* exp(-2 * k.^2 * lam^2));
q = min(max(q, 0), 1);
end

% -------------------------------------------------------------------------

function q = kuiperQ(lam)
%kuiperQ Upper tail of the limiting Kuiper distribution
if lam < 0.4
    q = 1;
    return
end
k = 1:100;
q = 2 * sum((4 * k.^2 * lam^2 - 1) .* exp(-2 * k.^2 * lam^2));
q = min(max(q, 0), 1);
end

% -------------------------------------------------------------------------

function c = kolmogorovCritical(alpha, n)
%kolmogorovCritical Critical value of the KS statistic (Stephens' correction)
lam = fzero(@(l) kolmogorovQ(l) - alpha, [0.3 6]);
c = lam / (sqrt(n) + 0.12 + 0.11 / sqrt(n));
end

% -------------------------------------------------------------------------

function c = kuiperCritical(alpha, n)
%kuiperCritical Critical value of the Kuiper statistic (Stephens' correction)
lam = fzero(@(l) kuiperQ(l) - alpha, [0.4 6]);
c = lam / (sqrt(n) + 0.155 + 0.24 / sqrt(n));
end

% -------------------------------------------------------------------------

function cm = divergingMap(m)
%divergingMap Blue-white-red colormap with m levels
t = linspace(0, 1, m).';
cm = [min(2 * t, 1), 1 - abs(2 * t - 1), min(2 * (1 - t), 1)];
end

% -------------------------------------------------------------------------
%% Internal dependencies: plots

function fig = plotDigitComparison(x, counts, expected, ttl, tickLabels, docked)
%plotDigitComparison Observed digit counts against Benford expectation
%
%   The first panel shows observed counts and expected counts.
%   The second panel shows cumulative observed and expected frequencies.

counts = counts(:);
expected = expected(:);
x = x(:);

fig = newFigure('Benford Digit Comparison', docked);

tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top panel: counts
nexttile(tl);
bar(x, counts, 0.6, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', [0.6 0.6 0.6]);
hold on
plot(x, expected, '-r', 'LineWidth', 2.5);
hold off
title(ttl)

if ~isempty(tickLabels)
    xticks(x)
    xticklabels(tickLabels)
end
xlim([min(x) - 0.5, max(x) + 0.5])
ylabel('Count')

% Bottom panel: cumulative frequencies
nexttile(tl);

obsCum = cumsum(counts);
expCum = cumsum(expected);

if obsCum(end) > 0
    obsCum = obsCum / obsCum(end);
end
if expCum(end) > 0
    expCum = expCum / expCum(end);
end

plot(x, obsCum, '-o', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.2 0.2 0.2]);
hold on
plot(x, expCum, '-r', 'LineWidth', 2.0);
hold off

if ~isempty(tickLabels)
    xticks(x)
    xticklabels(tickLabels)
end
xlim([min(x) - 0.5, max(x) + 0.5])
ylim([0 1])
ylabel('Cumulative frequency')
xlabel('Digit')
legend({'Observed', 'Expected'}, 'Location', 'best')
grid on
end

% -------------------------------------------------------------------------

function fig = plotDeviation(out, groups, alpha, docked)
%plotDeviation Observed minus expected proportions with simultaneous bands
%
%   For each digit group, bars show p_hat - p. The shaded region is the
%   Bonferroni band p_hat - p +/- z*sqrt(p(1-p)/n) with z the normal
%   quantile of level 1 - alpha/(2k), k being the number of cells, which
%   contains all the cells with probability about 1 - alpha under the
%   Benford hypothesis. Cells outside the band are in red.

nG = numel(groups);
fig = newFigure('Benford Deviations from Expected Proportions', docked);
tl = tiledlayout(fig, nG, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for g = 1:nG
    f = groups{g};
    c = out.counts.(f)(:);
    e = out.expected.(f)(:);
    p = e / sum(e);
    nn = sum(c);
    [xv, name] = groupAxis(f);

    dev = c / nn - p;
    k = numel(p);
    z = sqrt(2) * erfinv(1 - alpha / k);        % normal quantile of level 1-alpha/(2k)
    half = z * sqrt(p .* (1 - p) / nn);
    outside = abs(dev) > half;

    madv = mean(abs(dev));
    catg = madCategory(f, madv);

    ax = nexttile(tl);
    hold(ax, 'on')

    hnd = gobjects(0);
    leg = {};

    hband = fill(ax, [xv; flipud(xv)], [half; -flipud(half)], [0.80 0.88 1.00], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hnd(end + 1) = hband;
    leg{end + 1} = sprintf('%g%% simultaneous band', 100 * (1 - alpha));

    if any(~outside)
        hin = bar(ax, xv(~outside), dev(~outside), 0.7, ...
            'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
        hnd(end + 1) = hin;
        leg{end + 1} = 'Within band';
    end
    if any(outside)
        hout = bar(ax, xv(outside), dev(outside), 0.7, ...
            'FaceColor', [0.85 0.15 0.15], 'EdgeColor', 'none');
        hnd(end + 1) = hout;
        leg{end + 1} = 'Outside band';
    end
    yline(ax, 0, '-k');

    xlim(ax, [min(xv) - 0.5, max(xv) + 0.5]);
    if strcmp(f, 'firstTwo')
        xticks(ax, 10:10:90);
    else
        xticks(ax, xv);
    end
    ylabel(ax, 'Observed - expected proportion');
    title(ax, sprintf('%s: MAD = %.4f (%s)', name, madv, catg));
    legend(ax, hnd, leg, 'Location', 'best');
    grid(ax, 'on')
    if g == nG
        xlabel(ax, 'Digit');
    end
    hold(ax, 'off')
end
end

% -------------------------------------------------------------------------

function fig = plotHeatmap(out, alpha, docked)
%plotHeatmap Pearson residuals of the first-two-digit cells
%
%   Rows are the first digit (1-9), columns the second digit (0-9).

c = out.counts.firstTwo(:);
e = out.expected.firstTwo(:);
r = (c - e) ./ sqrt(e);

R = reshape(r, [10 9]).';              % row = first digit, column = second digit
k = numel(r);
zc = sqrt(2) * erfinv(1 - alpha / k);  % simultaneous (Bonferroni) cutoff
cl = max(4, min(max(abs(r)), 12));

fig = newFigure('Benford First-Two-Digit Residuals', docked);
ax = axes(fig);
imagesc(ax, 0:9, 1:9, R);
colormap(ax, divergingMap(256));
caxis(ax, [-cl cl]);
cb = colorbar(ax);
cb.Label.String = 'Pearson residual (observed - expected)/sqrt(expected)';

xticks(ax, 0:9);
yticks(ax, 1:9);
xlabel(ax, 'Second digit');
ylabel(ax, 'First digit');
title(ax, sprintf('First two digits: Pearson residuals (bold cells exceed the %g%% simultaneous cutoff)', ...
    100 * (1 - alpha)));

hold(ax, 'on')
for i = 1:9
    for j = 0:9
        d12 = 10 * i + j;
        flagged = abs(R(i, j + 1)) > zc;
        if flagged
            text(ax, j, i, num2str(d12), 'HorizontalAlignment', 'center', ...
                'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');
        else
            text(ax, j, i, num2str(d12), 'HorizontalAlignment', 'center', ...
                'FontSize', 7, 'Color', [0.45 0.45 0.45]);
        end
    end
end
hold(ax, 'off')
end

% -------------------------------------------------------------------------

function fig = plotMagnitude(x, docked)
%plotMagnitude Histogram of log10(x) and counts per order of magnitude

x = x(:);
x = x(isfinite(x) & x > 0);
lx = log10(x);

fig = newFigure('Benford Magnitude Overview', docked);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl);
histogram(ax1, lx, 'FaceColor', [0.239 0.149 0.659], 'EdgeColor', 'none');
xlabel(ax1, 'log_{10}(x)');
ylabel(ax1, 'Count');
title(ax1, sprintf('Range %.2f decades, std of log_{10}(x) = %.2f', ...
    max(lx) - min(lx), std(lx)));
grid(ax1, 'on')

ax2 = nexttile(tl);
d = floor(lx);
dv = min(d):max(d);
cnt = histcounts(d, [dv - 0.5, max(dv) + 0.5]);
bar(ax2, dv, cnt, 0.8, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
xlabel(ax2, 'Order of magnitude floor(log_{10}(x))');
ylabel(ax2, 'Count');
title(ax2, 'Observations per order of magnitude');
grid(ax2, 'on')

title(tl, 'Magnitude overview: a wide spread favours the Benford hypothesis');
end

% -------------------------------------------------------------------------

function fig = plotMantissa(mant, alpha, ngrid, docked)
%plotMantissa Mantissa histogram and ECDF difference with KS band and Kuiper strip
%
%   The ECDF difference F_n(u) - u is evaluated on a grid of ngrid+1
%   points, which does not require one marker per observation. The KS and
%   Kuiper statistics reported in the title are exact (computed from the
%   sorted mantissae). KS rejects when the curve leaves the band
%   +/- c_KS; Kuiper rejects when the range of the curve exceeds c_V, i.e.
%   when the curve does not fit within a horizontal strip of height c_V.

mant = mant(:);
mant = mant(isfinite(mant));
mant = mod(mant, 1);
n = numel(mant);

if n == 0
    fig = [];
    return
end

% ECDF difference on a grid
edges = linspace(0, 1, ngrid + 1);
cnt = histcounts(mant, edges);
u = edges(:);
D = [0; cumsum(cnt(:))] / n - u;

% Exact KS and Kuiper statistics
ms = sort(mant);
iv = (1:n).';
Dp = max(iv / n - ms);
Dm = max(ms - (iv - 1) / n);
Dks = max(Dp, Dm);
Vku = Dp + Dm;

cKS = kolmogorovCritical(alpha, n);
cKU = kuiperCritical(alpha, n);

fig = newFigure('Benford Mantissa Diagnostics', docked);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Histogram
ax1 = nexttile(tl);
nb = min(50, max(10, floor(sqrt(n) / 2)));
histogram(ax1, mant, 'BinLimits', [0 1], 'NumBins', nb, 'Normalization', 'pdf', ...
    'FaceColor', [0.239 0.149 0.659], 'EdgeColor', 'none');
hold(ax1, 'on')
yline(ax1, 1, '--r', 'LineWidth', 1.5);
hold(ax1, 'off')
xlim(ax1, [0 1]);
xlabel(ax1, 'Mantissa');
ylabel(ax1, 'Density');
title(ax1, 'Mantissa histogram');
grid(ax1, 'on')

% ECDF difference
ax2 = nexttile(tl);
hold(ax2, 'on')
mid = (max(D) + min(D)) / 2;
hstrip = patch(ax2, [0 1 1 0], mid + (cKU / 2) * [-1 -1 1 1], [0.85 0.85 0.85], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7);
hks = plot(ax2, [0 1], [cKS cKS], '--r', 'LineWidth', 1.3);
plot(ax2, [0 1], [-cKS -cKS], '--r', 'LineWidth', 1.3);
hcurve = plot(ax2, u, D, '-k', 'LineWidth', 1.6);
yline(ax2, 0, ':k');
hold(ax2, 'off')

ym = 1.1 * max([abs(D); cKS; abs(mid) + cKU / 2]);
ylim(ax2, [-ym ym]);
xlim(ax2, [0 1]);
xlabel(ax2, 'u');
ylabel(ax2, 'F_n(u) - u');
title(ax2, sprintf('KS D = %.4f (crit %.4f); Kuiper V = %.4f (crit %.4f)', ...
    Dks, cKS, Vku, cKU));
legend(ax2, [hcurve hks hstrip], ...
    {'Empirical - uniform cdf', sprintf('KS band (%g%%)', 100 * (1 - alpha)), ...
    'Kuiper strip (fits if V < crit)'}, 'Location', 'best');
grid(ax2, 'on')

title(tl, 'Mantissa diagnostics');
end

% -------------------------------------------------------------------------

function fig = plotDuplication(d10n, d0n, d0p, d10p, z, docked)
%plotDuplication Visualize duplicate excess counts for Benford analysis
%
%   fig = plotDuplication(d10n, d0n, d0p, d10p, z, docked) shows two
%   diagnostics:
%       1) how many duplicate observations fall in each value bucket
%       2) the corresponding standardized scores (z-scores)
%
%   The four duplicate buckets are:
%       d10n  duplicated negative values that are multiples of 10
%       d0n   duplicated negative values that are not multiples of 10
%       d0p   duplicated nonnegative values that are not multiples of 10
%       d10p  duplicated nonnegative values that are multiples of 10
%
%   Repeated rounded or coarse-grained values often indicate low
%   measurement precision, unit conversion effects, truncation, or other
%   preprocessing artifacts. Such patterns can also distort the digit
%   distribution used in Benford analysis.

counts = [d10n, d0n, d0p, d10p];
counts = counts(:).';

if nargin < 5 || isempty(z)
    z = [];
else
    z = z(:).';
end

if isempty(counts) || all(~isfinite(counts))
    fig = [];
    return
end

countLabels = { ...
    'Negative, x10', ...
    'Negative, other', ...
    'Nonnegative, other', ...
    'Nonnegative, x10'};

fig = newFigure('Duplicate Excess Diagnostics', docked);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: duplicate excess counts
nexttile(tl);
bar(counts, 'FaceColor', [0.239 0.149 0.659], 'EdgeColor', 'none');
xticks(1:numel(counts));
xticklabels(countLabels);
ylabel('Count');
title('Duplicate Excess by Value Type');
grid on

% Panel 2: standardized duplicate excess
nexttile(tl);

if isempty(z) || all(~isfinite(z))
    text(0.5, 0.5, 'No z-scores available', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
    axis off
else
    bar(z, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
    xticks(1:numel(z));
    if numel(z) == numel(countLabels)
        xticklabels(countLabels);
    end
    yline(0, '-k', 'LineWidth', 1);
    yline(1.96, '--r', 'LineWidth', 1);
    yline(-1.96, '--r', 'LineWidth', 1);
    ylabel('z-score');
    title('Standardized Duplicate Excess');
    grid on
end

title(tl, 'Duplication diagnostics');
end

% -------------------------------------------------------------------------

function fig = plotTwoStage(T, docked)
%plotTwoStage Null distributions, critical values and observed statistics
%
%   Top row: MAD. Bottom row: chi2-type statistic. First column: null
%   distribution of the statistic for H12 (first two digits). Second and
%   third columns: null distribution for H1 and H2 among the simulated
%   samples in which H12 is rejected, which is the distribution that the
%   critical values for H1 and H2 refer to.

lab = T.labels;
i12 = find(strcmp(lab, '12'), 1);
i1  = find(strcmp(lab, '1'), 1);
i2  = find(strcmp(lab, '2'), 1);
order = [i12 i1 i2];

colTitles = {'H_{12}: first two digits', ...
    'H_1 (among samples rejecting H_{12})', ...
    'H_2 (among samples rejecting H_{12})'};
rowNames = {'MAD', 'chi^2 statistic'};

sims = {T.madSim, T.chi2Sim};
cvs  = {T.cvMad(:), T.cvChi2(:)};
obss = {T.obsMad(:), T.obsChi2(:)};

fig = newFigure('Benford Two-Stage Test', docked);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:2
    S  = sims{r};
    cv = cvs{r};
    ob = obss{r};

    mask = S(:, i12) >= cv(i12);
    if ~any(mask)
        mask = true(size(mask));
    end

    for c = 1:3
        j = order(c);
        if c == 1
            v = S(:, j);
        else
            v = S(mask, j);
        end

        ax = nexttile(tl);
        hold(ax, 'on')
        hH = histogram(ax, v, 'Normalization', 'pdf', 'EdgeColor', 'none', ...
            'FaceColor', [0.6 0.6 0.6]);
        hC = xline(ax, cv(j), '--r', 'LineWidth', 1.5);
        hO = xline(ax, ob(j), '-k', 'LineWidth', 2);
        hold(ax, 'off')

        lo = min([v(:); ob(j); cv(j)]);
        hi = max([v(:); ob(j); cv(j)]);
        pad = 0.05 * (hi - lo);
        if pad > 0
            xlim(ax, [lo - pad, hi + pad]);
        end

        if r == 1
            title(ax, colTitles{c});
        end
        xlabel(ax, rowNames{r});
        if c == 1
            ylabel(ax, 'Null density');
        end
        if r == 1 && c == 1
            legend(ax, [hH hC hO], {'Null distribution', 'Critical value', 'Observed'}, ...
                'Location', 'best');
        end
        grid(ax, 'on')
    end
end

title(tl, 'Two-stage test: null distributions, critical values and observed statistics');
end

% -------------------------------------------------------------------------

function fig = plotSavvy(S, alpha, docked)
%plotSavvy Monte Carlo null sample of (Q1,Q_Delta), (KS2,Q_Delta), (KU2,Q_Delta)
%
%   S.nullStats has columns KS2, KU2, Q1, Q2 for the B null samples.
%   Q_Delta = Q2 - Q1. The red star is the observed sample. A sample which
%   mimics the first-digit law but not the fractional part of the
%   significand is well within the null range of Q1 but far in the tail of
%   Q_Delta.

N = S.nullStats;
B = size(N, 1);
QD = N(:, 4) - N(:, 3);

nulls = {N(:, 3), N(:, 1), N(:, 2)};
obsv  = [S.Q1.stat, S.KS2.stat, S.KU2.stat];
xlab  = {'Q_1 (first-digit Pearson)', 'KS_2', 'KU_2'};

critQD = sortedQuantile(QD, 1 - alpha);
sel = unique(round(linspace(1, B, min(B, 3000))));

fig = newFigure('Benford-Savvy Tests', docked);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for c = 1:3
    xn = nulls{c};
    ax = nexttile(tl);
    hold(ax, 'on')
    hN = scatter(ax, xn(sel), QD(sel), 6, [0.6 0.6 0.6], 'filled', ...
        'MarkerFaceAlpha', 0.4);
    hX = xline(ax, sortedQuantile(xn, 1 - alpha), '--b', 'LineWidth', 1.3);
    yline(ax, critQD, '--b', 'LineWidth', 1.3);
    hO = plot(ax, obsv(c), S.QDelta.stat, 'p', 'MarkerSize', 14, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    hold(ax, 'off')

    xl = [min([xn; obsv(c)]), max([xn; obsv(c)])];
    yl = [min([QD; S.QDelta.stat]), max([QD; S.QDelta.stat])];
    xp = 0.05 * (xl(2) - xl(1));
    yp = 0.05 * (yl(2) - yl(1));
    xlim(ax, [xl(1) - xp, xl(2) + xp]);
    ylim(ax, [yl(1) - yp, yl(2) + yp]);

    xlabel(ax, xlab{c});
    if c == 1
        ylabel(ax, 'Q_\Delta = Q_2 - Q_1');
        legend(ax, [hN hX hO], {'Null sample', ...
            sprintf('%g%% critical values', 100 * (1 - alpha)), 'Observed'}, ...
            'Location', 'best');
    end
    grid(ax, 'on')
end

title(tl, sprintf(['Benford-savvy tests (B = %d). Monte Carlo p-values: ' ...
    'G_{KS} = %.3g, G_{KU} = %.3g'], B, S.GKS.pMC, S.GKU.pMC));
end

% -------------------------------------------------------------------------

function fig = plotLacunary(L, ngrid, docked)
%plotLacunary p-value QQ plot and p-value histogram of the lacunary test
%
%   Under the null hypothesis the p-values are (about) uniform. The QQ
%   plot shows -log10 of the ordered p-values against -log10 of the
%   expected uniform order statistics, with pointwise 1-alpha bands from
%   the Beta(i, n-i+1) distribution of the i-th order statistic. Points
%   above the horizontal line are the flagged observations.

p = L.pValue(:);
n = numel(p);
a = L.alpha;
ps = sort(p);

idx = unique(round(logspace(0, log10(n), min(n, ngrid))));
idx = idx(idx >= 1 & idx <= n);
idx = idx(:);

xq = -log10((idx - 0.5) / n);
yq = -log10(max(ps(idx), realmin));
lo = betaincinv(a / 2, idx, n - idx + 1);
hi = betaincinv(1 - a / 2, idx, n - idx + 1);
ylo = -log10(max(hi, realmin));
yhi = -log10(max(lo, realmin));

fig = newFigure('Benford Lacunary Test', docked);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% QQ plot
ax1 = nexttile(tl);
hold(ax1, 'on')
hb = fill(ax1, [xq; flipud(xq)], [ylo; flipud(yhi)], [0.80 0.88 1.00], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8);
hd = plot(ax1, [0 max(xq)], [0 max(xq)], '-r', 'LineWidth', 1.3);
hp = plot(ax1, xq, yq, '.k', 'MarkerSize', 9);
ht = yline(ax1, -log10(a), '--k');
hold(ax1, 'off')
xlabel(ax1, 'Expected -log_{10}(p)');
ylabel(ax1, 'Observed -log_{10}(p)');
title(ax1, sprintf('QQ plot: %d of %d flagged at %g (about %.1f expected)', ...
    L.nReject, n, a, L.nExpected));
legend(ax1, [hp hd hb ht], {'Observations', 'Uniform p-values', ...
    sprintf('%g%% pointwise band', 100 * (1 - a)), 'Flag threshold'}, ...
    'Location', 'northwest');
grid(ax1, 'on')

% Histogram of the p-values
ax2 = nexttile(tl);
histogram(ax2, p, 'BinLimits', [0 1], 'NumBins', 20, ...
    'FaceColor', [0.239 0.149 0.659], 'EdgeColor', 'none');
hold(ax2, 'on')
yline(ax2, n / 20, '--r', 'LineWidth', 1.5);
hold(ax2, 'off')
xlabel(ax2, 'p-value');
ylabel(ax2, 'Count');
title(ax2, 'Histogram of the p-values (flat under H_0)');
grid(ax2, 'on')

title(tl, sprintf('Lacunary test, m = %d mantissa digits', L.m));
end

%FScategory:UTISTAT
