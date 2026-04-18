function out = mdMCARdistTest(Y, varargin)
% mdMCARdistTest
% Bootstrap test based on the change of Mahalanobis distances for complete rows
% when mu and Sigma are estimated:
%   (i) from complete rows only
%   (ii) from all rows using EM/TEM with missing values.
%
% The null is that the observed perturbation is compatible with MCAR.
%
% Required input
% --------------
% Y : n x p data matrix with NaNs.
%
% Optional input arguments (name/value pairs)
% -------------------------------------------
% 'alpha'      : trimming level. Default is 0.
% 'method'     : rescaling method used in mdTEM and mdPartialMD2full.
%                Default is 'pri'.
% 'nsimul'          : number of bootstrap replicates. Default is 499.
% 'conflev'    : confidence level for optional CI. Default 0.95.
% 'tol'        : tolerance passed to mdTEM. Default 1e-10.
% 'plots'      : if true, produce diagnostic plots. Default false.
%
% Output
% ------
% out is a structure containing:
%   out.pvalue
%   out.Tobs
%   out.Tboot
%   out.alpha
%   out.method
%   out.stat
%   out.n
%   out.p
%   out.nComplete
%   out.completeIdx
%   out.d2_cc
%   out.d2_all
%   out.logRatio
%   out.ciBoot
%   out.muHat
%   out.SigHat
%
% Notes
% -----
% The bootstrap is parametric under a Gaussian model:
%   Y*_full ~ N_p(muHat, SigHat)
% and the observed missingness mask is then imposed on Y*_full.
% In this way missingness is independent of the generated data, as under MCAR.
%
% Dependencies on path:
%   mdEM.m
%   mdTEM.m
%   mdPartialMD.m
%   mdPartialMD2full.m
%   mahalFS.m

% Default options
alpha   = 0;
method  = 'pri';
nsimul  = 499;
conflev = 0.95;
tol     = 1e-10;
plots  = false;

if mod(length(varargin),2) ~= 0
    error('Optional arguments must be supplied in name/value pairs.');
end

for i = 1:2:length(varargin)
    name = lower(varargin{i});
    val  = varargin{i+1};
    switch name
        case 'alpha'
            alpha = val;
        case 'method'
            method = val;
        case 'nsimul'
            nsimul = val;
        case 'conflev'
            conflev = val;
        case 'tol'
            tol = val;
        case 'plots'
            plots = val;
        otherwise
            error('Unknown option: %s', varargin{i});
    end
end


[n,p] = size(Y);
maskMiss = isnan(Y);
completeIdx = all(~maskMiss,2);
nComplete = sum(completeIdx);

if nComplete < p + 2
    error('Too few complete rows to compute the reference complete-case covariance.');
end

% Reference distances from complete rows only
Ycc = Y(completeIdx,:);
muCC = mean(Ycc,1);
SigCC = cov(Ycc);
d2_cc = mahalFS(Ycc, muCC, SigCC);

% Fit EM/TEM on all data
[d2_all_cc, muHat, SigHat] = local_fit_and_get_complete_distances(Y, completeIdx, alpha, method, tol);

% 4 observed test statistic
Tobs = local_statistic(d2_cc, d2_all_cc);

% Bootstrap under MCAR for the 4 tests
Tboot = nan(nsimul,4);

% Make Sigma numerically SPD if needed
SigGen = local_make_spd(SigCC);

% Cholesky once if possible
R = chol(SigGen, 'upper');

for j = 1:nsimul
    % Generate full data under fitted model
    YfullStar = randn(n,p) * R + muCC(:)';

    % Impose the original missingness pattern
    Ystar = YfullStar;
    Ystar(maskMiss) = NaN;

    % Complete-case reference in the bootstrap world:
    % use the truly complete generated rows corresponding to the same complete rows
    YccStar = YfullStar(completeIdx,:);
    muCCStar = mean(YccStar,1);
    SigCCStar = cov(YccStar);
    d2_cc_star = mahalFS(YccStar, muCCStar, SigCCStar);

    % All-data EM/TEM distances for the same complete rows
        [d2_all_cc_star] = local_fit_and_get_complete_distances(Ystar, completeIdx, alpha, method, tol);
       % Store the 4 statistics
        Tboot(j,:) = local_statistic(d2_cc_star, d2_all_cc_star);
end

% Remove failed bootstrap replicates
Tboot = rmmissing(Tboot);

if isempty(Tboot)
    error('All bootstrap replicates failed.');
end

% p-value for the 4 tests
        pvalue = (1 + sum(abs(Tboot) >= abs(Tobs))) / (length(Tboot) + 1);

% Bootstrap CI for the target parameter
alphaCI = 1 - conflev;
ciBoot = quantile(Tboot, [alphaCI/2 1-alphaCI/2]);

% Output
out = struct;
out.pvalue      = pvalue;
out.Tobs        = Tobs;
out.Tboot       = Tboot;
out.alpha       = alpha;
out.method      = method;
out.nComplete   = nComplete;
out.completeIdx = completeIdx;
out.d2_cc       = d2_cc;
out.d2_all      = d2_all_cc;
out.logRatio    = log(d2_all_cc ./ d2_cc);
out.ciBoot      = ciBoot;
out.muHat       = muHat;
out.SigHat      = SigHat;

% Optional plots
if plots
    statNames = {'median log-ratio', ...
                 'mean log-ratio', ...
                 'median difference', ...
                 'mean difference'};

    figure;

    % Create a 2x3 tiled layout
    tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

    % ---- LEFT: scatter plot spanning 2 rows ----
    nexttile([2 1])   % span 2 rows, 1 column
    scatter(d2_cc, d2_all_cc, 'o')
    refline(1,0)
    xlabel('Complete-case distance')
    ylabel('Distance from all-data EM/TEM')
    title(sprintf('alpha = %g', alpha))
    box on

    % ---- RIGHT: 4 bootstrap panels ----
    for j = 1:4
        nexttile
        histogram(Tboot(:,j))
        hold on
        xline(Tobs(j), 'r', 'LineWidth', 1.5)
        xlabel(statNames{j})
        title(sprintf('p = %.4g', pvalue(j)))
        box on
    end
end

end

% =======================================================================
function [d2_all_cc, muHat, SigHat] = local_fit_and_get_complete_distances(Y, completeIdx, alpha, method, tol)

p = size(Y,2); 

if alpha == 0
    outFit = mdEM(Y);
else
    outFit = mdTEM(Y, 'method', method, 'alpha', alpha, 'tol', tol);
end

muHat = outFit.loc;
SigHat = outFit.cov;

[d2_part, poss] = mdPartialMD(Y, muHat, SigHat);
d2_full = mdPartialMD2full(d2_part, p, poss, 'method', method);
d2_all_cc = d2_full(completeIdx);

end

% =======================================================================
function T = local_statistic(d2_cc, d2_all)
% Compute the 4 statistics
eps0 = 1e-12;
d2_cc  = max(d2_cc,  eps0);
d2_all = max(d2_all, eps0);
rat=log(d2_all ./ d2_cc);
T1 = median(rat);
T2 = mean(rat);
dif=d2_all - d2_cc;
T3 = median(dif);
T4 = mean(dif);
T=[T1 T2 T3 T4];
end

% =======================================================================
function Sspd = local_make_spd(S)

S = (S + S')/2;
[~,flag] = chol(S);

if flag == 0
    Sspd = S;
    return
end

% Simple ridge regularization if needed
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

error('Unable to regularize covariance matrix to positive definiteness.');
end