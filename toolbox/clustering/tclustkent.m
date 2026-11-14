function [out] = tclustkent(Y, k, alpha, restrfactor, varargin)
%tclustkent computes trimmed clustering for spherical data using mixtures of Kent distributions
%
%<a href="matlab: docsearchFS('tclustkent')">Link to the help page for this function</a>
%
%   tclustkent partitions the points in the n-by-3 data matrix Y, where
%   each row is assumed to be a unit vector on the unit sphere S^2 in R^3,
%   into k clusters modelled by Kent (Fisher-Bingham FB5) distributions.
%   It is the Kent-distribution counterpart of tclustvmf, in exactly the
%   same way that tclustvmf is the vMF counterpart of tclust.
%
%   The Kent distribution Kent(kappa,beta,Gamma) on S^2 has density
%
%     f(x; kappa,beta,Gamma) = (1/c(kappa,beta)) *
%         exp{ kappa*(gamma3'*x) + beta*[(gamma1'*x)^2 - (gamma2'*x)^2] }
%
%   where:
%     gamma3  is the mean direction (pole), a unit vector in R^3,
%     gamma1  is the major-axis direction (orthogonal to gamma3),
%     gamma2  is the minor-axis direction (= cross(gamma3,gamma1)),
%     kappa>=0 controls concentration (the vMF parameter),
%     beta>=0  controls ovalness of contours (0 <= beta < kappa/2),
%     c(kappa,beta) is the normalizing constant (see KentlogC).
%
%   beta=0 reduces the Kent density to the vMF density with parameter
%   kappa, and setting restrfactor=1 with beta=0 reduces tclustkent to
%   trimmed spherical k-means (like tclustvmf with restrfactor=1).
%
%   The Kent distribution is STRICTLY for S^2 (p=3 ambient dimension).
%   For higher-dimensional spherical data use tclustvmf instead.
%
%   Robustness is achieved via the same impartial trimming device used by
%   tclust and tclustvmf: at most n-h observations are left unassigned,
%   where h = fix(n*(1-alpha)). A constraint on max(kappa_j)/min(kappa_j)
%   prevents the classification likelihood from becoming unbounded as a
%   cluster collapses to a point (the vMF/Kent analogue of tclust's
%   eigenvalue restriction). An optional uniform noise component models
%   background contamination uniformly distributed on S^2.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-3 matrix of unit vectors on S^2.
%               Use option 'normalize' to rescale automatically.
%            k: Number of groups. Scalar. Positive integer.
%        alpha: Global trimming level. Scalar. 0<=alpha<1 (fraction) or
%               integer >=1 (count). alpha=0 means no trimming.
%  restrfactor: Restriction factor. Scalar >= 1. Constrains the ratio
%               max(kappa_j)/min(kappa_j). restrfactor=1 forces a common
%               kappa. See restrKentFS.
%
%  Optional input arguments:
%
%    equalweights : Cluster weights in the assignment. Logical. If false
%               (default) the mixing proportions pi_j are estimated and
%               enter the discriminant. If true all pi_j are fixed to 1/k
%               (or 1/kn with noise).
%                 Example - 'equalweights',true
%                 Data Types - logical
%       mixt  : Mixture modelling or crisp assignment. Scalar. mixt=0
%               (default) is the crisp CEM algorithm. mixt=2 is the EM
%               algorithm (soft trimming). mixt=1 is a hybrid (crisp
%               trimming, soft M-step).
%                 Example - 'mixt',2
%                 Data Types - double
%       nsamp : Number of subsamples. Scalar. Number of subsamples of
%               observations which are extracted to initialise the
%               algorithm. The default is min(300, nchoosek). Each
%               subsample contains 4 observations per cluster, which is
%               the minimum needed to seed a non-degenerate Kent
%               estimate in 3D.
%                 Example - 'nsamp',500
%                 Data Types - double
%     refsteps: Number of refining iterations. Scalar. Number of
%               concentration steps applied to each subsample. The
%               default is 15.
%                 Example - 'refsteps',10
%                 Data Types - double
%      reftol : Tolerance of the refining iterations. Scalar. Convergence
%               tolerance of the concentration steps. The default is
%               1e-06.
%                 Example - 'reftol',1e-08
%                 Data Types - double
%    normalize: Rescaling of the rows. Logical. If true each row of Y is
%               rescaled to unit norm. If false (default) an error is
%               thrown when the rows are not unit-norm.
%                 Example - 'normalize',true
%                 Data Types - logical
%     conflev : Confidence level of the confidence regions. Scalar. Value
%               in (0,1), the default is 0.95. The confidence region of a
%               group is the region of S^2 delimited by the contour of
%               the fitted Kent density which contains probability
%               conflev (an oval around the pole, elongated along the
%               major axis).
%                 Example - 'conflev',0.90
%                 Data Types - double
%  ellipsenpts: Points of the confidence regions. Scalar. Number of
%               points used to represent the boundary of each confidence
%               region. The default is 200. The points are ordered along
%               the contour (see More About).
%                 Example - 'ellipsenpts',300
%                 Data Types - double
%  markersize : Marker size. Scalar. Marker size for the scatter plots.
%               The default is 12.
%                 Example - 'markersize',20
%                 Data Types - double
%        msg  : Level of display. Scalar. The default is 1 (basic
%               progress). Use 0 for silent execution and 2 for more
%               detail.
%                 Example - 'msg',0
%                 Data Types - double
%      Ysave  : Save of the data. Logical. If true Y is stored in out.Y.
%               The default is false.
%                 Example - 'Ysave',true
%                 Data Types - logical
%      plots  : Plot on the sphere. Scalar. If 1 the classification is
%               shown on the unit sphere with function sphereplot: units
%               coloured by group, trimmed units as red crosses, units of
%               the noise component in gray, and (see plottype) the
%               confidence regions of the groups. The sphere can be
%               rotated with the mouse. The default is 0. Call sphereplot
%               directly, using out.idx, out.ellipse and out.gamma3, to
%               personalise the plot.
%                 Example - 'plots',1
%                 Data Types - double
%    plottype : Content of the plot. Character. 'points' (default) shows
%               the units only, 'ellipse' the confidence regions only and
%               'both' the units and the confidence regions.
%                 Example - 'plottype','both'
%                 Data Types - char
%      nocheck: Check of the data. Scalar. If 1 the check for NaN/Inf is
%               skipped. The default is 0.
%                 Example - 'nocheck',1
%                 Data Types - double
%       noise : Uniform noise component. Logical. If true an extra
%               (k+1)-th uniform-on-S^2 component is added to model
%               background noise. The default is false. The uniform
%               density on S^2 is 1/(4*pi), fixed (no free parameters).
%               Its mixing proportion is estimated jointly with the k
%               cluster proportions.
%                 Example - 'noise',true
%                 Data Types - logical
%  mstepiter : Iterations of the M-step. Scalar. Maximum number of
%               iterations of the numerical optimiser inside KentMstepFS
%               (kappa/beta estimation per M-step). The default is 50.
%               Reduce for speed, increase for accuracy in
%               high-concentration cases.
%                 Example - 'mstepiter',100
%                 Data Types - double
%
%  Output:
%
%  out :     A structure with the following fields:
%
%         out.idx  = n-by-1 vector of cluster assignments. Trimmed units
%                    have idx=0. Noise units have idx=k+1 (if noise=true).
%      out.gamma1  = k-by-3 matrix: estimated major-axis directions,
%                    one row per group.
%      out.gamma2  = k-by-3 matrix: estimated minor-axis directions.
%      out.gamma3  = k-by-3 matrix: estimated mean directions (poles).
%    out.kappaopt  = k-by-1 vector of estimated (restricted) kappa values.
%     out.betaopt  = k-by-1 vector of estimated beta values.
%    out.postprob  = n-by-kn posterior probability matrix (kn=k or k+1).
%     out.ellipse  = k-by-1 cell array; out.ellipse{j} is an
%                    ellipsenpts-by-3 matrix of points on the boundary of
%                    the conflev Kent confidence region of group j,
%                    ordered along the (closed) contour, so that it can
%                    be drawn with plot3.
%      out.conflev = confidence level used for out.ellipse.
%         out.siz  = tabulate(out.idx).
%           out.h  = number of untrimmed observations.
%         out.obj  = maximised trimmed log-likelihood.
%       out.NlogL  = -2 times classification log-likelihood of retained units.
%          out.bs  = elemental subset that gave the best solution.
%     out.fullsol  = nselected-by-1 vector of objective values per subsample.
%    out.notconver = fraction of subsamples without convergence.
% out.equalweights = logical, as input option.
%         out.mixt = scalar, as input option.
%        out.noise = logical, as input option.
%    out.noiseprop = (only if noise=true) estimated noise mixing proportion.
%      out.CLACLA  = (only if mixt==0) classification BIC: out.NlogL + nParam*log(h).
%      out.MIXMIX  = (only if mixt>0) mixture BIC: out.NlogLmixt + nParam*log(h).
%      out.MIXCLA  = (only if mixt>0) classification BIC at mixture estimates.
%            out.Y = (only if Ysave=true) the (possibly row-normalised) data.
%
% More About:
%
% CONCENTRATION STEPS (CEM/EM ALGORITHM ON S^2)
%
% Each subsample runs the following loop:
%  (i)  E-step: compute log-discriminants
%         d_{ij} = log(pi_j) - log c(kappa_j,beta_j)
%                  + kappa_j*(gamma3_j'*y_i)
%                  + beta_j*[(gamma1_j'*y_i)^2 - (gamma2_j'*y_i)^2]
%       for each unit i and group j; add the uniform noise term if noise=true.
%  (ii) Trimming: retain the h units with the largest discriminant (or
%       mixture log-density for mixt=2).
%  (iii)M-step: re-estimate (kappa_j, beta_j, Gamma_j) for each group
%       via KentMstepFS (closed-form Gamma, numerical kappa/beta).
%  (iv) Restriction: apply restrKentFS to enforce
%       max(kappa_j)/min(kappa_j) <= restrfactor.
%
% SUBSAMPLE SEEDING
%
% Each elemental subset draws k*4 observations (4 per cluster). With 4
% points on S^2 we can form a non-degenerate mean vector and second-moment
% matrix, so KentMstepFS can produce a valid seed. Using fewer points
% risks degenerate second-moment matrices with coincident eigenvalues.
%
% PARAMETER COUNT AND INFORMATION CRITERIA
%
% A Kent distribution on S^2 has 5 free parameters: kappa (1), beta (1),
% and 3 rotation angles for Gamma (5 total). The uniform noise component
% is parameter-free (density = 1/(4*pi) everywhere on S^2). With kn total
% components and estimated mixing proportions:
%
%   nParam = k*5 + (kn-1)   [equalweights=false]
%          = k*5             [equalweights=true]
%
% where kn=k (noise=false) or kn=k+1 (noise=true). The three information
% criteria are identical in structure to those of tclustvmf:
%   CLACLA = -2*CLA_CLA + nParam*log(h)   [mixt=0]
%   MIXMIX = -2*MIX_MIX + nParam*log(h)   [mixt>0]
%   MIXCLA = -2*MIX_CLA + nParam*log(h)   [mixt>0]
%
% KENT CONFIDENCE REGION (out.ellipse)
%
% The vMF confidence region is a spherical cap (a circle around the mean
% direction). The contours of the Kent density are instead ovals: the
% confidence region of group j is the set of points whose fitted density
% is larger than a threshold t_j, chosen in such a way that the region
% has probability conflev under Kent(kappa_j,beta_j,Gamma_j). In the
% canonical frame, with colatitude theta and longitude phi measured from
% the major axis, the log-density (up to the normalizing constant) is
%
%   kappa*cos(theta) + beta*sin(theta)^2*cos(2*phi).
%
% The boundary of the region is computed by function Kentcap, which
% can also be called directly. The angular semi-axes of the region are
% approximately proportional to 1/sqrt(kappa-2*beta) (major) and
% 1/sqrt(kappa+2*beta) (minor).
%
% RESTRICTION
%
% The restriction acts on the concentrations kappa_j. When restrKentFS
% modifies kappa_j, beta_j is rescaled so that beta_j/kappa_j (the shape
% of the contours of group j) does not change.
%
% References:
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
%   Mardia, K.V., Barber, S., Burdett, P.M., Kent, J.T. and Hamelryck, T.
%   (2021), Mixture models for spherical data with applications to protein
%   bioinformatics, arXiv:2104.13140.
%
%   Garcia-Escudero, L.A., Gordaliza, A., Matran, C. and Mayo-Iscar, A.
%   (2008), A General Trimming Approach to Robust Cluster Analysis,
%   Annals of Statistics, Vol. 36, pp. 1324-1345.
%
%   Cerioli, A., Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani, M.
%   (2017), Finding the Number of Groups in Model-Based Clustering via
%   Constrained Likelihoods, J. Computational and Graphical Statistics,
%   Vol. 27, pp. 404-416.
%
%   Dong, A. and Melnykov, V. (2024), Contaminated Kent mixture model for
%   clustering non-spherical directional data with heavy tails or scatter,
%   Statistics and Probability Letters, Vol. 208, 110058.
%
% See also: tclustvmf, tclust, sphereplot, Kentrnd, Kentcap, logKentpdfFS, KentMstepFS, restrKentFS, KentlogC
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustkent')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% tclustkent on two elliptical groups with uniform contamination.
    % The two groups are generated from Kent distributions (elliptical
    % caps): the first has pole [1 0 0] and is elongated along [0 0 1],
    % the second has pole [0 1 0] and is elongated along [1 0 0]. 40
    % units uniformly distributed on the sphere are added. Trimming level
    % is 10 per cent. The 95 per cent confidence regions are drawn.
    rng(1)
    n1=200; n2=200; nnoise=40;
    Y1 = Kentrnd(n1, 100, 40, [1 0 0], [0 0 1]);
    Y2 = Kentrnd(n2,  50, 15, [0 1 0], [1 0 0]);
    % uniform noise on the sphere (vMF with kappa=0)
    Ynoise = vMFrnd(nnoise, [0 0 1], 0);
    Y = [Y1; Y2; Ynoise];
    out = tclustkent(Y, 2, 0.1, 20, 'plots', 1, 'plottype', 'both');
    % Estimated parameters (true values: kappa 100 and 50, beta 40 and 15)
    disp(array2table([out.kappaopt out.betaopt], 'VariableNames', {'kappa','beta'}))
    disp('Estimated poles'); disp(out.gamma3)
    disp('Estimated major axes (defined up to sign)'); disp(out.gamma1)
%}
%
%{
    % The noise option: uniform background component.
    % Same two elliptical groups as above, with 60 units uniformly
    % distributed on the sphere. An extra component with uniform density
    % is added and a small trimming level (2 per cent) is used: most of
    % the contamination is captured by the noise component (out.idx==3,
    % in gray) rather than by trimming (out.idx==0, in red).
    rng(2)
    n1=200; n2=200; nnoise=60;
    Y1 = Kentrnd(n1, 100, 40, [1 0 0], [0 0 1]);
    Y2 = Kentrnd(n2,  50, 15, [0 1 0], [1 0 0]);
    % uniform noise on the sphere (vMF with kappa=0)
    Ynoise = vMFrnd(nnoise, [0 0 1], 0);
    Y = [Y1; Y2; Ynoise];
    out = tclustkent(Y, 2, 0.02, 20, 'noise', true, 'plots', 1, 'plottype', 'both');
    fprintf('Estimated noise proportion: %.1f%%\n', 100*out.noiseprop);
%}
%
%{
    % Comparison with tclustvmf when the groups are elongated.
    % Two strongly elongated groups, parallel to each other, whose poles
    % are only 20 degrees apart. The groups are well separated along
    % their minor axis but each of them is much longer than the distance
    % between the two. The vMF model can only describe circular caps and
    % cannot recover this structure, while the Kent model does.
    rng(3)
    n = 200;
    pole2 = [cosd(20) sind(20) 0];
    Y = [Kentrnd(n, 200, 90, [1 0 0], [0 0 1]);
         Kentrnd(n, 200, 90, pole2,   [0 0 1])];
    truelab = [ones(n,1); 2*ones(n,1)];
    outK = tclustkent(Y, 2, 0, 20, 'msg', 0);
    outV = tclustvmf(Y, 2, 0, 20, 'msg', 0);
    figure
    sphereplot(Y, 'group', outK.idx, 'cap', outK.ellipse, 'center', outK.gamma3, ...
        'h', subplot(1,2,1));
    title('tclustkent')
    sphereplot(Y, 'group', outV.idx, 'h', subplot(1,2,2), ...
        'cap', struct('mu', outV.muopt, 'kappa', outV.kappaopt));
    title('tclustvmf')
    fprintf('Adjusted Rand index: tclustkent %.3f, tclustvmf %.3f\n', ...
        RandIndexFS(truelab, outK.idx), RandIndexFS(truelab, outV.idx));
%}
%
%{
    % Effect of the restriction factor.
    % A very concentrated and a very dispersed elliptical group (true
    % kappa 200 and 10). With restrfactor=1 the two groups are forced to
    % have the same kappa; with restrfactor=50 the two concentrations are
    % free to differ. In both cases the ratios beta/kappa, i.e. the
    % shapes of the groups, are estimated separately for each group.
    rng(4)
    n = 300;
    Y = [Kentrnd(n, 200, 80, [1 0 0], [0 0 1]);
         Kentrnd(n,  10,  3, [0 1 0], [0 0 1])];
    out1  = tclustkent(Y, 2, 0.05,  1, 'msg', 0);
    out50 = tclustkent(Y, 2, 0.05, 50, 'msg', 0, 'plots', 1, 'plottype', 'both');
    disp('restrfactor=1:  kappa and beta');  disp([out1.kappaopt out1.betaopt])
    disp('restrfactor=50: kappa and beta');  disp([out50.kappaopt out50.betaopt])
%}
%
%{
    % Three elliptical groups with different orientation.
    % Posterior probabilities and group sizes.
    rng(5)
    n = 150;
    Y = [Kentrnd(n, 120, 50, [1 0 0], [0 1 1]);
         Kentrnd(n,  80, 30, [0 1 0], [1 0 0]);
         Kentrnd(n,  60, 10, [0 0 1], [1 1 0])];
    out = tclustkent(Y, 3, 0.05, 20, 'plots', 1, 'plottype', 'both');
    disp(out.siz)
%}
%
%{
    % The mixt option: classification versus mixture likelihood.
    % mixt=0 (default) is the crisp CEM algorithm; mixt=2 uses the
    % mixture likelihood both for trimming and for the M-step.
    rng(6)
    n = 200;
    Y = [Kentrnd(n, 100, 40, [1 0 0], [0 0 1]);
         Kentrnd(n,  50, 15, [0 1 0], [1 0 0])];
    out0 = tclustkent(Y, 2, 0.05, 20, 'mixt', 0, 'msg', 0);
    out2 = tclustkent(Y, 2, 0.05, 20, 'mixt', 2, 'msg', 0);
    fprintf('mixt=0: obj=%.2f, CLACLA=%.2f\n', out0.obj, out0.CLACLA);
    fprintf('mixt=2: obj=%.2f, MIXMIX=%.2f, MIXCLA=%.2f\n', out2.obj, out2.MIXMIX, out2.MIXCLA);
%}
%
%{
    % The confidence regions stored in out.ellipse.
    % out.ellipse{j} contains the points of the boundary of the region,
    % ordered along the contour. Here the regions at level 0.5 are
    % requested and the proportion of units of each group which fall
    % inside its region is computed.
    rng(7)
    n = 1000;
    Y = [Kentrnd(n, 100, 40, [1 0 0], [0 0 1]);
         Kentrnd(n,  50, 15, [0 1 0], [1 0 0])];
    out = tclustkent(Y, 2, 0, 20, 'conflev', 0.5, 'msg', 0, 'plots', 1, 'plottype', 'both');
    for j = 1:2
        % threshold of the log-density = its value on the boundary
        thr = mean(logKentpdfFS(out.ellipse{j}, out.gamma1(j,:), out.gamma2(j,:), ...
            out.gamma3(j,:), out.kappaopt(j), out.betaopt(j)));
        Yj = Y(out.idx==j,:);
        inside = logKentpdfFS(Yj, out.gamma1(j,:), out.gamma2(j,:), ...
            out.gamma3(j,:), out.kappaopt(j), out.betaopt(j)) >= thr;
        fprintf('Group %d: %.1f%% of the units inside the 50%% region\n', j, 100*mean(inside));
    end
%}

%% Beginning of code

nnargin = nargin;
vvarargin = varargin;
Y = aux.chkinputM(Y, nnargin, vvarargin);
[n, p] = size(Y);

if p ~= 3
    error('FSDA:tclustkent:WrongDim', ...
        ['tclustkent requires 3-column data (Kent distribution is defined on S^2 only).\n' ...
         'For higher-dimensional spherical data use tclustvmf instead.'])
end

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
    warning('FSDA:tclustkent:Wrongalpha', ...
        'You have not specified alpha: it is set to 0.05 by default');
end

if nargin < 4 || isempty(restrfactor)
    restrfactor = 20;
    warning('FSDA:tclustkent:Wrongrestrfact', ...
        'You have not specified restrfactor: it is set to 20 by default');
end

if restrfactor < 1
    error('FSDA:tclustkent:WrongRestrfactor', 'restrfactor must be >= 1')
end

if alpha < 0
    error('FSDA:tclustkent:WrongAlpha', 'alpha must be in [0,0.5) or a positive integer')
end
if alpha >= 1
    h = n - floor(alpha);
else
    h = fix(n * (1 - alpha));
end

% Elemental subset size: 4 observations per group (minimum for a
% non-degenerate Kent seed in 3D: need at least 3 non-collinear points
% for the mean vector and second-moment matrix, we use 4 for safety).
elemsize = k * 4;
ncomb = bc(n, elemsize);
nsampdef = min(300, ncomb);
refstepsdef = 15;
reftoldef = 1e-06;

options = struct('nsamp', nsampdef, 'plots', 0, 'plottype', 'points', ...
    'nocheck', 0, 'msg', 1, 'Ysave', 0, ...
    'refsteps', refstepsdef, 'reftol', reftoldef, 'equalweights', false, ...
    'mixt', 0, 'normalize', false, 'conflev', 0.95, 'ellipsenpts', 200, ...
    'markersize', 12, 'noise', false, 'mstepiter', 50);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions = varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:tclustkent:WrongInputOpt', ...
            'Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    inpchk = isfield(options, UserOptions);
    WrongOptions = UserOptions(inpchk == 0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:tclustkent:NonExistInputOpt', ...
            'In total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin > 4
    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
    if options.nsamp > ncomb
        disp('Number of subsets to extract greater than nchoosek. It is set to that value.');
        options.nsamp = 0;
    elseif options.nsamp < 0
        error('FSDA:tclustkent:WrongNsamp', ...
            'Number of subsets to extract must be 0 (all) or a positive number');
    end
end

nsamp       = options.nsamp;
refsteps    = options.refsteps;
reftol      = options.reftol;
equalweights = options.equalweights;
mixt        = options.mixt;
msg         = options.msg;
plots       = options.plots;
plottype    = options.plottype;
normalizeOpt = options.normalize;
conflev     = options.conflev;
ellipsenpts = options.ellipsenpts;
markersize  = options.markersize;
noise       = options.noise;
mstepiter   = options.mstepiter;

kn = k + double(noise);
lognoiseconst = -log(4*pi);   % log(1/(4*pi)) = uniform density on S^2

if conflev <= 0 || conflev >= 1
    error('FSDA:tclustkent:WrongConflev', 'conflev must be in (0,1)')
end
if ~any(strcmp(plottype, {'points','ellipse','both'}))
    error('FSDA:tclustkent:WrongPlottype', ...
        'plottype must be one of ''points'', ''ellipse'', ''both''')
end

% Check / enforce unit-norm rows
rownorms = sqrt(sum(Y.^2, 2));
if normalizeOpt
    Y = Y ./ rownorms;
else
    if max(abs(rownorms - 1)) > 1e-06
        error('FSDA:tclustkent:NotUnitNorm', ...
            'Rows of Y must be unit vectors (use the ''normalize'' option to rescale automatically)')
    end
end

tolrestr = 1e-08;

if coder.target('MATLAB')
    verLess2016b = verLessThanFS('9.1');
else
    verLess2016b = true;
end

%% Combinatorial part
[C, nselected] = subsets(nsamp, n, elemsize, ncomb, msg);

vopt = -1e+30;
fullsol = zeros(nselected, 1);
noconv = 0;

% Preallocate optimal parameter containers
gamma1opt = NaN(k, 3);
gamma2opt = NaN(k, 3);
gamma3opt = NaN(k, 3);
kappaopt  = NaN(k, 1);
betaopt   = NaN(k, 1);
nopt      = NaN(kn, 1);
bs        = [];

for i = 1:nselected

    index = C(i, :);

    niini   = zeros(kn, 1);
    g1ini   = zeros(k, 3);
    g2ini   = zeros(k, 3);
    g3ini   = zeros(k, 3);
    kapini  = zeros(k, 1);
    betaini = zeros(k, 1);
    % sufficient statistics of (kappa,beta), needed by restrKentFS
    rini    = zeros(k, 1);
    dqini   = zeros(k, 1);

    % --- Seed each group from an elemental subset of 4 points ---
    for j = 1:k
        ilow = (j-1)*4 + 1;
        iup  = j*4;
        selj = index(ilow:iup);
        Yselj = Y(selj, :);

        % Seed from elemental subset via KentMstepFS with unit weights
        wj = ones(4, 1);
        [g1j, g2j, g3j, kapj, betj, rj, dqj] = KentMstepFS(Yselj, wj, [], [], mstepiter);

        if isnan(kapj) || kapj <= 0
            % Degenerate seed: use the resultant as pole, set beta=0
            Rj = sum(Yselj, 1);
            nRj = norm(Rj);
            if nRj > 1e-12
                g3j = Rj / nRj;
            else
                g3j = Yselj(1, :);
            end
            % Build an arbitrary orthonormal frame
            [~, imin] = min(abs(g3j));
            evec = zeros(1, 3); evec(imin) = 1;
            g1j = evec - (evec * g3j') * g3j; g1j = g1j / norm(g1j);
            g2j = cross(g3j, g1j);
            kapj = 1; betj = 0;
            rj = nRj/4; dqj = 0;
        end

        niini(j) = 4;
        g1ini(j, :) = g1j;
        g2ini(j, :) = g2j;
        g3ini(j, :) = g3j;
        kapini(j)  = kapj;
        betaini(j) = betj;
        rini(j)    = rj;
        dqini(j)   = dqj;
    end

    if noise
        niini(k+1) = 4;
    end

    % Apply initial restriction
    [kapini, betaini] = restrKentFS(kapini, betaini, rini, dqini, niini(1:k), restrfactor, tolrestr);

    iter = 0;
    mudiff = 1e+15;
    ind = zeros(n, 1);
    postprob = zeros(n, kn);
    obj = 0;
    ll = zeros(n, kn);
    groupind = [];
    qq = [];
    qqunassigned = [];

    while (mudiff > reftol) && (iter < refsteps)
        iter = iter + 1;

        % E-step: compute log-discriminants
        for j = 1:k
            if isnan(kapini(j))
                ll(:, j) = -Inf;
            else
                lj = logKentpdfFS(Y, g1ini(j,:), g2ini(j,:), g3ini(j,:), kapini(j), betaini(j));
                if equalweights
                    ll(:, j) = lj;
                else
                    ll(:, j) = log(niini(j)/h) + lj;
                end
            end
        end
        if noise
            if equalweights
                ll(:, k+1) = lognoiseconst;
            else
                ll(:, k+1) = log(niini(k+1)/h) + lognoiseconst;
            end
        end

        if mixt == 2
            postprobold = postprob;
            [~, postprob, disc] = estepFS(ll, verLess2016b);
            [~, qq] = sort(disc, 'descend');
            qqunassigned = qq((h+1):n);
            qq = qq(1:h);
            postprob(qqunassigned, :) = 0;
            niini = sum(postprob)';
        else
            indold = ind;
            [disc, ind] = max(ll, [], 2);
            [~, qq] = sort(disc, 'descend');
            qqunassigned = qq((h+1):n);
            qq = qq(1:h);
            groupind = ind(qq);
            ind(qqunassigned) = 0;
        end

        if mixt == 1
            postprobold = postprob;
            [~, postprob] = estepFS(ll, verLess2016b);
            postprob(qqunassigned, :) = 0;
            niini = sum(postprob)';
        end

        % M-step: re-estimate Kent parameters for each group
        for j = 1:k
            if mixt >= 1
                wj = postprob(:, j);
                nj = niini(j);
            else
                niini(j) = sum(groupind == j);
                nj = niini(j);
                wj = double(ind == j);
                wj(qqunassigned) = 0;
            end

            if nj > 0
                [g1j, g2j, g3j, kapj, betj, rj, dqj] = KentMstepFS(Y, wj, kapini(j), betaini(j), mstepiter);
                if isnan(kapj)
                    g1ini(j,:) = NaN(1,3); g2ini(j,:) = NaN(1,3); g3ini(j,:) = NaN(1,3);
                    kapini(j) = NaN; betaini(j) = NaN; niini(j) = 0;
                else
                    g1ini(j,:) = g1j; g2ini(j,:) = g2j; g3ini(j,:) = g3j;
                    kapini(j) = kapj; betaini(j) = betj;
                    rini(j) = rj; dqini(j) = dqj;
                end
            else
                g1ini(j,:) = NaN(1,3); g2ini(j,:) = NaN(1,3); g3ini(j,:) = NaN(1,3);
                kapini(j) = NaN; betaini(j) = NaN; niini(j) = 0;
            end
        end

        if noise && mixt == 0
            niini(k+1) = sum(groupind == (k+1));
        end

        % Apply concentration restriction
        [kapini, betaini] = restrKentFS(kapini, betaini, rini, dqini, niini(1:k), restrfactor, tolrestr);

        % Compute objective function
        Ytri = Y(qq, :);
        obj = 0;
        if mixt >= 1
            log_lh = zeros(h, kn);
            for j = 1:k
                if ~isnan(kapini(j))
                    lj = logKentpdfFS(Ytri, g1ini(j,:), g2ini(j,:), g3ini(j,:), kapini(j), betaini(j));
                    log_lh(:, j) = log(niini(j)/h) + lj;
                else
                    log_lh(:, j) = -Inf;
                end
            end
            if noise
                if equalweights
                    log_lh(:, k+1) = lognoiseconst;
                else
                    log_lh(:, k+1) = log(niini(k+1)/h) + lognoiseconst;
                end
            end
            obj = estepFS(log_lh, verLess2016b);
        else
            for j = 1:k
                nj = niini(j);
                if nj > 0 && ~isnan(kapini(j))
                    Ytri_j = Y(qq(groupind == j), :);
                    lj = logKentpdfFS(Ytri_j, g1ini(j,:), g2ini(j,:), g3ini(j,:), kapini(j), betaini(j));
                    if equalweights
                        obj = obj + sum(lj);
                    else
                        obj = obj + nj*log(nj/h) + sum(lj);
                    end
                end
            end
            if noise && niini(k+1) > 0
                if equalweights
                    obj = obj + niini(k+1) * lognoiseconst;
                else
                    obj = obj + niini(k+1)*log(niini(k+1)/h) + niini(k+1)*lognoiseconst;
                end
            end
        end

        % Convergence check
        if mixt > 0
            mudiff = sum(sum(abs(postprob - postprobold))) / n;
        else
            mudiff = sum(abs(indold - ind) > 0) / n;
        end

        if iter == refsteps
            noconv = noconv + 1;
        end
    end % while loop

    fullsol(i) = obj;

    if obj >= vopt && ~any(isnan(kapini))
        vopt = obj;
        gamma1opt = g1ini;
        gamma2opt = g2ini;
        gamma3opt = g3ini;
        kappaopt  = kapini;
        betaopt   = betaini;
        nopt      = niini;
        bs        = index;
    end

end % subsample loop

if isempty(bs) || any(isnan(kappaopt))
    error('FSDA:tclustkent:NoConvergence', ...
        ['None of the %d subsamples produced a valid solution; ' ...
         'try increasing nsamp or refsteps, or reducing k'], nselected)
end

if nselected > 0 && noconv/nselected > 0.1
    disp('------------------------------')
    disp(['Warning: Number of subsets without convergence equal to ' ...
        num2str(100*noconv/nselected) '%'])
end

%% Final assignment with optimal parameters
ll = zeros(n, kn);
if equalweights
    for j = 1:k
        ll(:, j) = logKentpdfFS(Y, gamma1opt(j,:), gamma2opt(j,:), gamma3opt(j,:), kappaopt(j), betaopt(j));
    end
else
    for j = 1:k
        ll(:, j) = log(nopt(j)/h) + logKentpdfFS(Y, gamma1opt(j,:), gamma2opt(j,:), gamma3opt(j,:), kappaopt(j), betaopt(j));
    end
end
if noise
    if equalweights
        ll(:, k+1) = lognoiseconst;
    else
        ll(:, k+1) = log(nopt(k+1)/h) + lognoiseconst;
    end
end

[~, postprob, logpdf] = estepFS(ll, verLess2016b);
[disc, idx] = max(ll, [], 2);
[~, qq] = sort(disc, 'descend');

if mixt >= 1
    [~, qqmixt] = sort(logpdf, 'descend');
    unassignedmixt = qqmixt((h+1):n);
    assignedmixt = qqmixt(1:h);
    [~, idxmixt] = max(postprob, [], 2);
    idxmixt(unassignedmixt) = 0;
    postprob(unassignedmixt, :) = 0;
    out.idx = idxmixt;
    NlogLmixt = -estepFS(ll(assignedmixt,:), verLess2016b);
    out.NlogLmixt = 2*NlogLmixt;
else
    unassigned = qq((h+1):n);
    idx(unassigned) = 0;
    postprob(unassigned, :) = 0;
    out.idx = idx;
end

loglik = disc(qq(1:h));
NlogL  = -sum(loglik);

% Parameter count: k*5 (each Kent component has 5 free parameters:
% kappa, beta, and 3 rotation angles for Gamma on SO(3))
% plus kn-1 mixing proportion parameters when equalweights=false.
nParam = k * 5;
if ~equalweights
    nParam = nParam + (kn - 1);
end
logh = log(h);

if mixt > 0
    MIXMIX = 2*NlogLmixt + nParam*logh;
    MIXCLA = 2*NlogL + nParam*logh;
else
    CLACLA = 2*NlogL + nParam*logh;
end

out.gamma1   = gamma1opt;
out.gamma2   = gamma2opt;
out.gamma3   = gamma3opt;
out.kappaopt = kappaopt;
out.betaopt  = betaopt;
out.postprob = postprob;

% Boundary of the conflev confidence region of each group (see Kentcap)
out.ellipse = cell(k, 1);
for j = 1:k
    out.ellipse{j} = Kentcap(gamma1opt(j,:), gamma2opt(j,:), gamma3opt(j,:), ...
        kappaopt(j), betaopt(j), conflev, ellipsenpts);
end
out.conflev     = conflev;
out.siz         = tabulate(out.idx);
out.h           = h;
out.obj         = vopt;
out.NlogL       = 2*NlogL;
out.equalweights = equalweights;
out.mixt        = mixt;
out.noise       = noise;
if noise
    out.noiseprop = nopt(k+1)/h;
end
out.bs       = bs;
out.fullsol  = fullsol;
out.notconver = noconv / max(nselected, 1);
if mixt > 0
    out.MIXMIX = MIXMIX;
    out.MIXCLA = MIXCLA;
else
    out.CLACLA = CLACLA;
end
if options.Ysave
    out.Y = Y;
end

%% Plots
if plots == 1
    % The plot is produced by sphereplot
    if strcmp(plottype, 'points')
        capplot = [];
    else
        capplot = out.ellipse;
    end
    if noise
        noiseid = k+1;
        noisestr = sprintf(', noise prop.=%.1f%% (gray)', 100*out.noiseprop);
    else
        noiseid = [];
        noisestr = '';
    end
    hplot = sphereplot(Y, 'group', out.idx, 'noisegroup', noiseid, ...
        'cap', capplot, 'center', gamma3opt, 'markersize', markersize, ...
        'showpoints', ~strcmp(plottype, 'ellipse'));
    ngroups = numel(setdiff(unique(out.idx(:))', [0, k+1]));
    title(hplot, sprintf('%d groups found by tclustkent for alpha=%.2f, restrfactor=%.0f%s', ...
        ngroups, alpha, restrfactor, noisestr));
end

end
%FScategory:CLUS-RobClaMULT
