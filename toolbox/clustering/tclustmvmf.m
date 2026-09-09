function [out] = tclustmvmf(Y,k,alpha,restrfactor,varargin)
%tclustmvmf computes trimmed clustering for directional data using mixtures of von Mises-Fisher distributions
%
%<a href="matlab: docsearchFS('tclustmvmf')">Link to the help function</a>
%
%   tclustmvmf partitions the points in the n-by-p data matrix Y, where
%   each row is assumed to be a unit vector on the hypersphere S^(p-1) in
%   R^p, into k clusters. This is the directional-data counterpart of
%   tclust: the multivariate normal model N(mu_j,Sigma_j) is replaced by
%   the von Mises-Fisher model vMF(mu_j,kappa_j), where mu_j (unit
%   vector) is the mean direction and kappa_j >=0 is the concentration
%   parameter of group j (the vMF analogue of the inverse of Sigma_j: a
%   large kappa_j implies a tightly concentrated, "low variance" cluster
%   around mu_j, exactly as a large eigenvalue of Sigma_j^{-1} does in
%   the Gaussian case). Robustness is achieved via the same impartial
%   ("hard") trimming device used by tclust and tkmeans: at most
%   n-h observations are left unassigned, and a restriction on the
%   ratio max(kappa_j)/min(kappa_j) prevents the classification
%   likelihood from becoming unbounded as a cluster collapses onto a
%   single direction (the vMF analogue of tclust's eigenvalue
%   restriction, which prevents Sigma_j from collapsing to a singular
%   matrix).
%
%   Setting restrfactor=1 forces a common kappa across all groups, which
%   reduces the discriminant rule to cosine similarity: this is
%   "trimmed spherical k-means", the exact directional analogue of the
%   relationship between tclust (restrfactor=1) and tkmeans.
%
%   Optionally (see the 'noise' argument), an extra (k+1)-th component
%   modelling uniform "background" noise on the hypersphere can be added
%   alongside the k vMF clusters, to be used together with (typically a
%   smaller) trimming level: trimming remains a distribution-free
%   safeguard against contamination of unknown/arbitrary shape, while
%   the noise component explicitly captures the specific case of
%   contamination that is itself uniform on the sphere.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-p matrix of directional data:
%               each row of Y must be a unit vector (norm 1). Use the
%               'normalize' option to have tclustmvmf rescale the rows
%               for you; otherwise an error is thrown if rows are not
%               (numerically) unit-norm.
%            k: Number of groups. Scalar.
%        alpha: Global trimming level. Scalar. As in tclust, if
%               0<=alpha<1 clustering is based on h=fix(n*(1-alpha))
%               observations, else if alpha is an integer >=1 clustering
%               is based on h=n-floor(alpha).
%  restrfactor: Restriction factor on the concentration parameters.
%               Scalar >=1. Constrains max(kappa_j)/min(kappa_j) across
%               the k groups. restrfactor=1 forces a common kappa
%               (trimmed spherical k-means); larger values allow groups
%               of increasingly different concentration (tightness).
%
%  Optional input arguments:
%
%    equalweights : Cluster weights in the concentration/assignment
%               steps and in the likelihood. Logical. Default false
%               (mixing proportions pi_j=n_j/h are estimated and enter
%               the discriminant, exactly as in tclust).
%                 Example - 'equalweights',true
%       mixt  : Mixture modelling or crisp assignment. Scalar. Same
%               semantics as in tclust: mixt=0 (default) crisp
%               classification (CEM); mixt=1 crisp trimming criterion
%               with posterior-probability-weighted M-step; mixt=2 fully
%               soft trimming criterion based on the mixture density.
%                 Example - 'mixt',0
%       nsamp : Number of subsamples to extract. Scalar. Default is
%               min(300, nchoosek(n,k*(p+1))).
%                 Example - 'nsamp',500
%     refsteps: Number of concentration steps per subsample. Scalar.
%               Default 15.
%      reftol : Tolerance for the concentration steps. Scalar. Default 1e-06.
%    normalize: Rescale each row of Y to unit norm before clustering.
%               Logical. Default false (an error is thrown instead if
%               rows are not already unit-norm).
%                 Example - 'normalize',true
%  newtonsteps: Number of Newton-Raphson refinement steps used by
%               kappaFS to invert the mean resultant length into a
%               concentration parameter. Scalar. Default 2 (see kappaFS
%               for accuracy figures).
%     conflev : Confidence level of the group confidence caps returned in
%               out.ellipse (see vMFcap.m) and, when plotted, drawn as
%               the boundary circles. Scalar in (0,1). Default 0.95.
%                 Example - 'conflev',0.90
% ellipsenpts : Number of points used to represent each confidence cap
%               boundary in out.ellipse (evenly spaced around the exact
%               circle if p==3, a random sample of the boundary
%               (p-2)-sphere otherwise -- see vMFcap.m). Scalar. Default
%               100.
%                 Example - 'ellipsenpts',200
%  markersize : Marker size used for the points in the 'points'/'both'
%               plots. Scalar. Default 12.
%                 Example - 'markersize',20
%        msg  : Level of output to display. Scalar. Default 1.
%      Ysave  : Store original data matrix. Logical. Default false.
%      plots  : Plot the resulting classification. Scalar. Default 0.
%               If plots=1, produces a 3D scatter with points and/or
%               group confidence caps (see plottype), colored by final
%               assignment (trimmed units in red). For p==3 this is the
%               native space; for p~=3 (e.g. embeddings with hundreds of
%               dimensions) a 3D view must be supplied via plotdims.
%    plottype : What to draw when plots=1. String, one of 'points'
%               (default, as before), 'ellipse' (only the group
%               confidence caps from out.ellipse, outliers omitted since
%               they have no group), or 'both'.
%                 Example - 'plottype','both'
%    plotdims : 3D view of Y used for plotting when p~=3. Either a
%               vector of 3 distinct coordinate indices between 1 and p
%               (a simple choice of 3 original variables), or a p-by-3
%               matrix (e.g. the leading 3 loadings from your own PCA, or
%               any other linear projection you have decided on) used as
%               Y*plotdims. Required if plots=1 and p~=3: which 3
%               coordinates or which projection best represents
%               high-dimensional directional data is left to you --
%               tclustmvmf does not choose one. Ignored if p==3 unless
%               explicitly supplied (default uses the native 3
%               coordinates directly). When a projection is used, group
%               confidence caps are shown as an (approximate) projected
%               point cloud rather than an exact circle -- see "More
%               About" and vMFcap.m.
%                 Example - 'plotdims',[3 7 42]
%      nocheck: Skip the automated NaN/Inf check on Y. Scalar. Default 0.
%       noise : Add an extra, (k+1)-th mixture component modelling
%               uniform "background" noise on the hypersphere, alongside
%               the k vMF clusters. Logical. Default false. The idea
%               (see "More About") is that trimming and this noise
%               component play complementary roles: trimming remains a
%               distribution-free safeguard against an arbitrary
%               alpha-fraction of contamination of any shape, while the
%               noise component explicitly models the specific case of
%               contamination that is itself uniformly spread over the
%               sphere, without spending any of the trimming budget on
%               it. Units assigned to this component get out.idx=k+1,
%               distinct from both a real cluster (1,...,k) and a
%               trimmed unit (0); its estimated mixing proportion is
%               returned in out.noiseprop.
%                 Example - 'noise',true
%
%  Output:
%
%  out :     A structure containing the following fields:
%
%         out.idx  = n-by-1 vector containing the cluster assignment of
%                    each unit; trimmed units have idx=0. If the noise
%                    option is used, units assigned to the uniform
%                    background component instead have idx=k+1 (distinct
%                    from both a real cluster and a trimmed unit).
%       out.muopt  = k-by-p matrix containing the estimated mean
%                    directions (unit-norm rows).
%    out.kappaopt  = k-by-1 vector containing the estimated (restricted)
%                    concentration parameters.
%    out.postprob  = n-by-k matrix of posterior probabilities.
%     out.ellipse  = k-by-1 cell array; out.ellipse{j} is an
%                    ellipsenpts-by-p matrix of points on the boundary of
%                    the conflev confidence cap of group j, in the
%                    native p-dimensional space (see vMFcap.m): an exact
%                    circle if p==3, a Monte Carlo sample of the boundary
%                    otherwise. Project (e.g. via plotdims) before
%                    superimposing on a 3D view of Y.
%      out.conflev = confidence level used for out.ellipse, as the
%                    optional input argument.
%         out.siz  = matrix as returned by tabulate(out.idx): group id,
%                    size, percentage.
%           out.h  = number of observations used to compute the
%                    centroids (untrimmed units).
%         out.obj  = value of the maximized trimmed (classification or
%                    mixture, depending on mixt) log-likelihood.
%       out.NlogL  = -2 times the classification log-likelihood of the
%                    untrimmed units evaluated at the final estimates.
%           out.bs = indices of the elemental subset which gave rise to
%                    the optimal solution.
%     out.fullsol  = nselected-by-1 vector with the objective function
%                    value found in each of the nselected subsamples.
%    out.notconver = fraction of subsamples for which the concentration
%                    steps did not converge within refsteps iterations.
% out.equalweights = as the optional input argument.
%         out.mixt = as the optional input argument.
%        out.noise = as the optional input argument.
%    out.noiseprop = (only if noise=true) estimated mixing proportion of
%                    the uniform background component among the h
%                    retained units, nopt(k+1)/h.
%      out.CLACLA  = classification-likelihood BIC (only if mixt==0):
%                    2*out.NlogL + nParam*log(h).
%   out.MIXMIX/MIXCLA = mixture/classification BIC using mixture-based
%                    parameter estimates (only if mixt>0), mirroring
%                    tclust's New BIC / New ICL.
%            out.Y = original data matrix, stored only if Ysave is true.
%
% More About:
%
% The concentration steps alternate: (i) compute, for every unit and
% every group, the (weighted) vMF log-density log(pi_j)+log
% f(y;mu_j,kappa_j) via logvMFpdfFS; (ii) trim, keeping the h units with
% the largest value (over the criterion selected by mixt); (iii)
% re-estimate mu_j as the direction of the resultant vector of the units
% currently in group j, and kappa_j via kappaFS applied to the mean
% resultant length; (iv) restrict the kappa_j's so that
% max(kappa_j)/min(kappa_j) <= restrfactor.
%
% Note on the restriction step: tclust restricts the eigenvalues of
% Sigma_j via restreigen.m, which minimizes, over a shared floor m, the
% niini-weighted loss sum( (niini/n).*(log(e)+d./e) ) subject to
% e in [m, restrfactor*m] -- d being the raw per-direction sample
% eigenvalue and e the restricted one -- and exploits the fact that this
% specific loss is linear enough in e for the optimal m to have a closed
% algebraic form. That formula is tied to the Gaussian per-direction
% likelihood and is not a generic "restrict any positive scalars"
% routine, so tclustmvmf does not call restreigen directly on the
% kappa_j's. Instead restrKappaFS.m applies the same clip-to-[m,c*m]
% principle with the vMF-consistent loss -logCp(kappa)-kappa*r_j in
% place of log(e)+d/e (this loss is convex in kappa for the same reason
% the Gaussian one is, so the per-group constrained optimum for fixed m
% is still the raw estimate clipped to [m,c*m]); since logCp involves a
% Bessel function, the optimal m is found by direct 1-D numerical
% minimization rather than the closed-form algebra restreigen uses. See
% restrKappaFS.m for details and a derivation of the convexity claim.
%
% Note on plotting for p~=3 (e.g. embeddings, where p can be in the
% hundreds): the clustering itself is entirely dimension-agnostic, but
% visualizing it is not -- there is no canonical way to view directional
% data living on S^(p-1) in 2 or 3 dimensions. tclustmvmf therefore
% requires an explicit 'plotdims' (3 coordinates, or a projection matrix
% you supply, e.g. from your own PCA) rather than picking one for you;
% which view is "relevant" for a given high-dimensional problem is left
% open. Points are shown at their raw (generally shrunk, since a linear
% projection of a unit vector has norm <=1) projected coordinates, and
% group confidence caps -- exact circles only when p==3 -- are shown as
% an approximate projected point cloud for any other view (see
% vMFcap.m).
%
% Note on the noise option: adding a uniform "background" component
% alongside the real clusters is the vMF counterpart of the noise
% component introduced by Banfield and Raftery (1993) for Gaussian
% mixtures, where it is a uniform distribution over the convex hull of
% the data. That original proposal is known to have two issues (Hennig,
% 2004; Coretto and Hennig): it is not breakdown-robust, and because the
% convex hull depends on the data itself, the resulting density is not a
% proper, fixed statistical model, which causes the corresponding
% "maximum likelihood estimator" to lack standard asymptotic properties;
% Coretto and Hennig's own fix is to replace it with a fixed constant
% density that does not depend on the data. On the hypersphere this
% entire difficulty is sidestepped automatically: S^(p-1) is compact, so
% the uniform distribution on it is already a fixed, proper,
% parameter-free density (it is exactly the kappa=0 case of the vMF
% density used throughout, via logvMFpdfFS/vMFlogC) -- there is no
% data-dependent "convex hull" to define in the first place. This does
% not by itself establish breakdown-robustness of the resulting
% estimator (no such analysis is attempted here), but it does mean the
% noise component's density is well defined regardless of the sample.
% Combining an explicit noise component with impartial trimming, rather
% than using either mechanism alone, lets the two play complementary
% roles: trimming keeps its usual role as a distribution-free safeguard
% against an arbitrary alpha-fraction of contamination of unknown shape,
% while the noise component explicitly absorbs the specific case of
% contamination that is itself uniformly spread over the sphere, without
% spending trimming budget on it.
%
% References:
%
%   Banfield, J.D. and Raftery, A.E. (1993), Model-Based Gaussian and
%   Non-Gaussian Clustering, Biometrics, Vol. 49, 803-821.
%
%   Garcia-Escudero, L.A., Gordaliza, A., Matran, C. and Mayo-Iscar, A. (2008),
%   A General Trimming Approach to Robust Cluster Analysis. Annals
%   of Statistics, Vol. 36, 1324-1345.
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions, Journal of Machine Learning Research, Vol. 6.
%
% See also: tclust, tkmeans, logvMFpdfFS, kappaFS, restrKappaFS, vMFlogC, vMFcap
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustmvmf')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % tclustmvmf on a synthetic directional data set with two vMF
    % clusters (p=3) and 10% uniform background noise on the sphere.
    rng(1)
    n1=150; n2=150; nnoise=30;
    mu1=[1 0 0]; kappa1=30;
    mu2=[0 1 0]; kappa2=15;
    % crude vMF (von Mises-Fisher) sampler for illustration (not the
    % production sampler)
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF(mu1,kappa1,n1); simplevMF(mu2,kappa2,n2)];
    noise = normalize_rows(randn(nnoise,3));
    Y = [Y; noise];
    out = tclustmvmf(Y,2,0.1,50,'plots',1);                      % points only (default)
    out = tclustmvmf(Y,2,0.1,50,'plots',1,'plottype','both');    % points + 95% confidence caps
    out = tclustmvmf(Y,2,0.1,50,'plots',1,'plottype','ellipse'); % caps only, no points

    % Custom version of the 'both' plot, useful as a starting point if
    % you want more control than the built-in 'plots' option gives (e.g.
    % to add camlight/lighting, change the viewpoint, or export a
    % figure). Group confidence caps come straight from out.ellipse (see
    % vMFcap.m); no need to recompute them here.
    % A fixed discrete palette -- MATLAB's default axes colors (its
    % dark-red 7th swapped for dark gray, since red is reserved below)
    % plus a few extra well-separated RGB colors -- is cycled across the
    % groups found by the algorithm (there are typically few); red is
    % reserved exclusively for the trimmed/outlier group (idx==0).
    palette = [0.0000 0.4470 0.7410   % blue
               0.8500 0.3250 0.0980   % orange
               0.9290 0.6940 0.1250   % yellow
               0.4940 0.1840 0.5560   % purple
               0.4660 0.6740 0.1880   % green
               0.3010 0.7450 0.9330   % light blue
               0.2500 0.2500 0.2500   % dark gray
               0.0000 0.5000 0.5000   % teal
               0.6000 0.4000 0.2000   % brown
               0.5000 0.0000 0.5000   % violet
               0.5000 0.5000 0.0000   % olive
               0.0000 0.0000 0.5000   % navy
               1.0000 0.6000 0.7840   % pink
               0.0000 0.5000 0.0000]; % dark green
    outlierColor = [1 0 0];

    idx = out.idx;
    groups = setdiff(unique(idx(:))', 0);
    k = numel(groups);

    pointColors = zeros(size(Y,1),3);
    pointColors(idx==0,:) = repmat(outlierColor, sum(idx==0), 1);
    groupColors = zeros(k,3);
    for i = 1:k
        groupColors(i,:) = palette(mod(i-1,size(palette,1))+1, :);
        pointColors(idx==groups(i),:) = repmat(groupColors(i,:), sum(idx==groups(i)), 1);
    end

    [Yx,Yy,Yz] = sphere(60);
    figure;
    surf(Yx, Yy, Yz, 'FaceColor',[0.9 0.9 1], 'EdgeColor','none', 'FaceAlpha',0.25);
    hold on;
    scatter3(Y(:,1), Y(:,2), Y(:,3), 12, pointColors, 'filled');   % small points
    for i = 1:k
        P = out.ellipse{i};
        plot3(P(:,1), P(:,2), P(:,3), '-', 'Color', groupColors(i,:), 'LineWidth', 2);
    end
    axis equal; view(3); xlabel('X'); ylabel('Y'); zlabel('Z');
    %camlight headlight; lighting gouraud;
%}

%{
    % High-dimensional case (p=50, e.g. embeddings). The clustering
    % itself needs nothing special; only plotting requires choosing a 3D
    % view via 'plotdims', which is not picked automatically here.
    rng(2)
    p=50; n1=200; n2=200; nnoise=40;
    mu1=zeros(1,p); mu1(1)=1;   kappa1=25;
    mu2=zeros(1,p); mu2(2)=1;   kappa2=15;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,p)/sqrt(kappa));
    Y = [simplevMF(mu1,kappa1,n1); simplevMF(mu2,kappa2,n2)];
    noise = normalize_rows(randn(nnoise,p));
    Y = [Y; noise];

    % (i) view through 3 arbitrarily chosen original coordinates
    out = tclustmvmf(Y,2,0.1,50,'plots',1,'plottype','both','plotdims',[1 2 3]);

    % (ii) view through a projection onto a 3D subspace (here just the
    % top 3 principal components, purely as one reasonable choice among
    % many -- tclustmvmf leaves this choice to the user)
    [~,~,V] = svd(Y - mean(Y), 'econ');
    W = V(:,1:3);           % p-by-3 projection matrix
    out = tclustmvmf(Y,2,0.1,50,'plots',1,'plottype','both','plotdims',W);
%}

%{
    % The 'noise' option: an explicit uniform-background component
    % alongside the k real clusters, used together with (a typically
    % much smaller) trimming level. Here the contamination is genuinely
    % uniform on the sphere, exactly the case the noise component
    % targets, so most of it should be picked up by out.idx==k+1 (shown
    % in gray) rather than needing to be trimmed (out.idx==0, red).
    rng(3)
    n1=150; n2=150; nnoise=60;   % 20% uniform contamination
    mu1=[1 0 0]; kappa1=30;
    mu2=[0 1 0]; kappa2=15;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF(mu1,kappa1,n1); simplevMF(mu2,kappa2,n2)];
    noiseY = normalize_rows(randn(nnoise,3));
    Y = [Y; noiseY];

    out = tclustmvmf(Y,2,0.02,50,'noise',true,'plots',1,'plottype','both');
    fprintf('Estimated noise proportion: %.1f%% (true contamination: %.1f%%)\n', ...
        100*out.noiseprop, 100*nnoise/(n1+n2+nnoise));
%}

%% Beginning of code

nnargin=nargin;
vvarargin=varargin;
Y = aux.chkinputM(Y,nnargin,vvarargin);
[n, p]=size(Y);

if nargin<3 || isempty(alpha)
    alpha=0.05;
    warning('FSDA:tclustmvmf:Wrongalpha','You have not specified alpha: it is set to 0.05 by default');
end

if nargin<4 || isempty(restrfactor)
    restrfactor=100;
    warning('FSDA:tclustmvmf:Wrongrestrfact','You have not specified restrfactor: it is set to 100 by default');
end

if restrfactor<1
    error('FSDA:tclustmvmf:WrongRestrfactor','restrfactor must be a scalar >=1')
end

if alpha<0
    error('FSDA:tclustmvmf:WrongAlpha','alpha must be a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
end
if alpha>=1
    h=n-floor(alpha);
else
    h=fix(n*(1-alpha));
end

elemsize=k*(p+1);
ncomb=bc(n,elemsize);
nsampdef=min(300,ncomb);
refstepsdef=15;
reftoldef=1e-06;

options=struct('nsamp',nsampdef,'plots',0,'plottype','points','nocheck',0,'msg',1,'Ysave',0,...
    'refsteps',refstepsdef,'reftol',reftoldef,'equalweights',false,...
    'mixt',0,'normalize',false,'newtonsteps',2,'conflev',0.95,'ellipsenpts',100,...
    'markersize',12,'plotdims',[],'noise',false);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:tclustmvmf:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:tclustmvmf:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin>4
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    if options.nsamp>ncomb
        disp('Number of subsets to extract greater than nchoosek(n,k*(p+1)). It is set to that value.');
        options.nsamp=0;
    elseif options.nsamp<0
        error('FSDA:tclustmvmf:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
    end
end

nsamp=options.nsamp;
refsteps=options.refsteps;
reftol=options.reftol;
equalweights=options.equalweights;
mixt=options.mixt;
msg=options.msg;
plots=options.plots;
plottype=options.plottype;
normalizeOpt=options.normalize;
newtonsteps=options.newtonsteps;
conflev=options.conflev;
ellipsenpts=options.ellipsenpts;
markersize=options.markersize;
plotdims=options.plotdims;
noise=options.noise;

% kn is the total number of *modeled* components: the k real vMF
% clusters plus, if noise is requested, one extra uniform-on-the-sphere
% component (see "More About"). Used throughout to size the
% discriminant/posterior-probability matrices; muini/kappaini/rini stay
% k-long since the noise component has no mean direction or
% concentration to estimate.
kn=k+double(noise);
if noise
    lognoiseconst=vMFlogC(0,p);
end

if conflev<=0 || conflev>=1
    error('FSDA:tclustmvmf:WrongConflev','conflev must be a scalar in the interval (0,1)')
end
if ~any(strcmp(plottype,{'points','ellipse','both'}))
    error('FSDA:tclustmvmf:WrongPlottype','plottype must be one of ''points'', ''ellipse'', ''both''')
end

% Validate plotdims upfront (before the expensive clustering loop) if a
% plot has actually been requested. plotdims lets you view data/caps
% living in R^p (p possibly in the hundreds, e.g. embeddings) through
% either 3 chosen coordinates or a linear (e.g. PCA) projection down to
% 3 dimensions; tclustmvmf does not choose this view for you -- for
% p~=3 you must supply one.
if plots==1
    if isempty(plotdims)
        if p~=3
            error('FSDA:tclustmvmf:PlotdimsRequired',['Y has p=%d columns (not 3): plotting requires a 3D view of the ' ...
                'data, specified via the ''plotdims'' option, either as 3 coordinate indices (e.g. [1 2 3]) ' ...
                'or as a %d-by-3 projection matrix'],p,p);
        end
    elseif isvector(plotdims)
        if numel(plotdims)~=3 || any(plotdims<1) || any(plotdims>p) || any(plotdims~=fix(plotdims)) || numel(unique(plotdims))~=3
            error('FSDA:tclustmvmf:WrongPlotdims','As a vector, plotdims must contain 3 distinct integers between 1 and p')
        end
    elseif ismatrix(plotdims)
        if ~isequal(size(plotdims),[p 3])
            error('FSDA:tclustmvmf:WrongPlotdims','As a matrix, plotdims must be p-by-3 (p=%d)',p)
        end
    else
        error('FSDA:tclustmvmf:WrongPlotdims','plotdims must be either a vector of 3 coordinate indices or a p-by-3 projection matrix')
    end
end

% Check / enforce unit-norm rows
rownorms=sqrt(sum(Y.^2,2));
if normalizeOpt
    Y=Y./rownorms;
else
    if max(abs(rownorms-1))>1e-06
        error('FSDA:tclustmvmf:NotUnitNorm','Rows of Y must be unit vectors (use the ''normalize'' option to rescale automatically)')
    end
end

tolrestreigen=1e-08;

if coder.target('MATLAB')
    verLess2016b=verLessThanFS('9.1');
else
    verLess2016b=true;
end

%% Combinatorial part to extract the subsamples
[C,nselected] = subsets(nsamp,n,elemsize,ncomb,msg);

vopt=-1e+30;
fullsol=zeros(nselected,1);
noconv=0;
muopt=[];
kappaopt=[];
nopt=[];
bs=[];

for i=1:nselected

    index=C(i,:);

    niini=zeros(kn,1);
    muini=zeros(k,p);
    kappaini=zeros(k,1);
    rini=zeros(k,1);

    for j=1:k
        ilow=(j-1)*(p+1)+1;
        iup=j*(p+1);
        selj=index(ilow:iup);
        Yselj=Y(selj,:);
        Rj=sum(Yselj,1);
        normRj=norm(Rj);
        niini(j)=p+1;
        if normRj>1e-12
            muini(j,:)=Rj/normRj;
        else
            % near-antipodal cancellation in the seed: fall back to the
            % first unit of the elemental subset as the seed direction
            muini(j,:)=Yselj(1,:);
        end
        rini(j)=normRj/(p+1);
        kappaini(j)=kappaFS(rini(j),p,newtonsteps);
    end

    if noise
        % Seed the noise component's weight on the same footing as a
        % real group's elemental subsample (p+1 units): with no
        % dedicated seeding step of its own (there is nothing to seed --
        % a uniform density has no free parameters), starting it at
        % niini(k+1)=0 would give it log(0/h)=-Inf in the very first
        % discriminant and it could never be chosen (a permanent,
        % self-inflicted "empty component" from which hard/soft EM
        % cannot recover).
        niini(k+1)=p+1;
    end

    kappaini=restrKappaFS(kappaini,rini,niini(1:k),p,restrfactor,tolrestreigen);

    iter=0;
    mudiff=1e+15;
    ind=zeros(n,1);
    postprob=zeros(n,kn);
    obj=0;
    ll=zeros(n,kn);
    groupind=[];
    qq=[];
    qqunassigned=[];

    while (mudiff>reftol) && (iter<refsteps)
        iter=iter+1;

        for j=1:k
            if equalweights
                ll(:,j)=logvMFpdfFS(Y,muini(j,:),kappaini(j));
            else
                ll(:,j)=log(niini(j)/h)+logvMFpdfFS(Y,muini(j,:),kappaini(j));
            end
        end
        if noise
            if equalweights
                ll(:,k+1)=lognoiseconst;
            else
                ll(:,k+1)=log(niini(k+1)/h)+lognoiseconst;
            end
        end

        if mixt==2
            postprobold=postprob;
            [~,postprob,disc]=estepFS(ll,verLess2016b);
            [~,qq]=sort(disc,'descend');
            qqunassigned=qq((h+1):n);
            qq=qq(1:h);
            postprob(qqunassigned,:)=0;
            niini=(sum(postprob))';
        else
            indold=ind;
            [disc,ind]=max(ll,[],2);
            [~,qq]=sort(disc,'descend');
            qqunassigned=qq((h+1):n);
            qq=qq(1:h);
            groupind=ind(qq);
            ind(qqunassigned)=0;
        end

        if mixt==1
            postprobold=postprob;
            [~,postprob]=estepFS(ll,verLess2016b);
            postprob(qqunassigned,:)=0;
            niini=(sum(postprob))';
        end

        % M-step
        for j=1:k
            if mixt>=1
                % niini(j) was already updated above as sum(postprob(:,j))
                if niini(j)>0
                    Rj=sum(Y.*postprob(:,j),1);
                else
                    Rj=zeros(1,p);
                end
            else
                % niini(j) must be refreshed here: it is NOT updated
                % anywhere else in the crisp (mixt==0) path, unlike the
                % mixt>=1 branches which recompute it from postprob just
                % above. Without this line niini stays frozen at its
                % seed value, so rj=norm(Rj)/niini(j) below would divide
                % the resultant vector of the *current* (possibly much
                % larger) group by the *original* seed size -- silently
                % producing rj>>1 and, once clipped inside kappaFS, a
                % spuriously huge kappa. This is the single most
                % important line in the whole M-step.
                niini(j)=sum(groupind==j);
                if niini(j)>0
                    Rj=sum(Y(qq(groupind==j),:),1);
                else
                    Rj=zeros(1,p);
                end
            end

            % A group needs at least one untrimmed unit to update; with
            % niini(j)==0 it is marked empty (NaN), exactly as tclust
            % handles an empty component (cini becomes 0/0=NaN there;
            % here norm(Rj)==0 plays the same role). Unlike tclust,
            % nothing further is required from niini(j) alone: a group
            % with very few units can still push its raw kappaHat to an
            % enormous value (r landing close to 1 by chance, unlike a
            % Gaussian covariance from few points which stays finite),
            % but this is exactly what restrKappaFS is designed to rein
            % in via the shared-likelihood restriction below, mirroring
            % how tclust relies on restreigen alone (no extra minimum
            % group size) to tame a near-zero eigenvalue from a small
            % group.
            if niini(j)>0 && norm(Rj)>1e-12
                muini(j,:)=Rj/norm(Rj);
                rini(j)=norm(Rj)/niini(j);
                kappaini(j)=kappaFS(rini(j),p,newtonsteps);
            else
                muini(j,:)=NaN(1,p);
                kappaini(j)=NaN;
                rini(j)=NaN;
                niini(j)=0;
            end
        end

        if noise && mixt==0
            % For mixt>=1, niini(k+1) was already refreshed above (as
            % part of niini=(sum(postprob))', which spans all kn
            % columns); the crisp path is the only one that needs an
            % explicit update here, since its per-group M-step loop above
            % only ever touches j=1:k. Unlike a real group, niini(k+1)==0
            % is not treated as a failure requiring the whole candidate
            % to be discarded -- it simply means this particular
            % configuration finds no evidence of background noise, which
            % is a perfectly valid outcome, not a degenerate one.
            niini(k+1)=sum(groupind==(k+1));
        end

        kappaini=restrKappaFS(kappaini,rini,niini(1:k),p,restrfactor,tolrestreigen);

        % Objective function
        Ytri=Y(qq,:);
        obj=0;
        if mixt>=1
            log_lh=zeros(h,kn);
            for j=1:k
                if ~isnan(kappaini(j))
                    log_lh(:,j)=log(niini(j)/h)+logvMFpdfFS(Ytri,muini(j,:),kappaini(j));
                else
                    log_lh(:,j)=-Inf;
                end
            end
            if noise
                if equalweights
                    log_lh(:,k+1)=lognoiseconst;
                else
                    log_lh(:,k+1)=log(niini(k+1)/h)+lognoiseconst;
                end
            end
            obj=estepFS(log_lh,verLess2016b);
        else
            for j=1:k
                if niini(j)>0 && ~isnan(kappaini(j))
                    Ytri_j=Y(qq(groupind==j),:);
                    if equalweights
                        obj=obj+sum(logvMFpdfFS(Ytri_j,muini(j,:),kappaini(j)));
                    else
                        obj=obj+niini(j)*log(niini(j)/h)+sum(logvMFpdfFS(Ytri_j,muini(j,:),kappaini(j)));
                    end
                end
            end
            if noise && niini(k+1)>0
                % Every noise-assigned unit contributes the same constant
                % log-density, so the sum over niini(k+1) of them
                % collapses to a plain product.
                if equalweights
                    obj=obj+niini(k+1)*lognoiseconst;
                else
                    obj=obj+niini(k+1)*log(niini(k+1)/h)+niini(k+1)*lognoiseconst;
                end
            end
        end

        if mixt>0
            mudiff=sum(sum(abs(postprob-postprobold)))/n;
        else
            mudiff=sum(abs(indold-ind)>0)/n;
        end

        if iter==refsteps
            noconv=noconv+1;
        end
    end

    fullsol(i)=obj;

    if obj>=vopt && ~any(isnan(muini(:,1)))
        vopt=obj;
        muopt=muini;
        kappaopt=kappaini;
        nopt=niini;
        bs=index;
    end
end

if isempty(muopt)
    error('FSDA:tclustmvmf:NoConvergence',['None of the %d subsamples produced a valid (non-degenerate) ' ...
        'solution; try increasing nsamp or refsteps, or lowering k'],nselected)
end

if nselected>0 && noconv/nselected>0.1
    disp('------------------------------')
    disp(['Warning: Number of subsets without convergence equal to ' num2str(100*noconv/nselected) '%'])
end

%% Final assignment with the optimal parameters
ll=zeros(n,kn);
if equalweights
    for j=1:k
        ll(:,j)=logvMFpdfFS(Y,muopt(j,:),kappaopt(j));
    end
else
    for j=1:k
        ll(:,j)=log(nopt(j)/h)+logvMFpdfFS(Y,muopt(j,:),kappaopt(j));
    end
end
if noise
    if equalweights
        ll(:,k+1)=lognoiseconst;
    else
        ll(:,k+1)=log(nopt(k+1)/h)+lognoiseconst;
    end
end

[~,postprob,logpdf]=estepFS(ll,verLess2016b);
[disc,idx]=max(ll,[],2);
[~,qq]=sort(disc,'descend');

if mixt>=1
    [~,qqmixt]=sort(logpdf,'descend');
    unassignedmixt=qqmixt((h+1):n);
    assignedmixt=qqmixt(1:h);
    [~,idxmixt]=max(postprob,[],2);
    idxmixt(unassignedmixt)=0;
    postprob(unassignedmixt,:)=0;
    out.idx=idxmixt;
    NlogLmixt=-estepFS(ll(assignedmixt,:),verLess2016b);
    out.NlogLmixt=2*NlogLmixt;
else
    unassigned=qq((h+1):n);
    idx(unassigned)=0;
    postprob(unassigned,:)=0;
    out.idx=idx;
end

loglik=disc(qq(1:h));
NlogL=-sum(loglik);

% Number of estimated parameters: k mean directions, each with p-1 free
% parameters (a unit vector on S^(p-1) has p-1 degrees of freedom), plus
% k concentration parameters (the noise component, if present, has
% neither: a uniform density is parameter-free), plus kn-1 mixing
% proportions if these are estimated rather than fixed equal (kn=k+1
% proportions summing to 1 when noise is on, mirrors tclust's
% nParam/BIC logic).
nParam=k*(p-1)+k;
if ~equalweights
    nParam=nParam+(kn-1);
end
logh=log(h);

if mixt>0
    MIXMIX = 2*NlogLmixt + nParam*logh;
    MIXCLA = 2*NlogL + nParam*logh;
else
    CLACLA = 2*NlogL + nParam*logh;
end

out.muopt=muopt;
out.kappaopt=kappaopt;
out.postprob=postprob;

% out.ellipse{j} contains ellipsenpts-by-p points on the boundary of the
% conflev confidence cap (the directional analogue of a confidence
% ellipse) of group j, in the native p-dimensional space -- see
% vMFcap.m: exact for p==3, a Monte Carlo boundary sample otherwise.
% Trimmed units have no group and therefore no entry. Project these
% (together with Y) onto whatever 2D/3D view is of interest for p~=3
% (see 'plotdims' and the p>3 example above for how tclustmvmf's own
% plotting does this).
out.ellipse=cell(k,1);
for j=1:k
    out.ellipse{j}=vMFcap(muopt(j,:),kappaopt(j),conflev,ellipsenpts);
end
out.conflev=conflev;
out.siz=tabulate(out.idx);
out.h=h;
out.obj=vopt;
out.NlogL=2*NlogL;
out.equalweights=equalweights;
out.mixt=mixt;
out.noise=noise;
if noise
    % Diagnostic: estimated share of the retained h units explained by
    % the uniform background component rather than by any real cluster.
    out.noiseprop=nopt(k+1)/h;
end
out.bs=bs;
out.fullsol=fullsol;
out.notconver=noconv/max(nselected,1);
if mixt>0
    out.MIXMIX=MIXMIX;
    out.MIXCLA=MIXCLA;
else
    out.CLACLA=CLACLA;
end

if options.Ysave
    out.Y=Y;
end

%% Plots
if plots==1
    % Fixed discrete palette: MATLAB's default axes colors (its
    % dark-red 7th color swapped for dark gray) plus a few extra
    % well-separated RGB colors, cycled across the groups found by
    % the algorithm (there are typically few); red is reserved
    % exclusively for the trimmed/outlier group (idx==0).
    palette = [0.0000 0.4470 0.7410   % blue
               0.8500 0.3250 0.0980   % orange
               0.9290 0.6940 0.1250   % yellow
               0.4940 0.1840 0.5560   % purple
               0.4660 0.6740 0.1880   % green
               0.3010 0.7450 0.9330   % light blue
               0.2500 0.2500 0.2500   % dark gray
               0.0000 0.5000 0.5000   % teal
               0.6000 0.4000 0.2000   % brown
               0.5000 0.0000 0.5000   % violet
               0.5000 0.5000 0.0000   % olive
               0.0000 0.0000 0.5000   % navy
               1.0000 0.6000 0.7840   % pink
               0.0000 0.5000 0.0000]; % dark green
    outlierColor = [1 0 0];
    noiseColor = [0.6 0.6 0.6];

    idxplot=out.idx;
    if noise
        groupsplot=setdiff(unique(idxplot(:))',[0,k+1]);
    else
        groupsplot=setdiff(unique(idxplot(:))',0);
    end

    pointColors=zeros(n,3);
    pointColors(idxplot==0,:)=repmat(outlierColor,sum(idxplot==0),1);
    if noise
        pointColors(idxplot==(k+1),:)=repmat(noiseColor,sum(idxplot==(k+1)),1);
    end
    groupColors=zeros(numel(groupsplot),3);
    for j=1:numel(groupsplot)
        groupColors(j,:)=palette(mod(j-1,size(palette,1))+1,:);
        pointColors(idxplot==groupsplot(j),:)=repmat(groupColors(j,:),sum(idxplot==groupsplot(j)),1);
    end

    % View: identity for native p==3 with no explicit plotdims, coordinate
    % selection if plotdims is a 3-element index vector, or a linear
    % projection if plotdims is a p-by-3 matrix. tclustmvmf does not pick
    % this view for you when p~=3 (plotdims is then required, checked
    % above): which 3 coordinates or which projection best represents
    % high-dimensional directional data (e.g. embeddings) is left open.
    if isempty(plotdims)
        Yview=Y;
        capview=out.ellipse;
        isnativeview=true;
        axlabels={'X','Y','Z'};
    elseif isvector(plotdims)
        Yview=Y(:,plotdims);
        capview=cellfun(@(P) P(:,plotdims), out.ellipse, 'UniformOutput', false);
        isnativeview=(p==3);
        axlabels={sprintf('dim %d',plotdims(1)),sprintf('dim %d',plotdims(2)),sprintf('dim %d',plotdims(3))};
    else
        Yview=Y*plotdims;
        capview=cellfun(@(P) P*plotdims, out.ellipse, 'UniformOutput', false);
        isnativeview=false;
        axlabels={'proj 1','proj 2','proj 3'};
    end

    figure;
    [sx,sy,sz]=sphere(60);
    surf(sx,sy,sz,'FaceColor',[0.9 0.9 1],'EdgeColor','none','FaceAlpha',0.25);
    hold on
    if any(strcmp(plottype,{'points','both'}))
        scatter3(Yview(:,1),Yview(:,2),Yview(:,3),markersize,pointColors,'filled');
    end
    if any(strcmp(plottype,{'ellipse','both'}))
        for j=1:numel(groupsplot)
            Pcap=capview{j};
            if isnativeview
                % exact, ordered circle: draw as a closed curve
                plot3(Pcap(:,1),Pcap(:,2),Pcap(:,3),'-','Color',groupColors(j,:),'LineWidth',2);
            else
                % projected/high-dimensional case: Pcap is an unordered
                % Monte Carlo sample of the boundary, and/or the
                % projection has left the sphere, so it is shown as a
                % point cloud (its "shadow") rather than a connected
                % curve -- an approximate view of the confidence region,
                % exact only when isnativeview is true.
                scatter3(Pcap(:,1),Pcap(:,2),Pcap(:,3),6,groupColors(j,:),'filled');
            end
        end
    end
    axis equal
    view(3)
    xlabel(axlabels{1}); ylabel(axlabels{2}); zlabel(axlabels{3});
    if noise
        noisestr=sprintf(', noise prop.=%.1f%% (gray)',100*out.noiseprop);
    else
        noisestr='';
    end
    if isnativeview
        title(sprintf('%d groups found by tclustmvmf for alpha=%.2f, restrfactor=%.0f%s',numel(groupsplot),alpha,restrfactor,noisestr));
    else
        title(sprintf(['%d groups found by tclustmvmf for alpha=%.2f, restrfactor=%.0f%s\n' ...
            '(p=%d, showing a 3D view via ''plotdims''; caps are approximate/projected)'],numel(groupsplot),alpha,restrfactor,noisestr,p));
    end
    hold off
end

end
%FScategory:CLUS-RobClaMULT
