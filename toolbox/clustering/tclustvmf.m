function [out] = tclustvmf(Y,k,alpha,restrfactor,varargin)
%tclustvmf computes trimmed clustering for directional data using mixtures of von Mises-Fisher distributions
%
%<a href="matlab: docsearchFS('tclustvmf')">Link to the help page for this function</a>
%
%   tclustvmf partitions the points in the n-by-p data matrix Y, where
%   each row is assumed to be a unit vector on the hypersphere S^(p-1) in
%   R^p, into k clusters. This is the directional-data counterpart of
%   tclust: the multivariate normal model N(mu_j,Sigma_j) is replaced by
%   the von Mises-Fisher (vMF) model vMF(mu_j,kappa_j), where mu_j (a
%   unit vector) is the mean direction and kappa_j>=0 is the concentration
%   parameter of group j. kappa_j is the vMF analogue of the inverse of
%   Sigma_j: a large kappa_j implies a tightly concentrated, "low
%   variance" cluster around mu_j, exactly as a large eigenvalue of
%   Sigma_j^{-1} does in the Gaussian case. The vMF density on S^(p-1) is
%
%     f(y; mu, kappa) = C_p(kappa) * exp(kappa * mu' * y),
%
%   where C_p(kappa) is the normalizing constant (see vMFlogC.m) and y is
%   a unit vector. Robustness is achieved via the same impartial ("hard")
%   trimming device used by tclust and tkmeans: at most n-h observations
%   are left unassigned, where h = fix(n*(1-alpha)). A constraint on the
%   ratio max(kappa_j)/min(kappa_j) prevents the classification likelihood
%   from becoming unbounded as a cluster collapses onto a single direction
%   (the vMF analogue of tclust's eigenvalue restriction, which prevents
%   Sigma_j from collapsing to a singular matrix).
%
%   Setting restrfactor=1 forces a common kappa across all groups, which
%   reduces the discriminant rule to cosine similarity: this is "trimmed
%   spherical k-means", the exact directional analogue of the relationship
%   between tclust (restrfactor=1) and tkmeans.
%
%   Optionally (see option 'noise'), an extra (k+1)-th component modelling
%   uniform "background" noise on the hypersphere can be added alongside
%   the k vMF clusters. This is used together with a typically smaller
%   trimming level: trimming remains a distribution-free safeguard against
%   contamination of unknown/arbitrary shape, while the noise component
%   explicitly captures contamination that is itself uniform on the sphere.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-p matrix of directional data:
%               each row of Y must be a unit vector (norm 1). Use option
%               'normalize' to have tclustvmf rescale the rows for you;
%               otherwise an error is thrown if rows are not unit-norm.
%            k: Number of groups. Scalar. Positive integer specifying
%               how many vMF clusters to fit.
%        alpha: Global trimming level. Scalar. If 0<=alpha<1 clustering
%               is based on h=fix(n*(1-alpha)) observations; if alpha is
%               an integer >=1 then h=n-floor(alpha). alpha=0 means no
%               trimming.
%  restrfactor: Restriction factor on the concentration parameters.
%               Positive scalar, at least 1. Constrains the ratio
%               max(kappa_j)/min(kappa_j) across the k groups.
%               restrfactor=1 forces a common kappa (trimmed spherical
%               k-means); larger values allow groups of increasingly
%               different concentration. See restrKappaFS.m.
%
%  Optional input arguments:
%
%    equalweights : Cluster weights in the concentration/assignment
%               steps and in the log-likelihood. Logical. If false
%               (default) mixing proportions pi_j=n_j/h are estimated
%               and enter the discriminant, exactly as in tclust. If
%               true all pi_j are fixed to 1/k (or 1/kn with noise),
%               reducing to cosine-distance assignment within each group.
%                 Example - 'equalweights',true
%                 Data Types - logical
%       mixt  : Mixture modelling or crisp assignment. Scalar. The
%               semantics mirror tclust exactly. mixt=0 (default) is the
%               crisp classification likelihood, i.e. the CEM algorithm.
%               mixt=1 uses crisp trimming with a posterior-probability-
%               weighted M-step (a hybrid CEM/EM). mixt=2 is the fully
%               soft trimming criterion based on the mixture log-density,
%               i.e. the EM algorithm. The choice governs which
%               information criterion is returned: CLACLA for mixt=0,
%               MIXMIX and MIXCLA for mixt>0 (see Cerioli et al., 2017).
%                 Example - 'mixt',2
%                 Data Types - double
%       nsamp : Number of subsamples to extract. Scalar. Default is
%               min(300, nchoosek(n,k*(p+1))). Each subsample draws
%               k*(p+1) observations (p+1 per cluster), enough to
%               initialise a non-degenerate vMF estimate (direction from
%               the resultant vector, kappa from the mean resultant
%               length via kappaFS). Set nsamp=0 to enumerate all
%               subsets.
%                 Example - 'nsamp',500
%                 Data Types - double
%     refsteps: Number of concentration steps (CEM/EM iterations) run
%               on each subsample before comparing objective values.
%               Scalar. Default 15.
%                 Example - 'refsteps',10
%                 Data Types - double
%      reftol : Convergence tolerance for the concentration steps. Scalar.
%               The loop stops when the relative change in assignments
%               (mixt=0) or in posterior probabilities (mixt>0) falls
%               below this threshold. Default 1e-06.
%                 Example - 'reftol',1e-08
%                 Data Types - double
%    normalize: Rescale each row of Y to unit norm before clustering.
%               Logical. Default false (an error is thrown instead if
%               rows are not already unit-norm). Set to true when passing
%               raw embedding vectors or other non-normalised directional
%               data.
%                 Example - 'normalize',true
%                 Data Types - logical
%  newtonsteps: Number of Newton-Raphson steps used by kappaFS to invert
%               the mean resultant length r into the concentration
%               parameter kappa (solving A_p(kappa)=r, where A_p is the
%               ratio of Bessel functions I_{p/2}(kappa)/I_{p/2-1}(kappa)).
%               Scalar. Default 2 (see kappaFS for accuracy figures; 2
%               steps give about 1e-6 relative error for most r in (0,1)).
%                 Example - 'newtonsteps',5
%                 Data Types - double
%     conflev : Confidence level for the group caps returned in
%               out.ellipse and, when plots=1, drawn as boundary circles.
%               Scalar in (0,1). The cap is the region on S^(p-1) where
%               the vMF density exceeds its (1-conflev) quantile, i.e.
%               all directions within an angular radius theta_j of mu_j,
%               where theta_j = arccos(1 + log(1-conflev)/kappa_j).
%               Default 0.95. See vMFcap.m.
%                 Example - 'conflev',0.90
%                 Data Types - double
% ellipsenpts : Number of points used to represent each confidence cap
%               boundary in out.ellipse. Scalar. For p==3 these are
%               evenly spaced around the exact boundary circle; for p>3
%               they are a random sample of the boundary (p-2)-sphere
%               (see vMFcap.m). Default 100.
%                 Example - 'ellipsenpts',200
%                 Data Types - double
%  markersize : Marker size used for the scatter points in the 'points'
%               and 'both' plots. Scalar. Default 12.
%                 Example - 'markersize',20
%                 Data Types - double
%        msg  : Level of output displayed during the search. Scalar.
%               Default 1 (prints basic progress). Set to 0 to run
%               silently, or to 2 for more detail.
%                 Example - 'msg',0
%                 Data Types - double
%      Ysave  : Store the (possibly row-normalised) data matrix in the
%               output. Logical. Default false. If true, the matrix
%               used internally (after any normalisation) is returned as
%               out.Y.
%                 Example - 'Ysave',true
%                 Data Types - logical
%      plots  : Plot the resulting classification on the unit sphere.
%               Scalar. Default 0. If plots=1, produces a 3D scatter
%               with points and/or group confidence caps (see plottype),
%               colored by final assignment (trimmed units in red, noise
%               units in gray when noise=true). For p==3 this is the
%               native space; for p~=3 a 3D view must be supplied via
%               plotdims.
%                 Example - 'plots',1
%                 Data Types - double
%    plottype : What to draw when plots=1. String, one of 'points'
%               (default), 'ellipse' (confidence caps only, trimmed units
%               omitted since they carry no group), or 'both'.
%                 Example - 'plottype','both'
%                 Data Types - char|string
%    plotdims : 3D view of Y used for plotting when p~=3. Either a
%               vector of 3 distinct coordinate indices in 1:p (a simple
%               selection of original variables), or a p-by-3 matrix
%               (any linear projection, e.g. leading PCA loadings) used
%               as Y*plotdims. Required if plots=1 and p~=3: which 3
%               coordinates or which projection best represents the
%               high-dimensional sphere is left to the user. When a
%               projection is used, confidence caps are shown as an
%               approximate projected point cloud (see vMFcap.m and
%               "More About" below).
%                 Example - 'plotdims',[3 7 42]
%                 Data Types - double
%      nocheck: Skip the automated NaN/Inf check on Y. Scalar. Default 0.
%               Set to 1 for speed when the data are guaranteed clean.
%                 Example - 'nocheck',1
%                 Data Types - double
%       noise : Add an extra, (k+1)-th mixture component modelling
%               uniform noise on the hypersphere alongside the k vMF
%               clusters. Logical. Default false. The uniform density on
%               S^(p-1) is parameter-free and equals 1/surface_area(S^(p-1))
%               = C_p(0) everywhere (the kappa=0 special case of the vMF
%               density). Its mixing proportion pi_{k+1} is estimated
%               jointly with the k real cluster proportions. Units
%               assigned to this component get out.idx=k+1, distinct from
%               a real cluster (1,...,k) and a trimmed unit (0). The
%               estimated mixing proportion is in out.noiseprop. See
%               "More About" for the motivation and the relationship with
%               Banfield and Raftery (1993).
%                 Example - 'noise',true
%                 Data Types - logical
%
%  Output:
%
%  out :     A structure containing the following fields:
%
%         out.idx  = n-by-1 vector of cluster assignments. Trimmed units
%                    have idx=0. If noise=true, units assigned to the
%                    uniform background component have idx=k+1 (distinct
%                    from both a real cluster and a trimmed unit).
%       out.muopt  = k-by-p matrix of estimated mean directions; each
%                    row is a unit vector.
%    out.kappaopt  = k-by-1 vector of estimated (restricted) concentration
%                    parameters.
%    out.postprob  = n-by-kn matrix of posterior probabilities, where
%                    kn=k when noise=false and kn=k+1 when noise=true.
%                    Column j contains P(z_i=j|y_i) for j=1,...,kn.
%                    The last column (when noise=true) gives the posterior
%                    probability of belonging to the uniform background.
%                    Trimmed units have postprob set to zero.
%     out.ellipse  = k-by-1 cell array; out.ellipse{j} is an
%                    ellipsenpts-by-p matrix of points on the boundary of
%                    the conflev confidence cap of group j in the native
%                    p-dimensional space (see vMFcap.m): an exact circle
%                    if p==3, a Monte Carlo sample of the boundary
%                    (p-2)-sphere otherwise. Project via plotdims before
%                    superimposing on a 3D view of Y.
%      out.conflev = confidence level used for out.ellipse.
%         out.siz  = matrix as returned by tabulate(out.idx): group id,
%                    count, percentage. Rows for the noise group (k+1)
%                    and trimmed units (0) appear when present.
%           out.h  = number of untrimmed observations used to compute the
%                    final parameter estimates.
%         out.obj  = value of the maximized trimmed log-likelihood at the
%                    best solution: classification log-likelihood for
%                    mixt=0, mixture log-likelihood for mixt=2.
%       out.NlogL  = -2 times the classification log-likelihood of the h
%                    untrimmed units at the final estimates, computed as
%                    -2*sum_i max_j [log(pi_j)+log f(y_i;mu_j,kappa_j)]
%                    over i in the top-h retained units. This equals
%                    2*NlogL in the CLACLA and MIXCLA formulae. Available
%                    for all values of mixt.
%  out.NlogLmixt  = (only if mixt>0) -2 times the mixture log-likelihood
%                    of the h assigned units at the final estimates,
%                    computed as -2*sum_i log(sum_j pi_j*f(y_i;mu_j,kappa_j))
%                    over the h units retained by the mixture trimming
%                    criterion. This is the numerator of out.MIXMIX.
%           out.bs = indices of the elemental subset which gave rise to
%                    the best solution (the one maximising out.obj across
%                    all nselected subsamples).
%     out.fullsol  = nselected-by-1 vector of the objective function
%                    value found at the end of each subsample's
%                    concentration steps.
%    out.notconver = fraction of subsamples for which the concentration
%                    steps did not converge within refsteps iterations
%                    (mudiff > reftol at the last step).
% out.equalweights = logical, as the input option.
%         out.mixt = scalar, as the input option.
%        out.noise = logical, as the input option.
%    out.noiseprop = (only if noise=true) estimated mixing proportion
%                    pi_{k+1} of the uniform background component among
%                    the h retained units, taken from the M-step of the
%                    best subsample solution. This is a parametric
%                    estimate (the EM mixing proportion pi_{k+1} = n_{k+1}/h
%                    at convergence of the best subsample) and may differ
%                    from the empirical count sum(out.idx==k+1)/h from
%                    the final hard assignment.
%      out.CLACLA  = (only if mixt==0) classification-likelihood BIC:
%                    out.NlogL + nParam*log(h),
%                    where nParam = k*p + (k-1) [or k*p + k with noise]
%                    and the BIC penalty uses log(h) rather than log(n),
%                    following the FSDA convention for trimmed estimators
%                    (see Cerioli et al., 2017, and "More About" below).
%      out.MIXMIX  = (only if mixt>0) mixture-likelihood BIC:
%                    out.NlogLmixt + nParam*log(h). Uses the mixture
%                    log-likelihood (numerator) and mixture-based
%                    parameter estimates. Analogous to tclust's "New BIC"
%                    (Cerioli et al., 2017).
%      out.MIXCLA  = (only if mixt>0) classification-likelihood BIC at
%                    mixture estimates: out.NlogL + nParam*log(h). Uses
%                    the classification log-likelihood with the
%                    mixture-based parameter estimates. Analogous to
%                    tclust's "New ICL" (Cerioli et al., 2017).
%            out.Y = original (or row-normalised) data matrix, stored
%                    only if Ysave is true.
%
% More About:
%
% CONCENTRATION STEPS (CEM/EM algorithm on the hypersphere)
%
% The algorithm alternates four steps for each random subsample:
%  (i)  E-step: compute log-discriminants
%         d_{ij} = log(pi_j) + log C_p(kappa_j) + kappa_j*(mu_j'*y_i)
%       for each unit i and group j (and the uniform noise term if
%       noise=true); pi_j are the current mixing proportions.
%  (ii) Trimming: retain the h units with the largest d_{ij} value
%       (or largest log-mixture-density for mixt=2).
%  (iii)M-step: re-estimate mu_j as the direction of the resultant
%       vector R_j = sum_{i in group j} y_i, and kappa_j from the mean
%       resultant length r_j = ||R_j||/n_j via kappaFS.
%  (iv) Restriction: apply restrKappaFS to enforce
%       max(kappa_j)/min(kappa_j) <= restrfactor.
%
% PARAMETER COUNT AND INFORMATION CRITERIA
%
% A unit vector on S^(p-1) has p-1 free parameters (one constraint:
% ||mu||=1 from p Cartesian coordinates). The concentration parameter
% kappa_j adds 1 more. The uniform noise component is parameter-free
% (kappa=0 is fixed; its density is determined by the geometry of the
% sphere alone). With kn total components and estimated mixing proportions,
% there are kn-1 free mixing-proportion parameters. The total parameter
% count is therefore:
%
%   nParam = k*(p-1+1)  +  (kn-1)   [equalweights=false]
%          = k*p        +  (kn-1)
%
% where kn=k (noise=false) or kn=k+1 (noise=true).
%
% The three information criteria returned (Cerioli, Garcia-Escudero,
% Mayo-Iscar and Riani, 2017, JCGS) are:
%
%   CLACLA = -2*CLA_CLA + nParam*log(h)   [mixt=0]
%   MIXMIX = -2*MIX_MIX + nParam*log(h)   [mixt>0]
%   MIXCLA = -2*MIX_CLA + nParam*log(h)   [mixt>0]
%
% where CLA_CLA is the classification log-likelihood at CEM estimates,
% MIX_MIX is the mixture log-likelihood at EM estimates, and MIX_CLA is
% the classification log-likelihood evaluated at EM estimates. CLACLA
% and MIXMIX are BIC-type criteria; MIXCLA is an ICL-type criterion.
% NOTE: following the FSDA convention for trimmed estimators, the BIC
% penalty uses log(h) instead of log(n). This differs from the
% original CGMR 2017 paper (which uses log(n)) but is consistent with
% all FSDA functions built on tclust.
%
% RESTRICTION ON CONCENTRATION PARAMETERS
%
% tclust restricts the eigenvalues of Sigma_j via restreigen.m, which
% minimises a Gaussian log-likelihood-consistent weighted loss to find,
% analytically, the optimal common floor m. That formula is tied to the
% Gaussian likelihood and cannot be applied to kappa_j directly.
% restrKappaFS.m instead uses the vMF-consistent loss
%   L(kappa) = -logC_p(kappa) - kappa*r_j
% (the per-group vMF log-likelihood as a function of kappa for fixed r_j),
% which is convex in kappa for the same reason the Gaussian loss is
% convex in the eigenvalue. For fixed m, the constrained optimum of each
% kappa_j in [m, restrfactor*m] is again the unconstrained estimate
% clipped to that interval; the optimal m is found by 1-D minimisation
% (Bessel function involved, no closed form). See restrKappaFS.m.
%
% NOISE COMPONENT AND RELATIONSHIP WITH BANFIELD-RAFTERY (1993)
%
% For Gaussian mixtures, Banfield and Raftery (1993) proposed a uniform
% noise component over the convex hull of the data. That proposal has two
% known problems (Hennig 2004; Coretto and Hennig): (a) the convex-hull
% density depends on the data itself, so the model is not a fixed
% statistical model and the MLE lacks standard properties; (b) it is not
% breakdown-robust. On S^(p-1) both issues vanish automatically: the
% sphere is compact, the uniform distribution on it is a fixed proper
% density (C_p(0) = 1/omega_{p-1}, where omega_{p-1} is the surface area
% of S^(p-1)) that does not depend on the data at all, and is exactly the
% kappa=0 special case of the vMF density already in use. No data-dependent
% convex hull needs to be defined. Combined with impartial trimming, the
% noise component and the trimming mechanism play distinct roles: trimming
% absorbs an arbitrary alpha-fraction of contamination of any shape; the
% noise component explicitly models contamination that is itself uniformly
% spread over the sphere, without spending trimming budget on it.
%
% HIGH-DIMENSIONAL PLOTTING (p~=3)
%
% The clustering itself is entirely dimension-agnostic, but visualising
% S^(p-1) in 3D is not. tclustvmf requires an explicit 'plotdims' (3
% coordinate indices or a p-by-3 projection matrix) rather than picking
% a view automatically; which 3 coordinates or which projection is most
% informative is left to the user. When a linear projection is used,
% points appear at their projected (generally sub-unit-norm) coordinates
% and confidence caps are shown as approximate projected point clouds
% rather than exact circles (see vMFcap.m).
%
% References:
%
%   Banfield, J.D. and Raftery, A.E. (1993), Model-Based Gaussian and
%   Non-Gaussian Clustering, Biometrics, Vol. 49, 803-821.
%
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions, Journal of Machine Learning Research, Vol. 6,
%   pp. 1345-1382.
%
%   Cerioli, A., Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani, M.
%   (2017), Finding the Number of Groups in Model-Based Clustering via
%   Constrained Likelihoods, Journal of Computational and Graphical
%   Statistics, Vol. 27, pp. 404-416. doi:10.1080/10618600.2017.1390469
%   (source of CLACLA/MIXMIX/MIXCLA and their relationship to BIC/ICL)
%
%   Garcia-Escudero, L.A., Gordaliza, A., Matran, C. and Mayo-Iscar, A.
%   (2008), A General Trimming Approach to Robust Cluster Analysis,
%   Annals of Statistics, Vol. 36, pp. 1324-1345.
%   (tclust: the Euclidean parent of tclustvmf)
%
%   Hennig, C. (2004), Breakdown Points for Maximum Likelihood Estimators
%   of Location-Scale Mixtures, Annals of Statistics, Vol. 32, 1313-1340.
%   (limitations of the Banfield-Raftery noise component)
%
%   Mardia, K.V. and Jupp, P.E. (2000), Directional Statistics, Wiley,
%   Chichester. (vMF distribution, mean resultant length, kappa estimation)
%
% See also: tclust, tkmeans, logvMFpdfFS, kappaFS, restrKappaFS, vMFlogC, vMFcap
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustvmf')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% tclustvmf on a synthetic 3-D directional data set, default options.
    % Two vMF clusters (p=3) with 10% trimming and restriction factor 50.
    % Points and 95% confidence caps are drawn.
    rng(1)
    n1=150; n2=150; nnoise=30;
    mu1=[1 0 0]; kappa1=30;
    mu2=[0 1 0]; kappa2=15;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF(mu1,kappa1,n1); simplevMF(mu2,kappa2,n2)];
    noisePts = normalize_rows(randn(nnoise,3));
    Y = [Y; noisePts];
    out = tclustvmf(Y,2,0.1,50,'plots',1,'plottype','both');
%}
%
%{
    %% Effect of the equalweights option.
    % With equalweights=true the mixing proportions are fixed to 1/k
    % rather than estimated from the data.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],30,n1); simplevMF([0 1 0],15,n2)];
    out1 = tclustvmf(Y,2,0.05,20,'equalweights',false,'plots',0);
    out2 = tclustvmf(Y,2,0.05,20,'equalweights',true, 'plots',0);
    disp(out1.kappaopt); disp(out2.kappaopt);
%}
%
%{
    %% Effect of the mixt option: soft vs crisp assignment.
    % mixt=0 (default) is crisp CEM; mixt=2 uses the mixture log-density
    % for trimming and a soft M-step.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],20,n1); simplevMF([0 1 0],10,n2)];
    out0 = tclustvmf(Y,2,0.05,20,'mixt',0,'plots',0);
    out2 = tclustvmf(Y,2,0.05,20,'mixt',2,'plots',0);
    fprintf('mixt=0 obj=%.2f, mixt=2 obj=%.2f\n',out0.obj,out2.obj);
%}
%
%{
    %% Effect of nsamp: number of random subsamples.
    % Increasing nsamp improves solution quality at the cost of runtime.
    rng(1)
    n1=200; n2=200;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'nsamp',100,'plots',0);
%}
%
%{
    %% Effect of refsteps and reftol: concentration step control.
    % Fewer steps speed up the loop; tighter tolerance improves convergence.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'refsteps',10,'reftol',1e-08,'plots',0);
%}
%
%{
    %% The normalize option: rows need not be unit-norm on input.
    % When normalize=true tclustvmf rescales each row to unit norm before
    % clustering, so raw (non-normalised) data can be passed directly.
    rng(1)
    n1=150; n2=150;
    % Generate raw (un-normalised) directional data.
    Y = [randn(n1,3)+repmat([3 0 0],n1,1); randn(n2,3)+repmat([0 3 0],n2,1)];
    out = tclustvmf(Y,2,0.05,20,'normalize',true,'plots',0);
%}
%
%{
    %% The newtonsteps option: Newton-Raphson steps for kappa estimation.
    % More steps give a more accurate inversion of the mean resultant
    % length into a concentration parameter (see kappaFS).
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'newtonsteps',5,'plots',0);
%}
%
%{
    %% The conflev option: confidence level for the group caps.
    % out.ellipse{j} contains the boundary of the conflev cap of group j.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'conflev',0.90,'plots',1,'plottype','ellipse');
%}
%
%{
    %% The ellipsenpts option: boundary resolution of the confidence caps.
    % More points give a smoother circle when p==3.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'ellipsenpts',200,'plots',0);
    disp(size(out.ellipse{1}));
%}
%
%{
    %% The markersize option: point size in scatter plots.
    % Increase it for small samples or high-resolution displays.
    rng(1)
    n1=100; n2=100;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'plots',1,'markersize',20);
%}
%
%{
    %% The msg option: suppress or increase command-window output.
    % msg=0 runs silently; msg=1 (default) prints basic progress.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'msg',0,'plots',0);
%}
%
%{
    %% The Ysave option: store the (possibly normalised) data in out.Y.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'Ysave',true,'plots',0);
    disp(size(out.Y));
%}
%
%{
    %% The plots option: draw the classification scatter.
    % plots=1 opens a 3D sphere figure with the final assignment.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'plots',1);
%}
%
%{
    %% The plottype option: points, ellipse or both.
    % 'ellipse' draws only the confidence caps; 'both' adds the data points.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'plots',1,'plottype','both');
%}
%
%{
    %% The plotdims option: viewing high-dimensional data in 3D.
    % High-dimensional case (p=50, e.g. embeddings). The clustering
    % itself needs nothing special; only plotting requires choosing a 3D
    % view via 'plotdims', which is not picked automatically.
    rng(2)
    p=50; n1=200; n2=200; nnoise=40;
    mu1=zeros(1,p); mu1(1)=1;   kappa1=25;
    mu2=zeros(1,p); mu2(2)=1;   kappa2=15;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,p)/sqrt(kappa));
    Y = [simplevMF(mu1,kappa1,n1); simplevMF(mu2,kappa2,n2)];
    noisePts = normalize_rows(randn(nnoise,p));
    Y = [Y; noisePts];
    out = tclustvmf(Y,2,0.1,50,'plots',1,'plottype','both','plotdims',[1 2 3]);
%}
%
%{
    %% The nocheck option: skip the NaN/Inf input validation.
    % Useful for speed when the data are guaranteed to be clean.
    rng(1)
    n1=150; n2=150;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF([1 0 0],25,n1); simplevMF([0 1 0],12,n2)];
    out = tclustvmf(Y,2,0.05,20,'nocheck',1,'plots',0);
%}
%
%{
    %% The noise option: explicit uniform-background component.
    % An extra (k+1)-th component models uniform noise on the sphere,
    % used here with a small trimming level (2%). Most uniform
    % contamination is absorbed by the noise component (out.idx==k+1,
    % shown in gray) rather than by trimming (out.idx==0, red).
    rng(3)
    n1=150; n2=150; nnoise=60;
    mu1=[1 0 0]; kappa1=30;
    mu2=[0 1 0]; kappa2=15;
    normalize_rows = @(X) X ./ max(vecnorm(X,2,2), eps);
    simplevMF = @(mu,kappa,n) normalize_rows(mu + randn(n,3)/sqrt(kappa));
    Y = [simplevMF(mu1,kappa1,n1); simplevMF(mu2,kappa2,n2)];
    noiseY = normalize_rows(randn(nnoise,3));
    Y = [Y; noiseY];
    out = tclustvmf(Y,2,0.02,50,'noise',true,'plots',1,'plottype','both');
    fprintf('Estimated noise proportion: %.1f%%\n',100*out.noiseprop);
%}

%% Beginning of code

nnargin=nargin;
vvarargin=varargin;
Y = aux.chkinputM(Y,nnargin,vvarargin);
[n, p]=size(Y);

if nargin<3 || isempty(alpha)
    alpha=0.05;
    warning('FSDA:tclustvmf:Wrongalpha','You have not specified alpha: it is set to 0.05 by default');
end

if nargin<4 || isempty(restrfactor)
    restrfactor=100;
    warning('FSDA:tclustvmf:Wrongrestrfact','You have not specified restrfactor: it is set to 100 by default');
end

if restrfactor<1
    error('FSDA:tclustvmf:WrongRestrfactor','restrfactor must be a scalar >=1')
end

if alpha<0
    error('FSDA:tclustvmf:WrongAlpha','alpha must be a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
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
        error('FSDA:tclustvmf:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:tclustvmf:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
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
        error('FSDA:tclustvmf:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
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
    error('FSDA:tclustvmf:WrongConflev','conflev must be a scalar in the interval (0,1)')
end
if ~any(strcmp(plottype,{'points','ellipse','both'}))
    error('FSDA:tclustvmf:WrongPlottype','plottype must be one of ''points'', ''ellipse'', ''both''')
end

% Validate plotdims upfront (before the expensive clustering loop) if a
% plot has actually been requested. plotdims lets you view data/caps
% living in R^p (p possibly in the hundreds, e.g. embeddings) through
% either 3 chosen coordinates or a linear (e.g. PCA) projection down to
% 3 dimensions; tclustvmf does not choose this view for you -- for
% p~=3 you must supply one.
if plots==1
    if isempty(plotdims)
        if p~=3
            error('FSDA:tclustvmf:PlotdimsRequired',['Y has p=%d columns (not 3): plotting requires a 3D view of the ' ...
                'data, specified via the ''plotdims'' option, either as 3 coordinate indices (e.g. [1 2 3]) ' ...
                'or as a %d-by-3 projection matrix'],p,p);
        end
    elseif isvector(plotdims)
        if numel(plotdims)~=3 || any(plotdims<1) || any(plotdims>p) || any(plotdims~=fix(plotdims)) || numel(unique(plotdims))~=3
            error('FSDA:tclustvmf:WrongPlotdims','As a vector, plotdims must contain 3 distinct integers between 1 and p')
        end
    elseif ismatrix(plotdims)
        if ~isequal(size(plotdims),[p 3])
            error('FSDA:tclustvmf:WrongPlotdims','As a matrix, plotdims must be p-by-3 (p=%d)',p)
        end
    else
        error('FSDA:tclustvmf:WrongPlotdims','plotdims must be either a vector of 3 coordinate indices or a p-by-3 projection matrix')
    end
end

% Check / enforce unit-norm rows
rownorms=sqrt(sum(Y.^2,2));
if normalizeOpt
    Y=Y./rownorms;
else
    if max(abs(rownorms-1))>1e-06
        error('FSDA:tclustvmf:NotUnitNorm','Rows of Y must be unit vectors (use the ''normalize'' option to rescale automatically)')
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
    error('FSDA:tclustvmf:NoConvergence',['None of the %d subsamples produced a valid (non-degenerate) ' ...
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
% (see 'plotdims' and the p>3 example above for how tclustvmf's own
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
    % projection if plotdims is a p-by-3 matrix. tclustvmf does not pick
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
        title(sprintf('%d groups found by tclustvmf for alpha=%.2f, restrfactor=%.0f%s',numel(groupsplot),alpha,restrfactor,noisestr));
    else
        title(sprintf(['%d groups found by tclustvmf for alpha=%.2f, restrfactor=%.0f%s\n' ...
            '(p=%d, showing a 3D view via ''plotdims''; caps are approximate/projected)'],numel(groupsplot),alpha,restrfactor,noisestr,p));
    end
    hold off
end

end
%FScategory:CLUS-RobClaMULT
