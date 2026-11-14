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
%               subsets. With init='kmeans++' it is the number of random
%               starts (default 300), see init.
%                 Example - 'nsamp',500
%                 Data Types - double
%        init : Initialisation of the starting solutions. Character, matrix,
%               vector or cell array. 'auto' (default) uses 'subsets'
%               unless k*(p+1) is larger than n/2 (typical of high
%               dimensional data), in which case it uses 'kmeans++'.
%               'subsets' starts each of the nsamp trials from an
%               elemental subset of k*(p+1) random units, as in tclust.
%               'kmeans++' starts each trial from k units chosen with the
%               k-means++ rule on the cosine distance and from the
%               partition obtained by assigning every unit to the closest
%               of them. With 'kmeans++' nsamp must be positive and its
%               default is 300. A user-guided start can be given as a
%               k-by-p matrix of starting mean directions, or as a vector
%               of n labels (0 for units which are not used, 1,...,k for
%               the groups, k+1 for the noise units if noise=true); a
%               cell array contains several of them, one per trial. In
%               this case nsamp is ignored. See More About.
%                 Example - 'init','kmeans++'
%                 Data Types - char | double | cell
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
%               out.ellipse and, when plots=1, drawn in the plot.
%               Scalar in (0,1). The cap of group j is the region of
%               S^(p-1) where the vMF density exceeds a threshold, i.e.
%               the set of directions y with mu_j'*y>=t_j, whose
%               probability under vMF(mu_j,kappa_j) is conflev.
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
%               Scalar. Default 0. If plots=1 the units are shown
%               coloured by final assignment (trimmed units as red
%               crosses, noise units in gray when noise=true), with
%               and/or without the confidence caps of the groups (see
%               plottype). The plot is produced by sphereplot when p=2
%               (unit circle) or p=3 (unit sphere, which can be rotated
%               with the mouse) and by sphereprojplot when p>3 (the data
%               are projected on a 3-dimensional space, see plotdims).
%               Call these two functions directly, with out.idx,
%               out.muopt and out.kappaopt, to personalise the plot.
%                 Example - 'plots',1
%                 Data Types - double
%    plottype : What to draw when plots=1. String, one of 'points'
%               (default), 'ellipse' (confidence caps only, trimmed units
%               omitted since they carry no group), or 'both'.
%                 Example - 'plottype','both'
%                 Data Types - char|string
%    plotdims : 3-dimensional space used for the plot when p>3. It is
%               passed to sphereprojplot as option proj. If it is empty
%               (default) the space is chosen automatically: it is the
%               space which contains the estimated mean directions of the
%               groups, completed if k<3 with the principal directions of
%               the data orthogonal to them. Otherwise it is either a
%               vector of 3 distinct coordinate indices in 1:p (a simple
%               selection of original variables), or a p-by-3 matrix
%               which defines a linear projection, or one of the
%               character vectors 'means' and 'pca' (see sphereprojplot).
%               It is ignored when p<=3.
%                 Example - 'plotdims',[3 7 42]
%                 Data Types - double | char
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
%                    (p-2)-sphere if p>3. It can be passed to sphereplot
%                    or sphereprojplot as option cap.
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
%           out.bs = indices of the elemental subset (init='subsets') or of
%                    the k seed units (init='kmeans++') which gave rise to
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
% INITIALISATION (OPTION init)

% With init='auto' the choice is made on the number of units needed by an
% elemental subset: if k*(p+1)>n/2 init='kmeans++' is used, because such
% subsets are not a small random sample any more (and for k*(p+1)>n they
% cannot be extracted). A line is displayed when msg>0.
%
% K-MEANS++ (init='kmeans++')
%
% Each of the nsamp trials is started as follows (function vmfseedkmpp,
% at the end of this file):
%  (a) the first seed is a unit chosen uniformly at random; the j-th seed
%      is a unit chosen with probability proportional to d_i^2, where
%      d_i=1-y_i'*s is the cosine distance of unit i from the closest
%      seed s already chosen (seeds themselves have probability 0);
%  (b) every unit is assigned to the seed with the largest cosine, which
%      is the partition of one step of spherical k-means; the seeds
%      belong to their own group;
%  (c) mu_j is the direction of the resultant vector of group j and
%      r_j=||R_j||/n_j, as in any M-step. The initial weights are
%      n_j*h/n; with noise=true the noise component gets 10% of h and the
%      groups share the other 90%.
% USER-GUIDED START (init is a matrix, a vector or a cell array)
%
% If init is a k-by-p matrix each row is a starting mean direction (it is
% normalised): the units are assigned to the closest direction as in step
% (b) above, and (c) follows; a direction which attracts no unit starts
% with r_j=0 and weight 1. If init is a vector with n elements the groups
% are defined by the labels (0=ignored, j=1,...,k group j, k+1 noise
% unit if noise=true) and step (c) follows, with weights proportional to
% the number of units of each label (10% to the noise if no unit has
% label k+1). Each group needs at least one unit. If init is a cell array
% each element, of one of these two types, is a trial: nsamp is ignored
% and the best trial is kept. out.bs is empty.
%
% From this point the trial is identical to a trial started from an
% elemental subset (concentration steps, restriction, objective).
% The cost of (a)-(b) is k matrix-vector products with Y, that is
% O(n*p*k), against the O(k*p*p) work of the elemental subsets, and no
% unit has to be exchanged with another group to obtain a valid start
% even if p is of the order of n. Since the rule favours units which are
% far from the seeds, trials whose seeds fall on outliers or noise are
% possible: this is the role of the nsamp trials, and the best one is
% kept, as in the subsets approach.
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
% HIGH-DIMENSIONAL PLOTTING (p>3)
%
% The clustering itself is entirely dimension-agnostic, but visualising
% S^(p-1) in 3D is not. When p>3 the plot is produced by sphereprojplot,
% which projects the data orthogonally on a 3-dimensional space (by
% default the one which contains the estimated mean directions; see
% option plotdims). Differently from the case p<=3 this plot is not an
% exact representation: the projected units lie inside the unit sphere,
% the groups look tighter than they are because only 2 of the p-1
% directions of within-group variability are visible, and each cap is
% represented by the circle which delimits its projection. The title of
% the plot reports the fraction of the squared length of the units which
% is retained by the projection. See the help of sphereprojplot.
%
% HIGH-DIMENSIONAL COMPUTATIONS
%
% For p in the hundreds or thousands the Bessel functions which appear
% in the normalizing constant and in the equation for kappa cannot be
% represented in double precision. vMFlogC and kappaFS use logBesseliFS,
% which returns their logarithm in all cases.
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
% See also: tclust, tkmeans, tclustkent, sphereplot, sphereprojplot, vMFrnd, logvMFpdfFS, kappaFS, restrKappaFS, vMFlogC, vMFcap
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
    Y = [vMFrnd(n1,mu1,kappa1); vMFrnd(n2,mu2,kappa2)];
    % uniform noise on the sphere (vMF with kappa=0)
    noisePts = vMFrnd(nnoise,[0 0 1],0);
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
    Y = [vMFrnd(n1,[1 0 0],30); vMFrnd(n2,[0 1 0],15)];
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
    Y = [vMFrnd(n1,[1 0 0],20); vMFrnd(n2,[0 1 0],10)];
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
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'nsamp',100,'plots',0);
%}
%
%{
    %% High dimension: starting solutions chosen with the cosine distance.
    % With p=300 and k=3 an elemental subset would need k*(p+1)=903 units,
    % more than the n=490 units of the sample: the default init='subsets'
    % cannot be used, whereas init='kmeans++' needs k units per trial.
    % The data are 3 vMF groups plus uniform noise. The first 3 axes of
    % the space are used to show the result.
    rng(1)
    p=300;
    mu=eye(p);
    Y = [vMFrnd(150,mu(1,:),300); vMFrnd(150,mu(2,:),300); ...
         vMFrnd(150,mu(3,:),300); vMFrnd(40,mu(1,:),0)];
    out = tclustvmf(Y,3,0.05,20,'init','kmeans++','nsamp',50,'noise',true, ...
        'plots',1,'plotdims',[1 2 3]);
    disp(out.siz)
%}
%
%{
    %% User-guided start: starting directions and starting labels.
    % init can be a k-by-p matrix of starting mean directions (here the
    % first two axes of the space) or a vector with the n labels of an
    % initial partition (0 for units which are not used). A cell array
    % would contain several starts.
    rng(1)
    Y = [vMFrnd(200,[1 0 0],25); vMFrnd(200,[0 1 0],12)];
    out1 = tclustvmf(Y,2,0.05,20,'init',[1 0 0; 0 1 0],'plots',0);
    lab = zeros(400,1); lab(1:20)=1; lab(201:220)=2;
    out2 = tclustvmf(Y,2,0.05,20,'init',lab,'plots',0);
    fprintf('Objective: %.2f (directions), %.2f (labels)\n',out1.obj,out2.obj);
%}
%
%{
    %% Effect of refsteps and reftol: concentration step control.
    % Fewer steps speed up the loop; tighter tolerance improves convergence.
    rng(1)
    n1=150; n2=150;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
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
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'newtonsteps',5,'plots',0);
%}
%
%{
    %% The conflev option: confidence level for the group caps.
    % out.ellipse{j} contains the boundary of the conflev cap of group j.
    rng(1)
    n1=150; n2=150;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'conflev',0.90,'plots',1,'plottype','ellipse');
%}
%
%{
    %% The ellipsenpts option: boundary resolution of the confidence caps.
    % More points give a smoother circle when p==3.
    rng(1)
    n1=150; n2=150;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'ellipsenpts',200,'plots',0);
    disp(size(out.ellipse{1}));
%}
%
%{
    %% The markersize option: point size in scatter plots.
    % Increase it for small samples or high-resolution displays.
    rng(1)
    n1=100; n2=100;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'plots',1,'markersize',20);
%}
%
%{
    %% The msg option: suppress or increase command-window output.
    % msg=0 runs silently; msg=1 (default) prints basic progress.
    rng(1)
    n1=150; n2=150;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'msg',0,'plots',0);
%}
%
%{
    %% The Ysave option: store the (possibly normalised) data in out.Y.
    rng(1)
    n1=150; n2=150;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'Ysave',true,'plots',0);
    disp(size(out.Y));
%}
%
%{
    %% The plots option: draw the classification scatter.
    % plots=1 opens a 3D sphere figure with the final assignment.
    rng(1)
    n1=150; n2=150;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'plots',1);
%}
%
%{
    %% The plottype option: points, ellipse or both.
    % 'ellipse' draws only the confidence caps; 'both' adds the data points.
    rng(1)
    n1=150; n2=150;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
    out = tclustvmf(Y,2,0.05,20,'plots',1,'plottype','both');
%}
%
%{
    %% High-dimensional data and the plotdims option.
    % High-dimensional case (p=50, e.g. embeddings). The clustering
    % itself needs nothing special. The plot is produced by
    % sphereprojplot: by default the data are projected on the
    % 3-dimensional space which contains the estimated mean directions
    % of the groups; with option plotdims the user chooses 3 coordinates
    % or supplies a projection matrix.
    rng(2)
    p=50; n1=200; n2=200; nnoise=40;
    mu1=zeros(1,p); mu1(1)=1;   kappa1=100;
    mu2=zeros(1,p); mu2(2)=1;   kappa2=60;
    Y = [vMFrnd(n1,mu1,kappa1); vMFrnd(n2,mu2,kappa2)];
    % uniform noise on the hypersphere (vMF with kappa=0)
    noisePts = vMFrnd(nnoise,mu1,0);
    Y = [Y; noisePts];
    % default projection
    out = tclustvmf(Y,2,0.1,50,'plots',1,'plottype','both');
    % coordinates 1, 2 and 3
    out = tclustvmf(Y,2,0.1,50,'plots',1,'plottype','both','plotdims',[1 2 3]);
%}
%
%{
    %% The nocheck option: skip the NaN/Inf input validation.
    % Useful for speed when the data are guaranteed to be clean.
    rng(1)
    n1=150; n2=150;
    Y = [vMFrnd(n1,[1 0 0],25); vMFrnd(n2,[0 1 0],12)];
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
    Y = [vMFrnd(n1,mu1,kappa1); vMFrnd(n2,mu2,kappa2)];
    % uniform noise on the sphere (vMF with kappa=0)
    noiseY = vMFrnd(nnoise,[0 0 1],0);
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
    'markersize',12,'plotdims',[],'noise',false,'init','auto');

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
end

% Initialisation: userinit is true if the starts are given by the user
% (cell array startlist), kmpp is true if the starts are chosen with
% k-means++, otherwise they are elemental subsets
init=options.init;
userinit=~ischar(init);
kmpp=false;
if userinit
    if ~iscell(init)
        init={init};
    end
    startlist=init(:);
    for i=1:numel(startlist)
        vmfcheckstart(startlist{i},n,p,k,options.noise);
    end
    options.nsamp=numel(startlist);
else
    if ~any(strcmpi(init,{'auto','subsets','kmeans++'}))
        error('FSDA:tclustvmf:WrongInit','init must be ''auto'', ''subsets'', ''kmeans++'', a matrix, a vector or a cell array')
    end
    if strcmpi(init,'auto')
        kmpp = elemsize>n/2;
        if kmpp && options.msg>0
            fprintf('init=''kmeans++'' chosen because k*(p+1)=%d is larger than n/2=%g\n',elemsize,n/2);
        end
    else
        kmpp=strcmpi(init,'kmeans++');
    end
    if kmpp
        % ncomb is not relevant: nsamp is the number of random starts
        if ~any(strcmp(UserOptions,'nsamp'))
            options.nsamp=300;
        elseif ~(options.nsamp>0 && options.nsamp==fix(options.nsamp))
            error('FSDA:tclustvmf:WrongNsamp','With init=''kmeans++'' nsamp must be a positive integer');
        end
    elseif nargin>4
        if options.nsamp>ncomb
            disp('Number of subsets to extract greater than nchoosek(n,k*(p+1)). It is set to that value.');
            options.nsamp=0;
        elseif options.nsamp<0
            error('FSDA:tclustvmf:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
        end
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
% log-density of the uniform distribution on S^(p-1)
lognoiseconst=vMFlogC(0,p);

if conflev<=0 || conflev>=1
    error('FSDA:tclustvmf:WrongConflev','conflev must be a scalar in the interval (0,1)')
end
if ~any(strcmp(plottype,{'points','ellipse','both'}))
    error('FSDA:tclustvmf:WrongPlottype','plottype must be one of ''points'', ''ellipse'', ''both''')
end

% Validate plotdims upfront (before the expensive clustering loop) if a
% plot has actually been requested. plotdims is used only when p>3: it
% is the 3-dimensional space on which sphereprojplot projects the data.
% If it is empty the space is chosen by sphereprojplot.
if plots==1 && p>3 && ~isempty(plotdims)
    if ischar(plotdims)
        if ~any(strcmpi(plotdims,{'auto','means','pca'}))
            error('FSDA:tclustvmf:WrongPlotdims','As a character vector, plotdims must be ''auto'', ''means'' or ''pca''')
        end
    elseif isvector(plotdims)
        if numel(plotdims)~=3 || any(plotdims<1) || any(plotdims>p) || any(plotdims~=fix(plotdims)) || numel(unique(plotdims))~=3
            error('FSDA:tclustvmf:WrongPlotdims','As a vector, plotdims must contain 3 distinct integers between 1 and p')
        end
    elseif ~isequal(size(plotdims),[p 3])
        error('FSDA:tclustvmf:WrongPlotdims','As a matrix, plotdims must be p-by-3 (p=%d)',p)
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
if kmpp || userinit
    % nsamp random starts or starts given by the user: no elemental
    % subsets are needed
    C=[];
    nselected=nsamp;
else
    [C,nselected] = subsets(nsamp,n,elemsize,ncomb,msg);
end

vopt=-1e+30;
fullsol=zeros(nselected,1);
noconv=0;
muopt=[];
kappaopt=[];
nopt=[];
bs=[];

for i=1:nselected

    if userinit
        % starting groups given by the user (see vmfseeduser)
        [muini,rini,niini]=vmfseeduser(Y,startlist{i},k,h,logical(noise));
        index=[];
        kappaini=kappaFS(rini,p,newtonsteps);
    elseif kmpp
        % k-means++ seeds on the cosine distance, partition of the units
        % around them and its resultant vectors (see vmfseedkmpp)
        [muini,rini,niini,index]=vmfseedkmpp(Y,k,h,noise);
        kappaini=kappaFS(rini,p,newtonsteps);
    else
        index=C(i,:);

        % Seeds. The elemental subset is made of k blocks of p+1 units: the
        % resultant vectors of the k blocks are computed in one shot (row j
        % of Rk is the sum of the units of block j).
        Rk=reshape(sum(reshape(Y(index,:),p+1,k,p),1),k,p);
        normR=sqrt(sum(Rk.^2,2));
        muini=Rk./normR;
        degen=normR<=1e-12;
        if any(degen)
            % near-antipodal cancellation in the seed: fall back to the
            % first unit of the block as the seed direction
            muini(degen,:)=Y(index((find(degen)-1)*(p+1)+1),:);
        end
        rini=normR/(p+1);
        kappaini=kappaFS(rini,p,newtonsteps);

        % Each group starts with weight p+1. The noise component (if any) is
        % seeded on the same footing: with no dedicated seeding step of its
        % own (a uniform density has no free parameters), starting it at
        % niini(k+1)=0 would give it log(0/h)=-Inf in the very first
        % discriminant and it could never be chosen (a permanent,
        % self-inflicted "empty component" from which hard/soft EM cannot
        % recover).
        niini=repmat(p+1,kn,1);
    end

    kappaini=restrKappaFS(kappaini,rini,niini(1:k),p,restrfactor,tolrestreigen);

    iter=0;
    mudiff=1e+15;
    ind=zeros(n,1);
    postprob=zeros(n,kn);
    obj=0;
    % needll = true if the matrix of log-discriminants has to be computed
    % at the beginning of the iteration
    needll=true;

    while (mudiff>reftol) && (iter<refsteps)
        iter=iter+1;

        % E-step: n-by-kn matrix of log-discriminants, computed with a
        % single matrix product Y*muini'
        if needll
            ll=logdiscvmf(Y,muini,kappaini,niini,h,equalweights,noise,lognoiseconst);
        end

        if mixt==2
            postprobold=postprob;
            [~,postprob,disc]=estepFS(ll,verLess2016b);
        else
            indold=ind;
            [disc,ind]=max(ll,[],2);
        end

        % Trimming: the h units with the largest value of disc are kept
        % (no sorting is needed if there is no trimming)
        if h<n
            [~,qq]=sort(disc,'descend');
            qqunassigned=qq((h+1):n);
            qq=qq(1:h);
        else
            qq=(1:n)';
            qqunassigned=zeros(0,1);
        end

        if mixt==2
            postprob(qqunassigned,:)=0;
            niini=(sum(postprob,1))';
        else
            ind(qqunassigned)=0;
            if mixt==1
                postprobold=postprob;
                [~,postprob]=estepFS(ll,verLess2016b);
                postprob(qqunassigned,:)=0;
                niini=(sum(postprob,1))';
            end
        end

        % M-step: resultant vectors of all the groups (row j of Rk)
        if mixt>=1
            % niini was already updated above as sum(postprob,1)
            Rk=postprob(:,1:k)'*Y;
        else
            % Z is the n-by-kn sparse indicator matrix of the groups of
            % the untrimmed units. The sizes niini MUST be refreshed
            % here, because in the crisp (mixt==0) path they are not
            % updated anywhere else: with the sizes frozen at their seed
            % value, rini below would be the norm of the resultant of
            % the current group divided by the size of the seed.
            % niini(k+1), the size of the noise component, is refreshed
            % too. Unlike a real group, niini(k+1)==0 is not a failure:
            % it simply means that this configuration finds no evidence
            % of background noise.
            rows = qq(:);
            cols = ind(rows);
            if any(cols < 1 | cols > kn | isnan(cols))
                error('FSDA:tclustvmf:WrongGroupIndex', ...
                    'Cluster labels must be integers between 1 and %d', kn);
            end
            Z = sparse(rows, cols, 1, n, kn);
            niini=full(sum(Z,1))';
            Rk=full(Z(:,1:k)'*Y);
        end
        normR=sqrt(sum(Rk.^2,2));
        nik=niini(1:k);

        % A group needs at least one untrimmed unit to update; with
        % nik(j)==0 it is marked empty (NaN), exactly as tclust handles
        % an empty component. Nothing further is required from nik(j)
        % alone: a group with very few units can still push its raw
        % kappa to an enormous value (r close to 1 by chance), but this
        % is exactly what restrKappaFS is designed to rein in, mirroring
        % how tclust relies on restreigen alone to tame a near-zero
        % eigenvalue from a small group.
        okg=nik>0 & normR>1e-12;
        muini=NaN(k,p);
        rini=NaN(k,1);
        kappaini=NaN(k,1);
        muini(okg,:)=Rk(okg,:)./normR(okg);
        rini(okg)=normR(okg)./nik(okg);
        kappaini(okg)=kappaFS(rini(okg),p,newtonsteps);
        nik(~okg)=0;
        niini(1:k)=nik;

        kappaini=restrKappaFS(kappaini,rini,niini(1:k),p,restrfactor,tolrestreigen);

        % Objective function
        if mixt>=1
            % The log-discriminants with the new parameters are needed
            % for the objective function, and they are the same which
            % are needed by the E-step of the next iteration: they are
            % computed only once.
            ll=logdiscvmf(Y,muini,kappaini,niini,h,equalweights,noise,lognoiseconst);
            needll=false;
            llobj=ll(qq,:);
            if equalweights
                % as in tclust, the mixture likelihood of the k groups is
                % computed with the estimated sizes also when the
                % weights are not used in the E-step
                llobj(:,1:k)=llobj(:,1:k)+log(niini(1:k)'/h);
            end
            obj=estepFS(llobj,verLess2016b);
        else
            % Classification likelihood in closed form. mu_j is the
            % direction of the resultant vector R_j of group j, therefore
            % the sum over the units of the group of kappa_j*mu_j'*y_i is
            % kappa_j*norm(R_j), and the contribution of group j is
            % n_j*log(n_j/h) + n_j*log C_p(kappa_j) + kappa_j*norm(R_j).
            okg=okg & ~isnan(kappaini);
            obj=sum(nik(okg).*vMFlogC(kappaini(okg),p) + kappaini(okg).*normR(okg));
            if ~equalweights
                obj=obj+sum(nik(okg).*log(nik(okg)/h));
            end
            if noise && niini(k+1)>0
                % Every noise-assigned unit contributes the same constant
                % log-density
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
ll=logdiscvmf(Y,muopt,kappaopt,nopt,h,equalweights,noise,lognoiseconst);

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
    % The plot is produced by sphereplot (p<=3) or by sphereprojplot
    % (p>3). The caps are computed from the estimated parameters.
    if strcmp(plottype,'points')
        capplot=[];
    else
        capplot=struct('mu',muopt,'kappa',kappaopt);
    end
    if noise
        noiseid=k+1;
        noisestr=sprintf(', noise prop.=%.1f%% (gray)',100*out.noiseprop);
    else
        noiseid=[];
        noisestr='';
    end
    ngroups=numel(setdiff(unique(out.idx(:))',[0,k+1]));
    maintitle=sprintf('%d groups found by tclustvmf for alpha=%.2f, restrfactor=%.0f%s',ngroups,alpha,restrfactor,noisestr);
    plotopt={'group',out.idx,'noisegroup',noiseid,'cap',capplot,'conflev',conflev, ...
        'markersize',markersize,'showpoints',~strcmp(plottype,'ellipse')};

    if p<=3
        hplot=sphereplot(Y,plotopt{:});
        title(hplot,maintitle);
    else
        if isempty(plotdims)
            plotdims='auto';
        end
        [hplot,outproj]=sphereprojplot(Y,'proj',plotdims,'plottitle',false,plotopt{:});
        title(hplot,{maintitle, sprintf(['p=%d projected on 3 dimensions (%s): %.1f%% of the ' ...
            'squared length of the units retained'],p,outproj.type,100*outproj.retained)});
    end
end

end

function vmfcheckstart(start,n,p,k,noise)
%vmfcheckstart checks one user-guided start (see init in tclustvmf)
if isnumeric(start) && isequal(size(start),[k p])
    if any(~isfinite(start(:))) || any(sum(start.^2,2)==0)
        error('FSDA:tclustvmf:WrongInit','A matrix of starting directions must contain finite rows which are not null')
    end
elseif isnumeric(start) && isvector(start) && numel(start)==n
    kn=k+logical(noise);
    lab=start(:);
    if any(lab~=fix(lab)) || any(lab<0) || any(lab>kn)
        error('FSDA:tclustvmf:WrongInit','The labels of init must be integers between 0 and %d',kn)
    end
    cnt=accumarray(lab(lab>0),1,[kn 1]);
    if any(cnt(1:k)==0)
        error('FSDA:tclustvmf:WrongInit','Each of the groups 1,...,%d must contain at least one unit in init',k)
    end
else
    error('FSDA:tclustvmf:WrongInit','init must be a k-by-p matrix (k=%d, p=%d) or a vector with n=%d labels',k,p,n)
end
end

function [muini,rini,niini]=vmfseeduser(Y,start,k,h,noise)
%vmfseeduser starting values from a user-guided start (see init in
%tclustvmf): a k-by-p matrix of directions or a vector of n labels
n=size(Y,1);
if size(start,1)==k && size(start,2)==size(Y,2)
    % directions: assignment to the closest one (spherical k-means step)
    mu0=start./sqrt(sum(start.^2,2));
    [~,g]=max(Y*mu0',[],2);
    cnt=accumarray(g,1,[k 1]);
else
    % labels
    lab=start(:);
    g=lab;
    cnt=accumarray(lab(lab>0),1,[k+double(noise) 1]);
end
used=g>0 & g<=k;
Z=sparse(find(used),g(used),1,n,k);
Rk=full(Z'*Y);
normR=sqrt(sum(Rk.^2,2));
sz=full(sum(Z,1))';
muini=Rk./normR;
rini=normR./sz;
empty=sz==0;
if any(empty)
    % only possible with directions: the group starts from the direction
    % given by the user, with r=0 and weight 1
    muini(empty,:)=mu0(empty,:);
    rini(empty)=0;
    sz(empty)=1;
end
% initial weights (see More About, INITIALISATION)
if noise
    if numel(cnt)==k+1 && cnt(k+1)>0
        w=[sz; cnt(k+1)];
        w=w/sum(w);
    else
        w=[0.9*sz/sum(sz); 0.1];
    end
else
    w=sz/sum(sz);
end
niini=h*w;
end

function [muini,rini,niini,seedidx]=vmfseedkmpp(Y,k,h,noise)
%vmfseedkmpp starting values by k-means++ on the cosine distance
%
%  Y: n-by-p matrix of unit vectors; k: number of groups; h: number of
%  untrimmed units; noise: true if there is a uniform noise component.
%  muini: k-by-p mean directions; rini: k-by-1 mean resultant lengths;
%  niini: (k+noise)-by-1 initial weights; seedidx: k-by-1 seed units.
%  The steps (a)-(c) are those described in the help of tclustvmf
%  (More About, K-MEANS++ INITIALISATION).
n=size(Y,1);
seedidx=zeros(k,1);
seedidx(1)=randi(n);
% cosine distance of every unit from the closest seed, in [0,2]
dmin=max(1-Y*Y(seedidx(1),:)',0);
for j=2:k
    w=dmin.^2;
    w(seedidx(1:j-1))=0;
    if ~(sum(w)>0)
        % all the units coincide with the seeds: uniform choice among
        % the remaining ones
        w=ones(n,1);
        w(seedidx(1:j-1))=0;
    end
    cw=cumsum(w);
    seedidx(j)=min(find(cw>=rand*cw(end),1),n);
    dmin=min(dmin,max(1-Y*Y(seedidx(j),:)',0));
end
% (b) assignment to the seed with the largest cosine; seeds keep their
% own group
[~,g]=max(Y*Y(seedidx,:)',[],2);
g(seedidx)=(1:k)';
Z=sparse((1:n)',g,1,n,k);
sz=full(sum(Z,1))';
% (c) resultant vectors
Rk=full(Z'*Y);
normR=sqrt(sum(Rk.^2,2));
muini=Rk./normR;
rini=normR./sz;
if noise
    niini=[0.9*h*sz/n; 0.1*h];
else
    niini=h*sz/n;
end
end

function ll=logdiscvmf(Y,mu,kappa,ni,h,equalweights,noise,lognoiseconst)
% n-by-kn matrix of log-discriminants: column j (j=1,...,k) contains
% log(ni(j)/h)+log f(y_i; mu_j,kappa_j), without the first term if
% equalweights is true. Empty groups (kappa=NaN) have -Inf. If noise is
% true the last column refers to the uniform component.
k=numel(kappa);
ll=logvMFpdfFS(Y,mu,kappa);
if ~equalweights
    ll=ll+log(ni(1:k)'/h);
end
ll(:,isnan(kappa(:)') | any(isnan(mu),2)')=-Inf;
if noise
    if equalweights
        ll=[ll, repmat(lognoiseconst,size(Y,1),1)];
    else
        ll=[ll, repmat(log(ni(k+1)/h)+lognoiseconst,size(Y,1),1)];
    end
end
end
%FScategory:CLUS-RobClaMULT
