function [out, varargout] = tclustreg(y,X,k,restrfact,alphaLik,alphaX,varargin)
%tclustreg performs robust linear grouping analysis
%
%<a href="matlab: docsearchFS('tclustreg')">Link to the help function</a>
%
%  Required input arguments:
%
%         y : Response variable. Vector.
%             A vector with n elements that contains the response variable.
%             y can be either a row or a column vector.
%             Data Types - single|double
%
%         X : Explanatory variables (also called 'regressors'). Matrix.
%             Data matrix of dimension $(n \times p-1)$. Rows of X represent
%             observations, and columns represent variables. Missing values
%             (NaN's) and infinite values (Inf's) are allowed, since
%             observations (rows) with missing or infinite values will
%             automatically be excluded from the computations.
%             Data Types - single|double
%
%         k : Number of clusters. Scalar.
%             This is a guess on the number of data groups.
%             Data Types - single|double
%
% restrfact : restriction factor for regression residuals and covariance
%             matrices of the explanatory variables. Scalar or vector with two
%             elements. If restrfact is a scalar it controls the
%             differences among group scatters of the residuals. The value
%             1 is the strongest restriction. If restrfactor is a vector
%             with two elements the first element controls the differences
%             among group scatters of the residuals and the second the
%             differences among covariance matrices of the explanatory
%             variables. Note that restrfactor(2) is used just if
%             input option $alphaX=1$, that is if constrained weighted
%             model for X is assumed.
%            Data Types - single|double
%
%   alphaLik : Trimming level. Scalar.
%            alphaLik is a value between 0 and 0.5 or an  integer specifying
%            the number of observations which have to be trimmed. If
%            alphaLik=0 there is no trimming. More in detail, if 0<alphaLik<1
%            clustering is based on h=fix(n*(1-alphaLik)) observations.
%            Else if alphaLik is an integer greater than 1 clustering is
%            based on h=n-floor(alphaLik). More in detail, likelihood
%            contributions are sorted and the units associated with the
%            smallest n-h contributions are trimmed.
%            Data Types - single|double
%
%   alphaX : Second-level trimming or constrained weighted model for X. Scalar.
%            alphaX is a value in the interval [0 1].
%            - If alphaX=0 there is no second-level trimming.
%            - If alphaX is in the interval [0 0.5] it indicates the
%               fixed proportion of units subject to second level trimming.
%               In this case alphaX is usually smaller than alphaLik.
%               For further details see Garcia-Escudero et. al. (2010).
%            -  If alphaX is in the interval (0.5 1), it indicates a
%               Bonferronized confidence level to be used to identify the
%               units subject to second level trimming. In this case the
%               proportion of units subject to second level trimming is not
%               fixed a priori, but is determined adaptively.
%               For further details see Torti et al. (2018).
%            -  If alphaX=1, constrained weighted model for X is assumed
%               (Gershenfeld, 1997). The CWM estimator is able to
%               take into account different distributions for the explanatory
%               variables across groups, so overcoming an intrinsic limitation
%               of mixtures of regression, because they are implicitly
%               assumed equally distributed. Note that if alphaX=1 it is
%               also possible to apply using restrfactor(2) the constraints
%               on the cov matrices of the explanatory variables.
%               For further details about CWM see Garcia-Escudero et al.
%               (2017) or Torti et al. (2018).
%            Data Types - single|double
%
%
%  Optional input arguments:
%
%     intercept : Indicator for constant term. Scalar. If 1, a model with
%                constant term will be fitted (default), else no constant
%                term will be included.
%                Example - 'intercept',1
%                Data Types - double
%
%       mixt  : mixture modelling or crisp assignment. Scalar.
%               Option mixt specifies whether mixture modelling or crisp
%               assignment approach to model based clustering must be used.
%               In the case of mixture modelling parameter mixt also
%               controls which is the criterior to find the untrimmed units
%               in each step of the maximization
%               If mixt >=1 mixture modelling is assumed else crisp
%               assignment.
%                In mixture modelling the likelihood is given by
%                \[
%                \prod_{i=1}^n  \sum_{j=1}^k \pi_j \phi (y_i \; x_i' , \beta_j , \sigma_j),
%                \]
%               while in crisp assignment the likelihood is given by
%               \[
%               \prod_{j=1}^k   \prod _{i\in R_j} \pi_j  \phi (y_i \; x_i' , \beta_j , \sigma_j),
%               \]
%               where $R_j$ contains the indexes of the observations which
%               are assigned to group $j$,
%               Remark - if mixt>=1 previous parameter equalweights is
%               automatically set to 1.
%               Parameter mixt also controls the criterion to select the units to trim
%               if mixt == 2 the h units are those which give the largest
%               contribution to the likelihood that is the h largest
%               values of
%               \[
%                   \sum_{j=1}^k \pi_j \phi (y_i \; x_i' , \beta_j , \sigma_j)   \qquad
%                    i=1, 2, ..., n
%               \]
%               elseif mixt==1 the criterion to select the h units is
%               exactly the same as the one which is used in crisp
%               assignment. That is: the n units are allocated to a
%               cluster according to criterion
%               \[
%                \max_{j=1, \ldots, k} \hat \pi'_j \phi (y_i \; x_i' , \beta_j , \sigma_j)
%               \]
%               and then these n numbers are ordered and the units
%               associated with the largest h numbers are untrimmed.
%               Example - 'mixt',1
%               Data Types - single | double
%
%equalweights : cluster weights in the concentration and assignment steps.
%               Logical. A logical value specifying whether cluster weights
%               shall be considered in the concentration, assignment steps
%               and computation of the likelihood.
%               if equalweights = true we are (ideally) assuming equally
%               sized groups by maximizing the likelihood. Default value
%               false.
%                 Example - 'equalweights',true
%                 Data Types - Logical
%
%    nsamp : number of subsamples to extract.
%            Scalar or matrix with k*p columns.
%            If nsamp is a scalar it contains the number of subsamples
%            which will be extracted.
%            If nsamp=0 all subsets will be extracted.
%            If the number of all possible subset is <300 the
%            default is to extract all subsets, otherwise just 300.
%            If nsamp is a matrix it contains in the rows the indexes of
%            the subsets which have to be extracted. nsamp in this case can
%            be conveniently generated  by function subsets.
%            nsamp must have k*p columns. The first p columns are used to
%            estimate the regression coefficient of group 1... the last p
%            columns are used to estimate the regression coefficient of
%            group k
%             Example - 'nsamp',1000
%             Data Types - double
%
% refsteps:  Number of refining iterations. Scalar. Number of refining
%               iterations in each subsample.  Default is 10.
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'refsteps',15
%                 Data Types - single | double
%
%     reftol  : Tolerance for the refining steps. Scalar.
%               The default value is 1e-14;
%                 Example - 'reftol',1e-05
%                 Data Types - single | double
%
% commonslope  : Impose constraint of common slope regression coefficients.
%               Boolean.
%               If commonslope is true, the groups are forced to have the
%               same regression coefficients (apart from the intercepts).
%               The default value of commonslope is false;
%                 Example - 'commonslope',true
%                 Data Types - boolean
%
%    plots : Plot on the screen. Scalar. A flag to control the
%            generation of the plots.
%            If plots=1 a plot is showed on the screen with the
%            final allocation (and if size(X,2)==2 with the lines
%            associated to the groups)
%            Example - 'plots',1
%            Data Types - double
%
%   wtrim: Application of observation weights. Scalar or structure. If
%           wtrim is a scalar, a flag taking values
%          in [0, 1, 2, 3, 4], to control the application of weights on the
%          observations for betaestimation.
%          -  If \texttt{wtrim}=0 (no weights) and $\texttt{mixt}=0$, the
%             algorithm reduces to the standard tclustreg algorithm.
%          -  If \texttt{wtrim}=0 and \texttt{mixt}=2, the maximum posterior
%             probability $D_i$ of equation 7 of Garcia et al. 2010 is
%             computing by maximizing the log-likelihood contributions of
%             the mixture model of each observation.
%          -  If \texttt{wtrim} = 1, trimming is done by weighting the
%             observations using values specified in vector \texttt{we}.
%             In this case, vector \texttt{we} must be supplied by the
%             user. For instance, \texttt{we} = $X$.
%          -  If \texttt{wtrim} = 2, trimming is again done by weighting
%             the observations using values specified in vector \texttt{we}.
%             In this case, vector \texttt{we} is computed from the data as
%             a function of the density estimate $\mbox{pdfe}$.
%            Specifically, the weight of each observation is the
%            probability of retaining the observation, computed as
%            \[\mbox{pretain}_{i g} = 1 - \mbox{pdfe}_{ig}/\max_{ig}(\mbox{pdfe}_{ig})\]
%         -  If \texttt{wtrim} = 3, trimming is again done by weighting the
%            observations using values specified in vector \texttt{we}. In
%            this case, each element $we_i$ of vector \texttt{we} is a
%            Bernoulli random variable with probability of success
%            $\mbox{pdfe}_{ig}$. In the clustering framework this is done
%            under the constraint that no group is empty.
%         -  If \texttt{wtrim} = 4, trimming is done with the tandem approach
%            of Cerioli and Perrotta (2014).
%         -  If \texttt{wtrim} = 5 (TO BE IMPLEMENTED)
%          -  If \texttt{wtrim} = 6 (TO BE IMPLEMENTED)
%          If wtrim is a structure, it is composed by:
%         -  wtrim.wtype_beta: the weight for the beta estimation. It can be
%           0, 1, 2, 3, as in the case of wtrim scalar
%         -  wtrim.wtype_obj: the weight for the objective function. It can
%         be:
%             - '0': no weights in the objective function
%             - 'Z': Bernoulli random variable with probability of success
%            $\mbox{pdfe}_{ig}$
%             - 'w': a function of the density estimate $\mbox{pdfe}$.
%             - 'Zw': the product of the two above.
%             - 'user': user weights we.
%            Example - 'wtrim',1
%            Data Types - double
%
%      we: Vector of observation weights. Vector. A vector of size n-by-1
%          containing application-specific weights that the user needs to
%          apply to each observation. Default
%          value is  a vector of ones.
%            Example - 'we',[0.2 0.2 0.2 0.2 0.2]
%            Data Types - double
%
%        cup  :  pdf upper limit. Scalar. The upper limit for the pdf used
%                to compute the retantion probability. If cup = 1
%                (default), no upper limit is set.
%                Data Types - scalar
%                Example - cup, 0.8
%
%      pstar  :  thinning probability. Scalar. Probability with each a unit
%                enters in the thinning procedure. If pstar = 1 (default), all units
%                enter in the thinning procedure.
%                Data Types - scalar
%                Example - pstar, 0.95
%
%       k_dens_mixt: in the Poisson/Exponential mixture density function,
%                    number of clusters for density mixtures. Scalar.
%                    This is a guess on the number of data groups. Default
%                    value is 5.
%            Example - 'k_dens_mixt',6
%            Data Types - single|double
%
%   nsamp_dens_mixt: in the Poisson/Exponential mixture density function,
%                    number of subsamples to extract. Scalar. Default 300.
%                    Example - 'nsamp_dens_mixt',1000
%                    Data Types - double
%
%refsteps_dens_mixt: in the Poisson/Exponential mixture density function,
%                    number of refining iterations. Scalar. Number of refining
%                    iterations in each subsample.  Default is 10.
%                    Example - 'refsteps_dens_mixt',15
%                    Data Types - single | double
%
%  method_dens_mixt: in the Poisson/Exponential mixture density function,
%                    distribution to use. Character. If method_dens_mixt =
%                    'P', the Poisson distribution is used, with
%                    method_dens_mixt = 'E', the Exponential distribution
%                    is used. Default is 'P'.
%                    Example - 'method_dens_mixt','E'
%                    Data Types - char
%
%    msg  : Level of output to display. Scalar.
%           Scalar which controls whether to display or not messages
%           on the screen.
%           If msg==0 nothing is displayed on the screen.
%           If msg==1 (default) messages are displayed
%           on the screen about estimated time to compute the estimator
%           or the number of subsets in which there was no convergence.
%           If msg==2 detailed messages are displayed. For example the
%           information at iteration level.
%             Example - 'msg',1
%             Data Types - single | double
%
%RandNumbForNini: Pre-extracted random numbers to initialize proportions.
%                Matrix. Matrix with size k-by-size(nsamp,1) containing the
%                random numbers which are used to initialize the
%                proportions of the groups. This option is effective just
%                if nsamp is a matrix which contains pre-extracted
%                subsamples. The purpose of this option is to enable the
%                user to replicate the results in case routine tclust is
%                called using a parfor instruction (as it happens for
%                example in routine IC, where tclust is called through a
%                parfor for different values of the restriction factor).
%                The default value of RandNumbForNini is empty that is
%                random numbers from uniform are used.
%                   Example - 'RandNumbForNini',''
%                   Data Types - single | double
%
%      nocheck: Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               vector y and matrix X.
%               As default nocheck=0.
%                   Example - 'nocheck',1
%                   Data Types - single | double
%
%  Output:
%
%  out :  structure containing the following fields
%
%   out.bopt             = $(p+1) \times k$ matrix containing the regression
%                          parameters.
%
%   out.sigma2opt       = $k$ row vector containing the estimated group
%                          variances.
%
%   out.sigma2opt_corr    = $k$ row vector containing the estimated group
%                          variances corrected with  asymptotic consistency
%                          factor.
%
%         out.muXopt= k-by-p matrix containing cluster centroid
%                       locations. Robust estimate of final centroids of
%                       the groups. This output is present only if input
%                       option alphaX is 1.
%
%         out.sigmaXopt= p-by-p-by-k array containing estimated constrained
%                       covariance covariance matrices of the explanatory
%                       variables for the k groups. This output is present
%                       only if input option alphaX is 1.
%
%         out.cstepopt= scalar containing the concentration step where the
%                       objective function was the largest. This is useful
%                       when the objective function is not monotone (e.g.
%                       with second level trimming or with thinning).
%
%         out.subsetopt= scalar containing the subset id where the
%                       objective function was the largest.
%
%
%            out.idx  = n-by-1 vector containing assignment of each unit to
%                       each of the k groups. Cluster names are integer
%                       numbers from -2 to k.
%                       -1 indicates first level trimmed units.
%                       -2 indicates second level trimmed units.
%
%            out.siz  = Matrix of size k-by-3.
%                       1st col = sequence from -2 to k;
%                       2nd col = number of observations in each cluster;
%                       3rd col = percentage of observations in each
%                       cluster;
%                       Remark: 0 denotes thinned units (if the weights
%                       to find thinned units are 0 or 1, -1 indicates
%                       first level trimmed units and -2 indicates second
%                       level trimmed units).
%
%   out.postprobopt   = $n \times k$ matrix containing the final posterior
%                           probabilities. out.postprobopt(i,j) contains
%                           posterior probabilitiy of unit i from component
%                           (cluster) j. For the trimmed units posterior
%                           probabilities are 0. This output is always
%                           produced (independently of the value of mixt).
%
%          out.MIXMIX = BIC which uses parameters estimated using the
%                       mixture loglikelihood and the maximized mixture
%                       likelihood as goodness of fit measure.
%                       Remark: this output is present only if input option
%                       mixt is >0.
%
%          out.MIXCLA = BIC which uses the classification likelihood based on
%                       parameters estimated using the mixture likelihood
%                       (In some books this quantity is called ICL).
%                       This output is present only if input option
%                       mixt is >0.
%
%          out.CLACLA = BIC which uses the classification likelihood based on
%                       parameters estimated using the classification likelihood.
%                       Remark: this output is present only if input option
%                       mixt is =0.
%
%           out.obj   = scalar containing value of the objective function.
%
%          out.NlogL = Scalar. -2 log classification likelihood. In
%                       presence of full convergence -out.NlogL/2 is equal
%                       to out.obj.
%
%      out.NlogLmixt = Scalar. -2 log mixture likelihood. In
%                      presence of full convergence -out.NlogLmixt/2 is
%                      equal to out.obj. If input parameter mixt=0 then
%                      out.NlogLmixt is a missing value.
%
%               out.h = Scalar. Number of observations that have determined the
%                       regression coefficients (number of untrimmed units).
%
%          out.class = 'tclustreg'.
%
%  Optional Output:
%
%            C     : Indexes of extracted subsamples. Matrix.
%                    Matrix of size nsamp-by-k*p containing (in the rows)
%                    the indices of the subsamples extracted for computing
%                    the estimate.
%
% More About:
%
%  Computational issues to be addressed in future releases.
%  [1]
%  FSDA function "wthin" uses the MATLAB function ksdensity. The calls to
%  ksdensity have been optimized. The only possibility to further reduce
%  time execution is to replace ksdensity with a more efficient kernel
%  density estimator.
%  [2]
%  The weighted version of tclustreg requires weighted sampling. This is
%  now implemented in randsampleFS. A computaionally more efficient
%  algorithm, based on a binary tree approach introduced by
%      Wong, C.K. and M.C. Easton (1980) "An Efficient Method for Weighted
%      Sampling Without Replacement", SIAM Journal of Computing,
%      9(1):111-113.
%  is provided by recent releases of the MATLAB function datasample.
%  Unfortunately this function spends most of the self running time useless
%  parameter checking. To copy the function in the FSDA folder
%  FSDA/combinatorial, possibly rename it, and remove the option parameters
%  checks, is not sufficient, as datasample relies on a mex file wswor
%  which is platform dependent. The issue is usually referred to as "code
%  signature".
%  [3]
%  In the plots, the use of text to highlight the groups with their index
%  is terribly slow (more than 10 seconds to generate a scatter of 7000
%  units. ClickableMultiLegend and legend are also slow.
%  [4]
%  FSDA function restreigen could be improved. In some cases it is one of
%  the most expensive functions.
%
% REMARK: trimming vs thinning
% - trimming is the key feature of TCLUST, giving robustness to the model.
% - thinning is a new denoising feature introduced to mitigate the
%   distorting effect of very dense data areas. It is implemented via
%   observation weighting.
% For the sake of code readability, the relevant sections of the code are
% identified with a "TRIMMING" or "THINNING" tag.
%
% REMARK: the number of parameters to penalize the likelihood are given
% below:
% $k(p+1)$ = number of regression coefficients including the intercept.
% $k-1$ = number of proportions -1 (because their sum is 1).
% $(k-1)(1-1/restrfact(1))+1$ = constraints on the $\sigma^2_j$.
% To the above parameters, if $alphaX=1$ we must add:
% $k(p-1)p/2$ = rotation parameters for $\Sigma_X=cov(X)$.
% $(kp-1)(1-1/restrfactor(2))+1$ = eigenvalue parameters for
% $\Sigma_X=cov(X)$.
%
% See also: tclust, tkmeans, rlga
%
%
% References:
%
% Garcia-Escudero, L.A., Gordaliza A., Greselin F., Ingrassia S., and
% Mayo-Iscar A. (2016), The joint role of trimming and constraints in
% robust estimation for mixtures of gaussian factor analyzers,
% "Computational Statistics & Data Analysis", Vol. 99, pp. 131-147.
%
% Garcia-Escudero, L.A., Gordaliza, A., Greselin, F., Ingrassia, S. and
% Mayo-Iscar, A. (2017), Robust estimation of mixtures of regressions with
% random covariates, via trimming and constraints, "Statistics and
% Computing", Vol. 27, pp. 377-402.
%
% Garcia-Escudero, L.A., Gordaliza A., Mayo-Iscar A., and San Martin R.
% (2010), Robust clusterwise linear regression through trimming,
% "Computational Statistics and Data Analysis", Vol. 54, pp.3057-3069.
%
% Cerioli, A. and Perrotta, D. (2014). Robust Clustering Around Regression
% Lines with High Density Regions. Advances in Data Analysis and
% Classification, Vol. 8, pp. 5-26.
%
% Torti F., Perrotta D., Riani, M. and Cerioli A. (2018). Assessing Robust
% Methodologies for Clustering Linear Regression Data, "Advances in Data
% Analysis and Classification".
%
% Copyright 2008-2021.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustreg')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2019-11-20 18:33:09 #$: Date of the last commit
%
%
%
% Examples:
%
%{
    %% tclustreg of 'X data'.
    % The X data have been introduced by Gordaliza, Garcia-Escudero & Mayo-Iscar (2013).
    % The dataset presents two parallel components without contamination.
    X  = load('X.txt');
    y1 = X(:,end);
    X1 = X(:,1:end-1);

    k = 2 ;

    restrfact = 5; alpha1 = 0.05 ; alpha2 = 0.01;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2);

    restrfact = 2; alpha1 = 0.05 ; alpha2 = 0.01;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'mixt',2);

    % Comparison with robust linear grouping
    out = rlga(X,k,alpha2);

    cascade;
%}

%{
    % tclustreg of fishery data 1.
    clear all; close all;
    load fishery;
    X = fishery{:,:};
    % some jittering might be useful if there are many duplicated units
    X = X + 10^(-8) * abs(randn(677,2));

    %tclustreg on fishery data
    y1 = X(:,end);
    X1 = X(:,1:end-1);
    k = 3 ; restrfact = 50; alpha1 = 0.04 ; alpha2 = 0.01;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',0);
    title('TCLUST-REG');

    %lga and rlga on fishery data
    out=lga(X,3);

    alpha = 0.95;
    out=rlga(X,3,1-alpha);

    cascade;
%}

%{
    % tclustreg of fishery data 2.
    clear all; close all;
    load fishery;
    X = fishery{:,:};
    % some jittering is necessary because duplicated units are not treated
    % in tclustreg: this needs to be addressed
    X = X + 10^(-8) * abs(randn(677,2));
    y1 = X(:,end);
    X1 = X(:,1:end-1);
    
    % some arbitrary weights for the units
    we = sqrt(X1)/sum(sqrt(X1));
    
    % tclustreg required parameters
    k = 2; restrfact = 50; alpha1 = 0.04 ; alpha2 = 0.01;

    % now tclust is run on each combination of mixt and wtrim options

    disp('mixt = 0; wtrim = 0;');
    disp('standard tclustreg, with classification likelihood and without thinning' );
    mixt = 0; wtrim = 0;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 0;');
    disp('mixture likelihood, no thinning' );
    mixt = 2; wtrim = 0;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 0; wtrim = 1;');
    disp('classification likelihood, thinning based on user weights' );
    mixt = 0; wtrim = 1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'we',we,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 1;');
    disp('mixture likelihood, thinning based on user weights' );
    mixt = 2; wtrim = 1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'we',we,'wtrim',wtrim);

    disp('mixt = 0; wtrim = 2;');
    disp('classification likelihood, thinning based on retention probabilities' );
    mixt = 0; wtrim = 2; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 2;');
    disp('mixture likelihood, thinning based on retention probabilities' );
    mixt = 2; wtrim = 2; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 0; wtrim = 3;');
    disp('classification likelihood, thinning based on bernoulli weights' );
    mixt = 0; wtrim = 3; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 3;');
    disp('mixture likelihood, thinning based on bernoulli weights' );
    mixt = 2; wtrim = 3; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 0; wtrim = 4;');
    disp('classification likelihood, tandem thinning based on bernoulli weights' );
    mixt = 0; wtrim = 4; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 4;');
    disp('mixture likelihood, tandem thinning based on bernoulli weights' );
    mixt = 2; wtrim = 4; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 0; wtrim = struct; wtrim.wtype_beta=3;  wtrim.wtype_obj=''Z'';')
    disp('classification likelihood, componentwise thinning based on bernoulli weights, objective function based on bernoulli weights' );
    mixt = 0; wtrim = struct; wtrim.wtype_beta=3;  wtrim.wtype_obj='Z';
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 0; wtrim = struct; wtrim.wtype_beta=2;  wtrim.wtype_obj=''w'';')
    disp('classification likelihood, componentwise thinning based on retention probabilities, objective function based on retention probabilities' );
    mixt = 0; wtrim = struct; wtrim.wtype_beta=2;  wtrim.wtype_obj='w';
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    cascade
%}

%{
    % tclustreg of simulated data 1.
    % Generate mixture of regression using MixSimReg, with an average
    % overlapping at centroids = 0.01. Use all default options.
    rng(372,'twister');
    p=3;
    k=2;
    Q=MixSimreg(k,p,'BarOmega',0.001);
    n=500;
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);

    % plot the dataset
    yXplot(y,X);

    % run tclustreg
    out=tclustreg(y,X(:,2:end),k,50,0.01,0.01,'intercept',1);
%}

%{
    % tclustreg of simulated data 2.
    % Generate mixture of regression using MixSimReg, with an average
    % overlapping at centroids =0.01.
    rng(372,'twister');
    p=3;
    k=2;
    Q=MixSimreg(k,p,'BarOmega',0.001);
    n=200;
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);

    % plot the dataset
    yXplot(y,X);

    % Generate the elemental subsets used in tclustreg once and for all.
    nsamp  = 100;
    ncomb  = bc(n,p);
    method = [10*ones(n/2,1); ones(n/2,1)]; % weighted sampling using weights in "method"
    msg    = 0;
    for i=1:k
        C(:,(i-1)*p+1:i*p) = subsets(nsamp, n, p, ncomb, msg, method);
    end

    % tclustreg using samples in C
    out=tclustreg(y,X(:,2:end),k,50,0.01,0.01,'nsamp',C);
%}

%{
    % Example of nsamp passed as a matrix.
    X  = load('X.txt');
    y1 = X(:,end);
    X1 = X(:,1:end-1);
    n = 200;
    k = 2 ;

    restrfact = 5; alpha1 = 0.05 ; alpha2 = 0.01;
    nsamp=200;
    Cnsamp=subsets(nsamp,n,(size(X1,2)+1)*k);
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'nsamp',Cnsamp);
%}

%{
    %% CWM, adaptive and fixed 2nd level trimming.
    % Comparison among CWM, adaptive second level trimming and fixed second
    % level trimming
    close all
    % Generate two groups and a set of concentrated outliers
    rng('default')
    rng(100)
    p=1;
    n1=90;
    x1=(randn(n1,p)+2)*2.1-1;
    y1=-1.2*x1*1+randn(n1,1)*0.6+2.1;

    n2=90;
    x2=(randn(n1,p)+4)*1.5-4;
    y2=5+x2*0.7+randn(n2,p)*0.6;

    n3=20;
    x3=randn(n3,p)*0.1+12-4;
    y3=randn(n3,p)*0.1+16-2;
    %x3=0;
    %y3=0;

    X=[x1;x2;x3];
    y=[y1;y2;y3];

    n=n1+n2+n3;
    group=[ones(n1,1); 2*ones(n2,1); 3*ones(n3,1)];
    yXplot(y,X,group);
    legend('Location','best')
    % gscatter(X,y,group)
    % Run the 3 models and compare the results
    k=2;
    % Specify restriction factor for variance of residuals
    restrfact=5;
    % 10 per cent first level trimming.
    alphaLik=0.1;
    % number of subsamples to extract
    nsamp=1000;
    addnoisevariable = false;

    if addnoisevariable ==true
        X=[X randn(n,1)];
    end

    % In this case we use CWM model
    alphaXcwm=1;
    [outCWM]=tclustreg(y,X,k,restrfact,alphaLik,alphaXcwm,...
        'mixt',0,'plots',1,'msg',0,'nsamp',nsamp);
    title('CWM model')

    % In this case we use an adaptive second level of trimming.
    alphaXada=0.95;
    [outADA]=tclustreg(y,X,k,restrfact,alphaLik,alphaXada,...
        'mixt',0,'plots',1,'msg',0,'nsamp',nsamp);
    title('Adaptive second level trimming')

    % In this case we use a fixed level of trimming
    alphaXfixed=0.02;
    % Return in matrix X subsets used.
    [out,C]=tclustreg(y,X,k,restrfact,alphaLik,alphaXfixed,...
        'mixt',0,'plots',1,'msg',0,'nsamp',nsamp);
    title('Fixed second level trimming')
    disp('CWM model and adaptive second level trimming')
    disp('can recover the real structure of the data')
%}

%% Beginning of code
% Control variables, tolerances and internal flags
% warning('off');

% verLess2016b is true if current version is smaller than 2016b
verLess2016b=verLessThanFS(9.1);
%% Input parameters checking

nnargin   = nargin;
vvarargin = varargin;

[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% Check the presence of the intercept in matrix X
if min(max(X)-min(X))==0
    intercept = 1;
else
    intercept = 0;
end

% check restrfact option
if nargin < 4 || isempty(restrfact) || ~isnumeric(restrfact)
    restrfact = 12;         % This is the default in R
elseif min(restrfact)<1
    disp('Restriction factor smaller than 1. It is set to 1 (maximum constraint==>spherical groups)');
    restrfact(restrfact<1)=1;
else
end

% Check if restrfact is a scalar or a vector or length 2 (i.e. if the
% restriction is applied also on the explanatory variables)
if length(restrfact)>1
    % restriction factor among eigenvalues of covariance matrices of
    % explanatory variables
    restrfactX=restrfact(2);
    restrfact=restrfact(1);
else
    % No restriction of scatter matrices of explanatory variables
    restrfactX=Inf;
end

% checks on alpha1 (alphaLik) and alpha2 (alphaX)
if alphaLik < 0
    error('FSDA:tclustreg:error','error must a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
else
    % h is the number of observations not to be trimmed (used for fitting)
    if alphaLik < 1
        h = floor(n*(1-alphaLik));
    else
        h = n - floor(alphaLik);
    end
end

% checks on cwm, which decides if clusterwise regression has to be used
if alphaX < 0 || alphaX >1
    error('FSDA:tclustreg:WrongAlphaX','alphaX must a scalar in the interval [0 1]')
elseif alphaX==1
    cwm=1;
else
    cwm=0;
end

%% *Bivariate thinning* (_if wtrim == 4 and p == 2_)

% Bivariate thinning is applied once on the full dataset, at the start.
% This is done before setting the number of random samples nsamp.

if nargin>6
    % check if wtrim is among the user parameters
    chknwtrim = strcmp(varargin,'wtrim');
    if sum(chknwtrim)>0
        tmp = cell2mat(varargin(find(chknwtrim)+1));
        if ~isstruct(tmp)
            wtrimdef = 0;
            if cell2mat(varargin(find(chknwtrim)+1)) == 4
                interc = find(max(X,[],1)-min(X,[],1) == 0);
                if p == numel(interc) + 1
                    
                    % The bandwidth is chosen following Baddeley, as in the R
                    % spatstats package (density.ppp function).
                    bw = (range([X(:,numel(interc)+1),y]))/8;
                    
                    %in order to reproduce results comparable with the paper
                    %Cerioli and Perrotta (2013) the bandwidth is divided by 3
                    bw = bw/3;
                    % Another option to be considered follwing Baddeley is:
                    %bw = min(max(X),max(y))/8;
                    
                    % Thinning step
                    [Wt4,~] = wthin([X(:,numel(interc)+1),y], 'retainby','comp2one','bandwidth',bw);
                    id_unthinned = Wt4==1;
                    %id_thinned = Wt4==0;
                    
                    % save original data
                    Xori = X;
                    yori = y;
                    
                    % set retained data
                    X    = X(Wt4,:);
                    y    = y(Wt4);
                    
                    %recompute n on the retained data
                    n = size(y,1);
                end
            end
        else
            wtrimdef = struct;
        end
    else
        %no_wtrim = 1;
        wtrimdef = 0;
    end
    
    %% User options and their default values
    
    %%% - nsamp: the number of subsets to extract randomly, or the indexes of the initial subsets pre-specified by the User
    % Check whether option nsamp exists
    chknsamp = strcmp(varargin,'nsamp');
    if sum(chknsamp)>0
        nsamp=cell2mat(varargin(find(chknsamp)+1));
        
        % Check if options nsamp is a scalar
        if ~isscalar(nsamp)
            % if nsamp is not a scalar, it is a matrix containing in the
            % rows the indexes of the subsets which have to be extracted
            C=nsamp;
            [nsampdef,ncolC]=size(C);
            if ncolC ~= k*p
                disp('p is the total number of explanatory variables (including the constant if present)')
                error('FSDA:tclustreg:WrongInput','Input matrix C must contain k*p columns')
            end
            % The number of rows of nsamp (matrix C) is the number of
            % subsets which have to be extracted
            nselected=nsampdef;
            
            % Flag indicating if the user has selected a prior subset
            NoPriorSubsets=0;
            
            % In case of tandem thinning (Wtrim=4), the initial subset C
            % pre-specified by the user using the nsamp option might
            % include thinned units. In this case, we replace in C such
            % units with others that are close in terms of euclidean
            % distance. Verify the contiguity between the original and
            % replaced units with:
            %{
             figure;plot(Xori,yori,'.'); hold on ; text(Xori(nsamp(:)),yori(nsamp(:)),'X');
             figure;plot(Xori,yori,'.'); hold on ; text(X(C(:)),y(C(:)),'X');
            %}
            if sum(chknwtrim)>0 && ~isstruct(cell2mat(varargin(find(chknwtrim)+1)))
                if cell2mat(varargin(find(chknwtrim)+1))== 4
                    for f=1:size(C,1)*size(C,2)
                        if Wt4(C(f)) == 0
                            [~,C(f)] = min(pdist2([yori(C(f)),Xori(C(f))],[yori(Wt4),Xori(Wt4)]));
                        else
                            C(f) = sum(Wt4(1:C(f)));
                        end
                    end
                end
            end
            
        else
            % If nsamp is a scalar it simply contains the number of subsets
            % which have to be extracted. In this case NoPriorSubsets=1
            NoPriorSubsets=1;
        end
    else
        % If option nsamp is not supplied, then there are no prior subsets
        NoPriorSubsets=1;
    end
else
    % if nargin == 6, then the user has not supplied prior subsets
    NoPriorSubsets=1;
    wtrimdef = 0;
end

% If the user has not specified prior subsets (nsamp is not a scalar), then
% set the default number of samples to extract
if NoPriorSubsets == 1
    ncomb=bc(n,k*(p+intercept));
    nsampdef=min(300,ncomb);
end

%%% - Other user options

% default number of concentration steps
refstepsdef  = 10;

% default tolerance for comparing the classifications in two subsequent
% concentration steps
reftoldef=1e-5;

% default value for we: the observation weights
wedef = ones(n,1);

% default model: classification (mixt=0) or mixture likelihood (mixt=2)
mixtdef = 0;

% default choice for equalweight constraint
equalweightsdef = 0;

%seqk = sequence from 1 to the number of groups
seqk = 1:k;

%cup = pdf upper limit
cupdef = 1;

%pstar = thinning probability
pstardef = 1;

% commonslopedef = equal or different regression coefficients (excluding
% intercepts)
commonslopedef=false;

% automatic extraction of user options
options = struct('intercept',1,'mixt',mixtdef,...
    'nsamp',nsampdef,'refsteps',refstepsdef,...
    'reftol',reftoldef,'commonslope',commonslopedef,...
    'we',wedef,'wtrim',wtrimdef,...
    'equalweights',equalweightsdef,...
    'RandNumbForNini','','msg',1,'plots',1,...
    'nocheck',1,'k_dens_mixt',5,'nsamp_dens_mixt',nsampdef,...
    'refsteps_dens_mixt',refstepsdef,'method_dens_mixt','P','cup',cupdef,'pstar',pstardef);

if nargin > 6
    UserOptions = varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:tclustreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if all the specified optional arguments were present in
        % structure options. Remark: the nocheck option has already been dealt
        % by routine chkinputR.
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:tclustreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    % Write in structure 'options' the options chosen by the user
    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
    
    % Check if the number of subsamples to extract is reasonable
    if isscalar(options.nsamp) && options.nsamp>ncomb
        disp('Number of subsets to extract greater than (n k). It is set to (n k)');
        options.nsamp=0;
    elseif  options.nsamp<0
        error('FSDA:tclustreg:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
    end
end

% global variable controlling if messages are displayed in the console.
msg = options.msg;

% Graphs summarizing the results
plots = options.plots;

% Number of subsets to extract or matrix containing the subsets
nsamp = options.nsamp;

% Concentration steps
refsteps = options.refsteps;
reftol   = options.reftol;

% Common regression coefficients (excluding the intercepts)
commonslope=options.commonslope;

% Equalweights constraints
equalweights = options.equalweights;

% application-specific weights vector assigned by the user for beta
% estimation
we         = options.we;

%cup = pdf upper limit
cup = options.cup;

%pstar = thinning probability
pstar = options.pstar;

% Flag to control the type of thinning scheme for estimating beta
% (wtype_beta) and to compute obj function (wtype_obj)
if isstruct(options.wtrim)
    % Flag to control the type of thinning scheme for beta estimation
    wtype_beta      = options.wtrim.wtype_beta;
    % Flag to control the type of thinning scheme for obj function
    wtype_obj       = options.wtrim.wtype_obj;
else
    % if options.wtrim is a double it referes only to the beta estimation. No
    % weighting will be done in the obj function.
    wtype_beta      = options.wtrim;
    wtype_obj       ='0';
end
% Flag associated to the strategy for choosing the best refining step
% In the standard TCLUST the best refining step is granted to be the last
% one, because the objective funcion is monothonic. However, with second
% trimming level or componentwise thinning, the objective function may not
% be monothonic and a different strategy for choosing the best refining
% step can be considered.

zigzag = (alphaX > 0 && alphaX<1) || wtype_beta == 3 || wtype_beta == 2 || ~strcmp(wtype_obj, '0');

% Mixt option: type of membership of the observations to the sub-populations
% Control the mixture model to use (classification/mixture, likelihood or a
% combination of both):
%
% * mixt = 0: Classification likelihood
% * mixt = 1: Mixture likelihood, with crisp assignement
% * mixt = 2: Mixture likelihood
%
% $$ \prod_{j=1}^k  \prod_{i\in R_j} \phi (x_i;\theta_j) $$ $$ \quad $$
% $$ \prod_{j=1}^k  \prod_{i\in R_j} \pi_j \phi(x_i;\theta_j) $$ $$ \quad $$
% $$ \prod_{i=1}^n \left[ \sum_{j=1}^k \pi_j \phi (x_i;\theta_j)  \right] $$
mixt       = options.mixt;

if msg == 1
    switch mixt
        case 0
            % Classification likelihood.
            % To select the h untrimmed units, each unit is assigned to a
            % group and then we take the h best maxima
            disp('ClaLik with untrimmed units selected using crisp criterion');
        case 1
            % Mixture likelihood.
            % To select the h untrimmed units, each unit is assigned to a
            % group and then we take the h best maxima
            disp('MixLik with untrimmed units selected using crisp criterion');
        case 2
            % Mixture likelihood.
            % To select the h untrimmed units we take those with h largest
            % contributions to the likelihood
            disp('MixLik with untrimmed units selected using h largest lik contributions');
    end
end

% Initial mixing proportions $\pi_j$ can be user-defined or randomly generated
RandNumbForNini=options.RandNumbForNini;
if isempty(RandNumbForNini)
    NoPriorNini=1;
else
    NoPriorNini=0;
end

%% Initializations

%%% - Observation weights (we)

% Initialize we, a vector of application-specific weights associated to the
% observations, according to the trimming/thinning strategy.

switch wtype_beta
    case 0
        % standard case: weights must be all equal to 1. If user specifies
        % his weights vector differently, then the vector is rectified.
        if sum(we ~= wedef)>0
            disp('Warning: when "wtrim" is 0, "we" is set to a vector of ones');
            disp('         to give equal weights to all observations;');
            disp('         your vector "we" will not be considered.');
            we = wedef;
        end
    case 1
        % User specified weights: must be a column vector
        if sum(we == wedef)==n
            disp('Warning: when "wtrim" is 1, trimming is done by weighting');
            disp('         the observations using values specified in vector "we";');
            disp('         you left "we" to the default (i.e. a vector of ones,');
            disp('         giving equal weights to all observations);');
            disp('         we set them to a vector of 1/n, to sum to 1.');
        end
        % weights must be positive; if negative, values are translated so
        % that the minimum is 0
        if sum(we<0)>0
            we = we - min(we);
            disp('Warning: one or more of your weights are negative;');
            disp('         we added the minimum to all weights.');
        end
        % weights cannot be all equal to 0.
        if max(we) == 0
            we = wedef/n;
            disp('Warning: your weights are all zero;');
            disp('         we set them to a vector of 1/n, to sum to 1.');
        end
        
        % weights must be normalized so that to sum to 1
        we = we/nansum(we);
    case 2
        % weights will contain the componentwise univariate density probabilities.
        if sum(we ~= wedef)>0
            disp('Warning: when "wtrim" is 2, trimming is done by weighting');
            disp('         the observations according to the data density estimate;');
            disp('         your vector "we" will not be considered.');
            we = wedef;
        end
        
    case 3
        % weights will take {0,1} values, determined by componentwise univariate density probabilities.
        if sum(we ~= wedef)>0
            disp('Warning: when "wtrim" is 3, trimming is done by weighting');
            disp('         the observations with a Bernoulli random vector,');
            disp('         with probability of success depending on the data density estimate;');
            disp('         your vector "we" will not be considered.');
            we = wedef;
        end
        
    case 4
        % weights will take {0,1} values, determined by bivariate density probabilities.
        if length(we) ~= length(wedef)
            we = we(Wt4);
        end
        if sum(we ~= wedef)>0
            disp('Warning: when "wtrim" is 4, tclust is applied after thinning');
            disp('         observations with a Bernoulli random vector,');
            disp('         with probability of success depending on the data density estimate;');
            disp('         your vector "we" will not be considered.');
            we = wedef;
        end
end
%%% - Subsets extraction

%case with no prior subsets from the User
if NoPriorSubsets
    
    [C,nselected]=subsets(nsamp, n, (p+intercept)*k, ncomb, 0, we);
    
    % C = matrix which contains the indexes of the subsets to extract
    % for the k groups
    
    %nselected is set equal to the number of subsets:
    % - nselected = nsamp, if nsamp is a scalar;
    % - nselected = size(nsamp,1), if nsamp is a matrix (containing the initial subsets)
    
end

%%% - Output structures

% Store the initial subsets indices C
if nargout==2
    varargout={C};
end

% ll      = loglikelihood for each observations in each group
ll        = zeros(n,k);
% sigma2ini= variances deviation of each group
% sigma2ini = ones(1,k);


internationaltrade=false;
if internationaltrade==true
    wei=X(:,end).^2/var(X(:,end))+y.^2/var(y);
    weiForLikComputation=wei/max(wei);
else
    weiForLikComputation=1;
end

%%  RANDOM STARTS
% [bopt,sigma2opt,nopt,postprobopt,muXopt,sigmaXopt,cstepopt,vopt,subsetopt,idxopt,webeta,webetaopt]...
%     =tclustregcore(y,X,RandNumbForNini,reftol,refsteps,mixt,equalweights,h,nselected,k,restrfact,restrfactX,alphaLik,alphaX,...
%     seqk,zigzag,NoPriorNini,sigma2ini,msg,C,intercept,cwm,wtype_beta,we,wtype_obj);
[bopt,sigma2opt,nopt,postprobopt,muXopt,sigmaXopt,vopt,subsetopt,idxopt,webeta,cstepopt,webetaopt, Beta_all, obj_all] ...
    =tclustregcore(y,X,RandNumbForNini,reftol,refsteps,mixt,equalweights,h,nselected,k,restrfact,restrfactX,alphaLik,alphaX,...
    seqk,NoPriorNini,msg,C,intercept,cwm,commonslope,wtype_beta,we,wtype_obj,zigzag,weiForLikComputation,cup,pstar);

%%  END OF RANDOM STARTS



if isnan(sigma2opt)
    warning('FSDA:tclustreg:nocnvg','No convergence inside tclustreg')
    out=nan;
else
    %% Apply consistency factor based on the variance of the truncated normal distribution.
    
    % hh = number of non trimmed observations, after first and second level trimming
    hh = sum(nopt);
    
    % vt = variance of the truncated normal distribution
    % 1-hh/n is the trimming percentage
    vt = norminv(0.5*(1+hh/n));
    
    if hh<n
        factor = 1/sqrt(1-2*(n/hh)*vt.*normpdf(vt));
        % Apply the asymptotic consistency factor to the preliminary squared scale estimate
        sigma2opt_corr = sigma2opt*factor;
        % Apply small sample correction factor of Pison et al.
        sigma2opt_corr = sigma2opt_corr*corfactorRAW(1,n,hh/n);
    else
        sigma2opt_corr = sigma2opt;
    end
    
    %%  Set the output structure
    
    out                     = struct;
    out.class               = 'tclustreg';
    %cstepopt =             = cstep with maximum obj
    out.cstepopt            = cstepopt;
    out.subsetopt           = subsetopt;
    %   bopt                = regression parameters
    out.bopt                = bopt;
    %   sigmaopt0           = estimated group variances
    out.sigma2opt           = sigma2opt;
    %   sigma2opt_corr      = estimated group variances corrected with  asymptotic
    %                         consistency factor and small sample correction factor
    out.sigma2opt_corr      = sigma2opt_corr;
    %    out.wbetaopt            = webetaopt;
    
    if wtype_beta==5
        %TO BE IMPLEMENTED
        %out.w4trimopt_ovj_5 = w4trimopt_obj_5;
    end
    
    %CWM
    if cwm==1
        %out.muXopt= k-by-p matrix containing cluster centroid locations.
        out.muXopt           = muXopt';
        %out.sigmaXopt= p-by-p-by-k array containing estimated constrained
        %covariance covariance matrices of the explanatory variables for the k
        %groups.
        out.sigmaXopt        = sigmaXopt;
    end
    
    %   obj = value of the target function
    out.obj                  = vopt;
    
    if ~isempty(Beta_all)
        out.obj_all          = obj_all;
        out.Beta_all         = Beta_all;
    end
    %out.retained_id_all      = retained_id_all;
    
    out.C=C;
    
    %in tandem thinning it is necessary to save information about retained
    %units (idxopt) as well as thinned units.
    if wtype_beta == 4
        out.idx               = zeros(length(Xori),1);
        out.idx(id_unthinned) = idxopt;
    else
        out.idx               = idxopt;
    end
    % frequency distribution of the allocations
    out.siz=tabulateFS(idxopt(:,1));
    
    %postprobopt = posterior probabilities in the optimal cstep
    out.postprobopt     = postprobopt;
    
    % Store the indices in varargout
    if nargout==2
        varargout           = {C};
    end
    
    %% Compute INFORMATION CRITERIA
    
    %%% - Find NParam penalty term to use inside AIC and BIC
    
    % Add parameters referred to sigma2 restrictions
    % parameters associated to beta coefficients
    % npar=(p+(n-hh))*k;
    if commonslope == false
        npar=p*k;
    else
        npar=k+(p-intercept);
    end
    
    if equalweights==false      %to be generalized for equalweights==true
        npar=npar +(k-1);
    end
    nParam=npar+ (k-1)*(1-1/restrfact) +1;
    
    if cwm==1
        p1=p-intercept;
        nParam=nParam+ 0.5*p1*(p1-1)*k + (p1*k-1)*(1-1/restrfactX) +1;
    end
    
    % Specify which sigma2 to use in the final assignment
    % The two possibilities are siga2opt or sigma2opt_corr
    sigma2opt_corr=sigma2opt;
    
    % Discriminant functions for the assignments (use values of sigma2
    % corrected with Tallis
    if equalweights == 1
        for jj = 1:k
            ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*bopt(:,jj),0,sigma2opt_corr(jj));
            if cwm==1
                ll(:,jj)=  ll(:,jj)+ logmvnpdfFS(X(:,(intercept+1):end),muXopt(jj,:),sigmaXopt(:,:,jj));
            end
        end
    else
        for jj = 1:k
            ll(:,jj) = log((nopt(jj)/sum(nopt))) + logmvnpdfFS(y-X*bopt(:,jj),0,sigma2opt_corr(jj));
            if cwm==1
                ll(:,jj)=  ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muXopt(jj,:),sigmaXopt(:,:,jj));
            end
        end
    end
    
    % Now remove the rows which refer to first, or second level trimmed units
    % or thinned units
    
    delunits=false(n,1);
    delunits(idxopt(:,end)<0)=true;
    delunits(webeta==0)=true;
    
    ll(delunits,:)=[];
    
    if mixt>=1
        [NlogLmixt]=estepFS(ll,verLess2016b);
        % NlogLmixt is the negative of the maximized MIXTURE LOG-LIKELIHOOD
        % Note that if there was convergence NlogLmixt should be exactly equal to
        % -vopt
        NlogLmixt = -NlogLmixt;
    end
    
    loglik= max(ll,[],2);
    
    
    % NlogL is the negative of the CLASSIFICATION LOG-LIKELIHOOD  of the
    % untrimmed units
    % NlogL=-sum(max(ll(untrimmed units,[],2));
    % Note that if there was convergence NlogL should be exactly equal to
    % -vopt
    NlogL =-sum(loglik);
    
    
    
    logn=log(n);
    out.h=h;
    
    if mixt>0
        % MIXMIX = BIC which uses parameters estimated using the mixture loglikelihood
        % and the maximized mixture likelihood as goodness of fit measure (New BIC)
        MIXMIX  = 2*NlogLmixt +nParam*logn;
        
        % MIXCLA = BIC which uses the classification likelihood based on
        % parameters estimated using the mixture likelihood (New ICL)
        MIXCLA  = 2*NlogL +nParam*logn;
        
        out.MIXMIX=MIXMIX;
        out.MIXCLA=MIXCLA;
        
        % Store 2 times negative log likelihood for mixt and classification
        out.NlogL=2*NlogL;
        out.NlogLmixt=2*NlogLmixt;
    else
        % CLACLA = BIC which uses parameters estimated using the classification
        % likelihood and the maximized classification likelihood as goodness of fit
        % measure (New New)
        CLACLA  = 2*NlogL +nParam*logn;
        
        out.CLACLA=CLACLA;
        
        % Store 2 times negative log likelihood for mixt and classification
        out.NlogL=NlogL;
        out.NlogLmixt=[];
        
    end
    
    
    %% Generate plots
    
    if plots
        
        plot_type = out.class;
        
        % this is just for rotating colors in the plots
        clrdef = 'bkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcr';
        symdef = '+*d^v><phos+*d^v><phos+*d^v><phos+*d^v><phos';
        
        % The following plots are for the bi-variate case (i.e. v=1)
        if p-intercept < 2
            
            % initialize figure
            fh = figure('Name',[plot_type , ' plot'],'NumberTitle','off','Visible','on');
            gca(fh);
            hold on;
            
            for jj = 1:k
                group_label = ['Group ' num2str(jj)];
                
                % plot of the good units allocated to the current group.
                % Indices are taken after the second level trimming.
                % Trimmed points are not plotted by group.
                if wtype_beta==3
                    ucg = find(out.idx(:,end)==jj & webetaopt > 0);
                elseif wtype_beta==4
                    ucg = find(idxopt==jj);
                else
                    ucg = find(out.idx(:,1)==jj);
                end
                % misteriously text does not show a legend. This is why
                % we add a (ficticious) plot instruction with white symbols
                plot(X(ucg,end),y(ucg),'.w','DisplayName',[group_label ' (' num2str(length(ucg)) ' units)']);
                text(X(ucg,end),y(ucg),num2str(jj*ones(length(ucg),1)),...
                    'DisplayName',[group_label ' (' num2str(length(ucg)) ' units)'] , ...
                    'HorizontalAlignment','center','VerticalAlignment','middle',...
                    'Color',clrdef(jj), 'fontsize' , 12);
                
                % plot regression lines
                vv = [min(X(:,end)) max(X(:,end))];
                if intercept==1
                    plot(vv,bopt(1,jj)+bopt(2,jj)*vv,'DisplayName',[group_label ' fit' ],...
                        'Color',clrdef(jj));
                elseif intercept==0
                    plot(vv,bopt(:,jj)*vv,'DisplayName',[group_label ' fit' ],...
                        'Color',clrdef(jj));
                end
                
                if wtype_beta == 3
                    %plot the thinned (not trimmed) units
                    if jj == k
                        thinned_nt_trimmed = webetaopt;
                        thinned_nt_trimmed([ones(n,1) ;ones(n,1)]) = -12;
                        ucg = find(thinned_nt_trimmed == 0);
                        plot(X(ucg,end),y(ucg),symdef(jj),'color',clrdef(k+1),...
                            'DisplayName',['thinned units (' num2str(length(ucg)) ')' ]);
                    end
                end
            end
            
            % Plot the outliers (trimmed points)
            if wtype_beta==3
                ucg = find(out.idx(:,end)==-1);
            elseif wtype_beta==4
                ucg = find(idxopt==-1);
            else
                ucg = find(out.idx(:,1)==-1);
            end
            plot(X(ucg,end),y(ucg),'o','color','r','MarkerSize',8,...
                'DisplayName',['Trimmed units 1st (' num2str(length(ucg)) ')']);
            
            % Second level trimming points
            if wtype_beta==3
                ucg = find(out.idx(:,end)==-2);
            elseif wtype_beta==4
                ucg = find(idxopt==-2);
            else
                ucg = find(out.idx(:,1)==-2);
            end
            
            plot(X(ucg,end),y(ucg),'*','color','c',...
                'DisplayName',['Trimmed units 2nd (' num2str(length(ucg)) ')']);
            
            if wtype_beta == 4
                % in case of tandem thinning, plot the thinned points
                plot(Xori(~Wt4,end),yori(~Wt4),symdef(k+1),'color',clrdef(k+1),...
                    'DisplayName',['Thinned units (' num2str(length(Wt4) - sum(Wt4)) ')']);
            end
            
            % Position the legends and make them clickable. For some reason
            % clickableMultiLegend does not set properly the FontSize: to be fixed.
            lh=legend('show','Location','best');
            legstr = get(lh,'String');
            clickableMultiLegend(legstr,'Location','best','interpreter' , 'LaTex', 'fontsize' , 10);
            
            axis('manual');
            
            % control of the axis limits
            xmin = min(X(:,end)); xmax = max(X(:,end));
            ymin = min(y); ymax = max(y);
            deltax = (xmax - xmin) / 10;
            deltay = (ymax - ymin) / 10;
            
            xlim([xmin-deltax,xmax+deltax]);
            ylim([ymin-deltay,ymax+deltay]);
            
        else % In this case p > 2. A standard spmplot is used.
            
            if intercept
                YY = [X(:,2:end),y];
            else
                YY = [X,y];
            end
            
            % axis labels
            nameY = cellstr([repmat('X',size(YY,2)-1,1) , num2str((1:size(YY,2)-1)')]);
            nameY = [nameY ; 'y'];
            nameY = nameY';
            plo=struct;
            plo.nameY=nameY;
            plo.sym = [symdef(1:k) , 'o' ];
            plo.clr = [clrdef(1:k) , 'r' ];
            if sum(idxopt(:,end)==-2)>0
                plo.sym = [plo.sym , 's'];
                plo.clr = [plo.clr , 'r'];
            end
            
            % group names in the legend
            group = cell(n,1);
            group(idxopt(:,end)==-1) = {'Trimmed units'};
            group(idxopt(:,end)==-2) = {'Trimmed units level 2'};
            for iii = 1:k
                group(idxopt==iii) = {['Group ' num2str(iii)]};
            end
            
            % scatterplot
            spmplot(YY,group,plo,'hist');
            
        end
        
        % Title, reporting labels for beta coefficients
        betacoeff=sprintf('%0.3f; ',out.bopt(end,:));
        if wtype_beta ~= 3
            betacoeff=betacoeff(1:end-2); % remove the last
        end
        title({['$ wtrim_{beta}=' num2str(wtype_beta) ...
            '\quad wtrim_{obj}=' num2str(wtype_obj) ...
            '\quad mixt=' num2str(mixt) , ...
            '  \quad c='  num2str(restrfact) ...
            '\quad \alpha_{Lik}=' num2str(alphaLik) ...
            '\quad \alpha_X=' num2str(alphaX) '$'] , ...
            ['$ obj=' num2str(out.obj) '\quad b=(' betacoeff ') $']} , ...
            'interpreter' , 'LaTex', 'fontsize' , 14);
        
    end
end


%% _SUB-FUNCTIONS_

% corfactorRAW function
    function rawcorfac = corfactorRAW(p,n,alpha)
        
        if p > 2
            coeffqpkwad875=[-0.455179464070565,1.11192541278794,2;-0.294241208320834,1.09649329149811,3]';
            coeffqpkwad500=[-1.42764571687802,1.26263336932151,2;-1.06141115981725,1.28907991440387,3]';
            y1_500=1+(coeffqpkwad500(1,1)*1)/p^coeffqpkwad500(2,1);
            y2_500=1+(coeffqpkwad500(1,2)*1)/p^coeffqpkwad500(2,2);
            y1_875=1+(coeffqpkwad875(1,1)*1)/p^coeffqpkwad875(2,1);
            y2_875=1+(coeffqpkwad875(1,2)*1)/p^coeffqpkwad875(2,2);
            y1_500=log(1-y1_500);
            y2_500=log(1-y2_500);
            y_500=[y1_500;y2_500];
            A_500=[1,log(1/(coeffqpkwad500(3,1)*p^2));1,log(1/(coeffqpkwad500(3,2)*p^2))];
            coeffic_500=A_500\y_500;
            y1_875=log(1-y1_875);
            y2_875=log(1-y2_875);
            y_875=[y1_875;y2_875];
            A_875=[1,log(1/(coeffqpkwad875(3,1)*p^2));1,log(1/(coeffqpkwad875(3,2)*p^2))];
            coeffic_875=A_875\y_875;
            fp_500_n=1-(exp(coeffic_500(1))*1)/n^coeffic_500(2);
            fp_875_n=1-(exp(coeffic_875(1))*1)/n^coeffic_875(2);
        else
            if p == 2
                fp_500_n=1-(exp(0.673292623522027)*1)/n^0.691365864961895;
                fp_875_n=1-(exp(0.446537815635445)*1)/n^1.06690782995919;
            end
            if p == 1
                fp_500_n=1-(exp(0.262024211897096)*1)/n^0.604756680630497;
                fp_875_n=1-(exp(-0.351584646688712)*1)/n^1.01646567502486;
            end
        end
        if 0.5 <= alpha && alpha <= 0.875
            fp_alpha_n=fp_500_n+(fp_875_n-fp_500_n)/0.375*(alpha-0.5);
        end
        if 0.875 < alpha && alpha < 1
            fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
        end
        if alpha < 0.5
            fp_alpha_n = 1;
            if msg==1
                disp('Warning: problem in subfunction corfactorRAW')
                disp('alpha < 0.5')
            end
        end
        rawcorfac=1/fp_alpha_n;
        if rawcorfac <=0 || rawcorfac>50
            rawcorfac=1;
            if msg==1
                disp('Warning: problem in subfunction corfactorRAW')
                disp(['Correction factor for covariance matrix based on simulations found =' num2str(rawcorfac)])
                disp('Given that this value is clearly wrong we put it equal to 1 (no correction)')
                disp('This may happen when n is very small and p is large')
            end
        end
    end

end
%FScategory:CLUS-RobClaREG
