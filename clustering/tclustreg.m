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
%             matrices of explanatory variables. Scalar or vector with two
%             elements. If restrfact is a scalar it controls the
%             differences among group scatters of the residuals. The value
%             1 is the strongest restriction. If restrfactor is a vector
%             with two elements the first element controls the differences
%             among group scatters of the residuals and the second the
%             differences among covariance matrices of the explanatory
%             variables.
%            Data Types - single|double
%
%   alphaLik : Trimming level. Scalar.
%            alpha1 is a value between 0 and 0.5 or an  integer specifying
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
%            alphaX is a value between [0 and 1).
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
%               assumed equally distributed.
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
%equalweights : cluster weights in the concentration and assignment steps.
%               Logical. A logical value specifying whether cluster weights
%               shall be considered in the concentration, assignment steps
%               and computation of the likelihood.
%               if equalweights = true we are (ideally) assuming equally
%               sized groups by maximizing the likelihood
%                 Example - 'equalweights',true
%                 Data Types - Logical
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
% refsteps:  Number of refining iterations. Scalar. Number of refining
%               iterations in each subsample.  Default is 10.
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'refsteps',15
%                 Data Types - single | double
%    plots : Plot on the screen. Scalar. A flag to control the
%            generation of the plots.
%            If plots=1 a plot is showed on the screen with the
%            final allocation (and if size(X,2)==2 with the lines
%            associated to the groups)
%            Example - 'plots',1
%            Data Types - double
%   wtrim: Application of observation weights. Scalar. A flag taking values [0, 1, 2, 3, 4]
%          to control the application of weights on the observations.
%          -  If \texttt{wtrim}=0 (no weights) and \texttt{mixt}=0, the
%             algorithm reduces to the standard tclustreg algorithm.
%          -  If \texttt{wtrim}=0 and \texttt{mixt}=2, the maximum posterior
%             probability $D\_i$ of equation 7 of Garcia et al. 2010 is
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
%            Example - 'wtrim',1
%            Data Types - double
%      we: Vector of observation weights. Vector. A vector of size n-by-1
%          containing application-specific weights that the user needs to
%          apply to each observation. Default
%          value is  a vector of ones.
%            Example - 'we',[0.2 0.2 0.2 0.2 0.2]
%            Data Types - double
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
% RandNumbForNini: Pre-extracted random numbers to initialize proportions.
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
%  Output:
%
%  out :  structure containing the following fields
%
%   out.bopt             = $p \times k$ matrix containing the regression
%                          parameters.
%   out.sigma2opt       = $k$ row vector containing the estimated group
%                          variances.
%   out.sigma2opt_corr    = $k$ row vector containing the estimated group
%                          variances corrected with  asymptotic consistency
%                          factor.
%         out.muXopt= k-by-p matrix containing cluster centroid
%                       locations. Robust estimate of final centroids of
%                       the groups. This output is present just if input
%                       option alphaX is 1.
%         out.sigmaXopt= p-by-p-by-k array containing estimated constrained
%                       covariance covariance matrices of the explanatory
%                       variables for the k groups. This output is present
%                       just if input option alphaX is 1.
%            out.idx  = n-by-1 vector containing assignment of each unit to
%                       each of the k groups. Cluster names are integer
%                       numbers from -2 to k.
%                       -1 indicates first level trimmed units.
%                       -2 indicated second level trimmed units.
%            out.siz  = Matrix of size k-by-3.
%                       1st col = sequence from -2 to k;
%                       2nd col = number of observations in each cluster;
%                       3rd col = percentage of observations in each
%                       cluster;
%                       Remark: 0 denotes thinned units (if the weights are
%                       to find thinned units are 0 or 1, -1 indicates
%                       first level trimmed units and -2 indicates second
%                       level trimmed units).
% out.idx_before_tr = n-by-1 vector containing assignment of each unit to
%                       each of the k groups before applying first (and
%                       second level trimming). Cluster names are integer
%                       numbers from 1 to k. Note that while out.idx
%                       contains number which go from -2 to k,
%                       out.idx_before_tr only contains numbers which go
%                       from 1 to k.
%            out.post = n-by-k matrix containing posterior probabilities
%                       out.post(i,j) contains posterior probabilitiy of unit
%                       i from component (cluster) j. For the trimmed units
%                       posterior probabilities are 0.
%           out.vopt  = Scalar. The value of the target function.
%             out.we  = n-by-1 vector  containing the user-specific weigths
%                       of each observation, i.e. its contribution to the
%                       estimates.
%   out.postprobopt   = $n \times k$ matrix containing the final posterior
%                           probabilities. out.postprobopt(i,j) contains
%                           posterior probabilitiy of unit i from component
%                           (cluster) j. For the trimmed units posterior
%                           probabilities are 0.
%
%
%  Optional Output:
%
%            C     : Indexes of extracted subsamples. Matrix.
%                    Matrix of size nsamp-by-k*p containing (in the rows)
%                    the indices of the subsamples extracted for computing
%                    the estimate.
%
% More about:
%
%  Computational issues to be addressed in future releases.
%  1.
%  FSDA function "wthin" uses the MATLAB function ksdensity. The calls to
%  ksdensity have been optimized. The only possibility to further reduce
%  time execution is to replace ksdensity with a more efficient kernel
%  density estimator.
%  2.
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
%  3.
%  In the plots, the use of text to highlight the groups with their index
%  is terribly slow (more than 10 seconds to generate a scatter of 7000
%  units. ClickableMultiLegend and legend are also slow.
%  4.
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
% Copyright 2008-2018.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustreg')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
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

%}

%{
    % tclustreg of fishery data 1.
    clear all; close all;
    load fishery;
    X = fishery.data;
    % some jittering is necessary because duplicated units are not treated:
    % this needs to be addressed
    X = X + 10^(-8) * abs(randn(677,2));

    %tclustreg on fishery data
    y1 = X(:,end);
    X1 = X(:,1:end-1);
    k = 3 ; restrfact = 50; alpha1 = 0.04 ; alpha2 = 0.01;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',0);

    %lga and rlga on fishery data
    figure('name','RLGA');
    out=lga(X,3);
    clickableMultiLegend('1','2','3','data1','data2','data3');
    axis manual;

    alpha = 0.95;
    figure('name','LGA');
    out=rlga(X,3,1-alpha);
    clickableMultiLegend('0','1','2','3','data1','data2','data3');
    axis manual;
%}

%{
    % tclustreg of fishery data 2.
    clear all; close all;
    load fishery;
    X=fishery.data;
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

    k = 2 ;

    restrfact = 5; alpha1 = 0.05 ; alpha2 = 0.01;
    nsamp=200;
    Cnsamp=subsets(nsamp,n,(size(X1,2)+1)*k);
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'nsamp',Cnsamp);
%}


%% Internal Debugging Flags

warning('off');

%  csteps_stop == 0 is to monitor the stability of the objective
% function for a fixed number of loops (refsteps option). Otherwise, if it
% is set to 1, the loop stops when the classification does not change in 2
% consective c-steps.
csteps_stop = 0;

% assign_thinned_units == 1 is to assign thinned units to the most likely cluster.
assign_thinned_units = 1;

% Groups with less than skipthin_th units are not considered for thinning
skipthin_th = 50;

%% Initializations

% current and best objective function values
vopt = -1e+20;
obj = vopt;

% tolerance for restriction factor
tolrestreigen = 1e-08;

% index of the best concentration step
cstepopt = 0;

% use of repmat (from Release 8.2 repmat is faster than bsxfun)
if verLessThan('matlab','8.2.0') == 1
    userepmat=0;
else
    userepmat=1;
end

%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;

[~,p0] = size(X);
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% Intercept, yes/no
if p>p0
    intercept = 1;
else
    intercept = 0;
end

% check restrfact option
if nargin < 4 || isempty(restrfact) || ~isnumeric(restrfact)
    restrfact = 12;         % AGUSTIN: deafult 12 OK?
elseif min(restrfact)<1
    disp('Restriction factor smaller than 1. It is set to 1 (maximum contraint==>spherical groups)');
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


%% Bivariate thinning
%
% *If wtrim == 4 and p == 2*, bivariate thinning is applied
% once at the very beginning on the full dataset. This has to be done
% before setting the number of random samples nsamp.

if nargin>6
    chknwtrim = strcmp(varargin,'wtrim');
    if sum(chknwtrim)>0 && cell2mat(varargin(find(chknwtrim)+1)) == 4
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
end

%% Get user options: initial subsets pre-specified by the User.
%
% *nsamp* can contain the number of subsets to extract randomly, but also
% the indexes of the subsets which have to be extracted from the data
% matrix

if nargin>6
    
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
            if sum(chknwtrim)>0 && cell2mat(varargin(find(chknwtrim)+1)) == 4
                for f=1:size(C,1)*size(C,2)
                    if Wt4(C(f)) == 0
                        [~,C(f)] = min(pdist2([yori(C(f)),Xori(C(f))],[yori(Wt4),Xori(Wt4)]));
                    else
                        C(f) = sum(Wt4(1:C(f)));
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
end

% If the user has not specified prior subsets (nsamp is not a scalar), then
% set the default number of samples to extract
if NoPriorSubsets == 1
    ncomb=bc(n,k*(p+intercept));
    nsampdef=min(300,ncomb);
end


%% Get other user options

% default number of concentration steps
refstepsdef  = 10;

% default value for wtrim: the thinning strategy (wtrim can be 0,1,2,3,4)
wtrimdef = 0;

% default value for we: the observation weights
wedef = ones(n,1);

% default model: classification (mixt=0) or mixture likelihood (mixt=2)
mixtdef = 0;

% default choice for equalweight constraint
equalweightsdef = 1;

%seqk = sequence from 1 to the number of groups
seqk = 1:k;

% automatic extraction of user options
options = struct('intercept',1,'mixt',mixtdef,...
    'nsamp',nsampdef,'refsteps',refstepsdef,...
    'we',wedef,'wtrim',wtrimdef,...
    'equalweights',equalweightsdef,...
    'RandNumbForNini','','msg',1,'plots',1);

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

%% Set user's options

% global variable controlling if messages are displayed in the console.
msg = options.msg;

% Graphs summarizing the results
plots = options.plots;

% Number of subsets to extract or matrix containing the subsets
nsamp = options.nsamp;

% Concentration steps
refsteps = options.refsteps;

% Equalweights constraints
equalweights = options.equalweights;

% application-specific weights vector assigned by the user
we         = options.we;

% Flag to control the type of thinning scheme
wtrim      = options.wtrim;

%%% Flag associated to the strategy for choosing the best refining step
% In the standard TCLUST the best refining step is granted to be the last
% one, because the objective funcion is monothonic. However, with second
% trimming level or componentwise thinning, the objective function may not
% be monothonic and a different strategy for choosing the best refining
% step can be considered.

zigzag = (alphaX > 0 || wtrim == 3 || wtrim == 2);

%%% Mixt option: type of membership of the observations to the sub-populations
% Control the mixture model to use (classification/mixture, likelihood or a
% combination of both):
%
% * mixt = 0: Classification likelihood
% * mixt = 1: Mixture likelihood, with crisp assignement
% * mixt = 2: Mixture likelihood
%
% $$ \prod_{j=1}^k  \prod_{i\in R_j} \phi (x_i;\theta_j) $$ $$ \quad $$
% $$ \prod_{j=1}^k  \prod_{i\in R_j} \pi_j \phi(x_i;\mu_j,\Sigma_j) $$ $$ \quad $$
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


%% Initialize observation weights

% Initialize we, a vector of application-specific weights associated to the
% observations, according to the trimming/thinning strategy.
switch wtrim
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
        we = we(:);
        
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

%% Initial subsets extraction, if not pre-specified by the User.

%case with no prior subsets from the User
if NoPriorSubsets
    
    [C,nselected]=subsets(nsamp, n, (p+intercept)*k, ncomb, 0, we);
    
    % C = matrix which contains the indexes of the subsets to extract
    % for the k groups
    
    %nselected is set equal to the number of subsets:
    % - nselected = nsamp, if nsamp is a scalar;
    % - nselected = size(nsamp,1), if nsamp is a matrix (containing the initial subsets)
    
end

%% Initialize output structures

% timer to monitor performances.
tsampling = ceil(min(nselected/5 , 50));
time      = zeros(tsampling,1);

% Store the initial subsets indices C
if nargout==2
    varargout={C};
end

% ll      = loglikelihood for each observations in each group
ll        = zeros(n,k);
% sigma2ini= standard deviation of each group
sigma2ini = ones(1,k);
% bopt     = beta parameters of each group obtained in the best subset
bopt      = zeros(p,k);
% sigma2opt= sigma parameters of each group obtained in the best subset
sigma2opt = NaN(1,k);

%initializations for clusterwise regression
if cwm==1
    % muX = matrix of size (p-intercept)-by-k which contains centroids of
    % each extracted subset
    % sigmaX = matrix of size (p-intercept)-by-(p-intercept)-k which
    % contains covariance matrices of each extracted subset
    muX       = zeros(p-intercept,k);
    sigmaX    = zeros(p-intercept,p-intercept,k);
    
    % U = 3D array which contains the eigenvectors of covariance matrix of
    % the explanatory variables for each group
    U         = sigmaX;
    
    % Lambda_pk = matrix which contains the eigenvalues of sigmaX referred
    % to the gruops (first column is associated with group 1....)
    Lambda_pk = muX;
end

%w4trimopt = vector of {0,1} weights with 1 to identify units that are not
%thinned.  The vector is the one giving rise to the optimal solution.
%The intermediate vector, w4trim, is used in the first trimming step to
%compute the likelihood contribution on the units that are not thinned.
w4trimopt = ones(n,1);

%idxopt = (nx1) vector of values in {1, ... , k, -1 , -2} with cluster
%assignements (1,...,k) and trimmed units (-1 for first level trimming and
%-2 for second level trimming).
idxopt    = zeros(n,1);

%Remark. The assignement of the thinned units can be found with:
%A = [find(w4trimopt==0)' , idxopt(find(w4trimopt==0))']
%where the first column of A is the row index of the thinned unit and the
%second one is the cluster assignement.

% nopt = (kx1) vector containing the number of observations in each of
% the k groups in the optimal subset. The trimmed units are not counted.
nopt       = zeros(1,k);

% postprob = (nxk) matrix of posterior probabilities of the obtimal subset
postprobopt = zeros(n,k);

%obj_all: used to monitor the objective function during the refining steps.
obj_all = NaN(nselected,refsteps);
%indmax_all = NaN(n,refsteps,nselected);
%beta_all   = NaN(k,refsteps,nselected);

%%  RANDOM STARTS: beginning

for i =1:nselected
    
    switch msg
        case 1
            % monitor time execution
            if msg==1
                if i <= tsampling
                    tstart = tic;
                end
            end
        case 2
            % monitor iteration step
            disp(['Iteration ' num2str(i)]);
    end
    
    % ltkg = becomes 1 if a particular subset leads to less than k groups
    ltkg=0;
    
    % Beta = matrix of beta coefficients (j-col refers to j-th group)
    Beta = zeros(p,k);
    
    % to replicate results, fix the seed by uncommenting the line below
    % rng(1234);
    
    %%% -- Initialization of mixing proportions
    if NoPriorNini==1
        % if initial mixing proportions are not supplied by the user, they
        % are randomly assigned, making sure that:
        % a) minimum group size is positive
        % b) sum of group sizes is equal to h
        niin=0;
        itermax=0;
        while min(niin)==0 && itermax<1000
            randk=rand(k,1);
            niin=floor(h*randk/sum(randk));
            diffh=sum(niin)-h;
            [~,imin]=min(niin);
            niin(imin)=niin(imin)-diffh;
            itermax=itermax+1;
        end
        if itermax ==1000
            error('FSDA:tclustreg:WrongInput','Initialization of the group proportions failed')
        end
        niini=niin;
    else
        %check the initial mixing proportions supplied by the user
        randk=RandNumbForNini(:,i);
        itermax=0;
        
        % Make sure that the sum of niin is h
        niin=floor(h*randk/sum(randk));
        if sum(niin)<h
            niin(1)=niin(1)+h-sum(niin);
        end
        % Make sure that minimum group size is strictly positive
        while min(niin)==0 && itermax<1000
            ntoreplace=seqk(niin==0);
            niin(ntoreplace(1))=1;
            [~,posmax]=max(niin);
            niin(posmax)=niin(posmax)-1;
            itermax=itermax+1;
        end
        niini=niin;
    end
    
    %%% -- Initial regression parameters estimation
    
    % Estimate the regression parameters Beta and the Var-Cov matrix for
    % each of the k *initial* groups provided by the user or randomly
    % generated.
    
    %index: the i-th line of C, containing the i-th subset units
    index  = C(i,:);
    for j = 1:k
        %selj: the units of the current subset i, belonging to the group j
        ilow   = (j-1)*p+1;
        iup    = j*p;
        selj   = index(ilow:iup);
        
        %Xb and yb: X and y of the current subset i belonging to the group j
        Xb     = X(selj,:);
        yb     = y(selj,:);
        
        %If the model is without intercept and Xb is zero, the regression
        %cannot be computed. Some jittering is applied to Xb 
        if intercept == 0 && isscalar(Xb) && Xb == 0
            Xb = Xb + 0.0001*abs(randn(1,1));
        end
        
        %Beta and sigma2ini: estimation of regression parameters
        Beta(:,j) = Xb\yb;
        if length(yb)==1
            sigma2ini(j)= var(y);
        else
            sigma2ini(j) =var(yb);
        end
        
        % clusterwise regression: initialize X-distribution parameters
        if cwm==1
            
            % muX = matrix of size (p-intercept)-by-k which contains the k
            % centroids of the current subset i
            % sigmaX = matrix of size (p-intercept)-by-(p-intercept)-k which
            % contains covariance matrices of the current subset i
            muX(:,j)=mean(Xb(:,intercept+1:end),1);
            if length(yb)==1
                sigmaX(:,:,j)=var(X);
            else
                sigmaX(:,:,j)=cov(Xb(:,intercept+1:end));
            end
            
            % Lambdaj and Uj: eigenvalue eigenvector decomposition in the
            % current subset i, for group j
            [Uj,Lambdaj] = eig(sigmaX(:,:,j));
            % Store eigenvectors and eigenvalues of group j
            U(:,:,j)=Uj;
            Lambda_pk(:,j)=diag(Lambdaj);
        end
        
    end
    
    %%% -- Eigenvector-eigenvalue restriction
    % $$ \frac{ max_{g=1,\ldots,G} \sigma_g^2 \pi_g}{min_{g=1,\ldots,G} \sigma_g^2 \pi_g} < restrfact $$

    %restrict the eigenvalues according to the constraint specified in restrfact: 
    if equalweights==1
        sigma2ini= restreigen(sigma2ini,ones(k,1),restrfact,tolrestreigen,userepmat);
        if cwm==1
            autovalues= restreigen(Lambda_pk,ones(k,1),restrfactX,tolrestreigen,userepmat);
        end
    else
        sigma2ini= restreigen(sigma2ini,niini,restrfact,tolrestreigen,userepmat);
        if cwm==1
            autovalues= restreigen(Lambda_pk,niini,restrfactX,tolrestreigen,userepmat);
        end
    end
    
    % CWM: re-construct the covariance matrices keeping into account the
    % constraints on the eigenvalues
    if cwm==1
        for j=1:k
            sigmaX(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
            % Alternative code: more elegant but slower because diag is a
            % built in function
            % sigmaX(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
        end
    end
    
    %%% -- Log-likelihood of all observations based on the estimated regression parameters
    % equalweight: each group has the same weight, $1/k$.
    if equalweights == 1
        for jj = 1:k
            ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
            if cwm==1
                ll(:,jj)=ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(:,jj),sigmaX(:,:,jj));
            end
        end
    % the group weights (niini) are estimated
    else
        for jj = 1:k
            ll(:,jj) = log((niini(jj)/sum(niini))) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
            if cwm==1
                ll(:,jj)=ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(:,jj),sigmaX(:,:,jj));
            end
        end
    end
    
    %%% -- Estep: posterior probabilities of all observations
    
    % Mixture likelihood model
    if mixt > 0
        %E-step is run to compute the posterior probabilities of all
        %observations
        [~,postprob,~] = estepFS(ll);
        % idx: (nx1) vector indicating the group for which each observation
        % has the largest posterior probability. It takes values in  {1,
        % ... , k}; At the end it will take values in {1, ... , k,0, -1 ,
        % -2}, respectively for group assignement, thinned units, first and
        % second trimmed units.
        [~,idx]= max(postprob,[],2);
        
    %classification likelihood model
    else %  mixt == 0
        zeronk=zeros(n,k);
        % idx: (nx1) vector indicating the group for which each observation
        % has the largest likelihood. It takes values in  {1, ... , k};
        % At the end it will take values in {1, ... , k,0, -1 ,
        % -2}, respectively for group assignement, thinned units, first and
        % second trimmed units.
        [~,idx] = max(ll,[],2);
        postprob = zeronk;
        for j=1:k
            postprob(idx==j,j)=1;
        end
    end
    
    
    %% -- Concentration Steps
    
    indold = zeros(n,1)-1;
    for cstep = 1:refsteps
        
        %%% -- -- Componentwise Thinning (the wtrim option)
        
        % w4trim = vector of weights {0,1}; 1 identifies units that were not
        % thinned. In the first trimming step, it is used to compute the
        % likelihood contribution on the units that are potentially subject
        % to trimming.
        
        switch wtrim
            
            case 0
                % no observation weighting: w4trim is constant
                w4trim = ones(n,1);
                
            case 1
                % w4trim contains the user-specific weights
                w4trim = we;
                
            case 2
                % w4trim contains the retention probabilities on yhat,
                % estimated with function wthin 

                % small_group: kx1 vector identifying groups that are too
                % small to be thinned.
                % small_group_obs: nx1 matrix identifying observations in
                % small groups.
                small_group     = zeros(k,1);
                small_group_obs = zeros(n,1);
                
                w4trim = ones(n,1);
                for jj=1:k
                    % Boolean index of units forming group j
                    groupj=idx==jj;
                    
                    % ijj: indices of units in group jj
                    ijj = find(groupj);
                    
                    % update the weight vector w4trim: elements with
                    % indices of units in group jj are updated only if the
                    % group has more than skipthin_th units
                    if  length(ijj) > skipthin_th 
                        % pretain: the retention probabilities are based on
                        % the predicted values (yhat) estimated at the
                        % previous concentration step. REMARK: trimmed and
                        % non-trimmed units are both considered.
                        Xj = X(ijj,:);
                        yhat = Xj*Beta(:,jj);
                        
                        [~ , pretain] = wthin(yhat);
                        w4trim(ijj)   = pretain;
                        
                    else 
                        % The group is too small: skip thinning, flag
                        % with 1 vector small_group and with jj
                        % small_group_obs. If the group is completely
                        % empty, it is not necessary to flag the two
                        % vectors.
                        if ~isempty(ijj)
                            small_group(jj) = 1;
                            small_group_obs(:) = jj*groupj;
                        end
                    end
                    
                end
                
                % for too-small (but grater than zero) groups, where
                % thinning is not possible, w4trim is the median of the
                % weights of the other groups.
                if sum(small_group) > 0
                    %compute median on observations of enoughout-large
                    %groups
                    medianweights = nanmedian(w4trim(ismember(idx,find(small_group==0)))) ;
                     
                    id_small_groups = find(small_group==1);
                    for jj=1:sum(small_group)
                        w4trim(small_group_obs == id_small_groups(jj))=medianweights;
                        %it can happend that a group is completely empty
 %                       rep_median_we = repmat(medianweights,n_obs_current_gr,1);
 %                       w4trim(small_group_obs(:,pos_check_we_groups(jj)) == 1) = rep_median_we;
                    end
                end
                
            case 3
                
                % weights are the posterior probabilities multiplied by
                % the bernoulli weights
                
                ii = 0;
                w4trim=ones(n,1);
                for jj=1:k
                    % find indices of units in group jj
                    if cstep == 1
                        ijj = find(idx==jj);
                    else
                        ijj = find(idx_ne0==jj);
                    end
                    %ijj_ori = find(indall == jj);
                    % weight vector is updated only if the group has
                    % more than thinning_th observations abd if the
                    % beta of the group is not zero
                    if  numel(ijj)> skipthin_th
                        % Bernoulli weights based on density estimated
                        % on the component predicted values of the
                        % previous step. The retention probabilities
                        % (the second output argument of wthin, i.e.
                        % pretain) are not used.
                        Xj = X(ijj,:);
                        yhat = Xj*Beta(:,jj);
                        [Wt , ~] = wthin(yhat);
                        
                        %bw3 =  bwe(yhat,'robust');
                        %[Wt , ~] = wthin(yhat,'bandwidth',bw3);
                        
                        w4trim(ijj)=Wt;
                        
                        % count the thinned observations
                        nthinned = sum(Wt == 0);
                        ii = ii + nthinned;
                        %the next command should be descussed!!!!! FRT
                        niini(jj) = sum(Wt > 0);
                        %niini(jj) = niini(jj) - nthinned;
                        %eliminate thinned observations
                    end
                    
                end
                
            case 4
                % This is the bivariate thinning
                w4trim = Wt4;
        end
        
        % Mean of the weights must be 1
        % These weights will be used in first and second level trimming
        mean_w4trim = mean(w4trim);
        sum_w4trim = sum(w4trim);
        w4trim=w4trim/mean(w4trim);
        %%%%%%%%%%%%%%%%%%%END OF THINNING
        
        
        %not_empty_g = zeros(1,k);
        if sum(isnan(Beta(:)))>0
            break
        else
            
            
            %% 1) Log-likelihood (inside concentration steps)
            % Discriminant functions for the assignments
            if equalweights == 1
                for jj = 1:k
                    ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
                    if cwm==1
                        ll(:,jj)=  ll(:,jj)+ logmvnpdfFS(X(:,(intercept+1):end),muX(:,jj),sigmaX(:,:,jj));
                    end
                end
            else
                for jj = 1:k
                    ll(:,jj) = log((niini(jj)/sum(niini))) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
                    if cwm==1
                        ll(:,jj)=  ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(:,jj),sigmaX(:,:,jj));
                    end
                end
            end
            
            % idx is a (nx1) vector indicating in correspondence of which group
            % each observation has the largest posterior probability. It is a
            % vector of cluster assignements, taking values in  {1, ... , k};
            % At the end it will take values in {1, ... , k, 0, -1 , -2}.
            if mixt == 2
                %In the case of mixture modelling we need an estep to compute the log of the sum
                %of the log-likelihood (disc) and the normalized log-likelihood (postprob).
                [~,postprob,disc] = estepFS(ll);
                % indmax assigns each observation to the group with the largest posterior probability
                [~,idx]= max(postprob,[],2);
                
            else
                % we are in the case mixt=0 or mixt=1
                %We compute the maximum value of the log-likelihood
                %(disc) among the k groups, as in point 2.1 of appendix of Garcia Escudero et al.
                %(2010)
                
                
                [disc,idx] = max(ll,[],2);
                
                postprob=zeronk;
                for j=1:k
                    postprob(idx==j,j)=1;
                end
            end
            
            % Assign infinite weights to thinned units so that they can never be
            % trimmed
            disc(w4trim==0)=1e+20;
            
            
            if wtrim == 3
                % idx_ne0, computed only for wtrim = 3, is a (nx1) vector of
                % cluster assignements, taking values in  {1, ... , k}; tt the
                % end it will take values in {1, ... , k,  -1 , -2}.
                
                % idx is a (nx1) vector of cluster assignements, taking values
                % in  {1, ... , k, 0}; At the end it will take values in {1,
                % ... , k, 0, -1 , -2}.
                %idx could be always obtained as idx=idx_ne0.*w4trim
                
                idx_ne0=idx;
                idx(w4trim==0)=0;
            end
            % Sort the n likelihood contributions and save in qq the largest n*(1-alpha) likelihood
            % contributions
            [~,qq] = sort(disc,'ascend');
            
            
            %% 2) First level trimming (inside concentration steps)
            
            % Order the weights according to qq
            w4trimordered=w4trim(qq);
            
            % Find cumlative sum of weights
            cumsumww = cumsum(w4trimordered);
            
            % qqunassigned_small is a n-by-1 Boolean vector which
            % contains true in the first k positions (units which
            % have to be trimmed)
            qqunassigned_small = cumsumww <= n*alphaLik;
            
            % qqunassigned = indexes of units subject to first
            % level trimming
            qqunassigned = qq(qqunassigned_small);
            
            %indmax_before_tr should be saved in order to be able
            %to understand which groups the two trimmings affect more
            %indmax_before_tr = idxtt; TO DELETE
            
            
            % idx takes values in {1, ... , k, 0,  -1}. At the end it will
            % takes values in {1, ... , k, 0,  -1, -2}.
            
            idx(qqunassigned)=-1;
            if wtrim==3
                % idx_ne0 takes values in {1, ... , k,  -1}. At the end it
                % will takes values in {1, ... , k,  -1, -2}.
                idx_ne0(qqunassigned)=-1;
            end
            % Put equal to zero the posterior probability of unassigned
            % units
            postprob(qqunassigned,:)=0;
            
            % Now multiply posterior probabilities (n-by-k matrix) by
            % weights (n-by-1 vector) obtained using thinning
            if wtrim==3
                Z_all=postprob;
            end
            Z=bsxfun(@times,postprob,w4trim);
            
            %% 3) Second level of trimming (inside concentration steps) or cwm model
            % FS or MCD are used to find units to trim
            % Using untrimmed units find beta coefficients and
            % sigma2 using weighted regression (the weights are
            % based on thinning fixed once and for all before the
            % concentration steps)
            
            if alphaX<1
                for jj=1:k
                    % find indices of units belonging to groupj
                    groupjind = find(idx==jj);
                    
                    Xjnointercept  = X(groupjind,intercept+1:end);
                    Xj=X(groupjind,:);
                    njj=size(Xj,1);
                    
                    if alphaX>0.5
                        hj=njj;
                    else
                        % hjj = number of units to retan for group j
                        % This is simply h referred to group j
                        hj = floor(njj*(1-alphaX));
                    end
                    %check if a group is populated, that is check if size
                    %of group after second level trimming is not smaller
                    %than p and that number of units belonging to group
                    %before second level trimming has a size at least
                    %greater than p+2 (to run the FS)
                    if hj >=p && njj>p+2
                        
                        
                        if njj/(p-intercept)>10 && alphaX>0
                            
                            %The MCD is applied only when p=1, because in this case it is faster
                            %than the FS.
                            if p-intercept==1
                                
                                %R2016a has introduced robustcov, which could be used here as below.
                                %Remember however that mcd returns the squared distances, i.e. RAW.md = mah.^2.
                                %[~,~,mah,~,~] = robustcov(inliers,'Method','fmcd','NumTrials',nsampmcd,'OutlierFraction',alpha2b,'BiasCorrection',1); %
                                %plot(1:ni(jk),mah.^2,'-r',1:ni(jk),RAW.md,'-b');
                                %close;
                                if alphaX>0.5
                                    [~,REW]      = mcd(Xjnointercept,'msg',0,'conflev',1-(1-alphaX)/njj,'betathresh',1);
                                    if isfield(REW,'outliers')
                                        trimj=REW.outliers;
                                    else
                                        trimj = [];
                                    end
                                else
                                    [~,REW]      = mcd(Xjnointercept,'msg',0,'betathresh',1);
                                    
                                    if isfield(REW,'md')
                                        [~,indmdsor] = sort(REW.md);
                                        % Trimmed units (after second level
                                        % trimming)
                                        % indmax for second level trimmed units is
                                        % set to -1
                                        trimj=indmdsor(hj+1:end);
                                    else
                                        trimj = [];
                                    end
                                    
                                end
                            else
                                % Lines below are to find the second trimming points id_trim with
                                % the Forward Search rather than the MCD, using function FSMmmd
                                % Function FSMbsb is consderably faster than mcd when v>1.
                                % BBsel contains a NAN for the units
                                % not belonging to subset in step hj
                                if alphaX>0.5
                                    % 'init',round(njj*0.9)
                                    outj=FSM(Xjnointercept,'nocheck',1,'init',round(njj*0.9),'msg',0,'bonflev',alphaX,'plots',0);
                                    
                                    trimj=outj.outliers;
                                    % disp(trimj)
                                    if isnan(trimj)
                                        trimj=[];
                                    end
                                else
                                    [~,BBsel]=FSMbsb(Xjnointercept,0,'bsbsteps',hj,'init',hj,'nocheck',1,'msg',0);
                                    seqj=1:njj;
                                    trimj=seqj(isnan(BBsel));
                                end
                                
                            end
                        else
                            trimj=[];
                        end
                        
                        % idx takes values in {1, ... , k, 0,  -1,-2}.
                        idx(groupjind(trimj))=-2;
                        if wtrim==3
                            % idx_ne0 takes values in {1, ... , k, -1,-2}.
                            idx_ne0(groupjind(trimj))=-2;
                        end
                    else
                        ltkg = 1;
                        break
                    end
                    
                end
                
                % Set to zeros the rows of matrix Z which were trimmed based on
                % second level
                Z(idx==-2,:)=0;
                if wtrim==3
                    Z_all(idx==-2,:)=0;
                end
                % Set to zeros the rows of matrix postprob which were
                % trimmed based on second level
                postprob(idx==-2,:)=0;
                
            end
            
            if ltkg==1
                Beta=NaN;
                break
            end
            
            
            % Find beta coefficients and sigma2 using weighted
            % regression
            
            for jj=1:k
                
                %weights (for beta estimation) of observations
                %belonging to group iii, after second level trimming.
                sqweights = sqrt(Z(:,jj));
                % nj is the sum of the weights for gtoup j
                nj=sum(Z(:,jj));
                % ninini vector with the sum of the weights for each group
                niini(jj)=nj;
                
                if nj>p+1
                    Xw = bsxfun(@times, X, sqweights);
                    yw = y .* sqweights;
                    
                    % Estimate of beta from (re)weighted regression (RWLS)
                    breg = Xw\yw;
                    Beta(:,jj)=breg;
                    
                    % Find estimate of sigma2 after weighted regression
                    res2=(yw-Xw*breg).^2;
                    sigma2=sum(res2)/nj;
                    
                    sigma2ini(jj) = sigma2;
                    
                    if cwm ==1
                        muX(:,jj)=sum(bsxfun(@times, X(:,intercept+1:end), Z(:,jj)),1)/sum(Z(:,jj))';
                        sigmaX(:,:,jj)= (Xw(:,intercept+1:end)'*Xw(:,intercept+1:end))/sum(Z(:,jj))-muX(:,jj)*(muX(:,jj)');
                    end
                else
                    
                    ltkg=1;
                    break
                end
            end % loop on groups
            
            % get out of loop if you find less than k groups
            if ltkg==1
                Beta=NaN;
                break
            end
            
            %%%%%%%%%%%%%%%%%
            if cwm==1
                for j=1:k
                    % Eigenvalue eigenvector decomposition for group j
                    [Uj,Lambdaj] = eig(sigmaX(:,:,j));
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_pk(:,j)=diag(Lambdaj);
                end
                
            end
            
            % if equalweights =1 then equal proportions are supplied for group
            % sizes
            if equalweights==1
                sigma2ini= restreigen(sigma2ini,ones(k,1),restrfact,tolrestreigen,userepmat);
                
                if cwm==1
                    autovalues= restreigen(Lambda_pk,ones(k,1),restrfactX,tolrestreigen,userepmat);
                end
            else
                sigma2ini= restreigen(sigma2ini,niini,restrfact,tolrestreigen,userepmat);
                if cwm==1
                    autovalues= restreigen(Lambda_pk,niini,restrfactX,tolrestreigen,userepmat);
                end
                
            end
            
            
            
            if cwm==1
                % Covariance matrices are reconstructed keeping into account the
                % constraints on the eigenvalues
                for j=1:k
                    sigmaX(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
                    
                    % Alternative code: in principle more efficient but slower
                    % because diag is a built in function
                    % sigmaX(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
                end
            end
            
            %
            % Stop if two consecutive concentration steps have the same result
            if csteps_stop ==1
                if idx == indold
                    break
                else
                    indold = idx;
                end
            end
            %% 4) Computation of the value of the target function
            obj = 0;
            not_empty_g = seqk(~( niini <= p + 1 ));
            if mixt == 0
                
                for jj = not_empty_g
                    
                    ptermAIC = 2*p/mean_w4trim;
                    if equalweights ==1
                        %obj = obj + log(1/k) +...
                        %    sum(logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj)).*Z(:,jj));
                        %se nn si vogliono usare le thinnate nella obj
                        %function cancellare l'if e usare solo lo
                        %statement che c'e' nell'else
                        if wtrim==3 && assign_thinned_units ==1
                            obj = obj + log(1/k) +...
                                sum(logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj)).*Z_all(:,jj)) ;
                            
                        else
                            
                            obj = obj + log(1/k) +...
                                sum(logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj)).*Z(:,jj)) ;
                        end
                        
                    else
                        %obj = obj +    niini(jj)*log(niini(jj)/sum(niini))+...
                        %    sum(logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj)).*Z(:,jj));
                        obj = obj + niini(jj)*log(niini(jj)/sum(niini)) +...
                            sum(logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj)).*Z(:,jj));
                    end
                    
                    if cwm==1
                        obj=obj+sum(logmvnpdfFS(X(:,(intercept+1):end),muX(:,jj),sigmaX(:,:,jj)).*Z(:,jj));
                        %obj con AIC????????
                    end
                    
                end
                
                obj = 2*obj - ptermAIC;
                %TODO: check CWM and mixture likelihood
            else
                % compute mixture likelihood
                
                % Select all units not trimmed and not thinned
                % Target function is based on these units (not
                % trimmed and not thinned)
                
                log_lh=NaN(n,size(not_empty_g,2));
                
                for jj = 1:k
                    if equalweights ==1
                        log_lh(:,jj) = ...
                            log(1/k) + (logmvnpdfFS(y -X * Beta(:,jj),0,sigma2ini(jj) ) );
                    else
                        log_lh(:,jj) = ...
                            log(niini(jj)/sum(niini)) + (logmvnpdfFS(...
                            y - X*Beta(:,jj),0,sigma2ini(jj) ) );
                    end
                    if cwm==1
                        log_lh(:,jj)=log_lh(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(:,jj),sigmaX(:,:,jj));
                    end
                end
                log_lh(idx<=0,:)=[];
                obj = estepFS(log_lh);
                %obj con AIC????????
            end
            
            
        end
        % disp([cstep obj sum(indmax<=0)])
        obj_all(i,cstep) = obj;
        %        indmax_all(:,cstep,i) = indmax;
        %        beta_all(:,cstep,i) = Beta;
        
        
        %% Change the 'optimal' target value and 'optimal' parameters
        
        % The following if condition is done to check if an increase in the target value is achieved
        if zigzag == 1
            if obj >= vopt && sum(sum(isnan(Beta))) ==0
                if msg
                    disp(['sample = ' num2str(i) ', concentration step = ' num2str(cstep)])
                end
                cstepopt=cstep;
                vopt                    = obj;
                bopt                    = Beta;
                nopt                  = niini;
                sigma2opt                = sigma2ini;
                %w4trimopt = (nx1) vector of {0,1} for thinned and unthinned units respectively.
                w4trimopt               = w4trim;
                
                %idxopt = (nx1) vector of {-1,-2, 1, ..., k}
                if wtrim==3
                    idxopt                   = idx_ne0;
                else
                    idxopt                   = idx;
                end
                %idxopt_before_tropt     = indmax_before_tr; % TO DELETE
                if mixt ==2
                    postprobopt = postprob;
                end
                
                if cwm==1
                    muXopt=muX;
                    sigmaXopt=sigmaX;
                end
            end
        end
        
    end
    %%% Concentration steps ending
    
    
    %% Change the 'optimal' target value and 'optimal' parameters
    
    % The following if condition is done to check if an increase in the target value is achieved
    if zigzag == 0
        
        if obj >= vopt && sum(sum(isnan(Beta))) ==0
            cstepopt=cstep;
            vopt                    = obj;
            bopt                    = Beta;
            nopt                  = niini;
            sigma2opt                = sigma2ini;
            w4trimopt               = w4trim;
            %idxopt = (nx1) vector of {-1,-2, 1, ..., k}
            if wtrim==3
                idxopt                   = idx_ne0;
            else
                idxopt                   = idx;
            end
            
            if mixt ==2
                postprobopt = postprob;
            end
            
            if cwm==1
                muXopt=muX;
                sigmaXopt=sigmaX;
            end
        end
    end
    
    % monitor time execution
    if msg==1
        if i <= tsampling
            % sampling time until step tsampling
            time(i)=toc(tstart);
        elseif i==tsampling+1
            % stop sampling and print the estimated time
            fprintf('Total estimated time to complete tclustreg: %5.2f seconds \n', nselected*median(time));
        end
    end
end % end of loop over the nsamp subsets
%%%  Random starts: ending


%%  Compute quantities to be stored in the output structure or used in the plots

% Apply consistency factor based on the variance of the truncated normal distribution.

% number of non trimmed observations, after first and second level trimming
hh = sum(nopt);

% Compute variance of the truncated normal distribution
% Note that 1-hh/n is the trimming percentage
vt = norminv(0.5*(1+hh/n));

%factor=1/sqrt(1-(2*vt.*normpdf(vt))./(2*normcdf(vt)-1));
if hh<n
    factor = 1/sqrt(1-2*(n/hh)*vt.*normpdf(vt));
else
    factor=1;
end

% Apply the asymptotic consistency factor to the preliminary squared scale estimate
sigma2opt_corr=sigma2opt*factor;
% Apply small sample correction factor of Pison et al.
sigma2opt_corr=sigma2opt_corr*corfactorRAW(1,n,hh/n);

%%  Set the output structure

out                     = struct;
%vopt              = objective function
out.vopt                =vopt;
%cstepopt =        = cstep with maximum obj
out.cstepopt            =cstepopt;
%   bopt           = regression parameters
out.bopt                = bopt;
%   sigmaopt0      = estimated group variances
out.sigma2opt          = sigma2opt;
%   sigma2opt_corr  = estimated group variances corrected with  asymptotic
%   consistency factor and small sample correction factor
out.sigma2opt_corr       = sigma2opt_corr;

if cwm==1
    %            out.muXopt= k-by-p matrix containing cluster centroid
    %                       locations. Robust estimate of final centroids of
    %                       the groups.
    out.muXopt=muXopt;
    %         out.sigmaXopt= p-by-p-by-k array containing estimated constrained
    %                       covariance covariance matrices of the explanatory variables for the k groups.
    out.sigmaXopt=sigmaXopt;
end

%   obj           = value of the target function
out.obj                = vopt;
out.obj_all            = obj_all;
%out.indmax_all         = indmax_all;
%out.beta_all           = beta_all;
% out.idx = final allocation vector
% out.idx==1 for units allocated to group 1
% out.idx==2 for units allocated to group 2
%.....
% out.idx==k for units allocated to group k
% out.idx==-1 for first  level trimmed units
% out.idx==-2 for second level trimmed units

%in tandem thinning it is necessary to save information about not only retained units (idxopt), but also about thinned units.
if wtrim == 4
    out.idx = zeros(length(Xori),1);
    out.idx(id_unthinned)=idxopt;
else
    out.idx   = idxopt;
end
% frequency distribution of the allocations
out.siz=tabulateFS(idxopt(:,1));

%out.idx_before_tr   = idxopt_before_tropt; % TO DELETE

out.we                 =we;

if mixt == 2
    out.postprobopt     = postprobopt;
end
out.restrfact           = restrfact;
% Store the indices in varargout
if nargout==2
    varargout={C};
end

%% Compute INFORMATION CRITERIA

% Discriminant functions for the assignments
if equalweights == 1
    for jj = 1:k
        ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*bopt(:,jj),0,sigma2opt(jj));
        if cwm==1
            ll(:,jj)=  ll(:,jj)+ logmvnpdfFS(X(:,(intercept+1):end),muXopt(:,jj),sigmaXopt(:,:,jj));
        end
    end
else
    for jj = 1:k
        ll(:,jj) = log((nopt(jj)/sum(nopt))) + logmvnpdfFS(y-X*bopt(:,jj),0,sigma2opt(jj));
        if cwm==1
            ll(:,jj)=  ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muXopt(:,jj),sigmaXopt(:,:,jj));
        end
    end
end

% Now remove the rows which refer to first, or second level trimmed units
% or thinned units

delunits=false(n,1);
delunits(idxopt(:,end)<0)=true;
delunits(w4trim==0)=true;

ll(delunits,:)=[];

if mixt>=1
    [NlogLmixt]=estepFS(ll);
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


npar=p*k;

if equalweights==false
    npar=npar +(k-1);
end

% Add paramters referred to sigma2 restrictions
nParam=npar+ (k-1)*((1-1/restrfact)^(1-1/k)) +1;

if cwm==1
    nParam=nParam+ 0.5*p*(p-1)*k + (p*k-1)*((1-1/restrfactX)^(1-1/(p*k))) +1;
end

logh=log(h);

if mixt>0
    % MIXMIX = BIC which uses parameters estimated using the mixture loglikelihood
    % and the maximized mixture likelihood as goodness of fit measure (New BIC)
    MIXMIX  = 2*NlogLmixt +nParam*logh;
    
    % MIXCLA = BIC which uses the classification likelihood based on
    % parameters estimated using the mixture likelihood (New ICL)
    MIXCLA  = 2*NlogL +nParam*logh;
    
    out.MIXMIX=MIXMIX;
    out.MIXCLA=MIXCLA;
else
    % CLACLA = BIC which uses parameters estimated using the classification
    % likelihood and the maximized classification likelihood as goodness of fit
    % measure (New New)
    CLACLA  = 2*NlogL +nParam*logh;
    
    out.CLACLA=CLACLA;
end


%% Generate plots

if plots
    
    % this is just for rotating colors in the plots
    clrdef = 'bkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcr';
    symdef = '+*d^v><phos+*d^v><phos+*d^v><phos+*d^v><phos';
    
    % The following plots are for the bi-variate case (i.e. v=1)
    if p-intercept < 2
        
        % initialize figure
        fh = figure('Name','TclustReg plot','NumberTitle','off','Visible','on');
        gca(fh);
        hold on;
        
        title({['$ wtrim=' num2str(wtrim) '\quad mixt=' num2str(mixt) , '  \quad c=' num2str(restrfact) '\quad \alpha_1=' num2str(alphaLik) '\quad \alpha_2=' num2str(alphaX) '$'] , ...
            ['$ obj=' num2str(out.obj) '\quad b=(' sprintf('%0.3f ;',out.bopt(end,:)) ') $']} , ...
            'interpreter' , 'LaTex', 'fontsize' , 14);
        
        for jj = 1:k
            group_label = ['Group ' num2str(jj)];
            
            % plot of the good units allocated to the current group.
            % Indices are taken after the second level trimming.
            % Trimmed points are not plotted by group.
            if wtrim==4
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
            
            
            
        end
        
        % Plot the outliers (trimmed points)
        if wtrim==4
            ucg = find(idxopt==-1);
        else
            ucg = find(out.idx(:,1)==-1);
        end
        plot(X(ucg,end),y(ucg),'o','color','r','MarkerSize',8,...
            'DisplayName',['Trimmed units 1st (' num2str(length(ucg)) ')']);
        
        % Second level trimming points
        ucg = find(out.idx(:,1)==-2);
        plot(X(ucg,end),y(ucg),'*','color','c',...
            'DisplayName',['Trimmed units 2nd (' num2str(length(ucg)) ')']);
        
        if wtrim == 4
            % in case of tandem thinning, plot the thinned points
            plot(Xori(~Wt4,end),yori(~Wt4),symdef(k+1),'color',clrdef(k+1),...
                'DisplayName',['Thinned units (' num2str(length(Wt4) - sum(Wt4)) ')']);
        end
        
        % Position the legends and make them clickable. For some reason
        % clickableMultiLegend does not set properly the FontSize: to be fixed.
        lh=legend('show');
        legstr = get(lh,'String');
        clickableMultiLegend(legstr,'Location','northwest','interpreter' , 'LaTex', 'fontsize' , 12);
        
        axis('manual');
        
        % control of the axis limits
        xmin = min(X(:,end)); xmax = max(X(:,end));
        ymin = min(y); ymax = max(y);
        deltax = (xmax - xmin) / 10;
        deltay = (ymax - ymin) / 10;
        
        xlim([xmin-deltax,xmax+deltax]);
        ylim([ymin-deltay,ymax+deltay]);
        
        %%
        % initialize figure
        
        if wtrim == 3
            fh = figure('Name','TclustReg plot','NumberTitle','off','Visible','on');
            gca(fh);
            hold on;
            
            title({['$ wtrim=' num2str(wtrim) '\quad mixt=' num2str(mixt) , '  \quad c=' num2str(restrfact) '\quad \alpha_1=' num2str(alphaLik) '\quad \alpha_2=' num2str(alphaX) '$'] , ...
                ['$ obj=' num2str(out.obj) '\quad b=(' sprintf('%0.3f ;',out.bopt(end,:)) ') $']} , ...
                'interpreter' , 'LaTex', 'fontsize' , 14);
            
            for jj = 1:k
                group_label = ['Group ' num2str(jj)];
                
                % plot of the good units allocated to the current group.
                % Indices are taken after the second level trimming.
                % Trimmed points are not plotted by group.
                ucg = find(out.idx(:,end)==jj);
                
                
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
                
                %plot the thinned (not trimmed) units
                if jj == k
                    thinned_nt_trimmed = w4trimopt;
                    thinned_nt_trimmed([ones(n,1) ;ones(n,1)]) = -12;
                    ucg = find(thinned_nt_trimmed == 0);
                    plot(X(ucg,end),y(ucg),symdef(jj),'color',clrdef(k+1),...
                        'DisplayName',['thinned units (' num2str(length(ucg)) ')' ]);
                end
                
            end
            
            % Plot the outliers (trimmed points)
            ucg = find(out.idx(:,end)==-1);
            plot(X(ucg,end),y(ucg),'o','color','r','MarkerSize',8,...
                'DisplayName',['Trimmed units 1st (' num2str(length(ucg)) ')']);
            
            % Second level trimming points
            ucg = find(out.idx(:,end)==-2);
            plot(X(ucg,end),y(ucg),'*','color','c',...
                'DisplayName',['Trimmed units 2nd (' num2str(length(ucg)) ')']);
            
            if wtrim == 4
                % in case of tandem thinning, plot the thinned points
                plot(Xori(~Wt4,end),yori(~Wt4),symdef(k+1),'color',clrdef(k+1),...
                    'DisplayName',['Thinned units (' num2str(length(Wt4) - sum(Wt4)) ')']);
            end
            
            % Position the legends and make them clickable. For some reason
            % clickableMultiLegend does not set properly the FontSize: to be fixed.
            lh=legend('show');
            legstr = get(lh,'String');
            clickableMultiLegend(legstr,'Location','northwest','interpreter' , 'LaTex', 'fontsize' , 12);
            
            axis('manual');
            
            % control of the axis limits
            xmin = min(X(:,end)); xmax = max(X(:,end));
            ymin = min(y); ymax = max(y);
            deltax = (xmax - xmin) / 10;
            deltay = (ymax - ymin) / 10;
            
            xlim([xmin-deltax,xmax+deltax]);
            ylim([ymin-deltay,ymax+deltay]);
        end
    else
        
        % In this case p > 2. A standard spmplot is used.
        
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
        %%
        
        
        
    end
    
end


%% Subfunctions

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
