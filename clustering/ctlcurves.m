function out  = ctlcurves(Y, varargin)
%ctlcurves computes Classification Trimmed Likelihood Curves
%
%<a href="matlab: docsearchFS('ctlcurves')">Link to the help function</a>
%
%   ctlcurves applies tclust several times on a given dataset while
%   parameters alpha and k are altered. The resulting object gives an idea
%   of the optimal trimming level and number of clusters considering a
%   particular dataset.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. Data matrix containing n observations
%               on v variables, Rows of Y represent observations, and
%               columns represent variables. Observations (rows) with
%               missing (NaN) or or infinite (Inf) values will
%               automatically be excluded from the computations.
%                 Data Types -  double
%
%  Optional input arguments:
%
%        alpha: trimming level to monitor. Vector. Vector which specifies the
%               values of trimming levels which have to be considered.
%               For example if alpha=[0 0.05 0.1] ctlcurves considers these 3
%               values of trimming level.
%               The default for alpha is vector 0:0.02:0.10;
%                 Example - 'alpha',[0 0.05 0.1]
%                 Data Types -  double
%
%       bands  : confidence bands for the curves. boolean or struct. If
%               bands is a scalar boolean equal to true (default), 50 per
%               cent confidence bands are computed (and are shown on the
%               screen if plots=1), likelihood ratio tests are also
%               computed and the solutions found by these two methods are
%               given. If bands is a scalar boolean equal to false, bands
%               and likelihood ratio tests are not computed.
%             If bands is a struct bands are computed and the structure may
%             contain the following fields:
%           bands.conflev = scalar in the interval (0 1) which contains
%                   the confidence level of the bands (default is 0.5).
%           bands.nsamp =  Number of subsamples to extract in the
%                   bootstrap replicates. If this field is not present and
%                   or it is empty the number of subsamples which is used
%                   in the bootstrap replicates is equal to the one used
%                   for real data (input option nsamp).
%           bands.nsimul = number of replicates to use to create the
%                   confidence bands and likelihood ratio tests.
%                   The default value of bands.nsimul is
%                   60 in order to provide the output in a reasonal time.
%                   Note that for stable results we recommentd a value of
%                   bands.nsimul equal to 100.
%            bands.valSolution   = boolean which specifies if it is
%                   necessary to perform an outlier detection procedure on
%                   the components which have been found using optimalK and
%                   optimal alpha. If bands.valSolution is true then units
%                   detected as outliers in each component are assigned to
%                   the noise group. If bands.valSolution is false
%                   (default) or it is not present nothing is done.
%            bands.usepriorSol= Initialization of the EM for a particular subset.
%                   boolean. If bands.usepriorSol is true, we initialize
%                   the EM algorithm for one of the bands.nsamp subsets
%                   with the solution which has been found with real data,
%                   else if it is false (default) no prior information is
%                   used in any of the extracted subsets.
%   bands.outliersFromUniform = way of generating the outliers in the
%                   bootstrap replicates. Boolean. If outliersFromUniform
%                   is true (default) the outliers are generated using the
%                   uniform distribution in such a way that their squared
%                   Mahalanobis distance from the centroids of each
%                   existing group is larger then the quantile 1-0.999 of
%                   the Chi^2 distribution with p degrees of freedom. For
%                   additional details see input option noiseunits of
%                   simdataset. If outliersFromUniform is false the
%                   outliers are the units which have been trimmed after
%                   applying tclust for the particular combination of
%                   values of k and alpha.
%            bands.nsimulExtra =  number of replicates to use to compute the
%                   empirical pvalue of LRtest if multiple solutions are
%                   found given a value of alpha (default is 50).
%            bands.nsampExtra =  number of subsamples to extract in each
%                   replicate to compute the empirical pvalue of LRtest if
%                   multiple solutions are found given a value of alpha
%                   (default is 2000).
%           bands.usepriorSolExtra= Initialization of the EM for particular subset.
%                   boolean. If bands.usepriorSolExtra is true, we
%                   initialize the EM algorithm for one of the
%                   bands.nsampExtra subsets with the solution which has
%                   been found with real data, else if it is false
%                   (default) no prior information is used in any of the
%                   extracted subsets.
%                 Example - 'bands',true
%                 Data Types - logical | struct
%
%           kk: number of mixture components. Integer vector. Integer
%               vector specifying the number of mixture components
%               (clusters) for which trimmed likelihoods are calculated.
%               Vector. The default value of kk is 1:5.
%                 Example - 'kk',1:4
%                 Data Types - int16 | int32 | single | double
%
%       cshape : constraint to apply to each of the shape matrices.
%                Scalar greater or equal than 1. This options only works
%                is 'restrtype' is 'deter'.
%               When restrtype is deter the default value of the "shape"
%               constraint (as defined below) applied to each group is
%               fixed to $c_{shape}=10^{10}$, to ensure the procedure is
%               (virtually) affine equivariant. In other words, the
%               decomposition or the $j$-th scatter matrix $\Sigma_j$ is
%               \[
%               \Sigma_j=\lambda_j^{1/v} \Omega_j \Gamma_j \Omega_j'
%               \]
%               where $\Omega_j$ is an orthogonal matrix of eigenvectors, $\Gamma_j$ is a
%               diagonal matrix with $|\Gamma_j|=1$ and with elements
%               $\{\gamma_{j1},...,\gamma_{jv}\}$ in its diagonal (proportional to
%               the eigenvalues of the $\Sigma_j$ matrix) and
%               $|\Sigma_j|=\lambda_j$. The $\Gamma_j$ matrices are commonly
%               known as "shape" matrices, because they determine the shape of the
%               fitted cluster components. The following $k$
%               constraints are then imposed on the shape matrices:
%               \[
%               \frac{\max_{l=1,...,v} \gamma_{jl}}{\min_{l=1,...,v} \gamma_{jl}}\leq
%                   c_{shape}, \text{ for } j=1,...,k,
%               \]
%               In particular, if we are ideally searching for spherical
%               clusters it is necessary to set  $c_{shape}=1$. Models with
%               variable volume and spherical clusters are handled with
%               'restrtype' 'deter', $1<restrfactor<\infty$ and $cshape=1$.
%               The $restrfactor=cshape=1$ case yields a very constrained
%               parametrization because it implies spherical clusters with
%               equal volumes.
%                 Example - 'cshape',10
%                 Data Types - single | double
%
%equalweights : cluster weights in the concentration and assignment steps.
%               Logical. A logical value specifying whether cluster weights
%               shall be considered in the concentration, assignment steps
%               and computation of the likelihood.
%                 Example - 'equalweights',true
%                 Data Types - Logical
%
%       msg  :  Message on the screen. Scalar. Scalar which controls whether
%               to display or not messages about code execution.
%                 Example - 'msg',1
%                 Data Types - single | double
%
%       mixt  : Mixture modelling or crisp assignment. Scalar.
%               Option mixt specifies whether mixture modelling or crisp
%               assignment approach to model based clustering must be used.
%               In the case of mixture modelling parameter mixt also
%               controls which is the criterion to find the untrimmed units
%               in each step of the maximization.
%               If mixt >=1 mixture modelling is assumed else crisp
%               assignment. The default value is mixt=0 (i.e. crisp assignment).
%               In mixture modelling the likelihood is given by:
%                \[
%                \prod_{i=1}^n  \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j),
%                \]
%               while in crisp assignment the likelihood is given by:
%               \[
%               \prod_{j=1}^k   \prod _{i\in R_j} \phi (y_i; \; \theta_j),
%               \]
%               where $R_j$ contains the indexes of the observations which
%               are assigned to group $j$.
%               Remark - if mixt>=1 previous parameter equalweights is
%               automatically set to 1.
%               Parameter mixt also controls the criterion to select the
%               units to trim,
%               if mixt = 2 the h units are those which give the largest
%               contribution to the likelihood that is the h largest
%               values of:
%               \[
%                   \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j)   \qquad
%                    i=1, 2, ..., n,
%               \]
%               else if mixt=1 the criterion to select the h units is
%               exactly the same as the one which is used in crisp
%               assignment. That is: the n units are allocated to a
%               cluster according to criterion:
%               \[
%                \max_{j=1, \ldots, k} \hat \pi'_j \phi (y_i; \; \hat \theta_j)
%               \]
%               and then these n numbers are ordered and the units
%               associated with the largest h numbers are untrimmed.
%                   Example - 'mixt',1
%                   Data Types - single | double
%
%      nocheck: Check input arguments. Scalar. If nocheck is equal to 1
%               no check is performed on matrix Y. As default nocheck=0.
%                 Example - 'nocheck',10
%                 Data Types - single | double
%
%       nsamp : number of subsamples to extract. Scalar or matrix.
%               If nsamp is a scalar it contains the number of subsamples
%               which will be extracted. If nsamp=0 all subsets will be
%               extracted.
%               Remark - if the number of all possible subset is <300 the
%               default is to extract all subsets, otherwise just 300.
%               - If nsamp is a matrix it contains in the rows the indexes
%                 of the subsets which have to be extracted. nsamp in this
%                 case can be conveniently generated  by function subsets.
%                 nsamp can have k columns or k*(v+1) columns. If nsamp has
%                 k columns the k initial centroids each iteration i are
%                 given by X(nsamp(i,:),:) and the covariance matrices are
%                 equal to the identity.
%               - If nsamp has k*(v+1) columns the initial centroids and
%                 covariance matrices in iteration i are computed as follows
%                 X1=X(nsamp(i,:),:)
%                 mean(X1(1:v+1,:)) contains the initial centroid for group 1
%                 cov(X1(1:v+1,:)) contains the initial cov matrix for group 1               1
%                 mean(X1(v+2:2*v+2,:)) contains the initial centroid for group 2
%                 cov((v+2:2*v+2,:)) contains the initial cov matrix for group 2               1
%                 ...
%                 mean(X1((k-1)*v+1:k*(v+1))) contains the initial centroids for group k
%                 cov(X1((k-1)*v+1:k*(v+1))) contains the initial cov matrix for group k
%               REMARK - if nsamp is not a scalar option option below
%               startv1 is ignored. More precisely if nsamp has k columns
%               startv1=0 elseif nsamp has k*(v+1) columns option startv1=1.
%                 Example - 'nsamp',1000
%                 Data Types - double
%
%    refsteps : Number of refining iterations. Scalar. Number of refining
%               iterations in subsample.  Default is 15. refsteps = 0 means
%               "raw-subsampling" without iterations.
%                 Example - 'refsteps',10
%                 Data Types - single | double
%
%     reftol  : scalar. Default value of tolerance for the refining steps.
%               The default value is 1e-14;
%                 Example - 'reftol',1e-05
%                 Data Types - single | double
%
%  restrfactor: Restriction factor. Scalar. Positive scalar which
%               constrains the allowed differences
%               among group scatters. Larger values imply larger differences of
%               group scatters. On the other hand a value of 1 specifies the
%               strongest restriction forcing all eigenvalues/determinants
%               to be equal and so the method looks for similarly scattered
%               (respectively spherical) clusters. The default is to apply
%               restrfactor to eigenvalues. In order to apply restrfactor
%               to determinants it is is necessary to use optional input
%               argument cshape.
%                 Example - 'restrfactor',12
%                 Data Types -  double
%
%     restrtype : type of restriction. Character. The type of restriction to
%               be applied on the cluster scatter
%               matrices. Valid values are 'eigen' (default), or 'deter'.
%               eigen implies restriction on the eigenvalues while deter
%               implies restriction on the determinants. If restrtype is
%               'deter' it is possible to control the constraints on the
%               shape matrices using optional input argument cshape.
%                 Example - 'restrtype','deter'
%                 Data Types - char
%
%      startv1: how to initialize centroids and cov matrices. Scalar.
%               If startv1 is 1 then initial
%               centroids and and covariance matrices are based on (v+1)
%               observations randomly chosen, else each centroid is
%               initialized taking a random row of input data matrix and
%               covariance matrices are initialized with identity matrices.
%               Remark 1- in order to start with a routine which is in the
%               required parameter space, eigenvalue restrictions are
%               immediately applied. The default value of startv1 is 1.
%               Remark 2 - option startv1 is used just if nsamp is a scalar
%               (see for more details the help associated with nsamp).
%                 Example - 'startv1',1
%                 Data Types - single | double
%
%
%  cleanpool :  clean pool. Scalar. cleanpool is 1 if the parallel pool has
%               to be cleaned after the execution of the routine. Otherwise
%               it is 0. The default value of cleanpool is 0. Clearly this
%               option has an effect just if previous option numpool is >
%               1.
%                 Example - 'cleanpool',1
%                 Data Types - single | double
%
%
%     numpool : number of pools for parellel computing. Scalar.
%               If numpool > 1, the routine automatically checks if
%               the Parallel Computing Toolbox is installed and distributes
%               the random starts over numpool parallel processes. If
%               numpool <= 1, the random starts are run sequentially. By
%               default, numpool is set equal to the number of physical
%               cores available in the CPU (this choice may be inconvenient
%               if other applications are running concurrently). The same
%               happens if the numpool value chosen by the user exceeds the
%               available number of cores.
%               REMARK 1: Unless you adjust the cluster profile, the
%               default maximum number of workers is the same as the
%               number of computational (physical) cores on the machine.
%               REMARK 2: In modern computers the number of logical cores
%               is larger than the number of physical cores. By default,
%               MATLAB is not using all logical cores because, normally,
%               hyper-threading is enabled and some cores are reserved to
%               this feature.
%               REMARK 3: It is because of Remarks 3 that we have chosen as
%               default value for numpool the number of physical cores
%               rather than the number of logical ones. The user can
%               increase the number of parallel pool workers allocated to
%               the multiple start monitoring by:
%               - setting the NumWorkers option in the local cluster profile
%                 settings to the number of logical cores (Remark 2). To do
%                 so go on the menu "Home|Parallel|Manage Cluster Profile"
%                 and set the desired "Number of workers to start on your
%                 local machine".
%               - setting numpool to the desired number of workers;
%               Therefore, *if a parallel pool is not already open*,
%               UserOption numpool (if set) overwrites the number of
%               workers set in the local/current profile. Similarly, the
%               number of workers in the local/current profile overwrites
%               default value of 'numpool' obtained as feature('numCores')
%               (i.e. the number of physical cores).
%                 Example - 'numpool',4
%                 Data Types - double
%
%       plots : Plot on the screen. Scalar. If plots = 1, a plot of the
%               CTLcurves is shown on the screen. If input option bands is
%               not empty confidence bands are also shown.
%                 Example - 'plots',1
%                 Data Types - single | double
%
%       Ysave : save input matrix. Boolean.
%               Boolan that is set to true to request that the input matrix Y
%               is saved into the output structure out. Default is 1, that
%               is  matrix Y is saved inside output structure out.
%                 Example - 'Ysave',false
%                 Data Types - logical
%
%
%
%  Output:
%
%         out:   structure which contains the following fields:
%
%                out.Mu = cell of size length(kk)-by-length(alpha)
%                       containing the estimate of the centroids for each
%                       value of k and each value of alpha. More precisely,
%                       suppose kk=1:4 and alpha=[0 0.05 0.1], out.Mu{2,3}
%                       is a matrix with two rows and v columns containing
%                       the estimates of the centroids obtained when
%                       alpha=0.1.
%            out.Sigma = cell of size length(kk)-by-length(alpha)
%                       containing the estimate of the covariance matrices
%                       for each value of k and each value of alpha. More
%                       precisely, suppose kk=1:4 and alpha=[0 0.05 0.1],
%                       out.Sigma{2,3} is a 3D  array of size v-by-v-by-2
%                       containing the estimates of the covariance matrices
%                       obtained when alpha=0.1.
%            out.Pi   = cell of size length(kk)-by-length(alpha)
%                       containing the estimate of the group proportions
%                       for each value of k and each value of alpha. More
%                       precisely, suppose kk=1:4 and alpha=[0 0.05 0.1],
%                       out.Pi{2,3} is a 3D  array of size v-by-v-by-2
%                       containing the estimates of the covariance matrices
%                       obtained when alpha=0.1.
%            out.IDX   = cell of size length(kk)-by-length(alpha)
%                       containing the final assignment for each value of k
%                       and each value of alpha. More precisely, suppose
%                       kk=1:4 and alpha=[0 0.05 0.1], out.IDX{2,3} is a
%                       vector of length(n) containing the containinig the
%                       assignment of each unit obtained when alpha=0.1.
%                       Elements equal to zero denote unassigned units.
%           out.CTL    = matrix of size length(kk)-by-length(alpha)
%                       containing the values of the trimmed likelihood
%                       curves for each value of k and each value of alpha.
%                       This output (and all the other output below which
%                       start with CTL) are present only if input option
%                       bands is true or is a struct. All the other fields
%                       belo
%      out.CTLbands    = 3D array of size
%                       length(kk)-by-length(alpha)-by-nsimul containing
%                       the nsimul replicates of out.CTL.
%         out.CTLlikLB    =  matrix of size length(kk)-by-length(alpha)
%                       containing the lower confidence bands of the
%                       trimmed likelihood curves for each value of k and
%                       each value of alpha.
%         out.CTLlikUB    =  matrix of size length(kk)-by-length(alpha)
%                       containing the upper confidence bands of the
%                       trimmed likelihood curves for each value of k and
%                       each value of alpha. T
%         out.CTLlik050    =  matrix of size length(kk)-by-length(alpha)
%                       containing the central confidence bands of the
%                       trimmed likelihood curves for each value of k and
%                       each value of alpha.
%      out.CTLtentSol  = matrix with size m-by-4. Details of the ordered
%                          solutions where there was intersection between
%                          two consecutive trimmed likelihood curves. First
%                          column contains the value of k, second column
%                          the value of alpha, third column the index
%                          associated to the best value of alpha, fourth
%                          colum index associated with the best value of
%                          kk.
%       out.CTLoptimalAlpha = scalar, optimal value of trimming.
%       out.CTLoptimalK = scalar, optimal number of clusters, stored
%                        as a positive integer value.
%       out.CTLoptimalIDX  = n-by-1 vector containing assignment of each unit to
%                       each of the k groups in correspodence of
%                       OptimalAlpha and OptimalK. Cluster names are
%                       integer numbers from 1 to k. 0 indicates trimmed
%                       observations. The fields which follow which start
%                       with LRT refer to the likilhood ratio test
%
%        out.LRTpval =  table with size length(kk)-1-times-length(alpha)
%                           which stores the relative frequency in which
%                           the Likelihood ratio test is greater than the
%                           corresponding bootstrap test.
%                        as a positive integer value.
%        out.LRTtentSol  = matrix with size m-by-8. Details of the ordered
%                          solutions using the likelihood ratio test. First
%                          column (index): the index number of the
%                          solution. Second column (k): the value of k.
%                          Third column (alpha): the value of alpha. Fourth
%                          column (Truesol) contains 1 if the p-value
%                          beyond the threshold is always above $k$ shown
%                          in the second column for all $k^* >k$ else it
%                          contains 0. Fifth and sixth columns (kindex and
%                          alphaindex) contain the index numbers of optimal
%                          input values kk and alpha. Seventh column
%                          (kbestGivenalpha) contains 1 in correspondence
%                          of the best solution for each value of k (given
%                          alpha). Eight column (ProperSize) contains 1 if
%                          the solution which has been found has a minimum
%                          group size which is greater than n*max(alpha).
%      out.LRTtentSolt  = table with size m-by-5 containing the same
%                        information of array  out.TentSolLR in table format.
%      out.LRTtentSolIDX = matrix with size n-by-size(out.LRTtentSol,1) with
%                        the allocation associated with the tentative
%                        solutions found in out.LRTtentSol.
%                        First column refers to solution in row 1 of out.LRTtentSol ...
%       out.LRToptimalAlpha = scalar, optimal value of trimming.
%       out.LRToptimalK = scalar, optimal number of clusters, stored
%                        as a positive integer value.
%       out.LRToptimalIDX  = n-by-1 vector containing assignment of each unit to
%                       each of the k groups in correspodence of
%                       OptimalAlpha and OptimalK. Cluster names are
%                       integer numbers from 1 to k. 0 indicates trimmed
%                       observations. The fields which follow which start
%                       with LRT refer to the likelihood ratio test. The
%                       optimal solution is the first for which
%                       out.LRTtentSolt.kbestGivenalpha is 1 and
%                       out.LRTtentSolt.ProperSize is 1.


%                out.kk = vector containing the values of k (number of
%                       components) which have been considered. This  vector
%                       is equal to input optional argument kk if kk had been
%                       specified else it is equal to 1:5.
%                out.alpha = vector containing the values of the trimming
%                       level which have been considered. This
%                       vector is equal to input optional argument alpha.
%         out.restrfactor = scalar containing the restriction factor
%                       which has been used to compute tclust.
%                out.Y  = Original data matrix Y. The field is present if
%                       option Ysave is set to 1.
%
% More About:
%
% These curves show the values of the trimmed classification
% (log)-likelihoods when altering the trimming proportion alpha and the
% number of clusters k. The careful examination of these curves provides
% valuable information for choosing these parameters in a clustering
% problem. For instance, an appropriate k to be chosen is one that we do
% not observe a clear increase in the trimmed classification likelihood
% curve for k with respect to the k+1 curve for almost all the range of
% alpha values. Moreover, an appropriate choice of parameter alpha may be
% derived by determining where an initial fast increase of the trimmed
% classification likelihood curve stops for the final chosen k (Garcia
% Escudero et al. 2011). This routine adds confidence bands in order to
% provide an automatic criterion in order to choose k and alpha. Optimal
% trimming level is chosen as the smallest value of alpha such that the
% bands for two consecutive values of k intersect and computes a bootstrap
% test for two consecutive values of k given alpha.
%
%
% See also ctlcurves
%
%
% References:
%
% Garcia-Escudero, L.A.; Gordaliza, A.; Matran, C. and Mayo-Iscar A., (2011),
% "Exploring the number of groups in robust model-based
% clustering." Statistics and Computing, Vol. 21, pp. 585-599.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('ctlcurves')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of use of ctlcurves with all default options.
    % Load the geyser data
    rng('default')
    rng(10)
    Y=load('geyser2.txt');
    % Just use a very small sumber of subsets for speed reasons.
    nsamp=10;
    out=ctlcurves(Y,'nsamp',nsamp);
    % Show the automatic classification
    spmplot(Y,out.LRToptimalIDX);
%}

%{
    % Example of the use of options alpha, kk and bands.
    % Personalized levels of trimming.
    alpha=0:0.05:0.15;
    % bands is passed as a false boolean. No bands are shown.
    bands=false;
    Y=load('geyser2.txt');
    % Just use a very small sumber of subsets for speed reasons.
    rng(100)
    nsamp=10;
    out=ctlcurves(Y,'alpha',alpha,'kk',2:4,'bands',bands,'nsamp',nsamp);
%}

%{
    %% option bands passed as struct.
    % load the M5 data.
    Y=load('M5data.txt');
    Y=Y(:,1:2);
    bands=struct;
    % just 20 simulations to construct the bands.
    bands.nsimul=20;
    % Exclude the units detected as outliers for each group.
    bands.valSolution=true;
    % One of the extracted subsets is based on solutions found with real data 
    bands.usepriorSol=true;
    % Just use a very small sumber of subsets for speed reasons.
    nsamp=20;
    rng(100)
    out=ctlcurves(Y,'bands',bands,'kk',2:4,'alpha',0:0.02:0.1,'nsamp',nsamp);
    % Show final classification using intersection of confidence bands.
    spmplot(Y,out.CTLoptimalIDX);
%}

%{
    %% An example with many groups.
    % Simulate data with 5 groups in 3 dimensions
    rng(1)
    k=5; v=3;
    n=200;
    Y = MixSim(k, v, 'MaxOmega',0.01);
    [Y]=simdataset(n, Y.Pi, Y.Mu, Y.S, 'noiseunits', 10);
    % Just use a very small sumber of subsets for speed reasons.
    nsamp=5;
    out=ctlcurves(Y,'plots',0,'kk',4:6,'nsamp',nsamp);
    spmplot(Y,out.CTLoptimalIDX);
%}


%% Beginning of code

nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);
[n, v]=size(Y);

% check how many physical cores are available in the computer (warning:
% function 'feature' is undocumented; however, FSDA is automatically
% monitored for errors and other inconsistencies at each new MATLAB
% release).
numpool = feature('numCores');

% User options
% startv1def = default value of startv1 =1, initialization using covariance
% matrices based on v+1 units
startv1=1;

refsteps=15;
reftol=1e-5;

equalweights=false;
restr='eigen';

plots=1;
nsamp=500;
kk=1:5;
msg=1;
alphaTrim=0:0.02:0.10;

cleanpool=false;
% cshape. Constraint on the shape matrices inside each group which works only if restrtype is 'deter'
cshape=10^10;
restrfactor=100;
mixt=0;
bands=true;
outliersFromUniform=true;

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)

    options=struct('kk',kk,'mixt',mixt,'alpha',alphaTrim,...
        'nsamp',nsamp,'plots',plots,'nocheck',0,'bands',bands, ...
        'restrfactor', restrfactor, ...
        'msg',msg,'Ysave',1,'refsteps',refsteps,'equalweights',equalweights,...
        'reftol',reftol,'startv1',startv1,'restrtype',restr,'cshape',cshape,...
        'numpool',numpool, 'cleanpool', cleanpool);


    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:ctlcurves:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end

    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:ctlcurves:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end


    restr=options.restrtype;
    alphaTrim=options.alpha;
    kk=options.kk;
    nsamp=options.nsamp;        % Number of subsets to extract
    plots=options.plots;        % Plot the ctl curves
    equalweights=options.equalweights;    % Specify if assignment must take into account the size of the groups

    refsteps=options.refsteps;
    reftol=options.reftol;
    msg=options.msg;            % Scalar which controls the messages displayed on the screen

    mixt=options.mixt;
    cleanpool=options.cleanpool;
    numpool=options.numpool;
    cshape=options.cshape;
    restrfactor=options.restrfactor;
    bands=options.bands;

    % Make sure vectors kk and alphaTtrim are sorted.
    kk=sort(kk);
    alphaTrim=sort(alphaTrim);
end


lkk=length(kk);
lalpha=length(alphaTrim);

MuVal = cell(lkk,lalpha);
SigmaVal = MuVal;
PiVal = MuVal;
IDX=MuVal;
CTLVal=zeros(lkk,lalpha);

%% Preapare the pool (if required)

CnsampAll=cell(lkk,1);
gRandNumbForNiniAll=CnsampAll;

if isstruct(bands) || (islogical(bands) && bands ==true)
    ComputeBands = true;
else
    ComputeBands = false;
end

for k=1:lkk  % loop for different values of k (number of groups)

    seqk=kk(k);

    % Cnsamp=subsets(nsamp,n,(v+1)*seqk);
    %seqk = number of groups to consider
    if isscalar(nsamp)
        % For each value of seqk extract subsamples once and for all
        Cnsamp=subsets(nsamp,n,(v+1)*seqk);
    else
        Cnsamp=nsamp;
    end

    % For each value of k extract random numbers to initialize proportions
    % once and for all
    gRandNumbForNini=rand(seqk,nsamp);

    if ComputeBands==true
        CnsampAll{seqk}=Cnsamp;
        gRandNumbForNiniAll{seqk}=gRandNumbForNini;
    end

    parfor (j=1:lalpha , numpool)

        alphaTrimj=alphaTrim(j);

        outtc=tclust(Y,seqk,alphaTrimj,restrfactor,'nsamp',Cnsamp,'plots',0,'msg',0,'mixt',mixt, ...
            'restrtype',restr,'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
            'reftol',reftol,'RandNumbForNini',gRandNumbForNini,'cshape',cshape);

        % Store proportions
        [indice] = find(outtc.siz(:,1) >0);
        PiVal{k,j} = outtc.siz(indice,2)/sum(outtc.siz(indice,2));

        % columns (j) = values of alpha
        % rows (k) = number of groups

        % Store centroids
        MuVal{k,j} = outtc.muopt;
        % Store covariance matrices
        SigmaVal{k,j} = outtc.sigmaopt;

        % Store classification
        IDX{k,j} = outtc.idx;

        % Store objective function
        CTLVal(k,j) = outtc.obj;

    end
    if msg==1
        disp(['k=' num2str(seqk)])
    end
end

% gamma is the width of the confidence interval divided by 2
gamma=0.25;
% nsimul = default number of simulations to create the bands
nsimul=60;
% do not validate the components using outlier detection procedures.
valSolution=false;
% Initialization for number or subsets to extract when tclust is applied to
% simulated data.
nsampSimData=[];

out=struct;
out.Mu=MuVal;
out.Sigma=SigmaVal;
out.Pi=PiVal;
out.IDX=IDX;
out.CTL=CTLVal;

if ComputeBands==true
    % Suppress warning on All Workers
    parfevalOnAll(@warning,0,'off','all');
    if isstruct(bands)

        % make sure that the fields which are given are admissible
        bandsdef=struct;
        bandsdef.conflev=[];
        bandsdef.nsimul=[];
        bandsdef.valSolution=[];
        bandsdef.outliersFromUniform=[];
        bandsdef.nsampSimData=[];
        bandsdef.usepriorSol=[];
        bandsdef.usepriorSolExtra=[];
        bandsdef.nsimulExtra=[];
        bandsdef.nsampExtra=[];

        chkoptions(bandsdef,fieldnames(bands))

        if isfield(bands,'conflev')
            gamma=(1-bands.conflev)/2;
        end

        if isfield(bands,'nsimul')
            nsimul=bands.nsimul;
        end

        if isfield(bands,'valSolution')
            valSolution=bands.valSolution;
        end


        if isfield(bands,'outliersFromUniform')
            outliersFromUniform=bands.outliersFromUniform;
        end

        if isfield(bands,'nsamp')
            nsampSimData=bands.nsamp;
        else
            nsampSimData=[];
        end

        if isfield(bands,'usepriorSol')
            usepriorSol=bands.usepriorSol;
        else
            usepriorSol=false;
        end

        if isfield(bands,'usepriorSolExtra')
            usepriorSolExtra=bands.usepriorSolExtra;
        else
            usepriorSolExtra=false;
        end

        if isfield(bands,'nsimulExtra')
            nsimulExtra=bands.nsimulExtra;
        else
            nsimulExtra=50;
        end


        if isfield(bands,'nsampExtra')
            nsampExtra=bands.nsampExtra;
        else
            nsampExtra=2000;
        end
    else
        usepriorSol=false;
    end

    % CTLbands is a 3D array which will contain the replicates for the
    % solutions
    CTLbands=zeros(lkk,lalpha,nsimul);
    BandsCTLtest=zeros(lkk-1,lalpha,nsimul);
    % maxk = maximum allowed vector for k
    maxk=kk(end);

    for k=1:lkk  % loop for different values of k (number of groups)

        % Extract what depends just on k
        seqk=kk(k);
        CnsampAllk=CnsampAll{seqk};
        gRandNumbForNiniAllk=gRandNumbForNiniAll{seqk};

        if ~isempty(nsampSimData)
            CnsampAllkSimData=CnsampAllk(1:nsampSimData,:);
            gRandNumbForNiniAllkSimData=gRandNumbForNiniAllk(:,1:nsampSimData);
        else
            CnsampAllkSimData=CnsampAllk;
            gRandNumbForNiniAllkSimData=gRandNumbForNiniAllk;
        end

        if seqk<maxk
            CnsampAllkplus1=CnsampAll{seqk+1};
            gRandNumbForNiniAllkplus1=gRandNumbForNiniAll{seqk+1};
            if ~isempty(nsampSimData)
                CnsampAllkplus1SimData=CnsampAllkplus1(1:nsampSimData,:);
                gRandNumbForNiniAllkplus1SimData=gRandNumbForNiniAllkplus1(:,1:nsampSimData);
            else
                CnsampAllkplus1SimData=CnsampAllkplus1;
                gRandNumbForNiniAllkplus1SimData=gRandNumbForNiniAllkplus1;
            end

        else
            % parfor  needs that these variables are initialized
            CnsampAllkSimData=1;
            gRandNumbForNiniAllkplus1SimData=1;
            CnsampAllkplus1SimData=1;

        end

        if msg==1
            disp(['Bands k=' num2str(seqk)])
        end

        for j=1:lalpha

            ktrue = length(PiVal{k, j});
            Mutrue = MuVal{k, j};
            Mutrue=Mutrue(1:ktrue,:);
            Sigmatrue = SigmaVal{k, j};
            Sigmatrue = Sigmatrue(:,:,1:ktrue);
            Pitrue=PiVal{k, j};
            alphaTrimj=alphaTrim(j);
            ngood=round(n*(1-alphaTrimj));
            nout=n-ngood;
            idxkj=IDX{k,j};


            % if outliersFromUniform is false the outliers in the replicate
            % samples are the units which have been trimmed
            if outliersFromUniform == false
                Yadd=Y(idxkj==0,:);
            else
                Yadd=[];
            end

            if usepriorSol ==true
                if seqk<maxk
                    idxkplus1j=IDX{k+1,j};
                else
                    idxkplus1j=[];
                end
            else
                idxkj=[];
                idxkplus1j=[];
            end

            parfor (zz = 1:nsimul, numpool)
                %   for zz = 1:nsimul
                if outliersFromUniform == true
                    [Ysim]=simdataset(ngood, Pitrue, Mutrue, Sigmatrue,'noiseunits', nout);
                    if size(Ysim,1)<n
                        Yadd1=repmat(Ysim(end,:),n-size(Ysim,1),1);
                        Ysim=[Ysim;Yadd1];
                    end
                else
                    [Ysim]=simdataset(ngood, Pitrue, Mutrue, Sigmatrue,'noiseunits', 0);
                    Ysim=[Ysim;Yadd];
                end
                outtcSIM=tclust(Ysim,seqk,alphaTrimj,restrfactor,'nsamp',CnsampAllkSimData,'plots',0,'msg',0,'mixt',mixt, ...
                    'restrtype',restr,'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                    'reftol',reftol,'RandNumbForNini',gRandNumbForNiniAllkSimData,'cshape',cshape,...
                    'priorSol',idxkj);

                if seqk<maxk
                    outtcSIMkplus1=tclust(Ysim,seqk+1,alphaTrimj,restrfactor,'nsamp',CnsampAllkplus1SimData,'plots',0,'msg',0,'mixt',mixt, ...
                        'restrtype',restr,'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                        'reftol',reftol,'RandNumbForNini',gRandNumbForNiniAllkplus1SimData,'cshape',cshape, ...
                        'priorSol',idxkplus1j);
                    BandsCTLtest(k,j,zz)=outtcSIMkplus1.obj-outtcSIM.obj;
                end

                CTLbands(k, j, zz) = outtcSIM.obj;
            end
        end
    end
    out.CTLbands=CTLbands;

    CTLlikLB = zeros(lkk,lalpha);
    CTLlik050 = CTLlikLB;
    CRLlikUB =CTLlikLB;
    for k=1:lkk  % loop for different values of k (number of groups)
        parfor j=1:lalpha
            CTLlikLB(k,j) = quantile(CTLbands(k,j,:), gamma);
            CTLlik050(k,j) = median(CTLbands(k,j,:));
            CRLlikUB(k,j) = quantile(CTLbands(k,j,:), 1- gamma);
        end
    end

    out.CTLlikLB=CTLlikLB;
    out.CTLlik050=CTLlik050;
    out.CTLlikUB=CRLlikUB;

    % Call routine which computes the best tentative solutions using method
    % based on the intersection of confidence bands.
    [CTLtentSol,optimalK,alphafin,CTLidxOptimal]=findOptimalSolutions(CRLlikUB,CTLlikLB,CTLlik050,IDX,alphaTrim,lkk,kk);

    if valSolution == true
        % Validate the groups in correspondence of the best solution
        seq=1:n;
        ExtraZeros=false;
        for jj=1:optimalK
            seqjj=seq(CTLidxOptimal==jj);
            VALjj=FSM(Y(seqjj,:),'msg',0,'plots',0);
            if ~isnan(VALjj.outliers)
                % Set to 0 the units declared as outliers in each
                % component
                CTLidxOptimal(seqjj(VALjj.outliers))=0;
                ExtraZeros=true;
            end
        end
        if  ExtraZeros==true
            alphafin=sum(CTLidxOptimal==0)/n;
        end
    end


    out.CTLoptimalIDX=CTLidxOptimal;
    out.CTLoptimalAlpha=alphafin;
    out.CTLoptimalK=optimalK;
    out.CTLtentSol=CTLtentSol;
    % Store best classification

    % PART BELOW IS ASSOCIATED WITH LIKELIHOOD RATIO TEST (LRT)

    % Values of the  difference between target function using k+1 and k
    % using real data
    tobs= out.CTL(2:end,:)-out.CTL(1:end-1,:);

    % tboot = values of the difference between target function using k+1
    % and k groups for data simulated assuming k groups
    tboot=BandsCTLtest;
    tbootGTtobs=sum(tboot>tobs,3)/nsimul;

    varnam=strcat('alpha=',string(alphaTrim'));
    rownam=strcat('k=',string((kk(1:end-1)')),'_vs_k=',string((kk(2:end)')));
    out.LRTpval=array2table(tbootGTtobs,'VariableNames',varnam,'RowNames',rownam);

    LRTtentSol=zeros(lkk-1,8);
    LRTtentSolIDX=zeros(n,lkk-1);

    crit= 0.02;
    ij=1;

    for j=1:lalpha
        % For each value of alpha soli and indexSpuriousSolution are reset
        indexSpuriousSolution=true;
        soli=0;
        increment=1;
        while indexSpuriousSolution==true

            soli=find(tbootGTtobs((soli+increment):end,j)>=crit,1)+soli+increment-1;
            if ~isempty(soli)

                % Store solution number, value of k, value of alpha
                idxij=IDX{soli,j};
                LRTtentSolIDX(:,ij)=idxij;
                tabij=tabulate(idxij(idxij>0));
                if min(tabij(:,2))<n*max(alphaTrim)
                    properSize=false;
                else
                    properSize=true;
                end
                LRTtentSol(ij,[1:3 5:6 8])=[ij kk(soli) alphaTrim(j) soli j properSize];

                if soli<lkk

                    % sum(tbootGTtobs(soli+1:end,j)>=crit)==lkk-soli
                    if all(tbootGTtobs(soli+1:end,j)>=crit)
                        LRTtentSol(ij,4)=1;
                        indexSpuriousSolution=false;
                    else
                        LRTtentSol(ij,4)=0;
                        increment=find(tbootGTtobs(soli+1:end,j)<crit,1);

                    end
                    % Check whether this solution had been previously found
                    % that is if for the same value of k  we had already
                    % obtained the same solution
                    if ij>1
                        if LRTtentSol(ij,2) == LRTtentSol(ij-1,2) && LRTtentSol(ij,4)==LRTtentSol(ij-1,4)
                            LRTtentSol(ij,:)=[];
                            LRTtentSolIDX(:,ij)=[];
                            ij=ij-1;
                        end
                    end
                else
                    LRTtentSol(ij,4)=1;
                    indexSpuriousSolution=false;
                end

                % idxLR(:,ij)=IDX{i,soli};
                ij=ij+1;
            else
                indexSpuriousSolution=false;
            end
        end
    end

    if size(LRTtentSol,1)>=1
        LRTtentSol=LRTtentSol(1:ij-1,:);

        numsol=(1:size(LRTtentSol,1))';
        nsoleti="Sol"+numsol;

        % FOR EACH VALUE OF ALPHA COMPUTE ADDITIONAL LIK RATIO TESTS TO FIND
        % UNIQUE BEST SOLUTION
        % For example suppose that for a particular alpha we found that there
        % are the potential solutions k=3, k=6 and k=10 then LR test between
        % k=3 and k=6 is performed, and the best among these two values of k
        % is stored and is later tested against k=10
        % The best value among these 3 candidates gets a value of 1 in the last
        % column of matrix LRTtentSolt
        %         nsimulExtra=50;
        %         nsampExtra=5000;
        %         usepriorSolExtra=false;

        for j=1:lalpha
            alphaTrimj=alphaTrim(j);
            MultSolalpha=find(LRTtentSol(:,3)==alphaTrimj);
            kTentative=LRTtentSol(MultSolalpha,2);
            kTentativepos=LRTtentSol(MultSolalpha,5);

            if length(MultSolalpha)>1

                indbestk=1;
                for jjj=2:length(MultSolalpha)

                    % kpos = position of first value of k to use in the additional LR test
                    kpos=kTentativepos(jjj-1);
                    % kH1pos = position larger value of k to use in the additional LR test
                    kH1pos=kTentativepos(jjj);
                    % k and kH1= true values of k  to use in the additional
                    % LR test
                    k=kTentative(jjj-1);
                    kH1=kTentative(jjj);

                    ktrue = length(PiVal{kpos, j});
                    Mutrue = MuVal{kpos, j};
                    Mutrue=Mutrue(1:ktrue,:);
                    Sigmatrue = SigmaVal{kpos, j};
                    Sigmatrue = Sigmatrue(:,:,1:ktrue);
                    Pitrue=PiVal{kpos, j};
                    alphaTrimj=alphaTrim(j);
                    ngood=round(n*(1-alphaTrimj));
                    nout=n-ngood;
                    idxkj=IDX{kpos,j};

                    % if outliersFromUniform is false the outliers in the replicate
                    % samples are the units which have been trimmed
                    if outliersFromUniform == false
                        Yadd=Y(idxkj==0,:);
                    else
                        Yadd=[];
                    end
                    tboot=zeros(nsimulExtra,1);


                    if usepriorSolExtra ==true
                        idxkH1=IDX{kH1pos,j};
                    else
                        idxkj=[];
                        idxkH1=[];
                    end
                    parfor (zz = 1:nsimulExtra, numpool)
                        % Alternative instruction using traditional for
                        % instead of parfor
                        %   for zz = 1:nsimul
                        if outliersFromUniform == true
                            [Ysim]=simdataset(ngood, Pitrue, Mutrue, Sigmatrue,'noiseunits', nout);
                            if size(Ysim,1)<n
                                Yadd1=repmat(Ysim(end,:),n-size(Ysim,1),1);
                                Ysim=[Ysim;Yadd1];
                            end
                        else
                            [Ysim]=simdataset(ngood, Pitrue, Mutrue, Sigmatrue,'noiseunits', 0);
                            Ysim=[Ysim;Yadd];
                        end

                        outtcSIM=tclust(Ysim,k,alphaTrimj,restrfactor,'plots',0,'msg',0,'mixt',mixt, ...
                            'restrtype',restr,'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                            'reftol',reftol,'cshape',cshape,...
                            'priorSol',idxkj,'nsamp',nsampExtra);

                        outtcSIMkH1=tclust(Ysim,kH1,alphaTrimj,restrfactor,'plots',0,'msg',0,'mixt',mixt, ...
                            'restrtype',restr,'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                            'reftol',reftol,'cshape',cshape, ...
                            'priorSol',idxkH1,'nsamp',nsampExtra);

                        tboot(zz)=outtcSIMkH1.obj-outtcSIM.obj;

                    end


                    tobs= out.CTL(kH1pos,j)-out.CTL(kpos,j);

                    % tboot = values of the difference between target function using kH1
                    % and k groups for data simulated assuming k groups
                    addLRtest=sum(tboot>tobs)/nsimulExtra;

                    if addLRtest>crit
                        kTentative(jjj)=kTentative(jjj-1);
                    else
                        indbestk=indbestk+1;
                    end
                end
                LRTtentSol(MultSolalpha(indbestk),7)=1;
            elseif length(MultSolalpha) ==1
                LRTtentSol(MultSolalpha,7)=1;
            else
            end
        end

        LRTtentSolt=array2table(LRTtentSol,'RowNames',nsoleti, ...
            'VariableNames',{'index' 'k' 'alpha' 'Truesol' 'kindex' 'alphaindex' 'kbestGivenalpha' 'ProperSize'});

    else
        LRTtentSol=[];
        LRTtentSolt=[];
        LRTtentSolIDX=[];
    end

    out.LRTtentSol=LRTtentSol;
    out.LRTtentSolt=LRTtentSolt;
    out.LRTtentSolIDX=LRTtentSolIDX;

    % Find  LRToptimalK, LRToptimalAlpha and LRToptimalIDX
    if ~isempty(LRTtentSol)
        % Just in case there are multiple solutions take as optimal the one
        % with the smallest alpha associated with a minimum group size
        if size(LRTtentSol,1)>1
            indexBestSolLRT=find(LRTtentSol(:,7)==1 & LRTtentSol(:,8)==1);
        else
            indexBestSolLRT=1;
        end
        LRToptimalK=LRTtentSol(indexBestSolLRT,2);
        LRToptimalAlpha=LRTtentSol(indexBestSolLRT,3);
        LRToptimalIDX=LRTtentSolIDX(:,indexBestSolLRT);
    else
        LRToptimalK=[];
        LRToptimalAlpha=[];
        LRToptimalIDX=[];
    end
    out.LRToptimalK=LRToptimalK;
    out.LRToptimalAlpha=LRToptimalAlpha;
    out.LRToptimalIDX=LRToptimalIDX;
end

% Close pool and show messages if required
if cleanpool==true
    delete(gcp);
end


% thresh=0.05;
% for j=1:lalpha
%     nextj=false;
%     for i=1:lkk-1
%         if tbootGTtobs(i,j)>thresh
%             tbootGTtobsf(i+1:end,j)=NaN;
%             nextj=true;
%         end
%         if nextj==true
%             break
%         end
%     end
%
% end

out.kk=kk;
out.alpha=alphaTrim;
out.restrfactor=restrfactor;
out.Y=Y;

%% Plot section

if plots==1
    figure
    if ComputeBands==1
        linetype1 = repmat({'-.','-.','-.','-.','-.'},1,10);
        color = repmat({'r','g','b','c','k'},1,10);
        LineWidth = 1;
        hold('on')
        for i = 1:length(kk)
            plot(alphaTrim,CTLlikLB(i,:), 'LineStyle',linetype1{i}, 'Color', color{i}, 'LineWidth', LineWidth)
            % plot(alphaTrim,lik050(i,:), 'LineStyle',linetype{i}, 'Color', color{i}, 'LineWidth', LineWidth)
            text(alphaTrim(end),CTLlik050(i,end),[' k = ' num2str(kk(i))],'FontSize',16, 'Color', color{i})
            plot(alphaTrim,CRLlikUB(i,:), 'LineStyle',linetype1{i}, 'Color', color{i}, 'LineWidth', LineWidth)
        end
    else
        plot(alphaTrim', CTLVal')
        one=ones(length(alphaTrim),1);
        for i = 1:length(kk)
            text(alphaTrim',CTLVal(i,:)',num2str(kk(i)*one),'FontSize', 14)
        end
    end
    xlabel('Trimming level alpha')
    ylabel('Log likelihood')
    set(gca,'XTick',alphaTrim);
end

end

function [TentSol,kfin,alphafin,idxOptimal]=findOptimalSolutions(likUB,likLB,lik050,IDX,alphaTrim,lkk,kk)
conv=0;
% First column of TentSol will contain the value of k while the second
% column the associated trimming level
% Third column is the index of the vector alpha containing the best
% optimal trimming level
% Fourth column is the index of the vector kk containing the best number of
% groups
TentSol=zeros(lkk-1,4);
jj=1;

for j = 1:lkk-1
    alphaBest='';
    for jalpha =1:size(likUB,2)
        if (likUB(j,jalpha) > likLB(j+1,jalpha) &&  lik050(j+1,jalpha)> lik050(j,jalpha)) || ...
                (likLB(j,jalpha) < likUB(j+1,jalpha) &&  lik050(j+1,jalpha) < lik050(j,jalpha))
            % Find best trimming level
            alphaBest = alphaTrim(jalpha);
            break
        else
        end

    end

    if ~isempty(alphaBest)
        TentSol(jj,:)=[kk(j) alphaBest, jalpha j];
        jj=jj+1;
    end
end
if ~isempty(TentSol) && TentSol(1,1)>0
    conv=1;
    TentSol=TentSol(1:jj-1,:);
    kfin=TentSol(1,1);
    jalpha=TentSol(1,3); % index of alphaTrim associated to best solution
    jk=TentSol(1,4); % index of kk associated to best solution
    alphafin=TentSol(1,2);
else
    TentSol=NaN;
end

if conv == 1
    idxOptimal=IDX{jk,jalpha};
else
    disp('No intersection among the curves has been found for the selected trimming levels and number of groups')
    disp('Please increase k or alpha')
    alphafin=max(alphaTrim);
    kfin=max(kk);
    idxOptimal = IDX{end, end};
end
end

%FScategory:CLUS-RobClaMULT

