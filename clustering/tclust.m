function [out , varargout]  = tclust(Y,k,alpha,restrfactor,varargin)
%tclust computes trimmed clustering with scatter restrictions
%
%<a href="matlab: docsearchFS('tclust')">Link to the help function</a>
%
%   tclust partitions the points in the n-by-v data matrix Y
%   into k clusters.  This partition minimizes the trimmed sum, over all
%   clusters, of the within-cluster sums of point-to-cluster-centroid
%   distances.  Rows of Y correspond to points, columns correspond to
%   variables. tclust returns inside structure out an n-by-1 vector idx
%   containing the cluster indices of each point.  By default, tclust uses
%   (squared), possibly constrained, Mahalanobis distances.
%
%
%  Required input arguments:
%
%            Y: Input data. Matrix. Data matrix containing n observations on v variables.
%               Rows of Y represent observations, and columns
%               represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values will
%               automatically be excluded from the computations.
%
%            k: Number of groups. Scalar.
%               Scalar which specifies the number of groups.
%
%        alpha: Global trimming level. Scalar. alpha is a scalar between 0 and 0.5
%               or an integer specifying the number of observations which have to
%               be trimmed. If alpha=0 tclust reduces to traditional model
%               based or mixture clustering (mclust): see Matlab function
%               gmdistribution.
%               More in detail, if 0< alpha <1 clustering is based on
%                h=fix(n*(1-alpha)) observations,
%               else if alpha is an integer greater than 1 clustering is
%               based on h=n-floor(alpha);
%
%  restrfactor: Restriction factor. Scalar. Positive scalar which
%               constrains the allowed differences
%               among group scatters. Larger values imply larger differences of
%               group scatters. On the other hand a value of 1 specifies the
%               strongest restriction forcing all eigenvalues/determinants
%               to be equal and so the method looks for similarly scattered
%               (respectively spherical) clusters. The default is to apply
%               restrfactor to eigenvalues. In order to apply restrfactor
%               to determinants it is necessary to use optional input
%               argument restr.
%
%  Optional input arguments:
%
%       nsamp : Number of subsamples to extract.
%               Scalar or matrix.
%               If nsamp is a scalar it contains the number of subsamples
%               which will be extracted. If nsamp=0
%               all subsets will be extracted.
%               Remark - if the number of all possible subset is <300 the
%               default is to extract all subsets, otherwise just 300
%               If nsamp is a matrix it contains in the rows the indexes of
%               the subsets which have to be extracted. nsamp in this case
%               can be conveniently generated  by function subsets. nsamp can
%               have k columns or k*(v+1) columns. If nsamp has k columns
%               the k initial centroids each iteration i are given by
%               X(nsamp(i,:),:) and the covariance matrices are equal to the
%               identity.
%               If nsamp has k*(v+1) columns the initial centroids and covariance
%               matrices in iteration i are computed as follows:
%               X1=X(nsamp(i,:),:);
%               mean(X1(1:v+1,:)) contains the initial centroid for group  1;
%               cov(X1(1:v+1,:)) contains the initial cov matrix for group 1;
%               mean(X1(v+2:2*v+2,:)) contains the initial centroid for group 2;
%               cov((v+2:2*v+2,:)) contains the initial cov matrix for group 2;
%               ...;
%               mean(X1((k-1)*v+1:k*(v+1))) contains the initial centroids for group k;
%               cov(X1((k-1)*v+1:k*(v+1))) contains the initial cov matrix for group k.
%               REMARK - if nsamp is not a scalar option option below
%               startv1 is ignored. More precisely, if nsamp has k columns
%               startv1=0 elseif nsamp has k*(v+1) columns option startv1
%               =1.
%                 Example - 'nsamp',1000
%                 Data Types - double
%    refsteps : Number of refining iterations. Scalar. Number of refining
%               iterations in each subsample. Default is 15.
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'refsteps',10
%                 Data Types - single | double
%     reftol  : Tolerance for the refining steps. Scalar.
%               The default value is 1e-14;
%                 Example - 'reftol',1e-05
%                 Data Types - single | double
%equalweights : Cluster weights in the concentration and assignment steps.
%               Logical. A logical value specifying whether cluster weights
%               shall be considered in the concentration, assignment steps
%               and computation of the likelihood.
%               if equalweights = true we are (ideally) assuming equally
%               sized groups by maximizing:
%
%                 \[
%                 \sum_{j=1}^k \sum_{ x_i \in group_j } \log f(x_i; m_j , S_j)
%                 \]
%
%               else if equalweights = false (default) we allow for
%               different group weights by maximizing:
%
%                 \[
%                     \sum_{j=1}^k  \sum_{ x_i \in group_j }  \log \left[ \frac{n_j}{n}  f(x_i; m_j , S_j) \right]=
%                 \]
%                 \[
%                   = \sum_{j=1}^k n_j \log n_j/n + \sum_{j=1}^k \sum_{ x_i \in group_j} \log f(x_i; m_j , S_j) .
%                 \]
%
%               Remark: $\sum_{j=1}^k n_j \log n_j/n$ is the so called entropy
%               term
%                 Example - 'equalweights',true
%                 Data Types - Logical
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
% plots    :    Plot on the screen. Scalar or character or structure.
%
%               Case 1: plots option used as scalar.
%               - If plots=0 (default), plots are not generated.
%               - If plots=1, a plot with the classification is shown on
%                 the screen (using the spmplot function). The plot can be:
%                   * for v=1, an histogram of the univariate data.
%                   * for v=2, a bivariate scatterplot.
%                   * for v>2, a scatterplot matrix generated by spmplot.
%               When v>=2 plots offers the following additional features
%               (for v=1 the behaviour is forced to be as for plots=1):
%
%               Case 2: plots option used as character array.
%               - plots='contourf' adds in the background of the bivariate
%                 scatterplots a filled contour plot. The colormap of the
%                 filled contour is based on grey levels as default.
%                 This argument may also be inserted in a field named 'type'
%                 of a structure. In the latter case it is possible to
%                 specify the additional field 'cmap', which changes the
%                 default colors of the color map used. The field 'cmap'
%                 may be a three-column matrix of values in the range [0,1]
%                 where each row is an RGB triplet that defines one color.
%                 Check the colormap function for additional informations.
%               - plots='contour' adds in the background of the bivariate
%                 scatterplots a contour plot. The colormap of the contour
%                 is based on grey levels as default. This argument may
%                 also be inserted in a field named 'type' of a structure.
%                 In the latter case it is possible to specify the additional
%                 field 'cmap', which changes the default colors of the
%                 color map used. The field 'cmap' may be a three-column
%                 matrix of values in the range [0,1] where each row is an
%                 RGB triplet that defines one color.
%                 Check the colormap function for additional informations.
%               - plots='ellipse' superimposes confidence ellipses to
%                 each group in the bivariate scatterplots. The size of the
%                 ellipse is chi2inv(0.95,2), i.e. the confidence level used
%                 by default is 95%. This argument may also be inserted in
%                 a field named 'type' of a structure. In the latter case it
%                 is possible to specify the additional field 'conflev',
%                 which specifies the confidence level to use and it is a
%                 value between 0 and 1.
%               - plots='boxplotb' superimposes on the bivariate scatterplots
%                 the bivariate boxplots for each group, using the boxplotb
%                 function. This argument may also be inserted in a field
%                 named 'type' of a structure.
%
%               Case 3: plots option used as struct.
%                 If plots is a structure it may contain the following fields:
%                 plots.type = in this case the 'type' supplied
%                 is used to set the type of plot as when plots option is
%                 a character array. Therefore, plots.type can be:
%                 'contourf', 'contour', 'ellipse' or 'boxplotb'.
%                 plots.cmap = this field is used to set a colormap
%                 for the plot type. For example, plots.cmap = 'autumn'.
%                 See the MATLAB help of colormap for a list of colormap
%                 possiblilites.
%                 plots.conflev = this is the confidence level
%                 for the confidence ellipses. It must me a scalar between
%                 0 and 1. For example, one can set:
%                 plots.type = 'ellipse';
%                 plots.conflev = 0.5;
%
%               REMARK - The labels=0 are automatically excluded from the
%                          overlaying phase, considering them as outliers.
%                   Example - 'plots', 1
%                   Data Types - single | double | character | struct
%        msg  : Level of output to display. Scalar.
%               Scalar which controls whether to display or not messages
%               on the screen.
%               If msg=0 nothing is displayed on the screen.
%               If msg=1 (default) messages are displayed
%               on the screen about estimated time to compute the estimator
%               or the number of subsets in which there was no convergence.
%               If msg=2 detailed messages are displayed. For example the
%               information at iteration level.
%                   Example - 'msg',1
%                   Data Types - single | double
%      nocheck: Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix Y.
%               As default nocheck=0.
%                   Example - 'nocheck',1
%                   Data Types - single | double
%      startv1: How to initialize centroids and covariance matrices. Scalar.
%               If startv1 is 1 then initial centroids and covariance
%               matrices are based on (v+1) observations randomly chosen,
%               else each centroid is initialized taking a random row of
%               input data matrix and covariance matrices are initialized
%               with identity matrices. The default value of startv1 is 1.
%               Remark 1 - in order to start with a routine which is in the
%               required parameter space, eigenvalue restrictions are
%               immediately applied.
%               Remark 2 - option startv1 is used just if nsamp is a scalar
%               (see for more details the help associated with nsamp).
%                   Example - 'startv1',1
%                   Data Types - single | double
% RandNumbForNini: Pre-extracted random numbers to initialize proportions.
%                Matrix. Matrix with size k-by-size(nsamp,1) containing the
%                random numbers which are used to initialize the
%                proportions of the groups. This option is effective just
%                if nsamp is a matrix which contains pre-extracted
%                subsamples. The purpose of this option is to enable to
%                user to replicate the results in case routine tclust is
%                called using a parfor instruction (as it happens for
%                example in routine IC, where tclust is called through a
%                parfor for different values of the restriction factor).
%                The default value of RandNumbForNini is empty that is
%                random numbers from uniform are used.
%                   Example - 'RandNumbForNini',''
%                   Data Types - single | double
%     restrtype : type of restriction. Character. The type of restriction to
%               be applied on the cluster scatter
%               matrices. Valid values are 'eigen' (default), or 'deter'.
%               eigen implies restriction on the eigenvalues while deter
%               implies restriction on the determinants. If restrtype is
%               'deter' it is also possible to specify  through optional
%               parameter cshape the constraint to apply to the shape
%               matrices.
%                 Example - 'restrtype','deter'
%                 Data Types - char
%       cshape  : constraint to apply to the shape matrices. Scalar greater or
%               equal 1. This options only works is 'restrtype' is
%               'deter'.
%               When restrtype is deter the default value of the "shape" constraint (as
%               defined below) applied to each group is fixed to
%               $c_{shape}=10^{10}$, to ensure the procedure is (virtually)
%               affine equivariant. In other words, the decomposition or the
%               $j$-th scatter matrix $\Sigma_j$ is
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
%       Ysave : Save original input matrix. Scalar. Set Ysave to 1 to
%               request that the input matrix Y
%               is saved into the output structure out. Default is 0, id
%               est no saving is done.
%                 Example - 'Ysave',1
%                 Data Types - single | double
%
%
%  Output:
%
%         out:   structure which contains the following fields
%
%            out.muopt= k-by-v matrix containing cluster centroid
%                       locations. Robust estimate of final centroids of
%                       the groups.
%         out.sigmaopt= v-by-v-by-k array containing estimated constrained
%                       covariance for the k groups.
%            out.idx  = n-by-1 vector containing assignment of each unit to
%                       each of the k groups. Cluster names are integer
%                       numbers from 1 to k. 0 indicates trimmed
%                       observations.
%            out.siz  = Matrix of size (k+1)-by-3.
%                       1st col = sequence from 0 to k;
%                       2nd col = number of observations in each cluster;
%                       3rd col = percentage of observations in each
%                       cluster;
%                       Remark: 0 denotes unassigned units.
%         out.postprob = n-by-k matrix containing posterior probabilities
%                       out.postprob(i,j) contains posterior probabilitiy of unit
%                       i from component (cluster) j. For the trimmed units
%                       posterior probabilities are 0.
%             out.emp = "Empirical" statistics computed on final classification.
%                       Scalar or structure. When convergence is reached,
%                       out.emp=0. When convergence is not obtained, this
%                       field is a structure which contains the statistics
%                       of interest: idxemp (ordered from 0 to k*, k* being
%                       the number of groups with at least one observation
%                       and 0 representing the possible group of outliers),
%                       muemp, sigmaemp and sizemp, which are the empirical
%                       counterparts of idx, muopt, sigmaopt and siz.
%          out.MIXMIX = BIC which uses parameters estimated using the
%                       mixture loglikelihood and the maximized mixture
%                       likelihood as goodness of fit measure.
%                       Remark: this output is present just if input option
%                       mixt is >0.
%          out.MIXCLA = BIC which uses the classification likelihood based on
%                       parameters estimated using the mixture likelihood
%                       (In some books this quantity is called ICL).
%                       Remark: this output is present just if input option
%                       mixt is >0.
%          out.CLACLA = BIC which uses the classification likelihood based on
%                       parameters estimated using the classification likelihood.
%                       Remark: this output is present just if input option
%                       mixt is =0.
%       out.notconver = Scalar. Number of subsets without convergence
%              out.bs = k-by-1 vector containing the units forming initial
%                       subset associated with muopt.
%             out.obj = Scalar. Value of the objective function which is minimized
%                       (value of the best returned solution).
%                       If input option mixt >1 the likelihood which is
%                       maximized is a mixture likelihood as follows:
%                       \[
%                       \prod_{i=1}^h  \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j),
%                       \]
%                       else the likelihood which is maximized is a
%                       classification likelihood of the the form:
%                       \[
%                       \prod_{j=1}^k   \prod _{i\in R_j} \pi_j' \phi (y_i; \; \theta_j),
%                       \]
%                       where $R_j$ contains the indexes of the observations which are assigned to group $j$
%                       with the constraint that $\# \bigcup_{j=1}^k
%                       R_j=h$. In the classification likelihood if input
%                       option equalweights is set to 0, then $\pi_j'=1$, $j=1, ...,
%                       k$.
%   out.equalweights  = Logical. It is true if in the clustering procedure
%                       we (ideally) assumed equal cluster weights
%                       else it is false if we allowed for different
%                       cluster sizes.
%               out.h = Scalar. Number of observations that have determined the
%                       centroids (number of untrimmed units).
%          out.fullsol= Column vector of size nsamp which contains the
%                       value of the objective function at the end of the
%                       iterative process for each extracted subsample.
%              out.Y  = Original data matrix Y. The field is present if option
%                       Ysave is set to 1.
%
%  Optional Output:
%
%            C     : Indexes of extracted subsamples. Matrix.
%                    Matrix of size nsamp-by-(v+1)*k containing (in the
%                    rows) the indices of the subsamples extracted for
%                    computing the estimate.
%
% More About:
%
% This iterative algorithm initializes k clusters randomly and performs
% concentration steps in order to improve the current cluster assignment.
% The number of maximum concentration steps to be performed is given by
% input parameter refsteps. For approximately obtaining the global optimum,
% the system is initialized nsamp times and concentration steps are
% performed until convergence or refsteps is reached. When processing more
% complex data sets higher values of nsamp and refsteps have to be
% specified (obviously implying extra computation time). However, if more
% then 10 per cent of the iterations do not converge, a warning message is
% issued, indicating that nsamp has to be increased.
%
% See also: tkmeans, tclustIC, tclusteda
%
% References:
%
%   Garcia-Escudero, L.A., Gordaliza, A., Matran, C. and Mayo-Iscar, A. (2008),
%   A General Trimming Approach to Robust Cluster Analysis. Annals
%   of Statistics, Vol. 36, 1324-1345. [Technical Report available at:
%   www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf]
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('tclust')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % tclust of geyser data using k=3, alpha=0.1 and  restrfactor=10000.
    Y=load('geyser2.txt');
    out=tclust(Y,3,0.1,10000);
%}

%{
    %% Use of 'plots' option as a struct, to produce more complex plots.
    Y=load('geyser2.txt');
    close all
    out=tclust(Y,3,0.1,10000,'plots',1);
    title('plot with all default options','interpreter','LaTex','FontSize',18);

    % default confidence ellipses.
    out=tclust(Y,3,0.1,10000,'plots','ellipse');
    title('default confidence ellipses','interpreter','LaTex','FontSize',18);

    % confidence ellipses specified by the user
    plots.type = 'ellipse';

    plots.conflev = 0.5;
    out=tclust(Y,3,0.1,10000,'plots',plots);
    title('confidence ellipses set to 0.5','interpreter','LaTex','FontSize',18);

    plots.conflev = 0.9;
    out=tclust(Y,3,0.1,10000,'plots',plots);
    title('confidence ellipses set to 0.9','interpreter','LaTex','FontSize',18);

    % contour plots.
    out=tclust(Y,3,0.1,10000,'plots','contour');
    title('contour plot','interpreter','LaTex','FontSize',18);

    % filled contour plots with additional options
    plots.type = 'contourf';
    plots.cmap = autumn;
    out=tclust(Y,3,0.1,10000,'plots',plots);
    title('contourf plot with autumn colormap','interpreter','LaTex','FontSize',18);

    cascade
%}

%{
    % tclust of geyser with varargout.
    Y=load('geyser2.txt');
    nsamp=20;
    [out,MatrixContSubsets]=tclust(Y,3,0.1,10000,'nsamp',nsamp);
    % MatrixContSubsets is a matrix containing in the rows the indexes of
    % the nsamp subsets which have been extracted
%}

%{
    %% tclust of geyser data (output comparison).
    % We compare the output using three different values of
    % restriction factor.
    close all
    Y=load('geyser2.txt');
    restrfactor=10000;
    % nsamp = number of subsamples which will be extracted
    nsamp=500;
    out=tclust(Y,3,0.1,restrfactor,'nsamp',nsamp,'plots','ellipse');
    title(['Restriction factor =' num2str(restrfactor)],'interpreter','LaTex','FontSize',18)
    restrfactor=10;
    out=tclust(Y,3,0.1,restrfactor,'nsamp',nsamp,'refsteps',10,'plots','ellipse');
    title(['Restriction factor =' num2str(restrfactor)],'interpreter','LaTex','FontSize',18)
    % trimmed k-means solution restrfactor=1
    restrfactor=1;
    out=tclust(Y,3,0.1,restrfactor,'nsamp',nsamp,'refsteps',10,'plots','ellipse');
    title(['Restriction factor =' num2str(restrfactor) '. Trimmed k-means solution'],'interpreter','LaTex','FontSize',18)
    cascade
%}

%{
    %  tclust applied to the M5data.
    %  A bivariate data set obtained from three normal bivariate
    %  distributions with different scales and proportions 1:2:2. One of
    %  the components is very overlapped with another one. A 10 per cent
    %  background noise is added uniformly distributed in a rectangle
    %  containing the three normal components and not very overlapped with
    %  the three mixture components. A precise description of the M5 data
    %  set can be found in Garcia-Escudero et al. (2008).
    close all
    Y=load('M5data.txt');
    % plot(Y(:,1),Y(:,2),'o')
    % Scatter plot matrix with univariate boxplot on the main diagonal
    spmplot(Y(:,1:2),Y(:,3),[],'box')

    out=tclust(Y(:,1:2),3,0,1000,'nsamp',100,'plots',1)
    out=tclust(Y(:,1:2),3,0,10,'nsamp',100,'plots',1)
    out=tclust(Y(:,1:2),3,0.1,1,'nsamp',1000,'plots',1,'equalweights',1)
    out=tclust(Y(:,1:2),3,0.1,1000,'nsamp',100,'plots',1)

    cascade
%}

%{
    % tclust in presence of structured noise.
    % The data have been generated using the following R instructions
    %    set.seed (0)
    %    library(MASS)         
    %    v <- runif (100, -2 * pi, 2 * pi)
    %    noise <- cbind (100 + 25 * sin (v), 10 + 5 * v)
    %
    %    x <- rbind (
    %        rmvnorm (360, c (0.0,  0), matrix (c (1,  0,  0, 1), ncol = 2)),
    %        rmvnorm (540, c (5.0, 10), matrix (c (6, -2, -2, 6), ncol = 2)),
    %        noise)
    
    close all
    Y=load('structurednoise.txt');
    out=tclust(Y(:,1:2),2,0.1,100,'plots',1)
    out=tclust(Y(:,1:2),5,0.15,1,'plots',1)
    cascade
%}

%{
    % tclust applied to mixture100 data.
    % The data have been generated using the following R instructions
    %     set.seed (100)
    %     mixt <- rbind (rmvnorm (360, c (  0,  0), matrix (c (1,  0,  0,  1), ncol = 2)),
    %                rmvnorm (540, c (  5, 10), matrix (c (6, -2, -2,  6), ncol = 2)),
    %                rmvnorm (100, c (2.5,  5), matrix (c (50, 0,  0, 50), ncol = 2)))
    close all
    Y=load('mixture100.txt');
    out=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1)
    out=tclust(Y(:,1:2),3,0.05,1,'refsteps',20,'plots',1)
    cascade
%}

%{
    % tclust applied to mixture100 data, comparison of different options.
    close all
    Y=load('mixture100.txt');
    % Traditional tclust
    out1=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1)
    title('Traditional tclust','interpreter','LaTex','FontSize',18);
    % tclust with mixture models (selection of untrimmed units according to
    % likelihood contributions
    out2=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1,'mixt',1)
    title('tclust with mixture models (likelihood contributions)','interpreter','LaTex','FontSize',18);
    % Tclust with mixture models (selection of untrimmed units according to
    % densities weighted by estimates of the probability of the components)
    out3=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1,'mixt',2)
    title('tclust with mixture models (probability of the components)','interpreter','LaTex','FontSize',18);
    cascade
%}

%{
    % tclust using simulated data.
    % 5 groups and 5 variables
    n1=100;
    n2=80;
    n3=50;
    n4=80;
    n5=70;
    v=5;
    Y1=randn(n1,v)+5;
    Y2=randn(n2,v)+3;
    Y3=rand(n3,v)-2;
    Y4=rand(n4,v)+2;
    Y5=rand(n5,v);

    group=ones(n1+n2+n3+n4+n5,1);
    group(n1+1:n1+n2)=2;
    group(n1+n2+1:n1+n2+n3)=3;
    group(n1+n2+n3+1:n1+n2+n3+n4)=4;
    group(n1+n2+n3+n4+1:n1+n2+n3+n4+n5)=5;

    close all
    Y=[Y1;Y2;Y3;Y4;Y5];
    out=tclust(Y,5,0.05,1.3,'refsteps',20,'plots',1)
%}

%{
    %% Automatic choice of the best number of groups for geyser data.
    % close all
    Y=load('geyser2.txt');
    maxk=6;
    CLACLA=[(1:maxk)' zeros(maxk,1)];
    for j=1:maxk
        out=tclust(Y,j,0.1,5,'msg',0);
        CLACLA(j,2)=out.CLACLA;
    end
 
    MIXCLA=[(1:maxk)' zeros(maxk,1)];
    MIXMIX=MIXCLA;
    for j=1:maxk
        out=tclust(Y,j,0.1,5,'mixt',2,'msg',0);
        MIXMIX(j,2)=out.MIXMIX;
        MIXCLA(j,2)=out.MIXCLA;
    end
    
    subplot(1,3,1)
    plot(CLACLA(:,1),CLACLA(:,2))
    xlim([1 maxk])
    xlabel('Number of groups')
    ylabel('CLACLA')

    subplot(1,3,2)
    plot(MIXMIX(:,1),MIXMIX(:,2))
    xlabel('Number of groups')
    ylabel('MIXMIX')
    xlim([1 maxk])
    
    subplot(1,3,3)
    plot(MIXCLA(:,1),MIXCLA(:,2))
    xlabel('Number of groups')
    ylabel('MIXCLA (ICL)')
    xlim([1 maxk])
%}

%{
    % Automatic choice of the best number of groups for simulated data with
    % k=5 and v=5.
    close all
    n1=100;     % Generate 5 groups in 5 dimensions
    n2=80;
    n3=50;
    n4=80;
    n5=70;
    v=5;
    Y1=randn(n1,v)+5;
    Y2=randn(n2,v)+3;
    Y3=rand(n3,v)-2;
    Y4=rand(n4,v)+2;
    Y5=rand(n5,v);

    group=ones(n1+n2+n3+n4+n5,1);
    group(n1+1:n1+n2)=2;
    group(n1+n2+1:n1+n2+n3)=3;
    group(n1+n2+n3+1:n1+n2+n3+n4)=4;
    group(n1+n2+n3+n4+1:n1+n2+n3+n4+n5)=5;

    Y=[Y1;Y2;Y3;Y4;Y5];
    restrfactor=5;
    maxk=7;
    CLACLA=[(1:maxk)' zeros(maxk,1)];
    for j=1:maxk
        out=tclust(Y,j,0.1,restrfactor);
        CLACLA(j,2)=out.CLACLA;
    end
 
    MIXCLA=[(1:maxk)' zeros(maxk,1)];
    MIXMIX=MIXCLA;
    for j=1:maxk
        out=tclust(Y,j,0.1,restrfactor,'mixt',2);
        MIXMIX(j,2)=out.MIXMIX;
        MIXCLA(j,2)=out.MIXCLA;
    end
    
    subplot(1,3,1)
    plot(CLACLA(:,1),CLACLA(:,2))
    xlabel('Number of groups')
    ylabel('CLACLA')

    subplot(1,3,2)
    plot(MIXMIX(:,1),MIXMIX(:,2))
    xlabel('Number of groups')
    ylabel('MIXMIX')

    subplot(1,3,3)
    plot(MIXCLA(:,1),MIXCLA(:,2))
    xlabel('Number of groups')
    ylabel('MIXCLA (ICL)')
%}

%{
    % tclust applied to Swiss banknotes imposing determinant restriciton.
    close all
    load('swiss_banknotes');
    Y=swiss_banknotes.data;
    out=tclust(Y,3,0.01,20,'restrtype','deter','refsteps',20,'plots',1)
%}

%{
    % tclust applied to the Geyser data imposing determinant restriciton.
    close all
    Y=load('geyser2.txt');
    out=tclust(Y,4,0.1,10,'restrtype','deter','refsteps',20,'plots',1)
%}

%% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);
[n, v]=size(Y);

% User options
% startv1def = default value of startv1 =1, initialization using covariance
% matrices based on v+1 units
startv1def=1;

if nargin>4
    % Check whether option nsamp exists
    chknsamp = strcmp(varargin,'nsamp');
    
    % if the sum below is greater than 0 option nsamp exists
    if sum(chknsamp)>0
        nsamp=cell2mat(varargin(find(chknsamp)+1));
        
        % Check if options nsamp is a scalar
        if ~isscalar(nsamp)
            % if nsamp is not a scalar, it is a matrix which contains in
            % the rows the indexes of the subsets which have to be
            % extracted
            C=nsamp;
            [nsampdef,ncolC]=size(C);
            % The number of rows of nsamp (matrix C) is the number of
            % subsets which have to be extracted
            nselected=nsampdef;
            % If the number of columns of nsamp (matrix C) is equal to v
            % then the procedure is initialized using identity matrices
            % else using covariance matrices based on the (v+1)*k units
            if ncolC==v
                startv1=0;
            elseif ncolC==k*(v+1)
                startv1=1;
            else
                disp('If nsamp is not a scalar it must have v or k*(v+1) columns')
                disp('Please generate nsamp using')
                disp('nsamp=subsets(number_desired_subsets,n,k) or')
                disp('nsamp=subsets(number_desired_subsets,n,(v+1)*k)')
                error('FSDA:tclust:WrongNsamp','Wrong number of columns in matrix nsamp')
            end
            NoPriorSubsets=0;
        else
            % If nsamp is a scalar it simply contains the number of subsets
            % which have to be extracted. In this case NoPriorSubsets=1
            NoPriorSubsets=1;
            
            % In this case (nsamp is a scalar) we check whether the user has supplied option
            % startv1
            chkstartv1 = strcmp(varargin,'startv1');
            if sum(chkstartv1)>0
                startv1= cell2mat(varargin(find(chkstartv1)+1));
            else
                startv1=startv1def;
            end
        end
    else
        % If option nsamp is not supplied then for sure there are no prior
        % subsets
        NoPriorSubsets=1;
        
        % In this case (options nsamp does not exist) we check whether the
        % user has supplied option startv1
        chkstartv1 = strcmp(varargin,'startv1');
        if sum(chkstartv1)>0
            startv1= cell2mat(varargin(find(chkstartv1)+1));
        else
            startv1=startv1def;
        end
    end
else
    % if nargin ==4 for use the user has not supplied prior subsets.
    % Default value of startv1 is used
    NoPriorSubsets=1;
    startv1=startv1def;
end

% If the user has not specified prior subsets (nsamp is not a scalar) than
% according the value of startv1 we have a different value of ncomb
if NoPriorSubsets ==1
    % Remark: startv1 must be immediately checked because the calculation of
    % ncomb is immediately affected.
    
    if startv1 ==1
        ncomb=bc(n,k*(v+1));
    else
        % If the number of all possible subsets is <300 the default is to extract
        % all subsets otherwise just 300.
        % Notice that we use bc, a fast version of nchoosek. One may also use the
        % approximation floor(exp(gammaln(n+1)-gammaln(n-p+1)-gammaln(p+1))+0.5)
        ncomb=bc(n,k);
    end
    nsampdef=min(300,ncomb);
end

refstepsdef=15;
reftoldef=1e-5;

% tolrestreigen = tolerance to use in function restreigen
tolrestreigen=1e-08;

% Default
if nargin<3
    alpha=0.05;
    warning('FSDA:tclust:Wrongalpha','You have not specified alpha: it is set to 0.05 by default');
else
    if isempty(alpha)
        alpha=0.05;
        warning('FSDA:tclust:Wrongalpha','You have not specified alpha: it is set to 0.05 by default');
    end
end

% Fix alpha equal to the trimming size
% h = number of observations which is used to compute the centroids

if alpha<0
    error('FSDA:tclust:WrongAlpha','alpha must a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
end

if nargin<4
    restrfactor=12;
    warning('FSDA:tclust:Wrongrestrfact','You have not specified restrfactor: it is set to 12 by default');
end

% h = number of untrimmed units
if alpha>=1
    h=n-floor(alpha);
else
    h=fix(n*(1-alpha));
end

% restrnum=1 implies eigenvalue restriction
restrnum=1;
% cshape. Constraint on the shape matrices inside each group which works only if restrtype is 'deter'
cshape=10^10;

options=struct('nsamp',nsampdef,'RandNumbForNini','','plots',0,'nocheck',0,...
    'msg',1,'Ysave',0,'refsteps',refstepsdef,'equalweights',false,...
    'reftol',reftoldef,'mixt',0,'startv1',startv1def,'restrtype','eigen','cshape',cshape);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:tclust:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:tclust:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin > 4
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    % And check if the optional user parameters are reasonable.
    
    % Check number of subsamples to extract
    if isscalar(options.nsamp) && options.nsamp>ncomb
        disp('Number of subsets to extract greater than (n k). It is set to (n k)');
        options.nsamp=0;
    elseif  options.nsamp<0
        error('FSDA:tclust:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
    end
    
    % Check restriction factor
    if restrfactor<1
        disp('Restriction factor smaller than 1. It is set to 1 (maximum contraint==>spherical groups)');
        restrfactor=1;
    else
    end
    
    % Default restriction is the one based on the eigenvalues
    % restrnum=1 ==> restriction on the eigenvalues
    % restrnum=2 ==> restriction on the determinants
    % restrnum=3 ==> restriction on both
    
    restr=options.restrtype;
    if strcmp(restr,'eigen')
        restrnum=1;
    elseif strcmp(restr,'deter')==1
        restrnum=2;
        cshape=options.cshape;
        restrfactor=[restrfactor cshape];
        
    else
        error('FSDA:tclust:Wrongrestr','Wrong restriction');
    end
    
end

% Default values for the optional
% parameters are set inside structure 'options'

plots=options.plots;        % Plot of the resulting classification
nsamp=options.nsamp;        % Number of subsets to extract
equalweights=options.equalweights;    % Specify if assignment must take into account the size of the groups

refsteps=options.refsteps;
reftol=options.reftol;

RandNumbForNini=options.RandNumbForNini;
if isempty(RandNumbForNini)
    NoPriorNini=1;
else
    NoPriorNini=0;
end

%Initialize the objective function (trimmed variance) by a
%large  value
vopt=-1e+30;

msg=options.msg;            % Scalar which controls the messages displayed on the screen

mixt=options.mixt;         % if options.mixt==1 mixture model is assumed

if mixt>=1 && equalweights == true
    warning('FSDA:tclust:WrongEqualWeights','option equalweights must be different from 1 if mixture model approach is assumed')
    warning('FSDA:tclust:WrongEqualWeights','options equalweights is reset to 0')
end

%% Combinatorial part to extract the subsamples (if not already supplied by the user)
if NoPriorSubsets ==1
    if startv1 ==1 && k*(v+1) < n
        [C,nselected] = subsets(nsamp,n,k*(v+1),ncomb,msg);
    else
        [C,nselected] = subsets(nsamp,n,k,ncomb,msg);
        niinistart=repmat(floor(h/k),k,1);
    end
end

% Store the indices in varargout
if nargout==2
    varargout={C};
end

% ll = matrix of loglikelihoods for each unit from each cluster
% rows of ll are associated to units
% Columns of ll are associated to clusters
ll=zeros(n,k);
obj=1e+14;

% Create an identity matrix which will be used in fucntion logmvnpdfFS
eyev=eye(v);
% Create a copy of matrix Y which will be used in fucntion logmvnpdfFS
Y0tmp=zeros(n,v);

if mixt>=1
    % log_lh = h-by-k matrix containing the log of component conditional
    %          density weighted by the component probability.
    % log_lh = log( \pi_j \phi (y_i; \; \theta_j))
    log_lh=zeros(h,k);
end

% noconv = scalar linked to the number of times in which there was no
% convergence
noconv=0;

% The covariances are given initially by k identity matrices
ey=eye(v,v);
eyk=repmat(ey,[1 1 k]);

onev1=ones(v,1);

% Lambda_vk = matrix which will contain in column j the v (unrestricted)
% eigevalues of covariance matrix of group j (j=1, ..., k)
Lambda_vk=ones(v,k);

% sigmaopt = 3 dimensional array which will contain the estimates of the
% covariance matrices for the best solution
sigmaopt=zeros(v,v,k);

% Initialise and start timer.
tsampling = ceil(min(nselected/5 , 10));
time=zeros(tsampling,1);


% Lambda will contain the matrix of eigenvalues in each iteration for
% all groups. Lambda is a 3D array of size v-by-v-by-k
% Lambda=sigmaini;
% U will contain the eigenvectors of the cov matrices in each iteration
% for all groups (or the shape matrices if restrnum=2). U is a 3D array of size v-by-v-by-k
sigmaini=zeros(v,v,k);
U=sigmaini;

% fullsol = vector which stores value of the objective function in each
% iteration
fullsol=zeros(nselected,1);

verMatlab=verLessThan('matlab','8.2.0');

if verMatlab ==1
    userepmat=0;
else
    userepmat=1;
end

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
            disp('MixLik with untrimmed units selected using h largest likelihood contributions');
    end
end



%%  Random starts
for i=1:nselected
    
    if msg==1
        if i <= tsampling
            tstart = tic;
        end
    end
    
    if msg == 2
        disp(['Iteration ' num2str(i)])
    end
    
    if startv1 ==1
        if NoPriorNini==1
            randk=rand(k,1);
        else
            randk=RandNumbForNini(:,i);
        end
        
        niini=floor(h*randk/sum(randk));
        
        cini=zeros(k,v);
        for j=1:k
            ilow=(j-1)*(v+1)+1;
            iup=j*(v+1);
            index=C(i,:);
            selj=index(ilow:iup);
            % cini(j,:)=mean(Y(selj,:));
            Yselj=Y(selj,:);
            cini(j,:)=sum(Yselj)/(v+1);
            
            Yseljc = bsxfun(@minus,Yselj,cini(j,:));
            sigmaini(:,:,j) = (Yseljc' * Yseljc) / (v+1);
            % lines above are a faster solution for instruction below
            % sigmaini(:,:,j)=cov(Y(selj,:));
            
            % Eigenvalue eigenvector decomposition for group j
            [Uj,Lambdaj] = eig(sigmaini(:,:,j));
            % Store eigenvectors and eigenvalues of group j
            U(:,:,j)=Uj;
            Lambda_vk(:,j)=diag(Lambdaj);
        end
        
        Lambda_vk(Lambda_vk<0)=0;
        
        if restrnum==1
            
            % Restriction on the eigenvalue
            autovalues=restreigen(Lambda_vk,niini,restrfactor,tolrestreigen,userepmat);
            
        else
            % Restrictions on the determinants
            autovalues=restrdeter(Lambda_vk,niini,restrfactor,tolrestreigen,userepmat);
        end
        
        
        % Covariance matrices are reconstructed keeping into account the
        % constraints on the eigenvalues
        for j=1:k
            sigmaini(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
            
            % Alternative code: in principle more efficient but slower
            % because diag is a built in function
            % sigmaini(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
        end
        
    else
        
        % initialization of niini with equal proportions
        niini=niinistart;
        
        % extract a subset of size v
        index = C(i,:);
        
        % cini will contain the centroids in each iteration
        cini=Y(index,:);
        % sigmaini will contain the covariance matrices in each iteration
        sigmaini=eyk;
    end
    
    % sigmaopt will be final estimate of the covariance matrices
    % sigmaopt=sigmaini;
    
    iter=0;
    mudiff=1e+15;
    
    postprob=0;
    ind=0;
    
    % refsteps "concentration" steps will be carried out
    while ( (mudiff > reftol) && (iter < refsteps) )
        iter = iter + 1;
        if equalweights
            % In this case we are (ideally) assuming equally sized groups
            for j=1:k
                ll(:,j)= logmvnpdfFS(Y,cini(j,:),sigmaini(:,:,j),Y0tmp,eyev,n,v,0);
            end
        else
            
            % In this case we allow for different group weights or we are
            % assuming a mixture model
            for j=1:k
                % REMARK: we use log(niini(j)) instead of log(niini(j)/h)
                % because h is constant
                ll(:,j)= log(niini(j)/h) +  logmvnpdfFS(Y,cini(j,:),sigmaini(:,:,j),Y0tmp,eyev,n,v,0);
                % Line above is faster but equivalent to
                % ll(:,j)= (niini(j)/h)*mvnpdf(Y,cini(j,:),sigmaini(:,:,j));
            end
            
        end
        
        if mixt==2
            
            postprobold=postprob;
            
            [~,postprob,disc]=estepFS(ll);
            
            % Sort the n likelihood contributions
            % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
            [~,qq]=sort(disc,'descend');
            
            % qq = vector of size h which contains the indexes associated with the largest n(1-alpha)
            % (weighted) likelihood contributions
            qqunassigned=qq((h+1):n);
            qq=qq(1:h);
            
            % Ytri = n(1-alpha)-by-v matrix associated with the units
            % which have the largest n(1-alpha) likelihood contributions
            Ytri=Y(qq,:);
            
            postprob(qqunassigned,:)=0;
            
            % M-step update of niini
            % niini = numerator of component probabilities
            niini=(sum(postprob))';
            
        else
            indold=ind;
            
            % In this part we select the untrimmed units.
            % They are those which have the n(1-alpha) largest values among the
            % maxima of each row of matrix ll.
            % vector disc of length(n) contains the (weighted) contribution of
            % each unit to the log likelihood.
            [disc,ind]= max(ll,[],2);
            
            % Sort the n likelihood contributions
            % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
            [~,qq]=sort(disc,'descend');
            
            
            % qq = vector of size h which contains the indexes associated with the largest n(1-alpha)
            % (weighted) likelihood contributions
            qqunassigned=qq(h+1:end);
            qq=qq(1:h);
            
            % Ytri = n(1-alpha)-by-v matrix associated with the units
            % which have the largest n(1-alpha) likelihood contributions
            Ytri=Y(qq,:);
            % Ytriind = grouping indicator vector (of size n(1-alpha))
            % associated to Ytri
            groupind=ind(qq);
            
            % ind is the identifier vector
            % trimmed units have a value of ind=0
            ind(qqunassigned)=0;
        end
        
        
        if mixt == 1
            %  expll=exp(ll(qq,:));
            %  sumll=sum(expll,2);
            %  postprob=bsxfun(@rdivide,expll,sumll);
            
            % E-step: computation of posterior probabilities for untrimmed
            % units. In the context of mixture models posterior
            % probabilities will be used to estimate new component
            % probabilities of the mixtures, new centroids and new
            % covariance matrices
            
            postprobold=postprob;
            [~,postprob]=estepFS(ll);
            
            postprob(qqunassigned,:)=0;
            
            % M-step update of niini
            % niini = numerator of component probabilities
            niini=(sum(postprob))';
            
        end
        
        % M-step: parameters are updated
        % Matrix cini contains estimates of the new k centroids
        % Array sigmaini contains estimates of the new covariance matrices
        
        for j=1:k
            
            if mixt>=1
                % Matrix cini is updated using weighted means. The weights
                % are given by the posterior probabilities.
                % Note that Y is used instead of Ytri because posterior
                % probabilities for unassigned units are 0.
                cini(j,:)= sum(bsxfun(@times, Y, postprob(:,j)),1)/niini(j);
                
                if niini(j)>0
                    Ytric = bsxfun(@minus,Y,cini(j,:));
                    
                    sqweights = postprob(:,j).^(1/2);
                    
                    % Ytric = [X(:,1).*sqweights   X(:,2).*sqweights ...   X(:,end).*sqweights]
                    Ytric = bsxfun(@times, Ytric, sqweights);
                    
                    sigmaini(:,:,j) = (Ytric' * Ytric) / niini(j);
                    
                    % Eigenvalue eigenvector decomposition for group j
                    [Uj,Lambdaj] = eig(sigmaini(:,:,j));
                    
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_vk(:,j)=diag(Lambdaj);
                else
                    sigmaini(:,:,j)=ey;
                    U(:,:,j)=ey;
                    Lambda_vk(:,j)=onev1;
                end
                
            else  % This is the "crisp assignment" setting
                
                % Boolean index of units forming group j
                groupj=groupind==j;
                
                % Size of group j
                niini(j)=sum(groupj);
                
                % Group j values
                Ytrij=Ytri(groupj,:);
                
                % Means of group j
                cini(j,:)=sum(Ytrij)/niini(j);
                
                % niini=sum(Ytri(:,v+1)==j);
                if niini(j)>0
                    % Covariance of group j:
                    % sigmaini(:,:,j)=cov(Ytrij);
                    % cov would recompute the sample means; code below is
                    % more efficient
                    
                    % Important remark: DfM is a mex file with the compiled
                    % code of an efficient method to compute the following
                    % element by element operation:
                    % Ytrij = bsxfun(@minus,Ytrij,cini(j,:));
                    % The mex has been compiled for the following
                    % platforms: 32 and 64 bit MS Windows, 64 bit Linux
                    % and 64 bit MacOS. However, if you experience an error
                    % in correspondence to the DfM execution, you should
                    % comment the DfM line below and uncomment the bsxfun
                    % instruction above. In contexts where this is called
                    % many times, this solution is much more performant.
                    DfM(Ytrij,cini(j,:),Ytrij,niini(j),v);
                    
                    sigmaini(:,:,j) = (Ytrij' * Ytrij) / niini(j);
                    
                    % Eigenvalue eigenvector decomposition for group j
                    [Uj,Lambdaj] = eig(sigmaini(:,:,j));
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_vk(:,j)=diag(Lambdaj);
                    
                else
                    sigmaini(:,:,j)=ey;
                    U(:,:,j)=ey;
                    if restrnum==1
                        Lambda_vk(:,j)=onev1;
                    else
                        Lambda_vk(j)=1;
                    end
                end
                
            end
            
        end
        
        % Lambda_vk is a v-by-k  matrix whose jth column contains the
        % unrestricted eigenvalues of cov. matrix of group j   j=1, ..., k
        % The row below is just to avoid numerical problems
        Lambda_vk(Lambda_vk<0)=0;
        if restrnum==1
            autovalues=restreigen(Lambda_vk,niini,restrfactor,tolrestreigen,userepmat);
            
            
        elseif restrnum==2
            % Restriction on the determinants
            autovalues=restrdeter(Lambda_vk,niini,restrfactor,tolrestreigen,userepmat);
        end
        
        % Covariance matrices are reconstructed keeping into account the
        % constraints of the eigenvalues
        for j=1:k
            sigmaini(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
            
            % Alternative code: in principle more efficient but slower
            % because diag is a built in function
            % sigmaini(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
        end
        
        % BELOW THERE IS AN ALTERNATIVE WAY OF FINDING sigmaini without the loop
        % Note that the implementation below uses mex function mtimes
        %         autov=(autovalues(:)');
        %         if userepmat==1
        %             autov1=repmat(autov,v,1,1);
        %             ULambda=U.*reshape(autov1,v,v,k);
        %         else
        %             ULambda=reshape(bsxfun(@times,reshape(U,v,v*k),autov),v,v,k);
        %         end
        %         sigmainichk=mtimesx(ULambda,U,'T');
        
        
        % Alternative code based on gpuArrary
        % sigmainichk1=pagefun(@mtimes, gpuArray(sigmainichk), gpuArray(Ut));
        
        % Calculus of the objective function (E-step)
        % oldobj=obj;
        obj = 0;
        
        if mixt>=1
            
            % Likelihood for mixture modelling
            
            %   log_lh is a h-by-k matrix where k is the number of Gaussian components of the mixture
            %   log_lh is the log of component conditional density weighted by the component
            %   probability.  The probability of j-th component is niini(j)/h
            %   log_lh(i,j) is log (Pr(point i|component j) * Prob( component j))
            
            for j=1:k
                log_lh(:,j)=  log(niini(j)/h)+logmvnpdfFS(Ytri,cini(j,:),sigmaini(:,:,j),Y0tmp(1:h,:),eyev,h,v,0);
            end
            
            % obj contains the value of the log likelihood for mixture models
            obj=estepFS(log_lh);
            
        else
            
            % Likelihood for  crisp clustering
            for j=1:k
                % disp(ni)
                if niini(j)>0
                    if equalweights
                        % we simply sum the log of the densities for the untrimmed
                        % units
                        obj=obj+ sum(logmvnpdfFS(Ytri(groupind==j,:),cini(j,:),sigmaini(:,:,j)));
                    else
                        % niini(j)*log(niini(j)/h) is the so called entropy
                        % term which allows for different group weights
                        
                        niinij=niini(j);
                        obj=obj+ niini(j)*log(niinij/h)+sum(logmvnpdfFS(Ytri(groupind==j,:),cini(j,:),sigmaini(:,:,j),Y0tmp(1:niinij,:),eyev,niinij,v,0));
                    end
                end
            end
            
        end
        
        if mixt>0
            % if mixt >0 stopping criterion is referred to postprob
            mudiff=sum(sum(abs(postprob-postprobold)))/n;
            % disp(mudiff)
        else
            % if mixt=0 stopping criterion is referred to no modiification in the classification
            mudiff=sum(abs(indold-ind)>0)/n;
            % disp(mudiff)
        end
        
        %disp(num2str(iter))
        %disp(mudiff)
        %
        % Alternative stopping criterion was based  on the relative
        % modification of the objective function.
        %                  mudiff =abs(oldobj-obj)/abs(obj);
        %                  disp(['Iteration ' num2str(iter)])
        %                  disp([oldobj-obj]/abs(obj))
        %                  disp('monit')
        
        if iter==refsteps
            noconv=noconv+1;
        end
        
    end
    
    % Store value of the objective function for iteration i
    fullsol(i)=obj;
    
    % Store the centroids and the value of the objective function
    if obj>=vopt
        % vopt = value of the objective function in correspondence of the
        % best centroids
        vopt=obj;
        % muopt = matrix containing best centroids
        muopt=cini;
        % nopt = vector containing sizes of the groups
        nopt=niini;
        % format long;
        %disp(index)
        %disp(obj)
        
        % sigmaopt
        sigmaopt=sigmaini;
        
        % store the indexes of the subset which gave rise to the
        % optimal solution
        bs=index;
    end
    
    
    if msg==1
        if i <= tsampling
            % sampling time until step tsampling
            time(i)=toc(tstart);
        elseif i==tsampling+1
            % stop sampling and print the estimated time
            fprintf('Total estimated time to complete tclust: %5.2f seconds \n', nselected*median(time));
        end
    end
    
end
notconver=noconv/nselected;
if msg==1
    if notconver>0.1
        disp('------------------------------')
        disp(['Warning: Number of subsets without convergence equal to ' num2str(100*notconver) '%'])
    end
end

%% Store quantities in out structure
%exist('muopt')==0

% Procedure to order the non-empty components
if any(any(isnan(muopt)))
    
    % restore apropriate order of the components
    NanGroups = isnan(muopt(:,1)); % missing components
    
    % order of the components in nopt, muopt and sigmaopt
    nopt = [nopt(~NanGroups); nopt(NanGroups)];
    muopt = [muopt(~NanGroups,:); muopt(NanGroups,:)];
    sigmaopt(:,:,NanGroups) = NaN; % assign NaN on the empty clusters
    sigmaopt = cat(3, sigmaopt(:,:,~NanGroups), sigmaopt(:,:,NanGroups));
end

% With the best obtained values for the parameters, we compute the final
% assignments and parameters

% construct the  log of component conditional density weighted by the
% component probability.
% ll = log( \pi_j \phi (y_i; \; \theta_j))
% Get the likelihood for each point with each component
% ll is a n by k matrix,
% if equalweights is false
% ll(i,j) is log( (n_j/h) * f(x_i|\theta_j))
% else if  equalweights is true
% ll(i,j) is log( f(x_i|\theta_j))
% f(x_i|\theta_j) is multivariate normal with theta_j =(mu_j, \Sigma_j)
if equalweights
    for j=1:k
        if any(~isnan(muopt(j,:)))
            ll(:,j) = logmvnpdfFS(Y,muopt(j,:),sigmaopt(:,:,j),Y0tmp,eyev,n,v,0);
        else
            % avoid the computation for empty components and assign NaN
            ll(:,j) = NaN;
        end
    end
else
    for j=1:k
        if any(~isnan(muopt(j,:)))
            ll(:,j) = log(nopt(j)/h) + logmvnpdfFS(Y,muopt(j,:),sigmaopt(:,:,j),Y0tmp,eyev,n,v,0);
        else
            % avoid the computation for empty components and assign NaN
            ll(:,j) = NaN;
        end
    end
end

% matrix ll forms the input to compute both the MIXTURE and the
% CLASSIFICATION LIKELIHOOD

% postprob n x k containing posterior probabilities
% logpdf n x 1 vector containg the n contributions to the log
% likelihood of mixture models
[~,postprob,logpdf]=estepFS(ll);

% %%%%%
% In this part we select the untrimmed units
% They are those which have the n(1-alpha) largest values among the
% maxima of each row of matrix ll
% vector disc of length(n) contains the (weighted) contribution of
% each unit to the log likelihood
% idx = n x 1 vector containing the final assignments
% disc = n x 1 vector which contains the likelihood of each unit to
% the closest cluster

[disc,idx]= max(ll,[],2);

% Sort the n likelihood contributions
% qq contains the orderd (weighted) likelihood contributions
[~,qq]=sort(disc,'descend');

% %%%%%

% Find final trimmed and untrimmed units for final classification
if mixt>=1 % ==2
    
    % Sort the n likelihood contributions
    % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
    [~,qqmixt]=sort(logpdf,'descend');
    
    unassignedmixt=qqmixt(h+1:n);
    assignedmixt=qqmixt(1:h);
    
    % Store in vector idx the cluster associated to the highest posterior
    % probability
    [~,idxmixt]=max(postprob,[],2);
    idxmixt(unassignedmixt)=0;
    
    postprob(unassignedmixt,:)=0;
    % Remark:
    % If there was full convergence sum(logpdf(assigned)) = vopt
else
    
    unassigned=qq(h+1:n);
    % Assign observations to clusters and assign a 0 value to trimmed ones
    idx(unassigned)=0;
    postprob(unassigned,:)=0;
end


% Compute AIC and BIC

if mixt>=1
    % Compute value of the maximized MiXTURE log likelihood
    [NlogLmixt]=estepFS(ll(assignedmixt,:));
    
    % NlogLmixt is the negative of the maximized MIXTURE LOG-LIKELIHOOD
    % Note that if there was convergence NlogL should be exactly equal to
    % -vopt
    NlogLmixt = -NlogLmixt;
end

% Note that disc(qq(1:h)) is the contribution to the CLASSIFICATION
% loglikelihood of the untrimmed units
loglik=disc(qq(1:h));

% NlogL is the negative of the CLASSIFICATION LOG-LIKELIHOOD  of the
% untrimmed units
% NlogL=-sum(max(ll(untrimmed units,[],2));
% Note that if there was convergence NlogL should be exactly equal to
% -vopt
NlogL =-sum(loglik);

% Store robust estimate of final centroids of the groups
out.muopt=muopt;

% Store robust estimate of final covariance matrix of the groups
out.sigmaopt=sigmaopt;

% Store the assignments in matrix out. Unassigned units have an assignment
% equal to 0
if mixt>=1
    out.idx=idxmixt;
else
    out.idx=idx;
end

% siz = matrix of size k x 3,
% 1st col = sequence from 0 to k
% 2nd col = number of observations in each cluster
% 3rd col = percentage of observations in each cluster
% sum(out.siz(:,2))=n
% sum(out.siz(:,3))=100
siz=tabulate(out.idx);
out.siz=siz;

% Store n x k matrix containing posterior probability
% of each row from each component (cluster)
out.postprob=postprob;

% Number of estimated parameters
% k centroids of size v
% 0.5*v*(v+1) estimates for each of the k covariance matrices
npar=v*k;

% if equalweights = false the k-1 mixture proportions parameters must be added
if equalweights==false
    npar=npar +(k-1);
end

% Find number of constraints which are imposed on the covariance matrices
covunrestr=zeros(v,v,k);
% Find number of restricted eigenvalues for each group
constr=zeros(k,1);
for j=1:k
    
    selj=idx==j;
    nopt(j) = sum(selj);
    
    if nopt(j)>v
        covj=cov(Y(selj,:));
        covunrestr(:,:,j)=covj;
        % disp(['group' num2str(j) num2str(covj)])
        %         try
        [~,values] = eigs(covj);
        %         catch
        %             jj=100;
        %         end
        
        if v>1
            values=sort(diag(values));
            eigun=values./values(1);
            % maximum value of r is v-1
            r=sum(eigun(2:end)>=restrfactor(1));
            %
            constr(j)=r;
        else
            constr(j)=0;
        end
        
    else
    end
    
    %     eigunrestr((v*(j-1)+1):j*v)=diag(values);
    %     [~,values] = eigs(sigmaopt(:,:,j));
    %     eigrestr((v*(j-1)+1):j*v)=diag(values);
    %
end



%% Empirical quantities stored when there is no convergence

% unique ID found which are not outliers
UniqID = unique(idx(idx>0));

if length(UniqID) ~= k
    % Compute the empirical statistics when the algorithm does not reach
    % convergence (i.e. some cluster are missing)
    
    % initialize muemp and sigmaemp
    muemp    = nan(size(muopt));
    sigmaemp = nan(v,v,k);
    
    % iterate for each cluster found
    for j=1:length(UniqID)
        
        % assign the jj-th cluster ID
        jj = UniqID(j);
        
        if sum(idx==jj)>1
            % when more than one unit is in the jj-th cluster
            muemp(jj,:)      = mean(Y(idx==jj,:));
            sigmaemp(:,:,jj) = cov(Y(idx==jj,:));
        else
            % when one unit is in the jj-th cluster
            muemp(jj,:)      = Y(idx==jj,:);
            sigmaemp(:,:,jj) = 0;
        end
    end
    
    % restore apropriate order of the ID
    realID = 1:k;                           % searched clusters
    NanGroups = ~ismember(realID, UniqID);  % missing clusters
    % initialize idxemp
    idxemp = idx;
    % order the clusters found
    for i = 1:k-sum(NanGroups)
        posChang = idx==UniqID(i);
        UniqID(i) = UniqID(i) - sum(NanGroups(1:UniqID(i)));
        idxemp(posChang) = UniqID(i);
    end
    
    % put NaN at the end of muopt, sigmaopt, siz
    muemp = [muemp(~NanGroups,:); muemp(NanGroups,:)];
    sigmaemp = cat(3, sigmaemp(:,:,~NanGroups), sigmaemp(:,:,NanGroups));
    sizemp=tabulate(idxemp);
    misSiz = k-length(UniqID) ; % missing rows in siz
    sizemp = [sizemp; nan(misSiz, 3)];
    
    % Store empirical centroids, covariance matrices, mixing proportions
    % and ID
    emp = struct;
    emp.idxemp = idxemp;
    emp.muemp = muemp;
    emp.sigmaemp = sigmaemp;
    emp.sizemp = sizemp;
    
    % save the structure in the structure out
    out.emp = emp;
    
else
    % assign zero to the field when convergence is obtained
    emp = 0;
    out.emp = emp;
end


%% Compute INFORMATION CRITERIA

% nParam=npar+ 0.5*v*(v-1)*k + (v*k-1)*((1-1/restrfactor(1))^(1-1/(v*k))) +1;
nParam=npar+ 0.5*v*(v-1)*k + (v*k-1)*(1-1/restrfactor(1)) +1;

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


% Store the fraction of subsamples without convergence.
out.notconver=notconver;
alp=alpha>0;
if msg==1
    if size(siz,1)<k+alp
        disp(['Number of supplied clusters =' num2str(k)])
        disp(['Number of estimated clusters =' num2str(size(siz,1)-alp)])
        warning('FSDA:tclust:WrongKObtained','The total number of estimated clusters is smaller than the number supplied')
    end
end

if min(out.siz((1+alp):end,2)< n*0.02)
    warning('FSDA:tclust:TooSmallNj','Clusters with size < n * 0.02 found - try reducing k')
end

% Store units forming initial subset which gave rise to the optimal
% solution
out.bs=bs;

% Store value of the objective function (maximized trimmed log likelihood)
out.obj=vopt;

if out.obj==-1e+25
    warning('FSDA:tclust:NoConvergence','The result is artificially constrained due to restr.fact = 1')
end


out.equalweights=equalweights;

% Store the number of observations that have not been trimmed in the
% computation of the centroids
out.h=h;

out.fullsol=fullsol;

if options.Ysave
    % Store original data matrix
    out.Y=Y;
end


%% Create plots

% Plot the groups. Depending on v (univariate, bivariate, multivariate),
% we generate different plot types.
if  isstruct(plots) || (~iscell(plots) && isscalar(plots) && plots==1) || ... % integer equal to 1 or structure
        ((ischar(plots) || iscell(plots)) && max(strcmp(plots,{'contourf','contour','surf','mesh','ellipse','boxplotb'}))) || ... % char or cell of one of the specified string
        (iscell(plots) && isstruct(cell2mat(plots))) % cell containing a structure
    
    % The following statement is necessary because if mixt>0
    % idx was called idxmixt;
    idx=out.idx;
    
    % change the ID used if empirical values are evaluated
    if isstruct(emp)
        idx = idxemp;
    end
    
    if v==1
        
        % Univariate case: plot the histogram
        figure;
        histFS(Y,length(Y),idx);
        
    elseif v>=2
        
        % Bivariate plot, optionally with confidence ellipses, density
        % countours or bivariate boxplot
        
        % extract char or struct fro the cell if needed
        if iscell(plots)
            plots = cell2mat(plots);
        end
        
        % define what to superimpose on the plot
        if ischar(plots)
            overlay.type = plots;
        elseif isstruct(plots)
            overlay = plots;
        elseif plots==1
            % if plots=1 do not add anything to the scatter plot
            overlay ='';
        end
        
        % exclude outliers if present (when plots is char or struct)
        if any(idx<=0) && ~isempty(overlay)
            overlay.include = true(length(unique(idx)), 1);
            overlay.include(unique(idx)<=0) = false;
        end
        
        % differentiate for bivariate and multivariate data
        if v==2
            undock = [2 1];
        else
            undock = '';
        end
        
        % show axes label
        plo.labeladd=1;
        
        % use black color for outliers (i.e. k is the first color, used for group 0)
        if any(idx==0)
            plo.clr = 'kbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcy';
            plo.clr = plo.clr(1:length(unique(idx)));
        end
        
        % Add Labels Names
        % [To Enhance if used we lose the labels' order, it needs to
        % change on the spmplot call idx with id]
        % id=cellstr(num2str(idx));
        % id(idx==0)=cellstr('Trimmed units');
        
        % bivariate scatter
        figure;
        spmplot(Y, 'group', idx, 'plo', plo, 'undock', undock, 'overlay', overlay);
        
    end
    
    % add title
    str = sprintf('%d groups found by tclust for %s=%.2f and %s= %0.f', sum(unique(idx)>0),'$\alpha$',alpha,'$c$',restrfactor(1));
    title(str,'Interpreter','Latex'); % , 'fontsize', 14
    
elseif isscalar(plots) && plots == 0
    % does anything
else
    warning('FSDA:tclust:WrongInp','The parameter ''plots'' is not valid.');
end


end
%FScategory:CLUS-RobClaMULT