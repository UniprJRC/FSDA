function [OUTcell, out]  = tclusteda(Y,k,alpha,restrfactor,varargin)
%tclust computes trimmed clustering with scatter restrictions
%
%<a href="matlab: docsearchFS('tclust')">Link to the help function</a>
%
%   tclusteda performs tclust for a series of values of the trimming factor
%   alpha given k (number of groups) and given c (restriction factor). In
%   order to increase the speed of the computations, parellal for is used.
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
%        alpha: trimming level to monitor. vector. vector which specifies the
%               values of trimming levels which have to be considered.
%               alpha is a scalar between 0 and 0.5.
%               For example is alpha=[0.1 0.05 0] tclust considers these 3
%               values of trmming level.
%               If alpha=0 tclust reduces to traditional model
%               based or mixture clustering (mclust): see Matlab function
%               gmdistribution. The default value of alpha is [0.1 0.05 0]
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
% plots    :    Plot on the screen. Scalar, character, cell or struct.
%               - If plots=0 (default), plots are not generated.
%               - If plot=1, a plot with the classification is shown on
%                 the screen (using the spmplot function). The plot can be:
%                   * for v=1, an histogram of the univariate data.
%                   * for v=2, a bivariate scatterplot.
%                   * for v>2, a scatterplot matrix generated by spmplot.
%               When v>=2 plots offers the following additional features
%               (for v=1 the behaviour is forced to be as for plots=1):
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
%               REMARK - The labels=0 are automatically excluded from the
%                          overlaying phase, considering them as outliers.
%                   Example - 'plots', 1
%                   Data Types - single | double | string
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
%               implies restriction on the determinants.
%                 Example - 'restrtype','deter'
%                 Data Types - char
%       Ysave : Save original input matrix. Scalar. Set Ysave to 1 to request that the input matrix Y
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
%            out.muopt= k-by-v-length(alpha) array containing cluster centroid
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
%
% This iterative algorithm initializes k clusters randomly and
% performs "concentration steps" in order to improve the current cluster
% assignment. The number of maximum concentration steps to be performed is
% given by input parameter refsteps. For approximately obtaining the global
% optimum, the system is initialized nsamp times and concentration steps
% are performed until convergence or refsteps is reached. When processing
% more complex data sets higher values of nsamp and refsteps have to be
% specified (obviously implying extra computation time). However, if more
% then 10 per cent of the iterations do not converge, a warning message is
% issued, indicating that nsamp has to be increased.
%
% See also: tkmeans, estepFS
%
% References:
%
% Garcia-Escudero, L.A., Gordaliza, A., Matran, C. and Mayo-Iscar, A. (2008),
% A General Trimming Approach to Robust Cluster Analysis. Annals
% of Statistics, Vol.36, 1324-1345.
% Technical Report available at:
% http://www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf
%
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('tclust')">Link to the help function</a>
%
%$LastChangedDate:: 2018-01-29 18:52:14 #$: Date of the last commit

% Examples:

%{
    % tclust of geyser data using k=3, alpha=0.1 and  restrfactor=10000.
    Y=load('geyser2.txt');
    out=tclusteda(Y,3,[0.1 0.05 0],1000);
%}

%{
    %% tclust of geyser with classification plot.
    Y=load('geyser2.txt');
    close all
    out=tclust(Y,3,0.1,10000,'plots',1);

    % default confidence ellipses.
    out=tclust(Y,3,0.1,10000,'plots','ellipse');

    % confidence ellipses specified by the user
    plots.type = 'ellipse';
    plots.conflev = 0.5;
    out=tclust(Y,3,0.1,10000,'plots',plots);

    % contour plots.
    out=tclust(Y,3,0.1,10000,'plots','contour');

    % filled contour plots with additional options
    plots.type = 'contourf';
    plots.cmap = autumn;
    out=tclust(Y,3,0.1,10000,'plots',plots);

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
    title(['Restriction factor =' num2str(restrfactor)])
    restrfactor=10;
    out=tclust(Y,3,0.1,restrfactor,'nsamp',nsamp,'refsteps',10,'plots','ellipse');
    title(['Restriction factor =' num2str(restrfactor)])
    % trimmed k-means solution restrfactor=1
    restrfactor=1;
    out=tclust(Y,3,0.1,restrfactor,'nsamp',nsamp,'refsteps',10,'plots','ellipse');
    title(['Restriction factor =' num2str(restrfactor) '. Trimmed k-means solution'])
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
    title('Traditional tclust');
    % tclust with mixture models (selection of untrimmed units according to
    % likelihood contributions
    out2=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1,'mixt',1)
    title('tclust with mixture models (likelihood contributions)');
    % Tclust with mixture models (selection of untrimmed units according to
    % densities weighted by estimates of the probability of the components)
    out3=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1,'mixt',2)
    title('tclust with mixture models (probability of the components)');
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
    warning('FSDA:tclust:Wrongalpha','You have not specified alpha: it is set to [0.10 0.05 0] by default');
else
    if isempty(alpha)
        alpha=[0.10 0.05 0];
        warning('FSDA:tclust:Wrongalpha','You have not specified alpha: it is set to [0.10 0.05 0] by default');
    end
end

% Fix alpha equal to the trimming size
% h = number of observations which is used to compute the centroids

if min(alpha)<0
    error('FSDA:tclust:WrongAlpha','alpha must a vector with numbers in the interval [0 0.5]')
end

if nargin<4
    restrfactor=12;
    warning('FSDA:tclust:Wrongrestrfact','You have not specified restrfactor: it is set to 12 by default');
end

% hh= vector containing number of untrimmed units for each value of alpha
% h = number of untrimmed units
hh=fix(n*(1-alpha));

% restrnum=1 implies eigenvalue restriction
restrnum=1;

% check how many physical cores are available in the computer (warning:
% function 'feature' is undocumented; however, FSDA is automatically
% monitored for errors and other inconsistencies at each new MATLAB
% release).
numpool = feature('numCores');
cleanpool=false;
UnitsSameGroup='';

UserOptions=varargin(1:2:length(varargin));

if ~isempty(UserOptions)
    
    options=struct('nsamp',nsampdef,'RandNumbForNini','','plots',1,'nocheck',0,...
        'msg',1,'Ysave',0,'refsteps',refstepsdef,'equalweights',false,...
        'reftol',reftoldef,'mixt',0,'startv1',startv1def,'restrtype','eigen',...
        'UnitsSameGroup',UnitsSameGroup,...
        'numpool',numpool, 'cleanpool', cleanpool);
    
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
    else
        error('FSDA:tclust:Wrongrestr','Wrong restriction');
    end
    
    
    % Default values for the optional
    % parameters are set inside structure 'options'
    
    plots=options.plots;        % Plot of the resulting classification
    nsamp=options.nsamp;        % Number of subsets to extract
    equalweights=options.equalweights;    % Specify if assignment must take into account the size of the groups
    
    refsteps=options.refsteps;
    reftol=options.reftol;
    
    numpool=options.numpool;
    UnitsSameGroup=options.UnitsSameGroup;
    RandNumbForNini=options.RandNumbForNini;
    
    msg=options.msg;            % Scalar which controls the messages displayed on the screen
    
    mixt=options.mixt;         % if options.mixt==1 mixture model is assumed
    
end


if isempty(RandNumbForNini)
    NoPriorNini=1;
else
    NoPriorNini=0;
end



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
        niinistart=repmat(floor(n/k),k,1);
    end
end

% Store the indices in varargout
if nargout==2
    varargout={C};
end




% The covariances are given initially by k identity matrices
ey=eye(v,v);
eyk=repmat(ey,[1 1 k]);


% Lambda_vk = matrix which will contain in column j the v (unrestricted)
% eigevalues of covariance matrix of group j (j=1, ..., k)
Lambda_vk=ones(v,k);


% Lambda will contain the matrix of eigenvalues in each iteration for
% all groups. Lambda is a 3D array of size v-by-v-by-k
% Lambda=sigmaini;
% U will contain the eigenvectors of the cov matrices in each iteration
% for all groups (or the shape matrices if restrnum=2). U is a 3D array of size v-by-v-by-k
sigmaini=zeros(v,v,k);
U=sigmaini;

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


Cini=cell(nselected,1);
Sigmaini=Cini;
Niini=Cini;

%%  Random starts
for i=1:nselected
    
    
    if startv1 ==1
        if NoPriorNini==1
            randk=rand(k,1);
        else
            randk=RandNumbForNini(:,i);
        end
        
        niini=floor(n*randk/sum(randk));
        
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
    Cini{i}=cini;
    Sigmaini{i}=sigmaini;
    Niini{i}=niini;
end

% Number of estimated parameters
% k centroids of size v
% 0.5*v*(v+1) estimates for each of the k covariance matrices
npar=v*k;

% if equalweights = false the k-1 mixture proportions parameters must be added
if equalweights==false
    npar=npar +(k-1);
end
nParam=npar+ 0.5*v*(v-1)*k + (v*k-1)*((1-1/restrfactor)^(1-1/(v*k))) +1;

progbar = ProgressBar(length(hh));

OUTcell=cell(length(hh),1);
parfor (jjj=1:length(hh), numpool)
    out  = tclustcore(Y,Cini,Sigmaini,Niini,reftol,refsteps,mixt, ...
        equalweights,hh(jjj),nselected,k,restrnum,restrfactor,userepmat,nParam);
    OUTcell{jjj}=out;
    if msg == 1
        progbar.progress;  %#ok<PFBNS>
    end
end

%% Monitor the difference in centroids and covariance matrices


noisecluster=0;
c1=OUTcell{1}.idx;
IDX=zeros(n,length(alpha));
IDX(:,1)=c1;


ARmon=[alpha(2:end)' zeros(length(alpha)-1,3)];

for j=2:length(alpha)
    
    c2=OUTcell{j}.idx;
    IDX(:,j)=c2;
    
    [ARI]=RandIndexFS(c1,c2,noisecluster);
    
    ARmon(j-1,2)=ARI;
    c1=c2;
    % disp(j)
end

IDXold=IDX;

MU=zeros(k,v,length(alpha));
SIGMA=cell(length(alpha),1);
for j=1:length(alpha)
    MU(:,:,j)=OUTcell{j}.muopt;
    SIGMA(j)={OUTcell{j}.sigmaopt};
end
maxdist=zeros(length(alpha),1);
seqk=(1:k)';
for j=2:length(alpha)
    newlab=zeros(k,1);
    mindist=newlab;
    for ii=1:k
        % centroid of group ii for previous alpha value
        muii=MU(ii,:,j-1);
        % MU(:,:,j) =matrix of centroids for current alpha value
        [mind,indmin]=min(sum((muii-MU(:,:,j)).^2,2));
        newlab(ii)=indmin;
        mindist(ii)=mind;
    end
    % Store maximum among minimum distances
    [maxmindist,indmaxdist] =max(mindist);
    maxdist(j)=maxmindist;
    
    if isequal(sort(newlab),seqk)
        MU(:,:,j)=MU(newlab,:,j);
        SIGMA(j)= {SIGMA{j}(:,:,newlab)};
        for r=1:k
            IDX(IDXold(:,j)==newlab(r),j)=r;
        end
    else
        newlab(indmaxdist)=setdiff(sqk,newlab);
        disp('Preliminary relabelling not possible')
        if isequal(sort(newlab),seqk)
            MU(:,:,j)=MU(newlab,:,j);
            SIGMA(j)= {SIGMA{j}(:,:,newlab)};
            for r=1:k
                IDX(IDXold(:,j)==newlab(r),j)=r;
            end
        else
            disp('Relabelling not possible')
        end
    end
end

out=struct;
out.IDX=IDX;
out.MU=MU;
out.SIGMA=SIGMA;


for j=2:length(alpha)
    ARmon(j-1,3)=sum(sum( (MU(:,:,j)-MU(:,:,j-1)).^2, 2));
    dxdiag=0;
    for i=1:k
        dxdiag=dxdiag+sum((diag(SIGMA{j}(:,:,i))-diag(SIGMA{j-1}(:,:,i))).^2);
    end
    ARmon(j-1,4)=sum(sum(sum((SIGMA{j}-SIGMA{j-1}).^2,2)));
    % sumdistCOVonlydiag(j)=dxdiag;
end


if plots==1
    
    % ARI between two consecutive values of alpha
    subplot(3,1,1)
    plot(ARmon(:,1),ARmon(:,2))
    set(gca,'XDir','reverse');
    xlabel('\alpha')
    ylabel('ARI index')
    set(gca,'FontSize',16)
    
    % Monitoring of centroid changes
    subplot(3,1,2)
    plot(ARmon(:,1),ARmon(:,3))
    set(gca,'XDir','reverse');
    xlabel('\alpha')
    ylabel('Centroids')
    set(gca,'FontSize',16)
    
    % Monitoring of covariance matrices change
    subplot(3,1,3)
    plot(ARmon(:,1),ARmon(:,3))
    set(gca,'XDir','reverse');
    xlabel('\alpha')
    ylabel('Covariance')
    set(gca,'FontSize',16)
    
    
    %% Monitoring of allocation as alpha varies
    figure
    colord='brkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcy';
    symdef={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h'};
    
    %   legend = @() clickableMultiLegend();
    
    if v>2
        [~,Ypca,~,~,explained]=pca(zscore(Y),'NumComponents',2);
    else
        Ypca=Y;
    end
    for j=1:length(alpha)
        % IDX(:,j)=out{j}.idx;
        subplot(4,3,j)
        if alpha(j)~=0
            hh=gscatter(Ypca(:,1),Ypca(:,2),IDX(:,j),colord,[symdef{1:k+1}]);
        else
            hh=gscatter(Ypca(:,1),Ypca(:,2),IDX(:,j),colord(2:k+1),[symdef{2:k+1}]);
        end
        if v>2
            xlabel(['PCA1 - ' num2str(explained(1)) '%'])
            ylabel(['PCA2 - ' num2str(explained(2)) '%'])
        else
            xlabel('y1')
            ylabel('y2')
        end
        clickableMultiLegend(hh)
        title(['$\alpha=$' num2str(alpha(j))],'Interpreter','Latex')
    end
end

end
%FScategory:CLUS-RobClaMULT