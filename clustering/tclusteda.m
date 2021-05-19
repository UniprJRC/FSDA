function [out,varargout]  = tclusteda(Y,k,alpha,restrfactor,varargin)
%tclusteda computes tclust for a series of values of the trimming factor
%
%<a href="matlab: docsearchFS('tclusteda')">Link to the help function</a>
%
%   tclusteda performs tclust for a series of values of the trimming factor
%   alpha given k (number of groups) and given c (restriction factor). In
%   order to increase the speed of the computations, parfor is used.
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
%        alpha: trimming level to monitor. Vector. Vector which specifies the
%               values of trimming levels which have to be considered.
%               alpha is a vector which contains decreasing elements which
%               lie in the interval 0 and 0.5.
%               For example if alpha=[0.1 0.05 0] tclusteda considers these 3
%               values of trimming level.
%               If alpha=0 tclusteda reduces to traditional model
%               based or mixture clustering (mclust): see Matlab function
%               gmdistribution. The default for alpha is vector [0.1 0.05
%               0]. The sequence is forced to be monotonically decreasing.
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
%
%    refsteps : Number of refining iterations. Scalar. Number of refining
%               iterations in each subsample. Default is 15.
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'refsteps',10
%                 Data Types - single | double
%
%     reftol  : Tolerance for the refining steps. Scalar.
%               The default value is 1e-14;
%                 Example - 'reftol',1e-05
%                 Data Types - single | double
%
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
% plots    :    Plot on the screen. Scalar structure.
%
%               Case 1: plots option used as scalar.
%               - If plots=0,  plots are not generated.
%               - If plots=1 (default), two plots are shown on the screen.
%                 The first plot ("monitor plot") shows three panels
%                 monitoring between two consecutive values of alpha the
%                 change in classification using ARI index (top panel), the
%                 change in centroids using squared euclidean distances
%                 (central panel), the change in covariance matrices using
%                 squared euclidean distance (bottom panel).
%                 The second plot ("gscatter plot") shows a series of
%                 subplots which monitor the classification for each value
%                 of alpha. In order to make sure that consistent labels
%                 are used for the groups, between two consecutive values
%                 of alpha, we assign label r to a group if this group
%                 shows the smallest distance with group r for the previous
%                 value of alpha. The type of plot which is used to monitor
%                 the stability of the classification depends on the value
%                 of v.
%                   * for v=1, we use histograms of the univariate data
%                   (function histFS is called).
%                   * for v=2, we use the scatter plot of the two
%                   variables (function gscatter is called).
%                   * for v>2, we use the scatter plot of the first two
%                   principal components (function gscatter is called and
%                   we show on the axes titles the percentage of variance
%                   explained by the first two principal components).
%
%               Case 2: plots option used as struct.
%                 If plots is a structure it may contain the following fields:
%                 plots.name = cell array of strings which enables to
%                   specify which plot to display. plots.name = {'gscatter'}
%                   produces a figure with a series of subplots which show the
%                   classification for each value of alpha. plots.name = {'monitor'}
%                   shows a figure with 3 panels which monitor between two
%                   consecutive values of alpha the change in classification
%                   using ARI index (top panel), the change in centroids
%                   using squared euclidean distances (central panel), the
%                   change in covariance matrices using squared euclidean
%                   distance (bottom panel). If this field is
%                   not specified plots.name={'gscatter' 'monitor'} and
%                   both figures will be shown.
%                 plots.alphasel = numeric vector which speciies for which
%                   values of alpha it is possible to see the classification.
%                   For example if plots.alphasel =[ 0.05 0.02], the
%                   classification will be shown just for alpha=0.05 and
%                   alpha=0.02; If this field is
%                   not specified plots.alphasel=alpha and therefore the
%                   classification is shown for each value of alpha.
%                 plots.ylimy = 2D array of size 3-by 2 which specifies the
%                   lower and upper limits for the monitoring plots. The
%                   first row refers the ARI index (top panel), the second
%                   row refers to the the change in centroids using squared
%                   euclidean distances (central panel), the third row is
%                   associated with the change in covariance matrices using
%                   squared euclidean distance (bottom panel).
%                   Example - 'plots', 1
%                   Data Types - single | double | struct
%
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
%
%      nocheck: Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix Y.
%               As default nocheck=0.
%                   Example - 'nocheck',1
%                   Data Types - single | double
%
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
%
% RandNumbForNini: Pre-extracted random numbers to initialize proportions.
%                Matrix. Matrix with size k-by-size(nsamp,1) containing the
%                random numbers which are used to initialize the
%                proportions of the groups. This option is effective just
%                if nsamp is a matrix which contains pre-extracted
%                subsamples. The purpose of this option is to enable to
%                user to replicate the results in case routine tclust is
%                called using a parfor instruction (as it happens for
%                example in routine tclustIC, where tclust is called through a
%                parfor for different values of the restriction factor).
%                The default value of RandNumbForNini is empty that is
%                random numbers from uniform are used.
%                   Example - 'RandNumbForNini',''
%                   Data Types - single | double
%
%     restrtype : type of restriction. Character. The type of restriction to
%               be applied on the cluster scatter
%               matrices. Valid values are 'eigen' (default), or 'deter'.
%               eigen implies restriction on the eigenvalues while deter
%               implies restriction on the determinants.
%                 Example - 'restrtype','deter'
%                 Data Types - char
%
%       cshape :    constraint to apply to each of the shape matrices.
%                   Scalar greater or equal than 1. This options only works
%                   is 'restrtype' is 'deter'.
%               When restrtype is deter the default value of the "shape"
%               constraint (as defined below) applied to each group is fixed to
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
%
%   UnitsSameGroup :  list of the units which must (whenever possible)
%                   have a particular label. Numeric vector.  For example if
%                   UnitsSameGroup=[20 26], means that group which contains
%                   unit 20 is always labelled with number 1. Similarly,
%                   the group which contains unit 26 is always labelled
%                   with number 2, (unless it is found that unit 26 already
%                   belongs to group 1). In general, group which contains
%                   unit UnitsSameGroup(r) where r=2, ...length(kk)-1 is
%                   labelled with number r (unless it is found that unit
%                   UnitsSameGroup(r) has already been assigned to groups
%                   1, 2, ..., r-1).
%                 Example - 'UnitsSameGroup',[20 34]
%                 Data Types -  integer vector
%
%      numpool:     The number of parallel sessions to open. Integer. If
%                   numpool is not defined, then it is set equal to the
%                   number of physical cores in the computer.
%                 Example - 'numpool',4
%                 Data Types -  integer vector
%
%      cleanpool:   Function name. Scalar {0,1}. Indicated if the open pool
%                   must be closed or not. It is useful to leave it open if
%                   there are subsequent parallel sessions to execute, so
%                   that to save the time required to open a new pool.
%                 Example - 'cleanpool',true
%                   Data Types - integer | logical
%
%
%  Output:
%
%         out:   structure which contains the following fields
%
%            out.IDX  = n-by-length(alpha) vector containing assignment of each unit to
%                       each of the k groups. Cluster names are integer
%                       numbers from 1 to k. 0 indicates trimmed
%                       observations. First column refers of out.IDX refers
%                       to alpha(1), second column of out.IDX refers to
%                       alpha(2), ..., last column refers to alpha(end).
%
%            out.MU  =  3D array of size k-by-v-by-length(alpha) containing
%                       the monitoring of the centroid for each value of
%                       alpha. out.MU(:,:,1), refers to alpha(1) ...,
%                       out.MU(:,:,end) refers to alpha(end). First row in
%                       each slice refers to group 1, second row refers to
%                       group 2 ...
%
%         out.SIGMA  =  cell of length length(alpha) containing in element
%                       j, with j=1, 2, ...,  length(alpha), the 3D array
%                       of size v-by-v-k containing the k (constrained)
%                       estimated covariance matrices associated with
%                       alpha(j).
%
%         out.Amon  =  Amon stands for alpha monitoring. Matrix of size
%                      (length(alpha)-1)-by-4 which contains for two
%                       consecutive values of alpha the monitoring of three
%                       quantities (change in classification, change in
%                       centroid location, change in covariance location).
%                       1st col = value of alpha.
%                       2nd col = ARI index.
%                       3rd col = squared Euclidean distance between
%                           centroids.
%                       4th col = squared Euclidean distance between
%                           covariance matrices.
%
%              out.Y  = Original data matrix Y.
%
%  Optional Output:
%
%           outcell : cell of length length(alpha) which contains in jth
%                     position the structure which comes out from procedure
%                     tclust applied to alpha(j), with j =1, 2, ...,
%                     length(alpha).
%
% More About:
%
%
% This procedure extends to tclust the so called monitoring
% approach. The phylosophy is to investigate how the results change as the
% trimming proportion alpha reduces. This function enables us to monitor
% the change in classification (measured by the ARI index) and the change
% in centroids and covariance matrices (measured by the squared euclidean
% distances). In order to make sure that consistent labels are used for the
% groups, between two consecutive values of alpha, we assign label r to a
% group if this group shows the smallest distance with group r for the
% previous value of alpha.
%
% See also: tclust, tclustIC
%
% References:
%
%   Garcia-Escudero, L.A., Gordaliza, A., Matran, C. and Mayo-Iscar, A. (2008),
%   A General Trimming Approach to Robust Cluster Analysis. Annals
%   of Statistics, Vol. 36, 1324-1345. [Technical Report available at:
%   www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf]
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('tclusteda')">Link to the help function</a>
%
%$LastChangedDate:: 2018-01-29 18:52:14 #$: Date of the last commit

% Examples:

%{
    %% Monitoring using geyser data (all default options).
    close all
    Y=load('geyser2.txt');
    % alpha and restriction factor are not specified therefore for alpha
    % vector [0.10 0.05 0] is used while for the restriction factor, value c=12
    % is used
    k=3;
    [out]=tclusteda(Y,k);
%}

%{
    % Monitoring using geyser data with alpha and c specified.
    Y=load('geyser2.txt');
    close all
    % alphavec= vector which contains the trimming levels to consider
    alphavec=0.10:-0.01:0;
    % c = restriction factor to use
    c=100;
    % k= number of groups
    k=3;
    [out]=tclusteda(Y,k,alphavec,c);
%}

%{
    %% Monitoring using geyser data with option plots supplied as structure.
    Y=load('geyser2.txt');
    close all
    % alphavec= vector which contains the trimming levels to consider
    % in this case 31 values of alpha are considered
    alphavec=0.30:-0.01:0;
    % c = restriction factor to use
    c=100;
    % k= number of groups
    k=3;
    % The monitoring plot of allocation will shows just these four values of
    % alpha
    plots=struct;
    plots.alphasel=[0.2 0.10 0.05 0.01];
    [out]=tclusteda(Y,k,alphavec,c,'plots',plots);
%}

%{
    %% Monitoring geyser data with option UnitsSameGroup.
    Y=load('geyser2.txt');
    close all
    % alphavec= vector which contains the trimming levels to consider
    alphavec=0.30:-0.10:0;
    % c = restriction factor to use
    c=100;
    % k= number of groups
    k=3;
    % Make sure that group containing unit 10 is in a group which is labelled
    % group 1 and group containing unit 12 is in group which is labelled group 2
    UnitsSameGroup=[10 12];
    % Mixture model is used
    mixt=2;
    [out]=tclusteda(Y,k,alphavec,1000,'mixt',2,'UnitsSameGroup',UnitsSameGroup);
%}

%{
    %% tclusteda with M5 data.
    close all
    Y=load('M5data.txt');

    % alphavec= vector which contains the trimming levels to consider
    alphavec=0.10:-0.02:0;
    out=tclusteda(Y(:,1:2),3,alphavec,1000,'nsamp',1000,'plots',1);
%}

%{
    % Structured noise data ex1.
    close all
    Y=load('structurednoise.txt');
    alphavec=0.20:-0.01:0;
    out=tclusteda(Y,2,alphavec,100,'plots',1);
%}

%{
    % Structured noise data ex2.
    close all
    Y=load('structurednoise.txt');
    alphavec=0.20:-0.01:0;
    % just show the monitoring plot
    plots=struct;
    plots.name = {'monitor'};
    out=tclusteda(Y,2,alphavec,100,'plots',plots);
%}

%{
    % mixture100 data.
    close all
    Y=load('mixture100.txt');
    % Traditional tclust
    alphavec=0.20:-0.01:0;
    % just show the allocation plot
    plots=struct;
    plots.name = {'gscatter'};
    out=tclusteda(Y,2,alphavec,100,'plots',plots);
%}

%{
    %% tclusteda using simulated data.
    % 5 groups and 5 variables
    rng(100,'twister')
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
    n=size(Y,1);
    % Set number of groups
    k=5;

    % Example of the subsets precalculated
    nsamp=2000;
    nsampscalar=nsamp;
    nsamp=subsets(nsamp,n,(v+1)*k);
    % Random numbers to compute proportions computed once and for all
    RandNumbForNini=rand(k,nsampscalar);
    % The allocation is shown on the space of the first two principal
    % components
    out=tclusteda(Y,k,[],6,'plots',1,'RandNumbForNini',RandNumbForNini,'nsamp',nsamp);
%}

%{
    % tclusteda using determinant constraint.
    % Search for spherical clusters.
    % 5 groups and 5 variables
    rng(100,'twister')
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
    n=size(Y,1);
    % Set number of groups
    k=5;
    cshape=1
    out=tclusteda(Y,k,[],1000,'plots',1,'restrtype','deter','cshape',cshape);
%}

%{
    % An example of use of plots as a structure with field ylimy.
    load('swiss_banknotes');
    Y=swiss_banknotes{:,:};
    [n,v]=size(Y);
    alphavec=0.15:-0.01:0;
    % alphavec=0.12:-0.005:0;

    % c = restriction factor to use
    c=100;
    % k= number of groups
    k=2;
    % restriction on the determinants is imposed
    restrtype='deter';
    % Specify lower and upper limits for the monitoring plot
    plots=struct;
    % ylimits for monitoring of ARI index
    ylimARI=[0.95 1];
    % ylimits for change in centroids
    ylimCENT=[0 0.02];
    % ylimits for change in cov matrices
    ylimCOV=[0 0.01];
    ylimy=[ylimARI;ylimCENT;ylimCOV];
    plots.ylimy=ylimy;
    [outDet]=tclusteda(Y,k,alphavec,c,'restrtype',restrtype,'plots',plots,'nsamp',10000);
%}

%% Beginning of code

% Input parameters checking
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
                error('FSDA:tclusteda:WrongNsamp','Wrong number of columns in matrix nsamp')
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
    
    ncomb=bc(n,k*(v+1));
    nsamp=min(300,ncomb);
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
    nsampdef=min(2000,ncomb);
end

refstepsdef=15;
reftoldef=1e-5;

% tolrestreigen = tolerance to use in function restreigen
tolrestreigen=1e-08;

% Default
if nargin<3 || isempty(alpha)
    alpha=[0.10 0.05 0];
    warning('FSDA:tclusteda:Wrongalpha','You have not specified alpha: it is set to [0.10 0.05 0] by default');
end

% Fix alpha equal to the trimming size
% h = number of observations which is used to compute the centroids

if min(alpha)<0
    error('FSDA:tclusteda:WrongAlpha','alpha must a vector with numbers in the interval [0 0.5]')
end

if nargin<4 || isempty(restrfactor)
    restrfactor=12;
    warning('FSDA:tclusteda:Wrongrestrfact','You have not specified restrfactor: it is set to 12 by default');
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
RandNumbForNini='';
% cshape. Constraint on the shape matrices inside each group which works only if restrtype is 'deter'
cshape=10^10;

UserOptions=varargin(1:2:length(varargin));

if ~isempty(UserOptions)
    
    options=struct('nsamp',nsampdef,'RandNumbForNini','','plots',1,'nocheck',0,...
        'msg',1,'refsteps',refstepsdef,'equalweights',false,...
        'reftol',reftoldef,'mixt',0,'startv1',startv1def,'restrtype','eigen',...
        'UnitsSameGroup',UnitsSameGroup,...
        'numpool',numpool, 'cleanpool', cleanpool,'cshape',cshape);
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:tclusteda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDAeda:tclusteda:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
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
        error('FSDAeda:tclusteda:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
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
        error('FSDAeda:tclusteda:Wrongrestr','Wrong restriction');
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
    mixt=options.mixt;         % if options.mixt>0 mixture model is assumed
else
    mixt=0;
    msg=0;
    equalweights=false;
    reftol=reftoldef;
    refsteps=refstepsdef;
    plots=1;
end


if isempty(RandNumbForNini)
    NoPriorNini=1;
else
    NoPriorNini=0;
end

if mixt>=1 && equalweights == true
    warning('FSDA:tclusteda:WrongEqualWeights','option equalweights must be different from 1 if mixture model approach is assumed')
    warning('FSDA:tclusteda:WrongEqualWeights','options equalweights is reset to 0')
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

% verLess2016b is true if current version is smaller than 2016b
verLess2016b=verLessThanFS(9.1);

if verLess2016b ==1
    userepmat=1;
else
    userepmat=2;
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

%%  Loop through all subsets of nselected and store relevant quantitites
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
            autovalues=restrdeter(Lambda_vk,niini,restrfactor ,tolrestreigen,userepmat);
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

nnargout=nargout;

% Number of estimated parameters
% k centroids of size v
% 0.5*v*(v+1) estimates for each of the k covariance matrices
npar=v*k;

% if equalweights = false the k-1 mixture proportions parameters must be added
if equalweights==false
    npar=npar +(k-1);
end
nParam=npar+ 0.5*v*(v-1)*k + (v*k-1)*((1-1/restrfactor(1))^(1-1/(v*k))) +1;

lalpha=length(alpha);
if msg == 1
    progbar = ProgressBar(lalpha);
else
    progbar=[];
end

IDX=zeros(n,lalpha);
outcell=cell(lalpha,1);

MU=zeros(k,v,lalpha);
SIGMA=cell(lalpha,1);

parfor (j=1:lalpha, numpool)
    outj  = tclustcore(Y,Cini,Sigmaini,Niini,reftol,refsteps,mixt, ...
        equalweights,hh(j),nselected,k,restrnum,restrfactor,userepmat,nParam);
    
    if nnargout==2
        outcell{j}=outj;
    end
    
    IDX(:,j)=outj.idx;
    
    MU(:,:,j)=outj.muopt;
    SIGMA{j}=outj.sigmaopt;
    
    if msg == 1
        progbar.progress;  %#ok<PFBNS>
    end
end

if msg == 1
    progbar.stop;
end


if ~isempty(UnitsSameGroup)
    
    [IDXnew1, OldNewIndexes]=ClusterRelabel({IDX(:,1)}, UnitsSameGroup);
    
    MUold1=MU(:,:,1);
    SIGMAold1= SIGMA{1};
    
    MUnew1=MUold1;
    SIGMAnew1=SIGMAold1;
    for jj=1:size(OldNewIndexes,1)
        MUnew1(OldNewIndexes(jj,1),:)= MUold1(OldNewIndexes(jj,2),:);
        MUnew1(OldNewIndexes(jj,2),:)= MUold1(OldNewIndexes(jj,1),:);
        MUold1=MUnew1;
        
        SIGMAnew1(:,:,OldNewIndexes(jj,1))=SIGMAold1(:,:,OldNewIndexes(jj,2));
        SIGMAnew1(:,:,OldNewIndexes(jj,2))=SIGMAold1(:,:,OldNewIndexes(jj,1));
        SIGMAold1=SIGMAnew1;
    end
    IDX(:,1)=IDXnew1{:};
    MU(:,:,1)=MUold1;
    SIGMA{1}=SIGMAold1;
end


%% Monitor the difference in classification, centroids and covariance matrices

% Amon stands for alpha monitoring.
% Amon is the matrix of size lenght(alpha)-1-by 4 which contains for two
% consecutive values of alpha the monitoring of three quantities.
% 1st col = value of alpha
% 2nd col = ARI index
% 3rd col = squared Euclidean distance between consecutive centroids
% 4th col = squared Euclidean distance between consecutive covariance matrices
Amon=[alpha(2:end)' zeros(lalpha-1,3)];

noisecluster=0;

IDXold=IDX;

maxdist=zeros(lalpha,1);
seqk=(1:k)';

%verMatlab=verLessThanFS('9.2');

for j=2:lalpha
    newlab=zeros(k,1);
    mindist=newlab;
    for ii=1:k
        % centroid of group ii for previous alpha value
        muii=MU(ii,:,j-1);
        % MU(:,:,j) =matrix of centroids for current alpha value
        
        %if verMatlab==true;
        muij=bsxfun(@minus,muii,MU(:,:,j));
        %else
        %    muij=muii-MU(:,:,j);
        %end
        
        [mind,indmin]=min(sum(muij.^2,2));
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
        % In this case new1 contains the labels which never appeared inside
        % newlab. To this label we assign the maximum distance and check is
        % this time sequal(sort(newlab),seqk), that is we check whether
        % vector sort(newlab) of length k contain the numbers 1, 2, ..., k
        % if length(newl) two labels do not have the correspondence
        % therefore automatic relabelling is not possible.
        newl=setdiff(seqk,newlab);
        notinseqk=setdiff(seqk,newl);
        if length(newl)==1 && length(notinseqk)==1
            newlab(indmaxdist)=notinseqk;
            if isequal(sort(newlab),seqk)
                MU(:,:,j)=MU(newlab,:,j);
                SIGMA(j)= {SIGMA{j}(:,:,newlab)};
                for r=1:k
                    IDX(IDXold(:,j)==newlab(r),j)=r;
                end
            else
                disp(['Automatic relabelling not possible when alpha=' num2str(alpha(j))])
            end
        else
            disp(['Automatic relabelling not possible when alpha=' num2str(alpha(j))])
        end
    end
end

for j=2:lalpha
    
    % Compute ARI index between two consecutive alpha values
    [ARI]=RandIndexFS(IDX(:,j-1),IDX(:,j),noisecluster);
    % Store in the second column the ARI index
    Amon(j-1,2)=ARI;
    
    % Compute and store squared euclidean distance between consecutive
    % centroids
    Amon(j-1,3)=sum(sum( (MU(:,:,j)-MU(:,:,j-1)).^2, 2));
    
    % Compute and store squared euclidean distance between consecutive
    % covariance matrices (all elements of cov matrices are considered)
    dxdiag=0;
    for i=1:k
        dxdiag=dxdiag+sum((diag(SIGMA{j}(:,:,i))-diag(SIGMA{j-1}(:,:,i))).^2);
    end
    Amon(j-1,4)=sum(sum(sum((SIGMA{j}-SIGMA{j-1}).^2,2)));
    % sumdistCOVonlydiag(j)=dxdiag;
end

out=struct;

% Store classification
out.IDX=IDX;
% Store centroids
out.MU=MU;
% Store covariance matrices
out.SIGMA=SIGMA;
% Store ARI index, variation in centroid location and
% variation in covariance.
out.Amon=Amon;
% Store Y
out.Y=Y;

% Store the indices in varargout
if nnargout==2
    varargout=outcell;
end



%% Plotting part
if isstruct(plots)
    fplots=fieldnames(plots);
    
    d=find(strcmp('name',fplots));
    if d>0
        name=plots.name;
        if ~iscell(name)
            error('FSDA:tclusteda:Wronginput','plots.name must be a cell')
        end
    else
        name={'gscatter' 'monitor'};
    end
    
    d=find(strcmp('alphasel',fplots));
    if d>0
        alphasel=plots.alphasel;
    else
        alphasel=alpha;
    end
    
    
    d=find(strcmp('ylimy',fplots));
    if d>0
        ylimy=plots.ylimy;
        [nylim,vylim]=size(ylimy);
        if nylim~=3
            error('FSDA:tclusteda:Wronginput','plots.ylimy must be a matrix with 3 rows')
        end
        if vylim~=2
            error('FSDA:tclusteda:Wronginput','plots.ylimy must be a matrix with 2 columns')
        end
    else
        name={'gscatter' 'monitor'};
        alphasel=alpha;
        ylimy='';
    end
elseif plots==1
    name={'gscatter' 'monitor'};
    alphasel=alpha;
    ylimy='';
else
    name = '';
end

d=find(strcmp('monitor',name));

if d>0
    % ARI between two consecutive values of alpha
    subplot(3,1,1)
    plot(Amon(:,1),Amon(:,2))
    set(gca,'XDir','reverse');
    xlabel('\alpha')
    ylabel('ARI index')
    set(gca,'FontSize',16)
    if ~isempty(ylimy)
        ylim(ylimy(1,:))
    end
    
    % Monitoring of centroid changes
    subplot(3,1,2)
    plot(Amon(:,1),Amon(:,3))
    set(gca,'XDir','reverse');
    xlabel('\alpha')
    ylabel('Centroids')
    set(gca,'FontSize',16)
    if ~isempty(ylimy)
        ylim(ylimy(2,:))
    end
    
    % Monitoring of covariance matrices change
    subplot(3,1,3)
    plot(Amon(:,1),Amon(:,3))
    set(gca,'XDir','reverse');
    xlabel('\alpha')
    ylabel('Covariance')
    set(gca,'FontSize',16)
    if ~isempty(ylimy)
        ylim(ylimy(3,:))
    end
    
end

d=find(strcmp('gscatter',name));
if d>0
    
    % alphasel contains the indexes of the columns of matrix IDX which have
    % to be plotted
    
    % We use round(alpha*1e+7)/1e+7 to guarantee compatibility with old
    % versions of MATLAB. For the new versions the instruction would have
    % been:
    % [~,alphasel]=intersect(round(alpha,9),alphasel,'stable');
    [~,alphasel]=intersect(round(alpha*1e+7)/1e+7,round(alphasel*1e+7)/1e+7,'stable');
    lalphasel=length(alphasel);
    
    %% Monitoring of allocation
    if  lalphasel==1
        nr=1;
        nc=1;
    elseif lalphasel==2
        nr=2;
        nc=1;
    elseif lalphasel<=4
        nr=2;
        nc=2;
    elseif lalphasel<=6
        nr=3;
        nc=2;
    elseif lalphasel<=9
        nr=3;
        nc=3;
    elseif lalphasel<=12
        nr=3;
        nc=4;
    else
        nr=4;
        nc=4;
    end
    
    resup=1;
    figure('Name',['Monitoring allocation #' int2str(resup)])
    
    colord='brkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcy';
    symdef={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h'};
    
    % Plot first two principal components in presence of more than two
    % variables
    if v>2
        Yst=zscore(Y);
        [V,D]=eig(cov(Yst));
        [Dsort,ind] = sort(diag(D),'descend');
        Ypca=Yst*V(:,ind(1:2));
        explained=100*Dsort(1:2)/sum(Dsort);
        % Note that the rows above are just for retrocompatibility
        % Those who have a release >=2012B can use
        % [~,Ypca,~,~,explained]=pca(zscore(Y),'NumComponents',2);
    else
        Ypca=Y;
    end
    
    jk=1;
    for j=1:lalphasel
        
        % The monitoring must contain a maximum of 16 panels
        % If length(alpha) is greater than 16 a new set of 16 subpanels is
        if jk>16
            jk=1;
            resup=resup+1;
            figure('Name',['Monitoring allocation #' int2str(resup)])
            subplot(nr,nc,jk)
        else
            subplot(nr,nc,jk)
        end
        jk=jk+1;
        
        
        if v>=2
            if alpha(alphasel(j))~=0
                hh=gscatter(Ypca(:,1),Ypca(:,2),IDX(:,alphasel(j)),colord,[symdef{1:k+1}]);
            else
                hh=gscatter(Ypca(:,1),Ypca(:,2),IDX(:,alphasel(j)),colord(2:k+1),[symdef{2:k+1}]);
            end
            
            if v>2
                xlabel(['PCA1 - ' num2str(explained(1)) '%'])
                ylabel(['PCA2 - ' num2str(explained(2)) '%'])
            else
                xlabel('y1')
                ylabel('y2')
            end
            
            clickableMultiLegend(hh)
            if jk>2
                legend hide
            end
            axis manual
        else
            % Univariate case: plot the histogram
            if alpha(alphasel(j))~=0
                histFS(Y,10,IDX(:,alphasel(j)),[],[],colord)
            else
                histFS(Y,10,IDX(:,alphasel(j)),[],[],colord(2:k+1))
            end
        end
        title(['$\alpha=$' num2str(alpha(alphasel(j)))],'Interpreter','Latex')
    end
end

% Clear all persistent variables inside tclustcore
clear tclustcore

end
%FScategory:CLUS-RobClaMULT
