function out  = tclustICgpcm(Y, varargin)
%tclustICgpcm computes tclust for different number of groups k and restr. factors $c_{det}$ and $c_{shw}$
%
%<a href="matlab: docsearchFS('tclustICgpcm')">Link to the help function</a>
%
%   tclustICgpcm (where the two letters IC  stand for 'Information
%   Criterion') and gpcm stands for Gaussian parimonious clustering models
%   computes the values of BIC (MIXMIX), ICL (MIXCLA) or CLA (CLACLA), for
%   different values of k (number of groups) and different values of
%   $c_{det}$ (restriction factor for determinants), $c_{shw}$ (restriction
%   factor for within group shape element) and for a prespecified level of
%   trimming. If Parallel Computing toolbox is installed, parfor is used to
%   compute tclust for different values of restriction factors. In order to
%   minimize randomness, given k, the same subsets are used for each value
%   of the restriction factors.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. Data matrix containing n observations on v variables
%               Rows of Y represent observations, and columns represent
%               variables. Observations (rows) with missing (NaN) or or
%               infinite (Inf) values will automatically be excluded from
%               the computations.
%                 Data Types -  double
%
%  Optional input arguments:
%
%      pa  : Constraints to apply and model specification. Structure.
%            Structure containing the following fields:
%             pa.pars= type of Gaussian Parsimonious Clustering Model. Character.
%               A 3 letter word in the set:
%               'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',
%               'EVI','VEI','EEI','VII','EII'.
%               If pa is not supplied VVV model is used.
%             pa.cdet = scalar o vector containing the values of the
%               restriction factors for ratio of determinants which have to
%               be tested. If pa.cdet=1 all determinants are forced to be
%               equal. If this field is not present the values [1, 2, 4, 8,
%               16, 32, 64, 128] are used.
%               See section More About for additional details.
%             pa.shw = scalar o vector containing the values of the
%               restriction factors to the elements of
%               the shape matrices inside each group. If pa.shw=1 all diagonal
%               elements of the shape matrix of cluster j (with j=1, ...,
%               k) will be equal. If this field is not present the values [1, 2, 4, 8,
%               16, 32, 64, 128] are used.
%             pa.shb = scalar in the interval [1 Inf) which specifies the
%               the restriction which has to be applied to the elements of
%               the shape matrices across each group. If this field is not
%               present pa.shb=128 is used.
%             pa.maxiterS = positive integer which specifies the maximum
%               number of iterations to obtain the restricted shape matrix.
%               This parameter is used by routine restrshapepars. The
%               default value of pa.maxiterS is 5.
%             pa.maxiterR = positive integer which specifies the maximum
%               number of iterations to obtain the common rotation matrix
%               in presence of varying shape.
%               This parameter is used by routine cpcV. The
%               default value of pa.maxiterR is 20.
%          pa.maxiterDSR = positive integer which specifies the maximum
%               number of iterations to obtain the requested restricted
%               determinants, shape matrices and rotation. For all
%               parametrizations  pa.maxiterDSR is set to 1 apart from for
%               the specifications 'VVE', 'EVE' and 'VEE'. The default
%               value of pa.maxiterDSR is 20.
%           pa.tolS=tolerance to use to exit the iterative procedure for
%               estimating the shape. Scalar. The
%               iterative procedures stops when the relative difference of
%               a certain output matrix is smaller than itertol in two consecutive
%               iterations. The default value of pa.tol is 1e-12.
%      pa.zerotol = tolerance value to declare all input values equal to 0
%               in the eigenvalues restriction routine (file restreigen.m)
%               or in the final reconstruction of covariance matrices.
%               The default value of zerotol is 1e-10.
%          pa.msg = boolean which if set equal to true enables to monitor
%               the relative change of the estimates of lambda Gamma and
%               Omega in each iteration. The defaul value of pa.msg is
%               false, that is nothing is displayed in each iteration.
%   pa.userepmat  = scalar, which specifies whether to use implicit
%                   expansion or bsxfun.  pa.userepmat =2 implies implicit
%                   expansion, pa.userepmat=1 implies use of bsxfun. The
%                   default is to use implicit expansion (faster)
%                   if verLessThanFS(9.1) is false and bsxfun if MATLAB is
%                   older than 2016b.
%               Data Types - struct
%               Example - pa=struct; pa.cdet=10;
%
%           kk: number of mixture components. Integer vector. Integer
%               vector specifying the number of mixture components
%               (clusters) for which the BIC is to be calculated.
%               Vector. The default value of kk is 1:5.
%                 Example - 'kk',1:4
%                 Data Types - int16 | int32 | single | double
%
%
%      whichIC: type of information criterion. Character.
%               Character which specifies which information criteria must
%               be computed for each k (number of groups) and each value of
%               the restriction factors $c_{det}$ and $c_{shw}$.
%               Possible values for whichIC are:
%               'MIXMIX'  = a mixture model is fitted and to
%                   compute the information criterion the mixture
%                   likelihood is used. This option corresponds to the use of
%                   the Bayesian Information criterion (BIC). In output
%                   structure out just the matrix out.MIXMIX is given. If
%                   this input option is not specified, MIXMIX is assumed.
%               'MIXCLA'  = a mixture model is fitted but to compute the
%                   information criterion the classification likelihood is
%                   used. This option corresponds to the use of the
%                   Integrated Complete Likelihood (ICL). In output
%                   structure out just the matrix out.MIXCLA is given.
%               'CLACLA' =  everything is based on the classification
%                   likelihood. This information criterion will be called
%                   CLA. In output structure out just the matrix out.CLACLA
%                   is given.
%                 Example - 'whichIC','CLACLA'
%                 Data Types - character
%
%        alpha: global trimming level. Fraction or number of observations
%               which have to be trimmed. alpha is a scalar between 0 and
%               0.5 or an integer specifying the number of observations to
%               be trimmed. If alpha = 0 all observations are considered.
%               More in detail, if 0 < alpha < 1 clustering is based on
%               h = fix(n*(1-alpha)) observations. Else if alpha is an
%               integer greater than 1 clustering is based on h=n-floor(alpha).
%               The default value of alpha which is used is 0.
%                 Example - 'alpha',0
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
% RandNumbForNini: Pre-extracted random numbers to initialize proportions.
%                Matrix. Matrix with size k-by-size(nsamp,1) containing the
%                random numbers which are used to initialize the
%                proportions of the groups. This option is effective just
%                if nsamp is a matrix which contains pre-extracted
%                subsamples and k is a scalar. The purpose of this option
%                is to enable the user to replicate the results.
%                The default value of RandNumbForNini is empty,
%                that is random numbers from uniform are used.
%                   Example - 'RandNumbForNini',''
%                   Data Types - single | double
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
%equalweights : cluster weights in the concentration and assignment steps.
%               Logical. A logical value specifying whether cluster weights
%               shall be considered in the concentration, assignment steps
%               and computation of the likelihood. The default value is false.
%                 Example - 'equalweights',true
%                 Data Types - Logical
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
%       plots : Plot on the screen. Scalar. If plots = 1, three plots are produced.
%               The first is two panels plot. Left panel contains the
%               Information criterion as function of k for each value of
%               $c_{det}$ given the best value of $c_{shw}$ to the best
%               found value. The kind of information criterion which is
%               dispalyed depends on input option whichIC. It can be BIC
%               for mixture likelihood (MIXMIX) or ICL (MIXCLA) or BIC
%               for classification likelihood (CLACLA).
%               The second plot is a heatmap of the information criterion
%               for the values of $c_{det}$ and $c_{shw}$ which have been
%               considered.
%               The third plot is a two panels plot. The panel on the left
%               contains the values of IC for $c_{shb}$ given the best
%               values of $c_{det}$ and $c_{shw}$. This plot enables to
%               choose the best value of $c_{shb}$. The panel on the right
%               contains the values of IC in correspondence of the cases of
%               no rotation, equal rotation and varying rotation.
%                 Example - 'plots',1
%                 Data Types - single | double
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
%               available number of cores. REMARK 1: up to R2013b, there
%               was a limitation on the maximum number of cores that could
%               be addressed by the parallel processing toolbox (8 and,
%               more recently, 12). From R2014a, it is possible to run a
%               local cluster of more than 12 workers.
%               REMARK 2: Unless you adjust the cluster profile, the
%               default maximum number of workers is the same as the
%               number of computational (physical) cores on the machine.
%               REMARK 3: In modern computers the number of logical cores
%               is larger than the number of physical cores. By default,
%               MATLAB is not using all logical cores because, normally,
%               hyper-threading is enabled and some cores are reserved to
%               this feature.
%               REMARK 4: It is because of Remarks 3 that we have chosen as
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
%  cleanpool :  clean pool. Scalar. cleanpool is 1 if the parallel pool has
%               to be cleaned after the execution of the routine. Otherwise
%               it is 0. The default value of cleanpool is 1. Clearly this
%               option has an effect just if previous option numpool is >
%               1.
%                 Example - 'cleanpool',1
%                 Data Types - single | double
%
%       msg  :  Message on the screen. Scalar. Scalar which controls whether
%               to display or not messages about code execution.
%                 Example - 'msg',1
%                 Data Types - single | double
%
%      nocheck: Check input arguments. Scalar. If nocheck is equal to 1
%               no check is performed on matrix Y. As default nocheck=0.
%                 Example - 'nocheck',10
%                 Data Types - single | double
%
%       Ysave : save input matrix. Scalar. Scalar that is set to 1 to
%               request that the input matrix Y is saved into the output
%               structure out. Default is 1, that is  matrix Y is saved
%               inside output structure out.
%                 Example - 'Ysave',1
%                 Data Types - single | double
%
%   UnitsSameGroup : Units with same labels. Numeric vector.
%                   List of the units which must (whenever possible)
%                   have the same label.  For example if
%                   UnitsSameGroup=[20 26], it means that group which contains
%                   unit 20 is always labelled with number 1. Similarly,
%                   the group which contains unit 26 is always labelled
%                   with number 2, (unless it is found that unit 26 already
%                   belongs to group 1). In general, group which contains
%                   unit UnitsSameGroup(r) where r=2, ...length(kk)-1 is
%                   labelled with number r (unless it is found that unit
%                   UnitsSameGroup(r) has already been assigned to groups
%                   1, 2, ..., r-1). The default value of UnitsSameGroup is
%                   '' that is consistent labels are not imposed.
%                 Example - 'UnitsSameGroup',[12 20]
%                 Data Types - single | double
%
%
%       Remark: The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%
%
%  Output:
%
%         out:   structure which contains the following fields:
%
%                out.CLACLA = 3D array of size 5-times-8--times-8 if kk and
%                   pa.cdet and pa.cshw are not specififed else it is a 3D
%                   array of size
%                   length(kk)-times-length(pa.cdet)-times-length(pa.shw)
%                   containinig the value of the penalized classification
%                   likelihood. This output is present only if 'whichIC' is
%                   'CLACLA'.
%
%                out.IDXCLA = 3D array of size 5-times-8--times-8 if kk and
%                   pa.cdet and pa.cshw are not specififed else it is a 3D
%                   array of size
%                   length(kk)-times-length(pa.cdet)-times-length(pa.shw)
%                   containinig the assignment of each unit using
%                   the classification model. This output is present only
%                   if 'whichIC' is 'CLACLA'.
%
%                out.MIXMIX = 3D array of size 5-times-8-times-8 if kk and
%                   pa.cdet and pa.cshw are not specififed else it is a 3D
%                   array of size
%                   length(kk)-times-length(pa.cdet)-times-length(pa.shw)
%                   containinig the value of the penalized mixture
%                   likelihood. This output is present only if 'whichIC' is
%                   'MIXMIX'.
%
%                out.MIXCLA =D array of size 5-times-8-times-8 if kk and
%                   pa.cdet and pa.cshw are not specififed else it is a 3D
%                   array of size
%                   length(kk)-times-length(pa.cdet)-times-length(pa.shw)
%                   containinig the value of ICL. This
%                   output is present only if 'whichIC' is 'MIXCLA'.
%
%
%                out.IDXMIX = 3D array of size 5-times-8--times-8 if kk and
%                   pa.cdet and pa.cshw are not specififed else it is a 3D
%                   array of size
%                   length(kk)-times-length(pa.cdet)-times-length(pa.shw)
%                   containinig the assignment of each unit using using
%                   the mixture model. This output is present only if
%                   'whichIC' is 'MIXMIX'.
%
%                out.kk = vector containing the values of k (number of
%                   components) which have been considered. This  vector
%                   is equal to input optional argument kk if kk had been
%                   specified else it is equal to 1:5.
%
%                out.ccdet = vector containing the values of $c_det$ (values of the
%                   restriction factor for ratio of the determinant) which
%                   have been considered. This vector is equal to input
%                   argument pa.cdet else it is equal to [1, 2, 4, 8, 16,
%                   32, 64, 128].
%
%                out.kbest = scalar containing optimal value of k.
%

%                out.ccshw = vector containing the values of $c_det$ (values of the
%                   restriction factor for ratio of the shape elements
%                   inside each group) which have been considered. This
%                   vector is equal to input argument pa.cdet else it is
%                   equal to [1, 2, 4, 8, 16, 32, 64, 128].
%
%                out.cdetbest = scalar containing optimal value of
%                   restriction among determinants.
%
%                out.cshwbest = scalar containing optimal value of
%                   restriction among shape elements inside each group.
%
%                out.cshbbest = scalar containing optimal value of
%                   restriction among ordered shape elements across groups.
%
%                out.BICbest = scalar containing optimal value of
%                   BIC (smallest value of BIC).
%
%                out.modelbest = character identifying the best
%                   specification among the 14 GPCMs.
%
%                out.alpha = scalar containing the value of the trimming

%                out.modelbest = character identifying the best
%                   specification among the 14 GPCMs.
%                   level which has been used.
%
%               out.idx  = n-by-1 vector containing assignment of each unit to
%                       each of the k groups. Cluster names are integer
%                       numbers from 1 to k. 0 indicates trimmed
%                       observations. Note that this is the assignment
%                       using cshb=128;
%
%               out.idxcshbbest  = n-by-1 vector containing assignment of
%                       each unit to each of the k groups. Cluster names
%                       are integer numbers from 1 to k. 0 indicates
%                       trimmed observations. Note that this is the
%                       assignment which used best value of cshb;
%
%               out.pa = strcture containing the the gpcm parameters
%                   (pa.pars, pa.shb, pa.shw, pa.cdet) for the best model.
%
%                out.Y  = Original data matrix Y. The field is present if
%                   option Ysave is set to 1.
%
%
% See also tclustIC, tclust, tclustICsol, tclustICplot, carbikeplot
%
% References:
%
% Cerioli, A., Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani M. (2017),
% Finding the Number of Groups in Model-Based Clustering via Constrained
% Likelihoods, "Journal of Computational and Graphical Statistics", pp. 404-416,
% https://doi.org/10.1080/10618600.2017.1390469
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('tclustICgpcm')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Automatic choice of k, cdet, cshw, cshb for geyser data.
    % Use a small number fo subsets in order to reduce execution time.
    rng(100)
    nsamp=20;
    Y=load('geyser2.txt');
    % If no trimming is used 4 groups are found.
    out=tclustICgpcm(Y,'nsamp',nsamp);
%}

%{
    % Example of use of input options alpha and typeIC.
    Y=load('geyser2.txt');
    pa=struct;
    pa.pars='VVV';
    rng(100)
    nsamp=20;
    alpha=0.1;
    whichIC='MIXMIX';
    outIC=tclustICgpcm(Y,'pa',pa,'alpha',alpha,'whichIC',whichIC);
%}

%{
    % Automatic choice of k in an example with 3 components and prefixed overlap.
    rng('default') % Reinitialize the random number generator to its startup configuration
    rng(20000);
    ktrue=3;
    % n = number of observations
    n=150;
    % v= number of dimensions
    v=2;
    % Imposed average overlap
    BarOmega=0.04;
    restrfact=5;

    outg=MixSim(ktrue,v,'BarOmega',BarOmega, 'restrfactor',restrfact);
    % data generation given centroids and cov matrices
    [Y,id]=simdataset(n, outg.Pi, outg.Mu, outg.S);

    % Number of subsamples to extract (option nsamp) is very small
    % therefore a great variability is allowed.
    outIC=tclustICgpcm(Y,'nsamp',50);
%}

%{
   % An example with input options kk pa.cdet and pa.shw.
    Y=load('geyser2.txt');
    nsamp=100;
    pa=struct;
    pa.cdet=[2 4];
    pa.shw=[8 16 32];
    kk=[2 3 4 6];
    out=tclustICgpcm(Y,'pa',pa,'cleanpool',false,'plots',0,'alpha',0.1,'whichIC','CLACLA','kk',kk,'nsamp',nsamp);
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

plots=1;
nsamp=500;
alpha=0;
kk=1:5;
whichIC='MIXMIX';
msg=1;
cc=[1 2 4 8 16 32 64 128];
cleanpool=false;
UnitsSameGroup='';
RandNumbForNini='';
pa=struct;
pa.pars='VVV';

options=struct('pa',pa,'kk',kk,'whichIC',whichIC,'alpha',alpha,'nsamp',nsamp,'plots',plots,'nocheck',0,...
    'msg',msg,'Ysave',1,'refsteps',refsteps,'equalweights',equalweights,...
    'reftol',reftol,'startv1',startv1,...
    'UnitsSameGroup',UnitsSameGroup,'RandNumbForNini',RandNumbForNini,...
    'numpool',numpool, 'cleanpool', cleanpool);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:tclustIC:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:tclustICgpcm:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end


if nargin > 2
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    pa=options.pa;
    alpha=options.alpha;
    kk=options.kk;
    nsamp=options.nsamp;        % Number of subsets to extract
    plots=options.plots;        % Plot of the resulting classification
    equalweights=options.equalweights;    % Specify if assignment must take into account the size of the groups
    
    refsteps=options.refsteps;
    reftol=options.reftol;
    msg=options.msg;            % Scalar which controls the messages displayed on the screen
    
    whichIC=options.whichIC;
    cleanpool=options.cleanpool;
    numpool=options.numpool;
    UnitsSameGroup=options.UnitsSameGroup;
    RandNumbForNini=options.RandNumbForNini;
end

if strcmp(whichIC,'MIXMIX')
    typeIC=2;
elseif strcmp(whichIC,'MIXCLA')
    typeIC=1;
elseif strcmp(whichIC,'CLACLA')
    typeIC=0;
else
    warning('FSDA:tclustICgpcm:WrongOpt','Supplied string for whichIC is not supported.')
    error('FSDA:tclustICgpcm:WrongIC','Specified information criterion is not supported: possible values are ''MIXMIX'' , ''MIXCLA'',  ''CLACLA''')
end


if isfield(pa,'cdet')
    ccdet=pa.cdet;
else
    ccdet=cc;
end

if isfield(pa,'shw')
    ccshw=pa.shw;
else
    ccshw=cc;
end

if ~isfield(pa,'pars')
    pa.pars='VVV';
end

lcdet=length(ccdet);
lcshw=length(ccshw);
lkk=length(kk);

if typeIC==2
    MIXMIX=zeros(lkk,lcdet,lcshw);
    IDXMIX=cell(lkk,lcdet,lcshw);
elseif typeIC==1
    MIXCLA=zeros(lkk,lcdet,lcshw);
    IDXMIX=cell(lkk,lcdet,lcshw);
else %  typeIC==0
    CLACLA=zeros(lkk,lcdet,lcshw);
    IDXCLA=cell(lkk,lcdet,lcshw);
end


% Prepare rownames and colsnames for table which will contain
% in the rows the number of groups and in the columsn the values of c
rownamesIC=strcat(cellstr(repmat('k=',lkk,1)), cellstr(num2str(kk')));
rownamesIC=regexprep(rownamesIC,' ','');
colnamesIC=strcat(cellstr(repmat('cdet_',lcdet,1)), cellstr(num2str(ccdet')));
colnamesIC=regexprep(colnamesIC,' ','');


%% Preapare the pool (if required)
Cnsampall=cell(lkk,1);
gRandNumbForNiniall=cell(lkk,1);
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
    Cnsampall{k}=Cnsamp;
    
    if isempty(RandNumbForNini)
        % For each value of k extract random numbers to initialize proportions
        % once and for all
        gRandNumbForNini=rand(seqk,nsamp);
    else
        gRandNumbForNini=RandNumbForNini;
    end
    gRandNumbForNiniall{k}=gRandNumbForNini;
    pasel=pa;
    pasel.shb=128;
    
    for cshw=1:lcshw
        pasel.shw=ccshw(cshw);
        
        parfor (cdet=1:lcdet , numpool)
            % for cdet=1:lcc
            % Select  value for restriction among determinants
            pasel1=pasel;
            pasel1.cdet=ccdet(cdet);
            
            
            % columns = restr factor c_det
            % Third dimension values of c_shw
            % rows = number of groups
            % tclust using mixtures
            if typeIC>0
                outMixt=tclust(Y,seqk,alpha,pasel1,'nsamp',Cnsamp,'plots',0,'msg',0,'mixt',2, ...
                    'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                    'reftol',reftol,'RandNumbForNini',gRandNumbForNini);
                
                IDXMIX{k,cdet,cshw}=outMixt.idx;
                if typeIC==2
                    MIXMIX(k,cdet,cshw)=outMixt.MIXMIX;
                end
                if typeIC==1
                    MIXCLA(k,cdet,cshw)=outMixt.MIXCLA;
                end
            end
            
            if typeIC==0
                % tclust using classification likelihood
                outCla=tclust(Y,seqk,alpha, pasel1,'nsamp',Cnsamp,'plots',0,'msg',0, ...
                    'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                    'reftol',reftol,'RandNumbForNini',gRandNumbForNini);
                CLACLA(k,cdet,cshw)=outCla.CLACLA;
                IDXCLA{k,cdet,cshw}=outCla.idx;
            end
        end
    end
    if msg==1
        disp(['k=' num2str(seqk)])
    end
end

%% Close pool and show messages if required
if cleanpool==true
    delete(gcp);
end

out=struct;

% Call tclustreg with **V **E **I to decide about best rotation
modelb=zeros(3,1);
idxb=zeros(n,3);
models=cell(3,1);
typerot={'I';'E';'V'};



if typeIC==0 % CLACLA
    out.CLACLA=CLACLA;
    verMatlab = true;
    if verMatlab == false
        % out.(IC) is also given in table format
        out.CLACLAtable=array2table(CLACLA,'RowNames',rownamesIC,...
            'VariableNames',matlab.lang.makeValidName(colnamesIC));
    end
    
    % Store whenever possible consistent labels
    if ~isempty(UnitsSameGroup)
        IDXCLA=ClusterRelabel(IDXCLA,UnitsSameGroup);
    end
    
    out.IDXCLA=IDXCLA;
    % Store vectors kk and cc inside output structure out
    out.kk=kk;
    
    out.ccdet=ccdet;
    out.ccshw=ccshw;
    
    selIC=CLACLA;
    nameselIC='CLACLA';
    
    [kbest, cdetbest, cshwbest, bestBIC]=selICplot(selIC,ccdet,ccshw,kk,nameselIC,plots);
    
    % cdetbest and cshwbest are using in the refiing step to decide about
    % the best type of rotation
    pasel.cdet=cdetbest;
    pasel.shw=cshwbest;
    
    modelroot=pasel.pars;
    if cdetbest==1
        modelroot(1)='E';
    end
    
    
    pasel.k=kbest;
    pasel.v=v;
    
    
    % Find best estimate of cshbbest
    candshb=cc(cc<=cshwbest^((v-1)/v));
    
    modelshb=zeros(length(candshb),1);
    for jshb=1:length(candshb)
        pasel.shb=candshb(jshb);
        
        outCla=tclust(Y,kbest,alpha,pasel,'nsamp',Cnsampall{kk==kbest},'plots',0,'msg',0,'mixt',0, ...
            'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
            'reftol',reftol,'RandNumbForNini',gRandNumbForNiniall{kk==kbest});
        
        modelshb(jshb)=outCla.CLACLA;
    end
    
    [~,posminbestcshb]=min(modelshb);
    cshbbest=candshb(posminbestcshb);
    % Set best values of constraint across groups
    pasel.shb=cshbbest;
    
    if cshbbest==1
        modelroot(2)='E';
    end
    
    models{1}=[modelroot(1:2) 'I'];
    models{2}=[modelroot(1:2) 'E'];
    models{3}=[modelroot(1:2) 'V'];
    
    
    for jrot=1:3
        pasel.pars=models{jrot};
        
        outCla=tclust(Y,kbest,alpha,pasel,'nsamp',Cnsampall{kk==kbest},'plots',0,'msg',0,'mixt',0, ...
            'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
            'reftol',reftol,'RandNumbForNini',gRandNumbForNiniall{kk==kbest});
        
        idxb(:,jrot)=outCla.idx;
        modelb(jrot)=outCla.CLACLA;
    end
else  % MIXMIX or MIXCLA store IDXMIX
    if ~isempty(UnitsSameGroup)
        IDXMIX=ClusterRelabel(IDXMIX,UnitsSameGroup);
    end
    out.IDXMIX=IDXMIX;
    if typeIC==2
        out.MIXMIX=MIXMIX;
        selIC=MIXMIX;
        nameselIC='MIXMIX';
    else
        out.MIXCLA=MIXCLA;
        selIC=MIXCLA;
        nameselIC='MIXCLA';
    end
    % Store vectors kk and cc inside output structure out
    out.kk=kk;
    out.ccdet=ccdet;
    out.ccshw=ccshw;
    
    [kbest, cdetbest, cshwbest, bestBIC]=selICplot(selIC,ccdet,ccshw,kk,nameselIC,plots);
    %    [kbest, cdetbest, cshwbest, bestBIC]=selICplot(selIC,cc,lcc,kk,lkk,nameselIC,xkk,LineWidth,slintyp,styp,legstr);
    
    pasel.cdet=cdetbest;
    pasel.shw=cshwbest;
    
    modelroot=pasel.pars;
    if cdetbest==1
        modelroot(1)='E';
    end
    
    
    pasel.k=kbest;
    pasel.v=v;
    
    % Find best estimate of cshbbest
    candshb=ccshw(ccshw<=cshwbest^((v-1)/v));
    
    modelshb=zeros(length(candshb),1);
    for jshb=1:length(candshb)
        pasel.shb=candshb(jshb);
        
        outMixt=tclust(Y,kbest,alpha,pasel,'nsamp',Cnsampall{kbest},'plots',0,'msg',0,'mixt',2, ...
            'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
            'reftol',reftol,'RandNumbForNini',gRandNumbForNiniall{kbest});
        if typeIC==2
            modelshb(jshb)=outMixt.MIXMIX;
        end
        if typeIC==1
            modelshb(jshb)=outMixt.MIXCLA;
        end
    end
    
    [~,posminbestcshb]=min(modelshb);
    cshbbest=candshb(posminbestcshb);
    % Set best values of constraint across groups
    pasel.shb=cshbbest;
    
    if cshbbest==1
        modelroot(2)='E';
    end
    
    models{1}=[modelroot(1:2) 'I'];
    models{2}=[modelroot(1:2) 'E'];
    models{3}=[modelroot(1:2) 'V'];
    
    
    for jrot=1:3
        pasel.pars=models{jrot};
        outMixt=tclust(Y,kbest,alpha,pasel,'nsamp',Cnsampall{kbest},'plots',0,'msg',0,'mixt',2, ...
            'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
            'reftol',reftol,'RandNumbForNini',gRandNumbForNiniall{kbest});
        idxb(:,jrot)=outMixt.idx;
        
        if typeIC==2
            modelb(jrot)=outMixt.MIXMIX;
        end
        if typeIC==1
            modelb(jrot)=outMixt.MIXCLA;
        end
    end
end

[~,indminrot]=min(modelb);
pasel.pars=models{indminrot};

% Stope optimal values of the parameters
out.kbest=kbest;
out.cdetbest=cdetbest;
out.cshwbest=cshwbest;
out.cshbbest=cshbbest;
out.BICbest=bestBIC;
out.modelbest=models{indminrot};

out.idxcshbbest=idxb(:,indminrot);
% Store best classification before refining step for rotation
if typeIC>0
    out.idx=IDXMIX{kk==kbest,ccdet==cdetbest,ccshw==cshwbest};
else
    out.idx=IDXCLA{kk==kbest,ccdet==cdetbest,ccshw==cshwbest};
end

% Store best gpcm parameters
out.pa=pasel;


% Show refinements plots
if plots==1
    figure
    if length(candshb)>1
        nr=1;
        nc=2;
        subplot(nr,nc,1)
        plot(candshb,modelshb)
        xlabel('c_{shb}')
        ylabel('BIC to select best c_{shb}')
        title(['Best c_{shb}=' num2str(cshbbest)])
        subplot(nr,nc,2)
    else
        % if cshw is <=2 there is just one point and the plot is not ahown
    end
    plot(1:3,modelb,'o')
    xlabel('Type of rotation')
    ylabel('BIC to select best type of rotation')
    title(['Best rot =' typerot{indminrot}])
    xlim([1 3])
    set(gca,'XTick',1:3)
    set(gca,'XTickLabel',models);
    
    % Make the main BIC plot the current figure
    plBIC=findobj('type','figure','Name','BIC');
    if ~isempty(plBIC)
        figure(plBIC(1))
    end
    
end

% Store trimming level which has been used
out.alpha=alpha;
% Store original matrix
out.Y=Y;
end



function [kbest,cdetbest,cshwbest,BICbest]=selICplot(selIC,cdet,cshw,kk,nameselIC,plots)

[valmin,posmin]=min(selIC,[],'all','linear');
[lkk,lcdet,lcshw]=size(selIC);

% bestk, bestcdet and best cshw are respectively best values for number of
% groups (bestk), restriction among determinants (bestcdet) and restriction
% among the shape elements inside each group (bestcshw)
[bestk,bestcdet,bestcshw]=ind2sub([lkk lcdet lcshw],posmin);
cdetbest=cdet(bestcdet);
cshwbest=cshw(bestcshw);
kbest=kk(bestk);
% Store value of the minimum BIC
BICbest=valmin;

if plots==1
    figure('Name','BIC')
    % set line width of the trajectories of BIC
    LineWidth=1;
    % Define marker type
    styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
    lcc=max([lcdet,lcshw]);
    styp=repmat(styp,ceil(200/lcc),1);
    % Define line type
    slintyp={'-';'--';':';'-.'};
    slintyp=repmat(slintyp,ceil(200/lcc),1);
    % Define legend entries
    xkk=0:(1/(length(kk)-1)):1;
    
    nr=1;
    nc=2;
    subplot(nr,nc,1)
    selIC2D=selIC(:,:,bestcshw);
    plot1=plot(kk',selIC2D,'LineWidth',LineWidth);
    title([nameselIC ' | c_{shw}=' num2str(cshw(bestcshw))])
    % Add labels for the best value of cdet for each k
    cmin=zeros(lcdet,1);
    for j=1:lkk
        [~,posj]=min(selIC2D(j,:));
        cmin(j)=cdet(posj);
        text(xkk(j),0.98,['c_{det}=' num2str(cmin(j))],'Units','Normalized')
    end
    
    % Set line type and markers
    set(plot1,{'LineStyle'},slintyp(1:lcdet));
    set(plot1,{'Marker'},styp(1:lcdet))
    xlabel('Number of groups')
    set(gca,'xtick',kk)
    
    a=cell(lcdet,1);
    a(:)={'c_{det}='};
    if isrow(cdet)
        legstrcdet=strcat(a, cellstr(num2str(cdet')));
    else
        legstrcdet=strcat(a, cellstr(num2str(cdet')));
    end
    legend(legstrcdet,'location','best');
    
    subplot(nr,nc,2)
    selIC2D=squeeze(selIC(:,bestcdet,:));
    plot1=plot(kk',selIC2D,'LineWidth',LineWidth);
    title([nameselIC ' | cdet=' num2str(cdet(bestcdet))])
    % Add labels for the best value of cshw for each k
    cmin=zeros(lcshw,1);
    for j=1:lkk
        [~,posj]=min(selIC2D(j,:));
        cmin(j)=cshw(posj);
        text(xkk(j),0.98,['c_{shw}=' num2str(cmin(j))],'Units','Normalized')
    end
    
    % Set line type and markers
    set(plot1,{'LineStyle'},slintyp(1:lcshw));
    set(plot1,{'Marker'},styp(1:lcshw))
    xlabel('Number of groups')
    set(gca,'xtick',kk)
    
    a=cell(lcshw,1);
    a(:)={'c_{shw}='};
    if isrow(cshw)
        legstrcshw=strcat(a, cellstr(num2str(cshw')));
    else
        legstrcshw=strcat(a, cellstr(num2str(cshw')));
    end
    legend(legstrcshw,'location','best');
    % set(plot1,'Tag','BIC')
    disp('The labels of in the top part of Figure named BIC denote the values of $c_{det}$ ($c_{shw}$) for which IC is minimum')
    
    figure
    selIC2Dbestk=squeeze(selIC(bestk,:,:));
    heatmap(cshw,cdet,selIC2Dbestk)
    xlabel('c_{shw}')
    ylabel('c_{det}')
    title(['Heatmap for k=' num2str(kk(bestk)) '. Best c_{shw}=' ...
        num2str(cshw(bestcshw))  ', best c_{det}=' num2str(cdet(bestcdet))...
        ' min(' nameselIC ')=' num2str(BICbest,'%1.1f') ])
    
end
end

%FScategory:CLUS-RobClaMULT