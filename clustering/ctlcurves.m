function out  = ctlcurves(Y, varargin)
%tclustIC computes Classification Trimmed Likelihood Curves
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
%            Y: Input data. Matrix. Data matrix containing n observations on v variables
%               Rows of Y represent observations, and columns represent
%               variables. Observations (rows) with missing (NaN) or or
%               infinite (Inf) values will automatically be excluded from
%               the computations.
%                 Data Types -  double
%
%  Optional input arguments:
%
%           kk: number of mixture components. Integer vector. Integer
%               vector specifying the number of mixture components
%               (clusters) for which trimmed likelihoods are calculated.
%               Vector. The default value of kk is 1:5.
%                 Example - 'kk',1:4
%                 Data Types - int16 | int32 | single | double
%        alpha: trimming level to monitor. Vector. Vector which specifies the
%               values of trimming levels which have to be considered.
%               For example is alpha=[0 0.05 0,1] ctlcurves considers these 3
%               values of trimming level.
%               The default for alpha is vector
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
%                subsamples and k is a scalat. The purpose of this option
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
%               and computation of the likelihood.
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
%       plots : Plot on the screen. Scalar. If plots = 1, a plot of the
%               CTLcurves is shown on the screen. If input option bands is
%               not empty confidence bands are also shown.
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
%               it is 0. The default value of cleanpool is 0. Clearly this
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
%       Ysave : save input matrix. Scalar.
%               Scalar that is set to 1 to request that the input matrix Y
%               is saved into the output structure out. Default is 1, that
%               is  matrix Y is saved inside output structure out.
%                 Example - 'Ysave',1
%                 Data Types - single | double
%
%
%       cshape :    constraint to apply to each of the shape matrices.
%                   Scalar greater or equal than 1. This options only works is 'restrtype' is
%                   'deter'.
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
%                out.CLACLA = matrix of size 5-times-8 if kk and cc are not
%                   specififed else it is a matrix of size length(kk)-times
%                   length(cc) containinig the value of the penalized
%                   classification likelihood. This output is present only
%                   if 'whichIC' is 'CLACLA' or 'whichIC' is 'ALL'.
%                out.CLACLAtable = same output of CLACLA but in MATLAB
%                   table format (this field is present only if your MATLAB
%                   version is not<2013b).
%
%                out.IDXCLA = cell of size 5-times-8 if kk and cc are not
%                   specififed else it is a cell of size length(kk)-times
%                   length(cc). Each element of the cell is a vector of
%                   length n containinig the assignment of each unit using
%                   the classification model. This output is present only
%                   if 'whichIC' is 'CLACLA' or 'whichIC' is 'ALL'.
%
%                out.MIXMIX = matrix of size 5-times-8 if kk and cc are not
%                   specififed else it is a matrix of size length(kk)-times
%                   length(cc) containinig the value of the penalized
%                   mixture likelihood. This output is present only if
%                   'whichIC' is 'MIXMIX' or 'whichIC' is 'ALL'.
%                out.MIXMIXtable = same output of MIXMIX but in MATLAB
%                   table format (this field is present only if your MATLAB
%                   version is not<2013b).
%
%                out.MIXCLA = matrix of size 5-times-8 if kk and cc are not
%                   specififed else it is a matrix of size length(kk)-times
%                   length(cc) containinig the value of the ICL. This
%                   output is present only if 'whichIC' is 'MIXCLA' or
%                   'whichIC' is 'ALL'.
%                out.MIXCLAtable = same output of MIXCLA but in MATLAB
%                   table format (this field is present only if your MATLAB
%                   version is not<2013b).
%
%                out.IDXMIX = cell of size 5-times-8 if kk and cc are not
%                   specififed else it is a cell of size length(kk)-times
%                   length(cc). Each element of the cell is a vector of
%                   length n containinig the assignment of each unit using
%                   the mixture model. This output is present only if
%                   'whichIC' is 'MIXMIX', 'MIXCLA' or 'ALL'.
%
%                out.kk = vector containing the values of k (number of
%                   components) which have been considered. This  vector
%                   is equal to input optional argument kk if kk had been
%                   specified else it is equal to 1:5.
%
%                out.alphavec = vector containing the values of c (values of the
%                   restriction factor) which have been considered. This
%                   vector is equal to input optional argument alphavec if alphavec had
%                   been specified else it is equal to [0, 0.01, ..., 0.10].
%
%                out.restractor = scalar containing the trimming level
%                   which has been used.
%
%                out.OptimalK = scalar, optimal number of clusters, stored
%                   as a positive integer value. This output is present
%                   only if optional input argument is true.
%
%                out.Optimalalpha = scalar, optimal value of trimming. This
%                   output is present only if optional input argument is true.
%
%                out.Y  = Original data matrix Y. The field is present if
%                   option Ysave is set to 1.
%
%
% More About:
%
% These curves show the values of the trimmed classification
% (log-)likelihoods when altering the trimming proportion alpha and the
% number of clusters k. The careful examination of these curves provides
% valuable information for choosing these parameters in a clustering
% problem. For instance, an appropriate k to be chosen is one that we do
% not observe a clear increase in the trimmed classification likelihood
% curve for k with respect to the k+1 curve for almost all the range of
% alpha values. Moreover, an appropriate choice of parameter alpha may be
% derived by determining where an initial fast increase of the trimmed
% classification likelihood curve stops for the final chosen k. A more
% detailed explanation can be found in Garc<ed>a-Escudero et al. (2010).
%
%
% See also tclust, tclustICsol, tclustICplot
%
%
% References:
%
% Garcia-Escudero, L.A.; Gordaliza, A.; Matran, C. and Mayo-Iscar, A.
% (2011), "Exploring the number of groups in robust model-based
% clustering." Statistics and Computing, Vol. 21, pp. 585–599.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('ctlcurves')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

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
%alphaTrim=0:0.04:0.20;

cleanpool=false;
RandNumbForNini='';
% cshape. Constraint on the shape matrices inside each group which works only if restrtype is 'deter'
cshape=10^10;
restrfactor=100;
mixt=0;
bands=true;

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    options=struct('kk',kk,'alphavec',alphaTrim,'whichIC',whichIC,'alpha',alphaTrim,'nsamp',nsamp,'plots',plots,'nocheck',0,...
        'msg',msg,'Ysave',1,'refsteps',refsteps,'equalweights',equalweights,...
        'reftol',reftol,'startv1',startv1,'restrtype',restr,'cshape',cshape,...
        'RandNumbForNini',RandNumbForNini,...
        'numpool',numpool, 'cleanpool', cleanpool);
    
    
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
        error('FSDA:ctlcurves:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    
    restr=options.restrtype;
    alphaTrim=options.alphavec;
    kk=options.kk;
    nsamp=options.nsamp;        % Number of subsets to extract
    plots=options.plots;        % Plot of the resulting classification
    equalweights=options.equalweights;    % Specify if assignment must take into account the size of the groups
    
    refsteps=options.refsteps;
    reftol=options.reftol;
    msg=options.msg;            % Scalar which controls the messages displayed on the screen
    
    mixt=options.mixt;
    cleanpool=options.cleanpool;
    numpool=options.numpool;
    RandNumbForNini=options.RandNumbForNini;
    cshape=options.cshape;
    restrfactor=options.restrfactor;
    bands=options.bands;
end

lkk=length(kk);
lalpha=length(alphaTrim);

MuVal = cell(lkk,lalpha);
SigmaVal = MuVal;
PiVal = MuVal;
IDX=MuVal;
CTLtab=zeros(lkk,lalpha);

%% Preapare the pool (if required)
pariter=0;
[numpool,tstart, progbar, usePCT, usematlabpool] = PoolPrepare(numpool, pariter, UserOptions);

CnsampAll=cell(lkk,1);
gRandNumbForNiniAll=CnsampAll;

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
    
    if isempty(RandNumbForNini)
        % For each value of k extract random numbers to initialize proportions
        % once and for all
        gRandNumbForNini=rand(seqk,nsamp);
    else
        gRandNumbForNini=RandNumbForNini;
    end
    
    if bands==true
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
        
        % columns = values of alpha
        % rows = values of k
        
        % Store centroids
        MuVal{k,j} = outtc.muopt;
        SigmaVal{k,j} = outtc.sigmaopt;
        
        % Store classification
        IDX{k,j} = outtc.idx;
        
        % Store objective function
        CTLtab(k,j) = outtc.obj;
    end
    if msg==1
        disp(['k=' num2str(seqk)])
    end
end

if bands==true
    % Simulation to create the bands
    nsimul=100;
    BandsCTL=zeros(lkk,lalpha,nsimul);
    for k=1:lkk  % loop for different values of k (number of groups)
        if msg==1
            disp(['Bands k=' num2str(k)])

        end
        
        seqk=kk(k);
        for j=1:lalpha
       
            ktrue = length(PiVal{seqk, j});
            Mutrue = MuVal{seqk, j};
            Mutrue=Mutrue(1:ktrue,:);
            Sigmatrue = SigmaVal{seqk, j};
            Sigmatrue = Sigmatrue(:,:,1:ktrue);
            Pitrue=PiVal{seqk, j};
            alphaTrimj=alphaTrim(j);
            ngood=round(n*(1-alphaTrimj));
            nout=n-ngood;
            CnsampAllk=CnsampAll{seqk};
            gRandNumbForNiniAllk=gRandNumbForNiniAll{seqk};
            parfor zz = 1:nsimul
                [Ysim]=simdataset(ngood, Pitrue, Mutrue, Sigmatrue,'noiseunits', nout);
                if size(Ysim,1)<n
                    Yadd=repmat(Ysim(end,:),n-size(Ysim,1),1);
                    Ysim=[Ysim;Yadd];
                end
                outtcSIM=tclust(Ysim,seqk,alphaTrimj,restrfactor,'nsamp',CnsampAllk,'plots',0,'msg',0,'mixt',mixt, ...
                    'restrtype',restr,'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                    'reftol',reftol,'RandNumbForNini',gRandNumbForNiniAllk,'cshape',cshape);
                
                BandsCTL(k, j, zz) = outtcSIM.obj;
            end
        end
    end
end

out=struct;
out.Mu=MuVal;
out.Sigma=SigmaVal;
out.PiVal=PiVal;
out.IDX=IDX;
out.CTLtabl=CTLtab;
out.BandsCTL=BandsCTL;
gamma=0.25;

likLB = zeros(length(kk),length(alphaTrim));
lik050 = likLB;
likUB =likLB;
for k=1:length(kk)  % loop for different values of k (number of groups)
    parfor j=1:length(alphaTrim)
        likLB(k,j) = quantile(BandsCTL(k,j,:), gamma);
        lik050(k,j) = median(BandsCTL(k,j,:));
        likUB(k,j) = quantile(BandsCTL(k,j,:), 1- gamma);
    end
end

out.likLB=likLB;
out.lik050=lik050;
out.likUB=likUB;

for j = 1:length(kk)-1
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
        alphafin = alphaBest;
        kfin = j;
        conv = 1;
        break
    end
end


if conv == 1
    parSelec = [kfin alphafin];
else
    parSelec = [999 999];
end

out.parSelec=parSelec;
out.Optimalalpha=alphafin;
out.OptimalK=kfin;
% Store best classification
out.idx=IDX{kfin,jalpha};

if plots==1
    if bands==1
        linetype1 = {'-.','-.','-.','-.','-.'};
        % linetype = {'-','-','-','-','-'};
        color = {'r','g','b','c','k'};
        LineWidth = 1;
        hold('on')
        for i = 1:length(kk)
            plot(alphaTrim,likLB(i,:), 'LineStyle',linetype1{i}, 'Color', color{i}, 'LineWidth', LineWidth)
            % plot(alphaTrim,lik050(i,:), 'LineStyle',linetype{i}, 'Color', color{i}, 'LineWidth', LineWidth)
            text(alphaTrim(end),lik050(i,end),[' k = ' num2str(i)],'FontSize',16, 'Color', color{i})
            plot(alphaTrim,likUB(i,:), 'LineStyle',linetype1{i}, 'Color', color{i}, 'LineWidth', LineWidth)
        end
    else
        plot(alphaTrim, CTLtab)
        for i = kk
            text(alphaTrim,CTLtab(i,:),num2str(i*ones(length(alphaTrim),1)),'FontSize', 14)
        end
    end
    xlabel('Trimming level alpha')
    ylabel('Log likelihood')
    set(gca,'XTick',alphaTrim);
end

%% Close pool and show messages if required
PoolClose(cleanpool, tstart, progbar, usePCT, usematlabpool);


end


