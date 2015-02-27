function [out]  = MixSim(k,v,varargin)
%MixSim generates k clusters in v dimensions with given overlap
%
%<a href="matlab: docsearchFS('mixsim')">Link to the help function</a>
%
%   MixSim(k,v) generates k groups in v dimensions. It is possible to
%   control the average and maximum or standard deviation of overlapping.
%
%  Background: Given two generic clusters i and j with i ne j =1, ..., k,
%  indexed by \phi(x; \mu_i,\Sigma_i) and \phi(x; \mu_j,\Sigma_j) with
%  probabilities of occurrence \pi_i and \pi_j, the misclassification
%  probability with respect to cluster i (that is conditionally on x
%  belonging to cluster i,  which is called  w_j|i) is defined as
%  Pr[ \pi_i \phi(x;\mu_i,\Sigma_i) < \pi_j \phi(x;\mu_j,\Sigma_j)].
%  The matrix containing the misclassification probabilities w_j|i is
%  called OmegaMap
%  The probability of overlapping between groups i and j is given by
%            w_j|i + w_i|j          i,j=1,2, ..., k
%  The diagonal elements of OmegaMap are equal to 1.
%  The average overlap (which in the code is called below BarOmega) is
%  defined as the sum of the off diagonal elements of OmegaMap (matrix of
%  misclassification probabilities) divided by 0.5*k*(k-1)
%  The maximum overlap (which in the code is called MaxOmega) is defined as
%  max(w_j|i + w_i|j)  i \ne j=1,2, ..., k
%  The probability of overlapping w_j|i is nothing but the cdf of a linear
%  combination of non central chi^2 distributions with 1 degree of freedom
%  + a linear combination of N(0,1) evaluated in a
%  point c.  The coefficients of the linear combinations of non central
%  chi^2 and N(0,1) depend on the eigenvalues and eigenvectors of matrix
%  \Sigma_j|i = \Sigma^{0.5}_i \Sigma^{-1}_j \Sigma^{0.5}_i.
%  Point c depends on the same eigenvalues and eigenvectors
%  of matrix \Sigma_j|i, the mixing proportions \pi_i and \pi_j and the
%  determinants |\Sigma_i| and |\Sigma_j|
%  This probability is computed using routine ncx2mixtcdf
%
%  Required input arguments:
%
%            k: scalar, number of groups (components)
%            v: scalar, number of dimensions (variables).
%
%  Optional input arguments:
%
%    BarOmega : scalar, value of desired average overlap. The default value is ''
%    MaxOmega : scalar, value of desired maximum overlap. If BarOmega is empty
%               the default value of MaxOmega is 0.15
%    StdOmega : scalar, value of desired standard deviation of overlap.
%               Remark1: The probability of overlapping between two
%               clusters i and j (i ne j =1, 2, ..., k), called pij, is defined as the
%               sum of the two misclassification probabilities
%               pij=w_{j|i} + w_{i|j}
%               Remark2: it is possible to specify up to two values among
%               BarOmega MaxOmega and StdOmega.
%         sph : scalar boolean which specifies covariance matrix structure
%               sph=false (default) ==> non-spherical,
%               sph=true            ==> spherical = const*I
%         hom : scalar boolean which specifies heterogeneous or homogeneous
%               clusters
%               hom=false (default) ==> heterogeneous,
%               hom=true            ==> homogeneous \Sigma_1 = ... = \Sigma_k
%         ecc : scalar in the interval (0, 1] which defines maximum eccentricity.
%               For example, if ecc=0.9 (default value), we require for
%               each group that sqrt(1 - minL / maxL) <= 0.9 where minL and
%               maxL are respectively the min and max eigenvalue of the cov
%               matrix of a particular group
%  restrfactor: scalar in the interval [1 \infty] which specifies the
%               maximum ratio to allow between the largest eigenvalue and
%               the smallest eigenvalue of the k covariance matrices which
%               are generated. The default value is ''. More in details if for example
%               restrfactor=10 after generating the covariance matrices we
%               check that the ratio
%               \[
%                 \frac{   \max_{l=1, \ldots, v} \max_{j=1, \ldots, k}  \lambda_l(\hat \Sigma_j)}{   \min_{l=1, \ldots, v} \min_{j=1, \ldots, k}  \lambda_l(\hat \Sigma_j)}.
%               \]
%               between the largest eigenvalue of the k cov matrices
%               and the smallest eigenvalue of the k cov matrices is not
%               larger than restrfactor. In order to apply this restriction
%               (which is typical of tclust.m) we call routine restreigen.m
%       PiLow : scalar, value of the smallest mixing proportion (if 'PiLow'
%               is not reachable with respect to k, equal proportions are
%               taken; PiLow = 1.0 implies equal proportions by default).
%               PiLow must be a number in the interval (0 1]. Default value
%               0.
%         int : mean vectors are simulated uniformly on a hypercube with
%               sides specified by int = [lower.bound, upper.bound].
%               The default value of int is [0 1]
%        resN : maximum number of mixture resimulations to find a
%               similation setting with prespecified level of overlapping.
%               The default value of resN is 100
%         tol : vector of length 2.
%               tol(1) (which will be called tolmap) specifies
%               the tolerance between the requested and empirical
%               misclassification probabilities (default is 1e-06)
%               tol(2) (which will be called tolnxc2) specifies the
%               tolerance to use in routine ncx2mixtcdf (which computes cdf
%               of linear combinations of non central chi2 distributions).
%               The default value of tol(2) 1e-06
%         lim : maximum number of integration terms to use inside routine
%               ncx2mixtcdf. Default is 1e06.
%               REMARK: Optional parameters tolncx2=tol(2) and lim will be
%               used by function ncx2mixtcdf.m which computes the cdf of a
%               linear combination of non central chi2 r.v.. This is the
%               probability of misclassification
%     Display : Level of display.
%               'off' displays no output;
%               'notify' (default) displays output if requested
%               overlap cannot be reached in a particular simulation
%               'iter' displays output at each iteration of each
%               simulation
%      R_seed : scalar > 0 for the seed to be used to generate random numbers
%               in a R instance. This is used to check consistency of the
%               results obtained with the R package MixSim. See file
%               Connect_Matlab_with_R_HELP.m to know how to connect MATLAB
%               with R.  This option requires the installation of the
%               R-(D)COM Interface. Default is 0, i.e. random numbers are
%               generated by matlab.
%
%       Remark: The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%       Remark: If 'BarOmega', 'MaxOmega' and 'StdOmega' are not specified, 
%               the function generates a mixture solely based on
%               'MaxOmega'=0.15. If both BarOmega, StdOmega and MaxOmega are
%               empty values as follows
%               out=MixSim(3,4,'MaxOmega','','BarOmega','') the following
%               message appears on the screen Error: At least one overlap
%               characteristic between BarOmega, MaxOmega and StdOmega
%               should be specified...
%
%  Output:
%
%  The output consists of a structure 'out' containing the following fields:
%
%            out.Pi  : vector of length k containing mixing proportions.
%                       sum(out.Pi)=1
%            out.Mu  : k-by-v matrix consisting of components' mean vectors
%                      Each row of this matrix is a centroid of a group
%             out.S  : v-by-v-by-k array containing covariances for the k
%                      groups
%       out.OmegaMap : matrix of misclassification probabilities (k-by-k);
%                      OmegaMap(i,j) = w_{j|i} is the probability that X
%                      coming from the i-th component (group) is classified
%                      to the j-th component.
%       out.BarOmega : scalar. Value of average overlap.
%                      BarOmega is computed as
%                      (sum(sum(OmegaMap))-k)/(0.5*k(k-1))
%       out.MaxOmega : scalar. Value of maximum overlap. MaxOmega is the
%                       maximum of OmegaMap(i,j)+OmegaMap(j,i)
%                       (i ~= j)=1, 2, ..., k. In other words MaxOmega=
%                      OmegaMap(rcMax(1),rcMax(2))+OmegaMap(rcMax(2),rcMax(1))
%       out.StdOmega : scalar. Value of standard deviation (std) of overlap.
%                      StdOmega is the standard deviation of k*(k-1)/2
%                      probabilities of overlapping
%         out.rcMax  : vector of length 2. It containes the row and column
%                      numbers associated with  the pair of components
%                      producing maximum overlap 'MaxOmega'
%              fail  : scalar, flag value. 0 represents successful mixture
%                      generation, 1 represents failure.
%
% See also tkmeans, tclust, tclustreg, lga, rlga, ncx2mixtcdf, restreigen
%
% References:
%
%   Maitra, R. and Melnykov, V. (2010) Simulating data to study performance
%   of finite mixture modeling and clustering algorithms, The Journal of
%   Computational and Graphical Statistics, 2:19, 354-376. (to refer to
%   this publication we will use "MM2010 JCGS")
%
%   Melnykov, V., Chen, W.-C., and Maitra, R. (2012) MixSim: An R Package
%   for Simulating Data to Study Performance of Clustering Algorithms,
%   Journal of Statistical Software, 51:12, 1-25.
%
%   Davies, R. (1980) The distribution of a linear combination of
%   chi-square random variables, Applied Statistics, 29, 323-333.
%
%   Garcia-Escudero, L.A.; Gordaliza, A.; Matran, C. and Mayo-Iscar, A.
%   (2008), A General Trimming Approach to Robust Cluster Analysis. Annals
%   of Statistics, Vol.36, 1324-1345. Technical Report available at
%   www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf
%
%   Reference below documents the problem of the ill-conditioning of the
%   eigenvalue-eigenvector computation.
%
%   Numerische Mathematik, 19. August 1969, Volume 13, Issue 4, pp 293-304
%   Balancing a matrix for calculation of eigenvalues and eigenvectors
%   Dr. B. N. Parlett, Dr. C. Reinsch
%
%  Parlett, B. N. and C. Reinsch, Balancing a Matrix for Calculation of
%  Eigenvalues and Eigenvectors, Handbook for Auto. Comp., Vol. II, Linear
%  Algebra, 1971,pp. 315-326.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mixsim')">Link to the help function</a>
% Last modified 06-Feb-2015
%

% Examples:
%
%{
	% Generate 3 groups in 4 dimensions using maximum overlap equal to 0.15
    rng(10,'twister')
    out=MixSim(3,4)
    n=200;
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S);
    spmplot(X,id)
%}
%
%{
    % Generate 4 groups in 5 dimensions using average overlap of 0.05 and
    % maximum overlap equal to 0.15
    k=4;
    v=5;
    BarOmega=0.05;
    out=MixSim(4,5,'BarOmega',BarOmega, 'MaxOmega',0.15)

	% Check a posteriori the average overlap
    disp('Posterior average overlap')
    disp((sum(sum(out.OmegaMap))-k)/(0.5*k*(k-1)))
    
    % Check a posteriori the maximum overlap
    % Extract elements above the diagonal and sum them with the transpose
    % of the elements below the diagonal. The maximum of all these numbers
    % must be very close to the required maximum overlap
    cand=triu(out.OmegaMap,1)+(tril(out.OmegaMap,-1))'
    disp('Posterior average overlap')
    max(cand(:))
%}

%{
    % Example of use of optional input option restrfactor. In the first case
    % restrfactor is 1.1 and the clusters are roughly homogeneous. In the
    % second case no constraint is imposed on the ratio of maximum and
    % minimum eigevalue among clusters so elliptical shape clusters are
    % allowed. In both cases the same random seed together with the same level
    % of average and maximum overlapping is used
    state1=2;
    randn('state', state1);
    rand('state', state1);
    out=MixSim(3,5,'BarOmega',0.1, 'MaxOmega',0.2, 'restrfactor',1.1);
    state1=2;
    randn('state', state1);
    rand('state', state1);
    out1=MixSim(3,5,'BarOmega',0.1, 'MaxOmega',0.2);

    n=200;
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S);
    spmplot(X,id,[],'box');
    set(gcf,'Name','restrfactor=1.2: almost homogeneous groups')
    title('restrfactor=1.2: almost homogeneous groups','fontsize',18);

    [X1,id1]=simdataset(n, out1.Pi, out1.Mu, out1.S);
    figure;
    spmplot(X1,id1,[],'box')
    set(gcf,'Name','Heterogeneous groups')
    title('Heterogeneous groups','fontsize',18)
    cascade
%}

%{
    % Control of average and standard deviation of overlap. Given an
    % average value of overlap, we explore the differences between imposing a
    % small or a large value of standard deviation of overlap.
    clc
    close all
    rng(10,'twister')
    k=4;
    v=5;
    n=200;
    BarOmega=0.10;
    StdOmega=0.15;
    out=MixSim(k,v,'BarOmega',BarOmega, 'StdOmega',StdOmega,'resN',10, 'Display', 'iter');
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S);

    rng(10,'twister')
    StdOmega1=0.05;
    out1=MixSim(k,v,'BarOmega',BarOmega, 'StdOmega',StdOmega1,'resN',10, 'Display', 'iter');
    [X1,id1]=simdataset(n, out1.Pi, out1.Mu, out1.S);
    disp('Comparison using OmegaMap')
    disp('When StdOmega is large in this example groups 3 are 4 do show a strong overlap')
    disp('When StdOmega is large in this example groups 1, 2, 3 are quite separate')
    disp(out.OmegaMap)
    disp('When StdOmega is small the probabilities of overlapping are much more similar')
    disp(out1.OmegaMap)

    disp('Comparison using interactive scatter plot matrices')
    spmplot(X,id,[],'box');
    set(gcf,'name',['BarOmega=' num2str(BarOmega) ' StdOmega=' num2str(StdOmega)])
    title(['BarOmega=' num2str(BarOmega) ' StdOmega=' num2str(StdOmega)])
    figure
    spmplot(X1,id1,[],'box');
    set(gcf,'name',['BarOmega=' num2str(BarOmega) ' StdOmega=' num2str(StdOmega1)])
    title(['BarOmega=' num2str(BarOmega) ' StdOmega=' num2str(StdOmega1)])
    cascade
%}

%% User options

% Default
if nargin<2;
    error('FSDA:MixSim:Missingv','k=number of components and v = number of dimensions must be specified');
end

if (v < 1)
    error('FSDA:MixSim:Wrongv','Wrong number of dimensions v')
end

if k<=1
        error('FSDA:MixSim:Wrongk','Wrong number of mixture components k')
end

Rseeddef = 0;
BarOmegadef = '';
MaxOmegadef = '';
StdOmegadef = '';
eccdef      = 0.9;
PiLowdef    = 0;
intdef      = [0 1];
resNdef     = 100;
toldef      = [1e-06; 1e-06];
limdef      = 1e06;
restrfactordef='';

options=struct('R_seed', Rseeddef, 'BarOmega',BarOmegadef,'MaxOmega',MaxOmegadef,...
    'StdOmega',StdOmegadef, 'Display', 'notify', ...
    'sph',false,'hom',false,'ecc',eccdef,'PiLow',PiLowdef,...
    'int',intdef,'resN',resNdef,'tol',toldef,'lim',limdef,'restrfactor',restrfactordef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:MixSim:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present in
    % structure options
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:MixSim:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin > 2
    
    % If the user inside user options has only specified BarOmega but not
    % MaxOmega then MaxOmega is initialied with an empty value
    checkBarOmega = strcmp(UserOptions,'BarOmega');
    checkMaxOmega = strcmp(UserOptions,'MaxOmega');
    
    if sum(checkBarOmega) && ~sum(checkMaxOmega);
        options.MaxOmega='';
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
    
end

% Default values for the optional parameters are set inside structure
% 'options'
int=options.int;

rcMax  = [0 0];
Lbound = int(1);
Ubound = int(2);

R_seed   = options.R_seed;
MaxOmega = options.MaxOmega;
BarOmega = options.BarOmega;
StdOmega = options.StdOmega;
tol=options.tol;
% tolmap = tolerance between the requested and empirical misclassification
% probabilities
tolmap    = options.tol(1);
% tolncx2 = tolerance to use in routine ncx2mixtcdf (which computes cdf of
% linear combinations of non central chi2 distributions)
tolncx2  = options.tol(2);
sph      = options.sph;
hom      = options.hom;
ecc      = options.ecc;
PiLow    = options.PiLow;
resN     = options.resN;
lim      = options.lim;
restrfactor= options.restrfactor;
Display  = options.Display;

if R_seed > 0
    
    % Check if a connection exists
    global R_lInK_hANdle %#ok<TLEV>
    if isempty(R_lInK_hANdle)
        examp=which('Connect_Matlab_with_R_HELP.m');
        % examp1=strrep(examp,'\','\\');
        examp1=strrep(examp,filesep,[filesep filesep]);
        
        disp('To run MixSim independently from R, option R_seed must be 0');
        disp('To ensure replicability of R examples contained in file demoMixSim.R');
        disp('i.e. to use R random number generators, first run openR');
        disp(['See instructions in file <a href="matlab:opentoline(''',examp1,''',5)">Connect_Matlab_with_R_HELP.m</a>']);
        error('FSDA:MixSim:NoRConn','--------------------------');
    end
    
    setseed = ['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')'];
    [~] = evalR(setseed);
end

if ~islogical(sph)
    error('FSDA:MixSim:Wrongsph','option sph must be a logical value')
end

if ~islogical(hom)
    error('FSDA:MixSim:Wronghom','option hom must be a logical value')
end

if ecc <= 0 || ecc > 1
    error('FSDA:MixSim:Wrongecc','ecc must be a scalar in the interval (0 1]')
end

if PiLow < 0 || PiLow > 1
    error('FSDA:MixSim:WrongPiLow','Option PiLow must be in interval [0 1]')
end

if int(1) >= int(2)
    error('FSDA:MixSim:Wrongint','Second element of int must be greater than first, that is: int(2) >int(1)')
end

if resN < 1
    error('FSDA:MixSim:WrongresN','Number of resimulations cannot be smaller than 1')
end

if (min(options.tol) <= 0)
    error('FSDA:MixSim:Wrongtol','Wrong value of tolerance, it must a scalar stricty greater than 0')
end

if lim < 1
    error('FSDA:MixSim:Wronglim','Wrong value of lim, it cannot be smaller than 1')
end

if isempty(MaxOmega) && isempty(StdOmega)  && ~isempty(BarOmega)
    % method =0 ==> just BarOmega has been specified
    method = 0;
    Omega = BarOmega;
elseif isempty(BarOmega) && isempty(StdOmega)
    % method =1 ==> just MaxOmega has been specified
    method = 1;
    if isempty(MaxOmega)
        Omega=0.15;
    else
        Omega = MaxOmega;
    end
elseif isempty(BarOmega) && isempty(MaxOmega) && ~isempty(StdOmega)
    % method =1.5 ==> Just StdOmega has been specified
    method = 1.5;
    Omega=StdOmega;
elseif ~isempty(BarOmega) && ~isempty(MaxOmega) && isempty(StdOmega)
    % method =2 ==> both BarOmega and MaxOmega have been specified
    method = 2;
elseif  ~isempty(BarOmega) && ~isempty(StdOmega) && isempty(MaxOmega)
    % method =3 ==> both BarOmega and StdOmega have been specified
    method = 3;
elseif isempty(BarOmega) && isempty(MaxOmega)
    % method =-1 ==> both BarOmega and MaxOmega have not been specified
    method = -1;
elseif ~isempty(BarOmega) && ~isempty(StdOmega) && ~isempty(MaxOmega)
    % method =4 ==> both BarOmega MaxOmega and StdOmega have been specified
    method = 4;
    
else
    method= 2;
end

% Get in vector indabovediag the linear indices of the elements above
% diagonal in a matrix of size k-by-k. This will be necessary to compute the
% starndard deviation of overlapping
indabovediag=triu2vec(k,1);

if method == 0 || method == 1 || method == 1.5
    emax=ecc;
    Q = OmegaClust(Omega, method, v, k, PiLow, Lbound, Ubound, ...
        emax, tol, lim, resN, sph, hom, restrfactor, Display);
    
elseif method == 2
    emax=ecc;
    
    % In this case both BarOmega and MaxOmega have been specified
    Q = OmegaBarOmegaMax(v, k, PiLow, Lbound, Ubound, ...
        emax, tol, lim, resN, sph, hom, BarOmega, MaxOmega, restrfactor, Display);
    
elseif method ==3
    % In this case both BarOmega and StdOmega have been specified
    emax=ecc;
    StdOmega=options.StdOmega;
    Q=OmegaBarOmegaStd(v, k, PiLow, Lbound, Ubound, ...
        emax, tol, lim, resN, sph, hom, BarOmega, StdOmega, restrfactor, Display);
elseif method ==4
    % In this case both MaxOmega, BarOmega and StdOmega have been specified
    error('FSDA:MixSim:TooManyConstr','It is not possible to specify both MaxOmega, BarOmega and StdOmega at the same time')
    
elseif method~=-1
    error('FSDA:MixSim:WrongMethod','Should never enter here')
else
    % isempty(BarOmega) && isempty(MaxOmega)
    error('FSDA:MixSim:ConstrRequired','At least one overlap characteristic between MaxOmega and BarOmega should be specified')
end

out = Q;

%% Beginning of inner functions

    function  Q = OmegaClust(Omega, method, v, k, PiLow, Lbound, Ubound, ...
            emax, tol, lim, resN, sph, hom, restrfactor, Display)
        % OmegaClust = procedure when average or maximum overlap or Std of
        % overlap is specified (not more than one overlapping measure)
        %
        %  INPUT parameters
        %
        % Omega     : scalar containing requested overlap value
        % method    : scalar which specifies whether average or maximum
        %             overlap is requested. 
        %             If method == 0 average overlap is requested 
        %             If method == 1 max overlap is requested
        %             If method == 1.5 std of overlap is requested
        % v         : scalar, dimensionality (number of variables)
        % k         : scalar, number of components (groups)
        % PiLow     : smallest mixing proportion allowed
        % Lbound    : lower bound for uniform hypercube at which mean vectors are
        %             simulated
        % Ubound    : upper bound for uniform hypercube at which mean vectors are
        %             simulated
        % emax      : maximum eccentricity
        % tol, lim  : parameters for ncx2mixtcdf.m which computes the probability
        %             of overlapping
        % resN      : scalar, number of resamplings allowed
        % sph       : scalar. If sph =1 spherical covariance matrices are
        %             generated, else covariance matrices with a
        %             prespecified maximum eccentricity are generated
        % hom       : scalar. If hom =1 equal covariance matrices are
        %             generated
        %
        % Optional input parameters
        %
        %restrfactor: scalar in the interval [1 \infty] which specifies the
        %             maximum ratio to allow between the largest eigenvalue and
        %             the smallest eigenvalue of the k covariance matrices which
        %             are generated.
        %    Display: Level of display.
        %             'off' displays no output;
        %             'notify' (default) displays output if requested
        %             overlap cannot be reached in a particular simulation
        %             'iter' displays output at each iteration of each
        %             simulation
        %
        %
        %  OUTPUT parameters
        %
        %   A structure Q containing the following fields
        %           Pi   : vector of length k containing mixing proportions
        %           Mu   : matrix of size k-by-v contains mean vectors of
        %                  the k groups
        %           S    : array of size v-by-v-by-k containing covariance matrices
        %       OmegaMap : k-by-k matrix containing misclassification probabilities
        %       BarOmega : scalar, average overlap of the groups which have been
        %                  generated
        %       MaxOmega : scalar, maximum overlap of the groups which have been
        %                  generated
        %       StdOmega : scalar, standard deviation of overlap of the groups which have been
        %                  generated
        %        rcMax   : vector of length 2 containing the pair of
        %                  components producing the highest overlap
        %           fail : flag indicating if the process failed (1). If everything went
        %                  well fail=0
        %
        
        
        if nargin< 15
            prnt = 1;
        else
            switch Display
                case {'none','off'}
                    prnt = 0;
                case {'notify','notify-detailed'}
                    prnt = 1;
                case {'iter','iter-detailed'}
                    prnt=2;
                otherwise
                    prnt = 1;
            end
        end
        
        fixcl=zeros(k,1);
        
        for isamp=1:resN
            
            fail=0;
            
            % Generate parameters
            % procedure genPi generates (mixture proportions) k numbers
            % whose sum is equal to 1 and the smallest value is not smaller
            % than PiLow
            Pigen=genPi(k,PiLow);
            % procedure genMu generates random centroids
            Mugen=genMu(v, k, Lbound, Ubound);
            
            % The last input parameter of genSigmaEcc and genSphSigma is a
            % boolean which specifies whether the matrices are equal (1) or
            % not (0)
            
            % Generate the covariance matrices
            if sph == 0
                % genSigmaEcc generates the covariance matrices with a
                % prespecified level of eccentricity
                Sgen=genSigmaEcc(v, k, emax, hom);
            else
                % genSphSigma generates spherical covariance matrices
                Sgen=genSphSigma(v, k, hom);
            end
            
            if nargin>13 && ~isempty(restrfactor)
                U=zeros(v,v,k);
                S05=Sgen;
                Sinv=Sgen;
                detS=zeros(k,1);
                Lambda_vk=zeros(v,k);
                for j=1:k
                    [Uj,Lambdaj] = eig(Sgen(:,:,j));
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_vk(:,j)=diag(Lambdaj);
                end
                
                % The line below should be unnecessary
                Lambda_vk(Lambda_vk<0)=0;
                
                % Apply the restrictions to matrix Lambda_vk
                autovalues=restreigen(Lambda_vk,Pigen,restrfactor);
                
                for j=1:k
                    %disp(j)
                    Sgen(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
                    S05(:,:,j) = U(:,:,j)*diag(sqrt(autovalues(:,j)))* (U(:,:,j)');
                    Sinv(:,:,j) = U(:,:,j)*diag(autovalues(:,j).^-1)* (U(:,:,j)');
                    detS(j)=prod(autovalues(:,j));
                    % Alternative code: in principle more efficient but slower
                    % because diag is a built in function
                    % sigmaini(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
                end
                [li, di, const1]=ComputePars(v, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
            else
                % prepare parameters:  row 953 of file libOverlap.c
                [li, di, const1]=ComputePars(v, k, Pigen, Mugen, Sgen);
            end
            
            % check if desired overlap (maximum, average or std) is reachable
            asympt = 1;
            c = 0.0;
            
            [OmegaMap, Balpha, Malpha, rcMax] = ...
                GetOmegaMap(c, v, k, li, di, const1, fixcl, tolncx2, lim, asympt);
            
            if method == 0
                diff = Balpha - Omega;
                if prnt>1
                    disp(['Average empirical overlap - Average requested overlap=' num2str(diff)])
                end
            elseif method==1
                diff = Malpha - Omega;
                if prnt>1
                    disp(['Max empirical overlap - Max requested overlap=' num2str(diff)])
                end
            elseif method==1.5
                % Compute StdOverlap
                % Extract elements above diagonal and compute Stdoverlap
                cand=triu(OmegaMap,1)+(tril(OmegaMap,-1))';
                overlapv=cand(:);
                %                 overlapv=overlapv(overlapv>0);
                %                 if length(overlapv)<combk2;
                %                     overlapc=[overlapv; zeros(combk2-length(overlapv),1)];
                %                 else
                %                     overlapc=overlapv;
                %                 end
                %                 % Compute standard deviation of overlap for current
                %                 % solution
                %                 Stdalpha=std(overlapc);
                
                Stdalpha=std(overlapv(indabovediag));
                
                diff = Stdalpha - Omega;
                if prnt>1
                    disp(['Std of empirical overlap - Std requested overlap=' num2str(diff)])
                end
            else
                
            end
            
            if (diff < -tolmap) % Prefixed overlapping is not reachable
                if prnt>=1
                    disp(['Warning: the desired overlap is too high and cannot be reached in simulation '  num2str(isamp)]);
                end
                fail = 1;
            else
                lower=0;
                upper=4;
                % c is a constant which is used to scale the covariance
                % matrices of the groups in order to obtain a prespecified
                % level of average or maximum overlap
                c=0;
                while c==0
                    if method==0 || method ==1
                        [c, OmegaMap, Balpha, Malpha] = FindC(lower, upper, Omega, ...
                            method, v, k, li, di, const1, fixcl, tol, lim);
                    else
                        [c, OmegaMap, Balpha, Malpha, Stdalpha] = FindCStd(lower, upper, Omega, ...
                            v, k, li, di, const1, fixcl, tol, lim);
                    end
                    lower =upper;
                    upper=upper^2;
                    if upper>100000
                        if prnt>=1
                            disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                        end
                        fail=1;
                        break
                    end
                end
            end
            
            if fail==0
                Sgen=c*Sgen;
                break  % this break enables to get out from the resampling loop
            end
        end
        
        if method == 0
            if prnt>1
                disp(['Average empirical overlap - Average requested overlap=' num2str(Balpha-Omega)])
            end
        elseif method ==1
            if prnt>1
                disp(['Max empirical overlap - Max requested overlap=' num2str(Malpha-Omega)])
            end
        elseif method ==1.5
            if prnt>1
                disp(['Std empirical overlap - Std requested overlap=' num2str(Stdalpha-Omega)])
            end
        else
        end
        
        
        if isamp == resN
            warning('FSDA:MixSim:OverlapNotReached',['The desired overlap has not been reached in ' num2str(resN) ' simulations']);
            warning('FSDA:MixSim:NsimulTooSmall','Please increase the number of simulations allowed (option resN) or change the value of overlap');
            fail = 1;
        end
        
        % Compute standard deviation of overlap
        cand=triu(OmegaMap,1)+(tril(OmegaMap,-1))';
        overlapv=cand(:);
        
        
        % Recompute rcMax (the indices of the two groups producing the
        % highest overlap)
        [~,indmaxoverlapv]=max(overlapv);
        [rcMax(1), rcMax(2)]=ind2sub([k k],indmaxoverlapv);
        
        %         overlapv=overlapv(overlapv>0);
        %         if length(overlapv)<0.5*k*(k-1);
        %             overlapc=[overlapv; zeros(0.5*k*(k-1)-length(overlapv),1)];
        %         else
        %             overlapc=overlapv;
        %         end
        %         % Compute standard deviation of overlap for current
        %         % solution
        %         stdoverlap=std(overlapc);
        stdoverlap=std(overlapv(indabovediag));
        
        
        
        Q = struct;
        Q.OmegaMap=OmegaMap;
        Q.BarOmega=Balpha;
        Q.MaxOmega=Malpha;
        Q.StdOmega=stdoverlap;
        Q.fail=fail;
        Q.Pi=Pigen;
        Q.Mu=Mugen;
        Q.S=Sgen;
        Q.rcMax=rcMax;
        
    end




    function Q=OmegaBarOmegaMax(p, k, PiLow, Lbound, Ubound, ...
            emax, tol, lim, resN, sph, hom, BarOmega, MaxOmega, restrfactor, Display)
        % OmegaBarOmegaMax = procedure when average and maximum overlap are both specified
        %
        %
        %  INPUT parameters
        %
        %       p  : scalar, dimensionality
        %       k  : scalar, number of components
        %    PiLow : scalar, smallest mixing proportion allowed
        %   Lbound : scalar, lower bound for uniform hypercube at which mean vectors at simulated
        %   Ubound : scalar, upper bound for uniform hypercube at which mean vectors at simulated
        %     emax : scalar, maximum eccentricity
        % tol, lim : parameters for ncx2mixtcdf function which computes the
        %            probabilities of overlapping
        %     resN : scalar, number of resamplings allowed
        %      sph : scalar, if sph =1 spherical covariance matrices are
        %            generated
        %      hom : scalar. If hom =1 equal covariance matrices are
        %            generated. REMARK: in MM2010 JCGS this routine is
        %            always called using hom=0
        % BarOmega : scalar, required average overlap
        % MaxOmega : scalar, required maximum overlap
        %
        % Optional input parameters
        %
        %restrfactor: scalar in the interval [1 \infty] which specifies the
        %             maximum ratio to allow between the largest eigenvalue and
        %             the smallest eigenvalue of the k covariance matrices which
        %             are generated
        %    Display: Level of display.
        %             'off' displays no output;
        %             'notify' (default) displays output if requested
        %             overlap cannot be reached in a particular simulation
        %             'iter' displays output at each iteration of each
        %             simulation
        %
        %  OUTPUT parameters
        %
        %   A structure Q containing the following fields
        %
        %             Pi : mixing proportions
        %           Mu   : matrix of size k-by-v containing mean vectors of
        %                  the k groups
        %           S    : array of size v-by-v-by-k containing covariance matrices
        %       OmegaMap : k-by-k matrix containing misclassification probabilities
        %       BarOmega : scalar, average overlap of the groups which have been
        %                  generated
        %       MaxOmega : scalar, maximum overlap of the groups which have been
        %                  generated
        %       StdOmega : scalar, standard deviation of overlap of the groups
        %                  which have been generated
        %        rcMax   : vector of length 2 containing the pair of
        %                  components producing the highest overlap
        %           fail : flag indicating if the process failed (1). If
        %                  everything went well fail=0
        Balpha=BarOmega;
        Malpha=MaxOmega;
        if Malpha<Balpha || Malpha>Balpha*k*(k-1)/2
            disp('Both conditions should hold:')
            disp('1. MaxOverlap > AverOverlap')
            disp('2.  MaxOverlap < AverOverlap * K (K - 1) / 2')
            error('FSDA:MixSim:WrongOverlapSupplied','Incorrect values of average and maximum overlaps...');
            
        else
            
            if nargin< 15
                prnt = 1;
            else
                switch Display
                    case {'none','off'}
                        prnt = 0;
                    case {'notify','notify-detailed'}
                        prnt = 1;
                    case {'iter','iter-detailed'}
                        prnt=2;
                    otherwise
                        prnt = 1;
                end
            end
            li2=zeros(2, 2, p);
            di2=li2;
            const12=zeros(2);
            
            fix2=zeros(2,1);
            
            for isamp=1:resN
                
                % generate parameters
                % procedure genPi generates (mixture proportions) k numbers
                % whose sum is equal to 1 and the smallest value is not
                % smaller than PiLow
                Pigen=genPi(k,PiLow);
                
                % procedure genMu generates random centroids
                Mugen=genMu(p, k, Lbound, Ubound);
                
                % The last input parameter of genSigmaEcc and genSphSigma is a
                % boolean which specifies whether the matrices are equal
                % (1) or not (0)
                
                % Generate the covariance matrices
                if sph == 0
                    % genSigmaEcc generates the covariance matrices with a
                    % prespecified level of eccentricity
                    % in MM2010 JCGS the procedure is always called using hom=0
                    Sgen=genSigmaEcc(p, k, emax, hom);
                else
                    % genSphSigma generates spherical covariance matrices
                    % in MM2010 JCGS the procedure is always called using hom=0
                    Sgen=genSphSigma(p, k, hom);
                end
                
                if nargin>13 && ~isempty(restrfactor)
                    U=zeros(v,v,k);
                    S05=Sgen;
                    Sinv=Sgen;
                    detS=zeros(k,1);
                    Lambda_vk=zeros(v,k);
                    for j=1:k
                        [Uj,Lambdaj] = eig(Sgen(:,:,j));
                        % Store eigenvectors and eigenvalues of group j
                        U(:,:,j)=Uj;
                        Lambda_vk(:,j)=diag(Lambdaj);
                    end
                    
                    % The line below should be unnecessary
                    Lambda_vk(Lambda_vk<0)=0;
                    
                    % Apply the restrictions to matrix Lambda_vk
                    autovalues=restreigen(Lambda_vk,Pigen,restrfactor);
                    
                    for j=1:k
                        %disp(j)
                        Sgen(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
                        S05(:,:,j) = U(:,:,j)*diag(sqrt(autovalues(:,j)))* (U(:,:,j)');
                        Sinv(:,:,j) = U(:,:,j)*diag(autovalues(:,j).^-1)* (U(:,:,j)');
                        detS(j)=prod(autovalues(:,j));
                        % Alternative code: in principle more efficient but slower
                        % because diag is a built in function
                        % sigmaini(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
                    end
                    % prepare parameters:  row 953 of file libOverlap.c
                    [li, di, const1]=ComputePars(p, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
                    
                else
                    % prepare parameters:  row 953 of file libOverlap.c
                    [li, di, const1]=ComputePars(p, k, Pigen, Mugen, Sgen);
                    
                end
                
                % Check if maximum overlap is reachable.  Maximum
                % overlapping takes place when c \rightarrow \infty that is
                % when asympt =1
                asympt = 1;
                % If asympt=1 c does not matter so any value of c will do
                c = 0;
                
                % Initialize fail. A flag which indicates
                % whether the required overlapped has been
                % reached. fail=1 means required overlap
                % not reached
                fail=1;
                
                % fixcl = vector which specifies what are the clusters which
                % participate to the process of inflation. In this stage
                % all clusters partecipate to the process of inflation, so
                % fixcl is a vector of zeroes
                fixcl=zeros(k,1);
                
                [OmegaMap, Balpha, Malpha, rcMax]=GetOmegaMap(c, p, k, li, di, const1, fixcl, tolncx2, lim, asympt);
                
                % Malpha is the maximum level of overlapping which can be
                % reached with asympt=1
                diff = Malpha - MaxOmega;
                
                
                if diff >= -tolmap %  maximum level of overlapping is reachable
                    
                    lower = 0.0;
                    upper = 2^10;
                    
                    while fail ~=0
                        
                        % find constant c for two clusters which show the
                        % highest ovelapping
                        
                        % Extract parameters for the two clusters with the
                        % highest overlapping
                        li2(1,2,:) = li(rcMax(1),rcMax(2),:);
                        di2(1,2,:) = di(rcMax(1),rcMax(2),:);
                        const12(1,2)=const1(rcMax(1),rcMax(2));
                        li2(2,1,:) = li(rcMax(2),rcMax(1),:);
                        di2(2,1,:) = di(rcMax(2),rcMax(1),:);
                        const12(2,1)=const1(rcMax(2),rcMax(1));
                        
                        Malpha = MaxOmega;
                        
                        % the fourth input element of FindC is method (in
                        % this case 1 is supplied because maximum overlap
                        % is requested)
                        % The sixth input element of findC is k. In this
                        % case k=2 (the two groups which show the highest
                        % overlap)
                        c=FindC(lower, upper, Malpha, 1, p, 2, li2, di2, const12, fix2, tol, lim);
                        
                        if c == 0 % abnormal termination
                            if prnt>=1
                                disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                            end
                            fail = 1;
                            break
                        end
                        
                        asympt = 0;
                        % Compute map of misclassification probabilities
                        % (OmegaMap)
                        % using the value of c which comes out from procedure
                        % findC
                        % Balpha = Average overlap using c from findC
                        % Malpha = Maximum overlap using c from findC
                        [OmegaMap, Balpha, Malpha, rcMax]=GetOmegaMap(c, p, k, li, di, const1, fixcl, tolncx2, lim, asympt);
                        upper = c;
                        
                        % We hope that Balpha (overall average overlapping)
                        % obtained using c associated with the two clusters
                        % which showed the highest overlapping is greater
                        % than BarOmega (average requested overlapping).
                        diff = Balpha - BarOmega;
                        
                        % If diff<tolmap the desired average overlap characteristic is
                        % possibly unattainable using the candidate \mu
                        % \Sigma and \pi and a new candidate is needed
                        if (diff < -tolmap) % BarOmega is not reachable
                            if prnt>=1
                                disp(['Warning: the desired average overlap cannot be reached in simulation '  num2str(isamp)]);
                            end
                            fail = 1;
                            break
                        end
                        
                        % Now we make sure that none of pairwise overlaps
                        % (that is  make sure that the maximum pairwise
                        % overlap (which is Malpha) does not exceed MaxOmega
                        % (the maximum requested overlap). If this is the
                        % case  do another iteration of the loop (while
                        % fail ~=0) using rcMax which has just been found.
                        
                        
                        
                        diff = Malpha - MaxOmega;
                        if prnt>1
                            disp(['Average empirical overlap - Average requested overlap=' num2str(Balpha - BarOmega)])
                            disp(['Max empirical overlap - Max requested overlap=' num2str(diff)])
                        end
                        
                        if (diff < tolmap) %  MaxOmega has been reached
                            fail = 0;
                            break
                        end
                        
                    end
                    
                end
                
                if fail == 0
                    %  OmegaMax is reached and OmegaBar is reachable
                    %  correct covariances by multiplier c
                    
                    if nargin>13 && ~isempty(restrfactor)
                        S05=(c^0.5)*S05;
                        Sinv=(1/c)*Sinv;
                        detS=(c^v)*detS;
                        [li, di, const1]=ComputePars(p, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
                    else
                        Sgen=c*Sgen;
                        [li, di, const1]=ComputePars(p, k, Pigen, Mugen, Sgen);
                    end
                    
                    % The two clusters which enabled to obtain the highest
                    % overlap are kept unchanged all the way through the
                    % termination of the algorithm
                    fixcl(rcMax(1)) = 1;
                    fixcl(rcMax(2)) = 1;
                    upper = 1;
                    
                    % Now find the c which guarrantees the average
                    % requested ovelapping BarOmega
                    % method =0 because average overlapping is requested
                    method = 0;
                    Malphain=Malpha;
                    % FindC(lower, upper, Balpha, method, p, K, li, di, const1, fix, pars, lim, &c, OmegaMap, &Balpha, &Malpha, rcMax);
                    [c, OmegaMap, Balpha, Malpha, rcMax]=FindC(lower, upper, BarOmega, method, p, k, li, di, const1, fixcl, tol, lim);
                    
                    % If c =0 max number of iterations has been reached
                    % inside findc therefore another simulation is
                    % requested
                    if c==0 || abs(Malphain-Malpha)>tolmap
                        if prnt>=1
                            disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                        end
                        fail=1;
                    else
                        % correct covariances by multiplier c
                        for jj=1:k
                            if fixcl(jj) == 0
                                Sgen(:,:,jj) = c * Sgen(:,:,jj);
                            end
                        end
                        
                        
                        % Inside if  fail==0 OmegaMax is reached and OmegaBar is reachable
                        break
                    end
                    
                end
                
            end
            
            if prnt>1
                disp(['Average empirical overlap - Average requested overlap=' num2str(Balpha - BarOmega)])
                disp(['Max empirical overlap - Max requested overlap=' num2str(Malpha - MaxOmega)])
            end
            
            if isamp == resN && resN>1
                warning('FSDA:MixSim:OverlapNotReached',['The desired overlap has not been reached in ' num2str(resN) ' simulations']);
                warning('FSDA:MixSim:NsimulTooSmall','Increase the number of simulations allowed (option resN) or change the value of overlap');
                
                fail = 1;
                
            end
            
            % Compute standard deviation of overlap
            cand=triu(OmegaMap,1)+(tril(OmegaMap,-1))';
            overlapv=cand(:);
            %             overlapv=overlapv(overlapv>0);
            %             if length(overlapv)<0.5*k*(k-1);
            %                 overlapc=[overlapv; zeros(0.5*k*(k-1)-length(overlapv),1)];
            %             else
            %                 overlapc=overlapv;
            %             end
            %             % Compute standard deviation of overlap for current
            %             % solution
            %             stdoverlap=std(overlapc);
            
            stdoverlap=std(overlapv(indabovediag));
            
            Q=struct;
            Q.OmegaMap=OmegaMap;
            Q.BarOmega=Balpha;
            Q.MaxOmega=Malpha;
            Q.StdOmega=stdoverlap;
            Q.fail=fail;
            Q.Pi=Pigen;
            Q.Mu=Mugen;
            Q.S=Sgen;
            Q.rcMax=rcMax;
            
        end
    end


    function Q=OmegaBarOmegaStd(p, k, PiLow, Lbound, Ubound, ...
            emax, tol, lim, resN, sph, hom, BarOmega, StdOmega, restrfactor,Display)
        % OmegaBarOmegaStd = procedure when average and std of overlap are both specified
        %
        %
        %  INPUT parameters
        %
        %       p  : scalar, dimensionality
        %       k  : scalar, number of components
        %    PiLow : scalar, smallest mixing proportion allowed
        %   Lbound : scalar, lower bound for uniform hypercube at which mean vectors at simulated
        %   Ubound : scalar, upper bound for uniform hypercube at which mean vectors at simulated
        %     emax : scalar, maximum eccentricity
        % tol, lim : parameters for ncx2mixtcdf function which computes the
        %            probabilities of overlapping
        %     resN : scalar, number of resamplings allowed
        %      sph : scalar, if sph =1 spherical covariance matrices are
        %            generated
        %      hom : scalar. If hom =1 equal covariance matrices are
        %            generated. REMARK: in MM2010 JCGS this routine is
        %            always called using hom=0
        % BarOmega : scalar, required average overlap
        % StdOmega : scalar, required standard deviation of overlap
        %
        % Optional input parameters
        %
        %restrfactor: scalar in the interval [1 \infty] which specifies the
        %             maximum ratio to allow between the largest eigenvalue and
        %             the smallest eigenvalue of the k covariance matrices which
        %             are generated
        %    Display: Level of display.
        %             'off' displays no output;
        %             'notify' (default) displays output if requested
        %             overlap cannot be reached in a particular simulation
        %             'iter' displays output at each iteration of each
        %             simulation
        %
        %  OUTPUT parameters
        %
        %   A structure Q containing the following fields
        %
        %             Pi : mixing proportions
        %           Mu   : matrix of size k-by-v containing mean vectors of
        %                  the k groups
        %           S    : array of size v-by-v-by-k containing covariance matrices
        %       OmegaMap : k-by-k matrix containing misclassification probabilities
        %       BarOmega : scalar, average overlap of the groups which have been
        %                  generated
        %       MaxOmega : scalar, maximum overlap of the groups which have been
        %                  generated
        %        rcMax   : vector of length 2 containing the pair of
        %                  components producing the highest overlap
        %           fail : flag indicating if the process failed (1). If
        %                  everything went well fail=0
        
        if k<=2
            error('FSDA:MixSim:WrongOverlapSupplied','Average and std of overlap can be both set when k>2');
        end
        
        if nargin< 15
            prnt = 1;
        else
            switch Display
                case {'none','off'}
                    prnt = 0;
                case {'notify','notify-detailed'}
                    prnt = 1;
                case {'iter','iter-detailed'}
                    prnt=2;
                otherwise
                    prnt = 1;
            end
        end
        
        li2=zeros(2, 2, p);
        di2=li2;
        const12=zeros(2);
        
        fix2=zeros(2,1);
        
        MaxOmegaloopini=BarOmega*1.1;
        
        eps=tolmap*10;
        stdoverlap=NaN;
        
        for isamp=1:resN
            
            % generate parameters
            % procedure genPi generates (mixture proportions) k numbers
            % whose sum is equal to 1 and the smallest value is not
            % smaller than PiLow
            Pigen=genPi(k,PiLow);
            
            % procedure genMu generates random centroids
            Mugen=genMu(p, k, Lbound, Ubound);
            
            % The last input parameter of genSigmaEcc and genSphSigma is a
            % boolean which specifies whether the matrices are equal
            % (1) or not (0)
            
            % Generate the covariance matrices
            if sph == 0
                % genSigmaEcc generates the covariance matrices with a
                % prespecified level of eccentricity
                % in MM2010 JCGS the procedure is always called using hom=0
                Sgen=genSigmaEcc(p, k, emax, hom);
            else
                % genSphSigma generates spherical covariance matrices
                % in MM2010 JCGS the procedure is always called using hom=0
                Sgen=genSphSigma(p, k, hom);
            end
            
            if nargin>13 && ~isempty(restrfactor)
                U=zeros(v,v,k);
                S05=Sgen;
                Sinv=Sgen;
                detS=zeros(k,1);
                Lambda_vk=zeros(v,k);
                for j=1:k
                    [Uj,Lambdaj] = eig(Sgen(:,:,j));
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_vk(:,j)=diag(Lambdaj);
                end
                
                % The line below should be unnecessary
                Lambda_vk(Lambda_vk<0)=0;
                
                % Apply the restrictions to matrix Lambda_vk
                autovalues=restreigen(Lambda_vk,Pigen,restrfactor);
                
                for j=1:k
                    %disp(j)
                    Sgen(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
                    S05(:,:,j) = U(:,:,j)*diag(sqrt(autovalues(:,j)))* (U(:,:,j)');
                    Sinv(:,:,j) = U(:,:,j)*diag(autovalues(:,j).^-1)* (U(:,:,j)');
                    detS(j)=prod(autovalues(:,j));
                    % Alternative code: in principle more efficient but slower
                    % because diag is a built in function
                    % sigmaini(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
                end
                % prepare parameters:  row 953 of file libOverlap.c
                [liini, diini, const1ini]=ComputePars(p, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
                
            else
                % prepare parameters:  row 953 of file libOverlap.c
                [liini, diini, const1ini]=ComputePars(p, k, Pigen, Mugen, Sgen);
                
            end
            
            % Check if maximum overlap is reachable.  Maximum
            % overlapping takes place when c \rightarrow \infty that is
            % when asympt =1
            asympt = 1;
            % If asympt=1 c does not matter so any value of c will do
            c = 0;
            
            % fixcl = vector which specifies what are the clusters which
            % participate to the process of inflation. In this stage
            % all clusters partecipate to the process of inflation, so
            % fixcl is a vector of zeroes
            fixclini=zeros(k,1);
            
            [OmegaMap, Balphaini, Malphaini, rcMaxini]=GetOmegaMap(c, p, k, liini, diini, const1ini, fixclini, tolncx2, lim, asympt);
            
            % Malpha is the maximum level of overlapping which can be
            % reached with asympt=1
            
            
            % Check if sigmamax is reachable:
            % sigmamax= sqrt((M-xmin)*(xmax-M))
            % Note that here xmin=0
            sigmamax=sqrt(BarOmega*(Malphaini-BarOmega));
            
            
            diff = sigmamax - StdOmega;
            
            % Initialize fail
            fail=1;
            
            
            if diff < -tolmap % Requested StdOmega of overlapping is reachable
                disp('Requested sigma of ovelap must be smaller than')
                disp('BarOmega*(MaxOmega(achievable) - BarOmega)')
                disp(['In simulation ' num2str(isamp)])
                disp(['MaxOmega(achievable)=' num2str(Malphaini) ' and'])
                disp(['Maximum achievable sigma is=' num2str(sigmamax)]);
                Balpha=Balphaini;
                Malpha=Malphaini;
                
            else
                
                % Now we loop for different values of MaxOmega in order to
                % find the one which guarantees the required Std of
                % overlap
                if nargin>13 && ~isempty(restrfactor)
                    S05ini=S05;
                    Sinvini=Sinv;
                    detSini=detS;
                end
                Sgenini=Sgen;
                
                Erho1=10;
                step=0.03;
                step05=0;
                
                MaxOmegaloop=MaxOmegaloopini;
                iter=0;
                while abs(Erho1-1)>eps  
                    if step<1e-15
                        break
                    end
                    % if the value of MaxOmegaloop is greater than than the max
                    % overlap achivable (which is Malphaini) than requested
                    % standard deviation is too large and it is necessary to
                    % decrease it
                    if MaxOmegaloop> Malphaini
                        fail=1;
                        if prnt>=1
                            disp('Please decrease requested std of overlap')
                        end
                        break
                    end
                    
                    fail=1;
                    Sgen=Sgenini;
                    if nargin>13 && ~isempty(restrfactor)
                        S05=S05ini;
                        Sinv=Sinvini;
                        detS=detSini;
                    end
                    
                    Malpha=Malphaini;
                    Balpha=Balphaini;
                    % rcMaxini=[1;2];
                    rcMax=rcMaxini;
                    li=liini;
                    di=diini;
                    const1=const1ini;
                    fixcl=fixclini;
                    
                    lower = 0.0;
                    upper = 2^10;
                    
                    while fail ~=0
                        
                        % Extract parameters for the two clusters with the
                        % highest overlapping
                        li2(1,2,:) = li(rcMax(1),rcMax(2),:);
                        di2(1,2,:) = di(rcMax(1),rcMax(2),:);
                        const12(1,2)=const1(rcMax(1),rcMax(2));
                        li2(2,1,:) = li(rcMax(2),rcMax(1),:);
                        di2(2,1,:) = di(rcMax(2),rcMax(1),:);
                        const12(2,1)=const1(rcMax(2),rcMax(1));
                        
                        Malpha = MaxOmegaloop;
                        
                        % The fourth input element of FindC is method (in
                        % this case 1 is supplied because maximum overlap
                        % is requested)
                        % The sixth input element of findC is k. In this
                        % case k=2 (the two groups which show the highest
                        % overlap)
                        c=FindC(lower, upper, Malpha, 1, p, 2, li2, di2, const12, fix2, tol, lim);
                        
                        if c == 0 && prnt >=1 % abnormal termination
                            disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                            fail = 1;
                            break
                        end
                        
                        asympt = 0;
                        % Compute map of misclassification probabilities
                        % (OmegaMap)
                        % using the value of c which comes out from procedure
                        % findC
                        % Balpha = Average overlap using c from findC
                        % Malpha = Maximum overlap using c from findC
                        [OmegaMap, Balpha, Malpha, rcMax]=GetOmegaMap(c, p, k, li, di, const1, fixcl, tolncx2, lim, asympt);
                        upper = c;
                        
                        % We hope that Balpha (overall average overlapping)
                        % obtained using c associated with the two clusters
                        % which showed the highest overlapping is greater
                        % than BarOmega (average requested overlapping).
                        diff = Balpha - BarOmega;
                        % If diff<tolmap the desired average overlap characteristic is
                        % possibly unattainable using the candidate \mu
                        % \Sigma and \pi and a new candidate is needed
                        if (diff < -tolmap) % BarOmega is not reachable
                            if Erho1<1 && prnt>=1
                                disp(['Warning: sigma is too small in simulation '  num2str(isamp)]);
                            end
                            fail = 1;
                            break
                        end
                        
                        % Now we make sure that none of pairwise overlaps
                        % (that is  make sure that the maximum pairwise
                        % overlap (which is Malpha) does not exceed MaxOmega
                        % (the maximum requested overlap). If this is the
                        % case  do another iteration of the loop (while
                        % fail ~=0) using rcMax which has just been found.
                        
                        % TO CHECK REMOVED
                        % diff = Malpha - MaxOmegaloop;
                        %if (diff < tolmap) %  MaxOmega has been reached
                        fail = 0;
                        %  break
                        % end
                        
                    end
                    
                    
                    
                    if fail == 0
                        %  OmegaMax is reached and OmegaBar is reachable
                        %  correct covariances by multiplier C
                        
                        if nargin>13 && ~isempty(restrfactor)
                            S05=(c^0.5)*S05;
                            Sinv=(1/c)*Sinv;
                            detS=(c^v)*detS;
                            [li, di, const1]=ComputePars(p, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
                        else
                            Sgen=c*Sgen;
                            [li, di, const1]=ComputePars(p, k, Pigen, Mugen, Sgen);
                        end
                        
                        % The two clusters which enabled to obtain the highest
                        % overlap are kept unchanged all the way through the
                        % termination of the algorithm
                        fixcl(rcMax(1)) = 1;
                        fixcl(rcMax(2)) = 1;
                        upper = 1;
                        
                        % Now find the c which guarrantees the average
                        % requested ovelapping BarOmega
                        % method =0 because average overlapping is requested
                        method = 0;
                        % Malphain=Malpha;
                        % FindC(lower, upper, Balpha, method, p, K, li, di, const1, fix, pars, lim, &c, OmegaMap, &Balpha, &Malpha, rcMax);
                        [c, OmegaMap, Balpha, Malpha, rcMax]=FindC(lower, upper, BarOmega, method, p, k, li, di, const1, fixcl, tol, lim);
                        
                        % If c =0 max number of iterations has been reached
                        % inside findc therefore another simulation is
                        % requested
                        if c==0 % || abs(Malphain-Malpha)>1*tolmap
                            if Erho1<10 &&  prnt >=1
                                disp(['Warning: sigma is too large in simulation '  num2str(isamp)]);
                            end
                            fail=1;
                            break
                        else
                            % correct covariances by multiplier c
                            for jj=1:k
                                if fixcl(jj) == 0
                                    Sgen(:,:,jj) = c * Sgen(:,:,jj);
                                end
                            end
                            
                        end
                        
                    end
                    
                    % Compute standard deviation of overlap
                    cand=triu(OmegaMap,1)+(tril(OmegaMap,-1))';
                    overlapv=cand(:);
                    %                     overlapv=overlapv(overlapv>0);
                    %                     if length(overlapv)<combk2;
                    %                         overlapc=[overlapv; zeros(combk2-length(overlapv),1)];
                    %                     else
                    %                         overlapc=overlapv;
                    %                     end
                    %                     % Compute standard deviation of overlap for current
                    %                     % solution
                    %                     stdoverlap=std(overlapc);
                    stdoverlap=std(overlapv(indabovediag));
                    
                    Erho1old=Erho1;
                    Erho1=StdOmega/stdoverlap;
                    if prnt == 2
                        disp(['Iteration ' num2str(iter) ' in simulation ' num2str(isamp)])
                        disp(['Average overlap =' num2str(Balpha)])
                        disp('Ratio between std of required overlap and std of empirical overlap')
                        disp(Erho1)
                        disp(['Step=' num2str(step)])
                    end
                    
                    iter=iter+1;
                    if step05==1 || (Erho1old>1 && Erho1<1)
                        step=step*0.5;
                        step05=1;
                    else
                        step=step*0.95;
                    end
                    
                    if Erho1>1
                        MaxOmegaloop=MaxOmegaloop+step;
                    else
                        
                        if MaxOmegaloop-step<=BarOmega && iter >=5
                            if prnt >=1
                                disp(['Warning: sigma is too small in simulation '  num2str(isamp)]);
                            end
                            fail=1;
                            break
                        end
                        %MaxOmegaloop=max(MaxOmegaloop-step,MaxOmegaloopini);
                        MaxOmegaloop=max(MaxOmegaloop-step,BarOmega);
                        
                    end
                end
                if fail==0
                    %                     if prnt>=1
                    %                     disp(['average + sigma overlap not reached in simulation ' num2str(resN)])
                    %                     end
                    break
                end
            end
        end
        if isamp == resN && resN>1
                warning('FSDA:MixSim:OverlapNotReached',['The desired overlap has not been reached in ' num2str(resN) ' simulations']);
                warning('FSDA:MixSim:NsimulTooSmall','Increase the number of simulations allowed (option resN) or change the value of overlap');
            
            fail = 1;
            
        end
        
        
        
        Q=struct;
        Q.OmegaMap=OmegaMap;
        Q.BarOmega=Balpha;
        Q.MaxOmega=Malpha;
        Q.StdOmega=stdoverlap;
        Q.fail=fail;
        Q.Pi=Pigen;
        Q.Mu=Mugen;
        Q.S=Sgen;
        Q.rcMax=rcMax;
        
        
    end

    function Pigen=genPi(k,PiLow)
        % genPi generates vector of mixing proportions
        %
        %  Required input arguments:
        % 		k       : scalar, number of components
        % 		PiLow   : scalar, smallest possible mixing proportion
        %
        %  OUTPUT parameters
        %
        % 		Pigen : vector of length k containing mixing proportions
        %               Vector Pigen satisfies the following constraints:
        %               sum(Pigen)=1
        %               min(Pigen)>=PiLow
        %
        %  For example Pigen=genPi(4,0.24)
        %  produces Pigen=
        %    0.2440
        %    0.2517
        %    0.2565
        %    0.2478
        %
        flag = 0;
        
        if PiLow >= 1 || PiLow <= 0
            if PiLow < 0 || PiLow >= 1
                disp('Warning: PiLow is out of range... generated equal mixing proportions...');
            end
            Pigen=zeros(k,1);
            Pigen=Pigen+1/k;
            
        else
            
            if R_seed
                Pigen = evalR(['rgamma(' num2str(k) ',1)']);
                Pigen = Pigen';
            else
                Pigen = randg(1,k,1);
            end
            
            s=sum(Pigen);
            
            for j=1:k
                Pigen(j) = PiLow + Pigen(j) / s * (1 - k * PiLow);
                if (Pigen(j) < PiLow)
                    flag = 1;
                    break
                end
            end
            if (flag == 1)
                warning('FSDA:MixSim:WrongPiLow','PiLow is too high... generated equal mixing proportions...');
                Pigen=zeros(k,1)+1/k;
            end
        end
    end



    function Mugen=genMu(p, k, Lbound, Ubound)
        % genMu generates matrix of means of size k-by-p
        %
        %  Required input arguments:
        %               p - number of dimensions
        %               k - number of components
        %          Lbound - lower bound for the hypercube
        %          Ubound - upper bound for the hypercube
        %
        %  OUTPUT parameters
        % 		Mu - matrix of means of size k-by-p
        
        if R_seed
            % equivalent of 'rand(k,p)' in R is 'matrix(runif(k*p),k,p)'
            rn1s = ['matrix(runif(' num2str(k*p) '),' num2str(p) ',' num2str(k) ')'];
            rn1 = evalR(rn1s);
            rn1b = rn1';
            Mugen = Lbound + (Ubound-Lbound)*rn1b;
        else
            Mugen = Lbound + (Ubound-Lbound)*rand(k,p);
        end
    end


    function VC=genSigma(v)
        % genSigma generates covariance matrix based on (v + 1) observations
        %  Required input arguments:
        % 		v - number of dimensions
        %
        %  OUTPUT parameters
        %
        % 		VC - v-by-v covariance matrix based on v+1 observations
        % 		extracted from normal distribution
        n = v + 1;
        
        mu=zeros(v,1);
        x=zeros(n,v);
        
        
        for ii=1:n
            for jj=1:v
                if R_seed
                    % randn(1) in R is rnorm(1)
                    x(ii,jj)= evalR('rnorm(1)');
                    mu(jj) = mu(jj) + x(ii,jj);
                else
                    x(ii,jj)= randn(1);
                    mu(jj) = mu(jj) + x(ii,jj);
                end
            end
        end
        
        mu=mu/n;
        
        VC=zeros(v);
        
        for ii=1:n
            for jj=1:v
                for kk=1:v
                    VC(jj,kk) = VC(jj,kk) + (x(ii,jj) - mu(jj)) * (x(ii,kk) - mu(kk));
                end
            end
        end
        
        VC=VC/(n-1);
        % Note that VC is simply cov(X)
    end




    function S=genSigmaEcc(v, k, emax, hom)
        % genSigmaEcc generates covariance matrix with prespecified eccentricity
        %
        %  Required input arguments:
        %
        % 		v    - scalar, number of dimensions
        % 		k    - scalar, number of components
        % 		emax - scalar which defines maximum eccentricity
        %       hom  - scalar which specifies whether homogeneous (equal) or
        %              heterogeneous (different) covariance matrices must
        %              be generated. hom=1 implies homogeneous clusters
        %
        %  OUTPUT parameters
        %          S - array of size v-by-v-by-k containing the k
        %          variance-covariance matrices of the k groups
        %
        %    Remark: this procedure calls VC=genSigma and then checks
        %    whether the eigenvalues of VC satisfy a predetermined level of
        %    eccentricity
        
        % S = 3d array which contains the covariance matriced of the groups
        S=zeros(v,v,k);
        
        if hom == 0
            
            for kk=1:k
                
                VC=genSigma(v);
                S(:,:,kk)=VC;
                
                [V,L] = eig(VC);
                Eig=diag(L);
                minL=min(Eig);
                maxL=max(Eig);
                
                % e = eccentricity
                e = sqrt(1 - minL / maxL);
                
                if (e > emax)
                    L=zeros(v);
                    
                    for ii=1:v
                        Eig(ii) = maxL * (1 - emax * emax * (maxL - Eig(ii)) / (maxL - minL));
                        L(ii,ii)=Eig(ii);
                    end
                    R=V*L*V';
                    
                    % The two former instructions could be vectrorised as
                    % follows (but it does not seem necessary)
                    % Eig=maxL * (1 - emax^2* (maxL - Eig) / (maxL - minL));
                    % R=V*diag(Eig)*V';
                    S(:,:,kk)=R;
                end
            end
        else % homogeneous clusters
            
            VC=genSigma(v);
            for kk=1:k;
                S(:,:,kk)=VC;
            end
            
            [V,L] = eig(VC);
            Eig=diag(L);
            minL=min(Eig);
            maxL=max(Eig);
            
            
            % e = eccentricity
            e = sqrt(1 - minL / maxL);
            
            % If some simulated dispersion matrices have e>emax specified
            % by ecc all eigenvalues will be scaled in order to have
            % enew=emax;
            if (e > emax)
                L=zeros(v);
                
                for ii=1:v
                    Eig(ii) = maxL * (1 - (emax^2) * (maxL - Eig(ii)) / (maxL - minL));
                    L(ii,ii)=Eig(ii);
                end
                
                R=V*L*(V');
                for kk=1:k;
                    S(:,:,kk)=R;
                end
            end
        end
    end



    function S=genSphSigma(v,k,hom)
        %Generates spherical covariance matrices
        %
        %  Required input arguments:
        %
        % 		v    - scalar, number of dimensions
        % 		k    - scalar, number of components
        %       hom  - scalar which specifies whether homogeneous (equal) or
        %              heterogeneous (different) covariance matrices must
        %              be generated. hom=1 implies homogeneous clusters
        %
        %  OUTPUT parameters
        %
        % 		S :  3D array of size v-by-v-by-k containing the k
        % 		variance-covariance matrices
        %       if hom ==1
        %       S(:,:,j) = \sigma^2 I_v
        %       else
        %       S(:,:,j) = \sigma^2_j I_v
        %
        S=zeros(v,v,k);
        
        eyep=eye(v);
        if R_seed
            % rand(1) in R is runif(1)
            r = evalR('runif(1)');
        else
            r = rand(1);
        end
        
        for kkk=1:k
            if hom == 0
                if R_seed
                    % rand(1) in R is runif(1)
                    r = evalR('runif(1)');
                else
                    r = rand(1);
                end
            end
            LL=r*eyep;
            
            S(:,:,kkk)=LL;
        end
        
    end

    function  [c, OmegaMap2, BarOmega2, MaxOmega2, rcMax]=FindC(lower, upper, Omega, method, v, k, li, di, const1, fix, tol, lim)
        %find multiplier c to be applied to the cov matrices in the
        %interval [lower upper] in order to reach the required average or
        %maximum overlap
        %
        %  Required input arguments:
        %
        % lower : scalar - lower bound of the interval
        % upper : scalar - upper bound of the interval
        % Omega : scalar, associated with maximum or average overlapping requested
        % method : scalar which specifies whether average (method=0) or maximum
        %          overlap is requested
        %     v  : dimensionality
        %     k  : number of components
        % li, di, const1 : parameters needed for computing overlap,
        %          precalculated using routine ComputePars
        %    fix : vector of length k containing zeros or ones
        %          if fix(j) =1 cluster j does not participate to inflation
        %          or deflation. If fix=zeros(k,1) all clusters participate
        %          in inflation/deflation This parameter is used just if
        %          heterogeneous clusters are used
        %    tol : vector of length 2.
        %          tol(1) (which will be called tolmap) specifies
        %          the tolerance between the requested and empirical
        %          misclassification probabilities (default is 1e-06)
        %          tol(2) (which will be called tolnxc2) specifies the
        %          tolerance to use in routine ncx2mixtcdf (which computes cdf
        %          of linear combinations of non central chi2 distributions).
        %          The default value of tol(2) 1e-06
        %    lim : maximum number of integration terms default is 1e06.
        %          REMARK: Optional parameters tol and lim will be used by
        %          function ncx2mixtcdf.m which computes the cdf of a linear
        %          combination of non central chi2 r.v.. This is the
        %          probability of overlapping
        %
        %  OUTPUT parameters
        %
        %        c    : scalar inflation parameter
        %               c is a constant which is used to scale the covariance
        %               matrices of the groups in order to obtain a
        %               prespecified level of average or maximum overlap
        %   OmegaMap2 : k-by-k matrix containing map of misclassification
        %               probabilities
        %   BarOmega2 : scalar. Average overlap found using c
        %   MaxOmega2 : scalar. Maximum overlap found using c
        %      rcMax  : vector of length 2 containing the indexes associated
        %              to the pair of components producing the highest overlap
        %
        
        diff = Inf;
        stopIter = 200; % 500
        tolmap=tol(1);
        tolncx2=tol(2);
        
        sch = 0;
        asympt = 0;
        
        % Intervals which contain positive or negative powers of 2 are
        % considered. For example first interval is [0 1024] then if
        % MaxOmega2 (maximum overlap which has been found using c=512) is <
        % Omega (maximum required overlap), then the new interval becomes
        % [512 1024]  (c has to be increased and the new candidate c is
        % 0.5*(512+1024)) else the new interval becomes [0 512] (c has to
        % be decreased and the new candidate c is 0.5*(0+512)=256
        while abs(diff) > tolmap
            
            c = (lower + upper) / 2.0;
            
            [OmegaMap2, BarOmega2, MaxOmega2, rcMax]=GetOmegaMap(c, v, k, li, di, const1, fix, tolncx2, lim, asympt);
            
            if method == 0 % in this case average overlap is requested
                
                if BarOmega2 < Omega
                    % clusters are too far
                    lower = c;
                else
                    upper = c;
                end
                
                diff = BarOmega2 - Omega;
                
            else % in this case maximum overlap is requested
                
                if MaxOmega2 < Omega
                    % clusters are too far
                    lower = c;
                else
                    upper = c;
                end
                
                diff = MaxOmega2 - Omega;
                
            end
            
            sch = sch + 1;
            
            if sch == stopIter
                c = 0.0;
                disp(['Warning: required overlap was not reached in routine findC after ' num2str(stopIter) ' iterations...'])
                break
            end
        end
    end


    function  [c, OmegaMap2, BarOmega2, MaxOmega2, StdOmega2, rcMax]=FindCStd(lower, upper, Omega, v, k, li, di, const1, fix, tol, lim)
        %find multiplier c to be applied to the cov matrices in the
        %interval [lower upper] in order to reach the required std of
        %
        %  Required input arguments:
        %
        % lower : scalar - lower bound of the interval
        % upper : scalar - upper bound of the interval
        % Omega : scalar, associated with maximum or average overlapping requested
        %     v  : dimensionality
        %     k  : number of components
        % li, di, const1 : parameters needed for computing overlap,
        %          precalculated using routine ComputePars
        %    fix : vector of length k containing zeros or ones
        %          if fix(j) =1 cluster j does not participate to inflation
        %          or deflation. If fix=zeros(k,1) all clusters participate
        %          in inflation/deflation This parameter is used just if
        %          heterogeneous clusters are used
        %    tol : vector of length 2.
        %          tol(1) (which will be called tolmap) specifies
        %          the tolerance between the requested and empirical
        %          misclassification probabilities (default is 1e-06)
        %          tol(2) (which will be called tolnxc2) specifies the
        %          tolerance to use in routine ncx2mixtcdf (which computes cdf
        %          of linear combinations of non central chi2 distributions).
        %          The default value of tol(2) 1e-06
        %    lim : maximum number of integration terms default is 1e06.
        %          REMARK: Optional parameters tol and lim will be used by
        %          function ncx2mixtcdf.m which computes the cdf of a linear
        %          combination of non central chi2 r.v.. This is the
        %          probability of overlapping
        %
        %  OUTPUT parameters
        %
        %        c    : scalar inflation parameter
        %               c is a constant which is used to scale the covariance
        %               matrices of the groups in order to obtain a
        %               prespecified level of average or maximum overlap
        %   OmegaMap2 : k-by-k matrix containing map of misclassification
        %               probabilities
        %   BarOmega2 : scalar. Average overlap found using c
        %   MaxOmega2 : scalar. Maximum overlap found using c
        %   StdOmega2 : scalar. Std of overlap found using c
        %      rcMax  : vector of length 2 containing the indexes associated
        %              to the pair of components producing the highest overlap
        %
        
        diff = Inf;
        stopIter = 200; % 500
        tolmap=tol(1);
        tolncx2=tol(2);
        
        sch = 0;
        asympt = 0;
        
        % Intervals which contain positive or negative powers of 2 are
        % considered. For example first interval is [0 1024] then if
        % MaxOmega2 (maximum overlap which has been found using c=512) is <
        % Omega (maximum required overlap), then the new interval becomes
        % [512 1024]  (c has to be increased and the new candidate c is
        % 0.5*(512+1024)) else the new interval becomes [0 512] (c has to
        % be decreased and the new candidate c is 0.5*(0+512)=256
        while abs(diff) > tolmap
            
            c = (lower + upper) / 2.0;
            
            [OmegaMap2, BarOmega2, MaxOmega2, rcMax]=GetOmegaMap(c, v, k, li, di, const1, fix, tolncx2, lim, asympt);
            
            % Compute StdOverlap
            % Extract elements above diagonal and compute Stdoverlap
            cand=triu(OmegaMap2,1)+(tril(OmegaMap2,-1))';
            overlapv=cand(:);
            %             overlapv=overlapv(overlapv>0);
            %             if length(overlapv)<combk2;
            %                 overlapc=[overlapv; zeros(combk2-length(overlapv),1)];
            %             else
            %                 overlapc=overlapv;
            %             end
            %             % Compute standard deviation of overlap for current
            %             % solution
            %             StdOmega2=std(overlapc);
            StdOmega2=std(overlapv(indabovediag));
            
            
            if StdOmega2 < Omega
                % clusters are too far
                lower = c;
            else
                upper = c;
            end
            
            diff = StdOmega2 - Omega;
            
            sch = sch + 1;
            
            if sch == stopIter
                c = 0.0;
                disp(['Warning: required overlap was not reached in routine findC after ' num2str(stopIter) ' iterations...'])
                break
            end
        end
    end

end