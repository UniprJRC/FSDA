function [RAW,REW,varargout] = mcd(Y,varargin)
%mcd computes Minimum Covariance Determinant
%
%<a href="matlab: docsearchFS('mcd')">Link to the help function</a>
%
%  Required input arguments:
%
%    Y: Data matrix containing n observations on v variables.
%       Rows of Y represent observations, and columns represent variables.
%       Missing values (NaN's) and infinite values (Inf's) are allowed,
%       since observations (rows) with missing or infinite values will
%       be excluded from the computations.
%
%  Optional input arguments:
%
%      bdp    : Breakdown point. Scalar. (Number between 0
%               and 0.5). The default value is 0.5.
%               Example - 'bdp',1/4
%               Data Types - double
%
%
%      bestr  : Number of best solutions to store. Scalar. Number of "best locations"
%               to remember from the subsamples. These will be later iterated until
%               convergence (default=5)
%               Example - 'bestr',10
%               Data Types - double
%
%  betathresh : Distribution to use. Boolean. If betathresh = true the distribution
%               which is used to declare units as outliers is a mixture of Rocke
%               scaled F distribution and Beta else (default) traditional chi^2
%               distribution is used.
%               Example - 'betathresh',false
%               Data Types - logical
%
%
%     conflev : Confidence level. Scalar. Number between 0 and 1 containing
%               confidence level which is used to declare units as outliers
%               after reweighting.
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%               Example - 'conflev',0.99
%               Data Types - double
%
%  conflevrew : Confidence level for use for reweighting. Scalar. Number
%               between 0 and 1 containing confidence level which is used to do
%               the reweighting step. Default value is the one specified in
%               previous option conflev.
%               Example - 'conflevrew',0.99
%               Data Types - double
%
%        msg  : Display or not messages on the screen.
%               Scalar. If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the final
%               estimator else no message is displayed on the screen.
%               Example - 'msg',1
%               Data Types - double
%
%
%      nocheck: No check on input data. Scalar. If nocheck is equal to 1 no check
%               is performed on matrix Y. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
%      nsamp  : Number of subsamples. Scalar or matrix. 
%               If nsamp is a scalar, it contains the number of subsamples of size
%               v+1 which have to be extracted (if not given, default is
%               nsamp=1000). If nsamp=0 all subsets will be extracted. If
%               nsamp is a matrix it contains in the rows the indexes of
%               the subsets which have to be extracted. nsamp in this case
%               can be conveniently generated  by function subsets.
%               Example - 'nsamp',10000
%               Data Types - double
%
%    refsteps : Number of refining iterations. Scalar. Number of refining
%               iterations in each subsample (default = 3).
%               refsteps = 0 means "raw-subsampling" without iterations.
%               Example - 'refsteps',10
%               Data Types - double
%
%refstepsbestr: Number of refining iterations. Scalar. Number of refining iterations
%               for each best subset (default = 50).
%               Example - 'refstepsbestr',10
%               Data Types - double
%
%     reftol  : Refining steps tolerance. Scalar. Tolerance for the refining steps.
%               The default value is 1e-6;
%               Example - 'reftol',1e-8
%               Data Types - double
%
%
% reftolbestr : Tolerance for refining steps. Scalar. Value of tolerance for the
%               refining steps for each of the best subsets.
%               The default value is 1e-8;
%               Example - 'reftolbestr',1e-8
%               Data Types - double
%
%
%  restrfactor: Restriction factor. Scalar. 
%               Positive scalar greater or equal than 1 which constrains
%               the allowed differences of the eigenvalues of the scatter
%               matrix in each concentration step. The default value of
%               restrfactor is Inf, that is no constraint is imposed on the
%               eigenvalues of the covariace matrix in each subset.
%               Example - 'restrfactor',100
%               Data Types - double
%
%smallsamplecor: small sample correction factor. Boolean. Boolean which
%               defines whether to use or not small sample correction
%               factor to inflate the scale estimate.  If it is equal to
%               true the small sample correction factor is used. The
%               default value of smallsamplecor is true, that is the
%               correction is used. See
%               http://users.ugent.be/~svaelst/publications/corrections.pdf
%               for further details about the correction factor. 
%               Example - 'smallsamplecor',true 
%               Data Types - logical
%
%     tolMCD  : Tolerance to declare a subset as singular. Scalar. The
%               default value of tolMCD is exp(-50*v).
%               Example - 'tolMCD',1e-20
%               Data Types - double
%
%    ysaveRAW : save Y. boolean. Boolean that is set to true to request that the data
%               matrix Y is saved into the output structure RAW. This feature is
%               meant at simplifying the use of function malindexplot.
%               Default is false, i.e. no saving is done.
%               Example - 'ysaveRAW',true
%               Data Types - logical
%
%    ysaveREW : save Y. Boolean. Boolean that is set to true to request that the data
%               matrix Y is saved into the output structure REW. This feature is
%               meant at simplifying the use of function malindexplot.
%               Default is false, i.e. no saving is done.
%               Example - 'ysaveREW',true
%               Data Types - logical
%
%       plots : Plot on the screen. Scalar or structure.
%               If plots is a structure or scalar equal to 1, generates:
%               (1) two plots of Mahalanobis distances (raw and reweighted)
%               against index number. The confidence level used to draw the
%               confidence bands for the MD is given by the input option
%               conflev. If conflev is not specified a nominal 0.975
%               confidence interval will be used.
%               (2) two scatter plot matrices with the outliers (from raw
%               and reweighted mcd estimators) highlighted.
%               If plots is a structure it may contain the following fields
%                   plots.labeladd = if this option is '1', the outliers in the
%                       spm are labelled with their unit row index. The
%                       default value is labeladd='', i.e. no label is
%                       added.
%                   plots.nameY = cell array of strings containing the labels of
%                       the variables. As default value, the labels which
%                       are added are Y1, ...Yv.
%               Example - 'plots',1
%               Data Types - double or structure
%
%  Output:
%
%  The output consists of two structures RAW and REW. RAW refers to raw
%  MCD, on the other hand, REW refers to reweighted MCD
%
%         RAW:   structure which contains the following fields
%
%         RAW.h    = scalar. The number of observations that have
%                    determined the MCD estimator
%         RAW.loc  = 1 x v  vector containing raw MCD location of the data
%         RAW.cov  = robust MCD estimate of covariance matrix. 
%                    It is the raw MCD covariance matrix (multiplied by a
%                    finite sample correction factor and an asymptotic
%                    consistency factor).
%           RAW.cor= The raw MCD correlation matrix
%           RAW.obj= The determinant of the raw MCD covariance matrix.
%           RAW.bs = (v+1) x 1 vector containing the units forming best
%                    subset associated with MCD estimate of location.
%           RAW.md = n x 1 vector containing the estimates of the robust
%                    Mahalanobis distances (in squared units). This vector
%                    contains the distances of each observation from the
%                    raw MCD location of the data, relative to the raw MCD
%                    scatter matrix RAW.cov
%     RAW.outliers = A vector containing the list of the units declared as
%                    outliers using confidence level specified in input
%                    scalar conflev
%      RAW.conflev = Confidence level that was used to declare outliers
%      RAW.singsub = Number of subsets without full rank. Notice that
%                    out.singsub > 0.1*(number of subsamples) produces a
%                    warning
%      RAW.weights = n x 1 vector containing the estimates of the weights.
%                    Weights assume values 0 or 1. Weight is 1 if the
%                    associated observation has been used to compute
%                    centroid and covariance matrix. These weights
%                    determine which observations are used to compute the
%                    final MCD estimates. Unless there is a perfect fit
%                    sum(RAW.weights)=h
%        RAW.plane = In case of an exact fit, RAW.plane contains the
%                    coefficients of a (hyper)plane
%                    a_1(x_i1-m_1)+...+a_p(x_ip-m_p)=0
%                    containing at least h observations, where (m_1,...,m_p)
%                    is the MCD location of these observations.
%                    This field is present only if there is exact fit.
%            RAW.Y = Data matrix Y. This field is present only if option
%                    ysaveRAW was set to 1.
%        RAW.class = 'mcd'
%
%         REW : structure which contains the following fields:
%
%       REW.loc    = The robust location of the data, obtained after
%                    reweighting, if the raw MCD is not singular.
%                    Otherwise the raw MCD center is given here.
%       REW.cov    = The robust covariance matrix, obtained after
%                    reweighting and multiplying with a finite sample
%                    correction factor and an asymptotic consistency
%                    factor, if the raw MCD is not singular.  Otherwise the
%                    raw MCD covariance matrix is given here.
%       REW.cor    = The robust correlation matrix, obtained after reweighting
%       REW.md     = n x 1 vector containing the estimates of the robust
%                    Mahalanobis distances (in squared units). This vector
%                    contains the distances of each observation from the
%                    reweighted MCD location of the data, relative to the
%                    reweighted MCD scatter of the data These distances
%                    allow us to easily identify the outliers. If the
%                    reweighted MCD is singular, RAW.md is given here.
%     REW.outliers = A vector containing the list of the units declared as
%                    outliers after reweighting.
%      REW.weights = n x 1 vector containing the estimates of the weights.
%                    Weights assume values 0 or 1. Weight is 0 if the
%                    associated observation has been declared outlier.
%                    These weights determine which observations are used to
%                    compute the final MCD estimates.
%                    Remark: if the reweighted MCD is singular, RAW.weights
%                    is given here.
%       REW.method = In case of an exact fit, REW.method contains a
%                    character string containing information about the
%                    method and about singular subsamples (if any). This
%                    field is present only if there is exact fit.
%       REW.plane  = In case of an exact fit, REW.plane contains the
%                    coefficients of a (hyper)plane
%                    a_1(x_i1-m_1)+...+a_p(x_ip-m_p)=0
%                    containing at least h observations, where
%                    (m_1,...,m_p). This field is present only if there is
%                    exact fit.
%            REW.Y = Data matrix Y. The field is present only if option
%                    ysaveREW was set to 1.
%        REW.class = 'mcdr'.
%
%  Optional Output:
%
%            C     : matrix of size nsamp-by-v which contains the indices
%                    of the subsamples extracted for
%                    computing the estimate.
%
% More About:
%
% MCD computes the MCD estimator of a multivariate data set.  This
% estimator is given by the subset of h observations with smallest
% covariance determinant.  The MCD location estimate is then the mean of
% those h points, and the MCD scatter estimate is their covariance matrix.
% The default value of h is roughly 0.5n (where n is the total number of
% observations), but the user may choose each value between n/2 and n.
%
% The MCD method is intended for continuous Gaussian variables, and assumes
% that the number of observations n is at least 5 times the number of
% variables p. If p is too large relative to n, it is better to use options
% betathresh=1 (that is to use a threshold based on beta distribution (for
% the units which determine the centroid and the covariance matrix) and F
% distribution (for the units which are excluded from the computation of
% centroid and covariance matrix)
%
% The MCD method was introduced in:
%
%   Rousseeuw, P.J. (1984), Least Median of Squares Regression,
%   "Journal of the American Statistical Association", Vol. 79, pp. 871-881.
%
%   The program below uses the technique of concentration steps described in
%
%   Rousseeuw, P.J. and Van Driessen, K. (1999), "A Fast Algorithm for the
%   Minimum Covariance Determinant Estimator," Technometrics, 41, pp. 212-223.
%
% The MCD is a robust method in the sense that the estimates are not unduly
% influenced by outliers in the data, even if there are many outliers.
% Due to the MCD's robustness, we can detect outliers by their large
% robust distances. The latter are defined like the usual Mahalanobis
% distance, but based on the MCD location estimate and scatter matrix
% (instead of the nonrobust sample mean and covariance matrix).
%
%
% Remark: when more than h observations lie on a (hyper)plane, (perfect fit
% case) the program still yields the MCD location and scatter matrix, the
% latter being singular (as it should be), as well as the equation of the
% hyperplane.
%
% See also: mve.m
%
% References:
%
% Rousseeuw, P.J. and Van Driessen, K. 1999. A fast algorithm for the
% minimum covariance determinant estimator. Technometrics, 41:212?223.
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
% Acknowledgements:
%
% This function follows the lines of MATLAB/R code developed during the
% years by many authors. In particular, parts of the code rely on the LIBRA
% mcd implementation of Hubert and Verboven. For more details, see:
% http://wis.kuleuven.be/stat/robust/LIBRA.html,
% http://www.econ.kuleuven.be/public/NDBAE06/programs/
% and the R library Robustbase http://robustbase.r-forge.r-project.org/
% The core of our routines, e.g. the resampling approach, however, has been
% completely redesigned, with considerable increase of the computational
% performance. Note that, for the moment, FSDA does not adopt the 'divide
% and conquer' partitioning method proposed by Rousseeuw and Van Driessen
% to speed up computations for large datasets. This partitioning method is
% applied in the R and LIBRA implementations of the mcd.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mcd')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % mcd with all default options.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    RAW=mcd(Ycont);
%}

%{
    %% mcd with optional arguments.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    RAW=mcd(Ycont,'plots',1);
%}

%{
    % mcd monitoring the reweighted estimates.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [RAW,REW]=mcd(Ycont);
%}

%{
    % mcd monitoring the exctracted subsamples.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [RAW,REW,C]=mcd(Ycont);
%}

%{
    %% mcd applied to the aircraft data (default plots).
    % See Pison et al. 2002, Metrika.
    X = load('aircraft.txt');
    Y = X(:,1:end-1);
    [RAW,RES] = mcd(Y,'bdp',0.25,'plots',1);
%}

%{
    %% mcd applied to the aircraft data (plots using the scale of Pison et al).
    % See Pison et al. 2002, Metrika.
    X = load('aircraft.txt');
    Y = X(:,1:end-1);
    [RAW,REW] = mcd(Y,'bdp',0.25,'ysaveRAW',1);

    v=size(Y,2);
    % Compare the following figure with panel (b) of Fig. 8 of Pison et al.
    ylimy=[0 36];
    malindexplot(RAW.md,v,'conflev',0.975,'laby','robust distances','numlab',RAW.outliers,'ylimy',ylimy);
    title('Corrected MCD')

    % Compare the following figure with panel (4) of Fig. 8 of Pison et al.
    ylimy=[0 36];
    malindexplot(REW.md,v,'conflev',0.975,'laby','robust distances','numlab',REW.outliers,'ylimy',ylimy);
    title('Corrected reweighted MCD')

%}

%% Beginning of code

[n, v]=size(Y);

% default value of break down point
bdpdef=0.5;

% If the number of all possible subsets is <10000 the default is to extract
% all subsets otherwise just 1000.
% Notice that we use bc, a fast version of nchoosek. One may also use the
% approximation floor(exp(gammaln(n+1)-gammaln(n-p+1)-gammaln(p+1))+0.5)
ncomb=bc(n,v+1);
nsampdef=min(1000,ncomb);

% default value of number of refining iterations (C steps) for each
% extracted subset
refstepsdef=3;
% default value of tolerance for the refining steps convergence for  each
% extracted subset
reftoldef=1e-6;
% default value of number of best locs to remember
bestrdef=5;
% default value of number of refining iterations (C steps) for best subsets
refstepsbestrdef=50;
% default value of tolerance for the refining steps convergence for best
% subsets
reftolbestrdef=1e-8;
% other generic tolerance default
generictol = 1e-8;

% Tolerance to declare a subset as singular
% It was set to
% tolMCDdef=exp(-50*v);
% but this is useless, as the roundoff level is eps = 2^(-52)
tolMCDdef=eps('double');

% if smallsamplecor ==1 (then small sample correction factor is applied to the
% estimate of the scale)
smallsamplecor=true;

restrfactordef=Inf;

% store default values in the structure options
options=struct('nsamp',nsampdef,'refsteps',refstepsdef,'bestr',bestrdef,...
    'reftol',reftoldef,...
    'refstepsbestr',refstepsbestrdef,'reftolbestr',reftolbestrdef,...
    'bdp',bdpdef,'plots',0,'conflev',0.975,'conflevrew','',...
    'betathresh',0,'nocheck',0,'msg',1,'tolMCD',tolMCDdef,...
    'ysaveRAW',false,'ysaveREW',false,'smallsamplecor',smallsamplecor,...
    'restrfactor',restrfactordef);

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mcd:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

bdp = options.bdp;              % break down point
refsteps = options.refsteps;    % refining steps
bestr = options.bestr;          % best locs for refining steps till convergence
nsamp = options.nsamp;          % subsamples to extract
reftol = options.reftol;        % tolerance for refining steps

refstepsbestr=options.refstepsbestr;  % refining steps for the best subsets
reftolbestr=options.reftolbestr;      % tolerance for refining steps for the best subsets

% tolMCD threshold under which the determinant of the covariance matrix is
% thought to be singular
tolMCD=options.tolMCD;

% If msg =1 the total estimated time to compute MCD is printed on the screen
msg=options.msg;

restrfactor=options.restrfactor;

% Initialize the matrices which contain the best "bestr" estimates of
% location, indexes of subsets, cov matrices and objective function
bestlocs = zeros(bestr, v);

bestcovs=zeros(v,v,bestr);
bestobjs = Inf * ones(bestr,1);

% singsub = scalar which will contain the number of singular subsets which
% are extracted (that is the subsets of size p which are not full rank)
singsub=0;

% hmin is the minimum number of observations whose covariance determinant
% will be minimized.
hmin=floor(2*floor((n+v+1)/2)-n+2*(n-floor((n+v+1)/2))*(0.5));

h=floor(2*floor((n+v+1)/2)-n+2*(n-floor((n+v+1)/2))*(1-bdp));

if h < hmin
    error('FSDA:mcd:Wrongh',['The MCD must cover at least ' int2str(hmin) ' observations.'])
elseif h > n
    error('FSDA:mcd:Wrongh','h is greater than the number of non-missings and non-infinites.')
end

% bestsubset is the matrix which will contain the indexes of the bestr
% subsets. Each row refers to a subset.
% Remark: note that each subset should have v+1 elements. However, due to
% the fact that if the subset is singulat we continue adding randomly
% elements to it up to when it becomes non singular, it is possible that
% certain subsets have more than v+1 elements
bestsubset = zeros(bestr, h,'int8');

% write in structure RAW the value of h
RAW=struct;
RAW.h=h;
REW=struct;

z=zeros(1,v);
weights=zeros(1,n);
seq=1:n;
conflev = options.conflev;
thresh=chi2inv(conflev,v);

%% Standardization of the data with medians and mads

% The standardization of the data will now be performed.
med=median(Y);
madY = median(abs(bsxfun(@minus,Y, med)));

ii=find(madY < eps, 1 );
if ~isempty(ii)
    % The h-th order statistic is zero for the ii-th variable. The array plane contains
    % all the observations which have the same value for the ii-th variable.
    plane=find(abs(Y(:,ii)-med(ii)) < eps)';
    meanplane=mean(Y(plane,:));
    weights(plane)=1;
    if v==1
        out.weights = weights;
        out.loc     = meanplane;
        [out.cov,out.objective]=deal(0); %#ok<STRNU>
        REW.method=sprintf('\nUnivariate location and scale estimation.');
        REW.method=char(REW.method,sprintf('%g of the %g observations are identical.',length(plane),n));
        disp(REW.method);
    else
        z(ii)=1;
        REW.plane=z;
        covplane=cov(Y(plane,:));
        REW.method=sprintf('\nMinimum Covariance Determinant Estimator (fastMCD).');
        [RAW.center,RAW.cov,REW.center,REW.cov,RAW.objective,RAW.weights,REW.weights]=displRAW(3,length(plane),weights,n,v,meanplane,covplane,REW.method,z,h,ii);
    end
    return
end
% Standardization of the data with location estimate (median) and scale
% estimate (mad)
Y = bsxfun(@minus,Y, med);
Y = bsxfun(@rdivide, Y, madY);

% The standardized classical estimates are now computed.
clmean=mean(Y);
clcov=cov(Y);

%% If there is just one variable and h<n

% The univariate non-classical case is now handled.
if v==1 && h~=n
    
    REW.method=sprintf('\nUnivariate location and scale estimation.');
    
    [RAW.loc,RAW.cov]=mcduni(Y,h,1-bdp,smallsamplecor);
    
    % Assign weight=1 to the h units which show the smallest h squared
    % residuals
    residuals2=(Y-RAW.loc).^2;

    %{
    [sor,soridx]=sort(residuals2);    % deleted
    weights=zeros(n,1);               % deleted
    weights(soridx(1:h))=1;           % deleted
    RAW.weights=weights;              % deleted
    %}
    [hos , sor]= quickselectFS(residuals2,h);
    weights = residuals2<=hos;
    RAW.weights=weights;
    
    RAW.obj=1/(h-1)*sum(sor(1:h))*prod(madY)^2;
    md=residuals2/RAW.cov;
    
    weights=md<thresh;
    RAW.md=md;
    RAW.outliers=seq(md > thresh);
    
    % Reweighted part
    REW.loc=mean(Y(weights==1,:));
    REW.cov=cov(Y(weights==1,:));
    
    cfactor = consistencyfactor(sum(weights),n,v);
    if smallsamplecor == true
        cfactor = cfactor*corfactorREW(v,n,1-bdp);
    end
    
    REW.cov=cfactor*REW.cov;
    
    md=(Y-REW.loc).^2/REW.cov;
    
    REW.weights= md <= thresh;
    REW.outliers=seq(md > thresh);
    REW.md=md';
    
    % Transform back to the original scale
    [RAW.cov,RAW.loc]=trafo(RAW.cov,RAW.loc,med,madY);
    [REW.cov,REW.loc]=trafo(REW.cov,REW.loc,med,madY);
    
    if msg
        disp(REW.method);
    end
    
    return
end

%% If the determinant of cov(Y) is not full rank

if det(clcov) < tolMCD
    % all observations lie on a hyperplane.
    sigma = eigs_sigma(clcov);
    [z, ~]=eigs(clcov,1,sigma,struct('disp',0));
    weights=ones(n,1);
    
    correl=clcov./(sqrt(diag(clcov))*sqrt(diag(clcov))');
    
    [clcov,clmean]=trafo(clcov,clmean,med,madY);
    [RAW.center,RAW.cov,REW.center,REW.cov,RAW.objective,RAW.weights,...
        REW.weights] = displRAW(1,n,weights,n,v,clmean,clcov,...
        'Minimum Covariance Determinant Estimator',z./madY');
    
    [REW.cor,RAW.cor]=deal(correl);
    return
end

%% h==n is the case in which there is no trimming (classical case)

if h==n
    if msg
        disp(['The MCD estimates are equal to the classical estimates h=n=',num2str(h)]);
    end
    %  REW.method=char(REW.method,msgn);
    
    RAW.loc=clmean;
    RAW.cov=clcov;
    RAW.obj=det(clcov);
    md=mahalFS(Y,clmean,clcov);
    
    RAW.md=md;
    RAW.cor=RAW.cov./(sqrt(diag(RAW.cov))*sqrt(diag(RAW.cov))');
    
    weights=md<thresh;
    RAW.weights=weights;
    
    REW.loc=mean(Y(weights==1,:));
    REW.cov=cov(Y(weights==1,:));
    
    
    if det(REW.cov) < tolMCD
        [center,covar,z,correl,~,count]=exactfit(Y,NaN,med,madY,RAW.cor,REW.loc,REW.cov,n);
        
        REW.plane=z;
        
        REW.method=displREW(count,n,v,center,covar,REW.method,z,RAW.cor,correl);
        [RAW.cov,RAW.center]=trafo(RAW.cov,RAW.loc,med,madY);
        [REW.cov,REW.center]=trafo(REW.cov,REW.loc,med,madY);
        REW.md=RAW.md;
    else
        md=mahalFS(Y,REW.loc,REW.cov);
        weights=md<thresh;
        
        [RAW.cov,RAW.center]=trafo(RAW.cov,RAW.loc,med,madY);
        [REW.cov,REW.center]=trafo(REW.cov,REW.loc,med,madY);
        REW.md=md;
    end
    RAW.obj=RAW.obj*prod(madY)^2;
    REW.weights=weights;
    % disp(REW.method)
    
    return
end

%% Extract in the rows of matrix C the indexes of all required subsets

if isscalar(nsamp) % nsamp
    [C,nselected] = subsets(nsamp,n,v+1,ncomb,1);
else
    C=nsamp;
    if size(C,2) ~= v+1
        % check that the number of columns of C is v+1
        error('FSDA:mcd:WrongC','Matrix nsamp must have %.0f columns', v+1)
    end
    
    nselected=size(C,1);
end

% Store the indices in varargout
if nargout==3
    varargout={C};
end
% initialise and start timer.
tsampling = ceil(min(nselected/100 , 1000));
time=zeros(tsampling,1);

% ij is a scalar used to ensure that the best first bestr non singular
% subsets are stored
ij=1;
for i = 1:nselected
    if i <= tsampling, tic; end
    
    % extract a subset of size v+1 .
    index = C(i,:);
    
    Yj = Y(index,:);
    
    locj = mean(Yj);        % centroid of subset
    Sj = cov(Yj);           % covariance of subset
    
    % Check if the subset is in general position (rank<v)
    if det(Sj)< tolMCD
        singsub = singsub + 1;
        
        % The trial subsample is singular.
        % We distinguish two cases :
        %
        % 1. There are h or more observations in the subdataset that lie
        %    on the hyperplane and thus an exact fit.
        %
        % 2. There aren't h observations in the subdataset that lie on the
        %    hyperplane. We then extend the trial subsample until it isn't
        %    singular anymore.
        %
        % eigvct : contains the coefficients of the hyperplane.
        % If sigma is omitted, the eigenvalues largest in magnitude are
        % found. If sigma is zero, the eigenvalues smallest in magnitude
        % are found.
        
        sigma = eigs_sigma(Sj);
        [eigvct, ~]=eigs(Sj,1,sigma,struct('disp',0));
        
        % REMARK: Eigenvalues may be sensitive to perturbations. Normal
        % roundoff errors in the floating-point arithmetic computation have
        % the same effect of perturbing the original matrix. Thus, it would
        % be worth estimating the sensitivity of the eigenvalues using the
        % condition number of the matrix of eigenvectors. We are stydying
        % algorithms for the condition number more efficient than the
        % MATLAB function condest. Thus, lines below are for the moment
        % commented.
        %         perturb = condest(eigvct);
        %         if perturb > 10^3
        %             message = sprintf(' A perturbation in the covariance of subset could result \n in perturbations in the hyperplane coefficients \n that are %8.2f times as large.' , perturb);
        %             disp(message);
        %         end
        
        %dist=abs(sum((data-repmat(meanvct,n,1))'.*repmat(eigvct,1,n)));
        Ytilde  = bsxfun(@minus,Y, locj);
        dist    = abs(sum(bsxfun(@times,Ytilde, eigvct'),2));
        
        obsinplane = find(dist < generictol);
        % count : number of observations that lie on the hyperplane.
        count=length(obsinplane);
        
        if count >= h
            
            [center,covar,eigvct,correl]=exactfit(Y,obsinplane,med,madY,eigvct);
            REW.plane=eigvct;
            weights(obsinplane)=1;
            [RAW.center,RAW.cov,REW.center,REW.cov,RAW.objective,...
                RAW.weights,REW.weights]=displRAW(2,count,weights,n,v,center,covar,...
                'fastMCD',eigvct,correl);
            
            
            [REW.cor,RAW.cor]=deal(correl);
            
            return
            
        else
            % Find the indexes which do not belong to index
            setd=setdiff(seq,index);
            % Do a random permutation of these indexed
            % using utility of toolbox FSDA called shuffling
            setd=shuffling(int16(setd));
            % Remark: the instruction above is much more efficient than
            % setd=setd(randperm(length(setd)));
            
            jj=1;
            % Define a vector indexext which contains the original
            % indexes and then the random permutation of the reminaning
            % indexes
            indexext=[index setd];
            lindex=length(index);
            while det(Sj) < tolMCD
                % At each iteration extend the subset until the subset
                % becomes full rank
                index=indexext(1:lindex+jj);
                jj=jj+1;
                Sj=cov(Y(index,:));
            end
            locj=mean(Y(index,:));
            
        end
    end
    
    % Function IRWLSmult performs refsteps concentration steps of IRLS on elemental
    % start. Input:
    % - Y = datamatrix of dimension n x v
    % - locj = row vector containing (robust) centroid
    % - Sj = v x v covariance matrix
    % - refsteps = number of refining iterations
    % - reftol = tolerance for convergence of refining iterations
    outIRWLS = IRWLSmcd(Y, locj, Sj, h, refsteps, reftol, restrfactor);
    
    % If the value of the objective function is smaller than tolMCD
    % we have a perfect fit situation, that is there are h observations
    % that lie on the hyperplane.
    if outIRWLS.obj < tolMCD
        
        [center,covar,z,correl,obsinplane,count]=exactfit(Y,NaN,med,madY,NaN,...
            med,Sj,n);
        REW.plane=z;
        weights(obsinplane)=1;
        [RAW.center,RAW.cov,REW.center,REW.cov,RAW.objective,...
            RAW.weights,REW.weights]=displRAW(2,count,weights,n,v,center,covar,...
            'fastmcd',z);
        
        [REW.cor,RAW.cor]=deal(correl);
        
        return
        
    end
    
    % The output of IRWLSmult is a structure containing centroid, cov
    % matrix and estimate value of the objective function (which has
    % been minimized, that is |cov|
    locrw = outIRWLS.loc;
    covrw = outIRWLS.cov;
    objrw = outIRWLS.obj;
    
    % Compute Mahalanobis distances using locrw and covrw
    % mdrw = sqrt(mahalFS(Y,locrw,covrw));
    
    % to find s, save first the best bestr scales and shape matrices
    % (deriving from non singular subsets) and, from iteration bestr+1
    % (associated to another non singular subset), replace the worst scale
    % with a better one as follows
    if ij > bestr
        % from the second step check whether new loc and new shape belong
        % to the top best loc; if so keep loc and shape with
        % corresponding scale.
        
        
        if  objrw < max(bestobjs)
            % Find position of the maximum value of bestscale
            [~,ind] = max(bestobjs);
            bestobjs(ind) = objrw;
            bestlocs(ind,:) = locrw;
            bestcovs(:,:,ind) = covrw;
            % best subset associated with minimum value
            % of the objective function
            bestsubset(ind,1:length(index))=index;
        else
            
        end
    else
        bestobjs(ij) = objrw;
        bestlocs(ij,:) = locrw;
        bestcovs(:,:,ij) = covrw;
        bestsubset(ij,1:length(index)) = index;
        ij=ij+1;
    end
    
    
    % Write total estimation time to compute final estimate
    if i <= tsampling
        
        % sampling time until step tsampling
        time(i)=toc;
    elseif i==tsampling+1
        % stop sampling and print the estimated time
        if msg==1
            fprintf('Total estimated time to complete MCD: %5.2f seconds \n', nselected*median(time));
        end
    end
    
    
end
if singsub==nselected
    error('FSDA:mcd:NoFullRank','No subset had full rank. Please increase the number of subsets or check your design matrix X')
end

if singsub/nselected>0.1
    disp('------------------------------')
    disp(['Warning: Number of subsets without full rank equal to ' num2str(100*singsub/nsamp) '%'])
end

% perform C-steps on best 'bestr' solutions, till convergence or for a
% maximum of refstepsbestr steps using a convergence tolerance as specified
% by scalar reftolbestr

% this is to ensure that the condition tmp.scale < superbestscale in the
% next if statement is satisfied at least once
superbestobj = Inf;
for i=1:bestr
    tmp = IRWLSmcd(Y,bestlocs(i,:), bestcovs(:,:,i),h,refstepsbestr,reftolbestr,restrfactor);
    
    if tmp.obj < superbestobj
        superbestobj    = tmp.obj;
        superbestloc    = tmp.loc;
        superbestcov    = tmp.cov;
        superbestsubset = bestsubset(i,:);
        % weights = tmp.weights;
    end
end

% Remove the extra zeros in superbestsubset;
superbestsubset(superbestsubset==0)=[];

RAW.class   = 'mcd';
RAW.obj   = superbestobj;       % value of the objective function

RAW.bs=superbestsubset;

% cfactor: if we multiply the raw MCD covariance matrix by factor, we obtain
% consistency when the data come from a multivariate normal distribution.
cfactor = consistencyfactor(h,n,v);

% Apply small sample correction factor
if smallsamplecor == true
    cfactor = cfactor*corfactorRAW(v,n,1-bdp);
end

RAW.cov = cfactor*superbestcov;
RAW.obj = superbestobj*prod(madY)^2;

% Given that the data had been previously standardized
% it is necessary to find covariance and location in the original scale
[RAW.cov,RAW.loc]=trafo(RAW.cov,superbestloc,med,madY);

% Compute correlation matrix
RAW.cor=superbestcov./(sqrt(diag(superbestcov))*sqrt(diag(superbestcov))');

%Mahalanobis distances on standardized data
% Remember that MD are invariant under linear transformations of the data
md=mahalFS(Y,superbestloc,cfactor*superbestcov);

% Store vector of Mahalanobis distances (in squared units)
RAW.md = md;

% The first h smallest ordered Mahalanobis distances have weight equal to 1
%{
[~,soridx]=sort(md);       %deleted
weights=false(n,1);        %deleted
weights(soridx(1:h))=true; %deleted
RAW.weights=weights;       %deleted
%}
hos = quickselectFS(md,h);
weights = md<=hos;
RAW.weights=weights;

% Specify the distribution to use to compare Mahalanobis distances
% if betathresh==1 a mixture of scaled beta and F is used
% else traditional chi2 distribution is used
betathresh=options.betathresh;

conflevrew=options.conflevrew;
if isempty(conflevrew)
    conflevrew=conflev;
end

if betathresh==1
    % Compare md with the scaled F distribution of Hardin and Rocke (2005)
    [~,dfhr,~]=rockecs(n,v,h);
    xcor=(dfhr-v+1)/(v*dfhr);
    thresh=finv(conflevrew,v,dfhr-v+1)/xcor;
    weights=md<thresh;
else
    % Compare md with the chi2 distribution
    weights=md<chi2inv(conflevrew,v);
end

% RAW.weights=weights;
RAW.outliers=seq(md > thresh);

%  Store confidence level
RAW.conflev=conflev;

% Store total number of singular subsets
RAW.singsub=singsub;

% Compute reweighted estimates of location and covariance
% using standardized data
REW.loc=mean(Y(weights==1,:));
REW.cov=cov(Y(weights==1,:));

if betathresh==1
    % Consistency factor based on nominal trimming of conflevrew
    hh=floor(n*conflevrew);
    cfactor = consistencyfactor(hh,n,v);
else
    % Apply consistency factor to reweighted estimate of covariance
    hrew=sum(weights);
    if hrew<n
        cfactor=consistencyfactor(hrew,n,v);
    else
        cfactor=1;
    end
    
    if smallsamplecor==true
        % Apply small sample correction factor
        cfactor = cfactor*corfactorREW(v,n,1-bdp);
    end
    
end

REW.cov=cfactor*REW.cov;

% Find cov and location reweighted estimates on the original scale
[trcov,trcenter]=trafo(REW.cov,REW.loc,med,madY);

REW.cor=REW.cov./(sqrt(diag(REW.cov))*sqrt(diag(REW.cov))');

if det(trcov) < tolMCD
    
    [center,covar,z,~,~,count]=exactfit(Y,NaN,med,madY,z,REW.loc,REW.cov,n);
    REW.plane=z;
    
    correl=covar./(sqrt(diag(covar))*sqrt(diag(covar))');
    
    displREW(count,n,v,center,covar,'Minimum Covariance Determinant Estimato',z,correl);
    
    REW.weights=weights;
    REW.md=RAW.md;
else
    % Remark: given that Y has been standardized it is necessary to use
    % location and covariance estimates based on standardized data
    md=mahalFS(Y,REW.loc,REW.cov);
    
    if betathresh==1
        m=sum(weights);
        betain=((m-1)^2/m)*betainv(conflev,0.5*v,0.5*(m-v-1));
        % unitsin =units which have weight equal to 1
        
        md1=[md seq'];
        
        unitsin=seq(weights==1);
        mdsel=md1(unitsin,:);
        outbeta=mdsel(mdsel(:,1)>betain,2);
        
        fout=((m+1)/m)*(m-1)*v/(m-v)*finv(conflev,v,m-v);
        unitsout=seq(weights==0);
        
        mdsel=md1(unitsout,:);
        outf=mdsel(mdsel(:,1)>fout,2);
        
        outliers=[outbeta;outf];
        REW.outliers=outliers';
    else
        REW.outliers=seq(md > thresh);
    end
    
    REW.weights = weights;
    % Store reweighted Mahalanobis distances (in squared units)
    REW.md=md;
end

% Store reweighted estimates of location and covariance on the original
% scale
REW.cov=trcov;
REW.loc=trcenter;
REW.class='mcdr';

plo=options.plots;

% Plot Mahalanobis distances with outliers highlighted
if isstruct(plo) || (~isstruct(plo) && plo~=0)
    
    laby='Raw MCD Mahalanobis distances';
    malindexplot(RAW.md,v,'conflev',conflev,'laby',laby,'numlab',RAW.outliers,'tag','rawmcd');
    
    figure('Tag','pl_spm_outliers');
    group=ones(n,1);
    if ~isempty(RAW.outliers)
        group(RAW.outliers)=2;
    end
    spmplot(Y,group,plo);
    set(gcf,'Name',' Raw MCD: scatter plot matrix with outliers highlighted');
    
    laby='Reweighted MCD Mahalanobis distances';
    malindexplot(REW.md,v,'conflev',conflev,'laby',laby,'numlab',REW.outliers,'tag','rewmcd');
    
    figure('Tag','pl_spm_outliers');
    group=ones(n,1);
    if ~isempty(REW.outliers)
        group(REW.outliers)=2;
    end
    spmplot(Y,group,plo);
    set(gcf,'Name',' Reweighted MCD: scatter plot matrix with outliers highlighted');
    
end

% If requested, save data matrix Y into the output structure. This is used
% by function malindexplot.
if options.ysaveRAW
    RAW.Y = Y;
end
if options.ysaveREW
    REW.Y = Y;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subfunctions called by main function mcd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% subfunction IRWLSmcd
    function outIRWLS = IRWLSmcd(Y, initialloc, initialcov, h, refsteps, reftol, restrfactor)
        %IRWLSmult (iterative reweighted least squares) does refsteps refining steps from initialloc
        % for refsteps times or till convergence.
        %
        %  Required input arguments:
        %
        %    Y: Data matrix containining n observations on v variables.
        %       Rows of Y represent observations, and columns represent variables.
        % initialloc   : v x 1 vector containing initial estimate of location
        % initialcov   : v x v initial estimate of covariance matrix
        %   refsteps  : scalar, number of refining (IRLS) steps
        %   reftol    : relative convergence tolerance for the fully iterated
        %               best candidates. Deafult value is 1e-7
        %
        %  Optional input arguments:
        %
        %
        %  Output:
        %
        %  The output consists of a structure 'outIRWLS' containing:
        %      outIRWLS.loc     : v x 1 vector. Estimate of location after refsteps
        %                         refining steps.
        %      outIRWLS.cov     : v x v matrix. Estimate of the shape matrix after
        %                         refsteps refining steps.
        %      outIRWLS.obj     : scalar. Value of the objective function after refsteps refining
        %                         steps.
        %      outIRWLS.weights : n x 1 vector. Weights assigned to each
        %                         observation. h units have weights equal to 1 and
        %                         n-h units have weights equal to 0
        %
        % In the IRWLS procedure the value of loc and the value of the scale and
        % of the shape matrix are updated in each step
        
        loc = initialloc;
        % Mahalanobis distances from initialloc and Initialshape
        mahaldist = sqrt(mahalFS(Y, initialloc, initialcov));
        
        iter = 0;
        locdiff = 9999;
        if isfinite(restrfactor)
            userestrfactor=true;
        else
           userestrfactor=false;
        end

        while ( (locdiff > reftol) && (iter < refsteps) )
            iter = iter + 1;
            
            %{
            [~,sortdist]=sort(mahaldist);
            %obs_in_set = sort(sortdist(1:h)) ;  % sort removed: it is not necessary
            obs_in_set  = sortdist(1:h) ;
            %}
            
            hm = quickselectFS(mahaldist,h);
            obs_in_set = mahaldist<=hm;

            newloc      = mean(Y(obs_in_set,:));
            newcov      = cov(Y(obs_in_set,:));

            % Apply restriction factor to the covariance matrix of subset
    if userestrfactor==true
        [V,eige]=eig(newcov);
        eigenew=restreigen(diag(eige),h,restrfactor);
        newcov=V*diag(eigenew)*V';
    end
            
            obj         = det(newcov);
            
            % Compute MD
            mahaldist = sqrt(mahalFS(Y,newloc,newcov));
            
            % locdiff is linked to the tolerance
            locdiff = norm(newloc-loc,1)/norm(loc,1);
            loc = newloc;
            
        end
        weights=zeros(size(Y,1),1);
        weights(obs_in_set)=1;
        
        outIRWLS = struct('loc',newloc,'cov',newcov,'obj',obj,'weights',weights);
        % assignments below take twice the time of line above
        %outIRWLS.loc = newloc;
        %outIRWLS.cov = newcov;
        %outIRWLS.obj = obj;
        %outIRWLS.weights=weights;
    end

%% displRAW function
    function [raw_center,raw_cov,center,covar,raw_objective,raw_wt,mcd_wt]=displRAW(exactfit,...
            count,weights,n,p,center,covar,method,z,varargin)
        % exactfit = scalar which
        % exactfit=1 ==> The covariance matrix of the data is singular
        % exactfit=2 ==> The covariance matrix has become singular during the iterations of the MCD algorithm
        % exactfit=3 ==> The %g-th order statistic of the absolute deviation of variable %g is zero.
        
        % Determines some fields of the output argument REW for the exact fit situation.  It also
        % displays and writes the messages concerning the exact fit situation.  If the raw MCD
        % covariance matrix is not singular but the reweighted is, then the function displrw is
        % called instead of this function.
        
        [raw_center,center]=deal(center);
        [raw_cov,~]=deal(covar);
        raw_objective=0;
        mcd_wt=weights;
        raw_wt=weights;
        
        switch exactfit
            case 1
                msg='The covariance matrix of the data is singular.';
            case 2
                msg='The covariance matrix has become singular during the iterations of the MCD algorithm.';
            case 3
                msg=sprintf('The %g-th order statistic of the absolute deviation of variable %g is zero. ',varargin{1},varargin{2});
        end
        
        msg=sprintf([msg '\nThere are %g observations in the entire dataset of %g observations that lie on the \n'],count,n);
        switch p
            case 2
                msg=sprintf([msg 'line with equation %g(x_i1-m_1)%+g(x_i2-m_2)=0 \n'],z);
                msg=sprintf([msg 'where the mean (m_1,m_2) of these observations is the MCD location']);
            case 3
                msg=sprintf([msg 'plane with equation %g(x_i1-m_1)%+g(x_i2-m_2)%+g(x_i3-m_3)=0 \n'],z);
                msg=sprintf([msg 'where the mean (m_1,m_2,m_3) of these observations is the MCD location']);
            otherwise
                msg=sprintf([msg 'hyperplane with equation a_1 (x_i1-m_1) + ... + a_p (x_ip-m_p) = 0 \n']);
                msg=sprintf([msg 'with coefficients a_i equal to : \n\n']);
                msg=sprintf([msg sprintf('%g  ',z)]);
                msg=sprintf([msg '\n\nand where the mean (m_1,...,m_p) of these observations is the MCD location']);
        end
        
        % method=strvcat(method,[msg '.']);
        method=char(method,[msg '.']);
        
        disp(method)
        
    end

%% displREW function
    function method=displREW(count,n,p,center,covar,method,z,correl)
        % Displays and writes messages in the case the reweighted robust covariance
        % matrix is singular.
        
        msg=sprintf('The reweighted MCD scatter matrix is singular. \n');
        msg=sprintf([msg 'There are %g observations in the entire dataset of %g observations that lie on the\n'],count,n);
        
        switch p
            case 2
                msg=sprintf([msg 'line with equation %g(x_i1-m_1)%+g(x_i2-m_2)=0 \n\n'],z);
                msg=sprintf([msg 'where the mean (m_1,m_2) of these observations is : \n\n']);
            case 3
                msg=sprintf([msg 'plane with equation %g(x_i1-m_1)%+g(x_i2-m_2)%+g(x_i3-m_3)=0 \n\n'],z);
                msg=sprintf([msg 'where the mean (m_1,m_2,m_3) of these observations is : \n\n']);
            otherwise
                msg=sprintf([msg 'hyperplane with equation a_1 (x_i1-m_1) + ... + a_p (x_ip-m_p) = 0 \n']);
                msg=sprintf([msg 'with coefficients a_i equal to : \n\n']);
                msg=sprintf([msg sprintf('%g  ',z)]);
                msg=sprintf([msg '\n\nand where the mean (m_1,...,m_p) of these observations is : \n\n']);
        end
        
        msg=sprintf([msg sprintf('%g  ',center)]);
        msg=sprintf([msg '\n\nTheir covariance matrix equals : \n\n']);
        msg=sprintf([msg sprintf([repmat('% 13.5g ',1,p) '\n'],covar)]);
        
        msg=sprintf([msg '\n\nand their correlation matrix equals : \n\n']);
        msg=sprintf([msg sprintf([repmat('% 13.5g ',1,p) '\n'],correl)]);
        
        method=char(method,[msg '.']);
        disp(method)
        
    end

%% mcduni function
    function [initmean,initcov] = mcduni(y,h,alpha,smallsamplecorfactor)
        
        ncas=length(y);
        len=ncas-h+1;
        
        % The exact MCD algorithm for the univariate case.
        y=sort(y);
        
        % Initialize ay and aq
        ay=zeros(len,1);
        sq=ay;
        
        ay(1)=sum(y(1:h));
        
        for samp=2:len
            ay(samp)=ay(samp-1)-y(samp-1)+y(samp+h-1);
        end
        
        ay2=ay.^2/h;
        
        sq(1)=sum(y(1:h).^2)-ay2(1);
        
        for samp=2:len
            sq(samp)=sq(samp-1)-y(samp-1)^2+y(samp+h-1)^2-ay2(samp)+ay2(samp-1);
        end
        
        sqmin=min(sq);
        ijk=find(sq==sqmin);
        ndup=length(ijk);
        slutn(1:ndup)=ay(ijk);
        initmean=slutn(floor((ndup+1)/2))/h;
        
        if smallsamplecorfactor==true
            c1factor = corfactorRAW(1,ncas,alpha);
        end
        
        c1factor = c1factor * consistencyfactor(h,ncas,1);
        initcov  = c1factor * sqmin/(h-1);
    end

%% consistencyfactor function
    function rawconsfac = consistencyfactor(h,n,v,nu)
        % The consistency factor is used to take the effect of trimming
        % into account. 
        if nargin<4
            % This is the standard case, applied when uncontaminated data
            % are assumed to come from a multivariate Normal model.

            a=chi2inv(h/n,v);
            rawconsfac=(h/n)/(chi2cdf(a,v+2));
        else
            % This is the case of a heavy-tail scenario, when
            % uncontaminated data come from a multivariate Student-t
            % distribution. From Barabesi et al. (2023), Trimming
            % heavy-tailed multivariate data, submitted.

            alpha = (n-h)/n;
            integrand = @(u) 1 / (1 - betainv(u,v/2,nu/2));
            theintegral = integral(integrand,0,alpha);
            rawconsfac = ((nu-2) / (alpha*v) * theintegral - (nu - 2)/v)^(-1);
        end
    end

%% corfactorRAW function
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
        rawcorfac=1/fp_alpha_n;
        if rawcorfac <=0 || rawcorfac>50
            rawcorfac=1;
            disp('Warning: problem in subfunction corfactorRAW')
            disp(['Correction factor for covariance matrix based on simulations found =' num2str(rawcorfac)])
            disp('Given that this value is clearly wrong we put it equal to 1 (no correction)')
            disp('This may happen when n is very small and p is large')
        end
    end

%% corfactorREW function
    function rewcorfac=corfactorREW(p,n,alpha)
        
        if p > 2
            coeffrewqpkwad875=[-0.544482443573914,1.25994483222292,2;-0.343791072183285,1.25159004257133,3]';
            coeffrewqpkwad500=[-1.02842572724793,1.67659883081926,2;-0.26800273450853,1.35968562893582,3]';
            y1_500=1+(coeffrewqpkwad500(1,1)*1)/p^coeffrewqpkwad500(2,1);
            y2_500=1+(coeffrewqpkwad500(1,2)*1)/p^coeffrewqpkwad500(2,2);
            y1_875=1+(coeffrewqpkwad875(1,1)*1)/p^coeffrewqpkwad875(2,1);
            y2_875=1+(coeffrewqpkwad875(1,2)*1)/p^coeffrewqpkwad875(2,2);
            y1_500=log(1-y1_500);
            y2_500=log(1-y2_500);
            y_500=[y1_500;y2_500];
            A_500=[1,log(1/(coeffrewqpkwad500(3,1)*p^2));1,log(1/(coeffrewqpkwad500(3,2)*p^2))];
            coeffic_500=A_500\y_500;
            y1_875=log(1-y1_875);
            y2_875=log(1-y2_875);
            y_875=[y1_875;y2_875];
            A_875=[1,log(1/(coeffrewqpkwad875(3,1)*p^2));1,log(1/(coeffrewqpkwad875(3,2)*p^2))];
            coeffic_875=A_875\y_875;
            fp_500_n=1-(exp(coeffic_500(1))*1)/n^coeffic_500(2);
            fp_875_n=1-(exp(coeffic_875(1))*1)/n^coeffic_875(2);
        else
            if p == 2
                fp_500_n=1-(exp(3.11101712909049)*1)/n^1.91401056721863;
                fp_875_n=1-(exp(0.79473550581058)*1)/n^1.10081930350091;
            end
            if p == 1
                fp_500_n=1-(exp(1.11098143415027)*1)/n^1.5182890270453;
                fp_875_n=1-(exp(-0.66046776772861)*1)/n^0.88939595831888;
            end
        end
        if 0.5 <= alpha && alpha <= 0.875
            fp_alpha_n=fp_500_n+(fp_875_n-fp_500_n)/0.375*(alpha-0.5);
        end
        if 0.875 < alpha && alpha < 1
            fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
        end
        rewcorfac=1/fp_alpha_n;
        if rewcorfac <=0 || rewcorfac>50
            rewcorfac=1;
            disp('Warning: problem in subfunction corfactorREW');
            disp(['Correction factor for covariance matrix based on simulations found =' num2str(rewcorfac)]);
            disp('Given that this value is clearly wrong we put it equal to 1 (no correction)');
            disp('This may happen when n is very small and p is large');
        end
    end

%% trafo function
    function [covmat,meanvct]=trafo(covmat,meanvct,med,mad)
        
        % Transforms a mean vector and a covariance matrix to the original units.
        % Transform location
        meanvct=meanvct.*mad+med;
        
        %Transform covariance
        % ALTERNATIVE STATEMENT
        % covmat=covmat.*repmat(mad,size(covmat,1),1).*repmat(mad',1,size(covmat,1));
        covmat=bsxfun(@times, covmat, mad);
        covmat=bsxfun(@times, covmat, mad');
        
    end

%% exactfit function
    function [initmean,initcov,z,correl,varargout]=exactfit(Y,plane,med,mad,z,varargin)
        
        % This function is called in the case of an exact fit. It computes the
        % correlation matrix and transforms the coefficients of the hyperplane, the
        % mean, the covariance and the correlation matrix to the original units.
        
        if isnan(plane)
            [meanvct,covmat,n] = deal(varargin{:});
            sigma = eigs_sigma(covmat);
            [z, ~] = eigs(covmat,1,sigma,struct('disp',0));
            dist   = abs(sum((Y-repmat(meanvct,n,1))'.*repmat(z,1,n)));
            plane  = find(dist < generictol);
            varargout{1} = plane;
            varargout{2} = length(plane);
        end
        
        z=z./mad';
        [initcov,initmean]=trafo(cov(Y(plane,:)),mean(Y(plane,:)),med,mad);
        correl=initcov./(sqrt(diag(initcov))*sqrt(diag(initcov))');
    end

%% rockecs function
    function [fachr,dfhr,dfas]=rockecs(n,v,h)
        ratio  = h/n;
        
        qalpha = chi2inv(h/n,v);
        prob   = chi2cdf(qalpha,v+2);
        
        fachr  = prob/ratio;
        calpha = 1/fachr;
        
        alpha  = 1-ratio;
        c2     = -prob/2;
        prob   = chi2cdf(qalpha,v+4);
        c3     = -prob/2;
        c4     = c3*3;
        b1     = calpha*(c3-c4)/ratio;
        b2     = 0.5+calpha/ratio*(c3-qalpha/v*(c2+ratio/2));
        v1     = ratio*b1*b1*(alpha*(calpha*qalpha/v-1)^2-1)...
            -2*c3*calpha*calpha*(3*((b1-v*b2)^2)+(v+2)*b2*(2*b1-v*b2));
        v2     = n*calpha*calpha*(b1*(b1-v*b2)*ratio)^2;
        
        % Asymptotic degrees of freedom
        dfas   = 2/(calpha*calpha*v1/v2);
        xb1    = 0.725;
        xb2    = 0.00663;
        xb3    = 0.0780;
        dfhr   = dfas*exp(xb1-xb2*v-xb3*log(n));
    end
end
%FScategory:MULT-Multivariate
