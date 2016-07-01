function [RAW,REW,varargout] = mve(Y,varargin)
%mve computes Minimum volume ellipsoid
%
%<a href="matlab: docsearchFS('mve')">Link to the help function</a>
%
%  Required input arguments:
%
%    Y: Data matrix containining n observations on v variables.
%       Rows of Y represent observations, and columns represent variables.
%       Missing values (NaN's) and infinite values (Inf's) are allowed,
%       since observations (rows) with missing or infinite values will
%       be excluded from the computations.
%
%  Optional input arguments:
%
%      bdp    : scalar. Breakdown point. (Number between 0
%               and 0.5). The default value is 0.5.
%               Example - 'bdp',1/4 
%               Data Types - double
%      nsamp  : scalar. Number of subsamples of size v which have
%               to be extracted (if not given, default = 500).
%               Example - 'nsamp',10000 
%               Data Types - double
%    refsteps : Number of refining iterations. Scalar. Number of refining iterationsin each
%               subsample (default = 3).
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'refsteps',0 
%                 Data Types - single | double
%     reftol  : scalar. Default value of tolerance for the refining steps.
%               The default value is 1e-6;
%                 Example - 'reftol',1e-8 
%                 Data Types - single | double
%     conflev : Scalar. Number between 0 and 1 containing confidence level which is
%               used to declare units as outliers.
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%               Example - 'conflev',0.99
%               Data Types - double
%      nocheck: Scalar. If nocheck is equal to 1 no check is performed on
%               matrix Y. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%       plots : Plot on the screen. Scalar or structure.
%               If plots is a structure or scalar equal to 1, generates: 
%               (1) a plot of Mahalanobis distances against index number. The
%               confidence level used to draw the confidence bands for
%               the MD is given by the input option conflev. If conflev is
%               not specified a nominal 0.975 confidence interval will be
%               used.
%               (2) a scatter plot matrix with the outliers highlighted. 
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
%        msg  : scalar. Display or not messages
%               on the screen. If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the final
%               estimator else no message is displayed on the screen.
%               Example - 'msg',1
%               Data Types - double

%    ysaveRAW : scalar that is set to 1 to request that the data matrix Y
%               is saved into the output structure RAW. This feature is
%               meant at simplifying the use of function malindexplot.
%               Default is 0, i.e. no saving is done.
%               Example - 'ysaveRAW',1
%               Data Types - double
%    ysaveREW : scalar that is set to 1 to request that the data matrix Y
%               is saved into the output structure REW. This feature is
%               meant at simplifying the use of function malindexplot.
%               Default is 0, i.e. no saving is done.
%               Example - 'ysaveREW',1
%               Data Types - double
%  Output:
%
%  The output consists of two structures RAW and REW. RAW refers to raw
%  mve, on the other hand, REW refers to reweighted mve
%
%         RAW:   structure which contains the following fields
%         RAW.loc  = 1 x v  vector containing raw MCD location of the data
%         RAW.cov  = robust MCD estimate of
%                    covariance matrix. It is the raw MCD covariance matrix
%                    (multiplied by a finite sample correction factor and
%                    an asymptotic consistency factor).
%         RAW.cor  = The raw MVE correlation matrix
%           RAW.obj= The value of the objective function which has been minimized.
%           RAW.bs = (v+1) x 1 vector containing the units forming best subset
%                    associated with MVE estimate of location.
%          RAW.md  = n x 1 vector containing the estimates of the robust
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
%                    These weights determine which are the h observations are used to
%                    compute the final MVE estimates.
%            RAW.h = number of observations which have determined MVE.
%            RAW.Y = Data matrix Y. The field is present if option
%                    ysaveRAW is set to 1.
%        RAW.class = 'mve'
%
%
%         REW : structure which contains the following fields:
%      REW.loc     = The robust location of the data, obtained after reweighting, if
%                    the RAW.cov  is not singular.  Otherwise the raw MVE center is
%                    given here.
%       REW.cov    = The robust covariance matrix, obtained after reweighting and
%                    multiplying with a finite sample correction factor and an asymptotic
%                    consistency factor, if the raw MVE is not singular.  Otherwise the
%                    raw MVE covariance matrix is given here.
%       REW.cor    = The robust correlation matrix, obtained after reweighting
%       REW.md     = The distance of each observation to the final,
%                    reweighted MVE center of the data, relative to the
%                    reweighted MVE scatter of the data.  These distances allow
%                    us to easily identify the outliers. If the reweighted MVE
%                    is singular, RAW.md is given here.
%     REW.outliers = A vector containing the list of the units declared as
%                    outliers after reweighting.
%            REW.Y = Data matrix Y. The field is present if option
%                    ysaveRAW is set to 1.
%       REW.class = 'mve'
%            
%  Optional Output:
%
%            C     : matrix of size nsamp-by-v which contains the indices
%                    of the subsamples extracted for
%                    computing the estimate.
%
% More About:
%
% For each subset $J$ of $v+1$ observations
% $\mu_J$ and $C_J$ are the centroid and the covariance matrix based on
% subset $J$.
% 
%
% Rousseeuw and Leroy (RL) (eq. 1.25 chapter 7, p. 259) write the objective
% function for subset $J$ as
% \[
% RL_J=\left( med_{i=1, ..., n} \sqrt{ (y_i -\mu_J)' C_J^{-1} (y_i -\mu_J) } \right)^v \sqrt|C_J|
% \]
%
% Maronna Martin and Yohai (MMY), eq. (6.57), define $\Sigma_J = C_j /
% |C_j|^{1/v}$ and write the objective function for subset $J$ as
% \[
% MMY_J =  \hat \sigma \left( (y_i -\mu_J)' \Sigma_J^{-1} (y_i -\mu_J) \right) |C_J|^{1/v}
%       =  \hat \sigma \left( (y_i -\mu_J)' C_J^{-1} (y_i -\mu_J) \right) |C_J|^{1/v}
% \]
% where $\hat \sigma \left( (y_i -\mu_J)' C_J^{-1} (y_i -\mu_J) \right) = med_{i=1, ..., n}(y_i -\mu_J)' C_J^{-1} (y_i -\mu_J)$. 
% Note that $MMY_J= (RL)^{2/v}$.
%
%   To RAW.cov a consistency factor has been applied which is based on
%   chi2inv(1-bdp,v). On the other hand to REW.cov the usual asymptotic
%   consistency factor is applied. In this case we have used the empirical
%   percentage of trimming that is the ratio hemp/n where hemp is the
%   number of units which had a MD smaller than the cutoff level determined
%   by thresh=chi2inv(conflev,v); MD are computed using RAW.loc and
%   RAW.cov.
%
%
% The mve method is intended for continuous variables, and assumes that
% the number of observations n is at least 5 times the number of variables v.%
%
% See also: mcd.m
%
% References:
%
%   Rousseeuw, P.J. (1984), "Least Median of Squares Regression", Journal
%   of the American Statistical Association, Vol. 79, pp. 871-881.
%   Rousseeuw, P.J. and Leroy A.M. (1987), Robust regression and outlier
%   detection,  Wiley New York.
%
%  
%
% Acknowledgements:
%
% This function follows the lines of MATLAB/R code developed during the
% years by many authors.
% For more details see http://www.econ.kuleuven.be/public/NDBAE06/programs/
% and the R library robustbase http://robustbase.r-forge.r-project.org/
% The core of these routines, e.g. the resampling approach, however, has
% been completely redesigned, with considerable increase of the
% computational performance.
%
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mve')">Link to the help page for this function</a>
% Last modified 31-05-2016
%
% Examples:

%{
    % mve with all default options.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    RAW=mve(Ycont);
%}

%{
    %% mve with optional arguments.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    RAW=mve(Ycont,'plots',1);
%}

%{
    % mve monitoring the reweighted estimates.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [RAW,REW]=mve(Ycont);
%}

%{
    %% mve monitoring the extracted subsamples.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [RAW,REW,C]=mve(Ycont);
%}

%% Beginning of code

nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);
[n,v]=size(Y);
seq=1:n;

%% User options

% If the number of all possible subsets is <10000 the default is to extract
% all subsets otherwise just 10000.
% Notice that we use bc, a fast version of nchoosek. One may also use the
% approximation floor(exp(gammaln(n+1)-gammaln(n-p+1)-gammaln(p+1))+0.5)
ncomb=bc(n,v+1);
nsampdef=min(10000,ncomb);
% default value of tolerance for the refining steps convergence for  each
% extracted subset
reftoldef=1e-6;

bdpdef=0.5;

options=struct('nsamp',nsampdef,'bdp',bdpdef,...
    'plots',0,'nocheck',0,'conflev',0.975,'msg',1,...
    'ysaveRAW',0,'ysaveREW',0,'refsteps',0,'reftol',reftoldef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    % And check if the optional user parameters are reasonable.
    
    if options.bdp <0 || options.bdp >0.5
        error('FSDA:mve:WrongBdp','bdp must be a number greater than 0 and smaller or equal than 0.5')
    end
    
    
    % Check number of subsamples to extract
    if options.nsamp>ncomb
        disp('Number of subsets to extract greater than (n v+1). It is set to (n  +1)');
        options.nsamp=0;
    elseif  options.nsamp<0
        error('FSDA:mve:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
    end
end

% Default values for the optional
% parameters are set inside structure 'options'

bdp=options.bdp;
nsamp=options.nsamp;
h=floor(2*floor((n+v+1)/2)-n+2*(n-floor((n+v+1)/2))*(1-bdp));

% hmin is the minimum number of observations whose covariance determinant
% will be minimized.
hmin=floor(2*floor((n+v+1)/2)-n+2*(n-floor((n+v+1)/2))*(0.5));

if h < hmin
    error('FSDA:mve:Wrongh',['The MVE must cover at least ' int2str(hmin) ' observations.'])
elseif h > n
    error('FSDA:mve:Wrongh','h is greater than the number of non-missings and non-infinites.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% singsub = scalar which will contain the number of singular subsets which
% are extracted (that is the subsets of size v+1 which are not full rank)
singsub=0;

msg=options.msg;            % Scalar which controls the messages displayed on the screen

z=zeros(1,v);
weights=zeros(1,n);

conflev = options.conflev;
thresh=chi2inv(conflev,v);

refsteps = options.refsteps;    % refining steps
reftol = options.reftol;        % tolerance for refining steps


%% Standardization of the data with medians and mads
% The standardization of the data will now be performed.
med=median(Y);
% OLD mad=sort(abs(Y-repmat(med,n,1)));
mad = sort(abs(bsxfun(@minus,Y, med)));

mad=mad(h,:);
ii=find(mad < eps, 1 );
if ~isempty(ii)
    % The h-th order statistic is zero for the ii-th variable. The array plane contains
    % all the observations which have the same value for the ii-th variable.
    plane=find(abs(Y(:,ii)-med(ii)) < eps)';
    meanplane=mean(Y(plane,:));
    weights(plane)=1;
    if v==1
        out.weights=weights;
        out.loc=meanplane;
        [out.cov,out.objective]=deal(0);
        REW.method=sprintf('\nUnivariate location and scale estimation.');
        REW.method=char(REW.method,sprintf('%g of the %g observations are identical.',length(plane),n));
        disp(REW.method);
    else
        z(ii)=1;
        REW.plane=z;
        covplane=cov(Y(plane,:));
        [RAW.center,RAW.cov,REW.center,REW.cov,RAW.objective,RAW.weights,REW.weights,...
            REW.method]=displRAW(3,length(plane),weights,n,v,meanplane,covplane,'MVE',z,h,ii);
    end
    return
end
% Standardization of the data with location estimate (median) and scale
% estimate (mad)
Y = bsxfun(@minus,Y, med);
Y = bsxfun(@rdivide, Y, mad);


%% Combinatorial part to extract the subsamples
[C,nselected] = subsets(nsamp,n,v+1,ncomb,msg);
% Store the indices in varargout
if nargout==3
    varargout={C};
end

% volmin will contain the minimum value of objective funtion (volume of the ellipsoid)
volmin=Inf;

% initialise and start timer.
tsampling = ceil(min(nselected/100 , 1000));
time=zeros(tsampling,1);

for i=1:nselected
    if i <= tsampling, tic; end
    
    % Extract i-th row of matrix C
    s = C(i,:);
    
    A=Y(s,:);
    mA=mean(A);
    
    %Chapter 7, eq. 1.23
    Cj=cov(A);
    
    detCj=det(Cj);
    
    if  detCj> exp(-5*v)
        
        % Compute Mahalanobis distance
        md=mahalFS(Y,mA,Cj);
        
      
        if refsteps==0
            % Order them
            mdsqsor=sort(md);
        else
            
            iter = 0;
            locdiff = 9999;
            oldmA=mA;
            
            while ( (locdiff > reftol) && (iter < refsteps) )
                [~,sortdist]=sort(md);
                obs_in_set=sort(sortdist(1:h));
                mA=mean(Y(obs_in_set,:));
                Cj=cov(Y(obs_in_set,:));
                % Compute Mahalanobis distances
                md=mahalFS(Y,mA,Cj);
                % Order them
                mdsqsor=sort(md);
                
                % locdiff is linked to the tolerance
                locdiff = norm(mA-oldmA,1)/norm(oldmA,1);
                oldmA = mA;
                
                mj2=mdsqsor(h);
                mj=sqrt(mj2);
                vol=sqrt(detCj)*(mj^v);
                disp([iter vol])
                iter=iter+1;
            end
            
        end
        
        % Take the MD which occupies the h-ordered position
        mj2=mdsqsor(h);
        
        mj=sqrt(mj2);
        
        % Chapter 7, eq. 1.25 of Rousseeuw and Leroy (1987) - The objective function
        vol=sqrt(detCj)*(mj^v);
        
        
        if vol<volmin
            volmin=vol;
            
            % Store units forming best subset (which generated MVE)
            sbest=s;
            
        end
    else
        singsub=singsub+1;
    end
    
    if i <= tsampling
        
        % sampling time until step tsampling
        time(i)=toc;
    elseif i==tsampling+1
        % stop sampling and print the estimated time
        if msg==1
            fprintf('Total estimated time to complete MVE: %5.2f seconds \n', nselected*median(time));
        end
    end
    
end

if volmin==Inf
    error('FSDA:mve:NoFullRank','No subset had full rank. Please increase the number of subsets or check your design matrix X')
else
end

RAW=struct;
A=Y(sbest,:);

% Store number of observations which have determined MVE
RAW.h=h;

% Store best subset associated with MVE
RAW.bs=sbest;

mA=mean(A);

%Chapter 7, eq. 1.23 of Rousseeuw and Leroy (1977) p. 259
Cj=cov(A);

md=mahalFS(Y,mA,Cj);
mdsqsor=sort(md);
% Take the md which occupies the h-th ordered position
mj2=mdsqsor(h);

% weights vector will contain the units of the sample which have MD from
% RAW.loc and RAW.cov smaller than the h-th ordered distance
weights=md<=mj2;
RAW.weights=weights;

% Final estimate of covariance matrix
% See equation 1.26 of Rousseeuw and Leroy (1977) p. 259
RAW.cov=mj2*Cj/chi2inv(1-bdp,v);

% Compute final estimate of objective function in the original scale
RAW.obj=volmin*(prod(mad)^2);

% Recompute md using final estimate of cov matrix
% in the transformed scale
md=mahalFS(Y,mA,RAW.cov);

% Store raw Mahalanobis distances
RAW.md=md;

% Find final estimate of location and covariance in the original scale
[RAW.cov,RAW.loc]=trafo(RAW.cov,mA,med,mad);

% Store robust estimate of correlation matrix
RAW.cor=RAW.cov./(sqrt(diag(RAW.cov))*sqrt(diag(RAW.cov))');

% Find the outliers using md based on transformed data (Y), transformed
% location (mA) and transformed final estimate of covariance (RAW.cov before
% trafo)
RAW.outliers=seq(md > thresh);

% Store confidence level
RAW.conflev=conflev;

% Store total number of singular subsets
RAW.singsub=singsub;

% Store class
RAW.class='mve';

% Use the units not declared as outliers are used to form another estimate of location and covariance
weights=md<=thresh;
locrw=mean(Y(weights==1,:));
covrw=cov(Y(weights==1,:));
REW.loc=locrw;
REW.cov=covrw;

% Apply consistency factor to reweighted estimate of covariance
hrew=sum(weights);
if hrew<n
    factor=consistencyfactor(hrew,n,v);
else
    factor=1;
end
REW.cov=factor*REW.cov;

% Recompute md using final estimate of cov matrix
% in the transformed scale
md=mahalFS(Y,REW.loc,REW.cov);

% Store Mahalanobis distances from reweighted estimate of centroid and
% reweighted estimate of covariance matrix
REW.md=md;
REW.outliers=seq(md > thresh);

% Compute final reweighted estimate of covariance and location in the
% original scale
[REW.cov,REW.loc]=trafo(REW.cov,REW.loc,med,mad);

% Store final reweighted estimate of correlation matrix
REW.cor=REW.cov./(sqrt(diag(REW.cov))*sqrt(diag(REW.cov))');

REW.class='mver';

plo=options.plots;

% Plot Mahalanobis distances with outliers highlighted
if isstruct(plo) || (~isstruct(plo) && plo~=0)
    
    laby='Raw MVE Mahalanobis distances';
    malindexplot(RAW.md,v,'conflev',conflev,'laby',laby,'numlab',RAW.outliers);
    
    figure('Tag','pl_spm_outliers');
    group=ones(n,1);
    if ~isempty(RAW.outliers)
        group(RAW.outliers)=2;
    end
    spmplot(Y,group,plo);
    set(gcf,'Name',' Raw MVE: scatter plot matrix with outliers highlighted');
    
    laby='Reweighte MVE Mahalanobis distances';
    malindexplot(REW.md,v,'conflev',conflev,'laby',laby,'numlab',REW.outliers);
    
    figure('Tag','pl_spm_outliers');
    group=ones(n,1);
    if ~isempty(REW.outliers)
        group(REW.outliers)=2;
    end
    spmplot(Y,group,plo);
    set(gcf,'Name',' Reweighted MVE: scatter plot matrix with outliers highlighted');
end

% If requested, save data matrix Y into the output structure. This is used
% by function malindexplot.
if options.ysaveRAW
    RAW.Y = Y;
end
if options.ysaveREW
    REW.Y = Y;
end

end


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


function rawconsfac=consistencyfactor(h,n,v)
a=chi2inv(h/n,v);
rawconsfac=(h/n)/(chi2cdf(a,v+2));
end



%% displRAW function
function [raw_center,raw_cov,center,covar,raw_objective,raw_wt,mcd_wt,method]=displRAW(exactfit,...
    count,weights,n,p,center,covar,method,z,varargin)
% exactfit = scalar which
% exactfit=1 ==> The covariance matrix of the data is singular
% exactfit=2 ==> The covariance matrix has become singular during the iterations of the MVE algorithm
% exactfit=3 ==> The %g-th order statistic of the absolute deviation of variable %g is zero.

% Determines some fields of the output argument REW for the exact fit situation.  It also
% displays and writes the messages concerning the exact fit situation.  If
% the raw MVE covariance matrix is not singular but the reweighted is, then the function displrw is
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
        msg='The covariance matrix has become singular during the iterations of the MVE algorithm.';
    case 3
        msg=sprintf('The %g-th order statistic of the absolute deviation of variable %g is zero. ',varargin{1},varargin{2});
end

msg=sprintf([msg '\nThere are %g observations in the entire dataset of %g observations that lie on the \n'],count,n);
switch p
    case 2
        msg=sprintf([msg 'line with equation %g(x_i1-m_1)%+g(x_i2-m_2)=0 \n'],z);
        msg=sprintf([msg 'where the mean (m_1,m_2) of these observations is the MVE location']);
    case 3
        msg=sprintf([msg 'plane with equation %g(x_i1-m_1)%+g(x_i2-m_2)%+g(x_i3-m_3)=0 \n'],z);
        msg=sprintf([msg 'where the mean (m_1,m_2,m_3) of these observations is the MVE location']);
    otherwise
        msg=sprintf([msg 'hyperplane with equation a_1 (x_i1-m_1) + ... + a_p (x_ip-m_p) = 0 \n']);
        msg=sprintf([msg 'with coefficients a_i equal to : \n\n']);
        msg=sprintf([msg sprintf('%g  ',z)]);
        msg=sprintf([msg '\n\nand where the mean (m_1,...,m_p) of these observations is the MVE location']);
end

% method=strvcat(method,[msg '.']);
method=char(method,[msg '.']);

disp(method)

end

%FScategory:MULT-Multivariate