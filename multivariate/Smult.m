function [out , varargout] = Smult(Y,varargin)
%Smult computes S estimators in multivariate analysis
% Fast S algorithm for multivariate location estimation is based on Tukey's
% biweight function
%
%<a href="matlab: docsearchFS('Smult')">Link to the help function</a>
%
%  Required input arguments:
%
% Y :           Input data. Matrix. 
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
%  Optional input arguments:
%
%         bdp :  breakdown point. Scalar.
%               It measures the fraction of outliers
%               the algorithm should resist. In this case any value greater
%               than 0 but smaller or equal than 0.5 will do fine (default=0.5).
%               Note that given bdp nominal
%               efficiency is automatically determined.
%                 Example - 'bdp',0.4
%                 Data Types - double
%
%       nsamp   : Number of subsamples which will be extracted to find the
%                 robust estimator. Scalar. If nsamp=0 all subsets will be extracted.
%                 They will be (n choose p).
%                 If the number of all possible subset is <1000 the
%                 default is to extract all subsets otherwise just 1000.
%                 Example - 'nsamp',1000 
%                 Data Types - single | double
%
%    refsteps : Number of refining iterations. Scalar. Number of refining iterationsin each
%               subsample (default = 3).
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'refsteps',0 
%                 Data Types - single | double
%
%     reftol  : scalar. Default value of tolerance for the refining steps.
%               The default value is 1e-6;
%                 Example - 'reftol',1e-8 
%                 Data Types - single | double
%
%refstepsbestr: number of refining iterations for each best subset. Scalar.
%               Scalar defining number of refining iterations for each
%               best subset (default = 50).
%                 Example - 'refstepsbestr',10 
%                 Data Types - single | double
%
% reftolbestr : Tolerance for the refining steps. Scalar. 
%               Tolerance for the refining steps
%               for each of the best subsets
%               The default value is 1e-8;
%                 Example - 'reftolbestr',1e-10 
%                 Data Types - single | double
%
%     minsctol: tolerance for the iterative
%               procedure for finding the minimum value of the scale. Scalar. 
%               Value of tolerance for the iterative
%               procedure for finding the minimum value of the scale
%               for each subset and each of the best subsets
%               (It is used by subroutine minscale.m)
%               The default value is 1e-7;
%                 Example - 'minsctol',1e-7 
%                 Data Types - single | double
%
%      bestr  : number of "best betas" to remember. Scalar. Scalar defining
%               number of "best betas" to remember from the subsamples.
%               These will be later iterated until convergence (default=5)
%                 Example - 'bestr',10 
%                 Data Types - single | double
%
%     conflev :  Confidence level which is
%               used to declare units as outliers. Scalar. 
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%                 Example - 'conflev',0.99
%                 Data Types - double
%
%      nocheck : Check input arguments. Scalar. 
%               If nocheck is equal to 1 no check is performed on
%               matrix Y. As default nocheck=0.
%               Example - 'nocheck',1 
%               Data Types - double
%
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
%                       the variables. As default value, the labels which are
%                       added are Y1, ...Yv.
%                 Example - 'plots',0 
%                 Data Types - single | double
%
%        msg  : Level of output to display. Scalar. 
%               If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the final estimator
%               else no message is displayed on the screen
%                 Example - 'msg',0 
%                 Data Types - single | double
%
%       ysave : save input matrix Y. Scalar. 
%               Scalar that is set to 1 to request that the data matrix Y
%               is saved into the output structure out. This feature is
%               meant at simplifying the use of function malindexplot.
%               Default is 0, i.e. no saving is done. 
%               Example - 'ysave',1 
%               Data Types - double
%
%  Output:
%
%  out :     A structure containing the following fields
%
%         out.loc  = 1 x v  vector containing S estimate of location
%         out.shape= v x v matrix containing robust estimate of the shape
%                   matrix. Remark det|shape|=1
%         out.scale= scalar, robust estimate of the scale
%         out.cov  = out.scale^2 * out.shape = robust estimate of
%                   covariance matrix
%           out.bs = (v+1) x 1 vector containing the units forming best subset
%                    associated with S estimate of location.
%          out.md  = n x 1 vector containing the estimates of the robust
%                    Mahalanobis distances (in squared units)
%     out.outliers = A vector containing the list of the units declared as
%                   outliers using confidence level specified in input
%                   scalar conflev
%      out.conflev = Confidence level that was used to declare outliers
%      out.singsub = Number of subsets without full rank. Notice that
%                    out.singsub > 0.1*(number of subsamples) produces a
%                    warning
%      out.weights = n x 1 vector containing the estimates of the weights
%            out.Y = Data matrix Y. The field is present if option ysave
%                    is set to 1.
%        out.class = 'Smult'
%
%  Optional Output:
%
%            C        : matrix containing the indices of the subsamples 
%                       extracted for computing the estimate (the so called
%                       elemental sets).
%
%
% See also: MMmult
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
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
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('Smult')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Smult with all default options.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out]=Smult(Ycont);
%}

%{
    %% Smult with optional arguments.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out]=Smult(Ycont,'plots',1);
%}

%{
    % Smult with exctracted subsamples.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out,C]=Smult(Ycont);
%}

%% Beginning of code

% Input parameters checking
%chkinputM does not do any check if option nocheck=1
nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);

[n,v]=size(Y);

% default value of break down point
bdpdef=0.5;

% default values of subsamples to extract
ncomb=bc(n,v+1);
nsampdef=min(1000,ncomb);

% default value of number of refining iterations (C steps) for each extracted subset
refstepsdef=3;
% default value of tolerance for the refining steps convergence for  each extracted subset
reftoldef=1e-6;
% default value of number of best locs to remember
bestrdef=5;
% default value of number of refining iterations (C steps) for best subsets
refstepsbestrdef=50;
% default value of tolerance for the refining steps convergence for best subsets
reftolbestrdef=1e-8;
% default value of tolerance for finding the minimum value of the scale
% both for each extracted subset and each of the best subsets
minsctoldef=1e-7;

% store default values in the structure options
options=struct('nsamp',nsampdef,'refsteps',refstepsdef,'bestr',bestrdef,...
    'reftol',reftoldef,'minsctol',minsctoldef,...
    'refstepsbestr',refstepsbestrdef,'reftolbestr',reftolbestrdef,...
    'bdp',bdpdef,'plots',0,'conflev',0.975,'nocheck',0,'msg',1,'ysave',0);

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:Smult:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
minsctol = options.minsctol;    % tolerance for finding minimum value of the scale for each subset
msg=options.msg;                %

refstepsbestr=options.refstepsbestr;  % refining steps for the best subsets
reftolbestr=options.reftolbestr;      % tolerance for refining steps for the best subsets

% Find constant c linked to Tukey's biweight
% rho biweight is strictly increasing on [0 c] and constant on [c \infty)
% E(\rho) = kc = (c^2/6)*bdp, being kc the K of Rousseeuw and Leroy
c = TBbdp(bdp,v);
kc = (c^2/6)*bdp;

% Initialize the matrices which contain the best "bestr" estimates of
% location, index of subsets, shape matrices and scales
bestlocs = zeros(bestr, v);
bestsubset = zeros(bestr, v+1);
bestshapes=zeros(v,v,bestr);
bestscales = Inf * ones(bestr,1);
sworst = Inf;

% singsub = scalar which will contain the number of singular subsets which
% are extracted (that is the subsets of size p which are not full rank)
singsub=0;

%% Extract in the rows of matrix C the indexes of all required subsets
[C,nselected] = subsets(nsamp,n,v+1,ncomb,msg);
% Store the indices in varargout
if nargout==2
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
    
    % find a subset of size v+1 in general position (with rank v).
    index = C(i,:);
    
    Yj = Y(index,:);
    ranky = rank(Yj);
    
    if ranky==v
        
        locj = mean(Yj);        % centroid of subset
        Sj = cov(Yj);           % covariance of subset
        Gj = det(Sj)^(-1/v)*Sj; % shape of subset
        
        % Function IRWLSmult performs refsteps steps of IRLS on elemental
        % start. Input:
        % - Y = datamatrix of dimension n x v
        % - locj = row vector containing (robust) centroid
        % - Gj = v x v shape matrix
        % - scale = estimate of the scale (if scale=0). Scaled MAD of Mahalanobis
        %   distances using locj and Gj as centroid and shape matrix is used
        % - refsteps = number of refining iterations
        % - reftol = tolerance for convergence of refining iterations
        % - c = constant of Tukey's biweight linkted to breakdown point
        % - kc = (c^2/6)*bdp
        % Remark: in IRWLSmult the centroid and shape matrix are re-calculated
        % at each step; on the other hand in minscale they are kept fixed.
        outIRWLS = IRWLSmult(Y, locj, Gj, refsteps, reftol, c, kc);
        
        % The output of IRWLSmult is a structure containing centroid, shape
        % matrix and estimate of the scale
        locrw = outIRWLS.loc;
        shaperw = outIRWLS.shape;
        scalerw = outIRWLS.scale;
        
        % Compute Mahalanobis distances using locrw and shaperw
        mdrw = sqrt(mahalFS(Y,locrw,shaperw));
        
        % to find s, save first the best bestr scales and shape matrices
        % (deriving from non singular subsets) and, from iteration bestr+1
        % (associated to another non singular subset), replace the worst scale
        % with a better one as follows
        if ij > bestr
            % from the second step check whether new loc and new shape belong
            % to the top best loc; if so keep loc and shape with
            % corresponding scale.
            if  mean(TBrho(mdrw/sworst,c)) < kc % if >kc skip the sample
                % Find position of the maximum value of bestscale
                [~,ind] = max(bestscales);
                bestscales(ind) = minscale(mdrw,c,kc,scalerw,minsctol);
                bestlocs(ind,:) = locrw;
                bestshapes(:,:,ind) = shaperw;
                sworst = max(bestscales);
                % best subset sssociated with minimum value
                % of the scale
                bestsubset(ind,:)=index;
            end
        else
            bestscales(ij) = minscale(mdrw,c,kc,scalerw,minsctol);
            bestlocs(ij,:) = locrw;
            bestshapes(:,:,ij) = shaperw;
            bestsubset(ij,:) = index;
            ij=ij+1;
        end
    else
        singsub = singsub + 1;
    end
    
    % Write total estimate time to compute final estimate
    if i <= tsampling
        
        % sampling time until step tsampling
        time(i)=toc;
    elseif i==tsampling+1
        % stop sampling and print the estimated time
        if msg==1
            fprintf('Total estimated time to complete S estimator: %5.2f seconds \n', nselected*median(time));
        end
    end
    
end
if singsub==nsamp
    error('FSDA:Smult:NoFullRank','No subset had full rank. Please increase the number of subsets or check your design matrix X')
end

if singsub/nsamp>0.1
    disp('------------------------------')
    disp(['Warning: Number of subsets without full rank equal to ' num2str(100*singsub/nsamp) '%'])
end

% perform C-steps on best 'bestr' solutions, till convergence or for a
% maximum of refstepsbestr steps using a convergence tolerance as specified
% by scalar reftolbestr

% this is to ensure that the condition tmp.scale < superbestscale in the
% next if statement is satisfied at least once
superbestscale = Inf;
for i=1:bestr
    tmp = IRWLSmult(Y,bestlocs(i,:), bestshapes(:,:,i),refstepsbestr,reftolbestr,c,kc, bestscales(i));
    if tmp.scale < superbestscale
        superbestscale  = tmp.scale;
        superbestloc     = tmp.loc;
        superbestshape  = tmp.shape;
        superbestsubset = bestsubset(i,:);
        weights = tmp.weights;
    end
end

out.class   = 'S';
out.loc     = superbestloc;         % robust estimate of location
out.shape   = superbestshape;       % robust estimate of shape matrix
out.scale   = superbestscale;       % robust estimate of the scale
out.cov     = superbestscale^2*superbestshape; %robust estimate of covariance matrix
out.weights = weights;
out.md = mahalFS(Y,out.loc,out.cov);
out.bs=superbestsubset;             % Store units formin best subset

% Store in output structure the outliers found with confidence level conflev
conflev = options.conflev;
seq = 1:n;
out.outliers = seq(out.md > chi2inv(conflev,v) );
out.conflev = conflev;

plo=options.plots;

% Plot Mahalanobis distances with outliers highlighted
if isstruct(plo) || (~isstruct(plo) && plo~=0)
    
    laby='S Mahalanobis distances';
    malindexplot(out.md,v,'conflev',conflev,'laby',laby,'numlab',out.outliers);
    
    figure('Tag','pl_spm_outliers');
    group=ones(n,1);
    if ~isempty(out.outliers)
        group(out.outliers)=2;
    end
    spmplot(Y,group,plo);
    set(gcf,'Name',' S estimator: scatter plot matrix with outliers highlighted');

end

if options.ysave
    out.Y = Y;
end

end

% -------------------------------------------------------------------
% subfunction IRWLSmult
% -------------------------------------------------------------------

function outIRWLS=IRWLSmult(Y,initialloc, initialshape, refsteps, reftol, c, kc,initialscale)
%IRWLSmult (iterative reweighted least squares) does refsteps refining steps from initialloc
% for refsteps times or till convergence.
%
%  Required input arguments:
%
%    Y: Data matrix containining n observations on v variables.
%       Rows of Y represent observations, and columns represent variables.
% initialloc   : v x 1 vector containing initial estimate of location
% initialshape: v x v initial estimate of shape matrix
%   refsteps  : scalar, number of refining (IRLS) steps
%   reftol    : relative convergence tolerance for the fully iterated
%               best candidates. Deafult value is 1e-7
%    c        : scalar, tuning constant of the equation for Tukey biweight
%   kc        : scalar, tuning constant linked to Tukey's biweight
%
%  Optional input arguments:
%
% initialscale: scalar, initial estimate of the scale. If not defined,
%               scaled MAD of Mahalanobis distances from initialloc and
%               initialshape is used.
%
%
%  Output:
%
%  The output consists of a structure 'outIRWLS' containing:
%      outIRWLS.loc    : v x 1 vector. Estimate of location after refsteps
%                        refining steps.
%      outIRWLS.shape  : v x v matrix. Estimate of the shape matrix after
%                        refsteps refining steps.
%      outIRWLS.scale  : scalar. Estimate of scale after refsteps refining
%                        step.
%      outIRWLS.weights: n x 1 vector. Weights assigned to each observation
%
% In the IRWLS procedure the value of loc and the value of the scale and
% of the shape matrix are updated in each step

v = size(Y,2);
loc = initialloc;
% Mahalanobis distances from initialloc and Initialshape
mahaldist = sqrt(mahalFS(Y, initialloc, initialshape));


% The scaled MAD of Mahalanobis distances is the default for the initial scale
if (nargin < 8)
    initialscale = median(abs(mahaldist))/.6745;
end

scale = initialscale;

iter = 0;
locdiff = 9999;

while ( (locdiff > reftol) && (iter < refsteps) )
    iter = iter + 1;
    
    % Solve for the scale
    scale = scale* sqrt( mean(TBrho(mahaldist/scale,c))/kc);
    % mahaldist = vector of Mahalanobis distances from robust centroid and
    % robust shape matrix, which is changed in each step
    
    % compute w = n x 1 vector containing the weights (using TB)
    weights = TBwei(mahaldist/scale,c);
    
    % newloc = new estimate of location using the weights previously found
    % newloc = \sum_{i=1}^n y_i w(d_i) / \sum_{i=1}^n w(d_i)
    newloc = sum(bsxfun(@times,Y,weights),1)/sum(weights);
    
    % exit from the loop if the new loc has singular values. In such a
    % case, any intermediate estimate is not reliable and we can just
    % keep the initial loc and initial scale.
    if (any(isnan(newloc)))
        newloc = initialloc;
        newshape = initialshape;
        scale = initialscale;
        weights=NaN;
        break
    end
    
    % Res = n x v matrix which contains deviations from the robust estimate
    % of location
    Res = bsxfun(@minus,Y, newloc);
    
    % Multiplication of newshape by a constant (e.g. v) is unnecessary
    % because final value of newshape remains the same as det(newshape)=1.
    % For the same reason newshape remains the same if we use weights or
    % weights*(c^2/6)
    newshape= (Res')*bsxfun(@times,Res,weights);
    newshape = det(newshape)^(-1/v)*newshape;
    
    % Compute MD
    mahaldist = sqrt(mahalFS(Y,newloc,newshape));
    
    % locdiff is linked to the tolerance
    locdiff = norm(newloc-loc,1)/norm(loc,1);
    loc = newloc;
    
end

outIRWLS.loc = newloc;
outIRWLS.shape = newshape;
outIRWLS.scale = scale;
outIRWLS.weights=weights;
end
%FScategory:MULT-Multivariate
