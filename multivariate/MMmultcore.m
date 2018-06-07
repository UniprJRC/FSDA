function out = MMmultcore(Y,loc0,shape0,auxscale,varargin)
%MMmultcore computes multivariate MM estimators for a selected fixed scale
%
%
%<a href="matlab: docsearchFS('MMmultcore')">Link to the help function</a>
%
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
%        loc0 :  initial estimate of location. Vector.
%               Vector containing initial estimate of location (generally
%               an S estimate with high breakdown point, eg 0.5)
%     shape0 :  initial estimate of shape. Matrix. 
%               v x v matrix containing initial estimate of shape
%               (generally an S estimate with high breakdown point, eg 0.5)
%   auxscale :  initial estimate of scale. Scalar.
%               Scalar containing estimate of the scale (generally an S
%               estimate with high breakdown point).
%
%  Optional input arguments:
%
%      eff     : nominal efficiency. Scalar.
%                Scalar defining nominal efficiency (i.e. a number between
%                 0.5 and 0.99). The default value is 0.95
%                 Asymptotic nominal efficiency is:
%                 $(\int \psi' d\Phi)^2 / (\psi^2 d\Phi)$
%                 For example if eff=0.95 and v=3, c=5.49 in case of location
%                 efficiency or c=6.096 in case of shape efficiency
%                 Example - 'eff',0.99
%                 Data Types - double
%     effshape : locacation or scale effiicency. dummy scalar. 
%                If effshape=1 efficiency refers to shape 
%                efficiency else (default) efficiency refers to location
%                 Example - 'effshape',1
%                 Data Types - double
%     refsteps  : Maximum iterations. Scalar.
%                 Scalar defining maximum number of iterations in the MM
%                 loop. Default value is 100.
%                 Example - 'refsteps',10
%                 Data Types - double
%      reftol: Tolerance. Scalar.
%                 Scalar controlling tolerance in the MM loop.
%                 Default value is 1e-7
%                 Example - 'tol',1e-10
%                 Data Types - double
%     conflev :  Confidence level which is
%               used to declare units as outliers. Scalar.
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%                 Example - 'conflev',0.99
%                 Data Types - double
%       plots : Plot on the screen. Scalar or structure.
%               If plots = 1, generates a plot of Mahalanobis distances
%               against index number. The confidence level used to draw the
%               confidence bands for the MD is given by the input
%               option conflev. If conflev is not specified a nominal 0.975
%               confidence interval will be used.
%                 Example - 'plots',0 
%                 Data Types - single | double
%       nocheck : Check input arguments. Scalar. If nocheck is equal to 1
%                 no check is performed on
%                 matrix Y. As default nocheck=0.
%               Example - 'nocheck',1 
%               Data Types - double
%       ysave : input data matrix Y is saved into the output
%                structure out. Scalar. 
%               Scalar that is set to 1 to request that the data matrix Y
%               is saved into the output structure out. This feature is
%               meant at simplifying the use of function malindexplot.
%               Default is 0, i.e. no saving is done.
%               Example - 'ysave',1 
%
%  Output:
%
%  out :     A structure containing the following fields
%     out.loc =  v x 1 vector. Estimate of location after refsteps
%                  refining steps
%    out.shape=  v x v matrix. Estimate of the shape matrix after refsteps
%                  refining steps
%      out.cov=  v x v matrix. Estimate of the covariance matrix after
%                  refsteps refining steps
%       out.md=  n x 1 vector containing the estimates of the robust
%                  Mahalanobis distances (in squared units)
% out.outliers=  A vector containing the list of the units declared as
%                  outliers using confidence level specified in input
%                  scalar conflev
%  out.conflev= Confidence level that was used to declare outliers
%  out.weights=  n x 1 vector. Weights assigned to each observation
%    out.class= 'MM'
%       out.Y : Data matrix Y. The field is present if option ysave
%               is set to 1.
%
%
% More About:
%
% This routine does iterative reweighted least squares (IRWLS) steps from
% "initial location" (loc0) and shape matrix (shape0) keeping the estimate
% of the scale (auxscale) fixed.
%
% See also: Smult, MMmult
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
% and the R library robustbase http://robustbase.r-forge.r-project.org/.
% The core of these routines, e.g. the resampling approach, however, has
% been completely redesigned, with considerable increase of the
% computational performance.
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('MMmultcore')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % MMmultcore with all default options.
    % Determine, e.g., an S estimate and extract the required arguments for
    %  the MM estimate.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out]=Smult(Ycont);
    outMM=MMmultcore(Ycont,out.loc,out.shape,out.scale)
%}

%{
    %% MMmultcore with optional arguments.
    % Determine, e.g., an S estimate and extract the required arguments for
    % the MM estimate.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out]=Smult(Ycont,'plots',1);
    outMM=MMmultcore(Ycont,out.loc,out.shape,out.scale,'plots',1,'nocheck',1)
%}

%% Beginning of code

nnargin=nargin;
vvarargin=varargin;
[Y,n,v]= chkinputM(Y,nnargin,vvarargin);

% default nominal efficiency
effdef = 0.95;
% by default the nominal efficiency refers to location efficiency
effshapedef = 0;
% default value of number of maximum refining iterations
refstepsdef=50;
% default value of tolerance for the refining steps convergence
reftoldef=1e-6;

% store default values in the structure options
options=struct('refsteps',refstepsdef,'reftol',reftoldef,...
    'eff',effdef,'effshape',effshapedef,'conflev',0.975,...
    'plots',0,'nocheck',0,'ysave',0);

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:MMmultcore:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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

eff     = options.eff;      % nominal efficiency
effshape= options.effshape; % nominal efficiency refers to shape or location
refsteps= options.refsteps; % maximum refining iterations
reftol  = options.reftol;   % tolerance for refining iterations covergence

% constant for the Tukey's biweight function
if effshape==1
    c = TBeff(eff,v,1);
else
    c = TBeff(eff,v);
end

% Ytilde = deviations from centroid
Ytilde = bsxfun(@minus,Y, loc0);

% mahaldist = vector of Mahalanobis distances using as shape matrix shape0
mahaldist = sqrt(sum((Ytilde/shape0).*Ytilde,2));

% newobj = objective function which has to be minimized (keeping the value
% of the scale fixed)
origobj  = mean(TBrho(mahaldist/auxscale,c));
newobj   = origobj;

% compute M-estimate with auxiliary scale through IRLS steps, starting
% from S-estimate
iter   = 0;
oldobj = newobj + 1;
while ((oldobj - newobj) > reftol) && (iter < refsteps)
    iter = iter +1;
    
    % w = psi(x)/x
    weights = TBwei(mahaldist/auxscale,c);
    
    % Find new estimate of location using the weights previously found
    % newloc = \sum_{i=1}^n y_i weights(d_i) / \sum_{i=1}^n weights(d_i)
    newloc = sum(bsxfun(@times,Y,weights),1)/sum(weights);
    % exit from the loop if the new loc has singular values. In such a
    % case, any intermediate estimate is not reliable and we can just
    % keep the initial loc and initial scale.
    if (any(isnan(newloc)))
        newloc = loc0;
        newshape = shape0;
        weights=NaN;
        break
    end
    
    % Find new estimate of scaled covariance matrix using  the weights
    % previously found
    newshape = bsxfun(@times,Ytilde,weights)'*Ytilde;
    
    % newshape is a var cov matrix with determinant equal to 1
    newshape = det(newshape)^(-1/v)*newshape;
    
    % Compute Mahalanobis distances from centroid newloc and var newshape
    mahaldist=sqrt(mahalFS(Y,newloc,newshape));
    
    oldobj = newobj;
    newobj = mean(TBrho(mahaldist/auxscale,c));
end

out.class   = 'MM';
if newobj <= origobj
    out.loc = newloc;
    out.shape = newshape;
    out.cov = auxscale^2*newshape;
    out.weights = weights;
else % isn't supposed to happen
    warning(warnrank.state,'Initial solutions for location and shape parameters have been kept')
    warning(warnrank.state,'Because MM-loop does not produce better estimates');
    out.loc = loc0;
    out.shape = shape0;
    out.cov = auxscale^2*shape0;
    out.weights = NaN;
end
out.md = mahalFS(Y,out.loc,out.cov);

% Store in output structure the outliers found with confidence level conflev
conflev = options.conflev;
seq = 1:n;
out.outliers = seq(out.md > chi2inv(conflev,v) );
out.conflev = conflev;

if options.ysave
    out.Y = Y;
end

plo=options.plots;

% Plot Mahalanobis distances with outliers highlighted
if isstruct(plo) || (~isstruct(plo) && plo~=0)
    
    [n,v]=size(Y);
    
    laby='MM Mahalanobis distances';
    malindexplot(out.md,v,'conflev',conflev,'laby',laby,'numlab',out.outliers);
    
    figure('Tag','pl_spm_outliers');
    group=ones(n,1);
    if ~isempty(out.outliers)
        group(out.outliers)=2;
    end
    spmplot(Y,group,plo);
    set(gcf,'Name',' MM estimator: scatter plot matrix with outliers highlighted');   
    
end


end

%FScategory:MULT-Multivariate