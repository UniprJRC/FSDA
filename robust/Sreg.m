function [out , varargout] = Sreg(y,X,varargin)
%Sreg computes S estimators in linear regression
%
%<a href="matlab: docsearch('Sreg')">Link to the help function</a>
%
%  Required input arguments:
%
%    y:      A vector with n elements that contains the response variable.
%            It can be both a row or a column vector.
%    X :     Data matrix of explanatory variables (also called 'regressors')
%            of dimension (n x p-1). Rows of X represent observations, and
%            columns represent variables. 
%
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%
%  Optional input arguments:
%
%   intercept : If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%      bdp    : scalar defining breakdown point (i.e a number between 0 and 0.5)
%               The default value is 0.5
%      nsamp  : scalar defining number of subsamples of size p which have
%               to be extracted (if not given, default = 1000)
%    refsteps : scalar defining number of refining iterations in each
%               subsample (default = 3).
%               refsteps = 0 means "raw-subsampling" without iterations.
%     reftol  : scalar. Default value of tolerance for the refining steps
%               The default value is 1e-6;
%refstepsbestr: scalar defining number of refining iterations for each
%               best subset (default = 50).
% reftolbestr : scalar. Default value of tolerance for the refining steps
%               for each of the best subsets
%               The default value is 1e-8;
%     minsctol: scalar. Default value of tolerance for the iterative
%               procedure for finding the minimum value of the scale
%               for each subset and each of the best subsets
%               (It is used by subroutine minscale.m)
%               The default value is 1e-7;
%      bestr  : scalar defining number of "best betas" to remember from the
%               subsamples. These will be later iterated until convergence
%               (default=5)
%     conflev : Scalar between 0 and 1 containing confidence level which is
%               used to declare units as outliers.
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%        msg  : scalar which controls whether to display or not messages
%               on the screen If msg==1 (default) messages are displyed
%               on the screen about estimated time to compute the estimator
%               and the warnings about
%               'MATLAB:rankDeficientMatrix', 'MATLAB:singularMatrix' and
%               'MATLAB:nearlySingularMatrix' are set to off
%               else no message is displayed on the screen
%      nocheck: Scalar. If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0.
%       plots : Scalar or structure.
%               If plots = 1, generates a plot with the robust residuals
%               against index number. The confidence level used to draw the
%               confidence bands for the residuals is given by the input
%               option conflev. If conflev is not specified a nominal 0.975
%               confidence interval will be used.
%       yxsave : scalar that is set to 1 to request that the response 
%                vector y and data matrix X are saved into the output
%                structure out. Default is 0, i.e. no saving is done.
%  Output:
%
%  The output consists of a structure 'out' containing the following fields:
%
%            out.beta : vector containing the S estimator of regression
%                       coefficients
%            out.scale: scalar containing the estimate of the scale
%                       (sigma). This is the value of the objective function
%              out.bs : p x 1 vector containing the units forming best subset
%                       associated with S estimate of regression coefficient.
%        out.residuals: n x 1 vector containing the estimates of the robust
%                       scaled residuals
%        out.outliers : this output is present only if conflev has been
%                       specified. It is a vector containing the list of
%                       the units declared as outliers using confidence
%                       level specified in input scalar conflev
%         out.conflev : confidence level which is used to declare outliers.
%                       Remark: scalar out.conflev will be used to draw the
%                       horizontal line (confidence band) in the plot.
%         out.singsub : Number of subsets wihtout full rank. Notice that
%                       out.singsub > 0.1*(number of subsamples) produces a
%                       warning
%         out.weights : n x 1 vector containing the estimates of the weights
%           out.class : 'S'
%            out.y    : response vector Y. The field is present if option 
%                       yxsave is set to 1.
%            out.X    : data matrix X. The field is present if option 
%                       yxsave is set to 1.
%
%  Optional Output:
%
%            C     : matrix of the indices of the samples extracted for
%                    computing the estimate
%
% References:
%
% ``Robust Statistics, Theory and Methods'' by Maronna, Martin and Yohai;
% Wiley 2006.
%
% Acknowledgements
%
% This function follows the lines of MATLAB/R code developed during the
% years by many authors.
% For more details see http://www.econ.kuleuven.be/public/NDBAE06/programs/
% and the R library robustbase http://robustbase.r-forge.r-project.org/
% The core of these routines, e.g. the resampling approach, however, has
% been completely redesigned, with considerable increase of the
% computational performance.
%
% Copyright 2008-2013.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti
%
%
%<a href="matlab: docsearch('sreg')">Link to the help page for this function</a>
% Last modified 02-May-2013
%
% Examples:

%{
    % Sreg with all default options
    % Run this code to see the output shown in the help file
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    [out]=Sreg(ycont,X);
%}


%% Beginning of code
nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% default value of break down point
bdpdef=0.5;


% default values of subsamples to extract
ncomb=bc(n,p);
nsampdef=min(1000,ncomb);

% default value of number of refining iterations (C steps) for each extracted subset
refstepsdef=3;
% default value of tolerance for the refining steps convergence for  each extracted subset
reftoldef=1e-6;
% default value of number of best betas to remember
bestrdef=5;
% default value of number of refining iterations (C steps) for best subsets
refstepsbestrdef=50;
% default value of tolerance for the refining steps convergence for best subsets
reftolbestrdef=1e-8;
% default value of tolerance for finding the minimum value of the scale
% both for each extracted subset and each of the best subsets
minsctoldef=1e-7;

% store default values in the structure options
options=struct('intercept',1,'nsamp',nsampdef,'refsteps',refstepsdef,...
    'reftol',reftoldef,'refstepsbestr',refstepsbestrdef,'reftolbestr',reftolbestrdef,...
    'minsctol',minsctoldef,'bestr',bestrdef,'bdp',bdpdef,...
    'plots',0,'conflev',0.975,'nocheck',0,'msg',1,'yxsave',0);

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('Error:: number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end



% Get user values of warnings
warnrank=warning('query','MATLAB:rankDeficientMatrix');
warnsing=warning('query','MATLAB:singularMatrix');
warnnear=warning('query','MATLAB:nearlySingularMatrix');
% Set them to off inside this function
% At the end of the file they will be restored to previous values
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

bdp = options.bdp;              % break down point
refsteps = options.refsteps;    % refining steps
reftol = options.reftol;        % tolerance for refining steps
bestr = options.bestr;          % best betas for refining steps till convergence
nsamp = options.nsamp;          % subsamples to extract
minsctol = options.minsctol;    % tolerance for finding minimum value of the scale for each subset
refstepsbestr=options.refstepsbestr;  % refining steps for the best subsets
reftolbestr=options.reftolbestr;      % tolerance for refining steps for the best subsets
msg=options.msg;                % Scalar which controls the messages displayed on the screen

% Find constant c linked to Tukey's biweight
% rho biweight is strictly increasing on [0 c] and constant on [c \infty)
% E(\rho) = kc = (c^2/6)*bdp, being kc the K of Rousseeuw and Leroy
c  = TBbdp(bdp,1);
kc = bdp*(c^2/6);

bestbetas = zeros(bestr,p);
bestsubset = bestbetas;
bestscales = Inf * ones(bestr,1);
sworst = Inf;

% singsub = scalar which will contain the number of singular subsets which
% are extracted (that is the subsets of size p which are not full rank)
singsub=0;

% ij is a scalar used to ensure that the best first bestr non singular
% subsets are stored
ij=1;

%% Extract in the rows of matrix C the indexes of all required subsets
[C,nselected] = subsets(nsamp,n,p,ncomb,msg);
% Store the indices in varargout
if nargout==2
    varargout={C};
end

% initialise and start timer.
tsampling = ceil(min(nselected/100 , 1000));
time=zeros(tsampling,1);


for i = 1:nselected
    
    if i <= tsampling, tic; end
    
    % extract a subset of size p
    index = C(i,:);
    
    Xb = X(index,:);
    yb = y(index);
    
    % beta estimate
    beta = Xb\yb;
    
    if ~isnan(beta(1)) && ~isinf(beta(1))
        
        if refsteps > 0
            % do the refsteps refining steps
            % kc determines the breakdown point
            % c is linked to the biweight function
            tmp = IRWLSreg(y,X,beta,refsteps,reftol,c,kc);
            
            betarw = tmp.betarw;
            scalerw = tmp.scalerw;
            resrw = y - X * betarw;
        else
            % no refining steps
            betarw = beta;
            resrw = y - X * betarw;
            scalerw = median(abs(resrw))/.6745;
        end
        
        % to find s, save first the best bestr scales (deriving from non
        % singular subsets) and, from iteration bestr+1 (associated to
        % another non singular subset), replace the worst scale
        % with a better one as follows
        if ij > bestr
            % compute the objective function using current residuals and
            % the worst estimate of the scale among the bests previously
            % stored
            % scaletest = (1/n) \sum_i=1^n (u_i/(sworst*c))
            scaletest = mean(TBrho(resrw/sworst,c));
            if scaletest < kc
                % Find position of the maximum value of previously stored
                % best scales
                [~,ind] = max(bestscales);
                sbest = minscale(resrw,c,kc,scalerw,minsctol);
                % Store sbest, betarw and indexes of the units forming the
                % best subset associated with minimum value
                % of the scale
                bestscales(ind) = sbest;
                bestbetas(ind,:) = betarw';
                % best subset
                bestsubset(ind,:)=index;
                % sworst = the best scale among the bestr found up to now
                sworst = max(bestscales);
            end
        else
            bestscales(ij) = minscale(resrw,c,kc,scalerw,minsctol);
            bestbetas(ij,:) = betarw';
            bestsubset(ij,:) = index;
            ij=ij+1;
        end
    else
        singsub=singsub+1;
    end
    
    % Write total estimate time to compute final estimate
    if i <= tsampling
        % sampling time until step tsampling
        time(i)=toc;
    elseif i==tsampling+1
        % stop sampling and print the estimated time
        if msg==1
            fprintf('Total estimated time to complete S estimate: %5.2f seconds \n', nselected*median(time));
        end
    end
    
end

% perform C-steps on best 'bestr' solutions, till convergence or for a
% maximum of refstepsbestr steps using a convergence tolerance as specified
% by scalar reftolbestr

superbestscale = Inf;

for i=1:bestr
    tmp = IRWLSreg(y,X,bestbetas(i,:)',refstepsbestr,reftolbestr,c,kc,bestscales(i));
    
    if tmp.scalerw < superbestscale
        superbestscale = tmp.scalerw;
        superbestbeta = tmp.betarw;
        superbestsubset = bestsubset(i,:);
        weights = tmp.weights;
    end
end

% Store in output structure \beta, s, best subset and vector of S-weights
out.beta = superbestbeta;
out.scale = superbestscale;
out.bs = superbestsubset;
out.weights = weights;

% compute and store in output structure the S robust scaled residuals
out.residuals=(y-X*out.beta)/out.scale;

% Store in output structure the number of singular subsets
out.singsub=singsub;
if singsub/nselected>0.1;
    disp('------------------------------')
    disp(['Warning: Number of subsets without full rank equal to ' num2str(100*singsub/nselected) '%'])
end

% Restore the previous state of the warnings
warning(warnrank.state,'MATLAB:rankDeficientMatrix');
warning(warnsing.state,'MATLAB:singularMatrix');
warning(warnnear.state,'MATLAB:nearlySingularMatrix');


% Store in output structure the outliers found with confidence level conflev
conflev = options.conflev;
out.conflev = conflev;

conflev = (conflev+1)/2;
seq = 1:n;
out.outliers = seq( abs(out.residuals)>norminv(conflev) );

if options.yxsave
    if options.intercept==1;
        % Store X (without the column of ones if there is an intercept)
        out.X=X(:,2:end);
    else
        out.X=X;
    end
    % Store response
    out.y=y;
end

% Plot resindexplot with outliers highlighted
if options.plots==1;
    laby='Scaled S residuals';
    resindexplot(out.residuals,'conflev',out.conflev,'laby',laby,'numlab',out.outliers);
end

end

% -------------------------------------------------------------------
% subfunction IRWLSreg
% -------------------------------------------------------------------

function outIRWLS = IRWLSreg(y,X,initialbeta,refsteps,reftol,c,kc,initialscale)
%IRWLSreg (iterative reweighted least squares) does refsteps refining steps from initialbeta
%
%  Required input arguments:
%
%    y:         A vector with n elements that contains the response variable.
%               It can be both a row or column vector.
%    X :        Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p). Rows of X represent observations, and
%               columns represent variables.
% initialbeta : p x 1 vector containing initial estimate of beta
%   refsteps  : scalar, number of refining (IRLS) steps
%   reftol    : relative convergence tolerance
%               Default value is 1e-7
%    c        : scalar, tuning constant of the equation for Tukey biweight
%   kc        : scalar, tuning constant linked to Tukey's biweight
%
%  Optional input arguments:
%
% initialscale: scalar, initial estimate of the scale. If not defined,
%               scaled MAD of residuals is used.
%
%  Output:
%
%  The output consists of a structure 'outIRWLS' containing the following fields:
%      betarw  : p x 1 vector. Estimate of beta after refsteps refining steps
%     scalerw  : scalar. Estimate of scale after refsteps refining step
%     weights  : n x 1 vector. Weights assigned to each observation
%
% In the IRWLS procedure the value of beta and the value of the scale are
% updated in each step

% Residuals for the initialbeta
res = y - X * initialbeta;

% The scaled MAD of residuals is the initial scale estimate default value
if (nargin < 8)
    initialscale = median(abs(res))/.6745;
end

beta = initialbeta;
scale = initialscale;

iter = 0;
betadiff = 9999;

while ( (betadiff > reftol) && (iter < refsteps) )
    iter = iter + 1;
    
    % Solve for the scale
    scale = scale * sqrt(mean( TBrho(res/scale,c) ) / kc );
    
    % Compute n x 1 vector of weights (using TB)
    weights = TBwei(res/scale,c);
    sqweights = weights.^(1/2);
    
    % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
    Xw = bsxfun(@times, X, sqweights);
    yw = y .* sqweights;
    
    % estimate of beta from (re)weighted regression (RWLS)
    newbeta = Xw\yw;
    
    % exit from the loop if the new beta has singular values. In such a
    % case, any intermediate estimate is not reliable and we can just
    % keep the initialbeta and initial scale.
    if (any(isnan(newbeta)))
        newbeta = initialbeta;
        scale = initialscale;
        weights = NaN;
        break
    end
    
    % betadiff is linked to the tolerance (specified in scalar reftol)
    betadiff = norm(beta - newbeta,1) / norm(beta,1);
    
    % update residuals and beta
    res = y - X * newbeta;
    beta = newbeta;
    
end

% store final estimate of beta
outIRWLS.betarw = newbeta;
% store final estimate of scale
outIRWLS.scalerw = scale;
% store final estimate of the weights for each observation
outIRWLS.weights=weights;
end
