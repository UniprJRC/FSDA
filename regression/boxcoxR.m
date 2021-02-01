function out=boxcoxR(y,X, varargin)
%boxcoxR finds MLE of lambda in linear regression (and confidence interval) using Box Cox, YJ or extended YJ  transformation
%
%<a href="matlab: docsearchFS('boxcoxR')">Link to the help function</a>
%
% The function computes the profile log Likelihood for a range of values of
% the transformation parameter (lambda) and computes the MLE of lambda in
% the supplied range. Supported families are Box Cox, Yeo and Johnson and
% extended Yeo and Johnson (Atkinson et al. 2020).
%
%               The profile log-likelihood is computed as:
%               \[
%               -(n/2) \log( (y(\lambda)-X\beta(\lambda))'(y(\lambda)-X\beta(\lambda))/n) +\log J
%               \]
%               where  $y(\lambda)$ is the vector of transformed
%               observations using Box Cox family,  Yeo and Johnson or
%               extended Yao and Johnson family
%               \[
%               \beta(\lambda) = (X'X)^{-1} X' y(\lambda)
%               \]
%               and $J$ is the Jacobian of the transformation.
%
%  Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each element of y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, but observations (rows) with these values will
%               automatically be excluded from the computations.
%  X :          Predictor variables. Matrix. Matrix of explanatory
%               variables (also called 'regressors') of dimension n x (p-1)
%               where p denotes the number of explanatory variables
%               including the intercept.
%               Rows of X represent observations, and columns represent
%               variables. By default, there is a constant term in the
%               model, unless it is explicitly removed using input option
%               intercept; so, do not include a column of 1s in X. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with such values will
%               automatically be excluded from the computations.
%
%
% Optional input arguments:
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as comma-separated pair consisting of
%                 'Intercept' and either true or false, to respectively
%                 include or remove the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%    family :   parametric transformation to use. String. String which
%               identifies the family of transformations which must be
%               used. Character. Possible values are 'BoxCox' (default),
%               'YJ' (Yao and Yohnson) and 'YJpn' (extended Yeo and
%               Johnson).
%               The Box-Cox family of power transformations equals
%               $(y^{\lambda}-1)/\lambda$ for $\lambda$ not equal to zero,
%               and $\log(y)$ if $\lambda = 0$.
%               The Yeo-Johnson (YJ) transformation is the Box-Cox
%               transformation of $y+1$ for nonnegative values, and of
%               $|y|+1$ with parameter 2-lambda for y negative.
%               The extended Yeo-Johnson (YJpn) transformation is like
%               Yeo-Johnson, but admits two values of the transformation
%               parameters, respectively for positive and negative
%               observations.
%               Remark. BoxCox family can be used only if input y is
%               positive. Yeo-Johnson (and extended Yeo-Johnson family of
%               transformations do not have this limitation).
%               Example - 'family','YJ'
%               Data Types - char
%
%       nocheck : Check input arguments. Scalar. If nocheck is equal to 1
%                 no check is performed on vector y and matrix X. This
%                 means that y and X are left unchanged. Note also that the
%                 additional column of ones for the intercept is not added.
%                 As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
%        conflev : Confidence level for lambda. Scalar. The scalar is 
%                  between 0 and 1 and determines the confidence level for
%                  lambda, based on the asymptotic $chi1^2$ of twice the
%                  loglikelihood ratio. The default conflev value is 0.95;
%               Example - 'conflev',0.99
%               Data Types - double
%
%          laseq : Sequence of values of lambda to consider. Vector. Vector
%                  which contains the sequence of values of lambda for
%                  which the profile loglikelihood has to be computed if
%                  family is 'BoxCox' or 'YJ'. The default value of laseq
%                  is -2:0.001:2. This optional input is ignored if family
%                  is 'YJpn';
%               Example - 'laseq',[-1:0.001;0.7]
%               Data Types - double
%
%          laseqPos : Transformation for positive observations.
%                  Vector. Vector which contains the sequence of values of
%                  lambda which are used to transform positive observations
%                  when family 'YJpn'. The default value of laseqPos is
%                  -2:0.01:2. This optional input parameter is ignored if
%                  family is 'BoxCox' or 'YJ';
%               Example - 'laseqPos',[-1:0.001;0.7]
%               Data Types - double
%
%          laseqNeg : Transformation for negative observations.
%                  Vector. Vector which contains the sequence of values of
%                  lambda which are used to transform negative observations
%                  when family 'YJpn'. The default value of laseqNeg is
%                  -2:0.01:2. This optional input is ignored if family is
%                  'BoxCox' or 'YJ';
%               Example - 'laseqNeg',[-1:0.001;0.7]
%               Data Types - double
%
%   plots  :    Profile log likelihood for lambda. Boolean.
%               It specifies whether to show the profile log likelihood of
%               lambda. If plots is true, the plot of the profile
%               loglikelihood is produced together with the requested
%               confidence interval. The default value of prolik is false,
%               that is no plot is produced. If family is 'YJpn', a contour
%               plot is produced.
%               Example - 'plots',true
%               Data Types - boolean
%
%    usefmin :  use solver to find MLE of lambda. Boolean or struct.
%               This option applies only if family is YJpn. If usefmin is
%               true or usefmin is a struct, the maximum likelihood
%               estimates of $\lambda_P$ and $\lambda_N$ is computed using
%               the MATLAB solvers fminsearch or fminunc. The default value
%               of usefmin is false, that is the likelihood is evaluated at
%               the points laseqPos and laseqNeg without the solver.
%               If usefmin is a structure it may contain the following
%               fields:
%               usefmin.MaxIter = Maximum number of iterations (default is 1000).
%               usefmin.TolX   = Termination tolerance for the parameters
%                   (default is 1e-7).
%               usefmin.solver = name of the solver. Possible values are
%                   'fminsearch' (default) and 'fminunc'. fminunc needs the
%                   optimization toolbox.
%               usefmin.displayLevel = amount of information displayed by
%                   the algorithm. possible values are 'off' (displays no
%                   information, this is the default), 'final' (displays
%                   just the final output) and 'iter' (displays iterative
%                   output to the command window).
%               Example - 'usefmin',true
%               Data Types - boolean or struct
%
%
% Output:
%
%         out:   structure which contains the following fields
%
% out.lahat  =  best estimate of lambda. Scalar or vector of length 2.
%               out.lahat is a scalar if family is 'BoxCox' or 'YJ',
%               otherwise if family is 'YJpn' it is a vector of length 2
%               containing respectively MLE of transformation parameter for
%               positive and negative observations. This is the best
%               estimate among the values of lambda which are supplied in
%               input vector laseq if family is 'BoxCox' or 'YJ', or in
%               input vectors laseqPos and laseqNeg, if family is 'YJpn'.
% out.lahatci  = confidence intervals for MLE of lambda computed using Chi2
%               approximation using confidence level specified in input
%               option conflev. This argument is present only if family is
%               'BoxCox' or 'YJ'.
% out.LogLik   = matrix containing the value of the profile log-likelihood
%               for each value in laseq or laseqPos and laseqNeg. If family
%               is BoxCox' or 'YJ', the dimension of out.LogLik is
%               length(laseq)-by-2. In this case the first column contains
%               the values of laseq and the second column the values of the
%               profile log lik. If family is 'YJpn', the dimension of
%               out.LogLik is:
%               * length(laseqPos)-by-length(laseqNeg) if option usefmin is
%                 false (default) and therefore the maximization routine is
%                 not called.
%               * 9-by-9 matrix containing the value of the likelihood in
%                 correspondence of the points -2:0.5:2, if input option
%                 usefmin is true and boxcoxR is called for the first time.
%                 In order to find MLE of laPos and laNeg, the likelihood
%                 is computed for the 81 values in the meshgrid -2:0.5:2 to
%                 find a preliminary estimate of laPos and laNeg for the
%                 solver.
%               * empty value if input option usefmin is true and the solver
%                 has used the cached version of lambda to initialize the
%                 optimization routine. If you prefer that the maximization
%                 routine does not use the cached version of lambda,
%                 execute instruction clear boxcoxR before calling it.
%
% out.exitflag =  flag which informs about convergence. exitflag = 0
%                 implies normal convergence, else no convergence has been
%                 obtained. This ouptut is present only if input option
%                 usefmin is true or struct and family is YJpn.
%
%
% See also Score, FSRfan
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York. [see pp. 83-84]
% Box, G.E.P. and Cox, D.R. (1964), The analysis of transformations,
% Journal of the Royal Statistical Society, Vol. 26, pp. 211-252.
% Yeo, I.K and Johnson, R. (2000), A new family of power transformations to
% improve normality or symmetry, "Biometrika", Vol. 87, pp. 954-959.
% Atkinson, A.C., Riani, M. and Corbellini C. (2019), The Analysis of
% Transformations for Profit and Loss Data, "Journal of the Royal
% Statistical Society. Series C: Applied Statistics",
% https://doi.org/10.1111/rssc.12389 [ARC]
% Atkinson, A.C. Riani, M. and Corbellini A. (2020), The Box-Cox
% Transformation: Review and Extensions, "Statistical Science", in press.
%
% Acknowledgements:
%
% This function has been inspired by sumbmission
% https://www.mathworks.com/matlabcentral/fileexchange/10419-box-cox-power-transformation-for-linear-models
% in the file exchange written by Hovav Dror, hovav@hotmail.com, March 2006
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('boxcoxR')">Link to the help function</a>
%
%$LastChangedDate:: 2019-05-14 16:04:25 #$: Date of the last commit

% Examples:

%{
    %% boxcoxR with all default options.
    % Use the wool data.
    load('wool.txt','wool');
    y=wool(:,4);
    X=wool(:,1:3);
    out=boxcoxR(y,X);
    disp(['Estimate of lambda using Box Cox is =',num2str(out.lahat)])
%}

%{
    %% boxcoxR using YJ transformation.
    load('wool.txt','wool');
    y=wool(:,4);
    X=wool(:,1:3);
    out=boxcoxR(y,X,'family','YJ');
    disp(['Estimate of lambda using YJ family is =',num2str(out.lahat)])
%}

%{
    % Example of use of option confint.
    load('wool.txt','wool');
    y=wool(:,4);
    X=wool(:,1:3);
    % A 99.9 per cent confidence interval for lambda is used.
    out=boxcoxR(y,X,'conflev',0.999);
%}

%{
    % Example of use of option plots.
    load('wool.txt','wool');
    y=wool(:,4);
    X=wool(:,1:3);
    % Plot the profile loglikelihood
    out=boxcoxR(y,X,'plots',1);
%}

%{
    %% Example of use of option plots combined with  laseq.
    load('wool.txt','wool');
    y=wool(:,4);
    X=wool(:,1:3);
    % Plot the profile loglikelihood in the interval [-1 1]
    laseq=[-1:0.0001:1];
    out=boxcoxR(y,X,'plots',1,'laseq',laseq);
%}


%{
    %% Example of use of option family.
    close all
    YY=load('fondi_large.txt');
    y=YY(:,2);
    X=YY(:,[1 3]);
    out=boxcoxR(y,X,'family','YJpn','plots',1);
    % The contour plot suggestes that while positive observatios do not
    % have to transformed, negative observations have to be transformed using
    % lambda=0. For more details see Atkinson Riani and Corbellini (2020)
%}

%{
   %% Ex of the use of option usefmin.
    rng(500)
    % Generate regression data
    outsim=simulateLM(100,'R2',0.95);
    yori = outsim.y;
	X = outsim.X;
    % Transform in a different way positive and negative values
    laPos=0.2;
    laNeg=0.8;
    y=normYJpn(yori,[],[laPos laNeg],'inverse',true);
    % Use solver to find MLE of laPos and laNeg
    usefmin=struct;
    % specify maximum number of iterations
    usefmin.MaxIter=100;
    % Function boxcoxR initializes the optimization routine with the value
    % of lambda from the last call of the optmization routine. This trick
    % is very useful during the forward search when we use in step m+1 as
    % initial guess of laP and laN the final estimate of laP and laN in
    % step m.
    % The instruction below clear persisten variables in function boxcoxR
    clear boxcoxR
    out=boxcoxR(y,X,'family','YJpn','plots',1,'usefmin',usefmin);
    disp('MLE of laPos and laNeg')
    disp(out.lahat)
%}

%{
    %% Another example of the use of option usefmin.
    % In this example we specify as solver to use fminunc
    rng(1000)
    % Generate regression data
    outsim=simulateLM(100,'R2',0.6);
    yori = outsim.y;
	X = outsim.X;
    % Transform in a different way positive and negative values
    laPos=0.4;
    laNeg=-0.9;
    y=normYJpn(yori,[],[laPos laNeg],'inverse',true);
    % Use solver to find MLE of laPos and laNeg
    usefmin=struct;
    % specify maximum number of iterations
    usefmin.MaxIter=100;
    % Note that to specify as solver fminunc the optmization toolbox is
    % needed.
    usefmin.solver='fminunc';
    out=boxcoxR(y,X,'family','YJpn','plots',1,'usefmin',usefmin);
    disp('MLE of laPos and laNeg')
    disp(out.lahat)
    % Check that the values after the optmization are as expected
    assert(max(abs(out.lahat-[0.3166   -0.8825]))<1e-4,'Wrong values of laP and laN')
%}

%{
    %% Ex using simulated contaminated data.
    rng(10000)
    % Generate X and y data
    n=200;
    X=randn(n,3);
    beta=[ 1; 1; 1];
    sig=0.5;
    ytrue=X*beta+sig*randn(n,1);
    % Contaminate response
    ycont=ytrue;
    ycont(21:40)=ycont(21:40)+4;
    % Use two different values for laP and laN
    lapos=0.5;
    laneg=0;
    ytra=normYJpn(ytrue,[],[lapos laneg],'inverse',true,'Jacobian',false);
    yconttra=normYJpn(ycont,[],[lapos laneg],'inverse',true,'Jacobian',false);
    % In this example the true values of laP and laN are 0.5 and 0 however due
    % to contamination the MLE of lambda laP because very close to 0.
    % This wrongly suggests a unique value of lambda.
    subplot(2,1,1)
    out=boxcoxR(ytra,X,'family','YJpn','plots',1);
    subplot(2,1,2)
    outcont=boxcoxR(yconttra,X,'family','YJpn','plots',1);
%}

%% Beginning of code
persistent cachedlahatPreviousStep;

nnargin=nargin;
vvarargin=varargin;
[y,X,n] = chkinputR(y,X,nnargin,vvarargin);

if nargin<2
    error('FSDA:boxCoxR:missingInputs','It is necessary to supply both y and X');
end

% Specify the default values for the input option
laseq= -2:0.001:2;
laseqPos = -2:0.01:2;
laseqNeg = -2:0.01:2; % laseqPos;
plots=0;
family='BoxCox';
conflev=0.95;
usefmin=false;


% Write in structure 'options' the options chosen by the user
if nargin > 2
    
    options=struct('family',family,'plots',0,'laseq',laseq,'nocheck',0,...
        'intercept',1,'conflev',conflev,'laseqPos',laseqPos,...
        'laseqNeg',laseqNeg,'usefmin',usefmin);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:boxcoxR:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    laseq=options.laseq;
    plots=options.plots;
    family=options.family;
    conflev=options.conflev;
    laseqPos=options.laseqPos;
    laseqNeg=options.laseqNeg;
    usefmin=options.usefmin;
end


if strcmp(family,'BoxCox')
    BoxCoxTra=1;
    
elseif strcmp(family,'YJ')
    BoxCoxTra=2;
    
elseif strcmp(family,'YJpn')
    BoxCoxTra=3;
    
    MaxIter      = 1000;
    TolX         = 1e-6;
    typemin      = 1;
    displayLevel = 'off'; % no display during iterations
    
    if isstruct(usefmin)
        useOptim=true;
        if isfield(usefmin,'MaxIter') == true
            MaxIter=usefmin.MaxIter;
        end
        
        if isfield(usefmin,'TolX') == true
            TolX=usefmin.TolX;
        end
        
        if  isfield(usefmin,'solver') == true
            solver=usefmin.solver;
            if strcmp(solver,'fminsearch')==1
                typemin=1;
            elseif strcmp(solver,'fminunc')==1
                % Check if minimization toolbox is installed in current computer
                typemin=exist('fminunc','file');
                if typemin==0
                    warning('FSDA:boxcoxR:OptToolNotPresent','fminunc has been chosen as optmizer but optmization toolbox is not present')
                    warning('FSDA:boxcoxR:OptToolNotPresent','fminsearch will be used')
                    typemin=1;
                end
            else
                error('FSDA:boxcoxR:WrongInptOption',['solver in input structure' ...
                    ' usefmin can only be ''fminsearch'' or  ''fminunc'''])
            end
        end
        if isfield(usefmin,'displayLevel') == true
            displayLevel=usefmin.displayLevel;
        end
        
    elseif islogical(usefmin)
        useOptim=usefmin;
    else
        error('FSDA:boxcoxR:WrongType','Input option usefmin must be lofical or struct')
    end
else
    warning('FSDA:boxcoxR:WrongFamily','Transformation family which has been chosen is not supported')
    error('FSDA:boxcoxR:WrongFamily','Supported values are BoxCox or YeoJohnson')
end


[~, R] = qr(X,0);
E = X/R;
A = -E*E';
sel=1:n;
siz = size(A);
% Find linear indexes
% It is better to compute linind directly rather than calling sub2ind
% linind=sub2ind(siz,sel,sel);
linind = sel + (sel - 1).*siz(1);

% linind=sub2ind(size(A),sel,sel);
A(linind)=1+A(linind);
% Notice that:
% -E*E' = matrix -H = -X*inv(X'X)*X' computed through qr decomposition
% A = Matrix I - H


if BoxCoxTra == 1
    
    LogLik=[laseq(:) zeros(length(laseq),1)];
    ij=0;
    
    % Add the extra check on vector y
    if min(y)<0
        error('FSDA:boxcoxR:ynegative','BoxCox family cannot be computed because min(y) is smaller than 0. Please use Yeo-Johnson family')
    end
    
    SumLogY=sum(log(y));
    
    for la=laseq
        ij=ij+1;
        if la~=0
            ytra=(y.^la-1)/la;
        else
            ytra=log(y);
        end
        % Residual sum of squares using y(lambda) divided by n
        sigma2hat=ytra'*A*ytra/n;
        
        % logJ = log of the Jacobian
        logJ=(la-1)*SumLogY;
        
        
        %  Value of Profile Log likelihood
        LogLik(ij,2)=-0.5*n*log(sigma2hat)+logJ;
    end
    
elseif BoxCoxTra ==2
    
    
    LogLik=[laseq(:) zeros(length(laseq),1)];
    ij=0;
    
    nonnegs=y>=0;
    SumLogY=sum( log(   (1 + abs(y)).^(2 * nonnegs - 1)) );
    
    for la=laseq
        ij=ij+1;
        % Yeo and JOhnson transformation (without the Jacobian)
        ytra=normYJ(y,1,la,'Jacobian',false);
        
        % Residual sum of squares using y(lambda) divided by n
        sigma2hat=ytra'*A*ytra/n;
        
        % logJ = log of the Jacobian
        logJ=(la-1)*SumLogY;
        
        %  Value of Profile Log likelihood
        LogLik(ij,2)=-0.5*n*log(sigma2hat)+logJ;
    end
    
elseif BoxCoxTra ==3 % This is the case of two values of lambda
    
    nonnegs=y>=0;
    negs=~nonnegs;
    ynonnegs=y(nonnegs);
    ynegs=y(negs);
    logynonnegsp1=log(ynonnegs+1);
    log1mynegs=log(1-ynegs);
    SumLogYp=sum(logynonnegsp1);
    SumLogYn=sum(-log1mynegs);
    % cachedlahatPreviousStep=[];
    
    % If computeLogLikUsingGrid is true compute the likelihood using a grid of
    % values of laseqPos and laseqNeg
    OptimTruePreviousLambdaFalse=useOptim==true && isempty(cachedlahatPreviousStep);
    computeLogLikUsingGrid=OptimTruePreviousLambdaFalse || useOptim==false;
    
    if computeLogLikUsingGrid == true
        % if useOptim==true and the previous estimater of lambda  has not been cached
        % redefine laseqPos and laseqNeg in order to find
        % a rough initial estimate of lambda for the optmization
        if OptimTruePreviousLambdaFalse
            maxL=2;
            laseqPos=-maxL:0.5:maxL;
            laseqNeg=-maxL:0.5:maxL;
        end
        
        LogLik=zeros(length(laseqPos),length(laseqNeg));
        
        ijlaPos=0;
        for laPos=laseqPos
            ijlaPos=ijlaPos+1;
            ijlaNeg=0;
            ytra=y;
            % Yeo and Johnson transformation for positive observations (without the Jacobian)
            % YJ transformation is the Box-Cox transformation of
            % y+1 for nonnegative values of y
            % ytra(nonnegs)=normYJ(y(nonnegs),1,laPos,'Jacobian',false);
            if laPos ~=0
                % ytra(nonnegs)= ((ynonnegs+1).^laPos-1)/laPos;
                % ytra(nonnegs)= (laPos*log(ynonnegs+1) -1)/laPos;
                ytra(nonnegs)= (exp(laPos*logynonnegsp1)-1)/laPos;
            else
                % ytra(nonnegs)= log(ynonnegs+1);
                ytra(nonnegs)= logynonnegsp1;
            end
            
            
            for laNeg=laseqNeg
                ijlaNeg=ijlaNeg+1;
                
                % Yeo and Johnson transformation for negative values (without the Jacobian)
                % YJ transformation is the Box-Cox transformation of
                %  |y|+1 with parameter 2-lambda for y negative.
                % ytra(negs)=normYJ(y(negs),1,laNeg,'Jacobian',false);
                laNegm2=laNeg-2;
                if laNegm2~=0
                    % Slower version
                    % ytra(negs) = - ((-y(negs)+1).^(2-laNeg)-1)/(2-laNeg);
                    % Faster version
                    % ytra(negs)= -( exp( (2-laNeg)* log(-ynegs+1)) -1)/(2-laNeg);
                    % Even faster version
                    % ytra(negs)= -( exp((2-laNeg)* log1mynegs) -1 )/(2-laNeg);
                    % Even, even  faster version
                    ytra(negs)= ( exp(-laNegm2* log1mynegs) -1 )/laNegm2;
                else
                    % ytra(negs) = -log(-ynegs+1);
                    ytra(negs) = -log1mynegs;
                end
                
                % logJ = log of the Jacobian
                logJ=(laPos-1)*SumLogYp +(laNeg-1)*SumLogYn;
                
                % Residual sum of squares using y(lambda) divided by n
                sigma2hat=ytra'*A*ytra/n;
                % sigma2hat=sum(ytra.*(A*ytra))/n;
                
                %  Value of Profile Log likelihood
                LogLik(ijlaPos,ijlaNeg)=-0.5*n*log(sigma2hat)+logJ;
                
                % Alternative implementation using normalized transformed
                % observations (z(lambda))
                % ytra1=ytra/((exp(logJ))^(1/n));
                % sigma2hat1=ytra1'*A*ytra1/n;
                % LogLik(ijlaPos,ijlaNeg)=-0.5*n*log(sigma2hat1);
            end
        end
    else
        LogLik=[];
    end
end


out=struct;

if BoxCoxTra <=2
    % Find best Lambda:
    [maxLogLik,maxLoglikind]=max(LogLik(:,2));
    lahat=laseq(maxLoglikind);
    
    quant=chi2inv(conflev,1)/2;
    [maxLoglik,maxLoglikind]=max(LogLik(:,2));
    intersectPoint=maxLoglik-quant;
    indLow=find(LogLik(:,2)>intersectPoint,1,'first');
    lambdaLow=LogLik(indLow,1);
    indUp=find(LogLik(maxLoglikind:end,2)>intersectPoint,1,'last');
    lambdaUp=LogLik(indUp+maxLoglikind,1);
    lahatci=[lambdaLow  lambdaUp];
    
    % Store MLE of lambda
    out.lahat=lahat;
    % Confidence intervals for MLE of lambda
    out.lahatci=lahatci;
    
    % Plot of profile loglikelihood
    if plots==1
        plot(laseq,LogLik(:,2));
        AxisValues=axis;
        hold on
        % plot a dotted line showing the best lambda of the supplied range
        plot([lahat lahat],[AxisValues(3) maxLogLik],':');
        
        % Plot the confidence interval for lambda
        coo=axis;
        line(lambdaLow*ones(2,1),[coo(3) LogLik(indLow,2)],'Color','r'),
        line(lambdaUp*ones(2,1),  [coo(3) LogLik(indUp+maxLoglikind,2)],'Color','r'),
        title([num2str(conflev*100) ' per cent confidence interval for \lambda'])
        vdisp=(coo(4)-coo(3))/20;
        text(lambdaLow,coo(3)+vdisp,num2str(lambdaLow,'%2.2f'))
        text(lambdaUp,coo(3)+vdisp,num2str(lambdaUp,'%2.2f'))
        xlabel('\lambda');
        ylabel('Profile Log-Likelihood');
        stext=sprintf('%2.0f%% confidence interval for \\lambda',round(100*(conflev)));
        title(stext)
    end
else % This is the case of 2 values of lambda
    
    if computeLogLikUsingGrid == true
        % Find row and column index of max of LogLik
        [row, col] = find(ismember(LogLik, max(LogLik(:))));
        % (approximate) MLE of lambda
        lahat=[laseqPos(row(1)), laseqNeg(col(1))];
        
        % Show the contour plot of the loglikelihood
        if plots==1
            [Xlapos,Ylaneg] = meshgrid(laseqNeg,laseqPos);
            contour(Xlapos,Ylaneg,LogLik); % ,'ShowText','on')
            xlabel('\lambda_N')
            ylabel('\lambda_P')
            text(laseqNeg(col),laseqPos(row),'x')
            title(['$\hat \lambda_P=' num2str(lahat(1)) ', \hat \lambda_N=' num2str(lahat(2)) '$'] ,...
                'Interpreter','latex','FontSize',14)
        end
        
    else
        lahat=cachedlahatPreviousStep;
    end
    
    
    
    % Start the minimization
    if useOptim==true
        
        nlinfitOptions=statset('Display',displayLevel,'MaxIter',MaxIter,'TolX',TolX);
        
        % theta0=lahat+1e-3*randn(1,2);
        theta0=lahat;
        likfminOneParameter = @(laBoth)likfmin(laBoth, y, A, SumLogYp, SumLogYn, nonnegs, negs, n, logynonnegsp1, log1mynegs);
        if typemin == 1
            [lahatfinal,~,exitflag]  = ...
                fminsearch(likfminOneParameter,theta0,nlinfitOptions);
        else
            [lahatfinal,~,exitflag] = ...
                fminunc(likfminOneParameter,theta0,nlinfitOptions);
        end
        % Put in cache the value of lahatfinal and store it
        cachedlahatPreviousStep=lahatfinal;
        % Store MLE of lambda
        out.lahat=lahatfinal;
        out.exitflag=exitflag;
    else
        % Store MLE of lambda
        out.lahat=lahat;
    end
end

% Store profile Log lik
out.LogLik=LogLik;

end


% likfmin = Objective function to call with fminunc or fminsearch when useOptim==true
function objyhat=likfmin(laBoth, y, A, SumLogYp, SumLogYn, nonnegs, negs, n, logynonnegsp1, log1mynegs)

laPos=laBoth(1);
laNeg=laBoth(2);
ytra=y;
% Yeo and Johnson transformation for positive observations (without the Jacobian)
if laPos ~=0
    % ytra(nonnegs)= ((ynonnegs+1).^laPos-1)/laPos;
    % ytra(nonnegs)= (exp(laPos*log(ynonnegs+1))-1)/laPos;
    ytra(nonnegs)= (exp(laPos*logynonnegsp1)-1)/laPos;
else
    % ytra(nonnegs)= log(ynonnegs+1);
    ytra(nonnegs)= logynonnegsp1;
end

% Yeo and Johnson transformation for negative values (without the Jacobian)
if 2-laNeg~=0
    % ytra(negs)= -( exp( (2-laNeg)* log(-ynegs+1)) -1)/(2-laNeg);
    ytra(negs)= -( exp((2-laNeg)*log1mynegs) -1)/(2-laNeg);
else
    % ytra(negs) = -log(-ynegs+1);
    ytra(negs) = -log1mynegs;
end

% logJ = log of the Jacobian
logJ=(laPos-1)*SumLogYp +(laNeg-1)*SumLogYn;

% Residual sum of squares using y(lambda) divided by n
% sigma2hat=ytra'*A*ytra/n;
sigma2hat=sum(ytra.*(A*ytra))/n;
%  Value of objective function to minimize (negative log likelihood)
objyhat=0.5*n*log(sigma2hat)-logJ;
end

%FScategory:REG-Transformations
