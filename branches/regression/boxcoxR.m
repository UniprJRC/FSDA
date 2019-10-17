function out=boxcoxR(y,X, varargin)
%boxcoxR finds MLE of lambda in linear regression (and confidence interval) using Box Cox or YJ transformation
%
%<a href="matlab: docsearchFS('boxcoxR')">Link to the help function</a>
%
% It computes the profile log Likelihood for a range of values of the
% transforamtion parameter (lambda) and computes the MLE of lambda in the
% supplied range. Supported families are Box Cox and Yeo and Johnson.
%               The profile loglikelihood is computed as:
%               \[
%               -(n/2)*\log( (y(\lambda)-X\beta(\lambda))'(y(\lambda)-X\beta(\lambda))/n) +\log J
%               \]
%               where  $y(\lambda)$ is the vector of transformed
%               observations using Box Cox family or Yeo and Johnson family
%               \[
%               \beta(\lambda) = (X'X)^{-1} X' y(\lambda)
%               \]
%               and $J$ is the Jacobian of the transformation,
%
%  Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :          Predictor variables. Matrix. Matrix of explanatory
%               variables (also called 'regressors') of dimension n x (p-1)
%               where p denotes the number of explanatory variables
%               including the intercept.
%               Rows of X represent observations, and columns represent
%               variables. By default, there is a constant term in the
%               model, unless you explicitly remove it using input option
%               intercept, so do not include a column of 1s in X. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations.
%
%
%
% Optional input arguments:
%
%  intercept :  Indicator for constant term. Scalar.
%               If 1, a model with
%               constant term will be fitted (default), else no constant
%               term will be included.
%               Example - 'intercept',1
%               Data Types - double
%    family :   parametric transformation to use. String. String which
%               identifies the family of transformations which
%               must be used. Character. Possible values are 'BoxCox'
%               (default) or 'YJ'.
%               The Box-Cox family of power transformations equals
%               $(y^{\lambda}-1)/\lambda$ for $\lambda$ not equal to zero,
%               and $\log(y)$ if $\lambda = 0$.
%               The Yeo-Johnson (YJ) transformation is the Box-Cox
%               transformation of $y+1$ for nonnegative values, and of
%               $|y|+1$ with parameter 2-lambda for y negative.
%               The basic power transformation returns $y^{\lambda}$ if
%               $\lambda$ is not zero, and $\log(\lambda)$  otherwise.
%               Remark. BoxCox and the basic power family can be used just
%               if input y is positive. Yeo-Johnson family of
%               transformations does not have this limitation.
%               Example - 'family','YJ'
%               Data Types - char
%       nocheck : Check input arguments. Scalar. If nocheck is equal to 1
%                 no check is performed on
%                 vector y and matrix X. Notice that y and X are left
%                 unchanged. In other words the additional column of ones
%                 for the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
%               Example - 'laseq',-2:0.01:2
%               Data Types - double
%        conflev : Confidence level for lambda. Scalar.
%                  scalar between 0 and 1 determining
%                  confidence level for
%                  lambda based on the asymptotic chi1^2 of
%                  twice the loglikelihood ratio.
%                  The default confidence level is 0.95;
%               Example - 'conflev',0.99
%               Data Types - double
%          laseq : Sequence of values of lambda to consider. Vector. Vector
%                  which contains the sequence of values of lambda for
%                  which the profile loglikelihood has to be computed.
%                  The default value of laseq is -2:0.001:2;
%               Example - 'laseq',[-1:0.001;0.7]
%               Data Types - double
%   plots  :    Profile log likelihood for lambda. Boolean.
%               It specifies whether it is necessary to show the profile
%               log likelihood of lambda If plots is true, the plot of the
%               profile loglikelihood is produced together with the
%               requested confidence interval. The default value of prolik
%               is false, that is no plot is produced.
%               Example - 'plots',true
%               Data Types - boolean
%
% Output:
%
%         out:   structure which contains the following fields
%
% out.lahat  = scalar containing best estimate of lambda.
%               This is the best estimate among the values of lambda which
%               are supplied in input vector laseq.
% out.lahatci  = confidence intervals for MLE of lambda computed using Chi2
%               approximation using confidence level specified in input
%               option conflev.
% out.LogLik   = vector with the same length of input option laseq
%               containing the value of the profile loglikelihood for each
%               value in laseq.
%
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
%
% Acknowledgements:
%
% This function has been inspired by sumbmission
% https://www.mathworks.com/matlabcentral/fileexchange/10419-box-cox-power-transformation-for-linear-models
% in the file exchange written by Hovav Dror, hovav@hotmail.com, March 2006
%
% Copyright 2008-2019.
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

%% Beginning of code

nnargin=nargin;
vvarargin=varargin;
[y,X,n] = chkinputR(y,X,nnargin,vvarargin);

if nargin<2
    error('FSDA:boxCoxR:missingInputs','It is necessary to supply both y and X');
end

% Specify the defaule values for the input option
laseq= -2:0.001:2;
plots=0;
family='BoxCox';
conflev=0.95;

% Write in structure 'options' the options chosen by the user
if nargin > 2
    
    options=struct('family',family,'plots',0,'laseq',laseq,'nocheck',0,'intercept',1,'conflev',conflev);
    
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
end


if strcmp(family,'BoxCox')
    BoxCoxTra=true;
    
elseif strcmp(family,'YJ')
    BoxCoxTra=false;
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

LogLik=[laseq(:) zeros(length(laseq),1)];
ij=0;

if BoxCoxTra == true
    % Add the extra check on vector y
    if min(y)<0
        error('FSDA:Score:ynegative','BoxCox family cannot be computed because min(y) is smaller than 0. Please use Yeo-Johnson family')
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
    
else
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
    
end

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

out=struct;
% Store MLE of lambda
out.lahat=lahat;
% Confidence intervals for MLE of lambda
out.lahatci=lahatci;
% Store profile Log lik
out.LogLik=LogLik;
end





