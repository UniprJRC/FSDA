function [out]=avas(y,X,varargin)
%avas computes additivity and variance stabilization for regression
%
%<a href="matlab: docsearchFS('avas')">Link to the help page for this function</a>
%
%   This function differs from ace in that it uses a (nonparametric)
%   variance-stabilizing transformation for the response variable.
%
%
% Required input arguments:
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
%  Optional input arguments:
%
%   delrsq : termination threshold. Scalar. Iteration (in the outer loop)
%            stops when rsq changes less than delrsq in nterm. The default
%            value of delrsq is 0.01.
%           Example - 'delrsq',0.001
%           Data Types - double
%
%       l :  type of transformation. Vector. Vector of length p which
%           specifies the type of transformation for the explanatory
%           variables.
%           l(j)=1 => j-th variable assumes orderable values.
%           l(j)=2 => j-th variable assumes circular (periodic) values
%                 in the range (0.0,1.0) with period 1.0.
%           l(j)=3 => j-th variable transformation is to be monotone.
%           l(j)=4 => j-th variable transformation is to be linear.
%           l(j)=5 => j-th variable assumes categorical (unorderable) values.
%           j =1, 2, \ldots, p+1.
%           The default value of l is a vector of ones of length p,
%           that is the procedure assumes that both the explanatory
%           variables and the response have orderable values. Note that in
%           avas procedure the reponse is always transformed
%           Example - 'l',[3 3 1]
%           Data Types - double
%
%    maxit : maximum number of iterations for the outer loop. Scalar. The
%            default maximum number of iterations before exiting the outer
%            loop is 20.
%           Example - 'maxit',30
%           Data Types - double
%
%   nterm  : minimum number of consecutive iteration below the threshold
%           to terminate the outer loop. Positive scalar. This value
%           specifies how many consecutive iterations below the threshold
%           it is necesasry to have to declare convergence in the outer
%           loop. The default value of nterm is 3.
%           Example - 'nterm',5
%           Data Types - double
%
% PredictorOrderR2 : inclusion of the variables using R2. Boolean.
%           The default backfitting algorithm of avas (PredictorOrderR2=false)
%           just does one
%           iteration over the predictors, it may not find optimal
%           transformations and it will be dependent on the order of the
%           predictors in X. If PredictorOrderR2 in the backfitting
%           explanatory variables become candidate for transformation
%           depending on their R2 order of importance.
%           Example - 'PredictorOrderR2',true
%           Data Types - logical
%
%   rob  : Use outlier detection in each step fo the bakfitting procedure.
%           Boolean or struct.
%           If rob is not specified tyinitial is set to false and no
%           outlier detection routine is called.
%           If rob is true, FSR routine is called to detect outliers
%           before calling backfitting procedure and outliers are removed
%           If rob is a struct it is possible to specify in the
%           fields method, bdp and alphasim the values of the estiamtor to
%           use, the breakdown point and and confidence level to use to
%           declare the outliers.
%           More precisely:
%           rob.method = character which identifies the estiamtor to use.
%               Possible values are 'FS' (Forward search estimator', 'LTS'
%               (Least trimmed squares estimator), 'LMS' (Least median of
%               squares), 'S' (S estimator). If this field is not specified
%               'FS' is used.
%           rob.bdp = a scalar in the interval (0, 0.5] which specifies
%               breakdown point to use. If this field is not specified it is
%               set to 0.2.
%           rob.simalpha = simultaneous error rate one is willing to
%               tolerate to declare the outliers. If rob.simalpha is not
%               specified it is set to 0.01.
%           Example - 'rob',true
%           Data Types - boolean or struct
%
%   scail  : Inizialing values for the regressors. Boolean. If scail is
%           true (default value is false), a linear regression is done on
%           the independent variables $X_1$, ...$X_{p-1}$ centered, using y
%           standardized, getting the coefficients $b_1$, ..., $b_{p-1}$.
%           AVAS if scail is true takes the initial value of the transform
%           of $X_j$ to be $b_j(X_j-mean(X_j))$, $j=1, \ldots, p-1$. This
%           is done in order to limit the effect of the ordering of the
%           variables on the final results. The purpose of this initial
%           linear transformation is to let the the columns of matrix X
%           have the same weight when predicting $y$.
%           Note that this initialization option is always used inside the
%           ace routine but it is set to false in AVAS or compatibility
%           with the original fortran program of Tibshirani.
%           See the ADDED NOTE in the [Monotone Regression Splines in
%           Action]: Comment.  Leo Breiman Statistical Science, Vol. 3, No.
%           4 (Nov., 1988), pp. 442-445.
%           Example - 'scail',true
%           Data Types - logical
%
%   tyinitial  : Initial values for the transformed response. Boolean or struct.
%           If tyinitial is not specified tyinitial is set to false and the initial value for ty are
%           simply the standardized values.
%           If tyinitial is true, y is transformed using best BoxCox lambda
%           value among the five most common values of lambda (-1, -0.5, 0,
%           0.5, 1) and routines FSRfan.m and fanBIC.m are called.
%           If tyinitial is a struct it is possible to specify in the
%           fields la and family the values of lambda and the family to
%           use. More precisely:
%           tyinitial.la = values of the lambdas to consider inside FSRfan.
%           tyinitial.family= string which identifies the family of transformations which
%                   must be used. Character. Possible values are 'BoxCox'
%                   (default), 'YJ', or 'YJpn'.
%                   The Box-Cox family of power transformations equals
%                   $(y^{\lambda}-1)/\lambda$ for $\lambda$ not equal to zero,
%                   and $\log(y)$ if $\lambda = 0$.
%                   The Yeo-Johnson (YJ) transformation is the Box-Cox
%                   transformation of $y+1$ for nonnegative values, and of
%                   $|y|+1$ with parameter $2-\lambda$ for $y$ negative.
%                   Remember that BoxCox can be used just
%                   if input y is positive. Yeo-Johnson family of
%                   transformations does not have this limitation.
%                   If family is 'YJpn' Yeo-Johnson family is applied but in
%                   this case it is also possible to monitor (in the output
%                   arguments out.Scorep and out.Scoren) the score test
%                   respectively for positive and negative observations.
%           Example - 'tyinitial',true
%           Data Types - boolean or struct
%
%       w  : weights for the observations. Vector. Row or column vector of
%           length n containing the weights associated to each
%           observations. If w is not specified we assum $w=1$ for $i=1,
%           2, \ldots, n$.
%           Example - 'w',1:n
%           Data Types - double
%
%
% Output:
%
%         out:   structure which contains the following fields
%      out.ty  = n x 1 vector containing the transformed y values.
%      out.tX  = n x p matrix containing the transformed X matrix.
%     out.rsq  = the multiple R-squared value for the transformed values in
%               the last iteration of the outer loop.
%      out.y  = n x 1 vector containing the original y values.
%      out.X  = n x p matrix containing the original X matrix.
%      out.niter = scalar. Number of iterations which have been necessary
%               to achieve convergence.
%      out.outliers = k x 1 vector containing the units declared as outliers
%           when procedure is called with input option rob set to true. If
%           rob is false out.outliers=[].
%
% See also: ace.m, aceplot.m, smothr.m, supsmu.m, ctsub.m
%
% References:
%
%
% Tibshirani R. (1987), Estimating optimal transformations for regression,
% "Journal of the American Statistical Association", Vol. 83, 394-405.
% Wang D.  and Murphy M. (2005), Identifying nonlinear relationships
% regression using the ACE algorithm, "Journal of Applied Statistics",
% Vol. 32, pp. 243-258.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('avas')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:

%{
   %% Example of the use of avas based on the Wang and Murphy data.
   % In order to have the possibility of replicating the results in R using
   % library acepack function mtR is used to generate the random data.
    rng('default')
    seed=11;
    negstate=-30;
    n=200;
    X1 = mtR(n,0,seed)*2-1;
    X2 = mtR(n,0,negstate)*2-1;
    X3 = mtR(n,0,negstate)*2-1;
    X4 = mtR(n,0,negstate)*2-1;
    res=mtR(n,1,negstate);
    % Generate y
    y = log(4 + sin(3*X1) + abs(X2) + X3.^2 + X4 + .1*res );
    X = [X1 X2 X3 X4];
    % Apply the avas algorithm
    out= avas(y,X);
    % Show the output graphically using function aceplot
    aceplot(out)
%}

%{
    % Example 1 from TIB88: brain body weight data.
    % Comparison between ace and avas.
    YY=load('animals.txt');
    y=YY(1:62,2);
    X=YY(1:62,1);
    out=ace(y,X);
    aceplot(out)

    out=avas(y,X);
    aceplot(out)
    % https://vincentarelbundock.github.io/Rdatasets/doc/robustbase/Animals2.html
    % ## The `same' plot for Rousseeuw's subset:
    % data(Animals, package = "MASS")
    % brain <- Animals[c(1:24, 26:25, 27:28),]
    % plotbb(bbdat = brain)
%}

%{
    %% Example 3 from TIB88: unequal variances.
    n=200;
    rng('default')
    seed=100;
    negstate=-30;
    x = mtR(n,0,seed)*3;
    z = mtR(n,1,negstate);
    y=x+0.1*x.*z;
    X=x;
    nr=1;
    nc=2;
    outACE=ace(y,X);
    outAVAS=avas(y,X);
    subplot(nr,nc,1)
    yy=[y, outACE.ty outAVAS.ty];
    yysor=sortrows(yy,1);
    plot(yysor(:,1),yysor(:,2:3))

    % plot(y,[outACE.ty outAVAS.ty])
    title('Compare ACE and AVAS')
    xlabel('y')
    ylabel('ty')
    legend({'ACE' 'AVAS'},'Location','Best')

    subplot(nr,nc,2)
    XX=[X, outACE.tX outAVAS.tX];
    XXsor=sortrows(XX,1);
    plot(XXsor(:,1),XXsor(:,2:3))

    % plot(y,[outACE.ty outAVAS.ty])
    title('Compare ACE and AVAS')
    xlabel('x')
    ylabel('tX')
    legend({'ACE' 'AVAS'},'Location','Best')
    % For R users, the R code to reproduce the example above is given
    % below.
    % set.seed(100)
    % n=200
    % x <- runif(n)*3
    % z <- rnorm(n)
    % y=x+0.1*x*z;
    % out=ace(x,y)
    % out=avas(x,y)

%}

%{
    % Example 4 from TIB88: non constant underlying variance.
    close all
    negstate=-100;
    rng('default')
    seed=100;
    n=200;
    x1=mtR(n,1,seed);
    x2=mtR(n,1,negstate);
    z=mtR(n,1,negstate);
    % x1=randn(n,1);
    % x2=randn(n,1);
    % z=randn(n,1);
    absx1_x2=abs(x1-x2);
    y=x1+x2+(1/3)*absx1_x2.*z;
    X=[x1 x2];
    out=avas(y,X);
    tyhat=sum(out.tX,2);
    ty=out.ty;
    absty_tyhat=abs(ty-tyhat);
    % hold('on')
    plot(absx1_x2,absty_tyhat,'o');
    lsline
    % For R users, the R code to reproduce the example above is given
    % below.
    % # Example 4 non constant underlying variance
    %  set.seed(100)
    % n=200;
    % x1=rnorm(n);
    % x2=rnorm(n);
    % z=rnorm(n);
    % absx1_x2=abs(x1-x2);
    % y=x1+x2+(1/3)*absx1_x2*z;
    % X=cbind(x1,x2);
    % out=avas(X,y);
%}

%{
    %% Example 5 from TIB88: missing group variable.
    rng('default')
    seed=1000;
    n=100;
    x=mtR(n,0,seed)*10-5;
    z=mtR(n,1,-seed);
    %x=rand(n,1)*10-5;
    %z=randn(n,1);
    y1=3+5*x+z;
    y2=-3+5*x+z;
    y=[y1;y2];
    X=[x;x];
    out=avas(y,X);
    aceplot(out);
    % For R users, the R code to reproduce the example above is given
    % below.
    % # Example 5 missing group variable
    % set.seed(1000)
    % n=100;
    % x=runif(n)*10-5;
    % z=rnorm(n);
    % y1=3+5*x+z;
    % y2=-3+5*x+z;
    % y=c(y1,y2);
    % X=c(x,x)
    % out=avas(X,y);
%}

%{
    % Example 6 from TIB88: binary.
    seed=20;
    n=50;
    y1=exp(-0.5+0.5*mtR(n,1,seed));
    y2=exp(0.5+0.5*mtR(n,1,-seed));
    y=[y1;y2];
    X=[-0.5*ones(n,1); 0.5*ones(n,1)];
    out=avas(y,X);
    aceplot(out)
    % For R users, the R code to reproduce the example above is given
    % below.
    % Example 6 binary
    % set.seed(20)
    % n=50;
    % y1=exp(-0.5+0.5*rnorm(n));
    % y2=exp(0.5+0.5*rnorm(n));
    % y=c(y1,y2);
    % X=c(-0.5*rep(1,n), 0.5*rep(1,n));
    % out=avas(X,y);

%}

%{
    % Example 9 from TIB88: Nonmonotone function of X.
    n=200;
    x=rand(n,1)*2*pi;
    z=randn(n,1);
    y=exp(sin(x)+0.2*z);
    X=x;
    out=avas(y,X);
    aceplot(out)
%}


%{
    % Example of use of option tyinitial.
    % Generate the data.
    rng(2000)
    n=100;
    X=10*rand(n,1);
    sigma=0.1;
    a=2;
    b=0.3;
    y=a+b*X+sigma*randn(n,1);
    % The correct transformation is la=-0.5 (inverse square root)
    la=-0.5;
    y=normBoxCox(y,1,la,'inverse',true);
    % call of AVAS without option tyinitial
    subplot(2,1,1)
    outAVAS=avas(y,X,'l',4);
    tyfinal=outAVAS.ty;
    out=fitlm(X,tyfinal);
    plot(X,tyfinal,'o')
    lsline
    xlabel('Original x')
    ylabel('Without option tyinitial')
    title(['R2=' num2str(out.Rsquared.Ordinary)])
    % call of AVAS with option tinitial set to true
    subplot(2,1,2)
    outAVAStyini=avas(y,X,'l',4,'tyinitial',true);
    tyfinal=outAVAStyini.ty;
    out=fitlm(X,tyfinal);
    plot(X,tyfinal,'o')
    lsline
    xlabel('Original x')
    ylabel('With option tyinitial')
    title(['R2=' num2str(out.Rsquared.Ordinary)])
%}

%% Beginning of code

if nargin <2
    error('FSDA:ace:missingInputs','A required input argument is missing.')
end

[n,p]=size(X);

% l specifies how to transform the variables
% The first p values of l refer to the p explanatory variables.
% The last refers to the response
% 4 = linear transformation
% 3 = monotonic transformation
% ........
l=ones(p+1,1);
%  termination threshold for outer loop
delrsq=0.01;
% maxit = max. number of iterations for the outer loop
maxit=20;
% nterm : number of consecutive iterations for which
%         rsq must change less than delrsq for convergence.
nterm=3;
w=ones(n,1);

scail=false;
tyinitial=false;
rob=false;
PredictorOrderR2=false;

% c span, alpha : super smoother parameters.
% supermo=struct;


UserOptions=varargin(1:2:length(varargin));

if ~isempty(UserOptions)
    
    options=struct('l',l,'delrsq',delrsq,'nterm',nterm,...
        'w',w,'maxit',maxit,'scail',scail,'tyinitial',tyinitial,'rob',rob,...
        'PredictorOrderR2',PredictorOrderR2);
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:avas:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % We now overwrite inside structure options the default values with
    % those chosen by the user
    % Notice that in order to do this we use dynamic field names
    for j=1:2:length(varargin)
        options.(varargin{j})=varargin{j+1};
    end
    
    l=options.l;
    delrsq=options.delrsq;
    w=options.w;
    nterm=options.nterm;
    maxit=options.maxit;
    scail=options.scail;
    tyinitial=options.tyinitial;
    rob=options.rob;
end

if size(w,2)>1
    w=w';
end

if islogical(tyinitial)
    if tyinitial ==true
        callToFSRfan=true;
        la=[-1 -0.5 0 0.5 1];
        if min(y)>0
            family='BoxCox';
        else
            family='YJ';
        end
    else
        callToFSRfan=false;
    end
    
elseif isstruct(tyinitial)
    callToFSRfan=true;
    if isfield('tyinitial','la')
        la=tyinitial.la;
    else
        la=[-1 -0.5 0 0.5 1];
    end
    if isfield('tyinitial','family')
        family=tyinitial.family;
    else
        if min(y)>0
            family='BoxCox';
        else
            family='YJ';
        end
    end
else
    error('FSDA:avas:WrongInputOpt','tyinitial can only be a boolean or a struct.');
end

% Initial default value of break down point
bdpdef=0.2;

if islogical(rob)
    if rob ==true
        robustAVAS=true;
        % FS    ==> estimatorToUse=1 (default)
        % LTS   ==> estimatorToUse=2
        % LMS   ==> estimatorToUse=3
        % S     ==> estimatorToUse=4
        estimatorToUse=1;
        % In this case bdp is init (point to start FS)
        bdp=round(n*(1-bdpdef));
        % simalpa  1 per cent simultaneous
        simalpha=0.01;
    else
        robustAVAS=false;
    end
    
elseif isstruct(rob)
    robustAVAS=true;
    
    if isfield('rob','bdp')
        bdp=rob.bdp;
    else
        bdp=bdpdef;
    end
    if isfield('rob','simalpha')
        simalpha=rob.simalpha;
    else
        simalpha=0.01;
    end
    
    if isfield('rob','method')
        estimator=rob.method;
        if strcmp(estimator,'FS')
            estimatorToUse=1;
            % FS is defined using init and not bdp this is the reason of
            % the instruction below
            bdp=round(n*(1-bdp));
        elseif strcmp(estimator,'LTS')
            estimatorToUse=2;
        elseif strcmp(estimator,'LMS')
            estimatorToUse=3;
        elseif strcmp(estimator,'S')
            estimatorToUse=4;
        else
            error('FSDA:avas:WrongInputOpt',['rob.method can only be '...
                ' ''FS'', ''LTS'', ''LMS'' or ''S''']);
        end
        
    else
        estimatorToUse=1;
    end
else
    error('FSDA:avas:WrongInputOpt','rob can only be a boolean or a struct.');
end



if callToFSRfan ==true
    % FSRfan and fanplot with all default options
    [outFSR]=FSRfan(y,X,'msg',0,'la',la,'family',family,'plots',1);
    outLA=fanBIC(outFSR,'plots',0);
    % outLA.labest=-0.5;
    if strcmp(family,'BoxCox')
        y=normBoxCox(y,1,outLA.labest);
    elseif strcmp(family,'YJ')
        y=normYJ(y,1,outLA.labest);
    elseif strcmp(family,'YJpn')
        y=normYJpn(y,1,outLA.labest);
    end
end

seq=(1:n)';


% sw = sum of the weights
sw=sum(w);
% Center and standardize y
mey=sum(y.*w)/sw;
ty=y-mey;
sv=sum(ty.^2.*w)/sw;
ty=ty/sqrt(sv);
% ty is standardized
% var(ty,1) =1

% Center X matrix
verLess2016b=verLessThanFS(9.1);
if verLess2016b == false
    meX=sum(X.*w,1)/sw;
    tX=X-meX;
else
    meX=sum(bsxfun(@times,X,w),1)/sw;
    tX=bsxfun(@minus,X,meX);
end

if scail==true
    % Initial transformation for matrix X.
    % X is transformed so that its columns are equally weighted when predicting y.
    Xw=tX.*sqrt(w);
    yw=ty.*sqrt(w);
    b=Xw\yw;
    tX=tX.*(b');
end

% In the Fortran program vector b is found iteratively
% (procedure using subroutine scail)

% ct vector containing the
ct=zeros(maxit,1);
ct(1:nterm)=100;


% M matrix of dimension n-by-(p+1) containing the ranks
% (ordered positions) for each of the p explanatory variables
M=zeros(n,p);
for j=1:p
    [~,pos]=sort(X(:,j),1);
    M(:,j)=pos;
end


% \item initialize $iter=0$ (counter for number of iteration for outer
% loop), $maxit=20$ (maximum number of iterations for inner and outer
% loop), $rsq=0$ (initial value of $R^2$)

iter=0;
nt=0;
rsq=0;


if robustAVAS==true
    [bsb,outliers,ngood]=robAVAS(ty,tX,estimatorToUse,bdp,simalpha);
    sw=sum(w(bsb));
else
    outliers=[];
    ngood=n;
    bsb=seq;
end

% Call backfitting algorith to find transformed values of X
[tX,rsq]=backfitAVAS(ty,tX,X,w,M,l,rsq,maxit,sw,p,delrsq,bsb,outliers,PredictorOrderR2);

yspan=0;
lfinishOuterLoop=1;

while lfinishOuterLoop ==1 % Beginning of Outer Loop
    iter=iter+1;
    
    if robustAVAS==true && iter>1
        [bsb,outliers,ngood]=robAVAS(ty,tX,estimatorToUse,bdp,simalpha);
         sw=sum(w(bsb));
    end
    
    % (yhat contains fitted values and yhatord = fitted values sorted
    yhat=sum(tX,2);    % yhat is z10 in fortran program
    
    % tres = vector of residuals using transformed y and transformed X
    tres=ty-yhat;
    % If a squared residual is exactly zero leads to a problem when taking
    % logs therefore. Replace 0 with 1e-10
    tres(abs(tres)<1e-10)=1e-10;
    % logabsres = log of sqrt of transformed squared residuals
    % logabsres = log absolute values of residuals (z2 in the Fortran
    % program)
    logabsres=log(sqrt(tres.^2));
    
    if ~isempty(outliers)
        yhat(outliers)=Inf;
    end
    
    [yhatord,ordyhat]=sort(yhat);
    % ztar_sorted=log(sqrt((z2-z1).^2));
    logabsresOrdyhat=logabsres(ordyhat);
    % wOrdyat = weights using the ordering based on yhat
    wOrdyhat=w(ordyhat);
    
    % Now the residuals  are smoothed
    % smo=smothr(abs(l(pp1)),z(:,2),z(:,1),z(:,4));
    % x coord = fitted values based on tX (sorted)
    % y values ztar_sorted (log of |residuals| using the ordering of
    % fitted values)
    
    % Smooth the log of (the sample version of) the |residuals|
    % against fitted values. Use the ordering based on fitted values.
    % Smoothing is done for the units not declared as outliers.
    [smo,yspan]=rlsmo(yhatord(1:ngood),logabsresOrdyhat(1:ngood),wOrdyhat(1:ngood),yspan);
    
    % smoothresm1Ordyhat = 1/|e_{smoothed}|= v(u)^{-0.5}
    smoothresm1Ordyhat=exp(-smo);
    % sumlog=2*n*sum(smo.*w)/sw;
    
    %     tyOrdyhat=ty(ordyhat);
    %
    %     % Compute the variance stabilizing transformation
    %     % h(t)= \int_{\hat y_{(1)}}^t v(u)^{-0.5} du
    %     % yhatord = yhat sorted (z10 in original fortran code)
    %     % smoothresm1Ordyhat= reciprocal of smoothed |residuals|. z7 is v(u)^{-0.5}
    %     % tyOrdyhat = values of ty sorted using the ordering based on fitted
    %     % values (z8 n fortran code)
    %     % tynewOrdyhat (z9 in fortran code)
    %     tynewOrdyhat=ctsub(yhatord,smoothresm1Ordyhat,tyOrdyhat);
    %
    %     % Probably here there was a mistake in the original Fortran program
    %     % because z9 (which contains values ordered using yhat) has to be
    %     % weighted using z5=w(ordyhat) and not using w
    %     % 23033 continue
    %     %  call ctsub(n,z(1,10),z(1,7),z(1,8),z(1,9))
    %     %  sm=0
    %     %  do 23035 j=1,n
    %     %     sm=sm+w(j)*z(j,9)
    %     % 23035 continue
    %
    %     sm=sum(tynewOrdyhat.*w); % TODO replace w with wOrdyhat
    %     % Compute updated vector ty with mean removed
    %     ty(ordyhat)=tynewOrdyhat-sm/sw;
    
    % Compute updated transformed values
    ty=ctsub(yhatord(1:ngood),smoothresm1Ordyhat(1:ngood),ty);
    
    wbsb=w(bsb);
    sm=sum(ty(bsb).*wbsb);
    % sw=sum(wbsb);
    ty=ty-sm/sw;
    
    
    sv=sum(ty(bsb).^2.*wbsb)/sw;
    % Make sure that sv (variance of transformed y) is a positive number
    if sv<=0
        out=[];
        warning('FSDA:avas:negsv','Return a missing value')
        return
    end
    
    % TODO it is necessary to replace w with  wOrdyhat
    % sw=sum(wOrdyhat(1:ngood));
    svx=sum(yhatord(1:ngood).^2.*wOrdyhat(1:ngood))/sw;
    % ty is the new vector of transformed values standardized.
    ty=ty/sqrt(sv);
    % Note that each column of tX is standardized using sqrt of mean of
    % squared fitted values
    tX=tX/sqrt(svx);
    % Get the new transformation of X_j
    % that is backfit \hat g(y) on X_1, \ldots, X_p
    % to get new tX
    
    [tX,~]=backfitAVAS(ty,tX,X,w,M,l,rsq,maxit,sw,p,delrsq,bsb,outliers,PredictorOrderR2);
    % yhat contains fitted values (not sorted)
    yhat=sum(tX,2);    % z1 is z10 in AVAS
    
    rr=sum(((ty(bsb)-yhat(bsb)).^2).*wbsb)/sw;
    
    % rsq = new value of R2 (given that var(ty)=1)
    rsq=1-rr;
    
    % sumlog=sumlog+n*log(sv);
    % rnew=sumlog+rr;
    
    nt=mod(nt,nterm)+1;
    ct(nt)=rsq;
    cmn=100.0;
    cmx=-100.0;
    cmn=min([cmn; ct(1:nterm)]);
    cmx=max([cmx;ct(1:nterm)]);
    
    % Stopping condition for outer loop, for at least three consecutive
    % times the difference between two consecutive values of $R^2$ is
    % smaller than delrsq or $iter>=maxit$
    if (cmx-cmn <= delrsq  || iter>=maxit)
        % In this case go out of the OuterLoop
        lfinishOuterLoop=0;
    end
end

% Create output structure out
out=struct;
out.y=y;
out.ty=ty;
out.X=X;
out.tX=tX;
out.rsq=rsq;
out.niter=iter;
out.outliers=outliers;
end

function [bsb,outliers,ngood]=robAVAS(ty,tX,estimatorToUse,bdp,simalpha)
n=length(ty);
if estimatorToUse==1 % FS
    if simalpha ==0.01
        out=FSR(ty,tX,'plots',0,'msg',false,'init',bdp);
    else
        out=FSRr(ty,tX,'alpha',simalpha,'plots',0,'msg',false,'init',init);
    end
elseif estimatorToUse==2  % LTS
    out=LXS(ty,tX,'lms',2,'conflev',1-simalpha/n,'bdp',bdp);
elseif estimatorToUse ==3 % LMS
    out=LXS(ty,tX,'lms',1,'conflev',1-simalpha/length(ty),'bdp',bdp);
elseif estimatorToUse ==4 % S
    out=Sreg(ty,tX,'conflev',1-simalpha/length(ty),'bdp',bdp);
else
end
seq=(1:n)';
outl=out.outliers;
if any(isnan(outl)) || isempty(outl)
    ngood=n;
    bsb=seq;
    outliers=[];
else
    ngood=n-length(outl);
    bsb=setdiff(seq,outl);
    outliers=outl;
end
end

%FScategory:REG-Transformations