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
%       l :  type of transformation. Vector. Vector of length p which
%           specifies how the type of transformation for the explanatory
%           variables. 
%           l(j)=1 => j-th variable assumes orderable values.
%           l(j)=2 => j-th variable assumes circular (periodic) values
%                 in the range (0.0,1.0) with period 1.0.
%           l(j)=3 => j-th variable transformation is to be monotone.
%           l(j)=4 => j-th variable transformation is to be linear.
%           l(j)=5 => j-th variable assumes categorical (unorderable) values.
%           j =1, 2, \ldots, p+1.
%           The default value of l is is a vector of ones of length p+1,
%           that is the procedure assumes that both the explanatory
%           variables and the response have orderable values.
%           Example - 'l',[3 3 1]
%           Data Types - double
%       w  : weights for the observations. Vector. Row or column vector of
%           length n containing the weights associated to each
%           observations. If w is not specified we assum $w=1$ for $i=1,
%           2, \ldots, n$.
%           Example - 'w',1:n
%           Data Types - double
%   nterm  : minimum number of consecutive iteration below the threshold
%           to terminate the outer loop. Positive scalar. This value
%           specifies how many consecutive iterations below the threshold
%           it is necesasry to have to declare convergence in the outer
%           loop. The default value of nterm is 3.
%           Example - 'nterm',5
%           Data Types - double
%   delrsq : termination threshold. Scalar. Iteration (in the outer loop)
%            stops when rsq changes less than delrsq in nterm. The default
%            value of delrsq is 0.01.
%           Example - 'delrsq',0.001
%           Data Types - double
%    maxit : maximum number of iterations for the outer loop. Scalar. The
%            default maximum number of iterations before exiting the outer
%            loop is 20.
%           Example - 'maxit',30
%           Data Types - double
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
% Copyright 2008-2018.
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

% c span, alpha : super smoother parameters.
% supermo=struct;


UserOptions=varargin(1:2:length(varargin));

if ~isempty(UserOptions)
    
    options=struct('l',l,'delrsq',delrsq,'nterm',nterm,...
        'w',w,'maxit',maxit);
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSReda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
end

if size(w,2)>1
    w=w';
end

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
meX=sum(X.*w,1)/sw;
tX=X-meX;

% Initial transformation for matrix X.
% X is transformed so that its columns are equally weighted when predicting y.
% Xw=tX.*sqrt(w);
% yw=ty.*sqrt(w);
% b=Xw\yw;
% tX=tX.*(b');
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

% Call backfitting algorith to find transformed values of X
[tX,rsq]=backfitAVAS(ty,tX,X,w,M,l,rsq,maxit,sw,p,delrsq);

yspan=0;
lfinishOuterLoop=1;

while lfinishOuterLoop ==1 % Beginning of Outer Loop
    iter=iter+1;
    
    % z1 contains fitted values sorted
    z1=sum(tX,2);    % z1 is z10 in fortran program
    
    % tres = vector of residuals using transformed y and transformed X
    tres=ty-z1;
    % If a squared residual is exactly zero leads to a problem when taking
    % logs therefore. Replace 0 with 1e-10
    tres(abs(tres)<1e-10)=1e-10;
    % z2 = log of sqrt of transformed squared residuals
    % z2 = log absolute values of residuals
    z2=log(sqrt(tres.^2));
    
    [z1,ordyhat]=sort(z1);
    % ztar_sorted=log(sqrt((z2-z1).^2));
    ztar_sorted=z2(ordyhat);
    z5=w(ordyhat);
    
    % Now the residuals  are smoothed
    % smo=smothr(abs(l(pp1)),z(:,2),z(:,1),z(:,4));
    % x coord = fitted values based on tX (sorted)
    % y values ztar_sorted (log of |residuals| using the ordering of
    % fitted values)
    
    % Smooth the log of (the sample version of) the |residuals|
    % against fitted values. Use the ordering based on fitted values.
    [smo,yspan]=rlsmo(z1,ztar_sorted,z5,yspan);
    
    % z6=smo;
    z7=exp(-smo);
    sumlog=2*n*sum(smo.*w)/sw;
    z8=ty(ordyhat);
    
    % Compute the variance stabilizing transformation
    % h(t)= \int_{z1_1}^t v(u)^{-0.5} du
    % z1 = yhat sorted
    % z7 reciprocal of smoothed |residuals|. z7 is v(u)^{-0.5} 
    % z8 = values of ty sorted using the ordering of z1
    z9=ctsub(z1,z7,z8);
    
    sm=sum(z9.*w);
    ty(ordyhat)=z9-sm/sw;
    
    sv=sum(ty.^2.*w)/sw;
    % Make sure that sv (variance of transformed y) is a positive number
    if sv<=0
        out=[];
        warning('FSDA:avas:negsv','Return a missing value')
        return
    end
    
    svx=sum(z1.^2.*w)/sw;
    ty=ty/sqrt(sv);
    tX=tX/sqrt(svx);
    % Get the new transformation of X_j
    % that is backfit \hat g(y) on X_1, \ldots, X_p
    % to get new tX
    
    [tX,rsq]=backfitAVAS(ty,tX,X,w,M,l,rsq,maxit,sw,p,delrsq);
    % z1 contains fitted values sorted
    z1=sum(tX,2);    % z1 is z10 in AVAS
    
    rr=sum(((ty-z1).^2).*w)/sw;
    
    % New value of R2 
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
end

