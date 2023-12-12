function [out]=ace(y,X,varargin)
%ace computes alternative conditional expectation
%
%<a href="matlab: docsearchFS('ace')">Link to the help page for this function</a>
%
%   This function uses the alternating conditional expectation algorithm
%   to find the transformations of y and X that maximise the proportion of
%   variation in y explained by X. When X is a matrix, it is transformed so
%   that its columns are equally weighted when predicting y.
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
%       l :  type of transformation. Vector. Vector of length p+1 which
%           specifies how the type of transformation for the explanatory
%           variables and the response. The fitst p element of this vector
%           refer to the p explanatory variables, the last element refers to
%           the response.
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
%
%       w  : weights for the observations. Vector. Row or column vector of
%           length n containing the weights associated to each
%           observations. If w is not specified we assum $w=1$ for $i=1,
%           2, \ldots, n$.
%           Example - 'w',1:n
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
%   delrsq : termination threshold. Scalar. Iteration (in the outer loop)
%            stops when rsq changes less than delrsq in nterm. The default
%            value of delrsq is 0.01.
%           Example - 'delrsq',0.001
%           Data Types - double
%
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
%      out.outliers = k x 1 vector containing the units declared as outliers
%           when procedure is called with input option rob set to true. If
%           rob is false out.outliers=[].
%
% See also: aceplot.m, smothr.m, supsmu.m
%
% References:
%
%
% Breiman, L. and Friedman, J.H. (1985), Estimating optimal transformations
% for multiple regression and correlation, "Journal of the American
% Statistical Association", Vol. 80, pp. 580-597. 
% 
% Wang D.  and Murphy M. (2005), Identifying nonlinear relationships regression using the ACE
% algorithm, "Journal of Applied Statistics", Vol. 32, pp. 243-258.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('ace')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:

%{
   %% Example of the use of ace based on the Wang and Murphy data.
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
    % Apply the ace algorithm
    out= ace(y,X);
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


[varargin{:}] = convertStringsToChars(varargin{:});
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

sqrtw=sqrt(w);
yw=ty.*sqrtw;

% Center X matrix
% Initial transformation for matrix X.
% X is transformed so that its columns are equally weighted when predicting y.
verLess2016b=verLessThanFS('9.1');
if verLess2016b == false
    meX=sum(X.*w,1)/sw;
    tX=X-meX;
    Xw=tX.*sqrtw;
    b=Xw\yw;
    tX=tX.*(b');
else
    meX=sum(bsxfun(@times,X,w),1)/sw;
    tX=bsxfun(@minus,X,meX);
    Xw=bsxfun(@times,tX,sqrtw);
    b=Xw\yw;
    tX=bsxfun(@times,tX,(b'));
end

% In the Fortran program vector b is found iteratively
% (procedure using subroutine scail)
% [tXchk,bchk]=scail(w,sw,ty,tX,eps,maxit);

% ct vector containing the
ct=zeros(maxit,1);
ct(1:nterm)=100;


% M matrix of dimension n-by-(p+1) containing the ranks
% (ordered positions) for each of the p explanatory variables and for the
% response.
pp1=p+1;
M=zeros(n,p+1);
for j=1:p
    [~,pos]=sort(X(:,j),1);
    M(:,j)=pos;
end

[~,pos]=sort(y,1);
M(:,pp1)=pos;

% TODELETE
% M(:,1)=[1, 5, 4, 3, 2, 6, 7, 8, 9, 10];

% \item initialize $iter=0$ (counter for number of iteration for outer
% loop), $maxit=20$ (maximum number of iterations for inner and outer
% loop), $rsq=0$ (initial value of $R^2$)

iter=0;
nt=0;
rsq=0;

lfinishOuterLoop=1;

while lfinishOuterLoop ==1 % Beginning of Outer Loop
    iter=iter+1;
    
    % Call backfitting algorithm to find transformed values of X
    tX=backfit(ty,tX,X,w,M,l,rsq,maxit,sw,p,delrsq);
    
    ordy=M(:,pp1);
    z2=y(ordy);
    z4=w(ordy);
    z1=sum(tX(ordy,:),2);
    
    % Now the response is smoothed
    % smo=smothr(abs(l(pp1)),z(:,2),z(:,1),z(:,4));
    smo=smothr(l(pp1),z2,z1,z4);
    
    z3=smo;
    
    
    sm=sum(z3.*w(ordy))/sw;
    z2(ordy)=z1;
    
    z3=z3-sm;
    sv=sum((z3.^2).*z4)/sw;
    
    if sv<=0
        out=[];
        warning('FSDA:ace:negsv','Return a missing value')
        return
    else
        sv=1.0/sqrt(sv);
        
        ty(ordy)=z3*sv;
        
        sv=sum(((ty-z2).^2).*w);
        
        rsq=1.0-sv/sw;
        nt=mod(nt,nterm)+1;
        ct(nt)=rsq;
        cmn=100.0;
        cmx=-100.0;
        cmn=min([cmn; ct(1:nterm)]);
        cmx=max([cmx;ct(1:nterm)]);
    end
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
out.outliers=[];
end

%FScategory:REG-Transformations
