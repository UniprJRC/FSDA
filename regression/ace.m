function [out]=ace(y,X,varargin)
%ace computes alternative conditional expectation
% 
%<a href="matlab: docsearchFS('ace')">Link to the help page for this function</a>
%
%   This function uses the alternating conditional expectations algorithm
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
%           specifies how the type of transformation for the exaplanatory
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
% See also: aceplot.m
%
% References:
%
%
% Breiman, L. and Friedman, J.H. (1985), Estimating optimal
% transformations for multiple regression and correlation, "Journal of the
% American Statistical Association", Vol. 80, pp. 580-597.
% Wang D.  and Murphy M. (2005), Identifying nonlinear relationships
% regression using the ACE algorithm, "Journal of Applied Statistics",
% Vol. 32, pp. 243-258.
%
% Copyright 2008-2018.
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
Xw=tX.*sqrt(w);
yw=ty.*sqrt(w);
b=Xw\yw;
tX=tX.*(b');

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
    
    % ($nit$ = counter for number of iterations in the inner loop).
    nit=0;
    
    lfinishInnerLoop=1;
    
    % Iterate until e^2(\theta, \phi_1, ..., \phi_p) fails to decrease
    while lfinishInnerLoop==1 % Beginning of Inner Loop
        rsqi=rsq;
        nit=nit+1;
        
        % Find $Z_5 = ty- tX_1 \cdots - tX_p$.
        %  $Z_5$ is the $n \times 1$ vector which contains the residuals
        %  from transformed $y$ and transformed $X$. In what follows $Z_1,
        %  \ldots, Z_5$ are scratch variables (each of dimension $n \times
        %  1$).
        z5=ty-sum(tX,2);
        
        % Loop through all the explanatory variables
        % and transform them
        %   do 420 i=1,p (in Fortran code)
        for  j=1:p
            
            % sorted indexes of the j-th column of X
            ordXj=M(:,j);
            
            %   Find $Z_1= Z_5 + tX_i$ (note that the observations of all the matrices involved are ordered using the ith column of matrix $M$).
            %   Therefore to be precise   $Z_1= Z_5(M(:,i)) + tX_i(M(:,i))$.
            %   $Z_1$ is the $n \times 1$ vector
            %   which contains the residuals from transformed $y$ and transformed $X$ excluding
            %   variable $X_i$ with values ordered according to $X_i$.
            %   \[
            %   Z_1=  ty(M(:,i))- \sum_{j \ne i } tX_j(M(:,i)), \qquad j=1, 2, \ldots, p
            %   \]
            %    The first element of $Z_1$ is associated with the lowest value of $X_i$, the second element is associated with the second lowet value of $X_i$ ...,
            z1=z5(ordXj)+tX(ordXj,j);
            
            % z2 = original expl. variable X(:,j) sorted
            z2=X(ordXj,j);
            
            % z4 = corresponding weights
            z4=w(ordXj);
            
            %  smo=smothr(abs(l(i)),z(:,2),z(:,1),z(:,4));
            % Find smoothed values of X(:,i)
            smo=smothr(abs(l(j)),z2,z1,z4);
            
            % Weighted average of
            sm=sum(smo.*z4)/sw;
            % z3 = z3- mean(z3)
            z3=smo-sm;
            % sv=sum(((z(:,1)-z(:,3)).^2).*w(m(:,i)));
            
            % Compute residual sum of squares
            sv=sum(((z1-z3).^2).*z4);
            
            % Convert residual sum of squares into R2
            % Remember that y has been standardized
            sv=1-sv/sw;
            
            if sv <= rsq
                % lfinishInnerLoop=0;
            else
                % rsq = the multiple R-squared value for the transformed values
                rsq=sv;
                
                tX(ordXj,j)=z3;
                z5(ordXj)=z1-z3;
            end
            
        end       % End of the loop for the explanatory variables
        
        % Condition to exit the Inner Loop
        % There is just one explnatory variable ||
        % The change in R2 is smaller than delrsq
        % The maximum number of iterations has been achieved
        if (p == 1 || rsq-rsqi <= delrsq || nit == maxit)
            lfinishInnerLoop=0;
        end
    end % End of the Inner Loop
    
    
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
end

