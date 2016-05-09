function [y,X,id]=simdatasetreg(n, Pi, Beta, S, Xdistrib, varargin)
%simdatasetReg simulates regression dataset given the parameters of a mixture regression model
%
%
%<a href="matlab: docsearchFS('simdatasetreg')">Link to the help function</a>
%
%   [y,X,id]=simdatasetreg(n, Pi, Beta, S) generates a regression dataset
%   of size n from a mixture model with parameters 'Pi' (mixing
%   proportions), 'Beta' (matrix of regression coefficients), and 'S'
%   (vector of variances of the distributions of the points around each
%   regression hyperplane). Component sample sizes are produced as a
%   realization from a multinomial distribution with probabilities given by
%   mixing proportions. For example, if n=200, k=4 and Pi=(0.25, 0.25,
%   0.25, 0.25) function Nk1=mnrnd( n-k, Pi) is used to generate k integer
%   numbers (whose sum is n-k) from the multinominal distribution with
%   parameters n-k and Pi. The size of the groups is given by Nk1+1. The
%   first Nk1(1)+1  observations are generated using vector of regression
%   coefficients Beta(:,1) and variance S(1), ..., and the X simulated as
%   specified in structure Xdistrib, the last Nk1(k)+1 observations are
%   generated using using vector of regression coefficients Beta(:,k),
%   variance S(k) and the X simulated as specified in structure Xdistrib
%
%
%  Required input arguments:
%
%         n   : scalar, sample size  of the dataset
%        Pi   : vector of length k defining mixing proportions. \sum_{j=1}^k Pi=1
%      Beta   : p-by-k matrix containing (in the columns) regression
%               coefficients for the k groups
%         S   : vector of length k containing the variances of the k
%               regression hyperplanes
%    Xdistrib : structure which contains information about how to generate
%               each explanatory variable inside each group
%               The following options are admitted for Xdistrib
%                   Xdistrib.intercept = scalar equal to 1 if intercept is
%                       present. The default value of Xdistrib.intercept is 1.
%               The other fields of Xdistrib depend on the distribution
%               which is chosen.
%               NORMAL DISTRIBUTION N(mu, sigma)
%                   Xdistrib.type='normal';
%                   Xdistrib.mu = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters mu
%                       for each explanatory variable and each group. The
%                       default value of Xdistrib.mu is zeros(p-1, k).
%                   Xdistrib.sigma = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters
%                       sigma for each explanatory variable and each group.
%                       The default value of Xdistrib.sigma is ones(p-1, k)
%               UNIFORM DISTRIBUTION U(a, b)
%                   Xdistrib.type='uniform';
%                   Xdistrib.a = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters a
%                       for each explanatory variable and each group. The
%                       default value of Xdistrib.a is zeros(p-1, k).
%                   Xdistrib.b = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters b
%                       for each explanatory variable and each group. The
%                       default value of Xdistrib.b is ones(p-1, k).
%               HALF NORMAL DISTRIBUTION Half(sigma)= |N(0 sigma)|
%                   Xdistrib.type='halfnormal';
%                   Xdistrib.sigma = matrix of size (p-1)-by-k if
%                   (Xdistrib.intercept=1) or p-by-k if (Xdistrib.intercept=0)
%                   containing the parameters sigma for each explanatory variable
%                   and each group. The default value of Xdistrib.sigma is
%                   ones(p-1, k).
%TODO:simdatasetReg:OTHER_DISTRIB
%                   Xdistrib.type='user'.
%                   Xdistrib.X = matrix with at least n rows and p-1 (if
%                   intercept is present) or p (if intercept is not
%                   present) columns containing the values of the
%                   explanatory variables for the k groups.
%                   Xdistrib.id =identifier vector which labes the rows of
%                   matrix Xdistrib.X
%
%  Optional input arguments: (TODO TODO)
%
%TODO:simdatasetReg:INPUT_OPTIONS
%
%
%  Output:
%
%           y  : (n+nout)-by-1 vector containing the values of the
%                responses for the k groups
%           X  : matrix of size (n + nout)-by-(p + nnoise) containinng the
%                values of the explanatory variables for the k groups
%                Noise coordinates are provided in the last nnoise columns.
%           id : classification vector of length n + nout; 0 represents an
%                outlier.
%
%            REMARK: If nout outliers could not be generated a warning is
%                produced. In this case matrix X and vector id will have
%                just n rows.
%
%   DETAILS
% To make a dataset more challenging for clustering, a user might want to
% simulate noise variables or outliers. Parameter 'nnoise' specifies the
% desired number of noise variables. If an interval 'int' is specified,
% noise will be simulated from a Uniform distribution on the interval given
% by 'int'. Otherwise, noise will be simulated uniformly between the
% smallest and largest coordinates of mean vectors. 'nout' specifies the
% number of observations outside (1 - 'alpha') ellipsoidal contours for the
% weighted component distributions. Outliers are simulated on a hypercube
% specified by the interval 'int'. A user can apply an inverse Box-Cox
% transformation of y providing a coefficient 'lambda'. The value 1
% implies that no transformation is needed for the response
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('simdatasetreg')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:
%
%{
    % Example 1
    % Generate mixture of regression using an average overlapping at
    % centroids =0.01. Use all default options
    % 1) Beta is generated according to random normal for each group with
    % mu=0 and sigma=1
    % 2) X in each dimension and each group is generated according to U(0, 1)
    % 3) regression hyperplanes contain intercepts
    p=5;
    k=3;
    Q=MixSimreg(k,p,'BarOmega',0.01);
    n=200;
    % Q.Xdistrib.BarX in this case has dimension 5-by-3 and is equal to
    % 1.0000    1.0000    1.0000
    % 0.5000    0.5000    0.5000
    % 0.5000    0.5000    0.5000
    % 0.5000    0.5000    0.5000
    % 0.5000    0.5000    0.5000
    % Probabilities of overlapping are evaluated at
    % Q.Beta(:,1)'*Q.Xdistrib.BarX(:,1) ... Q.Beta(:,3)'*Q.Xdistrib.BarX(:,3)
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    spmplot([y X(:,2:end)],id)
%}


%% Beginning of code

if (n < 1)
    error('FSDA:simdatasetreg:Wrongn','Wrong sample size n...')
end

if sum(Pi <= 0)~=0 || sum(Pi >= 1) ~= 0
    error('FSDA:simdatasetreg:WrongPi','Wrong vector of mixing proportions Pi: the values must be in the interval (0 1)')
end

nnoisedef=0;
noutdef=0;
alphadef=0.001;
intdef='';
maxiterdef=1e+05;
lambdadef='';
Rseeddef = 0;

options=struct('nnoise',nnoisedef,'nout',noutdef,'alpha',alphadef,'int',intdef,...
    'maxiter',maxiterdef,'lambda',lambdadef,'R_seed', Rseeddef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:simdatasetreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:simdatasetreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the options chosen by the user
    for ii=1:2:length(varargin);
        options.(varargin{ii})=varargin{ii+1};
    end
    
end

R_seed   = options.R_seed;
nnoise=options.nnoise;
nout=options.nout;
alpha=options.alpha;
int=options.int;
maxiter=options.maxiter;
lambda=options.lambda;


if (nnoise < 0)
    error('FSDA:simdatasetreg:Wrongnnoise','Wrong value of nnoise: it cannot be smaller than 0')
end

if (nout < 0)
    error('FSDA:simdatasetreg:Wrongnout','Wrong value of nout: it cannot be smaller than 0')
end

if ((alpha >= 1) || (alpha <= 0))
    error('FSDA:simdatasetreg:WrongAlpha','Wrong value of alpha: it must be in the interval (0 1)')
end

if (maxiter < 1)
    error('FSDA:simdatasetreg:WrongMaxIter','Wrong value for maximum number of iterations: it cannot be <1')
end

[p,k]=size(Beta);

if (n >= k)
    mrr = mnrnd( n-k, Pi);
    
    % Nk contains the sizes of the clusters
    Nk = ones(1,k)+(mrr);
else
    error('FSDA:simdatasetreg:Wrongn','Sample size (n) cannot be less than the number of clusters (k)')
end


y=zeros(n,1);
id=zeros(n,1);

% Check if intercept is present
intercept=Xdistrib.intercept;
if intercept==1
    X=[ones(n,1) zeros(n,p-1)];
else
    X=zeros(n,p);
end

for j=1:k
    aa=sum(Nk(1:j-1))+1;
    bb=sum(Nk(1:j));
    id(aa:bb)=j*ones(Nk(j),1);
    
    % The dimension of a and b is p-1-by-k (if intercept is 1)
    % or p-by-k (if intercept is 0).
    if find(strcmp('Uniform',Xdistrib.type))
        if intercept==1
            Xab=bsxfun(@times, rand(Nk(j),p-1),(Xdistrib.b(:,j)-Xdistrib.a(:,j))');
            Xab=bsxfun(@plus,Xab,Xdistrib.a(:,j)');
            X(aa:bb,2:end)=Xab;
        else
            Xab=bsxfun(@times, rand(Nk(j),p),(Xdistrib.b(:,j)-Xdistrib.a(:,j))');
            Xab=bsxfun(@plus,Xab,Xdistrib.a(:,j)');
            X(aa:bb,:)=Xab;
        end
        
    elseif find(strcmp('Normal',Xdistrib.type))
        
        if intercept==1
            Xab=bsxfun(@times, randn(Nk(j),p-1),(Xdistrib.sigma(:,j))');
            Xab=bsxfun(@plus,Xab,Xdistrib.mu(:,j)');
            X(aa:bb,2:end)=Xab;
        else
            Xab=bsxfun(@times, randn(Nk(j),p),(Xdistrib.sigma(:,j))');
            Xab=bsxfun(@plus,Xab,Xdistrib.mu(:,j)');
            X(aa:bb,:)=Xab;
        end
        
    elseif find(strcmp('HalfNormal', Xdistrib.type))
        if intercept==1
            Xab=bsxfun(@times, abs(randn(Nk(j),p-1)),(Xdistrib.sigma(:,j))');
            X(aa:bb,2:end)=Xab;
        else
            Xab=bsxfun(@times, abs(randn(Nk(j),p)),(Xdistrib.sigma(:,j))');
            X(aa:bb,:)=Xab;
        end
    elseif find(strcmp('User',Xdistrib.type))
        d=find(strcmp(fieldnames(Xdistrib),'X'),1);
        if ~isempty(d);
            X= Xdistrib.X;
        else
            error('FSDA:simdatasetreg:MissingField','If string Xdistrib = ''User'' then the user must provide input matrix X')
        end
    else
        error('FSDA:simdatasetreg:Wrongbetadistrib','Possible values for option betadistrib are ''Normal'' ''Uniform'' ''HalfNormal'' and ''User'' ')
    end
    
    
    
    y(aa:bb)=X(aa:bb,:)*Beta(:,j)+sqrt(S(j))*randn(Nk(j),1);
    
    
    %    X(a:b,:) = rand(Nk(j),p)*BarX(j)*2;
    %     y(a:b)=X(a:b,:)   +sqrt(S(j))*randn(Nk(j),1);
    %        X(a:b,:) = randn(Nk(j),p) +Mu(j);
    
end

if nout ~= 0
    [Xout, fail] = getOutliers(nout, Beta, S, alpha, maxiter,int);
    if fail == 1
        warning('FSDA:simdatasetreg:Modifiedn',['Output matrix X will have just ' num2str(n) ...
            ' rows and not ' num2str(n+nout)])
        
    else
        X =[X;Xout];
        id =[id;zeros(nout,1)];
    end
end

if nnoise ~= 0
    if isempty(int)
        L = min(min(Beta));
        U = max(max(Beta));
    elseif strcmp('minmax',int)
        L = min(X);
        U = max(X);
    else
        L = int(1);
        U = int(2);
    end
    
    
    rrr = rand(n + nout,  nnoise);
    
    % Xnoise = (U-L)*rrr+L;
    if nnoise>length(U);
        U=repmat(U,10);
        L=repmat(L,10);
    end
    
    Xnoise=bsxfun(@times,rrr,(U(1:nnoise)-L(1:nnoise)));
    Xnoise=bsxfun(@plus,Xnoise,L(1:nnoise));
    X = [X, Xnoise];
    
end

if ~isempty(lambda)
    if (length(lambda) == p + nnoise)
        for j=1:(p + nnoise)
            X(:,j) = (lambda(j) * X(:, j) + 1).^(1/lambda(j)) - 1;
            if (sum(isnan(X(:,j))) ~= 0)
                warning('FSDA:simdatasetreg:NaNs','NaNs were produced during transformation')
            end
        end
    else
        error('FSDA:simdatasetreg:WrongLambda','The number of transformation coefficients lambda should be equal to ndimensions + nnoise')
    end
end
%% Inner functions
% Xout with nout rows which contains the outliers
% fail = scalar. If fail =1 than it was not possible to generate the
% outliers in the interval specified by input option int in maxiter trials
% else fail =0
    function [Xout,fail] = getOutliers(nout, Mu, S, alpha, maxiter,int)
        fail = 0;
        % maxiter = maximum number of iterations to generate outliers
        
        critval =chi2inv(1-alpha,p);
        
        Xout = zeros(nout,p);
        
        if isempty(int)
            L = min(min(Mu));
            U = max(max(Mu));
        elseif strcmp(int,'minmax')
            L = min(X);
            U = max(X);
            
        else
            L = int(1);
            U = int(2);
        end
        
        i = 1;
        
        iter=0;
        while (i <= nout  &&  iter<maxiter)
            iter=iter+1;
            if R_seed
                % equivalent of 'rand(k,p)' in R is 'matrix(runif(k*p),k,p)'
                rn1 = ['matrix(runif(' num2str(p) '),' num2str(1) ',' num2str(p) ')'];
                rr = evalR(rn1);
            else
                rr = rand(1,p);
            end
            
            Xout(i,:) = (U-L).*rr+L;
            
            ij=0;
            for jj=1:k
                if mahalFS(Xout(i,:),Mu(jj,:),S(:,:,jj)) <critval
                    ij=1;
                    break
                end
            end
            
            if ij==0
                i = i + 1;
            end
        end
        % If iter = maxiter than it was not possible  to generate nout
        % outliers in maxiter simulations.
        if iter== maxiter
            disp(['Warning: it was not possible to generate ' num2str(nout) ' outliers'])
            disp(['in ' num2str(maxiter) ' replicates in the interval [' num2str(L(1)) ...
                '--' num2str(U(1)) ']'])
            disp('Please modify the interval inside input option ''int'' ')
            disp('or increase input option ''alpha''')
            disp(['The values of int and alpha now are ' num2str(int) ' and ' num2str(alpha)]);
            
            % If max number of iteration has been reached fail is 1
            fail=1;
            Xout=Xout(1:i,:);
        end
    end
end
%FScategory:CLUS-MixSim