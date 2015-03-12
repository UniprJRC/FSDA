function [X,id]=simdataset(n, Pi, Mu, S,varargin)
%simdataset simulates a dataset given the parameters of finite mixture model with Gaussian components
%
%
%<a href="matlab: docsearchFS('simdataset')">Link to the help function</a>
%
%   simdataset(n, Pi, Mu, S) generates a matrix of size n-by-p containing n
%   observations p dimensions from k groups. In other words, this function
%   produces a dataset of n observations from a mixture model with
%   parameters 'Pi' (mixing proportions), 'Mu' (mean vectors), and 'S'
%   (covariance matrices). Mixture component sample sizes are produced as a
%   realization from a multinomial distribution with probabilities given by
%   mixing proportions. For example, if n=200, k=4 and Pi=(0.25, 0.25,
%   0.25, 0.25) function Nk1=mnrnd( n-k, Pi) is used to generate k integer
%   numbers (whose sum is n-k) from the multinominal distribution with
%   parameters n-k and Pi. The size of the groups is given by Nk1+1. The
%   first Nk1(1)+1  observations are generated using centroid Mu(1,:) and
%   covariance S(:,:,1), ..., the last Nk1(k)+1  observations are generated
%   using centroid Mu(k,:) and covariance S(:,:,k)
%
%  Required input arguments:
%
%         n   : scalar, sample size  of the dataset
%        Pi   : vector of length(k) defining mixing proportions. \sum_{j=1}^k Pi=1
%        Mu   : k-by-v matrix containing components' mean vectors
%         S   : v-by-v-by-k arrary containing components' covariance matrices
%
%  Optional input arguments:
%
%       nnoise : scalar, which specifies the number of noise variables (the
%               default value of nnoise is zero).
%         nout : scalar, which specifies the number of outlying observations.
%                The default value of nout is 0. If nout is (for example
%                10) than n+nout observations are generated.
%       alpha  : level for simulating outliers. The default value of alpha
%                is 0.001,
%       maxiter: maximum number of trials to simulate outliers. The default
%                value of maxout is 1e+05
%       int    : vector or string.
%                If int is a vector of length 2 it contains min and maximum
%                values of the interval in which noise has to be simulated
%                It int is empty (default) noise and outliers are simulated
%                uniformly between the smallest and largest coordinates of
%                mean vectors.
%                If int='minmax' noise and outliers are simulated uniformly
%                between the smallest and largest coordinates of simulated
%                data matrix X
%       lambda : vector of length p containing inverse Box-Cox
%                transformation coefficients. The value false (default)
%                implies that no transformation is applied to any variable.
%
%
%  Output:
%
%           X  : simulated dataset of size (n + nout)-by-(p + nnoise);
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
% transformation providing a vector of coefficients 'lambda'. The value 1
% implies that no transformation is needed for the corresponding
% coordinate
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('simdataset')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    out = MixSim(4,2,'BarOmega',0.01);
    n=60;
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S);
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'nout',10);
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'nout',10,'int','minmax');


    out = MixSim(4,3,'BarOmega',0.1);
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'nout',100);

%}


%{
    % Check if it is possible to generate the outliers given alpha  
    % in the interval specified by input option int and if it not possible
    % modify int by adding 0.1 until the outliers can be found
    out = MixSim(4,3,'BarOmega',0.1);
    % Point mass contamination of 30 observations in interval [0.4 0.4]
    nout=30;
    outint=[0.4 0.4];
    n=200;
    for j=1:10
        [X,id]=simdataset(n, out.Pi, out.Mu, out.S, 'nout', nout, 'alpha', 0.01, 'int', outint);
        if size(X,1)== n+nout
            break
        else
            outint=outint+0.1;
        end
    end
    spmplot(X,id)
%}

%% Beginning of code

if (n < 1)
    error('FSDA:simdataset:Wrongn','Wrong sample size n...')
end

if sum(Pi <= 0)~=0 || sum(Pi >= 1) ~= 0
    error('FSDA:simdataset:WrongPi','Wrong vector of mixing proportions Pi: the values must be in the interval (0 1)')
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
        error('FSDA:simdataset:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:simdataset:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
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
    error('FSDA:simdataset:Wrongnnoise','Wrong value of nnoise: it cannot be smaller than 0')
end

if (nout < 0)
    error('FSDA:simdataset:Wrongnout','Wrong value of nout: it cannot be smaller than 0')
end

if ((alpha >= 1) || (alpha <= 0))
    error('FSDA:simdataset:WrongAlpha','Wrong value of alpha: it must be in the interval (0 1)')
end

if (maxiter < 1)
    error('FSDA:simdataset:WrongMaxIter','Wrong value for maximum number of iterations: it cannot be <1')
end

[k,p]=size(Mu);

if (n >= k)
    
    if R_seed
        %mrr = rmultinom(1, n - K, Pi)
        nmk = n-k;
        putRdata('dd',nmk);
        putRdata('Pi',Pi);
        rn1m = 'rmultinom(1,dd,Pi)';
        mrr = evalR(rn1m);
        mrr = double(mrr');
        % Another solution, but more complicated
        % cPi = num2str(Pi');
        % cPi = regexprep(cPi, '[ ]{2,}', ',');
        % rn1m = ['rmultinom(1,' num2str(n - k) ', c(' cPi '))'];
        % mrr = evalR(rn1m);
        % mrr = double(mrr');
    else
        mrr = mnrnd( n-k, Pi);
    end
    % Nk contains the sizes of the clusters
    Nk = ones(1,k)+mrr;
else
    error('FSDA:simdataset:Wrongn','Sample size (n) cannot be less than the number of clusters (k)')
end

X=zeros(n,p);
id=zeros(n,1);
for j=1:k
    a=sum(Nk(1:j-1))+1;
    b=sum(Nk(1:j));
    id(a:b)=j*ones(Nk(j),1);
    
    if R_seed
        % mvrnorm(n = Nk[j], mu = Mu[j, ], Sigma = S[, , j])
        % mvrnorm(n = Nk[k], mu = Mu[k, ], Sigma = S[,, k])
        Muj = Mu(j,:);
        Sj  = S(:,:,j);
        Nkj = Nk(j);
        Xab = [];
        putRdata('Muj',Muj);
        putRdata('Sj',Sj);
        putRdata('Nkj',Nkj);
        putRdata('Xab',Xab);
        % Make sure that function mvrnorm is present
        exists=evalR('exists("mvrnorm",mode="function")');
        if ~exists
            % If exists is not true try to install library MASS which
            % includes function mvrnorm
            evalR('library(MASS)')
        end
        mvrnorms = 'Xab <- mvrnorm(n = Nkj, mu = Muj, Sigma = Sj)';
        evalR(mvrnorms);
        Xab = getRdata('Xab');
        if isempty(Xab)
            error('FSDA:simdataset:MissingRlibrary','Could not load library(MASS) in R, please install it')
        end
        X(a:b,:) = Xab;
    else
        X(a:b,:) = mvnrnd(Mu(j,:),S(:,:,j),Nk(j));
    end
end

if nout ~= 0
    [Xout, fail] = getOutliers(nout, Mu, S, alpha, maxiter,int);
    if fail == 1
        warning('FSDA:simdataset:Modifiedn',['Output matrix X will have just ' num2str(n) ...
            ' rows and not ' num2str(n+nout)])
    else
        X =[X;Xout];
        id =[id;zeros(nout,1)];
    end
end

if nnoise ~= 0
    if isempty(int)
        L = min(min(Mu));
        U = max(max(Mu));
    elseif strcmp('minmax',int)
        L = min(X);
        U = max(X);
    else
        L = int(1);
        U = int(2);
    end
    
    
    if R_seed
        % equivalent of 'rand(k,p)' in R is 'matrix(runif(k*p),k,p)'
        rn1s = ['matrix(runif(' num2str((n + nout)*nnoise) '),' num2str(n + nout) ',' num2str(nnoise) ')'];
        rrr = evalR(rn1s);
    else
        rrr = rand(n + nout,  nnoise);
    end
    
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
                warning('FSDA:simdataset:NaNs','NaNs were produced during transformation')
            end
        end
    else
        error('FSDA:simdataset:WrongLambda','The number of transformation coefficients lambda should be equal to ndimensions + nnoise')
    end
    % Remark: if lambda is very large MATLAB can create a complex X, so it
    % is necessary to have this further check
    X=real(X);
    
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