function [X,id]=simdataset(n, Pi, Mu, S, varargin)
%simdataset simulates a dataset given the parameters of finite mixture model with Gaussian components
%
%<a href="matlab: docsearchFS('simdataset')">Link to the help function</a>
%
%   simdataset(n, Pi, Mu, S) generates a matrix of size n-by-p containing n
%   observations p dimensions from k groups. More precisely, this function
%   produces a dataset of n observations from a mixture model with
%   parameters 'Pi' (mixing proportions), 'Mu' (mean vectors), and 'S'
%   (covariance matrices). Mixture component sample sizes are produced as a
%   realization from a multinomial distribution with probabilities given by
%   the mixing proportions. For example, if n=200, k=4 and Pi=(0.25, 0.25,
%   0.25, 0.25), function Nk1 = mnrnd( n-k, Pi) is used to generate k
%   integers (whose sum is n-k) from the multinominal distribution with
%   parameters n-k and Pi. The size of the groups is Nk1+1. The first
%   Nk1(1)+1  observations are generated using centroid Mu(1,:) and
%   covariance S(:,:,1), ... , the last Nk1(k)+1 observations are generated
%   using centroid Mu(k,:) and covariance S(:,:,k).
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
%    noisevars : empty value, scalar or structure.
%                - If noisevars is not specified or is an empty value
%                  (default) no noise variable is added to the matrix of
%                  simulated data.
%                - If noisevars is a scalar equal to r, then r new noise
%                  variables are added to the matrix of simulated data
%                  using the uniform distribution in the range [min(X)
%                  max(X)].
%                - If noisevars is a strcture it may contain the following
%                  fields:
%                  number: a scalar or a vector of length f. The sum of
%                       elements of vector 'number' is equal to the total
%                       number of noise variables to be addded.
%                  distribution: string or cell array of strings of length
%                       f which specifies the distribution to be used to
%                       simulate the noise variables.
%                       If field distribution is not present then the
%                       uniform distribution is used to simulate the noise
%                       variables.
%                       String 'distribution' can be one of the following
%                       values:
%                       * uniform = uniform distribution
%                       * normal  = normal distribution
%                       * t or T followed by a number which controls the
%                         degrees of freedom. For example, t6 specifies to
%                         generate the data according to a Student T with 6
%                         degrees of freedom.
%                       * chisquare followed by a number which controls the
%                         degreess of freedom. For example, chisquare8
%                         specifies to generate the data according to a Chi
%                         square distribution with 8 degrees of freedom.
%                  interval: string or vector of length 2 or matrix of length
%                         2-by-f which controls for each element of vector
%                         'number' or each element of cell 'distribution',
%                         the min and max of the generated data.
%                         If interval is empty (default), the noise variables
%                         are simulated uniformly between the smallest and
%                         the largest coordinates of mean vectors.
%                         If interval is 'minmax' the noise varaibles are
%                         simulated uniformly between the smallest and the
%                         largest coordinates of the simulated data matrix.
%                For example, the code:
%                   noisevars=struct;
%                   noisevars.number=[3 2];
%                   noisevars.distribution={'chisquare5' 'T3'};
%                   noisevars.interval='minmax';
%                adds 5 noise varibles, the first 3 generated using
%                the chi2 with 5 degrees of freedom and the last two
%                using the Student t with 3 degrees of freedom. Noise
%                variables are generated in the interval min(X) max(X).
%
%   noiseunits : missing value, scalar or structure, specifyin the number
%                and type of outlying observations. The default value of
%                noiseunits is 0.
%                - If noiseunits is a scalar t different from 0, then t
%                  units from the uniform distribution in the interval
%                  min(X) max(X) are generated in such a way that their
%                  squared Mahalanobis distance from the centroids of each
%                  existing group is larger then the quantile 1-0.999 of
%                  the chi^2 distribution with p degrees of freedom. In
%                  order to generate these these units the maximum number
%                  of attempts is equal to 10000.
%                - If noiseunits is a strcture it may contain the following
%                  fields:
%                  number : scalar, or vector of length f. The sum of the
%                       elements of vector 'number' is equal to the total
%                       number of outliers which are simulated.
%                  alpha : scalar or vector of legth f containing the
%                       level(s) of simulated outliers. The default value
%                       of alpha is 0.001.
%                  maxiter : maximum number of trials to simulate outliers.
%                       The default value of maxiter is 10000.
%                  typeout : list of length f containing the type of
%                       outliers which must be simulated. Possible values
%                       for typeout are:
%                       * unif (or uniform), if the outliers must be
%                         generated using the uniform distribution;
%                       * norm (or normal), if the outliers must be
%                         generated using the normal distribution;
%                       * Chisquarez, if the outliers must be generated
%                         using the Chi2 distribution with z degrees of
%                         freedom;
%                       * Tz or tz, if the outliers must be generated using
%                         the Student T distribution with z degrees of
%                         freedom;
%                       * pointmass, if the outliers are concentrated on a
%                         particular point;
%                       * componentwise, if the outliers must have the same
%                         coordinates of the existing rows of matrix X apart
%                         from a single coordinate (which will be to the min
%                         or max in that particular dimension).
%       lambda : vector of length p containing inverse Box-Cox
%                transformation coefficients. The value false (default)
%                implies that no transformation is applied to any variable.
%
%  Output:
% 
%           X  : simulated dataset of size (n + noiseunits)-by-(p + noisevars).
%                Noise coordinates are provided in the last noisevars columns.
%           id : classification vector of length n + noiseunits. Negative
%                numbers represents the groups associated to the
%                contaminated units.
%
%           REMARK: If noiseunits outliers could not be generated a warning
%                   is produced. In this case matrix X and vector id will
%                   have less than n + noiseunits rows.
%
%   DETAILS
%
% To make a dataset more challenging for clustering, a user might want to
% simulate noise variables or outliers. Parameter 'noisevars' specifies the
% desired number of noise variables. If an interval 'interval' is specified,
% noise will be simulated from a Uniform distribution on the interval given
% by 'interval'. Otherwise, noise will be simulated uniformly between the
% smallest and largest coordinates of mean vectors. 'noiseunits' specifies the
% number of observations outside (1 - 'alpha') ellipsoidal contours for the
% weighted component distributions. Outliers are simulated on a hypercube
% specified by the interval 'interval'. A user can apply an inverse Box-Cox
% transformation providing a vector of coefficients 'lambda'. The value 1
% implies that no transformation is needed for the corresponding coordinate
%
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
    spmplot(X,id);

    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noiseunits',10);
    spmplot(X,id);

    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noiseunits',10,'interval','minmax');
    spmplot(X,id);

    out = MixSim(4,3,'BarOmega',0.1);
    noiseunits=struct;
    noiseunits.typeout='uniform';
    noiseunits.number=10;
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noiseunits',noiseunits);
    spmplot(X,id);

    out = MixSim(2,2,'BarOmega',0.1);
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noiseunits',10,'noisevars',1);
    spmplot(X,id);

    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',1);
    spmplot(X,id);

    noiseunits=struct;
    noiseunits.number=10;
    noisevars=struct;
    noisevars.number=[1 1];
    noisevars.distribution={'chisquare5' 'T3'};
    noisevars.interval='minmax';
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noiseunits',noiseunits,'noisevars',noisevars);
    spmplot(X,id);
%}

%{
    % Check if it is possible to generate the outliers given alpha in the
    % interval specified by input option interval and if it not possible modify
    % interval by adding 0.1 until the outliers can be found
    out = MixSim(4,3,'BarOmega',0.1);
    % Point mass contamination of 30 observations in interval [0.4 0.4]
    nout=30;
    outint=[0.4 0.4];
    n=200;
    for j=1:10
        [X,id]=simdataset(n, out.Pi, out.Mu, out.S, 'noiseunits', nout, 'alpha', 0.01, 'interval', outint);
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

noisevarsdef    = '';
noiseunitsdef   = '';
alphadef        = 0.001;
intdef          = '';
maxiterdef      = 1e+05;
lambdadef       = '';
Rseeddef        = 0;

options=struct('noisevars',noisevarsdef,'noiseunits',noiseunitsdef,'alpha',alphadef,...
    'interval',intdef,'maxiter',maxiterdef,'lambda',lambdadef,'R_seed', Rseeddef);

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

R_seed     = options.R_seed;
noisevars  = options.noisevars;
noiseunits = options.noiseunits;
interval   = options.interval;
lambda     = options.lambda;

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
   
if (isstruct(noiseunits) || ~isempty(noiseunits)) || (isempty(noiseunits) && ~isempty(noisevars))
   
    if isstruct(noiseunits)  % Case where noiseunits is a structure
        fnoiseunits=fieldnames(noiseunits);

        d=find(strcmp('number',fnoiseunits));
        if d>0
            number=noiseunits.number;
            if (min(number) < 0)
                error('FSDA:simdataset:Wrongnnoise','Wrong value of outliers: it cannot be smaller than 0')
            end
        else
            number=1;
        end
        
        d=find(strcmp('typeout',fnoiseunits));
        if d>0
            typeout={noiseunits.typeout};
        else
            typeout={'uniform'};
        end
        
        d=find(strcmp('alpha',fnoiseunits));
        if d>0
            alpha=noiseunits.alpha;
        else
            alpha=0.001*ones(length(number),1);
        end
        
        d=find(strcmp('maxiter',fnoiseunits));
        if d>0
            maxiter=noiseunits.maxiter;
        else
            maxiter=10000;
        end
        
        if (number < 0)
            error('FSDA:simdataset:Wrongnout','Wrong value of number: it cannot be smaller than 0')
        end
        
        if ((max(alpha) >= 1) || (min(alpha) <= 0))
            error('FSDA:simdataset:WrongAlpha','Wrong value of alpha: it must be in the interval (0 1)')
        end
        
        if (maxiter < 1)
            error('FSDA:simdataset:WrongMaxIter','Wrong value for maximum number of iterations: it cannot be <1')
        end
        
        % noiseunits = total number of outliers  which has to be
        % simulated
        noiseunits=sum(number);
        
    else % Case where noiseunits is a scalar different from missing
        
        % noiseunits may be left empty when noisevars is specified, 
        % so this instruction is to generate one outlier by default        
        if isempty(noiseunits)
            noiseunits = 1;
        end
        
        number=noiseunits;
        %noiseunits=number; 
        typeout={'uniform'};
        alpha=0.001;
        maxiter=10000;
    end
    
    Xout=zeros(noiseunits, p);
    ni=0;
    idtmp=zeros(sum(number),1);
    for ii=1:length(number);
        typeouti=typeout{ii};
        
        [Xouti, ~] = getOutliers(number(ii), Mu, S, alpha(ii), maxiter, typeouti);
        Xout(ni+1:ni+size(Xouti,1),:)=Xouti;
        idtmp(ni+1:ni+size(Xouti,1))=-ii;
        ni=ni+size(Xouti,1);
    end
    
    if ni<noiseunits
        warning('FSDA:simdataset:Modifiedn',['Output matrix X will have just ' num2str(n+ni) ...
            ' rows and not ' num2str(n+noiseunits)])
    else
        X =[X;Xout];
    end
    id =[id;idtmp(1:ni)];
end

if isstruct(noisevars) || ~isempty(noisevars)
    if isstruct(noisevars)
        fnoisevars=fieldnames(noisevars);
        % labeladd option
        d=find(strcmp('number',fnoisevars));
        if d>0
            number=noisevars.number;
            if (min(number) < 0)
                error('FSDA:simdataset:Wrongnnoise','Wrong value of number of noisevars: it cannot be smaller than 0')
            end
        else
            number=1;
        end
        
        d=find(strcmp('distribution',fnoisevars));
        if d>0
            distribution=noisevars.distribution;
        else
            distribution='uniform';
        end
        
        d=find(strcmp('interval',fnoisevars));
        if d>0
            interval=noisevars.interval;
        else
            interval='';
        end
        
        % nvars = total number of noise variables to be simulated
        nvars=sum(number);
        
        %  if noisevars ~= 0
        if isempty(interval)
            L = min(min(Mu));
            U = max(max(Mu));
            L = L* ones(1,nvars);
            U = U* ones(1,nvars);
        elseif strcmp('minmax',interval)
            L = min(X);
            U = max(X);
        else
            L = interval(1,:);
            U = interval(2,:);
        end
    else % in this case noisevars is a scalar different from missing
        L = min(min(Mu));
        U = max(max(Mu));
        number=noisevars;
        nvars=number;
        distribution={'uniform'};
    end
    
    if R_seed
        % equivalent of 'rand(k,p)' in R is 'matrix(runif(k*p),k,p)'
        rn1s = ['matrix(runif(' num2str((n + noiseunits)*nvars) '),' num2str(n + noiseunits) ',' num2str(nvars) ')'];
        rrr = evalR(rn1s);
    else
        
        rrr=zeros(n + noiseunits, nvars);
        for ii=1:length(number);
            distributioni=distribution{ii};
            
            if strcmp(distributioni,'uniform')
                rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rand(n + noiseunits,  number(ii));
            elseif strcmp(distributioni,'norm') || strcmp(distributioni,'normal')
                % data generated from the normal distribution rescaled in the
                % interval [0 1]
                rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescale(randn(n + noiseunits,  number(ii)));
            elseif strcmp(distributioni(1),'T') || strcmp(distributioni(1),'t')
                nu=str2double(distributioni(2:end));
                rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescale(trnd(nu, n + noiseunits,  number(ii)));
            elseif strcmp(distributioni(1),'Chisquare')
                nu=str2double(distributioni(10:end));
                rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescale(chi2rnd(nu, n + noiseunits,  number(ii)));
            end
        end
    end
    
    % Values of noise variables were constrained to lie in the interval [0 1]
    % Now we rescale them to the interval [L U]
    Xnoise=bsxfun(@times,rrr,U-L);
    Xnoise=bsxfun(@plus,Xnoise,L);
    X = [X, Xnoise];  
end

if ~isempty(lambda)
    if (length(lambda) == p + noisevars)
        for j=1:(p + noisevars)
            X(:,j) = (lambda(j) * X(:, j) + 1).^(1/lambda(j)) - 1;
            if (sum(isnan(X(:,j))) ~= 0)
                warning('FSDA:simdataset:NaNs','NaNs were produced during transformation')
            end
        end
    else
        error('FSDA:simdataset:WrongLambda','The number of transformation coefficients lambda should be equal to ndimensions + noisevars')
    end
    % Remark: if lambda is very large MATLAB can create a complex X, so it
    % is necessary to have this further check
    X=real(X);
    
end
%% Nested functions
% Xout with nout rows which contains the outliers fail = scalar. If fail=1
% than it was not possible to generate the outliers in the interval
% specified by input option interval in maxiter trials else fail=0
    function [Xout,fail] = getOutliers(nout, Mu, S, alpha, maxiter, typeout)
        fail = 0;
        % maxiter = maximum number of iterations to generate outliers
        critval =chi2inv(1-alpha,p);
        
        Xout = zeros(nout,p);
        L = min(X);
        U = max(X);
        
        %         if isempty(interval)
        %             L = min(min(Mu));
        %             U = max(max(Mu));
        %         elseif strcmp(interval,'minmax')
        %             L = min(X);
        %             U = max(X);
        %
        %         else
        %             L = interval(1);
        %             U = interval(2);
        %         end
        
        i = 1;
        
        % Remark: maxiter1 must be much greater than nout
        if nout<2000
            maxiter1=20000;
        else
            maxiter1=nout*10;
            maxiter=maxiter1;
        end
        
        if strcmp(typeout(1:4),'unif')
            rrall = rand(maxiter1,p);
        elseif strcmp(typeout(1:4),'norm')
            rrall = rescale(randn(maxiter1,p));
        elseif strcmp(typeout(1),'T') || strcmp(typeout(1),'t')
            nuT=str2double(typeout(2:end));
            rrall = rescale(trnd(nuT, maxiter1, p));
        elseif strcmp(typeout(1:9),'Chisquare')
            nuC=str2double(typeout(10:end));
            rrall = rescale(chi2rnd(nuC,maxiter1,p));
        elseif   strcmp(typeout,'pointmass')
            rrall=rand(maxiter1,p);
        elseif   strcmp(typeout,'componentwise')
            % component wise contamination
        else
            error('FSDA:simdataset:WrongDistrib','Outlier distribution type not supported')
        end
        
        iter=0;
        while (i <= nout  &&  iter<maxiter)
            iter=iter+1;
            
            if strcmp(typeout,'componentwise')
                % extract one unit among those already extracted and
                % contaminate just a single random coordinate
                Xout(i,:)=X(randsample(1:n,1),:);
                contj=randsample(1:p,1);
                if rand(1,1)>0.5
                    Xout(i,contj)=U(contj);
                else
                    Xout(i,contj)=L(contj);
                end
            else
                
                if R_seed
                    % equivalent of 'rand(k,p)' in R is 'matrix(runif(k*p),k,p)'
                    rn1 = ['matrix(runif(' num2str(p) '),' num2str(1) ',' num2str(p) ')'];
                    rr = evalR(rn1);
                else
                    rindex=randsample(maxiter1,1);
                    rr=rrall(rindex,:);
                end
                
                Xout(i,:) = (U-L).*rr+L;
                
            end
            
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
            
            if ij==0 && strcmp(typeout,'pointmass')
                % Row of point mass contamination which has been found is
                % replicated nout times and stored inside matrix Xout
                % In other words, Xout will have nout equal rows
                Xalleq=repmat(Xout(1,:),nout,1);
                Xout=Xalleq;
                break
            end
        end
        
        % If iter = maxiter than it was not possible  to generate nout
        % outliers in maxiter simulations.
        if iter== maxiter
            disp(['Warning: it was not possible to generate ' num2str(nout) ' outliers'])
            disp(['in ' num2str(maxiter) ' replicates in the interval [' num2str(L(1)) ...
                '--' num2str(U(1)) ']'])
            disp('Please modify the interval inside input option ''interval'' ')
            disp('or increase input option ''alpha''')
            disp(['The values of interval and alpha now are ' num2str(interval) ' and ' num2str(alpha)]);
            
            % If max number of iteration has been reached fail is 1
            fail=1;
            Xout=Xout(1:i,:);
        end
    end
end