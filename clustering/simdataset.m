function [X,id]=simdataset(n, Pi, Mu, S, varargin)
%simdataset simulates a dataset given the parameters of finite mixture model with Gaussian components
%
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
%   0.25, 0.25) function Nk1=mnrnd( n-k, Pi) is used to generate k integers
%   (whose sum is n-k) from the multinominal distribution with parameters
%   n-k and Pi. The size of the groups is given by Nk1+1. The first
%   Nk1(1)+1  observations are generated using centroid Mu(1,:) and
%   covariance S(:,:,1), ..., the last Nk1(k)+1 observations are generated
%   using centroid Mu(k,:) and covariance S(:,:,k)
%
%   DETAILS
%
% To make a dataset more challenging for clustering, a user might want to
% simulate noise variables or outliers. The optional parameter 'noiseunits'
% controls the number and the type of outliers which must be added. The
% optional parameter 'noisevars' controls the number and the type of noise
% variables which must be added (it is possible to control the
% distribution, the interval and the number). Finally, the user can apply
% an inverse Box-Cox transformation providing a vector of coefficients
% 'lambda'. The value 1 implies that no transformation is needed for the
% corresponding coordinate.
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
%                   noisevars.distribution={'Chisquare5' 'T3'};
%                   noisevars.interval='minmax';
%                adds 5 noise varibles, the first 3 generated using
%                the Chi2 with 5 degrees of freedom and the last two
%                using the Student t with 3 degrees of freedom. Noise
%                variables are generated in the interval min(X) max(X).
%       lambda : vector of length v containing inverse Box-Cox
%                transformation coefficients. The value false (default)
%                implies that no transformation is applied to any variable.
%      R_seed : scalar > 0 for the seed to be used to generate random numbers
%               in a R instance. This is used to check consistency of the
%               results obtained with the R package MixSim. See file
%               Connect_Matlab_with_R_HELP.m to know how to connect MATLAB
%               with R.  This option requires the installation of the
%               R-(D)COM Interface. Default is 0, i.e. random numbers are
%               generated by matlab.
%
%
%  Output:
%
%           X  : simulated dataset of size (n + noiseunits)-by-(v + noisevars).
%                Noise coordinates are provided in the last noisevars columns.
%           id : classification vector of length n + noiseunits. Negative
%                numbers represents the groups associated to the
%                contaminated units.
%
%           REMARK: If noiseunits outliers could not be generated a warning
%                   is produced. In this case matrix X and vector id will
%                   have less than n + noiseunits rows.
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

    %  Simulate dataset with 10 outliers
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noiseunits',10);

    %  Simulate dataset with 100 outliers
    out = MixSim(4,3,'BarOmega',0.1);
    n=300;
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noiseunits',100);
    spmplot(X,id);

%}


%{
    rng('default')
    % Generate 4 groups in 2 dimensions
    rng(100)
    out = MixSim(4,2,'BarOmega',0.01);
    n=300;
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S);
    spmplot(X,id);
    title('4 groups without noise and outliers')
%}

%{
    %% Add outliers generated from uniform distribution
    n=300;
    noisevars=0;
    noiseunits=3000;
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with outliers from uniform')
%}

%{
    %% Add outliers generated from Chi2 with 5 degrees of freedom
    n=300;
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=3000;
    % Add asymmetric very concentrated noise
    noiseunits.typeout={'Chisquare5'};
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with outliers from $\chi^2_5$','Interpreter','Latex')
%}

%{
    %% Add outliers generated from Chi2 with 40 degrees of freedom
    n=300;
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=3000;
    % Add asymmetric concentrated noise
    noiseunits.typeout={'Chisquare40'};
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with outliers from $\chi^2_{40}$','Interpreter','Latex')
%}

%{
    %% Add outliers generated from normal distribution
    n=300;
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=3000;
    % Add normal noise
    noiseunits.typeout={'normal'};
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with outliers from normal distribution','Interpreter','Latex')
%}

%{
    %% Add outliers generated from Student T with 5 degrees of freedom
    n=300;
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=3000;
    % Add outliers from T5
    noiseunits.typeout={'T5'};
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with outliers from Student T with 5 degrees if freedom','Interpreter','Latex')
%}

%{
    %% Add componentwise contamination
    n=300;
    noisevars='';
    noiseunits=struct;
    noiseunits.number=3000;
    % Add asymmetric concentrated noise
    noiseunits.typeout={'componentwise'};
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with component wise outliers','Interpreter','Latex')
%}

%{
    %% Add outliers generated from Chisquare and T distribution
    n=300;
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=5000*ones(2,1);
    noiseunits.typeout={'Chisquare3','T20'};
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with outliers from $\chi^2_{3}$ and $T_{20}$','Interpreter','Latex')
%}

%{
    %% Add outliers from Chisquare and T distribution and use a personalized value of alpha
    n=300;
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=5000*ones(2,1);
    noiseunits.typeout={'Chisquare3','T20'};
    noiseunits.alpha=0.2;
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with outliers from $\chi^2_{3}$ and $T_{20}$ and $\alpha=0.2$','Interpreter','Latex')
%}

%{
    %% Add outliers from Chi2 + point mass contamination and add one noise variable
    noisevars=struct;
    noisevars.number=1;
    noiseunits=struct;
    noiseunits.number=[100 100];
    noiseunits.typeout={'pointmass' 'Chisquare5'};
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups with outliers from $\chi^2_{5}$ and point mass $+1$ noise var','Interpreter','Latex')
%}

%{
    %% Add 5 noise variables 
    n=300;
    noisevars=struct;
    noisevars.number=[2 3];
    noisevars.distribution={'Chisquare3','T20'};
    noiseunits='';
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id,[],'box');
    title('4 groups in 2 dims with 5 noise variables. First two from $\chi^2_{3}$ and last three from $T_{20}$','Interpreter','Latex')
%}

%{
    %% Add 3 noise variables 
    n=300;
    noisevars=struct;
    noisevars.number=[1 2];
    noisevars.distribution={'Chisquare3','T2'};
    noiseunits='';
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups in 2 dims with 3 noise variables. First from $\chi^2_{3}$ and last two from $T_{2}$','Interpreter','Latex')
%}

%{
    %% Add 3 noise variables and use a personalized interval
    n=300;
    noisevars=struct;
    noisevars.number=[1 2];
    noisevars.distribution={'Chisquare3','T20'};
    noisevars.int='minmax';
    noiseunits='';
    % In this example we supply min and max for each noise variable
    v1=sum(noisevars.number);
    noisevars.int=[3*ones(1,v1); 10*ones(1,v1)];
    [X,id]=simdataset(n, out.Pi, out.Mu, out.S,'noisevars',noisevars,'noiseunits',noiseunits);
    spmplot(X,id);
    title('4 groups in 2 dims with 3 noise variables with personalized interval','Interpreter','Latex')
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
lambdadef       = '';
Rseeddef        = 0;


options=struct('noisevars',noisevarsdef,'noiseunits',noiseunitsdef,...
        'lambda',lambdadef,'R_seed', Rseeddef);

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
lambda     = options.lambda;




[k,v]=size(Mu);

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

X=zeros(n,v);
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

if isstruct(noiseunits) || ~isempty(noiseunits)
    
    % Set all default values for outlier generation
    % One single outliers is generated from the uniform distribution.
    % Maximum 10000 iterations are used and this point must have a minimum
    % (squared) Mahalanobis distance from each centroid of existing groups greater
    % than 0.999 confidence level
    fnoiseunitsdef=struct;
    fnoiseunitsdef.number=1;
    fnoiseunitsdef.typeout={'uniform'};
    fnoiseunitsdef.alpha=0.001;
    fnoiseunitsdef.maxiter=10000;
    
    if isstruct(noiseunits)
        
        fnoiseunits=fieldnames(noiseunits);
        
        % Check if user options inside options.fnoiseunits are valid options
        chkoptions(fnoiseunitsdef,fnoiseunits)
        
        % labeladd option
        d=find(strcmp('number',fnoiseunits));
        if d>0
            number=noiseunits.number;
            if (min(number) < 0)
                error('FSDA:simdataset:Wrongnnoise','Wrong value of outliers: it cannot be smaller than 0')
            end
        else
            number=fnoiseunitsdef.number;
        end
        
        
        d=find(strcmp('typeout',fnoiseunits));
        if d>0
            typeout=noiseunits.typeout;
        else
            typeout=fnoiseunitsdef.typeout;
        end
        
        d=find(strcmp('alpha',fnoiseunits));
        if d>0
            alpha=noiseunits.alpha;
            % Just in case alpha is supplied as a scalar and length(number)
            % is greater than 1, then it is necessary to resize alpha to
            % make it have length equal to length(number)
            if isscalar(alpha) &&  length(number)>1
                alpha = alpha*ones(length(number),1);
            end
        else
            alpha = fnoiseunitsdef.alpha*ones(length(number),1);
        end
        
        d=find(strcmp('maxiter',fnoiseunits));
        if d>0
            maxiter=noiseunits.maxiter;
        else
            maxiter=fnoiseunitsdef.maxiter;
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
        
        % nout = total number of outliers  which has to be
        % simulated
        noiseunits=sum(number);
    else % in this case noisevars is a scalar different from missing
        number=noiseunits;
        noiseunits=number;
        typeout=fnoiseunitsdef.typeout;
        alpha=fnoiseunitsdef.alpha;
        maxiter=fnoiseunitsdef.maxiter;
    end
    
    Xout=zeros(noiseunits, v);
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
        noiseunits=ni;
        
    end
    X =[X;Xout(1:ni,:)];
    id =[id;idtmp(1:ni)];
else
    noiseunits=0;
end

if isstruct(noisevars) || ~isempty(noisevars)
    % Set all default values for noise variable generation
    % One single noise variable is generated from the uniform distribution.
    % The values of the noise variable range from min(min(Mu)) and max(max(Mu))
    noisevarsdef=struct;
    noisevarsdef.number=1;
    noisevarsdef.distribution={'uniform'};
    noisevarsdef.interval='';
    
    if isstruct(noisevars)
        fnoisevars=fieldnames(noisevars);
        
        % Check if user options inside options.fnoisevars are valid options
        chkoptions(noisevarsdef,fnoisevars)
        
        % number option
        d=find(strcmp('number',fnoisevars));
        if d>0
            number=noisevars.number;
            if (min(number) < 0)
                error('FSDA:simdataset:Wrongnnoise','Wrong value of number of noisevars: it cannot be smaller than 0')
            end
        else
            number=noisevarsdef.number;
        end
        
        % distribution
        d=find(strcmp('distribution',fnoisevars));
        if d>0
            distribution=noisevars.distribution;
        else
            distribution=noisevarsdef.distribution;
        end
        
        % interval
        d=find(strcmp('interval',fnoisevars));
        if d>0
            interval=noisevars.interval;
        else
            interval=noisevarsdef.interval;
        end
        
        % nvars = total number of noise variables which has to be
        % simulated
        nvars=sum(number);
        
        %  if noisevars ~= 0
        if isempty(interval)
            L = min(min(Mu));
            U = max(max(Mu));
            L = L* ones(1,nvars);
            U = U* ones(1,nvars);
        elseif strcmp('minmax',interval)
            L = min(min(X));
            U = max(max(X));
            L = L* ones(1,nvars);
            U = U* ones(1,nvars);
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
            
            if strcmp(distributioni(1),'T') || strcmp(distributioni(1),'t')
                nu=str2double(distributioni(2:end));
                rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescale(trnd(nu, n + noiseunits,  number(ii)));
            elseif strcmp(distributioni,'norm') || strcmp(distributioni,'normal')
                % data generated from the normal distribution rescaled in the
                % interval [0 1]
                rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescale(randn(n + noiseunits,  number(ii)));
            elseif strcmp(distributioni,'uniform')
                rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rand(n + noiseunits,  number(ii));
            elseif strcmp(distributioni(1:9),'Chisquare')
                nu=str2double(distributioni(10:end));
                rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescale(chi2rnd(nu, n + noiseunits,  number(ii)));
            else
                error('FSDA:simdataset:WrongDistrib','Variable distribution type not supported')
            end
        end
    end
    
    % Values of noise variables were constrained to lie in the interval [0 1]
    % Now we rescale them to the interval [L U]
    Xnoise=bsxfun(@times,rrr,U-L);
    Xnoise=bsxfun(@plus,Xnoise,L);
    X = [X, Xnoise];
else
    nvars=0;
end

if ~isempty(lambda)
    if (length(lambda) == v + nvars)
        for j=1:(v + nvars)
            X(:,j) = (lambda(j) * X(:, j) + 1).^(1/lambda(j)) - 1;
            if (sum(isnan(X(:,j))) ~= 0)
                warning('FSDA:simdataset:NaNs','NaNs were produced during transformation')
            end
        end
    else
        error('FSDA:simdataset:WrongLambda','The number of transformation coefficients lambda should be equal to v+ number of noise vars')
    end
    % Remark: if lambda is very large MATLAB can create a complex X, so it
    % is necessary to have this further check
    X=real(X);
    
end
%% Nested functions
% Xout with nout rows which contains the outliers
% fail = scalar. If fail =1 than it was not possible to generate the
% outliers in maxiter trials, using confidence interval based on alpha
% else fail =0
    function [Xout,fail] = getOutliers(nout, Mu, S, alpha, maxiter, typeout)
        fail = 0;
        % maxiter = maximum number of iterations to generate outliers
        critval =chi2inv(1-alpha,v);
        
        Xout = zeros(nout,v);
        L = min(X);
        U = max(X);
        
        i = 1;
        
        % Remark: maxiter1 must be much greater than nout
        if nout<2000
            maxiter1=20000;
        else
            maxiter1=nout*10;
            maxiter=maxiter1;
        end
        
        if strcmp(typeout(1),'T') || strcmp(typeout(1),'t')
            nuT=str2double(typeout(2:end));
            rrall = rescale(trnd(nuT, maxiter1, v));
        elseif strcmp(typeout(1:4),'norm')
            rrall = rescale(randn(maxiter1,v));
        elseif strcmp(typeout(1:4),'unif')
            rrall = rand(maxiter1,v);
        elseif strcmp(typeout(1:9),'Chisquare')
            nuC=str2double(typeout(10:end));
            rrall = rescale(chi2rnd(nuC,maxiter1,v));
        elseif   strcmp(typeout,'pointmass')
            rrall=rand(maxiter1,v);
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
                contj=randsample(1:v,1);
                if rand(1,1)>0.5
                    Xout(i,contj)=U(contj);
                else
                    Xout(i,contj)=L(contj);
                end
            else
                
                if R_seed
                    % equivalent of 'rand(k,p)' in R is 'matrix(runif(k*p),k,p)'
                    rn1 = ['matrix(runif(' num2str(v) '),' num2str(1) ',' num2str(v) ')'];
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
            disp(['Number of values which was possible to generate is equal to ' num2str(i)])
            disp('Please modify the type of outliers using option ''typeout'' ')
            disp('or increase input option ''alpha''')
            if isempty(interval)
                disp(['The values of int and alpha now are [0 1] and ' num2str(alpha)]);
            else
                disp(['The values of int and alpha now are ' num2str(interval) ' and ' num2str(alpha)]);
            end
            % If max number of iteration has been reached fail is 1
            fail=1;
            Xout=Xout(1:i,:);
        end
    end
end