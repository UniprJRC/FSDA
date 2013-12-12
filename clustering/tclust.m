function [out , varargout]  = tclust(Y,k,alpha,restrfactor,varargin)
%   FEEDBACK GIANLUCA TODO>
% Nel commento a riga 64, "plots : Scalar or structure" in realta' e' solo scalar.
% nocheck non e' implementato
% i due parametri nomes e msg non potrebbero ridursi a uno solo 
%

%tclust computes trimmed clustering
%
%<a href="matlab: docsearch('tclust')">Link to the help function</a>
%
%   tclust(Y, k, alpha) partitions the points in the n-by-v data matrix Y
%   into k clusters.  This partition minimizes the trimmed sum, over all
%   clusters, of the within-cluster sums of point-to-cluster-centroid
%   distances.  Rows of Y correspond to points, columns correspond to
%   variables. tclust returns inside structure out an n-by-1 vector idx
%   containing the cluster indices of each point.  By default, tclust uses
%   (squared), possibly constrained, Mahalanobis distances.
%
%
%  Required input arguments:
%
%            Y: Data matrix containining n observations on v variables
%               Rows of Y represent observations, and columns
%               represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values will
%               automatically be excluded from the computations.
%            k: scalar which specifies the number of groups
%        alpha: global trimming level. alpha is a scalar between 0 and 0.5
%               or an integer specifying the number of observations which have to
%               be trimmed. If alpha=0 tclust reduces to traditional model
%               based or mixture clustering (mclust): see Matlab function
%               gmdistribution.
%               More in detail, if 0< alpha <1 clustering is based on
%                h=fix(n*(1-alpha)) observations
%               Else if alpha is an integer greater than 1 clustering is
%               based on h=n-floor(alpha);
%    restrfact: positive scalar which constrains the allowed differences
%               among group scatters. Larger values imply larger differences of
%               group scatters. On the other hand a value of 1 specifies the
%               strongest restriction.
%
%  Optional input arguments:
%
%       nsamp : Number of subsamples which will be extracted to find the
%               partition. If nsamp=0 all subsets will be extracted.
%               They will be (n choose k).
%               Remark: if the number of all possible subset is <300 the
%               default is to extract all subsets, otherwise just 300
%    refsteps : scalar defining number of refining iterations in each
%               subsample (default = 15).
%     reftol  : scalar. Default value of tolerance for the refining steps
%               The default value is 1e-14;
%equalweights : a logical value specifying whether cluster weights
%               shall be considered in the concentration,
%               assignment steps and computation of the likelihood.
%               if equalweights = true we are (ideally) assuming equally sized groups by maximizing
%
%
%                 \[
%                 \sum_{j=1}^k \sum_{ x_i \in group_j } \log f(x_i; m_j , S_j)
%                 \]
%
%               else if equalweights = false (default) we allow for different group weights by maximizing
%
%                 \[
%                     \sum_{j=1}^k  \sum_{ x_i \in group_j }  \log \left[ \frac{n_j}{n}  f(x_i; m_j , S_j) \right]=
%                 \]
%                 \[
%                   = \sum_{j=1}^k n_j \log n_j/n + \sum_{j=1}^k \sum_{ x_i \in group_j} \log f(x_i; m_j , S_j)
%                 \]
%
%               Remark: \sum_{j=1}^k n_j \log n_j/n is the so called entropy
%               term
%       mixt  : scalar, which specifies whether mixture modelling or crisp
%               assignment approach to model based clustering must be used.
%               In the case of mixture modelling parameter mixt also
%               controls which is the criterior to find the untrimmed units
%               in each step of the maximization
%               If mixt >=1 mixture modelling is assumed else crisp
%               assignment.
%                In mixture modelling the likelihood is given by
%                \begin{equation}\label{mixlik}
%                \prod_{i=1}^n  \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j).
%                \end{equation}
%               while in crisp assignment the likelihood is given by
%               \begin{equation}\label{clalik}
%               \prod_{j=1}^k   \prod _{i\in R_j} \phi (y_i; \; \theta_j),
%               \end{equation}
%               where $R_j$ contains the indexes of the observations which
%               are assigned to group $j$,
%               Remark: if mixt>=1 previous parameter equalweights is
%               automatically set to 1
%               Parameter mixt also controls the criterion to select the units to trim
%               if mixt == 2 the h units are those which give the largest
%               contribution to the likelihood that is the h largest
%               values of
%                   \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j)
%                    i=1, 2, ..., n
%               elseif mixt==1 the criterior to select the h units is
%               exactly the same as the one which is used in crisp
%               assignment. That is: the n units are allocated to a cluster
%               cluster according to criterior
%                \max_{j=1, \ldots, k} \hat \pi'_j \phi (y_i; \; \hat \theta_j)
%               and then these n numbers are ordered and the units
%               associated with the largest h numbers are untrimmed.
%       plots : Scalar or structure.
%               If plots = 1, a plot with the classification is
%               shown on the screen.
%        msg  : scalar which controls whether to display or not messages
%               on the screen If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the estimator
%               else no message is displayed on the screen
%      nocheck: Scalar. If nocheck is equal to 1 no check is performed on
%               matrix Y.
%               As default nocheck=0.
%        nomes: Scalar. If nomes is equal to 1 no message about estimated
%               time to compute tclust is displayed, else if nomes is
%               equal to 0 (default), a message about estimated time is
%               displayed.
%      startv1: scalar. If startv is 1 than initial
%               centroids and and covariance matrices are based on (v+1)
%               observations randomly chosen, else each centroid is
%               initialized taking a random row of input data matrix and
%               covariance matrices are initialized with identity matrices.
%               Remark: in order to start with a routine which is in the
%               required parameter space, eigenvalue restrictions are
%               immediately applied. The defualt value of startv1 is 0.
%       Ysave : Scalar that is set to 1 to request that the input matrix Y
%               is saved into the output structure out. Default is 0, i.e.
%               no saving is done.
%
%
%       Remark: The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%
%
%  Output:
%
%  The output consists of a structure 'out' containing the following fields:
%
%            out.idx  : n-by-1 vector containing assignment of each unit to
%                       each of the k groups. Cluster names are integer
%                       numbers from 1 to k. 0 indicates trimmed
%                       observations.
%            out.muopt: k-by-v matrix containing cluster centroid locations.
%                       Robust estimate of final centroids of the groups.
%         out.sigmaopt: v-by-v-by-k array containing estimated constrained
%                       covariance for the k groups.
%              out.bs : k-by-1 vector containing the units forming initial
%                       subset associated with muopt.
%            out.post : n-by-k matrix containing posterior probabilities
%                       out.post(i,j) contains posterior probabilitiy of unit
%                       i from component (cluster) j. For the trimmed units
%                       posterior probabilities are 0
%            out.siz  : matrix of size k-by-3
%                       1st col = sequence from 0 to k
%                       2nd col = number of observations in each cluster
%                       3rd col = percentage of observations in each cluster
%                       Remark: 0 denotes unassigned units
%   out.equalweights  : logical. It is true if in the clustering procedure
%                       we (ideally) assumed equal cluster weights
%                       else it is false if we allowed for different
%                       cluster sizes
%               out.h : scalar. Number of observations that have determined the
%                       centroids (number of untrimmed units).
%             out.obj : scalar. Value of the objective function which is minimized
%                       (value of the best returned solution).
%                       If input option mixt >1 the likelihood which is
%                       maximized is a mixture likelihood as follows
%                       \prod_{i=1}^h  \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j).
%                       else the likelihood which is maximized is a classification likelihood of the the form
%                       \prod_{j=1}^k   \prod _{i\in R_j} \pi_j' \phi (y_i; \; \theta_j),
%                       where $R_j$ contains the indexes of the observations which are assigned to group $j$
%                       with the constraint that $\# \bigcup_{j=1}^k
%                       R_j=h$. In the classification likelihood is input
%                       option equalweights=0 then \pi_j'=1 j=1, ..., k
%       out.notconver : scalar. Number of subsets without convergence
%              out.Y  : original data matrix Y. The field is present if option
%                       Ysave is set to 1.
%            out.AIC  : AIC
%            out.BIC  : BIC
%
% See also tkmeans, estepFS.m
%
% References:
%
% Garcia-Escudero, L.A.; Gordaliza, A.; Matran, C. and Mayo-Iscar, A.
% (2008), "A General Trimming Approach to Robust Cluster Analysis". Annals
% of Statistics, Vol.36, 1324-1345. Technical Report available at
% www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf
%
%
% Copyright 2008-2014. FSDA toolbox
%
% DETAILS. This iterative algorithm initializes k clusters randomly and
% performs "concentration steps" in order to improve the current cluster
% assignment. The number of maximum concentration steps to be performed is
% given by input parameter refsteps. For approximately obtaining the global
% optimum, the system is initialized nsamp times and concentration steps
% are performed until convergence or refsteps is reached. When processing
% more complex data sets higher values of nsamp and refsteps have to be
% specified (obviously implying extra computation time). However, if more
% then 10% of the iterations do not converge, a warning message is issued,
% indicating that nsamp has to be increased.
%
%
%<a href="matlab: docsearch('tclust')">Link to the help function</a>
% Last modified 08-Dec-2013

% Examples:

%


%{
    % tclust using geyser data
    Y=load('geyser2.txt');
    out=tclust(Y,3,0.1,10000,'plots',1)
    out=tclust(Y,3,0.1,10,'nsamp',100,'refsteps',10,'plots',1)
    % k-means solution restrfactor=1
    out=tclust(Y,3,0.1,1,'nsamp',100,'refsteps',20,'plots',1)
%}

%{
    % M5data
    %  A bivariate data set obtained from three normal bivariate distributions
    %  with different scales and proportions 1:2:2. One of the components is very
    %  overlapped with another one. A 10% background noise is added uniformly
    %  distributed in a rectangle containing the three normal components and not
    %  very overlapped with the three mixture components. A precise description
    %  of the M5 data set can be found in García-Escudero et al. (2008).
    Y=load('M5data.txt');
    plot(Y(:,1),Y(:,2),'o')

    spmplot(Y(:,1:2),Y(:,3),[],'box')

    out=tclust(Y(:,1:2),3,0,1000,'nsamp',100,'plots',1)
    out=tclust(Y(:,1:2),3,0,10,'nsamp',100,'plots',1)
    out=tclust(Y(:,1:2),3,0.1,1,'nsamp',1000,'plots',1,'equalweights',1)
    out=tclust(Y(:,1:2),3,0.1,1000,'nsamp',100,'plots',1)

%}

%{
    % Trimmed k-means using structured noise
    % The data have been generated using the following R instructions
    %    set.seed (0)
    %    v <- runif (100, -2 * pi, 2 * pi)
    %    noise <- cbind (100 + 25 * sin (v), 10 + 5 * v)
    %
    %
    %    x <- rbind (
    %        rmvnorm (360, c (0.0,  0), matrix (c (1,  0,  0, 1), ncol = 2)),
    %        rmvnorm (540, c (5.0, 10), matrix (c (6, -2, -2, 6), ncol = 2)),
    %        noise)


    %
    Y=load('structurednoise.txt');
    out=tclust(Y(:,1:2),2,0.1,100,'plots',1)
    out=tclust(Y(:,1:2),5,0.15,1,'plots',1)

%}

%{
    % Trimmed k-means using mixture100 data
    % The data have been generated using the following R instructions
    %     set.seed (100)
    %     mixt <- rbind (rmvnorm (360, c (  0,  0), matrix (c (1,  0,  0,  1), ncol = 2)),
    %                rmvnorm (540, c (  5, 10), matrix (c (6, -2, -2,  6), ncol = 2)),
    %                rmvnorm (100, c (2.5,  5), matrix (c (50, 0,  0, 50), ncol = 2)))


    %
    Y=load('mixture100.txt');
    out=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1)
    out=tclust(Y(:,1:2),3,0.05,1,'refsteps',20,'plots',1)
%}

%{
    % Compare different algorithms
    Y=load('mixture100.txt');
    % Traditional Tclust
    out1=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1)
    % Tclust with mixture models (selection of untrimmed units according to
    % likelihood contributions
    out2=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1,'mixt',1)
    % Tclust with mixture models (selection of untrimmed units according to
    % densities weighted by estimates of the probability of the components)
    out3=tclust(Y(:,1:2),3,0.05,1000,'refsteps',20,'plots',1,'mixt',2)

%}

%{
    n1=100;
    n2=80;
    n3=50;
    n4=80;
    n5=70;
    v=5;
    Y1=randn(n1,v)+5;
    Y2=randn(n2,v)+3;
    Y3=rand(n3,v)-2;
    Y4=rand(n4,v)+2;
    Y5=rand(n5,v);

    group=ones(n1+n2+n3+n4+n5,1);
    group(n1+1:n1+n2)=2;
    group(n1+n2+1:n1+n2+n3)=3;
    group(n1+n2+n3+1:n1+n2+n3+n4)=4;
    group(n1+n2+n3+n4+1:n1+n2+n3+n4+n5)=5;


    Y=[Y1;Y2;Y3;Y4;Y5];
   out=tclust(Y,5,0.05,1.3,'refsteps',20,'plots',1)
   
%}


%% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);
[n, v]=size(Y);


%% User options
% Remark: startv1 must be immediately checked because the calculation of
% ncomb is immediately affected.
if nargin>4
    chkstartv1 = strcmp(varargin,'startv1');
    if sum(chkstartv1)>0
        startv1= cell2mat(varargin(find(chkstartv1)+1));
        % varargin(find(chkstartv1)+1);
        %     else
        %         startv1=1;
    else
        startv1=0;
    end
else
    startv1=0;
end

if startv1
    ncomb=bc(n,k*(v+1));
else
    % If the number of all possible subsets is <300 the default is to extract
    % all subsets otherwise just 300.
    % Notice that we use bc, a fast version of nchoosek. One may also use the
    % approximation floor(exp(gammaln(n+1)-gammaln(n-p+1)-gammaln(p+1))+0.5)
    ncomb=bc(n,k);
end
nsampdef=min(300,ncomb);
refstepsdef=15;
reftoldef=1e-5;

% Default
if nargin<3;
    alpha=0.05;
else
    if isempty(alpha)
        alpha=0.05;
    end
end

if nargin<4;
    restrfactor=12;
end


% Fix alpha equal to the trimming size
% h = number of observations which is used to compute the centroids

if alpha<0
    error('alpha must a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
end

if nargin<4
    restrfactor=12;
end


% h = number of untrimmed units
if alpha>=1
    h=n-floor(alpha);
else
    h=fix(n*(1-alpha));
end

options=struct('nsamp',nsampdef,'plots',0,'nocheck',0,'nomes',0,...
    'msg',1,'Ysave',0,'refsteps',refstepsdef,'equalweights',false,...
    'reftol',reftoldef,'mixt',0,'startv1',1);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('Error:: number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('Error:: in total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin > 4
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
    
    % And check if the optional user parameters are reasonable.
    
    % Check number of subsamples to extract
    if options.nsamp>ncomb;
        disp('Number of subsets to extract greater than (n k). It is set to (n k)');
        options.nsamp=0;
    elseif  options.nsamp<0;
        error('Number of subsets to extract must be 0 (all) or a positive number');
    end
    
    % Check restriction factor
    if restrfactor<1
        disp('Restriction factor smaller than 1. It is set to 1 (maximum contraint==>spherical groups)');
        restrfactor=1;
    else
    end
    
end

% Default values for the optional
% parameters are set inside structure 'options'

plots=options.plots;        % Plot of the resulting classification
nsamp=options.nsamp;        % Number of subsets to extract
equalweights=options.equalweights;    % Specify if assignment must take into account the size of the groups

refsteps=options.refsteps;
reftol=options.reftol;

%Initialize the objective function (trimmed variance) by a
%large  value
vopt=-1e+30;

msg=options.msg;            % Scalar which controls the messages displayed on the screen

nomes=options.nomes;        % if options.nomes==1 no message about estimated time to compute tclust is displayed

mixt=options.mixt;         % if options.mixt==1 mixture model is assumed

if mixt>=1 && equalweights == 1
    warning('options equalweights must be different from 1 if mixture model approach is assumed')
    warning('options equalweights is reset to 0')
end

%% Combinatorial part to extract the subsamples


if startv1
    [C,nselected] = subsets(nsamp,n,k*(v+1),ncomb,msg);
else
    [C,nselected] = subsets(nsamp,n,k,ncomb,msg);
    niinistart=repmat(floor(h/k),k,1);
end

% Store the indices in varargout
if nargout==2
    varargout={C};
end

% ll = matrix of loglikelihoods for each unit from each cluster
% rows of ll are associated to units
% Columns of ll are associated to clusters
ll=zeros(n,k);
obj=1e+14;


if mixt>=1
    % log_lh = h-by-k matrix containing the log of component conditional
    %               density weighted by the component probability.
    %               log_lh = log( \pi_j \phi (y_i; \; \theta_j))
    log_lh=zeros(h,k);
end

% noconv = scalar linked to the number of times in which there was no
% convergence
noconv=0;

% The covariances are given initially by k identity matrices
ey=eye(v,v);
eyk=repmat(ey,[1 1 k]);

onev1=ones(v,1);

% Lambda_vk = matrix which will contain in column j the v (unrestricted)
% eigevalues of covariance matrix of group j (j=1, ..., k)
Lambda_vk=ones(v,k);

% sigmaopt = 3 dimensional array which will contain the estimates of the
% covariance matrices for the best solution
sigmaopt=zeros(v,v,k);

% Initialise and start timer.
tsampling = ceil(min(nselected/5 , 10));
time=zeros(tsampling,1);


% Lambda will contain the matrix of eigenvalues in each iteration for
% all groups. Lambda is a 3D array of size v-by-v-by-k
% Lambda=sigmaini;
% U will contain the eigenvectors of the cov matrices in each iteration
% for all groups. U is a 3D array of size v-by-v-by-k
sigmaini=zeros(v,v,k);
U=sigmaini;

% fullsol = vector which stores value of the objective function in each
% iteration
fullsol=zeros(nselected,1);


%% Core of tclust function
for i=1:nselected
    if i <= tsampling
        tstart = tic;
    end
    
    
    
    if startv1
        
        % Initialize niini with with random numbers from uniform
        randk=rand(k,1);
        niini=floor(h*randk/sum(randk));
        
        cini=zeros(k,v);
        for j=1:k
            ilow=(j-1)*(v+1)+1;
            iup=j*(v+1);
            index=C(i,:);
            selj=index(ilow:iup);
            % cini(j,:)=mean(Y(selj,:));
            Yselj=Y(selj,:);
            cini(j,:)=sum(Yselj)/(v+1);
            
            Yseljc = bsxfun(@minus,Yselj,cini(j,:));
            sigmaini(:,:,j) = (Yseljc' * Yseljc) / (v+1);
            
            % sigmaini(:,:,j)=cov(Y(selj,:));
            
            % Eigenvalue eigenvector decomposition for group j
            [Uj,Lambdaj] = eig(sigmaini(:,:,j));
            % Store eigenvectors and eigenvalues of group j
            U(:,:,j)=Uj;
            Lambda_vk(:,j)=diag(Lambdaj);
        end
        
        Lambda_vk(Lambda_vk<0)=0;
        autovalues=restreigen(Lambda_vk,niini,restrfactor);
        
        % Covariance matrices are reconstructed keeping into account the
        % constraints of the eigenvalues
        for j=1:k
            %disp(j)
            sigmaini(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
            
            % Alternative code: in principle more efficient but slower
            % because diag is a built in function
            % sigmaini(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
            
        end
    else
        
        % initialization of niini with equal proportions
        niini=niinistart;
        
        % extract a subset of size v
        index = C(i,:);
        
        % cini will contain the centroids in each iteration
        cini=Y(index,:);
        % sigmaini will contain the covariance matrices in each iteration
        sigmaini=eyk;
    end
    
    % sigmaopt will be final estimate of the covariance matrices
    % sigmaopt=sigmaini;
    
    
    
    iter=0;
    mudiff=1e+15;
    
    postprob=0;
    ind=0;
    
    % refsteps "concentration" steps will be carried out
    while ( (mudiff > reftol) && (iter < refsteps) )
        iter = iter + 1;
        %disp(iter)
        if equalweights
            
            % In this case we are (ideally) assuming equally sized groups
            for j=1:k
                ll(:,j)= logmvnpdfFS(Y,cini(j,:),sigmaini(:,:,j));
            end
            
        else
            
            % In this case we allow for different group weights or we are
            % assuming a mixture model
            for j=1:k
                % REMARK: we use log(niini(j)) instead of log(niini(j)/h)
                % because h is constant
                ll(:,j)= log(niini(j)/h) +  logmvnpdfFS(Y,cini(j,:),sigmaini(:,:,j));
                % This is faster but equivalent to
                % ll(:,j)= (niini(j)/h)*mvnpdf(Y,cini(j,:),sigmaini(:,:,j));
            end
            
            
            
        end
        
        if mixt==2
            
            postprobold=postprob;
            
            [~,postprob,disc]=estepFS(ll);
            
            
            % Sort the n likelihood contributions
            % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
            [~,qq]=sort(disc,'descend');
            
            
            % qq = vector of size h which contains the indexes associated with the largest n(1-alpha)
            % (weighted) likelihood contributions
            qqunassigned=qq((h+1):n);
            qq=qq(1:h);
            
            % Ytri = n(1-alpha)-by-v matrix associated with the units
            % which have the largest n(1-alpha) likelihood contributions
            Ytri=Y(qq,:);
            
            
            postprob(qqunassigned,:)=0;
            
            
            % M-step update of niini
            % niini = numerator of component probabilities
            niini=(sum(postprob))';
            
        else
            indold=ind;
            
            % In this part we select the untrimmed units
            % They are those which have the n(1-alpha) largest values among the
            % maxima of each row of matrix ll
            % vector disc of length(n) contains the (weighted) contribution of
            % each unit to the log likelihood
            [disc,ind]= max(ll,[],2);
            
            % Sort the n likelihood contributions
            % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
            [~,qq]=sort(disc,'descend');
            
            
            % qq = vector of size h which contains the indexes associated with the largest n(1-alpha)
            % (weighted) likelihood contributions
            qqunassigned=qq(h+1:end);
            qq=qq(1:h);
            
            % Ytri = n(1-alpha)-by-v matrix associated with the units
            % which have the largest n(1-alpha) likelihood contributions
            Ytri=Y(qq,:);
            % Ytriind = grouping indicator vector (of size n(1-alpha))
            % associated to Ytri
            groupind=ind(qq);
            
            % ind is the identifier vector
            % trimmed units have a value of ind=0
            ind(qqunassigned)=0;
        end
        
        
        if mixt == 1
            %  expll=exp(ll(qq,:));
            %  sumll=sum(expll,2);
            %  postprob=bsxfun(@rdivide,expll,sumll);
            
            % E-step: computation of posterior probabilities for untrimmed
            % units. In the context of mixture models posterior
            % probabilities will be used to estimate new component
            % probabilities of the mixtures, new centroids and new
            % covariance matrices
            
            postprobold=postprob;
            % OLD             [~,postprob]=estepFS(ll(qq,:));
            [~,postprob]=estepFS(ll);
            
            postprob(qqunassigned,:)=0;
            
            % M-step update of niini
            % niini = numerator of component probabilities
            niini=(sum(postprob))';
            
        end
        
        % M-step: parameters are updated
        % Matrix cini contains estimates of the new k centroids
        % Array sigmaini contains estimates of the new covariance matrices
        
        for j=1:k
            
            if mixt>=1
                % Matrix cini is updated using weighted means. The weights
                % are given by the posterior probabilities
                % Note that Y is used instead of Ytri because posterior
                % probabilities for unassigned units are 0
                cini(j,:)= sum(bsxfun(@times, Y, postprob(:,j)),1)/niini(j);
                
                if niini(j)>0
                    Ytric = bsxfun(@minus,Y,cini(j,:));
                    
                    sqweights = postprob(:,j).^(1/2);
                    
                    % Ytric = [X(:,1).*sqweights   X(:,2).*sqweights ...   X(:,end).*sqweights]
                    Ytric = bsxfun(@times, Ytric, sqweights);
                    
                    sigmaini(:,:,j) = (Ytric' * Ytric) / niini(j);
                    
                    % Eigenvalue eigenvector decomposition for group j
                    [Uj,Lambdaj] = eig(sigmaini(:,:,j));
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_vk(:,j)=diag(Lambdaj);
                else
                    sigmaini(:,:,j)=ey;
                    U(:,:,j)=ey;
                    Lambda_vk(:,j)=onev1;
                end
                
            else  % crisp assignment
                % Boolean index of units forming group j
                groupj=groupind==j;
                
                % Size of group j
                niini(j)=sum(groupj);
                
                % Group j values
                Ytrij=Ytri(groupj,:);
                % Means of group j
                cini(j,:)=sum(Ytrij)/niini(j);
                
                % niini=sum(Ytri(:,v+1)==j);
                if niini(j)>0
                    % Covariance of group j
                    % sigmaini(:,:,j)=cov(Ytrij);
                    % cov would recompute the sample means; code below is more
                    % efficient
                    Ytrijc = bsxfun(@minus,Ytrij,cini(j,:));
                    sigmaini(:,:,j) = (Ytrijc' * Ytrijc) / niini(j);
                    
                    
                    % Eigenvalue eigenvector decomposition for group j
                    [Uj,Lambdaj] = eig(sigmaini(:,:,j));
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_vk(:,j)=diag(Lambdaj);
                else
                    sigmaini(:,:,j)=ey;
                    U(:,:,j)=ey;
                    Lambda_vk(:,j)=onev1;
                end
                
            end
            
        end
        %disp(['iteration' num2str(iter)])
        %disp(niini)
        
        
        % Lambda_vk is a v-by-k  matrix whose jth column contains the
        % unrestricted eigenvalues of cov. matrix of group j   j=1, ..., k
        % The row below is just to avoid numerical problems
        Lambda_vk(Lambda_vk<0)=0;
       
        autovalues=restreigen(Lambda_vk,niini,restrfactor);
        
        % Covariance matrices are reconstructed keeping into account the
        % constraints of the eigenvalues
        for j=1:k
            %disp(j)
            sigmaini(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
            
            % Alternative code: in principle more efficient but slower
            % because diag is a built in function
            % sigmaini(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
            
        end
        
        
        % Calculus of the objective function (E-step)
        % oldobj=obj;
        obj = 0;
        
        if mixt>=1
            
            % Likelihood for mixture modelling
            
            %   log_lh is a h-by-k matrix where k is the number of Gaussian components of the mixture
            %   log_lh is the log of component conditional density weighted by the component
            %   probability.  The probability of j-th component is niini(j)/h
            %   log_lh(i,j) is log (Pr(point i|component j) * Prob( component j))
            
            for j=1:k
                log_lh(:,j)=  log(niini(j)/h)+logmvnpdfFS(Ytri,cini(j,:),sigmaini(:,:,j));
            end
            
            
            % obj contains the value of the log likelihood for mixture models
            obj=estepFS(log_lh);
            
            
        else
            % Likelihood for  crisp clustering
            for j=1:k
                % disp(ni)
                if niini(j)>0
                    if equalweights
                        % we simply sum the log of the densities for the untrimmed
                        % units
                        obj=obj+ sum(logmvnpdfFS(Ytri(groupind==j,:),cini(j,:),sigmaini(:,:,j)));
                    else
                        % niini(j)*log(niini(j)/h) is the so called entropy term
                        % which allows for different group weights
                        obj=obj+ niini(j)*log(niini(j)/h)+sum(logmvnpdfFS(Ytri(groupind==j,:),cini(j,:),sigmaini(:,:,j)));
                    end
                end
            end
            
        end
        
        
        if mixt>0
            % if mixt >0 stopping criterion is referred to postprob
            mudiff=sum(sum(abs(postprob-postprobold)))/n;
            % disp(mudiff)
        else
            % if mixt=0 stopping criterior is referred to no modiification in the classification
            mudiff=sum(abs(indold-ind)>0)/n;
            % disp(mudiff)
        end
        
        %disp(num2str(iter))
        %disp(mudiff)
        
        % Alternative stopping criterion was based  on the relative
        % modification of the objective function.
        %                  mudiff =abs(oldobj-obj)/abs(obj);
        %                  disp(['Iteration ' num2str(iter)])
        %                  disp([oldobj-obj]/abs(obj))
        %                  disp('monit')
        
        if iter==refsteps;
            noconv=noconv+1;
        end
        
    end
    
    % Store value of the objective function for iteration i 
    fullsol(i)=obj;
    
    % Store the centroids and the value of the objective function
    if obj>=vopt,
        % vopt = value of the objective function in correspondence of the
        % best centroids
        vopt=obj;
        % muopt = matrix containing best centroids
        muopt=cini;
        % nopt = vector containing sizes of the groups
        nopt=niini;
        % format long;
        %disp(index)
        %disp(obj)
        
        % sigmaopt
        sigmaopt=sigmaini;
        
        % store the indexes of the subset which gave rise to the
        % optimal solution
        bs=index;
    end
    
    
    if ~nomes
        if i <= tsampling
            % sampling time until step tsampling
            time(i)=toc(tstart);
        elseif i==tsampling+1
            % stop sampling and print the estimated time
            if msg==1
                fprintf('Total estimated time to complete tclust: %5.2f seconds \n', nselected*median(time));
            end
        end
    end
    
end
notconver=noconv/nselected;
if notconver>0.1;
    disp('------------------------------')
    disp(['Warning: Number of subsets without convergence equal to ' num2str(100*notconver) '%'])
end

%% Store quantities in out structure
%exist('muopt')==0
% Store robust estimate of final centroids of the groups
out.muopt=muopt;

% Store robust estimate of final covariance matrix of the groups
out.sigmaopt=sigmaopt;

% Store units forming initial subset which gave rise to the optimal
% solution
out.bs=bs;


% With the best obtained values for the parameters, we compute the final
% assignments and parameters


% construct the  log of component conditional
% density weighted by the component probability.
% ll = log( \pi_j \phi (y_i; \; \theta_j))
% Get the likelihood for each point with each component
% ll is a n by k matrix,
% if equalweights is false
% ll(i,j) is log( (n_j/h) * f(x_i|\theta_j))
% else if  equalweights is true
% ll(i,j) is log( f(x_i|\theta_j))
% f(x_i|\theta_j) is multivariate normal with theta_j =(mu_j, \Sigma_j)

if equalweights
    for j=1:k
        ll(:,j)= logmvnpdfFS(Y,muopt(j,:),sigmaopt(:,:,j));
    end
else
    for j=1:k
        ll(:,j)= log(nopt(j)/h) + logmvnpdfFS(Y,muopt(j,:),sigmaopt(:,:,j));
    end
end


% postprob n x k containing posterior probabilities
% logpdf n x 1 vector containg the n contributions to the log
% likelihood of mixture models
[~,postprob,logpdf]=estepFS(ll);

% Find final trimmed and untrimmed units for final classification
if mixt==2
    
    
    
    % Sort the n likelihood contributions
    % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
    [~,qq]=sort(logpdf,'descend');
    
    unassigned=qq(h+1:n);
    assigned=qq(1:h);
    
    % Store in vector idx the cluster associated to the highest posterior
    % probability
    [~,idx]=max(postprob,[],2);
    idx(unassigned)=0;
    
    postprob(unassigned,:)=0;
    % Remark:
    % If there was full convergence sum(logpdf(assigned)) = vopt
else
    
    % In this part we select the untrimmed units
    % They are those which have the n(1-alpha) largest values among the
    % maxima of each row of matrix ll
    % vector disc of length(n) contains the (weighted) contribution of
    % each unit to the log likelihood
    % idx = n x 1 vector containing the final assignments
    % disc = n x 1 vector which contains the likelihood of each unit to
    % the closest cluster
    
    [disc,idx]= max(ll,[],2);
    
    % Sort the n likelihood contributions
    % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
    [~,qq]=sort(disc,'descend');
    
    unassigned=qq(h+1:n);
    assigned=qq(1:h);
    % Assign observations to clusters and assign a 0 value to trimmed ones
    idx(unassigned)=0;
    postprob(unassigned,:)=0;
end


% Compute AIC and BIC

if mixt>=1
    % Compute value of the maximized log likelihood
    [NlogLmixt]=estepFS(ll(assigned,:));
    
    
    % Note that is there was convergence NlogL should be exactly equal to
    % -vopt
    NlogL = -NlogLmixt;
else
    
    % Note that disc(qq(1:h)) is the contribution to the loglikelihood
    % of the untrimmed units
    loglik=disc(qq(1:h));
    
    % NlogL is the negative of the log-likelihood of the mixture of the
    % untrimmed units
    % NlogL=-sum(max(ll,[],2));
    % Note that if there was convergence NlogL should be exactly equal to
    % -vopt
    NlogL =-sum(loglik);
end


% Store the assignments in matrix out
% Unassigned units have an assignment equal to 0

out.idx=idx;


% siz = matrix of size k x 3,
% 1st col = sequence from 0 to k
% 2nd col = number of observations in each cluster
% 3rd col = percentage of observations in each cluster
% sum(out.siz(:,2))=n
% sum(out.siz(:,3))=100
siz=tabulate(out.idx);
out.siz=siz;


% Number of estimated parameters
% k centroids of size v
% 0.5*v*(v+1) estimates for each of the k covariance matrices
npar=v*k; % +0.5*v*(v+1)*k;

% if equalweights = 0 the k-1 mixture proportions parameters must be added
if equalweights==0
    npar=npar +(k-1);
end



% Find number of constraints which are imposed on the covariance matrices
covunrestr=zeros(v,v,k);
% Find number of restricted eigenvalues for each group
constr=zeros(k,1);
for j=1:k
    
    selj=idx==j;
    nopt(j) = sum(selj);
    
    if nopt(j)>v
        covj=cov(Y(selj,:));
        covunrestr(:,:,j)=covj;
        % disp(['group' num2str(j) num2str(covj)])
        %         try
        [~,values] = eigs(covj);
        %         catch
        %             jj=100;
        %         end
        
        if v>1;
            values=sort(diag(values));
            eigun=values./values(1);
            % maximum value of r is v-1
            r=sum(eigun(2:end)>=restrfactor);
            
            constr(j)=r;
            
        else
            constr(j)=0;
        end
        
    else
    end
    
    %     eigunrestr((v*(j-1)+1):j*v)=diag(values);
    %        [~,values] = eigs(sigmaopt(:,:,j));
    %    eigrestr((v*(j-1)+1):j*v)=diag(values);
    %
end



%% Compute AIC and BIC

% disp(nconstr)

nParam=npar+(0.5*v*(v+1)*k-1)*( (1-1/restrfactor)^(k*v-1) )+1;
% nParam=npar+0.5*v*(v+1)*k;

% NlogL = - maximized log likelihood (for untrimmed observations)
BIC = 2*NlogL + nParam*log(h); % Note log(h) instead of log(n)  h=untrimmed units
AIC = 2*NlogL + 2*nParam;


nParamOld=npar+0.5*v*(v+1)*k;
BICold = 2*NlogL + nParamOld*log(h); % Note log(h) instead of log(n)  h=untrimmed units
AICold = 2*NlogL + 2*nParamOld;


%
%
% if equalweights == 1
%     siz=tabulate(Ytri(:,end));
%
%     D=bsxfun(@minus,D,log(siz(:,1)'));
% end


% Store the fraction of subsamples without convergence.
out.notconver=notconver;
alp=alpha>0;
if size(siz,1)<k+alp;
    warning('The total number of estimated clusters is smaller than the number supplied')
    warning(['Number of supplied clusters =' num2str(k)])
    warning(['Number of estimated clusters =' num2str(size(siz,1)-alp)])
end


if min(out.siz((1+alp):end,2)< n*0.02)
    warning('FSDA:tclust','Clusters with size < n * 0.02 found - try reducing k')
end


% Store value of the objective function (maximized trimmed log likelihood)
out.obj=vopt;

if out.obj==-1e+25;
    warning('FSDA:tclust','The result is artificially constrained due to restr.fact = 1')
end


out.equalweights=equalweights;

% Store the number of observations that have not been trimmed in the
% computation of the centroids
out.h=h;

% Store n x k matrix containing posterior probability
% of each row from each component (cluster)
out.post=postprob;

out.BIC=BIC;
out.AIC=AIC;
out.BICold=BICold;
out.AICold=AICold;

out.fullsol=fullsol;

if options.Ysave
    % Store original data matrix
    out.Y=Y;
end


%% Create plots
% Plot the groups in the scatter plot matrix
if plots==1;
    if v==1
        histFS(Y,length(Y),idx)
    elseif v==2
        colors = 'brcmykgbrcmykgbrcmykg';
        figure
        hold('on')
        for j=1:k
            idxj=idx==j;
            if sum(idxj)>0
                plot(Y(idxj,1),Y(idxj,2),'o','color',colors(j));
                ellipse(muopt(j,:),sigmaopt(:,:,j))
            end
        end
        if alpha>0
            idxj=idx==0;
            plot(Y(idxj,1),Y(idxj,2),'x','color','k');
        end
        
        axis equal
        iidx=unique(idx);
        hall = findobj(gca, 'type', 'line');
        hellipses = findobj(gca, 'type', 'line','Marker','none');
        hpoints = setdiff(hall,hellipses);
        clickableMultiLegend(hpoints,cellstr(num2str(iidx)));
        axis manual;
    else
        id=cellstr(num2str(idx));
        id(idx==0)=cellstr('Trimmed units');
        spmplot(Y,id);
    end
end



end






