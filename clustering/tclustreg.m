function [out] = tclustreg(y,X,k,restrfact,alpha1,alpha2,varargin)
%tclustreg performs robust linear grouping analysis
%
%<a href="matlab: docsearchFS('tclustreg')">Link to the help function</a>
%
%  Required input arguments:
%
%         y : Response variable. Vector.
%             A vector with n elements that contains the response variable.
%             y can be either a row or a column vector.
%             Data Types - single|double
%
%         X : Explanatory variables (also called 'regressors'). Matrix.
%             Data matrix of dimension $(n \times p-1)$. Rows of X represent
%             observations, and columns represent variables. Missing values
%             (NaN's) and infinite values (Inf's) are allowed, since
%             observations (rows) with missing or infinite values will
%             automatically be excluded from the computations.
%             Data Types - single|double
%
%         k : Number of clusters. Scalar.
%             This is a guess on the number of data groups.
%             Data Types - single|double
%
% restrfact : Scatter constraint. Scalar.
%            This is a constant c controlling the differences among
%            group scatters. The value 1 is the strongest restriction.
%            Data Types - single|double
%
%   alpha1 : Trimming level. Scalar.
%            alpha1 is a value between 0 and 0.5 or an  integer specifying
%            the number of observations which have to be trimmed. If
%            alpha=0 there is no trimming. More in detail, if 0<alpha1<1
%            clustering is based on h=fix(n*(1-alpha1)) observations.
%            Else if alpha1 is an integer greater than 1 clustering is
%            based on h=n-floor(alpha1).
%            Data Types - single|double
%
%   alpha2 : Second-level trimming. Scalar.
%            alpha2 is a value between 0 and 0.5, usually smaller than
%            alpha1. If alpha2=0 there is no second-level trimming.
%            Data Types - single|double
%
%
%  Optional input arguments:
%
%intercept : Indicator for constant term. Scalar. If 1, a model with
%            constant term will be fitted (default), if 0, no constant
%            term will be included.
%            Example - 'intercept',1
%            Data Types - double
%
%    niter : Number of random starts. Scalar. An integer for the number
%            of iterations to attempt for convergence.
%            Example - niter = 20
%            Data Types - double
%   mixt   : mixture modelling or crisp assignmen. Scalar.
%            Option mixt specifies whether mixture modelling or crisp
%            assignment has to be used:
%            mixt = 2 is for mixture modelling;
%            mixt = 0 is for crisp assignment.
%            In mixture modelling, the likelihood is given by INSERT FORMULA DOME
%            In crisp assignment, the  likelihood is given by INSERT FORMULA DOME
%            Example - 'mixt',0
%            Data Types - single | double
%    nsamp : number of subsamples to extract.
%            Scalar or matrix.
%            If nsamp is a scalar it contains the number of subsamples
%            which will be extracted.
%            If nsamp=0 all subsets will be extracted.
%            Remark - if the number of all possible subset is <300 the
%            default is to extract all subsets, otherwise just 300.
%            If nsamp is a matrix it contains in the rows the indexes of
%            the subsets which have to be extracted. nsamp in this case can
%            be conveniently generated  by function subsets.
%            nsamp can have k columns or k*(v+1) columns. If nsamp has k
%            columns the k initial regression parameters in each iteration
%            i are given by X(nsamp(i,:),:) and the variances are equal to
%            the identity.
%            If nsamp has k*(v+1) columns the initial centroids and
%            covariance matrices in iteration i are computed as follows
%               X1=X(nsamp(i,:),:)
%               mean(X1(1:v+1,:)) contains the initial centroid for group 1
%               cov(X1(1:v+1,:)) contains the initial cov matrix for group 1               1
%               mean(X1(v+2:2*v+2,:)) contains the initial centroid for group 2
%               cov((v+2:2*v+2,:)) contains the initial cov matrix for group 2               1
%               ...
%               mean(X1((k-1)*v+1:k*(v+1))) contains the initial centroids for group k
%               cov(X1((k-1)*v+1:k*(v+1))) contains the initial cov matrix for group k
%               REMARK - if nsamp is not a scalar option option below
%               startv1 is ignored. More precisely if nsamp has k columns
%               startv1=0 elseif nsamp has k*(v+1) columns option startv1=1.
%             Example - 'nsamp',1000
%             Data Types - double
%  startv1: how to initialize regression parameters. Scalar.
%           If startv1 is 1 then initial regression parameters are based on
%           (v+1) observations randomly chosen, else each regression is
%           initialized taking a random row of input data matrix. Remark 1-
%           in order to start with a routine which is in the required
%           parameter space, eigenvalue restrictions are immediately
%           applied. The default value of startv1 is 1. Remark 2 - option
%           startv1 is used just if nsamp is a scalar (see for more details
%           the help associated with nsamp)
%           Example - 'startv1',1
%           Data Types - single | double
% Ksteps:  Number of refining iterations. Scalar. Number of refining
%               iterations in each subsample.  Default is 10.
%               Ksteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'Ksteps',15
%                 Data Types - single | double
%    plots : Plot on the screen. Scalar. A flag to control the
%            generation of the plots.
%            If plots=1 a plot is showed on the screen with the
%            final allocation (and if size(X,2)==2 with the lines
%            associated to the groups)
%            Example - 'plots',1
%            Data Types - double
%   wtrim: Application of weights. Scalar. A flag taking values [0, 1, 2, 3, 4]
%          to control the application of weights on the observations.
%          -  If \texttt{wtrim}=0 (no weights) and \texttt{mixt}=0, the
%             algorithm reduces to the standard tclustreg algorithm.
%          -  If \texttt{wtrim}=0 and \texttt{mixt}=2, the maximum posterior
%             probability $D\_i$ of equation 7 of Garcia et al. 2010 is
%             computing by maximizing the log-likelihood contributions of
%             the mixture model of each observation.
%          -  If \texttt{wtrim} = 1, trimming is done by weighting the
%             observations using values specified in vector \texttt{we}.
%             In this case, vector \texttt{we} must be supplied by the
%             user. For instance, \texttt{we} = $X$.
%          -  If \texttt{wtrim} = 2, trimming is again done by weighting
%             the observations using values specified in vector \texttt{we}.
%             In this case, vector \texttt{we} is computed from the data as
%             a function of the density estimate $\mbox{pdfe}$.
%            Specifically, the weight of each observation is the
%            probability of retaining the observation, computed as
%            \[\mbox{pretain}_{i g} = 1 - \mbox{pdfe}_{ig}/\max_{ig}(\mbox{pdfe}_{ig})\]
%         -  If \texttt{wtrim} = 3, trimming is again done by weighting the
%            observations using values specified in vector \texttt{we}. In
%            this case, each element $we_i$ of vector \texttt{we} is a
%            Bernoulli random variable with probability of success
%            $\mbox{pdfe}_{ig}$. In the clustering framework this is done
%            under the constraint that no group is empty.
%         -  If \texttt{wtrim} = 4, trimming is done with the tandem approach
%            of Cerioli and Perrotta (2014).
%            Example - 'wtrim',1
%            Data Types - double
%      we: Vector of observation weights. Vector. A vector of size nX1
%          containing the weights to apply to each observation. Default
%          value: vector of ones.
%            Example - 'we',[0.2 0.2 0.2 0.2 0.2]
%            Data Types - double
%eps_beta: minimum accepted difference between regression coefficients in
%           the initial subsets. Scalar. If the observation in the initial subsets are
%           collinear, it can happen that the number of groups identified is less than
%           p. To avoide this behavior, eps_beta>0 allows to start the refining steps
%           of the tclust algorithm from subsets chosen in a better way.
%           Default value: 0.
%            Example - 'eps_beta',0.01
%            Data Types - double
%        msg  : Level of output to display. Scalar.
%               Scalar which controls whether to display or not messages
%               on the screen. If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the estimator
%               or the number of subsets in which there was no convergence
%               else no message is displayed on the screen
%                 Example - 'msg',1
%                 Data Types - single | double
%
%  Output:
%
%  out :  structure containing the following fields
%
%   out.bopt           = $p-1 \times k$ matrix containing the regression
%                        parameters.
%   out.sigmaopt0      = $k$ row vector containing the estimated group
%                        variances.
%   out.sigmaopt_cons  = $k$ row vector containing the estimated group
%                        variances corrected with  asymptotic consistency factor
%   out.sigmaopt_pison = $k$ row vector containing the estimated group
%                            variances corrected with  asymptotic consistency factor
%                            and small sample correction factor of Pison et al.
%   out.numopt         = $k$ column vector containing the number of
%                        observations in each cluster
%                        after the second trimming.                                         .
%   out.vopt           = Scalar. The value of the target function.
%   out.asig1          = $n$ vector containing the cluster assigments after
%                        first trimming ('0' means a trimmed observation).
%   out.asig2          = $n$ vector containing the final cluster assigments
%                        after second trimming ('0' means a trimmed
%                        observation).
%   out.postprob       = $n$ vector containing the final posterior probability
%   out.count1_ng_lt_k = number of times that, after the first level of trimming, in a group there are not enought observations to compute the sigma
%   out.count1_eq_lt_k = number of times that, after the first level of trimming, in a group there are enought observations to compute the sigma
%   out.count2_ng_lt_k = number of times that, after the second level of trimming, in a group there are not enought observations to compute the sigma
%   out.count2_eq_lt_k = number of times that, after the second level of trimming, in a group there are enought observations to compute the sigma
%   out.nselected      = number of initial subsets actually
%                        extracted. If eps_beta is not specified or if it is set to
%                        zero, out.nselected = nsamp; otherwise out.nselected > nsamp
%  out.selj_all        = initial subsets extracted
%
% See also: tclust, tkmeans, estepFS
%
% References:
%
% Garcia-Escudero, L.A.; Gordaliza, A.; Matran, C. and Mayo-Iscar, A.
% (2008), "A General Trimming Approach to Robust Cluster Analysis". Annals
% of Statistics, Vol.36, 1324-1345. Technical Report available at
% www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf
%
% Cerioli, A. and Perrotta, D. (2014). "Robust Clustering Around Regression
% Lines with High Densoty Regions". Advances in Data Analysis and
% Classification, Volume 8, Issue 1, p. 5-26.
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustreg')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
%
%
% Examples:
%
%{
%% tclustreg of X data using number of groups k=2, restriction factor 50, alpha1 = 0.01, alpha2 = 0.01.
    X   = load('X.txt');
    out = lga(X,3);

    y1=X(:,end);
    X1=X(:,1:end-1);

    k = 3 ; restrfact = 5; alpha1 = 0.1 ; alpha2 = 0.1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2);

    k = 2 ; restrfact = 10; alpha1 = 0.005 ; alpha2 = 0.001;
    we = abs(X1/sum(X1));
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',1,'we',we,'wtrim',1,'mixt',2,'plots',0);
%}
%{
    load fishery;
    X=fishery.data;
    % some jittering is necessary because duplicated units are not treated
    % in tclustreg: this needs to be addressed
    X = X + 10^(-8) * abs(randn(677,2));

    out=lga(X,3);
    clickableMultiLegend('1','2','3','data1','data2','data3');
    axis manual;

    alpha = 0.95;
    out=rlga(X,3,1-alpha);
    clickableMultiLegend('0','1','2','3','data1','data2','data3');
    axis manual;


    y1 = X(:,end);
    X1 = X(:,1:end-1);
    k = 3 ; restrfact = 50; alpha1 = 0.04 ; alpha2 = 0.01;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',2);
%}

%{
    XX = X1;
    XX = sqrt(X1);
    XX = X1.^(1/3);
    we = XX/sum(XX);
    
    mixt = 0; wtrim = 1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'we',we,'wtrim',wtrim);

    mixt = 2; wtrim = 1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'we',we,'wtrim',wtrim);


    mixt = 0; wtrim = 2; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

%}

%{
    %% Generate mixture of regression using MixSimReg, with an average
    % overlapping at centroids =0.01. Use all default options.
    p=3;
    k=2;
    Q=MixSimreg(k,p,'BarOmega',0.001);
    n=400;
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    out=tclustreg(y,X,2,50,0.01,0.01,'intercept',1);

%}
%{
    % Generate mixture of regression using MixSimReg, with an average
    % overlapping at centroids =0.01. Use all default options.
    p=3;
    k=2;
    Q=MixSimreg(k,p,'BarOmega',0.001);
    n=400;
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
     we=X(:,2)/sum(X(:,2));
    out=tclustreg(y,X,2,50,0.01,0.01,'intercept',1,'we',we,'wtrim',1,'mixt',2);

%}


%% Check if optimization toolbox is installed in current computer
% to be done in next releases: introduce an optimizer

typemin = exist('fminunc','file');
if typemin ~=2
    error('FSDA:tclustreg:MissingOptToolbox','This function requires the optimization toolbox');
end


%% initializations

warning('off');

% internal, used for debugging purposes
debug_mode=0;

% 'oversamp' is a factor (which depends on the number of groups 'k') used
% to generate more samples in order to face the possibility that some
% subsets contain collinear obs.
% To obtain nsamp = 300 samples, 300*oversamp samples will be generated
oversamp = 10*k;

%number of times that, after the first level of trimming, in a group there
%are not enought observations to compute the sigma
count1_ng_lt_k = 0;

%number of times that, after the second level of trimming, in a group there
%are not enought observations to compute the sigma
count2_ng_lt_k= 0;

%number of times that, after the first level of trimming, in a group there
%are enought observations to compute the sigma
count1_ng_eq_k = 0;

%number of times that, after the second level of trimming, in a group there
%are enought observations to compute the sigma
count2_ng_eq_k= 0;

% tolerance for restriction factor
tolrestreigen = 1e-08;

%Initialization for the objective function (optimized during the random
%starts) through a very small value
vopt = -1e+20;

% this is just for rotating colors in the plots
clrdef = 'bkmgyrcbkmgyrcbkmgyrcbkmgyrcbkmgyrcbkmgyrcbkmgyrc';

% repmat from Release 8.2 is faster than bsxfun
verMatlab = verLessThan('matlab','8.2.0');
if verMatlab ==1
    userepmat=0;
else
    userepmat=1;
end

%% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% check restrfact option
if nargin < 4 || isempty(restrfact) || ~isnumeric(restrfact)
    restrfact = 12;
end

% checks on alpha1 and alpha2
if alpha1>0.5
    error('FSDA:tclustreg:WrongAlpha1','Error:: alpha1 must be in [0, 0.5]');
end

if alpha2>0.5
    error('FSDA:tclustreg:WrongAlpha2','Error:: alpha2 must be in [0, 0.5]');
end

% startv1def = default value of startv1 = 1
% initialization using covariance matrices based on v+1 units
startv1def = 1;

if nargin>6
    % Check whether option nsamp exists
    chknsamp = strcmp(varargin,'nsamp');
    
    % if the sum below is greater than 0 option nsamp exists
    if sum(chknsamp)>0
        nsamp=cell2mat(varargin(find(chknsamp)+1));
        
        % Check if options nsamp is a scalar
        if ~isscalar(nsamp)
            % if nsamp is not a scalar, it is a matrix which contains in
            % the rows the indexes of the subsets which have to be
            % extracted
            C=nsamp;
            [nsampdef,ncolC]=size(C);
            % The number of rows of nsamp (matrix C) is the number of
            % subsets which have to be extracted
            nselected=nsampdef;
            % If the number of columns of nsamp (matrix C) is equal to v
            % then the procedure is initialized using identity matrices
            % else using covariance matrices based on the (v+1)*k units
            if ncolC==p
                startv1=0;
            elseif ncolC==k*(p)%prima era p+1
                startv1=1;
            else
                disp('If nsamp is not a scalar it must have v or k*(v+1) columns')
                disp('Please generate nsamp using')
                disp('nsamp=subsets(number_desired_subsets,n,k) or')
                disp('nsamp=subsets(number_desired_subsets,n,(v+1)*k)')
                error('FSDA:tclust:WrongNsamp','Wrong number of columns in matrix nsamp')
            end
            NoPriorSubsets=0;
        else
            % If nsamp is a scalar it simply contains the number of subsets
            % which have to be extracted. In this case NoPriorSubsets=1
            NoPriorSubsets=1;
            
            % In this case (nsamp is a scalar) we check whether the user
            % has supplied option startv1
            chkstartv1 = strcmp(varargin,'startv1');
            if sum(chkstartv1)>0
                startv1= cell2mat(varargin(find(chkstartv1)+1));
            else
                startv1=startv1def;
            end
        end
    else
        % If option nsamp is no supplied then for sure there are no prior
        % subsets
        NoPriorSubsets=1;
        
        % In this case (options nsamp does not exist) we check whether the
        % user has supplied option startv1
        chkstartv1 = strcmp(varargin,'startv1');
        if sum(chkstartv1)>0
            startv1= cell2mat(varargin(find(chkstartv1)+1));
        else
            startv1=startv1def;
        end
    end
else
    % if nargin ==6 for use the user has not supplied prior subsets.
    % Default value of startv1 is used
    NoPriorSubsets=1;
    startv1=startv1def;
end


% If the user has not specified prior subsets (nsamp is not a scalar) than
% according the value of startv1 we have a different value of ncomb
if NoPriorSubsets ==1
    % Remark: startv1 must be immediately checked because the calculation of
    % ncomb is immediately affected.
    
    if startv1
        ncomb=bc(n,k*(p+1));
    else
        % If the number of all possible subsets is <300 the default is to
        % extract all subsets otherwise just 300.
        % Notice that we use bc, a fast version of nchoosek. One may also
        % use the approximation
        % floor(exp(gammaln(n+1)-gammaln(n-p+1)-gammaln(p+1))+0.5)
        ncomb=bc(n,k);
    end
    nsampdef=min(300,ncomb);
end
%% Defaults for optional arguments

% default number of random starts
niterdef = 20;

% default number of concentration starts
Kstepsdef  = 10;

%default value for wtrim
wtrimdef = 0;

%default value for we
wedef = ones(n,1);

%default model (mixture or classification likelihood)
mixtdef = 2;

% default for threshold controlling the distance between regression lines
% in the initialization phase. Zero threshold means that there is no
% control on the initial fits.
eps_beta_def = 0;

%% User options

options = struct('intercept',1,'mixt',mixtdef,...
    'nsamp',nsampdef,'niter',niterdef,'Ksteps',Kstepsdef,...
    'startv1',startv1def,'we',wedef,'wtrim',wtrimdef,'eps_beta',eps_beta_def,...
    'msg',0,'plots',1);

if nargin > 6
    
    UserOptions = varargin(1:2:length(varargin));
    
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:tclustreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if all the specified optional arguments were present in
        % structure options. Remark: the nocheck option has already been dealt
        % by routine chkinputR.
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:tclustreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    % Write in structure 'options' the options chosen by the user
    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
    
    % And check if the optional user parameters are reasonable.
    
    % Check number of subsamples to extract
    if isscalar(options.nsamp) && options.nsamp>ncomb
        disp('Number of subsets to extract greater than (n k). It is set to (n k)');
        options.nsamp=0;
    elseif  options.nsamp<0
        error('FSDA:tclust:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
    end
    
    % Check restriction factor
    if restrfact<1
        disp('Restriction factor smaller than 1. It is set to 1 (maximum contraint==>spherical groups)');
        restrfact=1;
    end
    
end

%% set user's options

msg = options.msg;

% Graph summarizing the results
plots = options.plots;

% Intercept, yes/no
intercept = options.intercept;

% Number of subsets to extract
nsamp = options.nsamp;

% Concentration steps
Ksteps = options.Ksteps;

% Threshold controlling the distance between regression lines in the
% initialization phase
eps_beta   = options.eps_beta;

% the weights vector
we         = options.we;

% flag to control the type of weighting scheme
wtrim      = options.wtrim;

switch wtrim
    case 0
        if sum(we ~= wedef)>0
            disp('Warning: when "wtrim" is 0, "we" is set to a vector of ones');
            disp('         to give equal weights to all observations;');
            disp('         your vector "we" will not be considered.');
            we = wedef;
        end
    case 1
        %we must be a column vector);
        we = we(:);
        
        if sum(we == wedef)==n
            disp('Warning: when "wtrim" is 1, trimming is done by weighting');
            disp('         the observations using values specified in vector "we";');
            disp('         you left "we" to the default (i.e. a vector of ones,');
            disp('         giving equal weights to all observations);');
            disp('         we set them to a vector of 1/n, to sum to 1.');
        end
        % weights must be positive; if negative, values are translated so
        % that the minimum is 0
        if sum(we<0)>0
            we = we - min(we);
            disp('Warning: one or more of your weights are negative;');
            disp('         we added the minimum to all weights.');
        end
        % weights cannot be all equal to 0.
        if max(we) == 0
            we = wedef/n;
            disp('Warning: your weights are all zero;');
            disp('         we set them to a vector of 1/n, to sum to 1.');
        end
        
        % weights must be normalized so that to sum to 1
        we = we/sum(we);
        
    case 2
        if sum(we ~= wedef)>0
            disp('Warning: when "wtrim" is 2, trimming is done by weighting');
            disp('         the observations according to the data density estimate;');
            disp('         your vector "we" will not be considered.');
            we = wedef;
        end
    case 3
        if sum(we ~= wedef)>0
            disp('Warning: when "wtrim" is 3, trimming is done by weighting');
            disp('         the observations with a Bernoulli random vector,');
            disp('         with probability of success depending on the data density estimate;');
            disp('         your vector "we" will not be considered.');
            we = wedef;
        end
    case 4
        if sum(we ~= wedef)>0
            disp('Warning: when "wtrim" is 4, tclust is applied after thinning');
            disp('         observations with a Bernoulli random vector,');
            disp('         with probability of success depending on the data density estimate;');
            disp('         your vector "we" will not be considered.');
            we = wedef;
        end
end

%option determining the model to use
mixt = options.mixt;
if msg == 1
    switch mixt
        case 0
            %ClaLik + Crisp
            disp('ClaLik + Crisp');
        case 1
            % each unit is assigned to a group and then we take the h
            % best maxima
            %ClaLik + PostProb
            disp('ClaLik + PostProb');
        case 2
            % we take the units with the h largest contributions to the
            % likelihood
            %MixLik + PostProb
            disp('MixLik + PostProb');
    end
end

%% Additional variables depending on user options

% Number of variables without considering the constant term. It is used for
% deciding the type of plot.
if intercept == 1
    v = p-1;
else
    v = p;
end

% First level trimming
if alpha1<1
    notrim = floor(n*(1-alpha1));
else
    notrim = n-floor(alpha1);
end
trimm = n-notrim;

% Total trimming after second trimming
if alpha1<1
    if alpha2<1
        trimm2 = floor(n*(1-alpha1)*(1-alpha2));
    elseif alpha2>=1
        trimm2 = floor(n*(1-alpha1)*(1-alpha2/n));
    end
elseif alpha1>=1
    if alpha2 >= 1
        trimm2 = floor(n-floor(alpha1)-floor(alpha2));
    elseif alpha2 < 1
        trimm2 = floor(n*(1-alpha1/n)*(1-alpha2));
    end
end

%% Combinatorial part to extract the subsamples (if not already supplied by the user)


%case with no prior subsets
if NoPriorSubsets
    %if stratv1 =1 the initial subsets are formed by k*(p) observations
    if startv1 && k*(v+1) < n
        % the number of initial subsets to be generated is nsamp*oversamp.
        % The input parameter nsamp is multiplied by a factor (oversamp) in
        % order to face the possibility that some subsets contain groups
        % which  are very closed one to the other and therefore have to be
        % eliminated and substituted with other subsets.
        for ns =1:nsamp*oversamp
            C(ns,:) = datasample(1:n,k*(p),'weights',we,'Replace',false); %was p+1
        end
        nselected = length(C)/oversamp;
        %if stratv1 =0 the initial subsets are formed by k observations
    else
        for ns =1:nsamp*oversamp
            % the number of initial subsets to be generated is nsamp*oversamp.
            % The input parameter nsamp is multiplied by a factor (oversamp) in
            % order to face the possibility that some subsets contain groups
            % which  are very closed one to the other and therefore have to be
            % eliminated and substituted with other subsets.
            C(ns,:) = datasample(1:n,k,'weights',we,'Replace',false);
        end
        nselected  = length(C)/oversamp;
        %niinistart = repmat(floor(notrim/k),k,1);
    end
end

%% Initialize structures

ll         = zeros(n,k);
ll_old     = zeros(n,k);
ni         = ones(1,k);
sigmaopt   = ni;
bopt       = zeros(p,k);
numopt     = 1:k;
fact3      = zeros(n,k);

%%  Random starts

count      = 0;
count_elim = 0;
iter       = 0;
sigmaini   = ones(1,k);
comb_beta  = combnk(1:k,2);
diff       = NaN(intercept+1,size(comb_beta,1));
selj_good_groups = NaN(nselected,k*p);
selj_elim_groups = NaN(nselected,k*p);

while iter < nselected
    iter  = iter+1;
    count = count+1;
    % display the count of the number of all selected samples. In
    % line 1096 it is diplayed the number of only valid samples.
    if msg == 1
        disp(['Iteration ' num2str(count)])
    end
    
    if startv1
        
        nameYY = zeros(p,k);
        % rng(1234);
        randk=rand(k,1);
        if alpha1<1
            niini=floor(fix(n*(1-alpha1))*randk/sum(randk));
        else
            niini=floor(fix(n - floor(alpha1))*randk/sum(randk));
        end
        while sum(niini == 0) >0
            randk=rand(k,1);
            % Initialize niini with with random numbers from uniform
            if alpha1<1
                niini=floor(fix(n*(1-alpha1))*randk/sum(randk));
            else
                niini=floor(fix(n -floor(alpha1))*randk/sum(randk));
            end
        end
        if debug_mode==1
            X_all=0;
            y_all=0;
            id_gr_all=0;
        end
        for j = 1:k
            ilow   = (j-1)*(p)+1; %was p+1
            iup    = j*(p);       %was p+1
            index  = C(iter,:);
            selj   = index(ilow:iup);
            selj_good_groups(count,ilow:iup) = selj;
            Xb     = X(selj,:);
            yb     = y(selj,:);
            ni(j)  = length(yb);
            nameYY(:,j) = Xb\yb;
            % now find residuals
            residuals = yb-Xb*nameYY(:,j);
            % Update sigmas through the mean square residuals
            if size(selj,2) > p
                sigmaini(j) = sum(residuals.^2)/ni(j);
            else
                sigmaini(j) =var(y);
            end
            if debug_mode == 1
                X_all     = [X_all ; Xb(:,intercept+1)]; %#ok<AGROW>
                y_all     = [y_all ; yb];                %#ok<AGROW>
                id_gr     = repmat(j,length(yb'),1);
                id_gr_all = [id_gr_all; id_gr];          %#ok<AGROW>
            end
        end
        sigmaini= restreigen(sigmaini,ni',restrfact,tolrestreigen,userepmat);
        if debug_mode==1
            figure;
            scatter(X(:,intercept+1),y);
            hold on;
            gscatter(X_all,y_all,id_gr_all);
        end
    else
        
        % initialization of niini with equal proportions
        % niini=niinistart;
        
        % extract a subset of size v
        index = C(count,:);
        Xb = X(index,:);
        yb = y(index);
        
        for j=1:k
            Xbj = Xb((1+(j-1)*p):(j*p),:);
            ybj = yb((1+(j-1)*p):(j*p));
            nameYY(:,j) = Xbj\ybj;
        end
        % sigmaini will contain the covariance matrices in each iteration
        sigmaini=ones(1,k);
    end
    
    %to be generalized to a generic k = 2 3 4 5 6 ...
    for par = 1:intercept + 1
        for gr = 1:size(comb_beta,1)
            diff(par,gr) = nameYY(par,comb_beta(gr,1)) - nameYY(par,comb_beta(gr,2));
        end
    end
    mindiff = min(abs(diff),[],2);
    %if both intercept and slope
    if sum(abs(mindiff) > eps_beta) >0
        good_initial_subs = 1;
    else
        count_elim = count_elim + 1;
        good_initial_subs = 0;
        selj_elim_groups(count_elim,:) = selj_good_groups(count,:);
        count = count-1;
        nselected = nselected+1;
    end
    
    if good_initial_subs == 1
        
        %% Concentration steps
        indold = zeros(n,1)-1;
        for t = 1:Ksteps
            % Discriminant functions for the assignments
            for jk = 1:k
                fact3(:,jk) = logmvnpdfFS(y-X*nameYY(:,jk),0,(sigmaini(jk)));%,zeros(n,v),eye(v),n,v,0);
                %alternative way to compute the normal probability density function
                %fact2(:,jk) = normpdf(y-X*nameYY(:,jk),0,sqrt(sigmaini(jk)));
            end
            %if there are extreme (contaminated) observations and the
            %groups are almost collinear (small sigmaini), it can happen
            %that the area computed by the normpdf is zero for all the k
            %groups. In this case we contaminate the k values close to zero
            %in such a way that these points are randomly assigned to one
            %of the k groups
            extreme_obs = find(sum(fact3,2)==0);
            for jk = 1:k
                fact3(extreme_obs,jk) = ...
                    fact3(extreme_obs,jk)+0.0000000001*abs(rand(length(extreme_obs),1));
                ll(:,jk) = log((niini(jk)/sum(niini))) + fact3(:,jk);
            end
            
            
            
            % trimming when there is no observation weighting
            if wtrim == 0  || wtrim == 2 || wtrim == 3 || wtrim == 4
                
                if mixt == 2
                    [~,postprob,disc] = estepFS(ll);
                elseif mixt ==0
                    [disc,indmax] = max(ll,[],2);
                end
                
                % Sort the n likelihood contributions
                % qq contains the largest n*(1-alpha)
                % likelihood contributions
                [~,qq] = sort(disc,'descend');
                
                % qq = vector of size h which contains the indexes
                % associated with the largest n(1-alpha)
                % likelihood contributions
                qqunassigned = qq((n-trimm+1):n);
                qq           = qq(1:n-trimm);
                
                % trimming with user weighting
            elseif wtrim ==1
                
                if mixt ==2
                    [~,postprob,disc] = estepFS(ll);
                    dsf = sum(ll,2);
                elseif mixt == 0
                    [dsf,indmax] = max(ll,[],2);
                end
                [ ~, id] =sort((dsf));
                cumsumyy = cumsum(we(id));
                if alpha1<1
                    qqunassigned_small = cumsumyy < alpha1*sum(we(id));
                else
                    qqunassigned_small = cumsumyy < alpha1/n*sum(we(id));
                end
                qqunassigned = id(qqunassigned_small);
                qq = setdiff((1:n)',qqunassigned);
            end
            
            % ytri and Xtri = n(1-alpha)-by-v matrix associated with
            % the units which have the largest n(1-alpha) likelihood
            % contributions
            Xtri = X(qq,:);
            ytri = y(qq,:);
            wetri = we(qq,:);
            
            if mixt == 2
                postprob(qqunassigned,:) = 0;
                postprobtri = postprob(qq,:);
                
                % M-step update of niini
                % niini = numerator of component probabilities
                niini=(nansum(postprob))';
                
                %the next lines are needed to assign each observation to
                %the group with the largest posterior probability
                [~,indmax]= max(postprob,[],2);
            end
            
            indtri=indmax(qq);
            xmod=[Xtri , ytri , indtri];
            
            
            % size of the groups
            for jj=1:k
                ni(jj) = sum(indtri==jj);
            end
            
            
            % Update of posterior probabilities via observation
            % weighting (option wtrim = 1,2,3). Needed to compute beta
            % coefficients
            
            % if wtrim == 0, no observation weighting is done
            if wtrim == 0
                if mixt == 2
                    weight = postprobtri;
                elseif mixt == 0
                    weight = repmat(wetri,1,k);
                end
            end
            
            % if wtrim == 1, the weights are the posterior
            % probabilities multiplied by the user weights
            if wtrim == 1
                if mixt == 2
                    weight = postprobtri .* repmat(wetri,1,k);
                elseif mixt == 0
                    weight =   repmat(wetri,1,k);
                end
            end
            
            % if wtrim == 2 || wtrim == 3, the weights are the
            % posterior probabilities multiplied by the kernel density
            % or bernoulli weights estimated on a component basis,
            % which requires an additional loop over groups
            if wtrim == 2 || wtrim == 3
                
                % initialize weight with posterior probabilities
                if mixt == 2
                    weight = postprobtri;
                elseif mixt == 0
                    weight =  repmat(wetri,1,k);
                end
                
                for jj=1:k
                    % find indices of units in group jj
                    ijj = find(indtri==jj);
                    
                    % weight vector is updated only if the group is not
                    % empty
                    if ~isempty(ijj)
                        % retention probabilities based on density
                        % estimated on the component predicted values
                        % of the previous step
                        yhattri = Xtri(ijj,:)*nameYY(:,jj);
                        [Wt , pretain] = wthin(yhattri);
                        if wtrim == 2
                            wetri = pretain;
                        else
                            wetri = Wt;
                        end
                        weight(ijj,jj) = weight(ijj,jj) .* wetri;
                    end
                end
            end
            
            % rescale weights so that the sum of elements in each raw
            % sum to one
            if wtrim > 0 && wtrim < 4
                weightsum = sum(weight,2);
                weight = bsxfun(@rdivide,weight,weightsum);
            end
            
            weightmod = [weight, indtri ];
            % Update the beta coefficients, possibly considering
            % the observation weighting vector we, according to the
            % wtrim option
            %                 for jj=1:k
            %                     nameYY(:,jj) = ...
            %                         bsxfun(@times,Xtri, sqrt(weight(:,jj))) \ ...
            %                         bsxfun(@times,ytri ,sqrt(weight(:,jj)));
            %                 end
            
            
            % If a cluster is empty or contains less than p+1 elements,
            % stop the concentration steps (not_empty_g != 0).
            
            % new k regression parameters and new sigmas
            xmodtemp    = zeros(n,p+2);
            indxmodtemp = 0;
            not_empty_g = ~( ni <= p + 1 );
            %count number of times the number of groups is lt k
            if sum(not_empty_g) == k
                count1_ng_eq_k = count1_ng_eq_k + 1;
            else
                count1_ng_lt_k = count1_ng_lt_k + 1;
            end
            
            %% second level of trimming
            jk = 0;
            for iii = not_empty_g
                jk = jk+1;
                
                %check if a group is populated
                if iii == 1
                    xmodj = xmod(xmod(:,end)==jk,:);
                    
                    weightmodj = weightmod(weightmod(:,end) == jk,:);
                    % Perform the second trimming
                    
                    % qqs contains contains the indexes of untrimmed units for
                    % group j
                    if alpha2 == 0
                        qqs = 1:ni(jk);
                    else
                        % Find the units with the smallest h distances.
                        % Apply mcd on the x space (without the intercept
                        % if present).
                        % REMARK: This is by far the computationally most
                        % expensive instruction of tclustreg. More
                        % precisely, the dominant expensive function inside
                        % mcd is IRWLSmcd.
                        if intercept
                            if alpha2 < 1
                                RAW = mcd(xmodj(:,2:p),'bdp',alpha2,'msg',0);
                            elseif alpha2 >= 1
                                RAW = mcd(xmodj(:,2:p),'bdp',alpha2/n,'msg',0);
                            end
                        else
                            if alpha2 < 1
                                RAW = mcd(xmodj(:,1:p),'bdp',alpha2,'msg',0);
                            elseif alpha2 >= 1
                                RAW = mcd(xmodj(:,1:p),'bdp',alpha2/n,'msg',0);
                            end
                        end
                        [~,indmdsor] = sort(RAW.md);
                        if alpha2 < 1
                            qqs = indmdsor(1:floor(ni(jk)*(1-alpha2)));
                        else
                            qqs = indmdsor(1:floor(ni(jk) - alpha2));
                        end
                    end
                    
                    %% Update betas through ordinary least squares regression
                    xxx = xmodj(qqs,1:p);
                    yyy = xmodj(qqs,p+1);
                    ni(jk) = length(yyy);
                    
                    weightmodj_jk = sqrt(weightmodj(qqs,jk));
                    breg =  (bsxfun(@times,xxx, weightmodj_jk)) \ (bsxfun(@times,yyy ,weightmodj_jk));
                    
                    nameYY(:,jk) = breg;
                    % now find residuals
                    residuals = yyy-xxx*breg;
                    % Update sigmas through the mean square residuals
                    sigmaini(jk) = sum((residuals .* weightmodj_jk).^2)/(sum((weightmodj_jk).^2));
                    
                    
                    xmodtemp((indxmodtemp+1):(indxmodtemp+ni(jk)),:) = xmodj(qqs,:);
                    indxmodtemp = indxmodtemp+ni(jk);
                    
                else
                    
                    % xmodj Data points in each group
                    xmodj = [];
                    if mixt == 2
                        postprobmodj = [];
                        weightmodj = [];
                    end
                    % Perform the second trimming
                    
                    % qqs contains contains the indexes of untrimmed units
                    % for group j
                    if alpha2 == 0
                        qqs = [];
                    else
                        qqs = [];
                    end
                    
                    % Update betas through ordinary least squares regression
                    xxx = [];
                    yyy = [];
                    ni(jk) = 0;
                    breg = NaN;
                    nameYY(:,jk) = breg;
                    xmodtemp((indxmodtemp+1):(indxmodtemp+ni(jk)),:) = xmodj(qqs,:);
                    indxmodtemp = indxmodtemp+ni(jk);
                    % New xmod
                    % If the scatters do not satisfy the restriction then a
                    % quadratic programming problem is solved
                    sigmaini(jk) = NaN;
                    %count the number of times in a group there are enough
                    %observations to compute the sigma
                    count1_ng_eq_k = count1_ng_eq_k + 1;
                end
                
            end
            
            sigmaini= restreigen(sigmaini,ni',restrfact,tolrestreigen,userepmat);
            if debug_mode == 1
                figure;
                gscatter(xmod(:,2),xmod(:,3),xmod(:,4));
                hold on;
                line([0 max(xmod(:,2))],[nameYY(1,1) , nameYY(1,1) + nameYY(2,1)*max(xmod(:,2))])
                hold on;line([0 max(xmod(:,2))],[nameYY(1,2) , nameYY(1,2) + nameYY(2,2)*max(xmod(:,2))])
                hold on;line([0 max(xmod(:,2))],[nameYY(1,3) , nameYY(1,3) + nameYY(2,3)*max(xmod(:,2))])
                
                %variance(jk) = var(residuals);
                % disp(sigmaini);
                
                hold off;
            end
            %if a group is emty, beta and sigma are computed as mean of the
            %other groups
            for j=1:k
                if isnan(sigmaini(j))
                    sigmaini(j) = nanmean(sigmaini);
                end
                if isnan(nameYY(:,j))
                    nameYY(:,j) = nanmean(nameYY,2);
                end
            end
            xmod = xmodtemp(1:indxmodtemp,:);
            
            % Stop if two consecutive concentration steps have the same result
            if indmax == indold
                break
            else
                indold = indmax;
            end
            
            %% Now compute the value of the target function
            obj = 0;
            not_empty_g = ~( ni <= p + 1 );
            if mixt == 0
                % Update weights
                niini(jk) = ni(jk);
                
                jk = 0;
                for iii = not_empty_g
                    jk = jk+1;
                    if iii ==1
                        yj = xmod(xmod(:,end) == jk,end-1);
                        Xj = xmod(xmod(:,end) == jk,1:end-2);
                        %the following command should be executed at the
                        %end of the for loop of the concentration steps.
                        %However here it is executed all steps, because of
                        %a break from the loop which is executed if two
                        %consecutive concentration steps have the same
                        %result
                        % if t == Ksteps
                        %                            obj_old = obj + niini(jk)*log(niini(jk)/trimm2) +...
                        %                            sum(log(normpdf(yj-Xj*nameYY(:,jk),0,sqrt(sigmaini(jk)))));
                        obj = obj + niini(jk)*log(niini(jk)/trimm2) +...
                            sum(logmvnpdfFS(yj-Xj*nameYY(:,jk),0,(sigmaini(jk))));
                        % end
                    else
                        %if a groupis missing, we do not compute the objective
                        %function for it.
                    end
                end
            elseif mixt == 2
                %the following command should be executed at the end of
                %the for loop of the concentration steps. However here
                %it is executed all steps, because of a break from the
                %loop which is executed if two consecutive
                %concentration steps have the same result
                %if t == Ksteps
                log_lh=NaN(size(xmod,1),size(not_empty_g,2));
                jk = 0;
                %log_lh = [];
                for iii = not_empty_g
                    jk = jk+1;
                    if iii ==1
                        log_lh(:,jk) = ...
                            log(niini(jk)/sum(niini)) + (logmvnpdfFS(...
                            xmod(:,end-1) - ...
                            xmod(:,1:(size(xmod,2)-2)) * ...
                            nameYY(:,jk),0,(sigmaini(jk)) ) );
                        
                        
                    else
                        %if a groupis missing, we do not compute the objective
                        %function for it.
                        %obj = obj-10^10;
                        log_lh(:,jk) = NaN(length(xmod),1);
                    end
                end
                
                group_missing = sum(isnan(log_lh),1)>0;
                log_lh(:,group_missing)=[];
                obj = estepFS(log_lh);
                if ~isempty(group_missing)
                    nameYY(:,group_missing) = NaN;
                    sigmaini(group_missing) = NaN;
                end
                %end
            end
            
        end % End of concentration steps
        
        %% Change the 'optimal' target value and 'optimal' parameters
        % This is done if an increase in the target value is achieved
        
        %if sum(not_empty_g ) == k   %this check is commented in order to estimate the effect of eps_beta
        if sum(sum(isnan(nameYY))) == 0
            if (obj >= vopt)
                vopt = obj;
                bopt = nameYY;
                numopt = niini;
                sigmaopt = sigmaini;
            end
        end
        %end
        
        % display the count of the number of valid selected samples. In
        % line 629 it is diplayed the number of alla samples, valid and not
        % valid.
        if msg == 1
            disp(['iter ' num2str(count)]);
        end
    end
end % end of loop over the nsamp subsets
if count < nsamp
    out = struct;
else
    %% Prepares the output structure and some variables for the plots
    
    % Assignment vectors:
    % - asig.1 will contain the clusters after the first trimming
    % - asig.2 will contain the clusters after after the second trimming
    asig1 = zeros(n,1);
    asig2 = asig1;
    
    % log-likelihoods for each unit and group
    not_empty_g = ~( numopt == 0 )';
    jk = 0;
    for iii = not_empty_g
        jk = jk+1;
        if iii == 1
            ll(:,jk) = log((numopt(jk)/sum(numopt)) )+ logmvnpdfFS(y-X*bopt(:,jk),0,(sigmaopt(jk)));
        else
            ll(:,jk) = NaN;
        end
    end
    %compute posterior probabilities. In principle it should be computed only
    %for mixt==2, but since it is among the output of the function, it is
    %computed also for mixt==1
    [~,postprob,~] = estepFS(ll);
    
    %% Determine observations to trim
    
    % boolean vectors indicating the good and outlying units
    [dist,indmax] = max(ll,[],2);
    % Sort the n likelihood contributions;
    [val,qq] = sort(dist,'descend');
    
    if wtrim ==0
        % qq is updated to be a vector of size h which contains the indexes
        % associated with the largest n(1-alpha) (weighted) likelihood
        % contributions
        qq  = qq(1:n-trimm);
        val = val(n-trimm);
        
    elseif wtrim ==1
        qq_acend = qq(end:-1:1);
        cumsumyy = cumsum(we(qq_acend));
        if alpha1 <1
            qqunassigned_small = cumsumyy < alpha1*sum(we(qq_acend));
        else
            qqunassigned_small = cumsumyy < alpha1/n*sum(we(qq_acend));
        end
        qqunassigned = qq_acend(qqunassigned_small);
        qq = setdiff((1:n)',qqunassigned);
        val = val(n-length(qqunassigned));
        %b_outl_unasigned = b_outl(qqunassigned_small==1);
        %b_outl_asigned = b_outl(qqunassigned_small==0);
    end
    b_good = (dist>=val);
    b_outl = (dist <val);
    
    % asig1: grouping variable for good units, with 0 for trimmed units
    for jk=1:k
        asig1((indmax == jk) & b_good) = jk;
    end
    
    % xmod: contains the good units and, in the last column, their group assignment
    xmod = [X(qq,:) y(qq) indmax(qq)];
    
    % used for collecting sencond level trimming points, used for plotting
    %xxx0_all = [];
    %yyy0_all = [];
    
    % IS THIS PART BELOW MADE JUST TO COUNT the number of times the number
    % of groups is lt k? IF YES, SHOULD THIS BE MOVED AT THE END OF THE
    % LOOP BELOW, AFTER SECOND LEVEL TRIMMING?
    % DOME BEGIN
    % go over the groups
    for jk = 1:k
        booljk = xmod(:,end) == jk;
        ni(jk) = sum(booljk);
    end
    
    not_empty_g = ~( ni <= p + 1 );
    %count the number of times the number of groups is lt k
    if sum(not_empty_g) == k
        count2_ng_eq_k = count2_ng_eq_k + 1;
    else
        count2_ng_lt_k = count2_ng_lt_k + 1;
    end
    % DOME END
    
    %% determine second level trimming points
    
    jk = 0;
    for iii = not_empty_g
        jk = jk+1;
        booljk = xmod(:,end) == jk;
        %b_outl_asigned_jk = b_outl_asigned(booljk);
        qqk = qq(booljk);
        ni(jk) = sum(booljk);
        %xmodjk = xmod(booljk & b_outl_asigned_jk,:);
        xmodjk = xmod(booljk ,:);
        if iii ==1
            if alpha2 == 0
                qqs = 1:ni(jk);
            else
                if intercept
                    if alpha2 <1
                        RAW = mcd(xmodjk(:,2:p),'bdp',alpha2,'msg',0);
                    elseif alpha2 >= 1
                        RAW = mcd(xmodjk(:,2:p),'bdp',alpha2/n,'msg',0);
                    end
                else
                    if alpha2 <1
                        RAW = mcd(xmodjk(:,1:p),'bdp',alpha2,'msg',0);
                    elseif alpha2 >= 1
                        RAW = mcd(xmodjk(:,1:p),'bdp',alpha2/n,'msg',0);
                    end
                end
                [~,indmdsor] = sort(RAW.md);
                qqs = indmdsor(1:floor(ni(jk)*(1-alpha2)));
            end
            
            % good units of the current group (used for plotting inside loop)
            %xxx = xmodjk(qqs,1:end-2);
            %yyy = xmodjk(qqs,end-1);
            %
            % second level trimming units of the current group
            %qqsn = setdiff(1:ni(jk),qqs);
            %xxx0 = xmodjk(qqsn,1:end-2);
            %yyy0 = xmodjk(qqsn,end-1);
            %
            % collect all second level trimming units in a same group, for plotting
            %xxx0_all = [xxx0_all ; xxx0(:,end)]; %#ok<AGROW>
            %yyy0_all = [yyy0_all ; yyy0];        %#ok<AGROW>
            
            qqf = qqk(qqs);
            asig2(qqf) = jk;
        end
        
        % plots    MOVED OUTSIDE THE LOOP, TOGETHER WITH OTHER PLOTS
        %
        %         % The following plots are for the bi-variate case (v=1)
        %         if (plots && v < 2)
        %
        %             % initialize figure
        %             if jk == 1
        %                 fh = figure('Name','TclustReg plot','NumberTitle','off','Visible','on');
        %                 gca(fh);
        %                 hold on;
        %                 xlabel('X');
        %                 ylabel('y');
        %                 title('TclustReg clustering','Fontsize',14);
        %                 %print(fh,[graphs 'PCvariance.png'],'-dpng');
        %                 %close(fh);
        %             end
        %
        %             % plot good units allocated to the current group
        %             group_label = ['Group ' num2str(jk)];
        %             plot(xxx(:,end),yyy,'.w','DisplayName',group_label);
        %             text(xxx(:,end),yyy,num2str(jk*ones(length(yyy),1)),...
        %                 'DisplayName',group_label , ...
        %                 'HorizontalAlignment','center',...
        %                 'VerticalAlignment','middle',...
        %                 'Color',clrdef(jk));
        %
        %             % second level trimming points are not plotted by group
        %             % plot(xxx0(:,end),yyy0,'*','color','c','DisplayName','Level-2 trim');
        %
        %             % plot regression lines
        %             vv = [min(X(:,end)) max(X(:,end))];
        %             if intercept==1
        %                 plot(vv,bopt(1,jk)+bopt(2,jk)*vv,...
        %                     'DisplayName',['fit of group '  num2str(jk)],...
        %                     'Color',clrdef(jk));
        %             elseif intercept==0
        %                 plot(vv,bopt(:,jk)*vv,...
        %                     'DisplayName',['fit of group '  num2str(jk)],...
        %                     'Color',clrdef(jk));
        %             end
        %
        %         end
        
    end
    
    %% Generate plots
    
    if plots
        
        % The following plots are for the bi-variate case (i.e. v=1)
        if v < 2
            
            % initialize figure
            fh = figure('Name','TclustReg plot','NumberTitle','off','Visible','on');
            gca(fh);
            hold on;
            xlabel('X');
            ylabel('y');
            title('TclustReg clustering','Fontsize',14);
            
            jk = 0;
            for iii = not_empty_g
                jk = jk+1;
                if iii>0
                    group_label = ['Group ' num2str(jk)];
                    
                    % plot of the good units allocated to the current group.
                    % Indices are taken after the second level trimming.
                    % Trimmed points are not plotted by group.
                    ucg = find(asig2==jk);
                    plot(X(ucg,end),y(ucg),'.w','DisplayName',group_label);
                    text(X(ucg,end),y(ucg),num2str(jk*ones(length(ucg),1)),...
                        'DisplayName',group_label , ...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle',...
                        'Color',clrdef(jk));
                    
                    % plot regression lines
                    vv = [min(X(:,end)) max(X(:,end))];
                    if intercept==1
                        plot(vv,bopt(1,jk)+bopt(2,jk)*vv,...
                            'DisplayName',[group_label ' fit'],...
                            'Color',clrdef(jk));
                    elseif intercept==0
                        plot(vv,bopt(:,jk)*vv,...
                            'DisplayName',[group_label ' fit'],...
                            'Color',clrdef(jk));
                    end
                end
            end
            
            % Plot the outliers (trimmed points)
            b_outl = (asig1==0);
            plot(X(b_outl,end),y(b_outl),'o','color','r',...
                'DisplayName','Trimmed units');
            
            % second level trimming points
            b_outl_2 = ~(asig1==asig2);
            xxx0_all = X(b_outl_2,end);
            yyy0_all = y(b_outl_2);
            plot(xxx0_all,yyy0_all,'*','color','c',...
                'DisplayName','L2 trimmed units');
            
            % position the legends and make them clickable
            lh=legend('show');
            %set(lh,'FontSize',14);
            axis('manual');
            legstr = get(lh,'String');
            clickableMultiLegend(legstr,'FontSize',14,'Location','northwest');
            %[hleg, hobj, hout, mout] = clickableMultiLegend(legstr,'FontSize',14,'Location','northwest');
            %Unfortunately custom markers for line objects are not possible in MATLAB
            %set(hobj(10),'Marker','1','Color','k');
            
        else
            
            % in this case p > 2 and a standard spmplot is used
            
            if intercept
                YY = [X(:,2:end),y];
            else
                YY = [X,y];
            end
            
            % axis labels
            nameYY = cellstr([repmat('X',size(YY,2)-1,1) , num2str((1:size(YY,2)-1)')]);
            nameYY = [nameYY ; 'y'];
            nameYY = nameYY';
            plo=struct;
            plo.nameY=nameYY;
            
            % group names in the legend
            group = cell(size(asig2,1),1);
            group(asig2==0) = {'Trimmed units'};
            for iii = 1:k
                group(asig2==iii) = {['Group ' num2str(iii)]};
            end
            
            % scatterplot
            spmplot(YY,group,plo,'hist');
            
            %group_l = cellstr([repmat('Group',k,1) , num2str((1:k)')]);
            %group_l = ['Trimmed units' ; group];
            %[hleg, hobj, hout, mout] =legend((out(1,end,:)));
        end
    end
    
    % If the scatters do not satisfy the restriction then a quadratic
    % programming problem is solved
    
    if sum(isnan(sigmaopt)) == 0
        sigmaopt_0 = restreigen(sigmaopt,ni',restrfact,tolrestreigen,userepmat);
        % the old (inefficient) approach was to use a quadratic programming
        % optimization, using function quadi
        %sigmaopt_0 = (quadi(sigmaopt.^(-1), restrfact)).^(-1);
    end
    
    %Apply consistency factor based on the variance of the truncated normal
    %distribution.
    % hh = sum(numopt) number of non trimmed observations, after first and
    % second level trimming.
    % 1-hh/n=trimming percentage
    hh = sum(numopt);
    % Compute variance of the truncated normal distribution.
    vt = norminv(0.5*(1+hh/n));
    %factor=1/sqrt(1-(2*vt.*normpdf(vt))./(2*normcdf(vt)-1));
    factor = 1/sqrt(1-2*(n/hh)*vt.*normpdf(vt));
    % Note that factor=sqrt(factor1)
    %     v=1;
    %     a=chi2inv(hh/n,1);
    %     factor1=(hh/n)/(chi2cdf(a,1+2));
    % Apply the asymptotic consistency factor to the preliminary scale estimate
    if ~isnan(factor)
        sigmaopt_cons=sigmaopt_0*factor;
        % Apply small sample correction factor of Pison et al.
        sigmaopt_pison=sigmaopt_cons*sqrt(corfactorRAW(1,n,hh/n));
    else
        sigmaopt_cons=sigmaopt_0;
        sigmaopt_pison=sigmaopt_cons;
    end
    
    %%  Set the output structure
    
    out                = struct;
    out.bopt           = bopt;
    out.sigmaopt_0     = sigmaopt_0;
    out.sigmaopt_cons  = sigmaopt_cons;
    out.sigmaopt_pison = sigmaopt_pison;
    out.numopt         = numopt;
    out.vopt           = vopt;
    out.asig1          = asig1;
    out.asig2          = asig2;
    out.postprob       = postprob;
    out.count1_ng_lt_k = count1_ng_lt_k;
    out.count2_ng_lt_k = count2_ng_lt_k;
    out.count1_ng_eq_k = count1_ng_eq_k;
    out.count2_ng_eq_k = count2_ng_eq_k;
    out.extra_inisubs  = nselected - nsampdef;
    out.selj_good      = selj_good_groups;
    out.selj_elim      = selj_elim_groups(1:count_elim,:);
    out.selj_all       = [selj_elim_groups(1:count_elim,:); selj_good_groups];
    
    %   bopt           = regression parameters
    %   sigmaopt0      = estimated group variances
    %   sigmaopt_cons  = estimated group variances corrected with  asymptotic consistency factor
    %   sigmaopt_pison = estimated group variances corrected with  asymptotic consistency factor and small sample correction factor of Pison et al.
    %   numopt         = number of observations in each cluster after the second trimming
    %   vopt           = value of the target function
    %   asig1          = cluster assigments after first trimming ('0' means a trimmed observation)
    %   asig2          = (-final-) cluster assigments after second trimming ('0' means a trimmed observation)
    %   postprob       = posterior probability
    %   count1_ng_lt_k = number of times that, after the first level of trimming, in a group there are not enought observations to compute the sigma
    %   count1_eq_lt_k = number of times that, after the first level of trimming, in a group there are enought observations to compute the sigma
    %   count2_ng_lt_k = number of times that, after the second level of trimming, in a group there are not enought observations to compute the sigma
    %   count2_eq_lt_k = number of times that, after the second level of trimming, in a group there are enought observations to compute the sigma
    %   extra_inisubs  = number of subsets generated above the number specified by the user (nsamp) because of small difference between pairwise regression parameters
    %   out.selj_good  = list of valid subsets and observations inside them
    %   out.selj_elim  =  list of not-valid subsets and observations inside them
    %   out.selj_all   =  list of all subsets (valid and not) and observations inside them
    
end

%% Subfunctions

% corfactorRAW function
    function rawcorfac = corfactorRAW(p,n,alpha)
        
        if p > 2
            coeffqpkwad875=[-0.455179464070565,1.11192541278794,2;-0.294241208320834,1.09649329149811,3]';
            coeffqpkwad500=[-1.42764571687802,1.26263336932151,2;-1.06141115981725,1.28907991440387,3]';
            y1_500=1+(coeffqpkwad500(1,1)*1)/p^coeffqpkwad500(2,1);
            y2_500=1+(coeffqpkwad500(1,2)*1)/p^coeffqpkwad500(2,2);
            y1_875=1+(coeffqpkwad875(1,1)*1)/p^coeffqpkwad875(2,1);
            y2_875=1+(coeffqpkwad875(1,2)*1)/p^coeffqpkwad875(2,2);
            y1_500=log(1-y1_500);
            y2_500=log(1-y2_500);
            y_500=[y1_500;y2_500];
            A_500=[1,log(1/(coeffqpkwad500(3,1)*p^2));1,log(1/(coeffqpkwad500(3,2)*p^2))];
            coeffic_500=A_500\y_500;
            y1_875=log(1-y1_875);
            y2_875=log(1-y2_875);
            y_875=[y1_875;y2_875];
            A_875=[1,log(1/(coeffqpkwad875(3,1)*p^2));1,log(1/(coeffqpkwad875(3,2)*p^2))];
            coeffic_875=A_875\y_875;
            fp_500_n=1-(exp(coeffic_500(1))*1)/n^coeffic_500(2);
            fp_875_n=1-(exp(coeffic_875(1))*1)/n^coeffic_875(2);
        else
            if p == 2
                fp_500_n=1-(exp(0.673292623522027)*1)/n^0.691365864961895;
                fp_875_n=1-(exp(0.446537815635445)*1)/n^1.06690782995919;
            end
            if p == 1
                fp_500_n=1-(exp(0.262024211897096)*1)/n^0.604756680630497;
                fp_875_n=1-(exp(-0.351584646688712)*1)/n^1.01646567502486;
            end
        end
        if 0.5 <= alpha && alpha <= 0.875
            fp_alpha_n=fp_500_n+(fp_875_n-fp_500_n)/0.375*(alpha-0.5);
        end
        if 0.875 < alpha && alpha < 1
            fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
        end
        if alpha < 0.5
            fp_alpha_n = 1;
            if msg==1
                disp('Warning: problem in subfunction corfactorRAW')
                disp('alpha < 0.5')
            end
        end
        rawcorfac=1/fp_alpha_n;
        if rawcorfac <=0 || rawcorfac>50
            rawcorfac=1;
            if msg==1
                disp('Warning: problem in subfunction corfactorRAW')
                disp(['Correction factor for covariance matrix based on simulations found =' num2str(rawcorfac)])
                disp('Given that this value is clearly wrong we put it equal to 1 (no correction)')
                disp('This may happen when n is very small and p is large')
            end
        end
    end

% Subfunction quadi prepares the quantities to call the matlab
% quadratic programming routine quadprog. It was used to constrain the
% scatters which do not satisfy the desired restriction. Now the
% function is replaced by the more efficient restreigen.m
%
%     function gnew = quadi(gg,factor)
%         if size(gg,1)>1
%             gg=gg';
%         end
%
%         % gnew will the new scatters
%         gnew = gg;
%
%         if (length(gg)>1)
%             %Sort scatters
%             [ggsor,ggsorind] = sort(gg);
%
%             % g(1) = smallest sigma
%             % ...
%             % g(end) = largest sigma
%             g = ggsor;
%
%             maximun = 10^5;
%
%             % Constant "c" defining the scatter constraint
%             factor = factor+0.0001;
%
%             % nscat is the number of scatter parameters
%             nscat = length(g);
%
%             Amat =zeros(nscat,nscat);
%             % rr = 1:nscat;
%
%             for ii =1:(nscat-1)
%                 Amat(ii,ii) = -1;
%                 Amat(ii,ii+1) =1;
%             end
%
%             % Definition of the quadratic problem
%             Amat(nscat,1) = factor;
%             Amat(nscat,nscat) = -1;
%             Vmat = diag([ones(nscat,1);zeros(nscat,1)]);
%             dvec = - [g,zeros(1,nscat)];
%             bvec = zeros(1,nscat);
%             uvecmax = maximun+zeros(1,2*nscat);
%             uvecmin = zeros(1,2*nscat);
%
%             Amat = [Amat,-1*eye(nscat)];
%
%             % Solve this quadratic problem
%             % a = quadprog(Vmat,dvec,[],[],Amat,bvec',uvecmin,uvecmax,g,'Algorithm','interior-point-convex');
%
%             % FSDATOAPP:tclustreg:DF
%             % Remark: for compatibilty with old version of MATLAB we use
%             % intruction optimset. However recent versions of Matlab accept
%             % function optimoptions as follows
%             % option = optimoptions('quadprog','algorithm','interior-point-convex','Display','off');
%             option = optimset('OutputFcn','quadprog','algorithm','interior-point-convex','Display','off');
%
%             a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],option);
%             %a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],'algorithm','interior-point-convex','Display','iter');
%             %a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],'algorithm','active-set');
%
%             gnew =a(1:nscat);
%
%             %Original order
%             gnew(ggsorind) = gnew;
%
%         end
%     end
% end

end

