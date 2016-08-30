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
%            In mixture modelling, the likelihood is given by
%           \prod_{i=1}^n \left[ \sum_{j=1}^k \pi_j \phi (x_i;\theta_j)  \right].
%            In crisp assignment, the  likelihood is given by
%           \prod_{j=1}^k  \prod_{i\in R_j} \phi (x_i;\theta_j)
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
%   wtrim: Application of observation weights. Scalar. A flag taking values [0, 1, 2, 3, 4]
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
    
    mixt = 0; wtrim = 0;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    mixt = 2; wtrim = 0;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    mixt = 0; wtrim = 1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'we',we,'wtrim',wtrim);

    mixt = 2; wtrim = 1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'we',we,'wtrim',wtrim);

    mixt = 0; wtrim = 2; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    mixt = 2; wtrim = 2; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    mixt = 0; wtrim = 3; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    mixt = 2; wtrim = 3; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    mixt = 0; wtrim = 4; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    mixt = 2; wtrim = 4; we = [];
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

%%  computational issues to be addressed in future releases

% MATLAB function datasample. For the moment, we have copied the function
% in folder FSDA/combinatorial, renamed it as datasampleFS, and removed the
% computationally expensive option parameters checks. Unfortunately, we had
% to copy in the same folder the mex file wswor, also renamed as wsworFS.
% The function should be re-written along the linels of
%      Wong, C.K. and M.C. Easton (1980) "An Efficient Method for Weighted
%      Sampling Without Replacement", SIAM Journal of Computing,
%      9(1):111-113.

% FSDA function wthin uses the MATLAB function ksdensity. The calls to
% ksdensity have been optimized. The only possibility to further reduce
% time execution is to replace ksdensity with a better kernel density estimator.

% In the plots, the use of text to highlight the groups with their index is
% terribly slow (more than 10 seconds to generate a scatter of 7000 units.
% ClickableMultiLegend and legend are also slow.

% FSDA function restreigen could be improved. In some cases it is one of
% the most expensive functions.

%% Check if optimization toolbox is installed in current computer
% to be done in next releases: introduce an optimizer

typemin = exist('fminunc','file');
if typemin ~=2
    error('FSDA:tclustreg:MissingOptToolbox','This function requires the optimization toolbox');
end


%% initializations

warning('off'); %#ok<WNOFF>

%Thinning threshold. A group with less than thinning_th units is not
%considered for thinning.
thinning_th = 50;

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
clrdef = 'bkmgrcbkmgrcbkmgrcbkmgrcbkmgrcbkmgrcbkmgrc';
symdef = '+*sd^v><pho*';
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

% check if wtrim option has been set by the user equal to 4. In such a case, the thinning step has
% to be done at the very beginning, in particular before setting the number of samples nsamp.
if nargin>6
    chknwtrim = strcmp(varargin,'wtrim');
    if sum(chknwtrim)>0
        if cell2mat(varargin(find(chknwtrim)+1)) == 4
            interc = find(max(X,[],1)-min(X,[],1) == 0);
            if p == numel(interc) +1
                
                %apply thinning on the full dataset if there is only one exploratory variable.
                [Wt4,~] = wthin([X(:,numel(interc)+1),y], 'retainby','comp2one');
                
                % save original data
                Xori = X;
                yori = y;
                % set retained data
                X     = X(Wt4,:);
                y      = y(Wt4);
                %recompute n
                n = size(y,1);
            end
        end
    end
end

% check restrfact option
if nargin < 4 || isempty(restrfact) || ~isnumeric(restrfact)
    restrfact = 12;
end

% checks on alpha1 and alpha2
if alpha1<0
    error('FSDA:tclust:WrongAlpha','alpha1 must a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
end

if alpha2<0
    error('FSDA:tclust:WrongAlpha','alpha2 must a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
end

% startv1def = default value of startv1 = 1
% initialization using covariance matrices based on v+1 units
startv1def = 1;

% nsamp option: is it is not set by the user, it has to be properly initialized
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

% 'oversamp' is a factor (which depends on the number of groups 'k') used
% to generate more samples in order to face the possibility that some
% subsets contain collinear obs and are therefore rejected.
% To obtain nsamp = 300 samples, 300*oversamp samples will be generated.
%If eps_beta = 0, i.e. selection of subsets is done indipendently from the
%values of the regression parameters, there is no reason to
%ovserample and therefore oversamp is equal to 1.
if eps_beta>0
    oversamp = 10*k;
else
    oversamp = 10;
end

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
        we = we/nansum(we);
        
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
        if length(we) ~= length(wedef)
            we = we(Wt4);
        end
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

% % Total trimming after second trimming
% if alpha1<1
%     if alpha2<1
%         trimm2 = floor(n*(1-alpha1)*(1-alpha2));
%     elseif alpha2>=1
%         trimm2 = floor(n*(1-alpha1)*(1-alpha2/n));
%     end
% elseif alpha1>=1
%     if alpha2 >= 1
%         trimm2 = floor(n-floor(alpha1)-floor(alpha2));
%     elseif alpha2 < 1
%         trimm2 = floor(n*(1-alpha1/n)*(1-alpha2));
%     end
% end

%% Combinatorial part to extract the subsamples (if not already supplied by the user)

% Note that we used function 'datasample' instead of our FSDA function
% 'subsets' because we need weighted sampling.
%      Wong, C.K. and M.C. Easton (1980) "An Efficient Method for Weighted
%      Sampling Without Replacement", SIAM Journal of Computing,
%      9(1):111-113.

%case with no prior subsets
if NoPriorSubsets
    
    %if stratv1 =1 the initial subsets are formed by k*(p) observations
    if startv1 && k*(v+1) < n
        %if stratv1 =1 the initial subsets are formed by k*p observations
        initial_subs_size = k*p;
    else
        %if stratv1 =0 the initial subsets are formed by k observations
        initial_subs_size = k;
    end
    % the number of initial subsets to be generated is nsamp*oversamp.
    % The input parameter nsamp is multiplied by a factor (oversamp) in
    % order to face the possibility that some subsets contain groups
    % which  are very closed one to the other and therefore have to be
    % eliminated and substituted with other subsets.
    for ns =1:nsamp*oversamp
        %C(ns,:) = datasample(1:n,initial_subs_size,'weights',we,'Replace',false);
        C(ns,:) = randsampleFS(n,initial_subs_size,we);
    end
    nselected  = length(C)/oversamp;
    %niinistart = repmat(floor(notrim/k),k,1);
    
end


check_obj_reached = 0;
while check_obj_reached == 0
    %% Initialize structures

    ll                                       = zeros(n,k);
    ni                                      = ones(1,k);
    sigmaopt                         = ni;
    bopt                                 = zeros(p,k);
    numopt                            = zeros(1,k);
    weopt                               = ones(n,1);
    trim1levelopt                   = ones(n,1);
    trim2levelopt                   = ones(n,1);
    indmaxopt                        = zeros(n,1);
    not_empty_gopt              = zeros(1,k);
    count1_ng_lt_kopt          = 0;
    count2_ng_lt_kopt          = 0;
    count1_ng_eq_kopt       = 0;
    count2_ng_eq_kopt       = 0;
    extra_inisubsopt             = nselected - nsampdef;
    selj_goodopt                   = zeros(nselected,k);
    selj_elimopt                     = zeros(nselected,k);
    selj_allopt                        = zeros(nselected,k);
    if mixt ==2
        postprobopt                 = zeros(n,k);
    end

    count                                = 0;
    count_elim                       = 0;
    iter                                    = 0;
    sigmaini                            = ones(1,k);
    comb_beta                       = combnk(1:k,2);
    diff                                     = NaN(intercept+1,size(comb_beta,1));
    selj_good_groups            = NaN(nselected,k*p);
    selj_elim_groups              = NaN(nselected,k*p);
    
    %%  Random starts
    while iter < nselected
        %iter =  iteration number including the steps where the subset is
        %refused because the regression lines are too closed one to the other
        iter  = iter+1;

        %count =  iteration number not including the steps where the subset is
        %refused because the regression lines are too closed one to the other
        count = count+1;

        if msg == 1
            disp(['Iteration ' num2str(count)])
        end

        if startv1

            nameYY = zeros(p,k);
            % in order to fix the seed decomment the following command
            % rng(1234);
            randk=rand(k,1);
            if alpha1<1
                niini=floor(fix(n*(1-alpha1))*randk/nansum(randk));
            else
                niini=floor(fix(n - floor(alpha1))*randk/nansum(randk));
            end
            while sum(niini == 0) >0
                randk=rand(k,1);
                % Initialize niini with with random numbers from uniform
                if alpha1<1
                    niini=floor(fix(n*(1-alpha1))*randk/nansum(randk));
                else
                    niini=floor(fix(n -floor(alpha1))*randk/nansum(randk));
                end
            end
            for j = 1:k
                ilow   = (j-1)*(p)+1;
                iup    = j*(p);
                index  = C(iter,:);
                selj   = index(ilow:iup);
                selj_good_groups(count,ilow:iup) = selj;
                Xb     = X(selj,:);
                yb     = y(selj,:);
                ni(j)  = length(yb);
                nameYY(:,j) = Xb\yb;
                %Update residuals
                residuals = yb-Xb*nameYY(:,j);
                % Update sigmas through the mean square residuals
                if size(selj,2) > p
                    sigmaini(j) = nansum(residuals.^2)/ni(j);
                else
                    sigmaini(j) =var(y);
                end
            end
            sigmaini= restreigen(sigmaini,ni',restrfact,tolrestreigen,userepmat);
        else

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

        %compute the differences between pairwise regression parameters
        for par = 1:intercept + 1
            for gr = 1:size(comb_beta,1)
                diff(par,gr) = nameYY(par,comb_beta(gr,1)) - nameYY(par,comb_beta(gr,2));
            end
        end

        %if the differences between regression parameters is lower than a
        %threshold eps_beta the subset will not be considered and substitute
        %with a new generated one
        mindiff = min(abs(diff),[],2);
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

            % CONCENTRATION STEPS

            indold = zeros(n,1)-1;
            check_beta = 0;
            for cstep = 1:Ksteps
                %all beta coefficient are NAN when from a precedent concentration step, almost all obs
                %have been assigned to the same group and the thinning was very strong, leaving few obs.
                %In this case we force to exit from the concentration step
                not_empty_g = zeros(1,k);
                if sum(isnan(nameYY)) == k
                    if check_beta == 0
                        count_elim = count_elim + 1;
                        selj_elim_groups(count_elim,:) = selj_good_groups(count,:);
                        count = count-1;
                        nselected = nselected+1;
                        check_beta = 1;
                    end

                else   

                    %% Log-likelihood

                    % Discriminant functions for the assignments
                    for jk = 1:k
                        ll(:,jk) = log((niini(jk)/sum(niini))) + logmvnpdfFS(y-X*nameYY(:,jk),0,(sigmaini(jk)));
                    end
                    %to compute the normal probability density function with
                    %normpdf(y-X*nameYY(:,jk),0,sqrt(sigmaini(jk))) instead of logmvnpdfFS leads to
                    %imprecise results in the queue of the distribution. For example, if there are extreme
                    %outliers and the groups are almost collinear (small sigmaini), it can happen that the
                    %area computed by the normpdf is zero for all the k groups. In this case we would have
                    %to perturb the k values close to zero in such a way that these points are randomly
                    %assigned to one of the k groups. This can be done with the following code:
                    %                 for jk = 1:k
                    %                    fact2(:,jk) = logmvnpdfFS(y-X*nameYY(:,jk),0,(sigmaini(jk)));
                    %                 end
                    %                extreme_obs = find(sum(fact2,2)==0);
                    %                 for jk = 1:k
                    %                    if ~isempty(extreme_obs)
                    %                        fact2(extreme_obs,jk) = fact2(extreme_obs,jk)+10^(-8)*abs(rand(length(extreme_obs),1));
                    %                    end
                    %                    ll(:,jk) = log((niini(jk)/sum(niini))) + fact2(:,jk);
                    %                 end

                    %In the case of crisp assignement we compute the maximum value of the log-likelihood
                    %(disc) among the k groups, as in point 2.1 of appendix of Garcia Escudero et al.
                    %(2010). In the case of mixture modelling we need an estep to compute the log of the sum
                    %of the log-likelihood (disc) and the normalized log-likelihood (postprob).
                    if mixt == 2
                        [~,postprob,disc] = estepFS(ll);
                    elseif mixt == 0
                        [disc,indmax] = max(ll,[],2);
                    end

                    % Sort the n likelihood contributions and save in qq the largest n*(1-alpha) likelihood
                    % contributions
                    [~,qq] = sort(disc,'ascend');


                    %% first level trimming
                    switch wtrim
                        case {0,4};
                            % standard case, without observation weighting

                            %[~,qq] = sort(disc,'descend');
                            %qqunassigned = qq((n-trimm+1):n);
                            %qq           = qq(1:n-trimm);

                            % qq = vector which contains the indexes associated with the largest n(1-alpha)
                            % likelihood contributions
                            % qqunassigned is vector which contains the indices of the remaining trimmed
                            % observations
                            qqunassigned = qq(1:trimm);
                            qqassigned   = qq(trimm + 1:n);

                        case {1,2,3};
                            % trimming with observation weighting, set by the user or density estimation
                            %if wtrim = 3 the weights should be re-initialized to one, to avoid the
                            %continuous elimination of observations from the estimates.
                            if wtrim == 3
                                we = wedef;
                            end
                            cumsumyy = cumsum(we(qq));
                            if alpha1<1
                                qqunassigned_small = cumsumyy < alpha1*nansum(we(qq));
                                qq_small = cumsumyy >= alpha1*nansum(we(qq));
                            else
                                qqunassigned_small = cumsumyy < alpha1/n*nansum(we(qq));
                                qq_small = cumsumyy >= alpha1/n*nansum(we(qq));
                            end
                            qqunassigned = qq(qqunassigned_small);
                            qqassigned = qq(qq_small);
                            % qq_small introduced because setdiff below is inefficient
                            %qq = setdiff((1:n)',qqunassigned);
                    end
                    qqunassigned = sort(qqunassigned);
                    qqassigned = sort(qqassigned);

                    %observations trimmed with the 1st level trimming in original scale
                    trim1level = qqunassigned;

                    % In case of mixture modeling:
                    if mixt == 2
                        % update the posterior probabilities
                        postprobuntri = postprob(qqassigned,:);
                        postprob(qqunassigned,:) = 0;

                        % M-step: update of niini, the numerator of component probabilities
                        niini=(nansum(postprob))';

                        % indmax assigns each observation to the group with the largest posterior probability
                        [~,indmax]= max(postprob,[],2);
                    end

                    % data and observation weights vectors associated with the units which have the
                    % largest n(1-alpha) likelihood contributions
                    Xuntri              = X(qqassigned,:);
                    Xtri                  = X(qqunassigned,:);
                    yuntri               = y(qqassigned,:);
                    weuntri            = we(qqassigned,:);
                    induntri            = indmax(qqassigned);
                    indtri                = indmax(qqunassigned);
                    xmoduntri        = [Xuntri , yuntri , induntri];

                    % size of the groups nj (could be named mixing proportions or group weights)
                    %histcount is more efficient but cannot be used because when one group is missing it
                    %rdoes not report "0" but missing. Therefore the result is a vector of length less than
                    %k.
                    %             if verLessThan('matlab','8.4')
                    %                 for jj=1:k
                    %                     ni(jj) = sum(indtri==jj);
                    %                 end
                    %             else
                    %                 ni = histcounts(indtri,k);
                    %             end
                    for jj=1:k
                        ni(jj) = sum(induntri==jj);
                    end
                    if mixt == 0
                        niini = ni;
                    end

                    %% Weights for the update of the model parameters

                    % The following switch statement is to compute the vector of weights that enter in the
                    % update of the model parameters (remember that in computing the value of the target
                    % function, we use the standard (un-weighted) version of the likelihood functions --
                    % mixture or classification -- ).

                    xmoduntri_unthinned = xmoduntri;
                    induntri_unthinned = induntri;
                    weuntri_unthinned = weuntri;
                    if mixt == 2
                        postprobuntri_unthinned = postprobuntri;
                    end
                    switch wtrim

                        case {0,4}
                            %no observation weighting, therefore:
                            if mixt == 2
                                % for mixture likelihood, the weights are the posterior probabilities
                                weightsuntri = postprobuntri;
                            elseif mixt == 0
                                % for crisp clustering the weights are a vector of ones
                                weightsuntri = repmat(weuntri,1,k);
                            end

                            weightsuntri_unthinned = weightsuntri;
                            
                        case 1
                            % user weights applied to each observation;   the weights are the posterior
                            % probabilities multiplied by the user weights

                            if mixt == 2
                                weightsuntri = postprobuntri .* repmat(weuntri,1,k);
                            elseif mixt == 0
                                weightsuntri =   repmat(weuntri,1,k);
                            end
                            
                            weightsuntri_unthinned = weightsuntri;

                        case 2
                            % weights are the posterior probabilities multiplied by the density

                            %initialize weight for trimming
                            %we = wedef ;

                            if mixt == 2
                                % for mixture likelihood, the weights are the posterior probabilities
                                weightsuntri = postprobuntri;
                            elseif mixt == 0
                                % for crisp clustering the weights are ...
                                weightsuntri =  repmat(weuntri,1,k);
                            end

                            % the next lines are written to assigne a weight to the trimmed observations.

                            %indall_good = vector of length n containing the id 1,...,k of the group the
                            %observation belongs to or "-1" if the observations was trimmed
                            indall = -ones(n,1);
                            indall(qqassigned) = induntri;

                            % indall_good_and_outl  = vector of length n containing the id 1,...,k of the group the
                            %observation belongs to, for trimmed and untrimmed units.
                            % indall_good  = vector of length n containing the id 1,...,k of the group the
                            %observation belongs to, for untrimmed units.
                            indall_good_and_outl = NaN(1,n);
                            indall_good_and_outl(qqassigned)   = induntri;
                            indall_good_and_outl(qqunassigned) = indtri;
                            %Xsort_ll is X sorted in ascending order of loglikelihood, as Xuntri etc.
                            Xsort_ll = ones(n,1+intercept);
                            Xsort_ll(1:length(qqassigned),1+intercept) = Xuntri;
                            Xsort_ll(length(qqassigned)+1:end,1+intercept) = Xtri; 
                            %if it is not possible to  compute wthin in a group
                            %(because there are less than thinning_th obs or
                            %because beta is zero), we are set to the median of the
                            %other weights. To do this we create a k-vector
                            %check_we_groups which contains k zeros, if all the k
                            %groups can be computed by wthin; otherwise it contains
                            %ones in correpondence of the groups which cannot be
                            %computed by wthin.
                            check_we_groups = zeros(k,1);
                            id_obs_in_groups_no_weight = zeros(n,k);
                            for jj=1:k
                                % ijj: indices of untrimmed units in group jj
                                ijj = find(induntri==jj);
                                %ijj_ori = vector of length n containing indices of untrimmed units in group jj
                                ijj_ori = indall == jj; % was find( indall == jj)

                                %ijj_ori_good_and_outl = vector of length n containing indices of trimmed and untrimmed units in group jj
                                ijj_ori_good_and_outl = indall_good_and_outl == jj; % was find( indall_good_and_outl == jj)
                                %ijj_ori_good = vector of length n containing indices of untrimmed units in group jj
                                %ijj_ori_good = indall_good == jj; % was find( indall_good_and_outl == jj)

                                % find indices of trimmed and untrimmed units in group jj
                                ijj_good_and_outl = find(ijj_ori_good_and_outl);


                                % weight vector is updated only if the group has
                                % more than thinning_th observations and if the
                                % beta of the group is not zero

                                %computation of weights for trimmed and non-trimmed units
                                if  numel(ijj_good_and_outl)>thinning_th && nameYY(end,jj) ~= 0
                                    % retention probabilities based on density
                                    % estimated on the component predicted values
                                    % of the previous step. The bernoulli weights
                                    % (the first output argument of wthin, i.e. Wt)
                                    % are not used.

                                    X_jj = Xsort_ll(ijj_good_and_outl,:);
                                    yhat = X_jj*nameYY(:,jj);

                                    [~ , pretain1] = wthin(yhat);

                                    we(ijj_ori_good_and_outl) = pretain1;

                                else
                                    check_we_groups(jj) = 1;
                                    id_obs_in_groups_no_weight(:,jj) = ijj_ori_good_and_outl';
                                end

                                %computation of weights only for non-trimmed units.
                                %It is necessary to compute weights only for
                                %untrimmed units.
                                if  numel(ijj)>thinning_th && nameYY(end,jj) ~= 0
                                    Xuntri_jj = Xuntri(ijj,:);
                                    yhatuntri = Xuntri_jj*nameYY(:,jj);

                                    [~ , pretain2] = wthin(yhatuntri);

                                    we(ijj_ori) = pretain2;

                                    %weights for the parameter estimation step. They
                                    %are computed only for non-trimmed units.
                                    weightsuntri(ijj,jj) = weightsuntri(ijj,jj) .* pretain2;
                                end

                                %                        %computation of weiights only for non-trimmed units
                                %                         if  numel(ijj)>thinning_th && nameYY(end,jj) ~= 0
                                %                             % retention probabilities based on density estimated on the component
                                %                             % predicted values of the previous step. The bernoulli weights
                                %                             % (the first output argument of wthin, i.e. Wt) are not used.
                                %                             Xtri_jj = Xtri(ijj,:);
                                %                             yhattri = Xtri_jj*nameYY(:,jj);
                                %
                                %                             [~ , pretain] = wthin(yhattri);
                                %
                                %                             we(ijj_ori) = pretain;
                                %                             weights(ijj,jj) = weights(ijj,jj) .* pretain;
                                %
                                %                         end


                            end
                            %if in some group it was not possible to run wthin,
                            %then we is fixed equal to the median of the weights
                            %(we) of the other groups. In order to simplify the
                            %code, the median is computed using the weights (we)
                            %different from one, which is the default value. In
                            %this way we underestimate the median weight (we), in
                            %case wthin assigns to some observations we=1
                            if sum(check_we_groups) > 0
                                median_we = nanmedian(we(we~=1)) ;
                                pos_check_we_groups = find(check_we_groups==1);
                                for jj=1:length(pos_check_we_groups)
                                    %if the number of observations in the group not
                                    %analysed with wthin is 0 there is necessity to
                                    %assign the median
                                    n_obs_current_gr = sum(id_obs_in_groups_no_weight(:,pos_check_we_groups(jj)));
                                    if n_obs_current_gr > 0
                                        rep_median_we = repmat(median_we,n_obs_current_gr,1);
                                        we(id_obs_in_groups_no_weight(:,pos_check_we_groups(jj)) == 1) = rep_median_we;
                                    end
                                end
                            end

                            weightsuntri_unthinned = weightsuntri;
                        case 3

                            % weights are the posterior probabilities multiplied by
                            % the bernoulli weights

                            %initialize weight for trimming
                            we = wedef ;

                            % initialize weight vector with posterior probabilities
                            if mixt == 2
                                weightsuntri = postprobuntri;
                            elseif mixt == 0
                                weightsuntri =  repmat(weuntri,1,k);
                            end
                            weightsuntri_unthinned = weightsuntri;
                            %indall = vector of length n containing the id 1,...,k
                            %of the group the observation belongs to or "-1" if the
                            %observations was trimmed
                            %indall = -ones(n,1);
                            %indall(qqassigned) = induntri;
                            ii = 0;


                            for jj=1:k
                                % find indices of units in group jj
                                ijj = find(induntri==jj);
                                %ijj_ori = find(indall == jj);
                                % weight vector is updated only if the group has
                                % more than thinning_th observations abd if the
                                % beta of the group is not zero
                                if  numel(ijj)> thinning_th && nameYY(end,jj) ~= 0
                                    % Bernoulli weights based on density estimated
                                    % on the component predicted values of the
                                    % previous step. The retention probabilities
                                    % (the second output argument of wthin, i.e.
                                    % pretain) are not used.
                                    Xuntri_jj = Xuntri(ijj,:);
                                    yhatuntri = Xuntri_jj*nameYY(:,jj);

                                    [Wt , ~] = wthin(yhatuntri);
        %                            figure;gscatter(Xtri_jj,ytri_jj,Wt)
                                    % the ids of the thinned observations. Values between [1 n_untrimmed]
                                    idWt0 = ijj(Wt == 0);
                                    %update of the we and wetri vector, necessary for doing
                                    %the trimming on all the observations in the
                                    %next step. wetri is n_untrimmed x 1 we is n x 1.
                                    weuntri(idWt0) = 0;
                                    weuntri_unthinned =  weuntri;

        %                            figure;gscatter(X,y,we)

                                    %update of the weights vector, necessary for
                                    %doing the regression parameter estimation on
                                    %the observations not trimmed and not thinned.
                                    %weight has size (n_not_trimmes x k)
                                    weightsuntri(ijj,jj) = weightsuntri(ijj,jj) .* Wt;

                                    % count the thinned observations
                                    nthinned = sum(Wt == 0);
                                    ii = ii + nthinned;
                                    %eliminate thinned observations
                                    ni(jj) = ni(jj) - nthinned;

                                end

                            end

                            we(qqassigned)=weuntri;
                            xmoduntri_unthinned(weuntri_unthinned == 0,:)=[];
                            induntri_unthinned(weuntri_unthinned == 0,:)= [];
                            weightsuntri_unthinned(weuntri_unthinned == 0,:)= [];
                            if mixt == 0
                                   niini =ni;
                            elseif mixt ==0
                                    postprob(we==0) = 0;
                                    postprobuntri_unthinned (weuntri==0) = 0;
                            end           

                    end

                    weightmoduntri_unthinned = [weightsuntri_unthinned, induntri_unthinned ];
                    % initializations of xmodtemp which is a working matrix used
                    % for creating xmod (which contains the results of the second
                    % trimming for all the observations) from xmodjj (which
                    % contains the results of the second trimming for the current
                    % group).
                    xmodtemp    = zeros(n,p+2);
                    % initializations of indxmodtemp which is a working scalar used
                    % to identify the rows of xmodtemp, where to append the
                    % following group results.
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
                    trim2level = [];
                    for iii = not_empty_g
                        jk = jk+1;

                        %check if a group is populated
                        if iii == 1
                            %extract x and y belonging to group iii

                            xmodjuntri_unthinned = xmoduntri_unthinned(xmoduntri_unthinned(:,end)==jk,:);
                            %extract the weights (for beta estimation) of
                            %observations belonging to group iii
                            weightmodjuntri_unthinned = weightmoduntri_unthinned(weightmoduntri_unthinned(:,end) == jk,:);


                            % qqs contains contains the indexes of untrimmed units
                            % (after 2nd level trimming)  for group j
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
                                        RAW = mcd(xmodjuntri_unthinned(:,2:p),'bdp',alpha2,'msg',0);
                                    elseif alpha2 >= 1
                                        RAW = mcd(xmodjuntri_unthinned(:,2:p),'bdp',alpha2/n,'msg',0);
                                    end
                                else
                                    if alpha2 < 1
                                        RAW = mcd(xmodjuntri_unthinned(:,1:p),'bdp',alpha2,'msg',0);
                                    elseif alpha2 >= 1
                                        RAW = mcd(xmodjuntri_unthinned(:,1:p),'bdp',alpha2/n,'msg',0);
                                    end
                                end
                                [~,indmdsor] = sort(RAW.md);

                                if alpha2 < 1
                                    qqs = indmdsor(1:floor(ni(jk)*(1-alpha2)));
                                else
                                    qqs = indmdsor(1:floor(ni(jk) - alpha2));
                                end

                                qqs_trim = setdiff(indmdsor, qqs);
                                %observations trimmed with the 2nd level trimming in original scale
                                trim2level = [trim2level ; xmodj(qqs_trim , size(xmodj,2)-1)];
                            end

                            %% new mixture parameters computed using OLS
                            %x and y belonging to group iii, after second level trimming.
                            xxx = xmodjuntri_unthinned(qqs,1:p);
                            yyy = xmodjuntri_unthinned(qqs,p+1);
                            %dimension of group iii, after second level trimming.
                            ni(jk) = length(yyy);
                            %weights (for beta estimation) of observations
                            %belonging to group iii, after second level trimming.
                            weightmodj_jk = sqrt(weightmodjuntri_unthinned(qqs,jk));
                            %weighted regression for group iii, after second level trimming.
                            breg =  (bsxfun(@times,xxx, weightmodj_jk)) \ (bsxfun(@times,yyy ,weightmodj_jk));
                            %store beta of the current group.
                            nameYY(:,jk) = breg;
                            % Update residuals
                            residuals = yyy-xxx*breg;
                            % Update sigmas through the mean square residuals
                            sigmaini(jk) = sum((residuals .* weightmodj_jk).^2)/(sum((weightmodj_jk).^2));
                            %xmodtemp is a working matrix necessary to concatenate
                            %the results of the second level trimming of the
                            %current group, with all the other groups.
                            xmodtemp((indxmodtemp+1):(indxmodtemp+ni(jk)),:) = xmodjuntri_unthinned(qqs,:);
                            indxmodtemp = indxmodtemp+ni(jk);

                        else

                            xmodjuntri_unthinned = xmoduntri_unthinned(xmoduntri_unthinned(:,end)==jk,:);

                            if alpha2 == 0
                                qqs = [];
                            else
                                qqs = [];
                            end

                            ni(jk) = 0;
                            breg = NaN;
                            nameYY(:,jk) = breg;
                            %xmodtemp is a working matrix necessary to concatenate the results of the second
                            %level trimming of the current group, with all the other groups.
                            xmodtemp((indxmodtemp+1):(indxmodtemp+ni(jk)),:) = xmodjuntri_unthinned(qqs,:);
                            indxmodtemp = indxmodtemp+ni(jk);
                            sigmaini(jk) = NaN;
                            %count the number of times in a group there are enough
                            %observations to compute the sigma
                            count1_ng_eq_k = count1_ng_eq_k + 1;
                        end

                    end

                    sigmaini= restreigen(sigmaini,ni',restrfact,tolrestreigen,userepmat);

                    %for computing the objective function, if a group is emty, beta and sigma are computed as mean of the
                    %other groups. In order to be passed to the next refining step, after having computed the objective function,
                    %they will be set at NaN.
                    for j=1:k
                        if isnan(sigmaini(j))
                            sigmaini(j) = nanmean(sigmaini);
                        end
                        if isnan(nameYY(:,j))
                            nameYY(:,j) = nanmean(nameYY,2);
                        end
                    end
                    xmoduntri_unthinned = xmodtemp(1:indxmodtemp,:);

                    % Stop if two consecutive concentration steps have the same result
                    if indmax == indold
                        break
                    else
                        indold = indmax;
                    end

                    %% Compute the value of the target function
                    obj = 0;
                    not_empty_g = ~( ni <= p + 1 );
                    if mixt == 0

                        jk = 0;
                        for iii = not_empty_g
                            jk = jk+1;
                            % Update weights
                            %niini(jk) = ni(jk);%bug
                            if iii ==1
                                yj = xmoduntri_unthinned(xmoduntri_unthinned(:,end) == jk,end-1);
                                Xj = xmoduntri_unthinned(xmoduntri_unthinned(:,end) == jk,1:end-2);

                                %The following line should be executed at the end
                                %of the concentration steps loop. It is executed at
                                %each step because of the break above, which is
                                %executed if two consecutive concentration steps
                                %have the same result

                                obj = obj + niini(jk)*log(niini(jk)/sum(niini)) +...
                                    sum(logmvnpdfFS(yj-Xj*nameYY(:,jk),0,(sigmaini(jk))));
                            else
                                %if a group is missing, we do not compute its
                                %objective function contribution.
                            end
                        end

                        %                 nameYY(:,~not_empty_g) = NaN;
                        %                 sigmaini(~not_empty_g) = NaN;

                    elseif mixt == 2
                        %the following command should be executed at the end of
                        %the for loop of the concentration steps. However here
                        %it is executed in all steps, because of the above break from the
                        %loop, which is executed if two consecutive
                        %concentration steps have the same result

                        log_lh=NaN(size(xmoduntri_unthinned,1),size(not_empty_g,2));

                        %log_lh = [];
                        jk = 0;
                        for iii = not_empty_g
                            jk = jk+1;
                            if iii ==1
                                log_lh(:,jk) = ...
                                    log(niini(jk)/sum(niini)) + (logmvnpdfFS(...
                                    xmoduntri_unthinned(:,end-1) - ...
                                    xmoduntri_unthinned(:,1:(size(xmoduntri_unthinned,2)-2)) * ...
                                    nameYY(:,jk),0,(sigmaini(jk)) ) );
                            else
                                %if a groupis missing, we do not compute the objective
                                %function for it.
                                log_lh(:,jk) = NaN(length(xmoduntri_unthinned),1);
                            end
                        end

                        group_missing = sum(isnan(log_lh),1)>0;
                        log_lh(:,group_missing)=[];
                        obj = estepFS(log_lh);
                        %                 if ~isempty(group_missing)
                        %                     nameYY(:,group_missing) = NaN;
                        %                     sigmaini(group_missing) = NaN;
                        %                 end
                    end
                end
            end % End of concentration steps

            %% Change the 'optimal' target value and 'optimal' parameters
            % This is done if an increase in the target value is achieved
            %this check has to be commented in order to estimate the effect of eps_beta
            if sum(not_empty_g ) == k
                if sum(sum(isnan(nameYY))) == 0
                    if (obj >= vopt)
                        check_obj_reached         = 1;
                        vopt                                   = obj;
                        bopt                                   = nameYY;
                        numopt                              = niini;
                        sigmaopt                           = sigmaini;
                        weopt                                = we;
                        trim1levelopt                     = trim1level;
                        trim2levelopt                     = trim2level;
                        indmaxopt                         = indmax;
                        not_empty_gopt              = not_empty_g;
                        count1_ng_lt_kopt          = count1_ng_lt_k;
                        count2_ng_lt_kopt          = count2_ng_lt_k;
                        count1_ng_eq_kopt       = count1_ng_eq_k;
                        count2_ng_eq_kopt       = count2_ng_eq_k;
                        extra_inisubsopt             = nselected - nsampdef;
                        selj_goodopt                  = selj_good_groups;
                        selj_elimopt                    = selj_elim_groups(1:count_elim,:);
                        selj_allopt                       = [selj_elim_groups(1:count_elim,:); selj_good_groups];

                        if mixt ==2
                            postprobopt = postprob;
                        end

                    end
                end
            end
        end
    end % end of loop over the nsamp subsets
    if check_obj_reached == 0
        restrfact = restrfact*2;
        disp(['--------- restrfact has been duplicated to ' num2str(restrfact)]);
    end
end
if count < nsamp
    out = struct;
else
    
    
    %% Generate plots
    
    if plots
        %count the number of obs in each group without the trimmed and the thinned
        trim1level_01opt = ones(n,1);
        trim2level_01opt = ones(n,1);
        trim1level_01opt(trim1levelopt) = 0;
        trim2level_01opt(trim2levelopt) = 0;

        % The following plots are for the bi-variate case (i.e. v=1)
        if v < 2
            
            % initialize figure
            fh = figure('Name','TclustReg plot','NumberTitle','off','Visible','on');
            gca(fh);
            hold on;
            xlabel('X');
            ylabel('y');
            title({ ['TclustReg clustering: ' 'mixt=' num2str(mixt) ' - wtrim=' num2str(wtrim)] },'Fontsize',12);
            
            jk = 0;
            for iii = not_empty_gopt
                jk = jk+1;
                if iii>0
                    group_label = ['Group ' num2str(jk)];
                    
                    % plot of the good units allocated to the current group.
                    % Indices are taken after the second level trimming.
                    % Trimmed points are not plotted by group.
                    if wtrim ==3
                        ucg = find(indmaxopt==jk & weopt == 1 & trim1level_01opt ==1 & trim2level_01opt == 1);
                    else
                        ucg = find(indmaxopt==jk &  trim1level_01opt ==1 & trim2level_01opt == 1);
                    end

                    % misteriously text does not show a legend. This is why
                    % we add a (ficticious) plot instruction with white symbols 
                    plot(X(ucg,end),y(ucg),'.w','DisplayName',[group_label ' (' num2str(length(ucg)) ' units)']);
                    text(X(ucg,end),y(ucg),num2str(jk*ones(length(ucg),1)),...
                        'DisplayName',[group_label ' (' num2str(length(ucg)) ' units)'] , ...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle',...
                        'Color',clrdef(jk));
                    
                    % plot regression lines
                    vv = [min(X(:,end)) max(X(:,end))];
                    if intercept==1
                        plot(vv,bopt(1,jk)+bopt(2,jk)*vv,...
                            'DisplayName',[group_label ' fit' ],...
                            'Color',clrdef(jk));
                    elseif intercept==0
                        plot(vv,bopt(:,jk)*vv,...
                            'DisplayName',[group_label ' fit' ],...
                            'Color',clrdef(jk));
                    end
                    
                    %plot the thinned units
                    if wtrim == 3
                        %                text(X(ucg,end),y(ucg),num2str(0*ones(length(ucg),1)),...
                        %                         'DisplayName','Thinned units' , ...
                        %                         'HorizontalAlignment','center',...
                        %                         'VerticalAlignment','middle',...
                        %                         'Color',clrdef(k+1));
                        ucg = find(indmaxopt==jk & weopt == 0);
                        plot(X(ucg,end),y(ucg),symdef(jk),'color',clrdef(k+1),...
                            'DisplayName',['Group ' num2str(jk) ': thinned units (' num2str(length(ucg)) ')' ]);
                    end
                    
                end
            end
            
            % Plot the outliers (trimmed points)
            plot(X(trim1levelopt,end),y(trim1levelopt),'o','color','r','MarkerSize',8,...
                'DisplayName',['Trimmed units (' num2str(length(y(trim1levelopt))) ')']);
            
            % Second level trimming points
            xxx0_all = X(trim2levelopt,end);
            yyy0_all = y(trim2levelopt);
            plot(xxx0_all,yyy0_all,'*','color','c',...
                'DisplayName',['L2 trimmed units (' num2str(length(yyy0_all)) ')']);
            
            if wtrim == 4
                % in case of tandem thinning, plot the thinned points
                plot(Xori(~Wt4,end),yori(~Wt4),symdef(k+1),'color',clrdef(k+1),...
                    'DisplayName',['Thinned units (' num2str(length(Wt4) - sum(Wt4)) ')']);
            end
            
            % Position the legends and make them clickable. For some reason
            % clickableMultiLegend does not set properly the FontSize: to be fixed. 
            lh=legend('show');
            legstr = get(lh,'String');
            clickableMultiLegend(legstr,'FontSize',14,'Location','northwest');
            
            axis('manual');
            
            % control of the axis limits
            xmin = min(X(:,end)); xmax = max(X(:,end));
            ymin = min(y); ymax = max(y);
            deltax = (xmax - xmin) / 10;
            deltay = (ymax - ymin) / 10;
            
            xlim([xmin-deltax,xmax+deltax]);
            ylim([ymin-deltay,ymax+deltay]);
            
        else
            
            % In this case p > 2. A standard spmplot is used.
            
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
    
    out                                   = struct;
    out.bopt                          = bopt;
    out.sigmaopt_0              = sigmaopt_0;
    out.sigmaopt_cons        = sigmaopt_cons;
    out.sigmaopt_pison       = sigmaopt_pison;
    out.numopt                      = numopt;
    out.vopt                            = vopt;
    out.asig_obs_to_group  = indmaxopt;
    out.weopt                         = weopt;
    out.trim1levelopt              = trim1level_01opt;
    out.trim2levelopt              = trim2level_01opt;
    out.count1_ng_lt_k          = count1_ng_lt_kopt;
    out.count2_ng_lt_k          = count2_ng_lt_kopt;
    out.count1_ng_eq_k       = count1_ng_eq_kopt;
    out.count2_ng_eq_k       = count2_ng_eq_kopt;
    out.extra_inisubs             = extra_inisubsopt;
    out.selj_good                  = selj_goodopt;
    out.selj_elim                    = selj_elimopt;
    out.selj_all                       = selj_allopt;
    if mixt == 2
        out.postprobopt          = postprobopt;
    end
    out.restrfact                     = restrfact;
    
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

