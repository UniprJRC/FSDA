function [out, varargout] = tclustreg(y,X,k,restrfact,alpha1,alpha2,varargin)
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
%equalweights : cluster weights in the concentration and assignment steps.
%               Logical. A logical value specifying whether cluster weights
%               shall be considered in the concentration, assignment steps
%               and computation of the likelihood.
%               if equalweights = true we are (ideally) assuming equally
%               sized groups by maximizing
%                 Example - 'equalweights',true
%                 Data Types - Logical
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
%      we: Vector of observation weights. Vector. A vector of size n-by-1
%          containing the weights to apply to each observation. Default
%          value: vector of ones.
%            Example - 'we',[0.2 0.2 0.2 0.2 0.2]
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
%   out.weopt       = $n$ vector containing the weigths of each
%                           observation, i.e. its contribution to the
%                           estimates.
%   out.asig_obs_to_group   = $n$ vector containing the  cluster assigments
%                       of all n observations (trimmed observations
%                       excluded).
% out.asig_obs_to_group_before_tr = $n$ vector containing the  cluster assigments
%                       of all n observations (trimmed observations
%                       included).
%   out.trim1levelopt       = $n$ vector containing the 1st level trimmed units (0) and
%                       1st level untrimmed (1) units.
%   out.trim2levelopt       = $n$ vector containing the 2nd level trimmed units (0) and
%                           2nd level untrimmed (1) units.
%   out.postprob       = $n$ vector containing the final posterior probability
%    out.C       = initial subsets extracted
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
    %% tclustreg of 'X data' (Gordaliza, Garcia-Escudero & Mayo-Iscar, 2013).
    % The dataset presents two parallel components without contamination.
    X  = load('X.txt');
    y1 = X(:,end);
    X1 = X(:,1:end-1);

    k = 2 ; restrfact = 5; alpha1 = 0.05 ; alpha2 = 0.01;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2);

    k = 2 ; restrfact = 2; alpha1 = 0.05 ; alpha2 = 0.01;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'mixt',2);

    % Comparison with robust linear grouping
    figure;
    out = rlga(X,k,alpha2);
    title('robust linear grouping on ''X data'' ');
%}

%{
    clear all; close all;
    load fishery;
    X = fishery.data;
    % some jittering is necessary because duplicated units are not treated:
    % this needs to be addressed
    X = X + 10^(-8) * abs(randn(677,2));

    %tclustreg on fishery data
    y1 = X(:,end);
    X1 = X(:,1:end-1);
    k = 3 ; restrfact = 50; alpha1 = 0.04 ; alpha2 = 0.01;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',0);

    %lga and rlga on fishery data
    figure('name','RLGA');
    out=lga(X,3);
    clickableMultiLegend('1','2','3','data1','data2','data3');
    axis manual;

    alpha = 0.95;
    figure('name','LGA');
    out=rlga(X,3,1-alpha);
    clickableMultiLegend('0','1','2','3','data1','data2','data3');
    axis manual;
%}

%{
    clear all; close all;
    load fishery;
    X=fishery.data;
    % some jittering is necessary because duplicated units are not treated
    % in tclustreg: this needs to be addressed
    X = X + 10^(-8) * abs(randn(677,2));
    y1 = X(:,end);
    X1 = X(:,1:end-1);
    
    % some arbitrary weights for the units
    we = sqrt(X1)/sum(sqrt(X1));
    
    % tclustreg required parameters 
    k = 3 ; restrfact = 50; alpha1 = 0.04 ; alpha2 = 0.01;

    % now tclust is run on each combination of mixt and wtrim options

    disp('mixt = 0; wtrim = 0;');
    disp('standard tclustreg, with classification likelihood and without thinning' );
    mixt = 0; wtrim = 0;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 0;');
    disp('mixture likelihood, no thinning' );
    mixt = 2; wtrim = 0;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 0; wtrim = 1;');
    disp('classification likelihood, thinning based on user meights' );
    mixt = 0; wtrim = 1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'we',we,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 1;');
    disp('mixture likelihood, thinning based on user meights' );
    mixt = 2; wtrim = 1;
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'we',we,'wtrim',wtrim);

    disp('mixt = 0; wtrim = 2;');
    disp('classification likelihood, thinning based on retention probabilities' );
    mixt = 0; wtrim = 2; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 2;');
    disp('mixture likelihood, thinning based on retention probabilities' );
    mixt = 2; wtrim = 2; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 0; wtrim = 3;');
    disp('classification likelihood, thinning based on bernoulli weights' );
    mixt = 0; wtrim = 3; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 3;');
    disp('mixture likelihood, thinning based on bernoulli weights' );
    mixt = 2; wtrim = 3; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 0; wtrim = 4;');
    disp('classification likelihood, tandem thinning based on bernoulli weights' );
    mixt = 0; wtrim = 4; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);

    disp('mixt = 2; wtrim = 4;');
    disp('mixture likelihood, tandem thinning based on bernoulli weights' );
    mixt = 2; wtrim = 4; we = [];
    out = tclustreg(y1,X1,k,restrfact,alpha1,alpha2,'intercept',0,'mixt',mixt,'wtrim',wtrim);
%}

%{
    %% Generate mixture of regression using MixSimReg, with an average
    % overlapping at centroids = 0.01. Use all default options.
    rng(372,'twister');
    p=3;
    k=2;
    Q=MixSimreg(k,p,'BarOmega',0.001);
    n=500;
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);

    % plot the dataset
    yXplot(y,X);

    % run tclustreg
    out=tclustreg(y,X,k,50,0.01,0.01,'intercept',1);
%}

%{
    % Generate mixture of regression using MixSimReg, with an average
    % overlapping at centroids =0.01. 
    p=3;
    k=2;
    Q=MixSimreg(k,p,'BarOmega',0.001);
    n=200;
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);

    % Generate the elemental subsets used in tclustreg once and for all.
    nsamp  = 100;
    ncomb  = bc(n,p);
    method = [10*ones(n/2,1); ones(n/2,1)]; % weighted sampling using weights in "method"
    msg    = 0;
    C      = subsets(nsamp, n, p, ncomb, msg, method);

    % tclustreg using samples in C
    out=tclustreg(y,X,k,50,0.01,0.01,'nsamp',C);
%}

%%  computational issues to be addressed in future releases

% FSDA function wthin uses the MATLAB function ksdensity. The calls to
% ksdensity have been optimized. The only possibility to further reduce
% time execution is to replace ksdensity with a more efficient kernel
% density estimator.

% The weighted version of tclustreg requires weighted sampling. This is
% now implemented in randsampleFS. A computaionally more efficient
% algorithm, based on a binary tree approach introduced by
%      Wong, C.K. and M.C. Easton (1980) "An Efficient Method for Weighted
%      Sampling Without Replacement", SIAM Journal of Computing,
%      9(1):111-113.
% is provided by recent releases of the MATLAB function datasample.
% Unfortunately this function spends most of the self running time useless
% parameter checking. To copy the function in the FSDA folder
% FSDA/combinatorial, possibly rename it, and remove the option parameters
% checks, is not sufficient, as datasample relies on a mex file wswor which
% is platform dependent. The issue is usually referred to as "code signature".

% In the plots, the use of text to highlight the groups with their index is
% terribly slow (more than 10 seconds to generate a scatter of 7000 units.
% ClickableMultiLegend and legend are also slow.

% FSDA function restreigen could be improved. In some cases it is one of
% the most expensive functions.

%% REMARK: trimming vs thinning
%
% - trimming is the key feature of TCLUST, giving robustness to the model.
% - thinning is a new denoising feature introduced to mitigate the distorting effect of
%   very dense data areas. It is implemented via observation weighting.
%
% For the sake of code readability, the relevant sections of the code will
% be identified with a "TRIMMING" or "THINNING" tag.


%% initializations

warning('off');

% Groups with less than thinning_th units are not considered for thinning
thinning_th = 50;

% tolerance for restriction factor
tolrestreigen = 1e-08;

% initial objective function value (optimized during the random starts)
vopt = -1e+20;

% this is just for rotating colors in the plots
clrdef = 'bkmgrcbkmgrcbkmgrcbkmgrcbkmgrcbkmgrcbkmgrc';
symdef = '+*sd^v><pho*';

% repmat from Release 8.2 is faster than bsxfun
if verLessThan('matlab','8.2.0') ==1
    userepmat=0;
else
    userepmat=1;
end

%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);


seq=1:n;

% THINNING: check if wtrim user option is equal to 4. In such a case,
% thinning is applied just at the very beginning and on the full dataset.
% This has to preceed the setting of number of random samples nsamp.
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
                X    = X(Wt4,:);
                y    = y(Wt4);
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

% nsamp option: if it is not set by the user, it has to be properly initialized
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

% default number of concentration starts
Kstepsdef  = 10;

%default value for wtrim
wtrimdef = 0;

%default value for we
wedef = ones(n,1);

%default model (mixture or classification likelihood)
mixtdef = 0;


%seqk = sequence from 1 to number of groups
seqk=1:k;

%% User options

options = struct('intercept',1,'mixt',mixtdef,...
    'nsamp',nsampdef,'Ksteps',Kstepsdef,...
    'startv1',startv1def,'we',wedef,'wtrim',wtrimdef,...
    'msg',0,'plots',1,'equalweights',1);

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

% Number of subsets to extract or matrix containing the subsets
nsamp = options.nsamp;

% Concentration steps
Ksteps = options.Ksteps;

% equalweights
equalweights = options.equalweights;


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

% % First level trimming
% if alpha1<1
%     notrim = floor(n*(1-alpha1));
% else
%     notrim = n-floor(alpha1);
% end
% trimm = n-notrim;



%% Combinatorial part to extract the subsamples (if not already supplied by the user)

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
        niinistart=repmat(floor(notrim/k),k,1);
    end
    % the number of initial subsets to be generated is nsamp*oversamp.
    % The input parameter nsamp is multiplied by a factor (oversamp) in
    % order to increase the number of generated subsets because in some subsets the regression lines
    % could be very closed one to the other and therefore have to be
    % substituted with other subsets.
    for ns =1:nsamp
        % Matlab function datasample adopts an efficient method for
        % weighted sampling based on binary trees. However, it consumes a
        % lot of of time in parameter checking. This is way we use our
        % function randsampleFS.
        %C(ns,:) = datasample(1:n,initial_subs_size,'weights',we,'Replace',false);
        C(ns,:) = randsampleFS(n,initial_subs_size,we);
    end
    %nselected is set equal to the number of subsets:
    % - nselected = nsamp, if nsamp is a scalar;
    % - nselected = size(nsamp,1), if nsamp is a matrix (containing the initial subsets)
    nselected  = length(C);
    
end

%initialization of iter: index which is used to extract the subset in the C
%matrix. Iter takes value in the interval 1:nselected, where nselected is
%the number of rows of C, i.e. nsamp if nsamp is a scalar, or
%size(nsamp,1) if nsamp is a matrix .
iter                = 0;

%this while cycle is necessary when at the end of an entire loop on the
%extracted subsets (while iter < nselected), the best subset has not been
%identified (i.e. the best subset does not exist) or if it has less than k
%groups. If it happens another loop on other extracted subsets (the
%following nselected rows in C matrix) is repeated.

    
    %% Initialize structures
    
    ll                  = zeros(n,k);
    sigmaini            = ones(1,k);
    %initialize the output of tclustreg, i.e. the parameters obtained in the best (optimum) subset
    bopt             = zeros(p,k);
    %numopt = (kX1) vector containing the number of observations in each of
    %the k groups in the optimal subset.
    numopt              = zeros(1,k);
    weopt               = ones(n,1);
    %trim1levelopt = (nx1) vector of 0 and 1 indicating respectively units
    %trimmed and not trimmed in the first trimming level in the optimal
    %subset
    trim1levelopt       = ones(n,1);
    %trim2levelopt = (nx1) vector of 0 and 1 indicating respectively units
    %trimmed and not trimmed in the second trimming level in the optimal
    %subset
    trim2levelopt       = ones(n,1);
    %indmaxopt = (nx1) vector of 1, ... , k, indicating the assignement of
    %all the n observations (trimmed observations included) to one of the k
    %groups  in the optimal subset
    indmaxopt           = zeros(n,1);
    if mixt ==2
        postprobopt     = zeros(n,k);
    end
    
    %%  Random starts
    while iter < nselected
        
        % lessthankgroups will be equal 1 if for a particlar subset less
        % than k groups are found
        lessthankgroups=0;
        
        %iter =  iteration number  of the loop on the subsets, including:
        % -subsets which are refused because the regression lines are too
        % closed one to the other; -an entire loop of subsets which does
        % not lead to any best subset or leads to a best subset with less
        % than k groups.
        iter  = iter+1;
        
        
        
        if msg == 1
            disp(['Iteration ' num2str(iter)])
        end
        %if startv1 = 1, covariance matrices is estimated in each initial subset composed by  (v+1)*k units
        if startv1
            % Beta = matrix of beta coefficients (j-col is referred to j-th
            % group)
            Beta = zeros(p,k);
            % in order to fix the seed decomment the following command
            % rng(1234);
            randk=rand(k,1);
            %the number of observations in each group is randomly assigned
            if alpha1<1
                niini=floor(fix(n*(1-alpha1))*randk/nansum(randk));
            else
                niini=floor(fix(n - floor(alpha1))*randk/nansum(randk));
            end
            %if one or more groups have been initialized with 0 observations, we repete
            %the initialization until all the k groups are populated.
            while sum(niini == 0) >0
                randk=rand(k,1);
                % Initialize niini with with random numbers from uniform
                if alpha1<1
                    niini=floor(fix(n*(1-alpha1))*randk/sum(randk));
                else
                    niini=floor(fix(n -floor(alpha1))*randk/sum(randk));
                end
            end
            %given the current subsets extracted observations (C(iter,:)),
            %X and y are identified for each group and the
            %residuals and covariance (sigmaini) are computed
            for j = 1:k
                ilow   = (j-1)*(p)+1;
                iup    = j*(p);
                index  = C(iter,:);
                selj   = index(ilow:iup);
                Xb     = X(selj,:);
                yb     = y(selj,:);
                niini(j)  = length(yb);
                Beta(:,j) = Xb\yb;
                %Update residuals
                residuals = yb-Xb*Beta(:,j);
                % Update sigmas through the mean square residuals
                if size(selj,2) > p
                    sigmaini(j) = nansum(residuals.^2)/niini(j);
                else
                    if length(yb)==1
                        sigmaini(j)= var(y);
                    else
                        sigmaini(j) =var(yb);
                    end
                end
            end
            sigmaini= restreigen(sigmaini,niini,restrfact,tolrestreigen,userepmat);
        else
            %given the current subsets extracted observations (C(iter,:)),
            %identification of X and y for each group, computation of
            %residuals and sigma
            % initialization of niini with equal proportions
  
            % extract a subset of size v
            index  = C(iter,:);
            Xb = X(index,:);
            yb = y(index);
            
            for j=1:k
                Xbj = Xb((1+(j-1)*p):(j*p),:);
                ybj = yb((1+(j-1)*p):(j*p));
                Beta(:,j) = Xbj\ybj;
                niinistart(j) =  length(ybj);
            end
            niini=niinistart;
            sigmaini=ones(1,k);
        end
        
        
            
            
            %% a) Log-likelihood based on all  n observations
            % Given starting values of parameters (coming from a subset)
            if equalweights == 1
                for jj = 1:k
                    ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*Beta(:,jj),0,sigmaini(jj));
                end
            else
                for jj = 1:k
                    ll(:,jj) = log((niini(jj)/sum(niini))) + logmvnpdfFS(y-X*Beta(:,jj),0,sigmaini(jj));
                end
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
                % indmax assigns each observation to the group with the largest posterior probability
                [~,indmax]= max(postprob,[],2);
                
            elseif mixt == 0
                postprob = ones(n,k);
                [disc,indmax] = max(ll,[],2);
            end
            
            
            %% b) THINNING ONCE AND FOR ALL
            % The following switch statement is to compute the vector of weights that enter in the
            % update of the model parameters (remember that in computing the value of the target
            % function, we use the standard (un-weighted) version of the likelihood functions --
            % mixture or classification -- ).
            
            % pretain is the same for all units
            w4trim=ones(n,1);
            
            switch wtrim
                
                case 0
                    %no observation weighting,
                    
                    
                case 1
                    % user weights applied to each observation;   the weights are the posterior
                    % probabilities multiplied by the user weights
                    
                    % pretain is the same for all units
                    w4trim=we;
                case 2
                    % weights are the posterior probabilities multiplied by the density
                    % the next lines are written to assigne a weight to the trimmed observations.
                    check_we_groups = zeros(k,1);
                    id_obs_in_groups_no_weight = zeros(n,k);
                    for jj=1:k
                        % Boolean index of units forming group j
                        groupj=indmax==jj;
                        % ijj: indices of units in group jj
                        ijj = find(groupj);
                        
                        
                        % weight vector is updated only if the group has
                        % more than thinning_th observations and if the
                        % beta of the group is not zero
                        
                        %computation of weights for trimmed and non-trimmed units
                        %????????????????????? TO CHECL
                        if  length(ijj)>thinning_th
                            % retention probabilities based on density
                            % estimated on the component predicted values
                            % of the previous step. The bernoulli weights
                            % (the first output argument of wthin, i.e. Wt)
                            % are not used.
                            
                            Xj = X(ijj,:);
                            yhat = Xj*Beta(:,jj);
                            
                            [~ , pretain] = wthin(yhat);
                            w4trim(ijj) = pretain;
                            
                        else
  %                          error('Qua Non ho capito')
                            check_we_groups(jj) = 1;
                            id_obs_in_groups_no_weight(:,jj) = groupj;
                        end
                        
                    end
                    %if in some group it was not possible to run wthin,
                    %then we is fixed equal to the median of the weights
                    %(we) of the other groups. In order to simplify the
                    %code, the median is computed using the weights (we)
                    %different from one, which is the default value. In
                    %this way we underestimate the median weight (we), in
                    %case wthin assigns to some observations we=1
                    if sum(check_we_groups) > 0
                        median_we = nanmedian(w4trim(w4trim~=1)) ;
                        pos_check_we_groups = find(check_we_groups==1);
                        for jj=1:length(pos_check_we_groups)
                            %if the number of observations in the group not
                            %analysed with wthin is 0 there is necessity to
                            %assign the median
                            n_obs_current_gr = sum(id_obs_in_groups_no_weight(:,pos_check_we_groups(jj)));
                            if n_obs_current_gr > 0
                                rep_median_we = repmat(median_we,n_obs_current_gr,1);
                                w4trim(id_obs_in_groups_no_weight(:,pos_check_we_groups(jj)) == 1) = rep_median_we;
                            end
                        end
                    end
                    
                    
                case 3
                    
                    % weights are the posterior probabilities multiplied by
                    % the bernoulli weights
                    
                    ii = 0;
                    
                    for jj=1:k
                        % find indices of units in group jj
                        ijj = find(indmax==jj);
                        %ijj_ori = find(indall == jj);
                        % weight vector is updated only if the group has
                        % more than thinning_th observations abd if the
                        % beta of the group is not zero
                        if  numel(ijj)> thinning_th
                            % Bernoulli weights based on density estimated
                            % on the component predicted values of the
                            % previous step. The retention probabilities
                            % (the second output argument of wthin, i.e.
                            % pretain) are not used.
                            Xj = X(ijj,:);
                            yhat = Xj*Beta(:,jj);
                            
                            [Wt , ~] = wthin(yhat);
                            
                            w4trim(ijj)=Wt;
                            
                            % count the thinned observations
                            nthinned = sum(Wt == 0);
                            ii = ii + nthinned;
                            %eliminate thinned observations
                            niini(jj) = niini(jj) - nthinned;
                        end
                        
                    end
                    
            end
            
            % Mean of the weights must be 1
            % These weights will be used in first and second level trimming
            w4trim=w4trim/mean(w4trim);
            
            % Now multiply posterior probabilities (n-by-k matrix) by weights (n-by-1 vector)
            w4beta=bsxfun(@times,postprob,w4trim);
            
            
%             % Note that w4trim is a vector while w4beta is a matrix
%             
%             not_empty_g = ~( niini <= p + 1 );
%             
%             %count number of times the number of groups is lt k
%             if sum(not_empty_g) == k
%                 n_refsteps_eq_k_gr = n_refsteps_eq_k_gr + 1;
%             else
%                 n_refsteps_lt_k_gr = n_refsteps_lt_k_gr + 1;
%             end
            
            
            
            %% CONCENTRATION STEPS
            
            indold = zeros(n,1)-1;
            for cstep = 1:Ksteps
                
                %not_empty_g = zeros(1,k);
                if all(isnan(Beta(:)))
                    break
                else
                    
                    
                    %% 1) Log-likelihood (inside concentration steps)
                    % Discriminant functions for the assignments
                    if equalweights == 1
                        for jj = 1:k
                            ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*Beta(:,jj),0,sigmaini(jj));
                        end
                    else
                        for jj = 1:k
                            ll(:,jj) = log((niini(jj)/sum(niini))) + logmvnpdfFS(y-X*Beta(:,jj),0,sigmaini(jj));
                        end
                    end
                    
                    
                    %In the case of crisp assignement we compute the maximum value of the log-likelihood
                    %(disc) among the k groups, as in point 2.1 of appendix of Garcia Escudero et al.
                    %(2010). In the case of mixture modelling we need an estep to compute the log of the sum
                    %of the log-likelihood (disc) and the normalized log-likelihood (postprob).
                    if mixt == 2
                        [~,postprob,disc] = estepFS(ll);
                        % indmax assigns each observation to the group with the largest posterior probability
                        [~,indmax]= max(postprob,[],2);
                        
                    elseif mixt == 0
                        [disc,indmax] = max(ll,[],2);
                    end
                    
                    % Sort the n likelihood contributions and save in qq the largest n*(1-alpha) likelihood
                    % contributions
                    [~,qq] = sort(disc,'ascend');
                    
                    
                    %% 2) First level trimming (inside concentration steps)
                    % Order the weights according to qq
                    w4trimordered=w4trim(qq);
                    
                    % Find cumlative sum of weights
                    cumsumww = cumsum(w4trimordered);
                    
                    % qqunassigned_small is a n-by-1 Boolean vector which
                    % contains true in the first k positions (units which
                    % have to be trimmed)
                    qqunassigned_small = cumsumww < n*alpha1;
                    
                    % qqunassigned = indexes of units subject to first
                    % level trimming
                    qqunassigned = qq(qqunassigned_small);
                    
                    %indmax_before_tr should be saved in order to be able
                    %to understand which groups the two trimmings affect more
                    indmax_before_tr = indmax;
                    
                    % indmax of the units subject to first level trimming
                    % is set equal to 0
                    indmax(qqunassigned)=0;
                    
                    %% 3) Second level of trimming (inside concentration steps)
                    % FS or MCD are used to find units to trim
                    % Using untrimmed units find beta coefficients and
                    % sigma2 using weighted regression (the weights are
                    % based on thinning fixed once and for all before the
                    % concentraion steps)
                    
                    for jj=1:k
                        % find indices of units belonging to groupj
                        groupjind = find(indmax==jj);
                        
                        %check if a group is populated
                        if length(groupjind) >p+2
                            
                            Xjnointercept  = X(groupjind,intercept+1:end);
                            Xj=X(groupjind,:);
                            
                            yj  = y(groupjind);
                            w4betaj=w4beta(groupjind,jj);
                            w4betaj=w4betaj/mean(w4betaj);
                            
                            njj=size(Xj,1);
                            
                            
                            % hjj = number of units to retan for group j
                            % This is simply h referred to group j
                            hj = floor(njj*(1-alpha2));
                            
                            
                            %The MCD is applied only when v=1, because in this case it is faster
                            %than the FS.
                            if v==1
                                
                                [~,REW]      = mcd(Xjnointercept,'bdp',alpha2,'msg',0);
                                %R2016a has introduced robustcov, which could be used here as below.
                                %Remember however that mcd returns the squared distances, i.e. RAW.md = mah.^2.
                                %[~,~,mah,~,~] = robustcov(inliers,'Method','fmcd','NumTrials',nsampmcd,'OutlierFraction',alpha2b,'BiasCorrection',1); %
                                %plot(1:ni(jk),mah.^2,'-r',1:ni(jk),RAW.md,'-b');
                                %close;
                                [~,indmdsor] = sort(REW.md);
                                % Trimmed units (after second level
                                % trimming)
                                % indmax for second level trimmed units is
                                % set to -1
                                trimj=indmdsor(hj+1:end);
                                
                            else
                                % Lines below are to find the second trimming points id_trim with
                                % the Forward Search rather than the MCD, using function FSMmmd
                                % Function FSMbsb is consderably faster than mcd when v>1.
                                % BBsel contains a NAN for the units
                                % not belonging to subset in step hj
                                [~,BBsel]=FSMbsb(Xjnointercept,0,'bsbsteps',hj,'init',hj,'nocheck',1,'msg',0);
                                trimj=seq(isnan(BBsel));
                            end
                            
                            % indmax for second level trimmed units is
                            % equal to -1
                            indmax(groupjind(trimj))=-1;
                            
                            % Redefine Xj and yj ater removing 2nd level
                            % trimmed units
                            Xj(trimj,:)=[];
                            yj(trimj)=[];
                            w4betaj(trimj)=[];
                            nj=length(yj);
                            
                            % Number of units belonging to group j after
                            % second level trimming
                            niini(jj)=nj;
                            
                            
                            % Find beta coefficients and sigma2 using weighted
                            % regression
                            
                            %weights (for beta estimation) of observations
                            %belonging to group iii, after second level trimming.
                            sqweights = sqrt(w4betaj);
                            
                            Xw = bsxfun(@times, Xj, sqweights);
                            yw = yj .* sqweights;
                            
                            % Estimate of beta from (re)weighted regression (RWLS)
                            breg = Xw\yw;
                            Beta(:,jj)=breg;
                            
                            % Find estimate of sigma2 after weighted regression
                            res2=(yw-Xw*breg).^2;
                            sigma2=sum(res2)/(nj-size(Xw,2));
                            
                            sigmaini(jj) = sigma2;
                            
                        else
                            lessthankgroups=1;
                            break
                        end
                    end % loop on groups
                    
                    % get out of loop if you find less than k groups
                    if lessthankgroups==1
                        Beta=NaN;
                        break
                    end
                    
                    % ?????????????????????????
                    %for computing the objective function, if a group is empty, beta and sigma are
                    %computed as mean of the other groups. In order to be passed to the next
                    %refining step, after having computed the objective function, they will be set
                    %at NaN.
                    for j=1:k
                        if isnan(sigmaini(j))
                            
                            disp('sigmaininan')
                            error('nansigma')
                        end
                        if isnan(Beta(:,j))
                            
                            disp('betanan')
                            error('betanan')
                        end
                    end
                    
                    % Apply eigenvalue restrictions on sigmas from
                    % regression
                    sigmaini= restreigen(sigmaini,niini,restrfact,tolrestreigen,userepmat);
                    
                    % Stop if two consecutive concentration steps have the same result
                    if indmax == indold
                        break
                    else
                        indold = indmax;
                    end
                    
                    %% 4) Computation of the value of the target function
                    obj = 0;
                    not_empty_g = seqk(~( niini <= p + 1 ));
                    if mixt == 0
                        
                        for jj = not_empty_g
                            
                            yj = y(indmax==jj);
                            Xj = X(indmax==jj,:);
                            
                            %The following line should be executed at the end
                            %of the concentration steps loop. It is executed at
                            %each step because of the break above, which is
                            %executed if two consecutive concentration steps
                            %have the same result
                            if equalweights ==1
                                obj = obj + log(1/k) +...
                                    sum(logmvnpdfFS(yj-Xj*Beta(:,jj),0,sigmaini(jj)));
                            else
                                obj = obj + niini(jj)*log(niini(jj)/sum(niini)) +...
                                    sum(logmvnpdfFS(yj-Xj*Beta(:,jj),0,sigmaini(jj)));
                            end
                            
                        end
                        
                        
                    elseif mixt == 2
                        
                        log_lh=NaN(sum(niini),size(not_empty_g,2));
                        % Select all units not trimmed and not thinned
                        % Target function is based on these units (not
                        % trimmed and not thinned)
                        ytri=y(indmax>0);
                        Xtri=X(indmax>0,:);
                        
                        for jj = 1:k
                            if equalweights ==1
                                log_lh(:,jj) = ...
                                    log(1/k) + (logmvnpdfFS(ytri -Xtri * Beta(:,jj),0,sigmaini(jj) ) );
                            else
                                log_lh(:,jj) = ...
                                    log(niini(jj)/sum(niini)) + (logmvnpdfFS(...
                                    ytri - Xtri*Beta(:,jj),0,sigmaini(jj) ) );
                            end
                        end
                        
                        obj = estepFS(log_lh);
                        %                 if ~isempty(group_missing)
                        %                     nameYY(:,group_missing) = NaN;
                        %                     sigmaini(group_missing) = NaN;
                        %                 end
                    end
                    
                    
                end
            end % End of concentration steps
            
            %% Change the 'optimal' target value and 'optimal' parameters
            
            %the following two checks have to be commented if one can accept an optimal
            %subset with less than k groups.
            %if  (length(not_empty_g )== k && sum(sum(isnan(Beta))) == 0)

                % The following if condition is done to check if an increase in the target value is achieved
                if obj >= vopt
                     vopt                    = obj;
                    bopt                    = Beta;
                    numopt                  = niini;
                    sigmaopt                = sigmaini;
                    w4trimopt                   = w4trim;
                    weopt                   = we;
                    trim1levelopt           = indmax(indmax==0);
                    trim2levelopt           = indmax(indmax==-1);
                    indmax_before_tropt     = indmax_before_tr;
                    indmaxopt               = indmax;
                    %                     not_empty_gopt          = not_empty_g;
                    %                     n_refsteps_lt_k_gropt       = n_refsteps_lt_k_gr;
                    %                     n_refsteps_eq_k_gropt       = n_refsteps_eq_k_gr;
                    %                     extra_inisubsopt        = nselected - nsampdef;
                    %                     good_subs_idopt            = good_subs_id;
                    %                     elim_subs_idopt            = elim_subs_id(1:count_elim,:);
                    %                     all_subs_idopt             = [elim_subs_id(1:count_elim,:); good_subs_id];
                    if mixt ==2
                        postprobopt = postprob;
                    end
                    
                end
        %end
    end % end of loop over the nsamp subsets


    
    %%  Compute quantities to be stored in the output structure or used in the plots
    
    % Apply restriction if the scatters do not satisfy the user condition.
    % Note that the old (inefficient) approach was to use a quadratic programming
    % optimization, using function quadi, i.e. by running
    % sigmaopt_0 = (quadi(sigmaopt.^(-1), restrfact)).^(-1);
    if sum(isnan(sigmaopt)) == 0
        sigmaopt_0 = restreigen(sigmaopt,numopt,restrfact,tolrestreigen,userepmat);
    end
    
    % Apply consistency factor based on the variance of the truncated normal distribution.
    
    % number of non trimmed observations, after first and second level trimming
    hh = sum(numopt);
    
    % Compute variance of the truncated normal distribution
    % Note that 1-hh/n is the trimming percentage
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
    
    %     % units classification: 0 if trimmed at the 1st step, 1 otherwise
    %     trim1level_01opt = ones(n,1);
    %     trim1level_01opt(trim1levelopt) = 0;
    %
    %     % units classification: 0 if trimmed at the 2st step, 1 otherwise
    %     trim2level_01opt = ones(n,1);
    %     trim2level_01opt(trim2levelopt) = 0;
    %
    %     % final cluster assigments after first trimming ('0' means a trimmed observation)
    %     asig1 = indmaxopt .* weopt .* trim1level_01opt ;
    %     % final cluster assigments after first and second trimming ('0' means a trimmed observation)
    %     asig2 = asig1 .* trim2level_01opt ;
    
    %%  Set the output structure
    
    out                     = struct;
    %   bopt           = regression parameters
    out.bopt                = bopt;
    %   sigmaopt0      = estimated group variances
    out.sigmaopt_0          = sigmaopt_0;
    %   sigmaopt_cons  = estimated group variances corrected with  asymptotic consistency factor
    out.sigmaopt_cons       = sigmaopt_cons;
    %   sigmaopt_pison = estimated group variances corrected with  asymptotic consistency factor and small sample correction factor of Pison et al.
    out.sigmaopt_pison      = sigmaopt_pison;
    %   numopt         = number of observations in each cluster after the second trimming
    out.numopt              = numopt;
    %   vopt           = value of the target function
    out.vopt                = vopt;
    out.asig_obs_to_group   = indmaxopt;
    out.asig_obs_to_group_before_tr   = indmax_before_tropt;
    out.trim1levelopt       = trim1levelopt;
    out.trim2levelopt           = trim2levelopt;
    out.weopt               =weopt;
    out.C                   =C;
    %     out.asig1               = asig1;
    %     out.asig2               = asig2;
    %     out.trim2levelopt       = trim2level_01opt;
    %     out.count1_ng_lt_k      = n_refsteps_lt_k_gropt;
    %     out.count1_ng_eq_k      = n_refsteps_eq_k_gropt;
    %     out.numb_opt_sample_with_lt_k_gr = ntclust_no_opt_subs;
    %     out.extra_inisubs       = extra_inisubsopt;
    %     out.selj_good           = good_subs_idopt;
    %     out.selj_elim           = elim_subs_idopt;
    %     out.selj_all            = all_subs_idopt;
    if mixt == 2
        out.postprobopt     = postprobopt;
    end
    out.restrfact           = restrfact;
    
    
    
    % Store the indices in varargout
if nargout==2
    varargout={C};
end


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
            title({ ['TclustReg clustering: ' 'mixt=' num2str(mixt) ' - wtrim=' num2str(wtrim)] },'Fontsize',12);
            
            
            for jj = 1:k
                group_label = ['Group ' num2str(jj)];
                
                % plot of the good units allocated to the current group.
                % Indices are taken after the second level trimming.
                % Trimmed points are not plotted by group.
                ucg = find(indmaxopt==jj);
                
                
                % misteriously text does not show a legend. This is why
                % we add a (ficticious) plot instruction with white symbols
                plot(X(ucg,end),y(ucg),'.w','DisplayName',[group_label ' (' num2str(length(ucg)) ' units)']);
                text(X(ucg,end),y(ucg),num2str(jj*ones(length(ucg),1)),...
                    'DisplayName',[group_label ' (' num2str(length(ucg)) ' units)'] , ...
                    'HorizontalAlignment','center','VerticalAlignment','middle',...
                    'Color',clrdef(jj));
                
                % plot regression lines
                vv = [min(X(:,end)) max(X(:,end))];
                if intercept==1
                    plot(vv,bopt(1,jj)+bopt(2,jj)*vv,'DisplayName',[group_label ' fit' ],...
                        'Color',clrdef(jj));
                elseif intercept==0
                    plot(vv,bopt(:,jj)*vv,'DisplayName',[group_label ' fit' ],...
                        'Color',clrdef(jj));
                end
                
                %plot the thinned units
                if wtrim == 3
                    ucg = find(w4trimopt == 0);
                    plot(X(ucg,end),y(ucg),symdef(jj),'color',clrdef(k+1),...
                        'DisplayName',['Group ' num2str(jj) ': thinned units (' num2str(length(ucg)) ')' ]);
                end
                
            end
            
            % Plot the outliers (trimmed points)
            ucg = find(indmaxopt==0);
            plot(X(ucg,end),y(ucg),'o','color','r','MarkerSize',8,...
                'DisplayName',['Trimmed units (' num2str(length(ucg)) ')']);
            
            % Second level trimming points
            ucg = find(indmaxopt==-1);
            plot(X(ucg,end),y(ucg),'*','color','c',...
                'DisplayName',['L2 trimmed units (' num2str(length(ucg)) ')']);
            
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
            nameY = cellstr([repmat('X',size(YY,2)-1,1) , num2str((1:size(YY,2)-1)')]);
            nameY = [nameY ; 'y'];
            nameY = nameY';
            plo=struct;
            plo.nameY=nameY;
            
            % group names in the legend
            group = cell(n,1);
            group(indmaxopt==0) = {'Trimmed units'};
            group(indmaxopt==-1) = {'Trimmed units level 2'};
            for iii = 1:k
                group(indmaxopt==iii) = {['Group ' num2str(iii)]};
            end
            
            % scatterplot
            spmplot(YY,group,plo,'hist');
            
        end
        
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
%
%     %Check if optimization toolbox is installed in current computer
%     typemin = exist('fminunc','file');
%     if typemin ~=2
%         error('FSDA:tclustreg:MissingOptToolbox','This function requires the optimization toolbox');
%     end
%
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

