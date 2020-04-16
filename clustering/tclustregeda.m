function [out, varargout] = tclustregeda(y,X,k,restrfact,alphaLik,alphaX,varargin)
%tclustregeda performs robust linear grouping analysis for a series of values of the trimming factor
%
%<a href="matlab: docsearchFS('tclustregeda')">Link to the help function</a>
%
%   tclustregeda performs tclustregreg for a series of values of the trimming
%   factor alpha given k (number of groups) and given restrfactor
%   (restriction factor) and alphaX (second level trimming or cluster
%   weighted model). In order to increase the speed of the computations,
%   parfor is used.
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
% restrfact : restriction factor for regression residuals and covariance
%             matrices of the explanatory variables. Scalar or vector with two
%             elements. If restrfact is a scalar it controls the
%             differences among group scatters of the residuals. The value
%             1 is the strongest restriction. If restrfactor is a vector
%             with two elements the first element controls the differences
%             among group scatters of the residuals and the second the
%             differences among covariance matrices of the explanatory
%             variables. Note that restrfactor(2) is used just if
%             input option $alphaX=1$, that is if constrained weighted
%             model for X is assumed.
%            Data Types - single|double
%
%     alphaLik: trimming level to monitor. Vector. Vector which specifies the
%               values of trimming levels which have to be considered.
%               alpha is a vector which contains decreasing elements which
%               lie in the interval 0 and 0.5.
%               For example if alpha=[0.1 0.05 0] tclustregeda considers these 3
%               values of trimming level.
%               If alphaLik=0 tclustregeda does not trim. The default for
%               alphaLik is vector [0.1 0.05 0]. The sequence is forced to be
%               monotonically decreasing.
%   alphaX : Second-level trimming or constrained weighted model for X. Scalar.
%            alphaX is a value in the interval [0 1].
%            - If alphaX=0 there is no second-level trimming.
%            - If alphaX is in the interval [0 0.5] it indicates the
%               fixed proportion of units subject to second level trimming.
%               In this case alphaX is usually smaller than alphaLik.
%               For further details see Garcia-Escudero et. al. (2010).
%            -  If alphaX is in the interval (0.5 1), it indicates a
%               Bonferronized confidence level to be used to identify the
%               units subject to second level trimming. In this case the
%               proportion of units subject to second level trimming is not
%               fixed a priori, but is determined adaptively.
%               For further details see Torti et al. (2018).
%            -  If alphaX=1, constrained weighted model for X is assumed
%               (Gershenfeld, 1997). The CWM estimator is able to
%               take into account different distributions for the explanatory
%               variables across groups, so overcoming an intrinsic limitation
%               of mixtures of regression, because they are implicitly
%               assumed equally distributed. Note that if alphaX=1 it is
%               also possible to apply using restrfactor(2) the constraints
%               on the cov matrices of the explanatory variables.
%               For further details about CWM see Garcia-Escudero et al.
%               (2017) or Torti et al. (2018).
%            Data Types - single|double
%
%  Optional input arguments:
%
%     intercept : Indicator for constant term. Scalar. If 1, a model with
%                constant term will be fitted (default), else no constant
%                term will be included.
%                Example - 'intercept',1
%                Data Types - double
%
%       mixt  : mixture modelling or crisp assignment. Scalar.
%               Option mixt specifies whether mixture modelling or crisp
%               assignment approach to model based clustering must be used.
%               In the case of mixture modelling parameter mixt also
%               controls which is the criterior to find the untrimmed units
%               in each step of the maximization
%               If mixt >=1 mixture modelling is assumed else crisp
%               assignment.
%                In mixture modelling the likelihood is given by
%                \[
%                \prod_{i=1}^n  \sum_{j=1}^k \pi_j \phi (y_i \; x_i' , \beta_j , \sigma_j),
%                \]
%               while in crisp assignment the likelihood is given by
%               \[
%               \prod_{j=1}^k   \prod _{i\in R_j} \pi_j  \phi (y_i \; x_i' , \beta_j , \sigma_j),
%               \]
%               where $R_j$ contains the indexes of the observations which
%               are assigned to group $j$,
%               Remark - if mixt>=1 previous parameter equalweights is
%               automatically set to 1.
%               Parameter mixt also controls the criterion to select the units to trim
%               if mixt == 2 the h units are those which give the largest
%               contribution to the likelihood that is the h largest
%               values of
%               \[
%                   \sum_{j=1}^k \pi_j \phi (y_i \; x_i' , \beta_j , \sigma_j)   \qquad
%                    i=1, 2, ..., n
%               \]
%               elseif mixt==1 the criterion to select the h units is
%               exactly the same as the one which is used in crisp
%               assignment. That is: the n units are allocated to a
%               cluster according to criterion
%               \[
%                \max_{j=1, \ldots, k} \hat \pi'_j \phi (y_i \; x_i' , \beta_j , \sigma_j)
%               \]
%               and then these n numbers are ordered and the units
%               associated with the largest h numbers are untrimmed.
%               Example - 'mixt',1
%               Data Types - single | double
%
%equalweights : cluster weights in the concentration and assignment steps.
%               Logical. A logical value specifying whether cluster weights
%               shall be considered in the concentration, assignment steps
%               and computation of the likelihood.
%               if equalweights = true we are (ideally) assuming equally
%               sized groups by maximizing the likelihood
%                 Example - 'equalweights',true
%                 Data Types - Logical
%
%    nsamp : number of subsamples to extract.
%            Scalar or matrix with k*p columns.
%            If nsamp is a scalar it contains the number of subsamples
%            which will be extracted.
%            If nsamp=0 all subsets will be extracted.
%            If the number of all possible subset is <300 the
%            default is to extract all subsets, otherwise just 300.
%            If nsamp is a matrix it contains in the rows the indexes of
%            the subsets which have to be extracted. nsamp in this case can
%            be conveniently generated  by function subsets.
%            nsamp must have k*p columns. The first p columns are used to
%            estimate the regression coefficient of group 1... the last p
%            columns are used to estimate the regression coefficient of
%            group k
%             Example - 'nsamp',1000
%             Data Types - double
%
% refsteps:  Number of refining iterations. Scalar. Number of refining
%               iterations in each subsample.  Default is 10.
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'refsteps',15
%                 Data Types - single | double
%
%     reftol  : Tolerance for the refining steps. Scalar.
%               The default value is 1e-14;
%                 Example - 'reftol',1e-05
%                 Data Types - single | double
%
%
%
% plots    :    Plot on the screen. Scalar structure.
%
%               Case 1: plots option used as scalar.
%               - If plots=0,  plots are not generated.
%               - If plots=1 (default), 6 plots are shown on the screen.
%                 The first plot ("monitor") shows 4 or 6 panels
%                 monitoring between two consecutive values of alpha the
%                 change in classification using ARI index (top left panel),
%                 the relative change in beta (top right panel)
%                 ($||\beta_{\alpha_1}-\beta_{\alpha_2}||^2/||beta_{\alpha_2}||^2
%                 the relative change in \sigma^2 (third panel) the
%                 relative change in \sigma^2 corrected (fourth panel) and
%                 if alphaX=1, the relative change in centroids (fifth
%                 panel), the relative change in covariance matrices using
%                 squared euclidean distance (sixth panel).
%                 The second plot ("UnitsTrmOrChgCla") is a bsb plot which
%                 shows the monitoring of the units which changed
%                 classification (shown in red) or were trimmed at least
%                 once.
%                 The third plot shows ("PostProb") is a 4 panel parallel
%                 coordinate plots which shows the monitoring of posterior
%                 probabilities. Four versions of parallel coordinateds
%                 plots are proposed. The first two use function
%                 parallelcoords and refer respectively to all the units
%                 (left panel) or the units which were trimmed or changed
%                 allocation (right panels). The two bottom panels
%                 are exactly equal as top panels but use function
%                 parallelplot.
%                 The fourth plot ("PostProb") is 4 panel parallel
%                 coordinate plots which shows the monitoring of posterior
%                 probabilities. Four versions of parallel coordinateds
%                 plots are proposed. The first two use function
%                 parallelcoords and refer respectively to all the units
%                 (left panel) or the units which were trimmed or changed
%                 allocation (right panels). The two bottom panels
%                 are exactly equal as top panels but use function
%                 parallelplot.
%                 The fifth plot is a scatter of y versus X (or a yXplot in
%                 case X ha more than one column) with all the regression
%                 lines for each value of alphaLik shown. Units trimmed are
%                 shown in correspondence of max(alphaLik).
%                 The sixth plot ("gscatter plot") shows a series of
%                 subplots which monitor the classification for each value
%                 of alpha. In order to make sure that consistent labels
%                 are used for the groups, between two consecutive values
%                 of alpha, we assign label r to a group if this group
%                 shows the smallest distance with group r for the previous
%                 value of alpha. The type of plot which is used to monitor
%                 the stability of the classification depends on the value
%                 of p (number of explanatory variables excluding the intercept).
%                   * for p=0, we use histograms of the univariate data
%                   (function histFS is called).
%                   * for p=1, we use the scatter plot of y
%                   against the unique explanatory variable (function gscatter is called).
%                   * for p>=2, we use partial least square regression (see
%                   function plsregress) and use the scatter plot of y
%                   against the predictor scores Xs, that is, the first PLS
%                   component that is linear combination of the variables
%                   in X.
%
%               Case 2: plots option used as struct.
%                 If plots is a structure it may contain the following fields:
%                 plots.name = cell array of strings which enables to
%                   specify which plot to display.
%                   plots.name={'monitor'; 'UnitsTrmOrChgCla'; 'PostProb'; 'Sigma';...
%                               'ScatterWithRegLines'; 'gscatter'};
%                   is exactly equivalent to plots=1
%                   For the explanation of the above plots see plots=1.
%                   If plots.name=={ 'monitor'; 'UnitsTrmOrChgCla'; 'PostProb'; 'Sigma';...
%                                   'ScatterWithRegLines'; 'gscatter'; ...
%                                   'Beta';'Siz'}; or
%                   plots.name={'all'};
%                   it is also possible to monitor the beta coefficients
%                   for each group ('Beta') and the size of the groups (Siz).
%                   Note that the trajectories of beta coefficients are
%                   standardized in order to have all of them on a
%                   comparable scale.
%                 plots.alphasel = numeric vector which speciies for which
%                   values of alpha it is possible to see the
%                   classification (in plot gscatter) or the
%                   superimposition of regression lines (in plot
%                   ScatterWithRegLines) .
%                   For example if plots.alphasel =[ 0.05 0.02], the
%                   classification in plot gscatter and the regression
%                   lines in plot ScatterWithRegLines will be shown just
%                   for alphaLik=0.05 and alphaLik=0.02; If this field is not
%                   specified plots.alphasel=alphaLik and therefore the
%                   classification is shown for each value of alphaLik.
%                   Example - 'plots', 1
%                   Data Types - single | double | struct
%
%        msg  : Level of output to display. Scalar.
%               Scalar which controls whether to display or not messages
%               on the screen.
%               If msg=0 nothing is displayed on the screen.
%               If msg=1 (default) messages are displayed
%               on the screen about estimated time to compute the estimator
%               or the number of subsets in which there was no convergence.
%               If msg=2 detailed messages are displayed. For example the
%               information at iteration level.
%                   Example - 'msg',1
%                   Data Types - single | double
%
%      nocheck: Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               vector y and matrix X.
%               As default nocheck=0.
%                   Example - 'nocheck',1
%                   Data Types - single | double
%
%
%RandNumbForNini: Pre-extracted random numbers to initialize proportions.
%                Matrix. Matrix with size k-by-size(nsamp,1) containing the
%                random numbers which are used to initialize the
%                proportions of the groups. This option is effective just
%                if nsamp is a matrix which contains pre-extracted
%                subsamples. The purpose of this option is to enable the
%                user to replicate the results in case routine tclust is
%                called using a parfor instruction (as it happens for
%                example in routine IC, where tclust is called through a
%                parfor for different values of the restriction factor).
%                The default value of RandNumbForNini is empty that is
%                random numbers from uniform are used.
%                   Example - 'RandNumbForNini',''
%                   Data Types - single | double
%
%
%   UnitsSameGroup :  list of the units which must (whenever possible)
%                   have a particular label. Numeric vector.  For example if
%                   UnitsSameGroup=[20 26], means that group which contains
%                   unit 20 is always labelled with number 1. Similarly,
%                   the group which contains unit 26 is always labelled
%                   with number 2, (unless it is found that unit 26 already
%                   belongs to group 1). In general, group which contains
%                   unit UnitsSameGroup(r) where r=2, ...length(kk)-1 is
%                   labelled with number r (unless it is found that unit
%                   UnitsSameGroup(r) has already been assigned to groups
%                   1, 2, ..., r-1).
%                 Example - 'UnitsSameGroup',[20 34]
%                 Data Types -  integer vector
%
%      numpool:     The number of parallel sessions to open. Integer. If
%                   numpool is not defined, then it is set equal to the
%                   number of physical cores in the computer.
%                 Example - 'numpool',4
%                 Data Types -  integer vector
%
%      cleanpool:   Function name. Scalar {0,1}. Indicated if the open pool
%                   must be closed or not. It is useful to leave it open if
%                   there are subsequent parallel sessions to execute, so
%                   that to save the time required to open a new pool.
%                 Example - 'cleanpool',true
%                   Data Types - integer | logical
%
%      we: Vector of observation weights. Vector. A vector of size n-by-1
%          containing application-specific weights that the user needs to
%          apply to each observation. Default
%          value is  a vector of ones.
%            Example - 'we',[0.2 0.2 0.2 0.2 0.2]
%            Data Types - double
%
%  Output:
%
%         out:   structure which contains the following fields
%
%            out.IDX  = n-by-length(alphaLik) matrix containing assignment of each unit to
%                       each of the k groups. Cluster names are integer
%                       numbers from 1 to k:
%                       -1 indicates first level trimmed observations;
%                       -2 indicates second level trimmed observations.
%                       observations. First column refers of out.IDX refers
%                       to alphaLik(1), second column of out.IDX refers to
%                       alphaLik(2), ..., last column refers to alphaLik(end).
%
%            out.Beta  =  3D array of size k-by-p-by-length(alphaLik) containing
%                       the monitoring of the regression coefficients for each value of
%                       alphaLik. out.Beta(:,:,1), refers to alphaLik(1) ...,
%                       out.Beta(:,:,end) refers to alphaLik(end). First row in
%                       each slice refers to group 1, second row refers to
%                       group 2 ...
%
%         out.Sigma2y  =  matrix of size k-by-length(alphaLik) containing in column
%                       j, with j=1, 2, ...,  length(alphaLik), the
%                       estimates of the k (constrained) variances of the
%                       regressions lines (hyperplanes) associated with alphaLik(j).
%
%         out.Sigma2yc  =  matrix of size k-by-length(alphaLik) containing in column
%                       j, with j=1, 2, ...,  length(alphaLik), the
%                       estimates of the k (constrained) unbiased variances of the
%                       regressions lines (hyperplanes) associated with alphaLik(j).
%                       In order to make the estimates of sigmas unbiased
%                       we apply Tallis correction factor.
%
%         out.Nopt  =  matrix of size k-by-length(alphaLik) containing in column
%                       j, with j=1, 2, ...,  length(alphaLik), the
%                       sizes of the of the k groups.
%
%         out.Vopt  =  column vector of length(alphaLik) containing the
%                      value of the target likelihod for each value of alphaLik.
%
%
%         out.Amon  =  Amon stands for alphaLik monitoring. Matrix of size
%                      (length(alphaLik)-1)-by-7 which contains for two
%                       consecutive values of alpha the monitoring of six
%                       quantities (change in classification, change in
%                       betas, sigmas, correted sigmas and if cwm also
%                       centroid and coariance in the X space.
%                       1st col = value of alphaLik.
%                       2nd col = ARI index.
%                       3rd col = relative squared Euclidean distance between
%                           two consecutive beta.
%                       4th col = relative squared Euclidean distance between
%                           two consecutive vector of variances of the
%                           errors of the k regressions.
%                       5th col = relative squared Euclidean distance between
%                           two consecutive vector of correct variances of the
%                           errors of the k regressions.
%                       6th col = relative squared Euclidean distance between
%                           two consecutive $\hat \mu_X$.
%                       7th col = relative squared Euclidean distance between
%                           two consecutive $\hat \Sigma_X$.
%
% out.unitplot = Matrix containing information about the
%                       units (n1) which were trimmed or changed classification
%                       at least once in the forward search. The size of
%                       out.UnitsTrmOrChgCla is n1-by-length(alphaLik)+1;
%                       1st col = list of the units which were trimmed or
%                       changed classification at least once.
%                       2nd col = allocation of the n1 units in step alphaLik(1)
%                       3rd col = allocation of the n1 units in step alphaLik(2)
%                       ...
%                       last col = allocation of the n1 units in step alphaLik(end)
%
%     out.units= structure containing the following fields:
%                       units.UnitsTrmOrChgCla=units trimmed at least onece
%                           or changed classification at least once.
%                       units.UnitsChgCla=units which changed
%                           classification at least once (i.e. from group 1
%                           to group 3 ...).
%                       units.UnitsTrm=units trimmed at least once.
%                       units.UnitsNeverAssigned=units never assigned
%                           (i.e. all those which have always been trimmed by
%                           first level or second level).
%
%              out.y  = original response y.
%
%              out.X  = original X matrix.
%
%  Optional Output:
%
%           outcell : cell of length length(alpha) which contains in jth
%                     position the structure which comes out from procedure
%                     tclustreg applied to alphaLik(j), with j =1, 2, ...,
%                     length(alphaLik).
%
% More About:
%
%
% This procedure extends to tclustreg the so called monitoring approach.
% The phylosophy is to investigate how the results change as the trimming
% proportion alpha reduces. This function enables us to monitor the change
% in classification (measured by the ARI index) and the change in
% regression coefficients and error variances (measured by the relative
% squared euclidean distances). In order to make sure that consistent
% labels are used for the groups, between two consecutive values of alpha,
% we assign label r to a group if this group shows the smallest distance
% with group r for the previous value of alpha. In otder to be consistent
% along the different runs it is possible to specify through option
% UnitsSameGroup the list of the units which must (whenever possible) have
% a particular label.
%
%
% See also: tclustreg, tclustIC, tclusteda
%
% References:
%
% Garcia-Escudero, L.A., Gordaliza A., Greselin F., Ingrassia S., and
% Mayo-Iscar A. (2016), The joint role of trimming and constraints in
% robust estimation for mixtures of gaussian factor analyzers,
% "Computational Statistics & Data Analysis", Vol. 99, pp. 131-147.
%
% Garcia-Escudero, L.A., Gordaliza, A., Greselin, F., Ingrassia, S. and
% Mayo-Iscar, A. (2017), Robust estimation of mixtures of regressions with
% random covariates, via trimming and constraints, "Statistics and
% Computing", Vol. 27, pp. 377-402.
%
% Garcia-Escudero, L.A., Gordaliza A., Mayo-Iscar A., and San Martin R.
% (2010), Robust clusterwise linear regression through trimming,
% "Computational Statistics and Data Analysis", Vol. 54, pp.3057-3069.
%
% Cerioli, A. and Perrotta, D. (2014). Robust Clustering Around Regression
% Lines with High Density Regions. Advances in Data Analysis and
% Classification, Vol. 8, pp. 5-26.
%
% Torti F., Perrotta D., Riani, M. and Cerioli A. (2018). Assessing Robust
% Methodologies for Clustering Linear Regression Data, "Advances in Data
% Analysis and Classification".
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('tclustregeda')">Link to the help function</a>
%
%$LastChangedDate:: 2018-01-29 18:52:14 #$: Date of the last commit

% Examples:

%{
    %% tclustreg of contaminated X data using all default options.
    % The X data have been introduced by Gordaliza, Garcia-Escudero & Mayo-Iscar (2013).
    % The dataset presents two parallel components without contamination.
    X  = load('X.txt');
    y = X(:,end);
    X =X(:,1:end-1);
    % Contaminate the first 4 units
    y(1:4)=y(1:4)+6;
    % Use 2 groups
    k = 2 ;
    % Set restriction factor
    restrfact = 5;
    % Value of trimming
    alphaLik = 0.10:-0.01:0;
    % cwm
    alphaX = 1;
    out = tclustregeda(y,X,k,restrfact,alphaLik,alphaX);
%}

%{
    %% tclustreg with a noise variable and personalized plots.
    % Use the X data of the previous example
    X  = load('X.txt');
    y = X(:,end);
    rng(100)
    X = [randn(length(y),1) X(:,1:end-1)];
    % Contaminate the first 4 units
    y(1:4)=y(1:4)+6;
    % Use 2 groups
    k = 2 ;

    restrfact = 5;
    % Value of trimming
    alphaLik = [0.10 0.06 0.03 0];
    % cwm
    alphaX = 1;
    % Personalize plots. Just show the gscatter plot.
    % In this case there is more than one explanatory variable therefore PLS
    % regression (adding the dummies for the classified units) is performed.
    % In the gscatter plots the percentage of variance explained by the first
    % linear combination of the X variables is given in the title.
    plots=struct;
    plots.name={'gscatter'};

    out = tclustregeda(y,X,k,restrfact,alphaLik,alphaX,'plots',plots);
%}
%% Beginning of code
% Control variables, tolerances and internal flags
warning('off');


%% Input parameters checking

nnargin   = nargin;
vvarargin = varargin;

[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% Check the presence of the intercept in matrix X
if min(max(X)-min(X))==0
    intercept = 1;
else
    intercept = 0;
end

% check restrfact option
if nargin < 4 || isempty(restrfact) || ~isnumeric(restrfact)
    restrfact = 12;         % This is the default in R
elseif min(restrfact)<1
    disp('Restriction factor smaller than 1. It is set to 1 (maximum constraint==>spherical groups)');
    restrfact(restrfact<1)=1;
else
end

% Check if restrfact is a scalar or a vector or length 2 (i.e. if the
% restriction is applied also on the explanatory variables)
if length(restrfact)>1
    % restriction factor among eigenvalues of covariance matrices of
    % explanatory variables
    restrfactX=restrfact(2);
    restrfact=restrfact(1);
else
    % No restriction of scatter matrices of explanatory variables
    restrfactX=Inf;
end

% checks on alpha1 (alphaLik) and alpha2 (alphaX)
if min(alphaLik) < 0 || max(alphaLik)>0.5
    error('FSDA:tclustreg:error','error must a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
end

% checks on cwm, which decides if clusterwise regression has to be used
if alphaX < 0 || alphaX >1
    error('FSDA:tclustreg:WrongAlphaX','alphaX must a scalar in the interval [0 1]')
elseif alphaX==1
    cwm=1;
else
    cwm=0;
end


%% User options and their default values

%%% - nsamp: the number of subsets to extract randomly, or the indexes of the initial subsets pre-specified by the User
numpool = feature('numCores');
cleanpool=false;

if nargin>6
    
    % Check whether option nsamp exists
    chknsamp = strcmp(varargin,'nsamp');
    if sum(chknsamp)>0
        nsamp=cell2mat(varargin(find(chknsamp)+1));
        
        % Check if options nsamp is a scalar
        if ~isscalar(nsamp)
            % if nsamp is not a scalar, it is a matrix containing in the
            % rows the indexes of the subsets which have to be extracted
            C=nsamp;
            [nsampdef,ncolC]=size(C);
            if ncolC ~= k*p
                disp('p is the total number of explanatory variables (including the constant if present)')
                error('FSDA:tclustreg:WrongInput','Input matrix C must contain k*p columns')
            end
            % The number of rows of nsamp (matrix C) is the number of
            % subsets which have to be extracted
            nselected=nsampdef;
            
            % Flag indicating if the user has selected a prior subset
            NoPriorSubsets=0;
            
            % In case of tandem thinning (Wtrim=4), the initial subset C
            % pre-specified by the user using the nsamp option might
            % include thinned units. In this case, we replace in C such
            % units with others that are close in terms of euclidean
            % distance. Verify the contiguity between the original and
            % replaced units with:
            %{
             figure;plot(Xori,yori,'.'); hold on ; text(Xori(nsamp(:)),yori(nsamp(:)),'X');
             figure;plot(Xori,yori,'.'); hold on ; text(X(C(:)),y(C(:)),'X');
            %}
            if sum(chknwtrim)>0 && ~isstruct(cell2mat(varargin(find(chknwtrim)+1)))
                if cell2mat(varargin(find(chknwtrim)+1))== 4
                    for f=1:size(C,1)*size(C,2)
                        if Wt4(C(f)) == 0
                            [~,C(f)] = min(pdist2([yori(C(f)),Xori(C(f))],[yori(Wt4),Xori(Wt4)]));
                        else
                            C(f) = sum(Wt4(1:C(f)));
                        end
                    end
                end
            end
            
        else
            % If nsamp is a scalar it simply contains the number of subsets
            % which have to be extracted. In this case NoPriorSubsets=1
            NoPriorSubsets=1;
        end
    else
        % If option nsamp is not supplied, then there are no prior subsets
        NoPriorSubsets=1;
    end
else
    % if nargin == 6, then the user has not supplied prior subsets
    NoPriorSubsets=1;
end
% If the user has not specified prior subsets (nsamp is not a scalar), then
% set the default number of samples to extract
if NoPriorSubsets == 1
    ncomb=bc(n,k*(p+intercept));
    nsampdef=min(300,ncomb);
end

%%% - Other user options

% default number of concentration steps
refstepsdef  = 10;

% default tolerance for comparing the classifications in two subsequent
% concentration steps
reftoldef=1e-5;

% default value for we: the observation weights
wedef = ones(n,1);

% default model: classification (mixt=0) or mixture likelihood (mixt=2)
mixtdef = 0;

% default choice for equalweight constraint
equalweightsdef = 1;

%seqk = sequence from 1 to the number of groups
seqk = 1:k;

plots=1;
UnitsSameGroup='';

% automatic extraction of user options
options = struct('intercept',1,'mixt',mixtdef,...
    'nsamp',nsampdef,'refsteps',refstepsdef,...
    'reftol',reftoldef,...
    'we',wedef,'numpool',numpool,'cleanpool', cleanpool,...
    'equalweights',equalweightsdef,...
    'RandNumbForNini','','msg',1,'plots',plots,...
    'nocheck',1,'UnitsSameGroup',UnitsSameGroup);

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
    
    % Check if the number of subsamples to extract is reasonable
    if isscalar(options.nsamp) && options.nsamp>ncomb
        disp('Number of subsets to extract greater than (n k). It is set to (n k)');
        options.nsamp=0;
    elseif  options.nsamp<0
        error('FSDA:tclustreg:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
    end
end
numpool=options.numpool;

% global variable controlling if messages are displayed in the console.
msg = options.msg;

% Graphs summarizing the results
plots = options.plots;

% Number of subsets to extract or matrix containing the subsets
nsamp = options.nsamp;

% Concentration steps
refsteps = options.refsteps;
reftol   = options.reftol;

% Equalweights constraints
equalweights = options.equalweights;


UnitsSameGroup=options.UnitsSameGroup;

% application-specific weights vector assigned by the user for beta
% estimation
we         = options.we;


% Mixt option: type of membership of the observations to the sub-populations
% Control the mixture model to use (classification/mixture, likelihood or a
% combination of both):
%
% * mixt = 0: Classification likelihood
% * mixt = 1: Mixture likelihood, with crisp assignement
% * mixt = 2: Mixture likelihood
%
% $$ \prod_{j=1}^k  \prod_{i\in R_j} \phi (x_i;\theta_j) $$ $$ \quad $$
% $$ \prod_{j=1}^k  \prod_{i\in R_j} \pi_j \phi(x_i;\theta_j) $$ $$ \quad $$
% $$ \prod_{i=1}^n \left[ \sum_{j=1}^k \pi_j \phi (x_i;\theta_j)  \right] $$
mixt       = options.mixt;

if msg == 1
    switch mixt
        case 0
            % Classification likelihood.
            % To select the h untrimmed units, each unit is assigned to a
            % group and then we take the h best maxima
            disp('ClaLik with untrimmed units selected using crisp criterion');
        case 1
            % Mixture likelihood.
            % To select the h untrimmed units, each unit is assigned to a
            % group and then we take the h best maxima
            disp('MixLik with untrimmed units selected using crisp criterion');
        case 2
            % Mixture likelihood.
            % To select the h untrimmed units we take those with h largest
            % contributions to the likelihood
            disp('MixLik with untrimmed units selected using h largest lik contributions');
    end
end

% Initial mixing proportions $\pi_j$ can be user-defined or randomly generated
RandNumbForNini=options.RandNumbForNini;
if isempty(RandNumbForNini)
    NoPriorNini=1;
else
    NoPriorNini=0;
end

%% Initializations

%%% - Subsets extraction

%case with no prior subsets from the User
if NoPriorSubsets
    
    [C,nselected]=subsets(nsamp, n, (p+intercept)*k, ncomb, 0, we);
    
    % C = matrix which contains the indexes of the subsets to extract
    % for the k groups
    
    %nselected is set equal to the number of subsets:
    % - nselected = nsamp, if nsamp is a scalar;
    % - nselected = size(nsamp,1), if nsamp is a matrix (containing the initial subsets)
    
end

%%% - Output structures

% Store the initial subsets indices C
if nargout==2
    varargout={C};
end
% Flag associated to the strategy for choosing the best refining step
% In the standard TCLUST the best refining step is granted to be the last
% one, because the objective funcion is monothonic. However, with second
% trimming level or componentwise thinning, the objective function may not
% be monothonic and a different strategy for choosing the best refining
% step can be considered.
wtype_beta=0;
wtype_obj=0;
zigzag = (alphaX > 0 && alphaX<1) || wtype_beta == 3 || wtype_beta == 2 || ~strcmp(wtype_obj, '0');

% Make sure alphaLik is a column vector;
alphaLik=alphaLik(:);

lalpha=length(alphaLik);
if msg == 1
    progbar = ProgressBar(lalpha);
else
    progbar=[];
end

IDX=zeros(n,lalpha);

MU=zeros(k,p-1,lalpha);
SIGMA=cell(lalpha,1);
Beta=zeros(k,p,lalpha);
Nopt=zeros(k,lalpha);
Sigma2y=zeros(k,lalpha);
Sigma2yc=Sigma2y;
Vopt=zeros(lalpha,1);
Postprob=zeros(n,k,lalpha);
% Do not show messages during each execution of tclustregcore
msgrs=0;
parfor (j=1:lalpha, numpool)
    % for j=1:lalpha
    h = floor(n*(1-alphaLik(j)));
    
    %%  RANDOM STARTS
    [bopt,sigma2opt,nopt,postprobopt,muXopt,sigmaXopt,vopt,~,idxopt]...
        =tclustregcore(y,X,RandNumbForNini,reftol,refsteps,mixt,equalweights,h,nselected,k,restrfact,restrfactX,alphaLik(j),alphaX,...
        seqk,NoPriorNini,msgrs,C,intercept,cwm,wtype_beta,we,wtype_obj,zigzag);
    
    %%  END OF RANDOM STARTS
    
    IDX(:,j)=idxopt;
    Beta(:,:,j)=bopt';
    Nopt(:,j)=nopt;
    MU(:,:,j)=muXopt;
    SIGMA{j}=sigmaXopt;
    Vopt(j)=vopt;
    Postprob(:,:,j)=postprobopt;
    
    % hh = number of non trimmed observations, after first and second level trimming
    hh = sum(nopt);
    
    % vt = variance of the truncated normal distribution
    % 1-hh/n is the trimming percentage
    vt = norminv(0.5*(1+hh/n));
    
    sigma2opt=sigma2opt';
    if hh<n
        factor = 1/sqrt(1-2*(n/hh)*vt.*normpdf(vt));
        % Apply the asymptotic consistency factor to the preliminary squared scale estimate
        sigma2opt_corr = sigma2opt*factor;
    else
        sigma2opt_corr = sigma2opt;
    end
    
    Sigma2y(:,j)=sigma2opt;
    Sigma2yc(:,j)=sigma2opt_corr;
    
    if msg == 1
        progbar.progress; %#ok<PFBNS>
    end
end

if msg == 1
    progbar.stop;
end


if ~isempty(UnitsSameGroup)
    % Note that this operation is just applied fo first column of IDX
    idx1=IDX(:,1);
    trimmed1=idx1==-1;
    trimmed2=idx1==-2;
    trimmed=idx1<0;
    idxtmp=idx1;
    idxtmp(trimmed)=0;
    [IDXnew1, OldNewIndexes]=ClusterRelabel({idxtmp}, UnitsSameGroup);
    
    % MUold1 is k-by-p (rows refer to groups)
    % Mu is k-by-p-length(alpha)
    MUold1=MU(:,:,1);
    % Sigmaold1 is pxpxk
    SIGMAold1= SIGMA{1};
    
    Betaold1=Beta(:,:,1);
    Sigma2yold1=Sigma2y(:,1);
    Sigma2ycold1=Sigma2yc(:,1);
    Postprobold1=Postprob(:,:,1);
    
    MUnew1=MUold1;
    SIGMAnew1=SIGMAold1;
    Betanew1=Betaold1;
    Sigma2ynew1=Sigma2yold1;
    Sigma2ycnew1=Sigma2ycold1;
    Postprobnew1=Postprobold1;
    
    for jj=1:size(OldNewIndexes,1)
        
        % in SIGMAnew1 k is in the third dimension
        SIGMAnew1(:,:,OldNewIndexes(jj,1))=SIGMAold1(:,:,OldNewIndexes(jj,2));
        SIGMAnew1(:,:,OldNewIndexes(jj,2))=SIGMAold1(:,:,OldNewIndexes(jj,1));
        SIGMAold1=SIGMAnew1;
        
        % in MUnew1 Betanew1, Sigma2ynew1, and Sigma2ycnew1 k is in the rows
        MUnew1(OldNewIndexes(jj,1),:)= MUold1(OldNewIndexes(jj,2),:);
        MUnew1(OldNewIndexes(jj,2),:)= MUold1(OldNewIndexes(jj,1),:);
        MUold1=MUnew1;
        
        Betanew1(OldNewIndexes(jj,1),:)=Betanew1(OldNewIndexes(jj,2),:);
        Betanew1(OldNewIndexes(jj,2),:)=Betanew1(OldNewIndexes(jj,1),:);
        Betaold1=Betanew1;
        
        Sigma2ynew1(OldNewIndexes(jj,1))=Sigma2ynew1(OldNewIndexes(jj,2));
        Sigma2ynew1(OldNewIndexes(jj,2))=Sigma2ynew1(OldNewIndexes(jj,1));
        Sigma2yold1=Sigma2ynew1;
        
        Sigma2ycnew1(OldNewIndexes(jj,1))=Sigma2ycnew1(OldNewIndexes(jj,2));
        Sigma2ycnew1(OldNewIndexes(jj,2))=Sigma2ycnew1(OldNewIndexes(jj,1));
        Sigma2ycold1=Sigma2ycnew1;
        
        % In Postprobnew1 k is in the columns
        Postprobnew1(:,OldNewIndexes(jj,1))=Postprobnew1(:,OldNewIndexes(jj,2));
        Postprobnew1(:,OldNewIndexes(jj,2))=Postprobnew1(:,OldNewIndexes(jj,1));
        Postprobold1=Postprobnew1;
        
    end
    IDX(:,1)=IDXnew1{:};
    IDX(trimmed1,1)=-1;
    IDX(trimmed2,1)=-2;
    
    MU(:,:,1)=MUold1;
    SIGMA{1}=SIGMAold1;
    Beta(:,:,1)=Betaold1;
    Sigma2y(:,1)=Sigma2yold1;
    Sigma2yc(:,1)=Sigma2ycold1;
    Postprob(:,:,1)=Postprobold1;
end

IDXold=IDX;
maxdist=zeros(lalpha,1);
seqk=(1:k)';

for j=2:lalpha
    newlab=zeros(k,1);
    mindist=newlab;
    for ii=1:k
        % centroid of group ii for previous alpha value
        muii=Beta(ii,:,j-1);
        % MU(:,:,j) =matrix of centroids for current alpha value
        
        %if verMatlab==true;
        muij=bsxfun(@minus,muii,Beta(:,:,j));
        %else
        %    muij=muii-MU(:,:,j);
        %end
        
        [mind,indmin]=min(sum(muij.^2,2));
        newlab(ii)=indmin;
        mindist(ii)=mind;
    end
    % Store maximum among minimum distances
    [maxmindist,indmaxdist] =max(mindist);
    maxdist(j)=maxmindist;
    
    if isequal(sort(newlab),seqk)
        MU(:,:,j)=MU(newlab,:,j);
        SIGMA(j)= {SIGMA{j}(:,:,newlab)};
        Beta(:,:,j)=Beta(newlab,:,j);
        Sigma2y(:,j)=Sigma2y(newlab,j);
        Sigma2yc(:,j)=Sigma2yc(newlab,j);
        Postprob(:,:,j)=Postprob(:,newlab,j);
        
        for r=1:k
            IDX(IDXold(:,j)==newlab(r),j)=r;
        end
    else
        newlab(indmaxdist)=setdiff(seqk,newlab);
        disp(['Preliminary relabelling not possible when alpha=' num2str(alphaLik(j))])
        if isequal(sort(newlab),seqk)
            MU(:,:,j)=MU(newlab,:,j);
            SIGMA(j)= {SIGMA{j}(:,:,newlab)};
            Beta(:,:,j)=Beta(newlab,:,j);
            Sigma2y(:,j)=Sigma2y(newlab,j);
            Sigma2yc(:,j)=Sigma2yc(newlab,j);
            Postprob(:,:,j)=Postprob(:,newlab,j);
            
            
            for r=1:k
                IDX(IDXold(:,j)==newlab(r),j)=r;
            end
        else
            disp(['Automatic relabelling not possible when alpha=' num2str(alphaLik(j))])
        end
    end
end


%%  Set the output structure

out       = struct;
out.IDX   = IDX;
out.Sigma2y=Sigma2y;
out.Sigma2yc=Sigma2yc;
out.Beta=Beta;
out.Nopt=Nopt;
out.MU=MU;
out.SIGMA=SIGMA;
out.Vopt=Vopt;
out.Postprob=Postprob;

% Store the indices in varargout
if nargout==2
    varargout={C};
end


%% Generate plots
% plotdef = list of the plots which are produced by default (is plots=1)
plotdef={'monitor'; 'UnitsTrmOrChgCla'; 'PostProb'; 'Sigma';...
    'ScatterWithRegLines'; 'gscatter'};
% plotall = list of all available plots
plotall={'monitor'; 'UnitsTrmOrChgCla'; 'PostProb'; 'Sigma';...
    'ScatterWithRegLines'; 'gscatter'; ...
    'Beta'; 'Siz'};

clrdef = 'bkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcr';
symdef = '+sd^v><phos+*d^v><phos+*d^v><phos+*d^v><phos';

col1stLevelTrimmedUnits='r';
sym1stLevelTrimmedUnits='o';

col2ndLevelTrimmedUnits='c';
sym2ndLevelTrimmedUnits='*';

if isstruct(plots)
    fplots=fieldnames(plots);
    
    d=find(strcmp('name',fplots));
    if d>0
        name=plots.name;
        if ~iscell(name)
            error('FSDA:tclustregeda:Wronginput','plots.name must be a cell')
        end
        if strcmp(name,{'all'})
            name=plotall;
        else
            % Check that the specified names is in the list of available names.
            chkoptions(cell2struct(plotall,plotall),name)
        end
    else
        name=plotdef;
    end
    
    d=find(strcmp('alphasel',fplots));
    if d>0
        alphasel=plots.alphasel;
    else
        alphasel=alphaLik;
    end
    
elseif plots==1
    name=plotdef;
    alphasel=alphaLik;
    % ylimy='';
end



%% Create a series of variables which will be used in more than one plot and
% therefore are computed once and for all
IDXmin0=IDX<=0;
IDXwithNaN=IDX;
IDXwithNaN(IDXmin0)=NaN;
seq=(1:n)';

% Find the rows which are not constant, that is the units which were trimmed
% or those which changed classification
UnitsTrmOrChgCla=seq(max(IDX,[],2)-min(IDX,[],2)~=0 | (max(IDX,[],2)<0) );

% Find units which changed classification at least once
UnitsChgCla=seq(max(IDXwithNaN,[],2)-min(IDXwithNaN,[],2)~=0 & (max(IDXwithNaN,[],2)>0));

% UnitsTrm= list of the units trimmed (first level or seond level) at least once.
UnitsTrm=setdiff(UnitsTrmOrChgCla,UnitsChgCla);

% Units which have never been assigned
UnitsNeverAssigned=seq(sum(isnan(IDXwithNaN),2)==lalpha);

units=struct;
units.UnitsTrmOrChgCla=UnitsTrmOrChgCla;
units.UnitsChgCla=UnitsChgCla;
units.UnitsTrm=UnitsTrm;
units.UnitsNeverAssigned=UnitsNeverAssigned;

out.units=units;


% String to include in the legends
legendGroups=[repmat('Group ',k,1) num2str((1:k)')];


% Amon stands for alpha monitoring.
% Amon is the matrix of size lenght(alpha)-1-by 4 which contains for two
% consecutive values of alpha the monitoring of three quantities.
% 1st col = value of alpha
% 2nd col = ARI index
% 3rd col = squared Euclidean distance between consecutive beta
% 4th col = squared Euclidean distance between consecutive sigma2
% 5th col = squared Euclidean distance between consecutive sigma2c
% 6th col = squared Euclidean distance between consecutive X centroids
% 7th col = squared Euclidean distance between consecutive X covariance matrices
Amon=[alphaLik(2:end) zeros(lalpha-1,6)];

noisecluster=0;
IDXmin0=IDX<=0;
IDXm=IDX;
IDXm(IDXmin0)=0;
for j=2:lalpha
    
    % Compute ARI index between two consecutive alpha values
    [ARI]=RandIndexFS(IDXm(:,j-1),IDXm(:,j),noisecluster);
    % Store in the second column the ARI index
    Amon(j-1,2)=ARI;
    
    % Compute and store squared euclidean distance between consecutive
    % centroids
    Amon(j-1,3)=sum(sum( (Beta(:,:,j)-Beta(:,:,j-1)).^2, 2)) / sum(sum( Beta(:,:,j-1)).^2, 2);
    
    % Compute and store squared euclidean distance between consecutive
    % sigma2
    Amon(j-1,4)=sum(sum( (Sigma2y(:,j)-Sigma2y(:,j-1)).^2, 2)) /sum(sum( (Sigma2y(:,j-1)).^2, 2));
    
    % Compute and store squared euclidean distance between consecutive
    % sigma2c
    Amon(j-1,5)=sum(sum( (Sigma2yc(:,j)-Sigma2yc(:,j-1)).^2, 2))/ sum(sum( (Sigma2yc(:,j-1)).^2, 2));
    
    % Compute and store squared euclidean distance between consecutive
    % centroids
    Amon(j-1,6)=sum(sum( (MU(:,:,j)-MU(:,:,j-1)).^2, 2)) / sum(sum( (MU(:,:,j-1)).^2, 2));
    
    % Compute and store squared euclidean distance between consecutive
    % covariance matrices (all elements of cov matrices are considered)
    Amon(j-1,7)=sum(sum(sum((SIGMA{j}-SIGMA{j-1}).^2,2)))/ sum(sum(sum((SIGMA{j-1}).^2,2)));
end
out.Amon=Amon;



% alphasel contains the indexes of the columns of matrix IDX which have
% to be plotted
% We use round(alpha*1e+7)/1e+7 to guarantee compatibility with old
% versions of MATLAB. For the new versions the instruction would have
% been:
% [~,alphasel]=intersect(round(alpha,9),alphasel,'stable');
[~,alphasel]=intersect(round(alphaLik*1e+7)/1e+7,round(alphasel*1e+7)/1e+7,'stable');
lalphasel=length(alphasel);


%% Produce all necessary calculations for UnitsTrmOrChgCla plot
IDs=[UnitsTrmOrChgCla IDX(UnitsTrmOrChgCla,:)];
[n1,k1]=size(IDs);
onex=ones(n1,1);
seqIDs=(1:n1)';

Alltrueselj=zeros(n1,1);
ij=0;
for j=1:k
    Nj=IDs(:,2:end)==j;
    sumNj=sum(Nj,2);
    selj=seqIDs(sumNj>0);
    %     if mod(j,2)~=0
    [~,indj]=sort(sumNj(selj),'descend');
    %     else
    %        [~,indj]=sort(sumNj(selj),'ascend');
    %     end
    trueselj=selj(indj);
    trueseljf=setdiff(trueselj,Alltrueselj,'stable');
    lt=length(trueseljf);
    if lt>0
        Alltrueselj(ij+1:ij+lt)=trueseljf;
        ij=ij+lt;
    end
end
% UnitsNeverAssigned

unitsNeverAssigned=setdiff(seqIDs,Alltrueselj);
if ~isempty(unitsNeverAssigned)
    Alltrueselj(ij+1:end)=unitsNeverAssigned;
end

% IDt is the same as IDs but the rows are rearranged in order to
% have units in group 1 at the bottom then units of group 2...
IDt=IDs(Alltrueselj,:);
% IDt(IDt==-1)=[];
out.UnitsTrmOrChgCla=IDt;

%% 1 Monitor change of statistics between two consecutive values of alphaLik
namej='monitor';
d=find(strcmp(namej,name));
if d>0
    figure('Name',namej,'Visible','on');
    
    plotsname={'ARI','$\hat \beta$','$\hat \sigma^2$','$\hat \sigma^2_c$'...
        '$\hat \mu_X$' '$\hat \Sigma_X$'};
    if alphaX==1
        nr=3;
        nc=2;
    else
        nr=2;
        nc=2;
    end
    
    for j=1:nr*nc
        subplot(nr,nc,j)
        plot(Amon(:,1),Amon(:,j+1))
        xlim([min(alphaLik),max(alphaLik)])
        % set(gca,'XTickLabel',num2str(alpha1'))
        
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),lalpha))
        set(gca,'XTickLabel',num2str(flipud(alphaLik)))
        
        set(gca,'XDir','reverse')
        xlabel('Level of trimmming')
        title(plotsname{j},'interpreter','latex')
    end
    sgtitle('Monitor changes between two consecutive alpha values')
end


%% 2 Monitoring stability of classification (units plot)
namej='UnitsTrmOrChgCla';
d=find(strcmp(namej,name));
if d>0
    
    figure('Name',namej,'Visible','on');
    xlim([0 k1+1])
    ylim([0 n1+1])
    set(gca,'XTick',0:k1+1);
    
    for j=1:k1
        strj=cellstr(num2str(IDt(:,j)));
        % Empty spaces for trimmed units
        strj(strcmp(strj,'-1'))={''};
        
        % set colors red for trimmed units blue for those which changed
        % classification
        % UnitsTrm
        h=text(j*onex,seqIDs,strj,'HorizontalAlignment','center');
        if j==1
            [~,indredcolor]=intersect(IDt(:,1),UnitsChgCla);
            col=repmat({'b'},n1,1);
            col(indredcolor)={'r'};
            set(h,{'Color'},col)
        end
    end
    alpha1str=num2str(alphaLik(:));
    newxtcklab=cell(k1+2,1);
    newxtcklab([1:2 k1+2])={''};
    newxtcklab(3:k1+1)=cellstr(alpha1str);
    set(gca,'xticklabels',newxtcklab)
    set(gca,'yticklabels','')
    
    xlabel('Trimming level')
    ylabel('Units trimmed at least once or which changed assignment')
    title('In red color units which changed classification')
    
    if ~isempty(unitsNeverAssigned)
        hline=refline(0,n1-length(unitsNeverAssigned)+0.5);
        hline.Color = 'm';
        title('In red units which changed classification. Above the hor. line unit never classified.')
    else
        title('In red color units which changed classification')
    end
    
end


%% 3 Monitoring posterior probabilities
namej='PostProb';
d=find(strcmp(namej,name));
if d>0
    figure('Name',namej,'Visible','on');
    Prob1=squeeze(Postprob(:,1,:));
    Prob1(IDXmin0)=NaN;
    
    group=cell(n,1);
    group(1:n)={'Units which never changed assignment'};
    group(UnitsChgCla)={'Units which changed assignment'};
    group(UnitsTrm)={'Trimmed units'};
    
    subplot(2,2,1)
    parallelcoords(Prob1,'Group',group, 'Labels',alpha1str)
    ylim([-0.05 1.05])
    xlabel('Level of trimming')
    ylabel('Post prob. group 1 all units')
    set(legend,'Location','best')
    
    subplot(2,2,2)
    Prob1sel=Prob1(UnitsTrmOrChgCla,:);
    groupsel=group(UnitsTrmOrChgCla);
    
    parallelcoords(Prob1sel,'Group',groupsel,...
        'Labels',alpha1str)
    ylim([-0.05 1.05])
    xlabel('Level of trimming')
    ylabel('Post prob. group 1 selected units')
    set(legend,'Location','best')
    % Add the label of the units whose final post prob is intermediate
    unitswithText=Prob1sel(:,end)>0.05 &  Prob1sel(:,end)<0.95;
    text(lalpha*ones(sum(unitswithText),1),Prob1sel(unitswithText,end),...
        cellstr(num2str(UnitsTrmOrChgCla(unitswithText))));
    
    subplot(2,2,3)
    % Prob1table=array2table(Prob1,'VariableNames',cellstr(num2str(alphaLik)));
    parallelplot(Prob1,'GroupData',group)
    set(gca,'CoordinateTickLabels',cellstr(num2str(alphaLik)))
    xlabel('Level of trimming')
    ylabel('Post prob. group 1 all units')
    legend('off')
    
    subplot(2,2,4)
    parallelplot(Prob1(UnitsTrmOrChgCla,:),'GroupData',group(UnitsTrmOrChgCla))
    set(gca,'CoordinateTickLabels',cellstr(num2str(alphaLik)))
    xlabel('Level of trimming')
    ylabel('Post prob. group 1 selected units')
    legend('off')
    sgtitle('Monitor posterior probabilities')
end


%% Monitoring  of sigma2 and sigma2corr
namej='Sigma';
d=find(strcmp(namej,name));
if d>0
    figure('Name',namej,'Visible','on');
    % Sigma2y is of dimension k-by-length(alpha)
    subplot(2,1,1)
    % Sigma2y is k-by-length(alphaLik)
    h=plot(alphaLik,Sigma2y');
    % set the colors using the order in clrdef
    set(h,{'Color'},cellstr(clrdef(1:k)'))
    
    xlim([min(alphaLik),max(alphaLik)])
    % set(gca,'XTickLabel',num2str(alpha1'))
    
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),lalpha))
    set(gca,'XTickLabel',num2str(flipud(alphaLik)))
    
    set(gca,'XDir','reverse')
    xlabel('Level of trimmming')
    ylabel('$\hat \sigma^2$','Interpreter','latex')
    legend(legendGroups)
    legend('hide')
    
    subplot(2,1,2)
    h=plot(alphaLik,Sigma2yc');
    % set the colors using the order in clrdef
    set(h,{'Color'},cellstr(clrdef(1:k)'))
    
    xlim([min(alphaLik),max(alphaLik)])
    % set(gca,'XTickLabel',num2str(alpha1'))
    
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),lalpha))
    set(gca,'XTickLabel',num2str(flipud(alphaLik)))
    
    set(gca,'XDir','reverse')
    xlabel('Level of trimmming')
    ylabel('$\hat \sigma^2_c$','Interpreter','latex')
    legend(legendGroups)
    legend('show')
    clickableMultiLegend(h)
    sgtitle('Monitor error variances')
end


%% Plot scatter with all regression lines (hyperplanes)
namej='ScatterWithRegLines';
d=find(strcmp(namej,name));
if d>0
    
    idx=IDX(:,1);
    % this is just for rotating colors in the plots
    
    intercept=1;
    % The following plots are for the bi-variate case (i.e. v=1)
    if p-intercept < 2
        
        % initialize figure
        figure('Name',namej,'Visible','on');
        hold on;
        alpha1range=['[' num2str(max(alphaLik)) ' \; ' num2str(min(alphaLik)) ']' ];
        title({['$\quad mixt=' num2str(mixt) , '  \quad c_{\hat \sigma^2}='...
            num2str(restrfact) '\quad \alpha_{Lik}=' alpha1range ...
            '\quad \alpha_{\Sigma_X}=' num2str(alphaX) '$']} , ...
            'interpreter' , 'LaTex', 'fontsize' , 14);
        
        % plot regression lines
        vv = [min(X(:,end)) max(X(:,end))];
        
        % hRegLines vector of graphic handles which will contain regression
        % lines
        hRegLines = gobjects(k,1);
        % hText vector of graphic handles which will contain the labels of
        % groups
        hText = gobjects(k,1);
        for jj = 1:k
            group_label = ['Group ' num2str(jj)];
            
            % jj refers to groupj 1,...k
            % jjj refers to trimming
            if intercept==1
                
                for jjj=1:lalpha
                    %   if jjj==1
                    gr=plot(vv,Beta(jj,1,jjj)+Beta(jj,2,jjj)*vv,'DisplayName',[group_label ' fit' ],...
                        'Color',clrdef(jj)); %#ok<NASGU>
                end
                eval(['hRegLines(' num2str(jj) ')=gr;']);
                
            elseif intercept==0
                for jjj=1:lalpha
                    gr=plot(vv,Beta(jj,1,jjj)*vv,'DisplayName',[group_label ' fit' ],...
                        'Color',clrdef(jj)); %#ok<NASGU>
                end
                eval(['hRegLines(' num2str(jj) ')=gr;']);
            end
            
            
            ucg = find(idx==jj);
            % we add a (ficticious) plot instruction with white symbols
            texth=plot(X(ucg,end),y(ucg),'.w','DisplayName',[group_label ' (' num2str(length(ucg)) ' units)']); %#ok<NASGU>
            text(X(ucg,end),y(ucg),num2str(jj*ones(length(ucg),1)),...
                'DisplayName',[group_label ' (' num2str(length(ucg)) ' units)'] , ...
                'HorizontalAlignment','center','VerticalAlignment','middle',...
                'Color',clrdef(jj), 'fontsize' , 12);
            eval(['hText(' num2str(jj) ')=texth;']);
        end
        
        % Plot the outliers (trimmed points)
        ucg = find(idx==-1);
        % hunitsMinus1 = graphical handle to the first level trimmed units
        hunitsMinus1=plot(X(ucg,end),y(ucg),sym1stLevelTrimmedUnits,'color',col1stLevelTrimmedUnits,'MarkerSize',8,...
            'DisplayName',['Trimmed units 1st (' num2str(length(ucg)) ')']);
        
        % Plot second level trimmed units (if there are)
        ucg = find(idx==-2);
        % hunitsMinus2 = graphical handle to second level trimmed units
        hunitsMinus2=plot(X(ucg,end),y(ucg),sym2ndLevelTrimmedUnits,'color',col2ndLevelTrimmedUnits,...
            'DisplayName',['Trimmed units 2nd (' num2str(length(ucg)) ')']);
        
        
        % Add clickable multilegend
        clickableMultiLegend([hRegLines; hText; hunitsMinus1; hunitsMinus2],...
            'Location','best','interpreter' , 'LaTex', 'fontsize' , 10) % ,'TextColor','r');
        axis('manual');
        
        % control of the axis limits
        xmin = min(X(:,end)); xmax = max(X(:,end));
        ymin = min(y); ymax = max(y);
        deltax = (xmax - xmin) / 10;
        deltay = (ymax - ymin) / 10;
        
        xlim([xmin-deltax,xmax+deltax]);
        ylim([ymin-deltay,ymax+deltay]);
        
    else
        % In this case p > 2. A standard yXplot is used.
        
        bunitsMinus1=sum(idx==-1)>0;
        bunitsMinus2=sum(idx==-2)>0;
        if bunitsMinus1 ==true && bunitsMinus2 && true
            clrdefj=[col2ndLevelTrimmedUnits col1stLevelTrimmedUnits clrdef];
            symdefj=[sym2ndLevelTrimmedUnits sym1stLevelTrimmedUnits symdef];
            numsym=k+2;
        elseif  bunitsMinus2 && true
            clrdefj=[col2ndLevelTrimmedUnits clrdef];
            symdefj=[sym2ndLevelTrimmedUnits symdef];
            numsym=k+1;
        elseif bunitsMinus1 ==true
            clrdefj=[col1stLevelTrimmedUnits clrdef];
            symdefj=[sym1stLevelTrimmedUnits symdef];
            numsym=k+1;
        else
            clrdefj=clrdef;
            symdefj=symdef;
            numsym=k;
        end
        
        plo=struct;
        plo.clr = clrdefj(1:numsym);
        plo.sym = symdefj(1:numsym);
        
        % group names in the legend
        group = cell(n,1);
        group(idx==-1) = {'Trimmed units'};
        group(idx==-2) = {'Trimmed units level 2'};
        for iii = 1:k
            group(idx==iii) = {['Group ' num2str(iii)]};
        end
        
        % yXplot
        % Remark: it is necessary to sort idx because in this way idx=-2 (if present) is
        % the first symbol, idx=-1 is the second ....
        [~,indsor]=sort(idx);
        [~,AX,~]=yXplot(y(indsor),X(indsor,:),group(indsor),plo);
        % Dimension of Beta is k-by-p-by-length(alpha)
        for j = 1:length(AX)
            % Make the axes of the panel with handle AX(i) the current axes.
            set(gcf,'CurrentAxes',AX(j));
            
            for i=1:k
                % refline slope and intercept
                % add a line for each value of alpha1
                for jj=1:length(alphaLik)
                    hline=refline([Beta(i,j+1,jj) Beta(i,1,jj)]);
                    hline.Color=clrdef(i);
                end
            end
            
        end
    end
end

%% Monitoring of allocation (using gscatter)
d=find(strcmp('gscatter',name));
if d>0
    
    % Monitoring of allocation
    if  lalphasel==1
        nr=1;
        nc=1;
    elseif lalphasel==2
        nr=2;
        nc=1;
    elseif lalphasel<=4
        nr=2;
        nc=2;
    elseif lalphasel<=6
        nr=3;
        nc=2;
    elseif lalphasel<=9
        nr=3;
        nc=3;
    elseif lalphasel<=12
        nr=3;
        nc=4;
    else
        nr=4;
        nc=4;
    end
    
    resup=1;
    figure('Name',['Monitoring allocation #' int2str(resup)])
    
    
    
    jk=1;
    for j=1:lalphasel
        
        % The monitoring must contain a maximum of 16 panels
        % If length(alpha) is greater than 16 a new set of 16 subpanels is
        if jk>16
            jk=1;
            resup=resup+1;
            figure('Name',['Monitoring allocation #' int2str(resup)])
            subplot(nr,nc,jk)
        else
            subplot(nr,nc,jk)
        end
        jk=jk+1;
        
        idxselj=IDX(:,alphasel(j));
        % Check if inside idxselj there are 1st or 2nd level trimmed units
        bunitsMinus1=sum(idxselj==-1)>0;
        bunitsMinus2=sum(idxselj==-2)>0;
        if bunitsMinus1 ==true && bunitsMinus2 && true
            clrdefj=[col2ndLevelTrimmedUnits col1stLevelTrimmedUnits clrdef];
            symdefj=[sym2ndLevelTrimmedUnits sym1stLevelTrimmedUnits symdef];
        elseif  bunitsMinus2 && true
            clrdefj=[col2ndLevelTrimmedUnits clrdef];
            symdefj=[sym2ndLevelTrimmedUnits symdef];
        elseif bunitsMinus1 ==true
            clrdefj=[col1stLevelTrimmedUnits clrdef];
            symdefj=[sym1stLevelTrimmedUnits symdef];
        else
            clrdefj=clrdef;
            symdefj=symdef;
        end
        
        if p>2
            % In this case PLS regression is used.
            % In order to take into account group structure k-1 dummy
            % variables are added to matrix X. Of course trimmed units for
            % that particular value of alphaLik are excluded, but included
            % in the gscatter plots.
            DUM=zeros(n,k-1);
            for jj=1:k-1
                DUM(idxselj==jj,jj)=1;
            end
            Xext=[X DUM];
            idxgt0=idxselj>0;
            % training set
            ysel0=y(idxgt0);
            Xsel0=Xext(idxgt0,:);
            Xsel1=Xext(~idxgt0,:);
            mXsel0=mean(Xsel0);
            [~,~,XS0,~,~,PCTVAR,~,stats] = plsregress(Xsel0,ysel0,1);
            % [XL,yl,XS0,YS0,beta,PCTVAR,MSE,stats] = plsregress(Xsel0,ysel0,1);
            
            % Find best predictor for non trimmed (XS0) and trimmed units
            XS1=zeros(n,1);
            XS1(idxgt0)=XS0;
            % REMARK
            % XS0chk=(Xsel0-mXsel0)*stats.W;
            XS1(~idxgt0)=(Xsel1-mXsel0)*stats.W;
            
            hh=gscatter(XS1,y,idxselj,clrdefj,symdefj);
            
            xlabel('PLS predictor')
            ylabel('y')
            
            clickableMultiLegend(hh)
            if jk>2
                legend hide
            end
            axis manual
            alphajtxt=num2str(alphaLik(alphasel(j)));
            title(['$\alpha=$' alphajtxt 'Var. \; explained=' num2str(100*PCTVAR(2,1),3) ],'Interpreter','Latex')
            
        elseif p==2
            hh=gscatter(X(:,end),y,idxselj,clrdefj,symdefj);
            
            xlabel('x1')
            ylabel('y')
            
            clickableMultiLegend(hh)
            if jk>2
                legend hide
            end
            axis manual
            title(['$\alpha=$' num2str(alphaLik(alphasel(j)))],'Interpreter','Latex')
        else
            % Univariate case: plot the histogram
            histFS(y,10,idxselj,[],[],clrdefj)
            title(['$\alpha=$' num2str(alphaLik(alphasel(j)))],'Interpreter','Latex')
        end
    end
end

%% Monitoring of beta regression coefficients (standardized)
namej='Beta';
d=find(strcmp(namej,name));
if d>0
    figure('Name',namej,'Visible','on');
    % Dimension of Beta is k-by-p-by-length(alpha)
    % If p is 2 first column contains interecepts and second slopes
    % Standardize Beta along the 3rd dimension
    Betast=zscore(Beta,0,3);
    if  p==1
        nr=1;
        nc=1;
    elseif p==2
        nr=2;
        nc=1;
    elseif p<=4
        nr=2;
        nc=2;
    elseif p<=6
        nr=3;
        nc=2;
    elseif p<=9
        nr=3;
        nc=3;
    elseif p<=12
        nr=3;
        nc=4;
    else
        nr=4;
        nc=4;
    end
    
    for j=1:p
        subplot(nr,nc,j)
        h=plot(alphaLik(:),squeeze(Betast(:,j,:))');
        % set the colors using the order in clrdef
        set(h,{'Color'},cellstr(clrdef(1:k)'))
        
        xlim([min(alphaLik),max(alphaLik)])
        % set(gca,'XTickLabel',num2str(alpha1'))
        
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),lalpha))
        set(gca,'XTickLabel',num2str(flipud(alphaLik)))
        
        set(gca,'XDir','reverse')
        xlabel('Level of trimmming')
        ylabel(['$\hat \beta_' num2str(j-1) '$'],'Interpreter','latex')
        legend(legendGroups)
        legend('hide')
        if j==p
            legend('show')
            clickableMultiLegend(h)
        end
    end
    sgtitle('Monitor estimated regression coefficients')
end

%% Monitor group size
namej='Siz';
d=find(strcmp(namej,name));
if d>0
    % Monitoring of group size
    figure('Name',namej,'Visible','on');
    h=plot(alphaLik(:),out.Nopt');
    % set the colors using the order in clrdef
    set(h,{'Color'},cellstr(clrdef(1:k)'))
    
    xlim([min(alphaLik),max(alphaLik)])
    % set(gca,'XTickLabel',num2str(alpha1'))
    
    lalpha=length(alphaLik);
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),lalpha))
    set(gca,'XTickLabel',num2str(flipud(alphaLik)))
    
    set(gca,'XDir','reverse')
    title('Monitor group size')
    xlabel('Level of trimmming')
    legend(legendGroups)
    legend('show')
    clickableMultiLegend(h)
    
end

end
%FScategory:CLUS-RobClaMULT