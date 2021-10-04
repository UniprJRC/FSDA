function [out, varargout] = tclustregeda(y,X,k,restrfact,alphaLik,alphaX,varargin)
%tclustregeda performs robust linear grouping analysis for a series of values of the trimming factor
%
%<a href="matlab: docsearchFS('tclustregeda')">Link to the help function</a>
%
%   tclustregeda performs tclustreg for a series of values of the trimming
%   factor alpha, for given k (number of groups), restrfactor (restriction
%   factor) and alphaX (second level trimming or cluster weighted model).
%   In order to increase the speed of the computations, parfor is used.
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
%               sized groups by maximizing the likelihood. Default value
%               false.
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
%            be conveniently generated by function subsets.
%            nsamp must have k*p columns. The first p columns are used to
%            estimate the regression coefficient of group 1... the last p
%            columns are used to estimate the regression coefficient of
%            group k.
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
% commonslope  : Impose constraint of common slope regression coefficients. Boolean.
%               If commonslope is true, the groups are forced to have the
%               same regression coefficients (apart from the intercepts).
%               The default value of commonslope is false;
%                 Example - 'commonslope',true
%                 Data Types - boolean
%
%
% plots    :    Plot on the screen. Scalar structure.
%
%               Case 1: plots option used as scalar.
%               - If plots=0, plots are not generated.
%               - If plots=1 (default), 6 plots are shown on the screen.
%                 The first plot ("monitor") shows 4 or 6 panels
%                 monitoring between two consecutive values of alpha the
%                 change in classification using ARI index (top left panel),
%                 the relative change in beta (top right panel)
%                 ($||\beta_{\alpha_1}-\beta_{\alpha_2}||^2/||\beta_{\alpha_2}||^2
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
%                           two consecutive vectors of variances of the
%                           errors of the k regressions.
%                       5th col = relative squared Euclidean distance between
%                           two consecutive vectors of correct variances of the
%                           errors of the k regressions.
%                       6th col = relative squared Euclidean distance between
%                           two consecutive $\hat \mu_X$.
%                       7th col = relative squared Euclidean distance between
%                           two consecutive $\hat \Sigma_X$.
%
%            out.MU  =  3D array of size k-by-(p-1)-by-length(alphaLik) containing
%                       the monitoring of the X centroids for each value of
%                       alphaLik. out.MU(:,:,1), refers to alphaLik(1) ...,
%                       out.MU(:,:,end) refers to alphaLik(end). First row in
%                       each slice refers to group 1, second row refers to
%                       group 2 ... This field is present only if input
%                       option alphaX is 1.
%
%         out.SIGMA  =  cell of length length(alphaLik) containing in element
%                       j, with j=1, 2, ...,  length(alphaLik), the 3D array
%                       of size (p-1)-by-(p-1)-by-k containing the k (constrained)
%                       estimated covariance matrices of X associated with
%                       alphaLik(j). This field is present only if input
%                       option alphaX is 1.
%
% out.UnitsTrmOrChgCla = Matrix containing information about the
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
%   out.Postprob      = Posterior probabilities. 3D array of size
%                       n-by-k-length(alphaLik) containing the monitoring
%                       of posterior probabilities.
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
% This procedure extends to tclustreg the so-called monitoring approach.
% The philosophy is to investigate how the results change as the trimming
% proportion alpha reduces. This function enables us to monitor the change
% in classification (measured by the ARI index) and the change in
% regression coefficients and error variances (measured by the relative
% squared Euclidean distances).
% Note that there is a non-identifiability of a finite mixture distribution
% caused by the invariance of the mixture density function to components
% relabeling. In order to make sure that consistent labels are used for the
% groups, between two consecutive values of alpha, we call
% we assign label r to a group if this group shows the smallest distance
% with group r for the previous value of alpha.
% More precisely, once the labelling is fixed for the largest value of the
% trimming factor supplied, if $\alpha_r$ and $\alpha_s$ denote two
% consecutive levels of trimming ($\alpha_r>\alpha_s$) and
% \beta_{j,\alpha_r}
%given \beta_{j,\alpha_r}
%
% In order to be consistent
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
% Copyright 2008-2021.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('tclustregeda')">Link to the help function</a>
%
%$LastChangedDate:: 2018-01-29 18:52:14 #$: Date of the last commit

% Examples:

%{
    % tclustregeda of contaminated X data using all default options.
    % The X data have been introduced by Gordaliza, Garcia-Escudero and Mayo-Iscar (2013).
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
    %% tclustregeda with a noise variable and personalized plots.
    % Use the X data of the previous example.
    %  set(0,'DefaultFigureWindowStyle','docked')
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
    % Personalize plots. Just show the gscatter plot. In this case given
    % that there is more than one explanatory variable  PLS regression
    % (adding the dummies for the classified units) is performed. In the
    % gscatter plots the percentage of variance explained by the first
    % linear combination of the X variables is given in the title of each
    % panel of the gscatter plot.
    plots=struct;
    plots.name={'gscatter'};
    out = tclustregeda(y,X,k,restrfact,alphaLik,alphaX,'plots',plots);
%}

%{
    % tclustregeda: example with multiple explanatory variables, with yXplot.
    X  = load('X.txt');
    y = X(:,end);
    rng(100)
    X = [randn(length(y),4) X(:,1:end-1)];
    k=2;
    alphaLik = [0.10:-0.01:0]' ;
    alphaX = 0;
    restrfact =1;
    mixt=2;
    plots=struct;
    % plots.name={'gscatter','postprob'};
    % plots.alphasel=[0.05 0.03 0];
    plots.name={'all'};
    % 'UnitsSameGroup',152,
    out = tclustregeda(y,X,k,restrfact,alphaLik,alphaX,'mixt',2,'msg',0,'plots',plots);
%}

%% Beginning of code
% Control variables, tolerances and internal flags
warning('off');

verbertotest = 9.5; % R2018b
vafter91=~verLessThanFS(verbertotest); % >=2018b

verbertotest = 9.6; % R2019a
vafter95=~verLessThanFS(verbertotest); % >=2019a

scrsz = get(groot,'ScreenSize');
[left , bottom] = deal(scrsz(1) , scrsz(2));%

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
    error('FSDA:tclustregeda:error','error must a scalar in the interval [0 0.5] or an integer specifying the number of units to trim')
end

% checks on cwm, which decides if clusterwise regression has to be used
if alphaX < 0 || alphaX >1
    error('FSDA:tclustregeda:WrongAlphaX','alphaX must a scalar in the interval [0 1]')
elseif alphaX==1
    cwm = 1;
else
    cwm = 0;
end

%% User options and their default values

%%% - nsamp: the number of subsets to extract randomly, or the indexes of the initial subsets pre-specified by the User
numpool   = feature('numCores');
cleanpool = false;

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
            nselected      = nsampdef;
            
            % Flag indicating if the user has selected a prior subset
            NoPriorSubsets = 0;
            
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
            %             if sum(chknwtrim)>0 && ~isstruct(cell2mat(varargin(find(chknwtrim)+1)))
            %                 if cell2mat(varargin(find(chknwtrim)+1))== 4
            %                     for f=1:size(C,1)*size(C,2)
            %                         if Wt4(C(f)) == 0
            %                             [~,C(f)] = min(pdist2([yori(C(f)),Xori(C(f))],[yori(Wt4),Xori(Wt4)]));
            %                         else
            %                             C(f) = sum(Wt4(1:C(f)));
            %                         end
            %                     end
            %                 end
            %             end
            
        else
            % If nsamp is a scalar it simply contains the number of subsets
            % which have to be extracted. In this case NoPriorSubsets=1
            NoPriorSubsets = 1;
        end
    else
        % If option nsamp is not supplied, then there are no prior subsets
        NoPriorSubsets = 1;
    end
else
    % if nargin == 6, then the user has not supplied prior subsets
    NoPriorSubsets = 1;
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
reftoldef    = 1e-5;

% default value for we: the observation weights
wedef        = ones(n,1);

% default model: classification (mixt=0) or mixture likelihood (mixt=2)
mixtdef      = 0;

% default choice for equalweight constraint
equalweightsdef = 0;

%seqk = sequence from 1 to the number of groups
seqk = 1:k;

plots          = 1;
UnitsSameGroup = '';
commonslopedef=false;

% automatic extraction of user options
options = struct('intercept',1,'mixt',mixtdef,...
    'nsamp',nsampdef,'refsteps',refstepsdef,...
    'reftol',reftoldef,'commonslope',commonslopedef,...
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
numpool = options.numpool;

% global variable controlling if messages are displayed in the console.
msg     = options.msg;

% Graphs summarizing the results
plots   = options.plots;

% Number of subsets to extract or matrix containing the subsets
nsamp   = options.nsamp;

% Concentration steps
refsteps = options.refsteps;
reftol   = options.reftol;

% Common slope constraint
commonslope=options.commonslope;


% Equalweights constraints
equalweights   = options.equalweights;

UnitsSameGroup = options.UnitsSameGroup;

% application-specific weights vector assigned by the user for beta
% estimation
we             = options.we;


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
    NoPriorNini = 1;
else
    NoPriorNini = 0;
end

%% Initializations

%%% - Subsets extraction

%case with no prior subsets from the User
if NoPriorSubsets
    
    [C,nselected] = subsets(nsamp, n, (p+intercept)*k, ncomb, 0, we);
    
    % C = matrix of indexes of the subsets to extract for the k groups
    
    %nselected is set equal to the number of subsets:
    % - nselected = nsamp, if nsamp is a scalar;
    % - nselected = size(nsamp,1), if nsamp is ther matrix of the initial subsets
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
wtype_beta = 0;
wtype_obj  = 0;
zigzag     = (alphaX > 0 && alphaX<1) || wtype_beta == 3 || ...
    wtype_beta == 2 || ~strcmp(wtype_obj, '0');

% Make sure alphaLik is a column vector;
alphaLik = alphaLik(:);

lalpha   = length(alphaLik);
if msg == 1
    progbar = ProgressBar(lalpha);
else
    progbar=[];
end

IDX     = zeros(n,lalpha);

% Check for the existence of the X space
% p is number of explanatory variables including intercept (if present)
if ((p-1 > 0 && intercept==0) || (p>1 && intercept==1)) && cwm==1
    existXspace=true;
else
    existXspace=false;
end

% Note that Mu and SIGMA are meaningful just if existXspace=true;
MU      = zeros(k,p-intercept,lalpha);
SIGMA   = cell(lalpha,1);

Beta    = zeros(k,p,lalpha);
Nopt    = zeros(k,lalpha);
Sigma2y = zeros(k,lalpha);
Sigma2yc= Sigma2y;
Vopt    = zeros(lalpha,1);
Postprob=zeros(n,k,lalpha);

% Do not show messages during each execution of tclustregcore
msgrs=0;

internationaltrade=false;
if internationaltrade==true
    wei=X(:,end).^2/var(X(:,end))+y.^2/var(y);
    weiForLikComputation=wei/max(wei);
else
    weiForLikComputation=1;
end

parfor (j=1:lalpha, numpool)
    % for j=1:lalpha
    
    h = floor(n*(1-alphaLik(j)));
    
    %%  RANDOM STARTS
    
    [bopt,sigma2opt,nopt,postprobopt,muXopt,sigmaXopt,vopt,~,idxopt]...
        = tclustregcore(y,X,RandNumbForNini,reftol,refsteps,mixt,...
        equalweights,h,nselected,k,restrfact,restrfactX,alphaLik(j),alphaX,...
        seqk,NoPriorNini,msgrs,C,intercept,cwm,commonslope,wtype_beta,we,wtype_obj,zigzag,weiForLikComputation);
    
    %%  END OF RANDOM STARTS
    
    IDX(:,j)    = idxopt;
    Beta(:,:,j) = bopt';
    Nopt(:,j)   = nopt;
    
    if existXspace == true
        MU(:,:,j)   = muXopt;
        SIGMA{j}    = sigmaXopt;
    end
    
    Vopt(j)     = vopt;
    Postprob(:,:,j) = postprobopt;
    
    % nnto = number of non trimmed observations, after first and second level trimming
    nnto = sum(nopt);
    
    % vt = variance of the truncated normal distribution
    % 1-nnto/n is the trimming percentage
    vt = norminv(0.5*(1+nnto/n));
    
    sigma2opt=sigma2opt';
    if nnto<n
        factor = 1/sqrt(1-2*(n/nnto)*vt.*normpdf(vt));
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
    idx1        = IDX(:,1);
    trimmed1    = idx1==-1;
    trimmed2    = idx1==-2;
    trimmed     = idx1<0;
    idxtmp      = idx1;
    idxtmp(trimmed) = 0;
    [IDXnew1, OldNewIndexes]=ClusterRelabel({idxtmp}, UnitsSameGroup);
    
    if existXspace == true
        % MUold1 is k-by-p (rows refer to groups)
        % Mu is k-by-p-length(alpha)
        MUold1      = MU(:,:,1);
        % Sigmaold1 is pxpxk
        SIGMAold1   = SIGMA{1};
    end
    
    Betaold1    = Beta(:,:,1);
    Sigma2yold1 = Sigma2y(:,1);
    Sigma2ycold1= Sigma2yc(:,1);
    Postprobold1= Postprob(:,:,1);
    
    MUnew1      = MUold1;
    SIGMAnew1   = SIGMAold1;
    Betanew1    = Betaold1;
    Sigma2ynew1 = Sigma2yold1;
    Sigma2ycnew1= Sigma2ycold1;
    Postprobnew1= Postprobold1;
    
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
    IDX(:,1)        = IDXnew1{:};
    IDX(trimmed1,1) = -1;
    IDX(trimmed2,1) = -2;
    
    if existXspace == true
        MU(:,:,1)       = MUold1;
        SIGMA{1}        = SIGMAold1;
    end
    
    Beta(:,:,1)     = Betaold1;
    Sigma2y(:,1)    = Sigma2yold1;
    Sigma2yc(:,1)   = Sigma2ycold1;
    Postprob(:,:,1) = Postprobold1;
end

IDXold = IDX;
maxdist= zeros(lalpha,1);
seqk   = (1:k)';

for j=2:lalpha
    newlab=zeros(k,1);
    mindist=newlab;
    for ii=1:k
        % muii= regression coefficients of group ii for previous alpha value
        muii=Beta(ii,:,j-1);
        
        % Beta(:,:,j) =matrix of regression coefficients for current alpha value
        muij=bsxfun(@minus,muii,Beta(:,:,j));
        
        % It is better to consider relative changes in the coefficients
        % mu_ij is a matrix of size kxp which will contain in row i, i=1, 2, ..., k
        % (\beta_previoustrimming-\beta_currenttrimming)/\beta_previoustrimming
        % We will then take the sum of squares of the elements of each row
        muij = bsxfun(@rdivide, muij, muii);
        
        % For example if indmin is equal to 2 it means that probabbly
        % previous labelled group ii now has ben labelled indmin
        [mind,indmin]=min(sum(muij.^2,2));
        newlab(ii)=indmin;
        mindist(ii)=mind;
    end
    % Store maximum among minimum distances
    [maxmindist,indmaxdist] =max(mindist);
    maxdist(j)=maxmindist;
    
    % If for example in the case k=3, newlab is 2 3 1, it means that
    % previous group 1 now has been labelled group 2, previous group 3 now
    % has been labelled group 3 and previous group 3 now has been labelled
    % group 1.
    % Therefore if isequal(sort(newlab),seqk) relabelling is easy
    if isequal(sort(newlab),seqk)
        
        if existXspace == true
            MU(:,:,j)=MU(newlab,:,j);
            SIGMA(j)= {SIGMA{j}(:,:,newlab)};
        end
        
        Beta(:,:,j)=Beta(newlab,:,j);
        Sigma2y(:,j)=Sigma2y(newlab,j);
        Sigma2yc(:,j)=Sigma2yc(newlab,j);
        Postprob(:,:,j)=Postprob(:,newlab,j);
        
        for r=1:k
            IDX(IDXold(:,j)==newlab(r),j)=r;
        end
    else
        % In this case new1 contains the labels which never appeared inside
        % newlab. To this label we assign the maximum distance and check is
        % this time sequal(sort(newlab),seqk), that is we check whether
        % vector sort(newlab) of length k contain the numbers 1, 2, ..., k
        % if length(newl) >1 two labels do not have the correspondence
        % therefore automatic relabelling is not possible.
        newl=setdiff(seqk,newlab);
        if length(newl)==1
            newlab(indmaxdist)=newl;
            
            if isequal(sort(newlab),seqk)
                if existXspace == true
                    MU(:,:,j)=MU(newlab,:,j);
                    SIGMA(j)= {SIGMA{j}(:,:,newlab)};
                end
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
        else
            disp(['Automatic relabelling not possible when alpha=' num2str(alphaLik(j))])
        end
    end
end


%%  Set the output structure

out         = struct;
out.IDX     = IDX;
out.Sigma2y = Sigma2y;
out.Sigma2yc= Sigma2yc;
out.Beta    = Beta;
out.Nopt    = Nopt;
out.MU      = MU;
out.SIGMA   = SIGMA;
out.Vopt    = Vopt;
out.Postprob= Postprob;

% Store the indices in varargout
if nargout==2
    varargout={C};
end


%% Generate plots

% default for width of lines;
figureResize        = 1.5;
plotLineWidth       = 1.5;
plotLineWidthGrad   = 0.5;
xyTickFontSize      = 12;
xxTickAngleVal      = 45;
xyLabelSize         = 16;
yLabelLatexSize     = 18;
legendSize          = 14;
titleSize           = 18;
subtitleSize        = 16;
matrixFont          = 8;

% Colormaps used to give a gradient to the lines obtained for different
% \alpha values. Colormap follows our color rotation standard:
% clrdef = 'bkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcr';
clrdefmap = zeros(lalpha,3,6);
clrdefmap(:,:,1) = flipud(cmapFS(FSColors.darkblue.RGB ,FSColors.blueish.RGB  ,lalpha));
clrdefmap(:,:,2) = flipud(cmapFS(FSColors.black.RGB    ,FSColors.greysh.RGB   ,lalpha));
clrdefmap(:,:,3) = flipud(cmapFS(FSColors.darkpurpl.RGB,FSColors.purplish.RGB ,lalpha));
clrdefmap(:,:,4) = flipud(cmapFS(FSColors.darkgreen.RGB,FSColors.greenish.RGB ,lalpha));
clrdefmap(:,:,5) = flipud(cmapFS(FSColors.darkcyan.RGB ,FSColors.lightblue.RGB,lalpha));
clrdefmap(:,:,6) = flipud(cmapFS(FSColors.darkred.RGB  ,FSColors.reddish.RGB  ,lalpha));
clrdefmap        = repmat(clrdefmap,1,1,4);

% plotdef = list of the plots which are produced by default (is plots=1)
plotdef={'monitor'; 'UnitsTrmOrChgCla'; 'PostProb'; 'Sigma';...
    'ScatterWithRegLines'; 'gscatter'};
% plotall = list of all available plots
plotall={'monitor'; 'UnitsTrmOrChgCla'; 'PostProb'; 'Sigma';...
    'ScatterWithRegLines'; 'gscatter'; ...
    'Beta'; 'Siz'};

clrdef = 'bkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcr';
symdef = '+sd^v><phos+*d^v><phos+*d^v><phos+*d^v><phos';
linedef = {'-','--',':','-.'};
linedef=repmat(linedef,1,5);

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
else
    % NO PLOT IS SHOWN ON THE SCREEEN
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
legendGroups=cellstr([repmat('Group ',k,1) num2str((1:k)')]);

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
Amon=[alphaLik(2:end) NaN(lalpha-1,6)];

noisecluster=0;
IDXmin0=IDX<=0;
IDXm=IDX;
IDXm(IDXmin0)=0;
for j=2:lalpha
    
    % Make sure there was convergence (i.e. no missing values in
    % IDXm(:,j-1) and IDXm(:,j)).
    if any(isnan(IDXm(:,j-1))) || any(isnan(IDXm(:,j)))
    else
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
        
        if existXspace == true
            % Compute and store squared euclidean distance between consecutive
            % centroids
            Amon(j-1,6)=sum(sum( (MU(:,:,j)-MU(:,:,j-1)).^2, 2)) / sum(sum( (MU(:,:,j-1)).^2, 2));
            
            % Compute and store squared euclidean distance between consecutive
            % covariance matrices (all elements of cov matrices are considered)
            Amon(j-1,7)=sum(sum(sum((SIGMA{j}-SIGMA{j-1}).^2,2)))/ sum(sum(sum((SIGMA{j-1}).^2,2)));
        else
            Amon(j-1,6:7)=[NaN NaN];
        end
    end
end
out.Amon=Amon;

if (isnumeric(plots) && plots~=0) || isstruct(plots)
    % alphasel contains the indexes of the columns of matrix IDX which have
    % to be plotted. We use round(alpha*1e+7)/1e+7 to guarantee
    % compatibility with old versions of MATLAB. For the new versions the
    % instruction would have been:
    % [~,alphasel]=intersect(round(alpha,9),alphasel,'stable');
    [~,alphasel]=intersect(round(alphaLik*1e+7)/1e+7,round(alphasel*1e+7)/1e+7,'stable');
    lalphasel=length(alphasel);
end

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

alpha1str=num2str(alphaLik(:));

%% Monitoring plots follow

if (isnumeric(plots) && plots~=0) || isstruct(plots)
    
    %% 1 Monitor change of statistics between two consecutive values of alphaLik
    
    namej = 'monitor';
    tit   = {'Tclustreg monitoring plot' , 'Changes between two consecutive \alpha-values'};
    tit2  = {'Tclustreg monitoring plot -- Changes between two consecutive $\alpha$-values' , ' '};
    d=find(strcmp(namej,name));
    if d>0
        hf1 = figure('Name',namej,'Visible','off');
        set(hf1,'Tag','tclusteda');
        plotsname={'ARI','$\hat \beta$','$\hat \sigma^2$','$\hat \sigma^2_c$'...
            '$\hat \mu_X$' '$\hat \Sigma_X$'};
        
        if existXspace == true
            nr=3;
            nc=2;
        else
            nr=2;
            nc=2;
        end
        
        for j=1:nr*nc
            subplot(nr,nc,j);
            plot(Amon(:,1),Amon(:,j+1),'LineWidth',plotLineWidth);
            set(gca,'xticklabel',[],'XGrid','on');
            xlim([min(alphaLik),max(alphaLik)]);
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),lalpha), 'FontSize' , xyTickFontSize);
            if j>nc*(nr-1)
                xtickangle(gca,xxTickAngleVal);
                set(gca,'XTickLabel',num2str(flipud(alphaLik)), 'FontSize' , xyTickFontSize);
                xlabel('Level of trimmming', 'FontSize' , xyLabelSize);
            end
            % sigma: plots 3,4,6
            set(gca,'XDir','reverse');
            title(plotsname{j},'interpreter','latex', 'FontSize' , subtitleSize);
            if vafter95
                axtoolbar('Visible','off');
            end
        end
        if vafter95 == true
            sgtitle(tit , 'FontSize' , titleSize);
        else
            a  = axes;
            t1 = title(tit2, 'FontSize' , titleSize, 'FontWeight', 'normal', 'Interpreter' , 'latex');
            a.Visible = 'off'; % set(a,'Visible','off');
            t1.Visible = 'on'; % set(t1,'Visible','on');
        end
        
    end
    
    %% 2 Monitoring stability of classification (units plot)
    
    namej='UnitsTrmOrChgCla';
    tit1 = {'Tclustreg monitoring plot' , 'Stability of classification - changes in purple' , 'Above the horizontal line unit never classified'};
    tit2 = {'Tclustreg monitoring plot' , 'Stability of classification - changes in purple'};
    
    symbseq = {'$\clubsuit$' , '$\diamondsuit$' , '$\heartsuit$' , ...
        '$\spadesuit$' , '$\circ$' , '$\bullet$' , ...
        '$\nabla$' , '$\o$' , '$\copyright$' , '$*$' , '$+$'};
    
    d=find(strcmp(namej,name));
    if d>0
        
        % permutation criterion
        criterion = 0; % default is 0
        switch criterion
            case 0
                IDtt = IDt;
            case 1
                IDtt = topkrows(IDt,size(IDt,1),2:size(IDt,2),'descend');
            case 2
                [~,indChgCla]=intersect(IDt(:,1),UnitsChgCla);
                IDta  = IDt(indChgCla,:);
                IDtb  = IDt(setdiff(1:size(IDt,1),indChgCla),:);
                IDtta = topkrows(IDta,size(IDta,1),2:size(IDta,2),'descend');
                IDttb = topkrows(IDtb,size(IDtb,1),2:size(IDtb,2),'descend');
                IDtt  = [IDtta ; IDttb];
            otherwise
                IDtt = IDt;
        end
        
        % Generate figure.
        hf2=figure('CreateFcn',{@movegui, 'south'},'Name',namej,'Visible','off');
        
        % If the reclassified units are more than thsz, resize the figure
        sz = size(IDtt);
        thsz = 50;
        if sz(1)>thsz
            current_pos = get(hf2,'Position');
            set(hf2,'Position',[left,bottom,figureResize,figureResize].*current_pos);
        end
        set(hf2,'Tag','tclusteda','Resize','off');
        
        % x and y limits
        xlim([0 k1+1]);
        ylim([0 n1+1]);
        set(gca,'XTick',0:k1+1, 'FontSize' , xyTickFontSize);
        
        % fill matrix
        for j=1:k1
            if j==1
                % unit number
                strj=cellstr(num2str(IDtt(:,j)));
                if sz(1)>thsz
                    xpos = onex-(mod(seqIDs,2)*0.5);
                else
                    xpos = onex;
                end
                h = text(xpos,seqIDs,strj,...
                    'HorizontalAlignment','center', ...
                    'FontSize' , matrixFont+4 , 'interpreter','none');
                
                [~,indredcolor]=intersect(IDtt(:,1),UnitsChgCla);
                %col=repmat({'b'},n1,1);
                %col(indredcolor)={'r'};
                %col=repmat({FSColors.darkgrey.RGB},n1,1);
                %col(indredcolor)={FSColors.purplish.RGB};
                col=repmat({[0.4660 0.6740 0.1880]},n1,1);
                col(indredcolor)={[0.6350 0.0780 0.1840]};
                st=repmat({'normal'},n1,1);
                st(indredcolor)={'normal'};
                
                set(h,{'Color'},col);
                set(h,{'FontWeight'},st);
            else
                % classes associated to units
                ipos = find(IDtt(:,j)>0);
                ineg = IDtt(:,j)<=0;
                strj = cell(size(strj));
                strj(ipos) = cellstr(symbseq(abs(IDtt(ipos,j))));
                % Empty spaces for trimmed units
                strj(ineg) = {''};
                hg = text(j*onex,seqIDs,strj,'HorizontalAlignment','center', 'FontSize' , matrixFont,'interpreter','latex');
                colorsG = cellstr(repmat('w',size(strj)));
                colorsG(ipos) = cellstr(clrdef(abs(IDtt(ipos,j)))');
                set(hg,{'Color'},colorsG);
            end
            
        end
        
        newxtcklab=cell(k1+2,1);
        newxtcklab([1:2 k1+2])={''};
        newxtcklab(3:k1+1)=cellstr(alpha1str);
        set(gca,'xticklabels',newxtcklab, 'FontSize' , xyTickFontSize)
        set(gca,'yticklabels','', 'FontSize' , xyTickFontSize)
        
        xlabel({'Trimming level'}, 'FontSize' , xyLabelSize);
        ylabel({'Units trimmed at least once' , 'or which changed assignment'}, 'FontSize' , xyLabelSize);
        
        % Group indication in title
        un=unique(IDtt(:,2:end));
        un=un(un>0);
        Gleg = strcat('G' , string(num2cell(un)) ,  '=' , symbseq(un)');
        tit1(size(tit1,2)+1) = cellstr(strjoin(Gleg));
        tit2(size(tit2,2)+1) = cellstr(strjoin(Gleg));
        
        % Add title
        if ~isempty(unitsNeverAssigned)
            hline=refline(0,n1-length(unitsNeverAssigned)+0.5);
            hline.Color = 'm';
            title(tit1 , 'FontSize' , titleSize , 'FontWeight', 'normal','interpreter','latex');
        else
            title(tit2 , 'FontSize' , titleSize , 'FontWeight', 'normal','interpreter','latex');
        end
        
        pan('off');
        if vafter95
            axtoolbar(gca,'Visible','off');
        end
        %movegui(gcf,'south');
    end
    
    %% 3 Monitoring  of sigma2 and sigma2corr
    namej = 'Sigma';
    tit   = {'Tclustreg monitoring plot' , ['Error variances for restriction factor c=' num2str(restrfact)]};
    
    d=find(strcmp(namej,name));
    if d>0
        hf4 = figure('Name',namej,'Visible','off');
        set(hf4,'Tag','tclusteda');
        
        % first subplot
        subplot(2,1,1);
        % Sigma2y is k-by-length(alphaLik)
        h1  = plot(alphaLik,Sigma2y','LineWidth',plotLineWidth);
        % set the colors and linestyle
        set(h1,{'Color'},cellstr(clrdef(1:k)'),{'LineStyle'},linedef(1:k)',{'DisplayName'},legendGroups);
        
        xlim([min(alphaLik),max(alphaLik)])
        % set(gca,'XTickLabel',num2str(alpha1'))
        
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),lalpha), 'FontSize' , xyTickFontSize);
        set(gca,'XTickLabel',num2str(flipud(alphaLik)), 'FontSize' , xyTickFontSize);
        set(gca,'XDir','reverse','XGrid','on');
        
        xlabel('Level of trimmming', 'FontSize' , xyLabelSize);
        ylabel('$\hat \sigma^2_j$','Interpreter','latex', 'FontSize' , yLabelLatexSize);
        axis('manual');
        if vafter95
            axtoolbar('Visible','off');
        end
        
        % second subplot
        subplot(2,1,2);
        h2  = plot(alphaLik,Sigma2yc','LineWidth',plotLineWidth);
        % set the colors and linestyle
        set(h2,{'Color'},cellstr(clrdef(1:k)'),{'LineStyle'},linedef(1:k)',{'DisplayName'},legendGroups);
        
        xlim([min(alphaLik),max(alphaLik)]);
        % set(gca,'XTickLabel',num2str(alpha1'))
        
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),lalpha), 'FontSize' , xyTickFontSize);
        set(gca,'XTickLabel',num2str(flipud(alphaLik)), 'FontSize' , xyTickFontSize);
        set(gca,'XDir','reverse','XGrid','on');
        
        xlabel('Level of trimmming', 'FontSize' , xyLabelSize);
        ylabel('$\hat \sigma^2_{cj}$','Interpreter','latex', 'FontSize' , yLabelLatexSize);
        
        axis('manual');
        if vafter95
            axtoolbar('Visible','off');
        end
        
        if vafter95 == true
            sgtitle(tit , 'FontSize' , titleSize, 'FontWeight', 'normal');
        else
            a  = axes;
            t1 = title(tit, 'FontSize' , titleSize, 'FontWeight', 'normal');
            a.Visible = 'off'; % set(a,'Visible','off');
            t1.Visible = 'on'; % set(t1,'Visible','on');
        end
        
        if vafter91
            clickableMultiLegend(h1, 'FontSize' , legendSize);
        else
            clickableMultiLegend(h1);
        end
        
    end
    
    % this follows the first plot with legend to avoid the updates of all
    % existing graphics, because of the subsequent calls to legend and drawnow.
    drawnow limitrate nocallbacks;
    
    %% 4 Plot scatter with all regression lines (hyperplanes)
    namej = 'ScatterWithRegLines';
    alpha1range = ['[' num2str(max(alphaLik)) ' \; ' num2str(min(alphaLik)) ']' ];
    tit = {'Tclustreg monitoring plot' , ['$\quad mixt=' num2str(mixt) , '  \quad c_{\hat \sigma^2}='...
        num2str(restrfact) '\quad \alpha_{Lik}=' alpha1range ...
        '\quad \alpha_{\Sigma_X}=' num2str(alphaX) '$']};
    
    d=find(strcmp(namej,name));
    if d>0
        
        idx = IDX(:,1); % used for rotating colors in the plots
        
        if p-intercept < 2
            % The following plots are for the bi-variate case (i.e. v=1)
            
            % initialize figure
            hf5 = figure('Name',namej,'Visible','off');
            set(hf5,'Tag','tclusteda');
            
            hold on;
            
            % plot regression lines
            vv = [min(X(:,end)) max(X(:,end))];
            
            % hRegLines vector of graphic handles containing regression lines
            hRegLines = gobjects(k,1);
            % hText vector of graphic handles containing the labels of groups
            hText = gobjects(k,1);
            for jj = 1:k
                group_label = ['Group ' num2str(jj)];
                
                % jj refers to groupj 1,...k
                % jjj refers to trimming
                if intercept==1
                    
                    for jjj=1:lalpha
                        gr=plot(vv,Beta(jj,1,jjj)+Beta(jj,2,jjj)*vv,...
                            'DisplayName',[group_label ' fit' ],...
                            'Color',clrdefmap(jjj,:,jj),...
                            'LineWidth',plotLineWidthGrad,...
                            'LineStyle',linedef{jj});  %#ok<NASGU> clrdef(jj)
                    end
                    eval(['hRegLines(' num2str(jj) ')=gr;']);
                    
                elseif intercept==0
                    
                    for jjj=1:lalpha
                        gr=plot(vv,Beta(jj,1,jjj)*vv,...
                            'DisplayName',[group_label ' fit' ],...
                            'Color',clrdefmap(jjj,:,jj),...
                            'LineWidth',plotLineWidthGrad,...
                            'LineStyle',linedef{jj});  %#ok<NASGU> clrdef(jj)
                    end
                    eval(['hRegLines(' num2str(jj) ')=gr;']);
                end
                
                ucg = find(idx==jj);
                % we add a (ficticious) plot instruction with white symbols
                texth=plot(X(ucg,end),y(ucg),'.w','DisplayName',[group_label ' (' num2str(length(ucg)) ' units)']); %#ok<NASGU>
                text(X(ucg,end),y(ucg),num2str(jj*ones(length(ucg),1)),...
                    'DisplayName',[group_label ' (' num2str(length(ucg)) ' units)'] , ...
                    'HorizontalAlignment','center','VerticalAlignment','middle',...
                    'Color',clrdef(jj), 'FontSize' , 12);
                eval(['hText(' num2str(jj) ')=texth;']);
            end
            
            % Plot the outliers (trimmed points)
            ucg = find(idx==-1);
            % hunitsMinus1 = graphical handle to the first level trimmed units
            hunitsMinus1=plot(X(ucg,end),y(ucg),...
                sym1stLevelTrimmedUnits,'color',col1stLevelTrimmedUnits,...
                'MarkerSize',8,...
                'DisplayName',['Trimmed units 1st (' num2str(length(ucg)) ')']);
            
            % Plot second level trimmed units (if there are)
            ucg = find(idx==-2);
            % hunitsMinus2 = graphical handle to second level trimmed units
            hunitsMinus2=plot(X(ucg,end),y(ucg),...
                sym2ndLevelTrimmedUnits,'color',col2ndLevelTrimmedUnits,...
                'DisplayName',['Trimmed units 2nd (' num2str(length(ucg)) ')']);
            
            
            % Add clickable multilegend
            if vafter91
                clickableMultiLegend([hRegLines; hText; hunitsMinus1; hunitsMinus2],...
                    'Location','best','interpreter' , 'LaTex', 'FontSize' , legendSize) % ,'TextColor','r');
            else
                legend([hRegLines; hText; hunitsMinus1; hunitsMinus2],...
                    'Location','best')
            end
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
            plo.labeladd = ''; %DDD
            
            % group names in the legend
            group = cell(n,1);
            for iii = 1:k
                group(idx==iii) = {['Group ' num2str(iii)]};
            end
            group(idx==-1) = {'Trimmed units'};
            group(idx==-2) = {'Trimmed units level 2'};
            
            % yXplot
            % Remark: it is necessary to sort idx because in this way idx=-2 (if present) is
            % the first symbol, idx=-1 is the second ....
            [~,indsor]=sort(idx);
            [~,AX,~]=yXplot(y(indsor),X(indsor,:),group(indsor),plo);
            set(gcf, 'Tag' , 'Trimmed units');
            % Dimension of Beta is k-by-p-by-length(alpha)
            for j = 1:length(AX)
                % Make the axes of the panel with handle AX(i) the current axes.
                set(gcf,'CurrentAxes',AX(j));
                for i=1:k
                    if intercept==1
                        indexesOtherCoef=[1 setdiff((2:p),j+1)];
                    else
                        indexesOtherCoef=setdiff((1:p),j);
                    end
                    
                    % add a refline for each value of alpha1
                    for jj=1:length(alphaLik)
                        
                        idxi=IDX(:,jj)==i;
                        % Beta = k-by-p-by-length(alphaLik)
                        meanOtherX=mean(X(idxi,indexesOtherCoef)*Beta(i,indexesOtherCoef,jj)');
                        % hline  = refline(Beta(i,j+intercept,jj),meanOtherX);
                        xlimits = get(AX(j),'Xlim');
                        hline = line(xlimits , meanOtherX + Beta(i,j+intercept,jj).*xlimits);
                        
                        hline.Color     = clrdefmap(jj,:,i);
                        hline.LineStyle = linedef{i};
                    end
                end
            end
        end
        
        if vafter95 == true
            sgtitle(tit, 'FontSize' , titleSize , 'FontWeight', 'normal', 'interpreter' , 'latex');
            axtoolbar('Visible','off');
        else
            a  = axes;
            t1 = title(tit, 'FontSize' , titleSize , 'FontWeight', 'normal', 'interpreter' , 'latex');
            a.Visible = 'off'; % set(a,'Visible','off');
            t1.Visible = 'on'; % set(t1,'Visible','on');
        end
        
    end
    
    
    %% 5 Monitoring of allocation (using gscatter)
    d=find(strcmp('gscatter',name));
    tit0 = {'Tclustreg monitoring plot -- allocation of units' , ...
        'Variance Explained $\cal{V}$, for different trimming levels $\alpha$'};
    
    if d>0
        
        biggerfig = false;
        % Monitoring of allocation
        switch lalphasel
            case 1
                nr=1;
                nc=1;
            case 2
                nr=2;
                nc=1;
            case {3,4}
                nr=2;
                nc=2;
            case {5,6}
                nr=3;
                nc=2;
            case {7,8,9}
                nr=3;
                nc=3;
                biggerfig = true;
            case {10,11,12}
                nr=3;
                nc=4;
                biggerfig = true;
            otherwise
                nr=4;
                nc=4;
                biggerfig = true;
        end
        
        resup=1;
        hf6 = figure('CreateFcn',{@movegui, 'south'},'Name',['Monitoring allocation #' int2str(resup)],'Visible','off');
        set(hf6,'Tag','tclusteda');
        
        if biggerfig
            newpos = get(hf6,'Position');
            set(hf6,'Position',[left,bottom,figureResize,figureResize].*newpos);
        end
        
        if p-intercept>1
            XLmon=zeros(p+k-1,lalphasel);
            XLmonW=zeros(p+k-1,lalphasel);
        end
        
        jk=1;
        for j=1:lalphasel
            
            % The monitoring must contain a maximum of 16 panels
            % If length(alpha) is greater than 16 a new set of 16 subpanels is
            if jk>16
                jk=1;
                resup=resup+1;
                hf6b=figure('Name',['Monitoring allocation #' int2str(resup)],'Visible','off');
                set(hf6b,'Tag','tclusteda');
                subplot(nr,nc,jk);
            else
                subplot(nr,nc,jk);
            end
            
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
            
            if p-intercept>1 % More than one explanatory variable (excluding intercept)
                % In this case PLS regression is used.
                % In order to take into account group structure k-1 dummy
                % variables are added to matrix X. Of course trimmed units for
                % that particular value of alphaLik are excluded, but included
                % in the gscatter plots.
                DUM    = zeros(n,k-1);
                for jj = 1:k-1
                    DUM(idxselj==jj,jj) = 1;
                end
                Xext    = [X DUM];
                idxgt0  = idxselj>0;
                % training set
                ysel0   = y(idxgt0);
                Xsel0   = Xext(idxgt0,:);
                Xsel1   = Xext(~idxgt0,:);
                mXsel0  = mean(Xsel0);
                [XL,~,XS0,~,~,PCTVAR,~,stats] = plsregress(Xsel0,ysel0,1);
                XLmon(:,j)=XL;
                XLmonW(:,j)=stats.W;
                
                % Find best predictor for non trimmed (XS0) and trimmed units
                XS1     = zeros(n,1);
                XS1(idxgt0) = XS0;
                % REMARK
                % XS0chk=(Xsel0-mXsel0)*stats.W;
                XS1(~idxgt0) = (Xsel1-mXsel0)*stats.W;
                
                hh = gscatter(XS1,y,idxselj,clrdefj,symdefj);
                
                if jk>nc*(nr-1)
                    xlabel('PLS predictor', 'FontSize' , xyLabelSize);
                else
                    xlabel(' ');
                end
                if ismember(jk,1:nc:nc*nr)
                    ylabel('y', 'FontSize' , xyLabelSize);
                else
                    ylabel(' ');
                end
                
                if vafter91
                    clickableMultiLegend(hh, 'FontSize' , legendSize);
                else
                    clickableMultiLegend(hh);
                end
                
                if jk>1
                    legend hide
                end
                axis manual
                alphajtxt=num2str(alphaLik(alphasel(j)));
                title(['$\alpha$=' alphajtxt ' - $\cal{V}$=' num2str(100*PCTVAR(2,1),3) ],'Interpreter','latex', 'FontSize' , subtitleSize)
                
            elseif p-intercept>0  % Just one explanatory variable (excluding intercept)
                if any(isnan(idxselj))
                    hh=scatter(X(:,end),y);
                else
                    hh  =  gscatter(X(:,end),y,idxselj,clrdefj,symdefj);
                end
                % To use line function would be much faster than gscatter
                %             iiii = unique(idxselj);
                %             for idx=1:numel(iiii)
                %                 props = {'LineStyle','none','Marker',symdefj(idx),'MarkerEdge',clrdefj(idx),'MarkerSize',6};
                %                 line([X(idxselj==iiii(idx),end),X(idxselj==iiii(idx),end)],[y(idxselj==iiii(idx)),y(idxselj==iiii(idx))],...
                %                 props{:});
                %             end
                %             hhh=findobj(gcf,'Type','Line');
                
                if jk>nc*(nr-1)
                    xlabel('x1', 'FontSize' , xyLabelSize);
                end
                if ~ismember(jk,1:nc:nc*nr)
                    ylabel('', 'FontSize' , xyLabelSize);
                end
                
                if vafter91
                    clickableMultiLegend(hh, 'FontSize' , legendSize);
                else
                    clickableMultiLegend(hh);
                end
                if jk>1
                    legend hide
                end
                axis manual
                title(['$\alpha=$' num2str(alphaLik(alphasel(j)))],'Interpreter','Latex', 'FontSize' , subtitleSize)
            else
                % Univariate case: plot the histogram
                histFS(y,10,idxselj,[],[],clrdefj)
                title(['$\alpha=$' num2str(alphaLik(alphasel(j)))],'Interpreter','Latex', 'FontSize' , subtitleSize)
            end
            jk=jk+1;
        end
        
        if vafter95 == true
            sgtitle(tit0, 'FontSize' , titleSize , 'FontWeight', 'normal', 'interpreter' , 'latex');
        else
            a  = axes;
            t1 = title(tit0, 'FontSize' , titleSize , 'FontWeight', 'normal', 'interpreter' , 'latex');
            a.Visible = 'off'; % set(a,'Visible','off');
            t1.Visible = 'on'; % set(t1,'Visible','on');
        end
        
        if p-intercept>1
            % Undocumented store variable importance in presence of more than
            % one explanatory variable
            out.XLmon=XLmon;
            out.XLmonW=XLmonW;
        end
        %movegui(gcf,'south');
        
    end
    
    
    %% 6 Monitoring of beta regression coefficients (standardized)
    namej='Beta';
    tit = 'Tclustreg monitoring plot -- Estimated regression coefficients';
    
    d=find(strcmp(namej,name));
    if d>0
        hf7 = figure('CreateFcn',{@movegui, 'south'},'Name',namej,'Visible','off');
        set(hf7,'Tag','tclusteda');
        
        % Dimension of Beta is k-by-p-by-length(alpha)
        % If p is 2 first column contains interecepts and second slopes
        % Standardize Beta along the 3rd dimension
        Betast=zscore(Beta,0,3);
        biggerfig = false;
        switch p
            case 1 %p==1
                nr=1;
                nc=1;
            case 2 % p==2
                nr=2;
                nc=1;
            case {3,4} % p<=4
                nr=2;
                nc=2;
            case {5,6} %p<=6
                nr=3;
                nc=2;
                biggerfig = true;
            case {7,8,9} % p<=9
                nr=3;
                nc=3;
                biggerfig = true;
            case {10,11,12} % p<=12
                nr=3;
                nc=4;
                biggerfig = true;
            otherwise
                nr=4;
                nc=4;
                biggerfig = true;
        end
        if biggerfig
            newpos = get(hf7,'Position');
            set(hf7,'Position',[left,bottom,1.5,1.5].*newpos);
        end
        
        for j=1:p
            subplot(nr,nc,j);
            h=plot(alphaLik(:),squeeze(Betast(:,j,:))','LineWidth',plotLineWidth);
            % set the colors using the order in clrdef
            set(h,{'Color'},cellstr(clrdef(1:k)'));
            set(h,{'LineStyle'},linedef(1:k)');
            
            xlim([min(alphaLik),max(alphaLik)]);
            % set(gca,'XTickLabel',num2str(alpha1'))
            
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),lalpha)    , 'FontSize' , xyTickFontSize);
            set(gca,'XTickLabel',num2str(flipud(alphaLik)), 'FontSize' , xyTickFontSize);
            set(gca,'XDir','reverse','XGrid','on');
            
            if j>nc*(nr-1)
                xlabel('Level of trimmming', 'FontSize' , xyLabelSize);
            end
            ylabel(['$\hat \beta_' num2str(j-1) '$'],'Interpreter','latex', 'FontSize' , yLabelLatexSize);
            xtickangle(xxTickAngleVal);
            
            legend(legendGroups, 'FontSize' , legendSize);
            legend('hide');
            if j==1
                legend('show');
                if vafter91
                    clickableMultiLegend(h, 'FontSize' , legendSize);
                else
                    clickableMultiLegend(h);
                end
            end
            if vafter95
                axtoolbar('Visible','off');
            end
            axis('manual');
        end
        
        if vafter95 == true
            sgtitle(tit , 'FontSize' , titleSize , 'FontWeight', 'normal');
        else
            a  = axes;
            t1 = title(tit , 'FontSize' , titleSize , 'FontWeight', 'normal');
            a.Visible = 'off'; % set(a,'Visible','off');
            t1.Visible = 'on'; % set(t1,'Visible','on');
        end
        %movegui(gcf,'south');
        
    end
    
    
    %% 7 Monitor group size
    namej='Siz';
    tit = {'Tclustreg monitoring plot' , 'Group size'};
    
    d=find(strcmp(namej,name));
    if d>0
        % Monitoring of group size
        hf8 = figure('Name',namej,'Visible','off');
        set(hf8,'Tag','tclusteda');
        
        h=plot(alphaLik(:),out.Nopt','LineWidth',plotLineWidth);
        % set the colors using the order in clrdef
        set(h,{'Color'},cellstr(clrdef(1:k)'));
        set(h,{'LineStyle'},linedef(1:k)');
        
        xlim([min(alphaLik),max(alphaLik)]);
        % set(gca,'XTickLabel',num2str(alpha1'))
        
        lalpha=length(alphaLik);
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),lalpha), 'FontSize' , xyTickFontSize);
        set(gca,'XTickLabel',num2str(flipud(alphaLik)), 'FontSize' , xyTickFontSize);
        
        set(gca,'XDir','reverse');
        xlabel('Level of trimmming', 'FontSize' , xyLabelSize);
        legend(legendGroups);
        legend('show');
        if vafter91
            clickableMultiLegend(h, 'FontSize' , legendSize);
        else
            clickableMultiLegend(h);
        end
        title(tit,'FontSize' , titleSize , 'FontWeight', 'normal');
        
    end
    
    
    %% 8 (ex3) Monitoring posterior probabilities
    namej = 'PostProb';
    tit   = {'Tclustreg monitoring plot -- Posterior probabilities'};
    lw = 1;  % line width of all plots
    
    d=find(strcmp(namej,name));
    if d>0
        hf3 = figure('CreateFcn',{@movegui, 'south'},'Name',namej,'Visible','off');
        set(hf3,'Tag','tclusteda');
        set(hf3,'Position',[left,bottom,figureResize,figureResize].*get(hf3,'Position'));
        
        Prob1=squeeze(Postprob(:,1,:));
        Prob1(IDXmin0)=NaN;
        
        group=cell(n,1);
        group(1:n)={'Units which never changed assignment'};
        group(UnitsChgCla)={'Units which changed assignment'};
        group(UnitsTrm)={'Trimmed units'};
        
        subplot(2,2,1);
        parallelcoords(Prob1,'Group',group, 'Labels',alpha1str,'LineWidth',lw);
        ylim([-0.05 1.05]); xlim manual;
        xlabel('Level of trimming','FontSize',xyLabelSize);
        ylabel('Post prob. group 1 all units','FontSize',xyLabelSize);
        %legend('off');
        %hlpc1 = legend('hide');
        set(legend,'Location','best');
        if vafter95
            axtoolbar('Visible','off');
        end
        
        subplot(2,2,2);
        Prob1sel=Prob1(UnitsTrmOrChgCla,:);
        groupsel=group(UnitsTrmOrChgCla);
        
        parallelcoords(Prob1sel,'Group',groupsel,'Labels',alpha1str,'LineWidth',lw);
        ylim([-0.05 1.05]); xlim manual;
        xlabel('Level of trimming','FontSize',xyLabelSize);
        ylabel('Post prob. group 1 selected units','FontSize',xyLabelSize);
        %legend('off');
        %hlpc2 = legend('hide');
        set(legend,'Location','best');
        
        % Add the label of the units whose final post prob is intermediate
        unitswithText=Prob1sel(:,end)>0.05 &  Prob1sel(:,end)<0.95;
        text(lalpha*ones(sum(unitswithText),1),Prob1sel(unitswithText,end),...
            cellstr(num2str(UnitsTrmOrChgCla(unitswithText))));
        if vafter95
            axtoolbar('Visible','off');
        end
        
        subplot(2,2,3)
        % Prob1table=array2table(Prob1,'VariableNames',cellstr(num2str(alphaLik)));
        if vafter95
            parallelplot(Prob1,'GroupData',group,'FontSize',12,'LineWidth',plotLineWidth);
            set(gca,'CoordinateTickLabels',cellstr(num2str(alphaLik)))
            xlabel('Level of trimming');
            ylabel('Post prob. group 1 all units');
            legend('off');
        end
        
        subplot(2,2,4)
        if vafter95
            parallelplot(Prob1(UnitsTrmOrChgCla,:),'GroupData',group(UnitsTrmOrChgCla),...
                'FontSize',12,'LineWidth',plotLineWidth,'LineAlpha',0.99);
            set(gca,'CoordinateTickLabels',cellstr(num2str(alphaLik)));
            xlabel('Level of trimming');
            ylabel('Post prob. group 1 selected units');
            legend('off');
        end
        
        if vafter95 == true
            sgtitle(tit , 'FontSize' , titleSize, 'FontWeight', 'normal');
        else
            a  = axes;
            t1 = title(tit, 'FontSize' , titleSize, 'FontWeight', 'normal');
            a.Visible = 'off'; % set(a,'Visible','off');
            t1.Visible = 'on'; % set(t1,'Visible','on');
        end
        %movegui(gcf,'south');
    end
    
    % make all figures visible again
    set(findobj('Tag','tclusteda'),'Visible','on');
    
    % make heavy legends visible again
    %set([hlpc1,hlpc2],'Location','best','Visible','on');
    
end

%% inner function
    function cm = cmapFS(cstart,cend,m)
        %cmapFS creates m RGB colors, with a gradient from cstart to cend
        % We use intensity colors in [0 1], not RGB values in [1 255]. So,
        % output cm is a matrix of m*3 elements containing the color
        % gradient.
        %
        % Example:
        %    c = cmapFS([1 0 0],[0.5 0.8 1],12);
        %    surf(peaks)
        %    colormap(c);
        
        % default number of colors in the gradient is 32.
        if nargin < 3
            m=32;
        end
        
        % conversion in case numbers are given in the [1 255] standard
        if max(cstart)>1 && max(cstart) <= 255
            cstart=cstart./255;
        end
        if max(cend)>1 && max(cend) <= 255
            cend=cend./255;
        end
        
        mm = (cend-cstart)/(m-1);
        iR  = (cstart(1):mm(1):cend(1))';
        iG  = (cstart(2):mm(2):cend(2))';
        iB  = (cstart(3):mm(3):cend(3))';
        cm  = [iR ,iG ,iB];
    end

end
%FScategory:CLUS-RobClaREG
