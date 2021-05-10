function [out, varargout] = LTSts(y,varargin)
%LTSts extends LTS estimator to time series
%
%<a href="matlab: docsearchFS('LTSts')">Link to the help function</a>
%
% It is possible to set a model with a trend (up to third order), a
% seasonality (constant or of varying amplitude and with a different number
% of harmonics) and a level shift (in this last case it is possible to
% specify the window in which level shift has to be searched for).
%
%  Required input arguments:
%
%    y:         Time series to analyze. Vector. A row or a column vector
%               with T elements, which contains the time series.
%
%
%  Optional input arguments:
%
%
%         bdp : breakdown point. Scalar. It measures the fraction of outliers
%               the algorithm should resist. In this case any value greater
%               than 0 but smaller or equal than 0.5 will do fine. Please
%               specify h or bdp, but not both.
%                 Example - 'bdp',0.4
%                 Data Types - double
%
%     conflev : Confidence level. Scalar. Scalar between 0 and 1 containing
%               Confidence level which is used to declare units as
%               outliers. Usually conflev=0.95, 0.975 0.99 (individual
%               alpha) or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous
%               alpha). Default value is 0.975.
%                 Example - 'conflev',0.99
%                 Data Types - double
%
%  dispresults : Display results of final fit. Boolean. If dispresults is
%               true,  labels of coefficients, estimated coefficients,
%               standard errors, tstat and p-values are shown on the
%               screen in a fully formatted way. The default value of
%               dispresults is false.
%               Example - 'dispresults',true
%               Data Types - logical
%
%
%           h : The number of observations that determined the least
%               trimmed squares estimator. Scalar. h is an integer greater
%               than p (number of columns of matrix X including the
%               intercept but smaller then n. If the purpose is outlier
%               detection than h does not have to be smaller than
%               [0.5*(T+p+1)]. The default value of h is [0.75*T]. Note
%               that if h is supplied input argument bdp is ignored.
%                 Example - 'h',round(n*0.75)
%                 Data Types - double
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
% lshiftlocref: Parameters for local shift refinement. Structure.
%               This option is used just if model.lshift is greater then 0.
%               In order to precisely identify level shift position it is
%               necessary to consider a local sum of squares varying the
%               position of the level shift around the first tentative
%               position keeping all the other parameters fixed. This
%               structure contains the following fields:
%               lshiftlocref.wlength = scalar greater than 0 which
%                   identifies the length of the window. The default value
%                   is 15, that is the tentative level shift position
%                   varies from tl-15, tl-15, ..., tl+14, tl+15, where tl is
%                   the best preliminary tentative level shift position.
%              lshiftlocref.typeres = scalar which identifies the type of
%                   residuals to consider. If typerres =1, the local
%                   residuals sum of squares is based on huberized (scaled)
%                   residuals (this is the default
%                   choice) else raw residuals are used.
%              lshiftlocref.huberc= tuning constant for Huber estimator just
%                   in case lshiftlocref.typeres=1. The default value is 2.
%               Example - 'lshiftlocref',lshiftlocref.typeres=2
%               Data Types - struct
%
%       lts   : structure which controls a set of options of the
%               maximization procedure. Structure. Structure with the
%               following fields:
%                   lts.bestr   = scalar defining number of "best betas" to
%                               remember from the subsamples. These will be
%                               later iterated until convergence.
%                               The default is 20 (10 of them are the best
%                               from previous iteration in case a level
%                               shift is present).
%                  lts.refsteps = scalar defining number of concentration
%                               steps (default = 2). refsteps = 0 means
%                               "raw-subsampling" without iterations.
%             lts.refstepsbestr = scalar defining maximum number of refining
%                               steps for each best subset (default=50).
%                   lts.reftol  = scalar. Default value of tolerance for
%                               the refining steps
%                               The default value is 1e-6;
%              lts.reftolbestr  = scalar. Default value of tolerance for
%                               the refining steps for each of the best
%                               subsets The default value is 1e-8.
%                Example - 'lts',lts
%                Data Types - struct
%              Remark: if lts is an empty value all default values of
%              structure lts will be used.
%
%      model :  model type. Structure. A structure which specifies the model
%               which will be used. The model structure contains the following
%               fields:
%               model.s = scalar (length of seasonal period). For monthly
%                         data s=12 (default), for quartely data s=4, ...
%               model.trend = scalar (order of the trend component).
%                       trend = 0 implies no trend;
%                       trend = 1 implies linear trend with intercept (default);
%                       trend = 2 implies quadratic trend;
%                       trend = 3 implies cubic trend.
%                       Admissible values for trend are, 0, 1, 2 and 3.
%                       In the paper RPRH to denote the order of the trend
%                       symbol A is used.
%               model.seasonal = scalar (integer specifying number of
%                        frequencies, i.e. harmonics, in the seasonal
%                        component. Possible values for seasonal are
%                        $0,1, 2, ..., [s/2]$, where $[s/2]=floor(s/2)$.
%                        If seasonal =0 we assume there is no seasonal
%                        component.
%                        If seasonal =1 (default) we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 sin ( 2 \pi t/s)$;
%                        if seasonal =2 we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 \sin ( 2 \pi t/s)
%                        + \beta_3 \cos(4 \pi t/s) + \beta_4 \sin (4 \pi t/s)$.
%                        Note that when $s$ is even the sine term disappears
%                        for $j=s/2$ and so the maximum number of
%                        trigonometric parameters is $s-1$.
%                        If seasonal is a number greater than 100 then it
%                        is possible to specify how the seasonal component
%                        grows over time.
%                        For example, seasonal = 101 implies a seasonal
%                        component which just uses one frequency
%                        which grows linearly over time as follows:
%                        $(1+\beta_3 t)\times ( \beta_1 cos( 2 \pi t/s) +
%                        \beta_2 \sin ( 2 \pi t/s))$.
%                        For example, seasonal =201 implies a seasonal
%                        component which just uses one frequency
%                        which grows in a quadratic way over time as
%                        follows:
%                        $(1+\beta_3 t + \beta_4  t^2)\times( \beta_1 \cos(
%                        2 \pi t/s) + \beta_2 \sin ( 2 \pi t/s))$.
%                        seasonal =0 implies a non seasonal model.
%                       In the paper RPRH to denote the number of
%                       frequencies of the seasonal component
%                       symbol B is used, while symbol G is used to denote
%                       the order of the trend of the seasonal component.
%                       Therefore, for example, model.seasonal=201
%                       corresponds to B=1 and G=2, while model.seasonal=3
%                       corresponds to B=3 and G=0;
%               model.X  =  matrix of size T-by-nexpl containing the
%                         values of nexpl extra covariates which are likely
%                         to affect y.
%               model.lshift = scalar or vector associated to level shift
%                       component. lshift=0 (default) implies no level
%                       shift component.
%                       If model.lshift is vector of positive integers,
%                         then it is associated to the positions of level
%                         shifts which have to be considered. The most
%                         significant one is included in the fitted model.
%                         For example if model.lshift =[13 20 36] a
%                         tentative level shift is imposed in position
%                         $t=13$, $t=20$ and $t=36$. The most significant
%                         among these positions in included in the final
%                         model. In other words, the following extra
%                         parameters are added to the final model:
%                         $\beta_{LS1}* I(t \geq \beta_{LS2})$ where
%                         $\beta_{LS1}$ is a real number (associated with
%                         the magnitude of the level shift) and
%                         $\beta_{LS2}$ is an integer which assumes values
%                         13, 20 or 36 and and $I$ denotes the indicator
%                         function.
%                         As a particular case, if model.lshift =13 then a
%                         level shift in position $t=13$ is added to the
%                         model. In other words the following additional
%                         parameters are added: $\beta_{LS1}* I(t \geq 13)$
%                         where $\beta_{LS1}$ is a real number and $I$
%                         denotes the indicator function.
%                       If lshift = -1 tentative level shifts are
%                         considered for positions $p+1,p+2, ..., T-p$ and
%                         the most significant one is included in the final
%                         model ($p$ is the total number of parameters in
%                         the fitted model). Note that lshift=-1 is not
%                         supported for C-coder translation.
%                       In the paper RPRH $\beta_{LS1}$ is denoted with
%                       symbol $\delta_1$, while, $\beta_{LS2}$ is denoted
%                       with symbol $\delta_2$.
%               model.ARp = scalar greater or equal than 0 which
%                         specifies the length of the autoregressive
%                         component. The default value of model.ARp is 0,
%                         that is there is no autoregressive component.
%                 Example - 'model', model
%                 Data Types - struct
%               Remark: the default model is for monthly data with a linear
%               trend (2 parameters) + seasonal component with just one
%               harmonic (2 parameters), no additional explanatory
%               variables and no level shift that is
%                               model=struct;
%                               model.s=12;
%                               model.trend=1;
%                               model.seasonal=1;
%                               model.X='';
%                               model.lshift=0;
%               Using the notation of the paper RPRH we have A=1, B=1; and
%               $\delta_1=0$.
%
%        msg  : Messages on the screen. Boolean.
%               Scalar which controls whether to display or not messages
%               on the screen. If msg==true (default) messages are displayed on
%               the screen about estimated time to compute the estimator
%               and the warnings about 'MATLAB:rankDeficientMatrix',
%               'MATLAB:singularMatrix' and 'MATLAB:nearlySingularMatrix'
%               are set to off else no message is displayed on the screen
%               Example - 'msg',true
%               Data Types - logical
%
%nbestindexes : position of the best solutions. Positive integer. For each
%               tentative level shift solution, it is interesenting to
%               understand whether best solutions of target function come
%               from subsets associated with current level shift solution
%               or from best solutions from previous tentative level shift
%               position.  The indexes from 1 to lts.bestr/2 are associated
%               with subsets just extracted. The indexes from lts.bestr/2+1
%               to lts.bestr are associated with best solutions from
%               previous tentative level shift. More precisely:
%               index lts.bestr/2+1 is associated with best solution from
%               previous tentative level shift;
%               index lts.bestr/2+2 is associated with second best solution
%               from previous tentative level shift;
%               ...
%               nbestindexes is an integer which specifies how many indexes
%               we want to store. The default value of nbestindexes  is 3.
%               Example - 'nbestindexes',5
%               Data Types - double
%
%      nocheck: Check input arguments. Boolean. If nocheck is equal to true no
%               check is performed on matrix y and matrix X. Notice that y
%               and X are left unchanged. In other words the additioanl
%               column of ones for the intercept is not added. As default
%               nocheck=false. The controls on h, bdp and nsamp still remain.
%               Example - 'nocheck',true
%               Data Types - boolean
%
%       nsamp : number of subsamples to extract. Scalar or vector of length 2.
%               Vector of length 1 or 2 which controls the number of
%               subsamples which will be extracted to find the robust
%               estimator. If lshift is not equal to 0 then nsamp(1)
%               controls the number of subsets which have to be extracted
%               to find the solution for t=lshift(1). nsamp(2) controls the
%               number of subsets which have to be extracted to find the
%               solution for t=lshift(2), lshift(3), ..., lshift(end).
%               Note that nsamp(2) is generally smaller than nsamp(1)
%               because in order to compute the best solution for
%               t=lshift(2), lshift(3), ..., lshift(end), we use the lts.bestr/2
%               best solutions from previous t (after shifting the
%               position of the level shift in the estimator of beta). If
%               lshift is a vector of positive integers the default value
%               of nsamp is (500 250). If
%               lshift is a vector of positive integers and nsamp is supplied as a scalar the default
%               is to extract [nsamp/2] subsamples for t=lshift(1),
%               lshift(2), ... Therefore, for example, in order to extract
%               600 subsamples for t=lshift(1) and 300 subsamples for t=
%               lshift(2) ... you can use nsamp =600 or nsamp=[600 300].
%               The default value of nsamp is 1000;
%                 Example - 'nsamp',500
%                 Data Types - double
%               Remark: if nsamp=0 all subsets will be extracted.
%               They will be (n choose p).
%
% refstepsALS :   Maximum iterations inside ALS. Scalar. Maximum number
%                 of iterations inside ALS routine. Default value of
%                 tolerance for the refining steps inside ALS routine. The
%                 default value is 50.
%                 Example - 'refstepsALS',20
%                 Data Types - double
%
%
%  reftolALS  :   Tolerance inside ALS. Scalar. Tolerance value of tolerance
%                 for the refining steps inside ALS routine. The default
%                 value is 1e-03.
%                 Example - 'reftolALS',1e-05
%                 Data Types - double
%
%SmallSampleCor:Small sample correction factor to control empirical size of
%               the test.  Scalar equal to 1 or 2 (default) or 3 or 4.
%               - If SmallSampleCor=1 in the reweighting step the nominal
%                 threshold based on $\chi^2_{0.99}$ is multiplied by the
%                 small sample correction factor which guarrantees that the
%                 empirical size of the test is equal to the nominal size.
%                 Given that the correction factors were obtained through
%                 simulation for a linear model, the number of explanatory
%                 which is used to compute the correction factor refers to
%                 all explanatory variables except the non linear components
%                 in the seasonal part of the model. For example, in a model
%                 with linear trend 4 seasonal harmonics + level shift and
%                 second order trend in the seasonal component the number of
%                 explanatory variables used is 11 = total number of
%                 variables -2 = 2 (linear trend) + 8 (4 seasonal harmonics)
%                 +1 (level shift).
%               - If SmallSampleCor =2 Gervini and Yohai procedure is called
%                 with 'iterating' false and 'alpha' 0.99 is invoked, that is:
%                 weights=GYfilt(stdres,'iterating',false,'alpha',0.99);
%               - If SmallSampleCor =3 Gervini and Yohai procedure  called
%                 with 'iterating' true and 'alpha' 0.99 is invoked, that is:
%                 weights=GYfilt(stdres,'iterating',true,'alpha',0.99);
%               - If SmallSampleCor =4  $\chi^2_{0.99}$ threshold is used that is:
%                 weights = abs(stdres)<=sqrt(chi2inv(0.99,1));
%                 Example - 'SmallSampleCor',3
%                 Data Types - double
%
%
%       yxsave : store X and y. Boolean. Scalar that is set to 1 to request
%                that the response vector y and data matrix X are saved
%                into the output structure out.
%                Default is 0, i.e. no saving is done.
%               Example - 'yxsave',1
%               Data Types - logical
%
%       plots : Plots on the screen. Scalar.
%               If plots = 1, a two panel plot will be shown on the screen.
%               The upper panel contains the orginal time series with
%               fitted values. The bottom panel will contain the plot
%               of robust residuals against index number. The confidence
%               level which is used to draw the horizontal lines associated
%               with the bands for the residuals is specified in input
%               option conflev. If conflev is missing a nominal 0.975
%               confidence interval will be used. If plots =2 the following
%               additional plots will be shown on the screen.
%               1) Boxplot of the distribution of the lts.bestr values of
%               the target function for each tentative level shift position;
%               2) A two panel plot which shows the values of the local sum
%               of squares varying the position of the level shift around
%               the first tentative position keeping all the other
%               parameters fixed. Top panel refers to Huberized residuals
%               sum of squares and bottom panel refers to residual sum of
%               squares.
%               3) A plot which shows the indexes of the best nbestindexes
%               solutions for each tentative level shift position.
%               4) A plot which shows the relative frequency of inclusion
%               of each unit in the best h-subset after lts.refsteps
%               refining steps.
%               5) A plot which shows the relative frequency of inclusion
%               of each unit inside the best nbestindexes subsets which are
%               brought to full convergence.
%               The default value of plot is 0 i.e. no plot is shown on the
%               screen.
%                 Example - 'plots',1
%                 Data Types - double
%
%       Remark: The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%
%  Output:
%
%  out :     A structure containing the following fields
%
%             out.B =   Matrix containing estimated beta coefficients,
%                       (including the intercept when options.intercept=true)
%                       standard errors, t-stat and p-values
%                       The content of matrix B is as follows:
%                       1st col = beta coefficients
%                        The order of the beta coefficients is as follows:
%                        1) trend elements (if present). If the trend is
%                        of order two there are r+1 coefficients if the
%                        intercept is present otherwise there are just r
%                        components;
%                        2) linear part of seasonal component 2, 4, 6, ...,
%                        s-2, s-1 coefficients (if present);
%                        3) coefficients associated with the matrix of
%                        explanatory variables which have a potential effect
%                        on the time series under study (X or
%                        autoregressive part); If model.ARp>0 the first
%                        model.ARp elements refer to the autoregressive
%                        component.
%                        4) non linear part of seasonal component, that is
%                        varying amplitude. If varying amplitude is of order
%                        k there are k coefficients (if present);
%                        5) level shift component (if present). In out.B it
%                        is shown just the real number which identifies the
%                        magnitude of the upward (downward) level shift.
%                        The integer which specifies the time in which
%                        level shift takes place is given in output
%                        out.posLS.
%                       2nd col = standard errors;
%                       3rd col = t-statistics;
%                       4th col = p values.
%          out.Btable = same thing out.B but n table format.
%               out.h = The number of observations that have determined the
%                       initial LTS estimator, i.e. the value of h.
%              out.bs = Vector containing the units with the smallest p+k
%                       squared residuals before the reweighting step,
%                       where p is the total number of the parameters in
%                       the model and p+k is smallest number of units such
%                       that the design matrix is full rank.
%                       out.bs can be used to initialize the forward
%                       search.
%         out.Hsubset = matrix of size T-by-r
%                       containing units forming best H subset for each
%                       tentative level shift which is considered. r is
%                       number of tentative level shift positions whicha re
%                       considered. For example if model.lshift=[13 21 40]
%                       r is equal to 3. Units belonging to subset are
%                       given with their row number, units not belonging to
%                       subset have missing values
%                       This output is present just if input option
%                       model.lshift is not equal to 0.
%           out.posLS = scalar associated with best tentative level shift
%                       position. This output is present just if input
%                       option model.lshift is not equal to 0.
%     out.numscale2 = matrix of size lts.bestr-by-(T-2*lshift) containing
%                       (in the columns) the values of the lts.bestr smallest
%                       values of the target function. Target function = truncated
%                       residuals sum of squares.
%     out.BestIndexes = matrix of size nbestindexes-by-(T-2*lshift)
%                       containing in each column the indexes
%                       associated with the best nbestindexes solutions.
%                       The indexes from lts.bestr/2+1 to lts.bestr are
%                       associated with best solutions from previous
%                       tentative level shift.
%                       More precisely:
%                       index lts.bestr/2+1 is associated with best solution
%                       from previous tentative level shift;
%                       index lts.bestr/2+2 is associated with best solution
%                       from previous tentative level shift.
%                       This output is present just if input option
%                       model.lshift is not equal to 0.
%         out.Likloc  = matrix of size (2*lshiftlocref.wlength+1)-by-3
%                       containing local sum of squares of residuals in
%                       order to decide best position of level shift:
%                       1st col = position of level shift;
%                       2nd col = local sum of squares of huberized residuals;
%                       3rd col = local sum of squares of raw residuals.
%                       This output is present just if input option
%                       model.lshift is not equal to 0.
%             out.RES = Matrix of size T-by-(T-lshift) containing scaled
%                       residuals for all the T units of the original time
%                       series monitored in steps lshift+1, lshift+2, ...,
%                       T-lshift, where lshift+1 is the first tentative
%                       level shift position, lshift +2 is the second level
%                       shift position, and so on. This output is present
%                       just if input option model.lshift is not equal to 0.
%            out.yhat = vector of fitted values after final (NLS=non linear
%                       least squares) step.
%                       $ (\hat \eta_1, \hat \eta_2, \ldots, \hat \eta_T)'$
%       out.residuals = Vector T-by-1 containing the scaled residuals from
%                       after final NLS step.
%         out.weights = Vector containing weights after adaptive
%                       reweighting. The elements of
%                       this vector are 0 or 1. These weights identify the
%                       observations which are used to compute the final
%                       NLS estimate.
%           out.scale = Final scale estimate of the residuals using final weights.
%                     \[
%                     \hat \sigma = cor \times \sum_{i \in S_m} [y_i- \eta(x_i,\hat \beta)]^2/(m-p)
%                     \]
%                     where $S_m$ is a set of cardinality $m$ which
%                     contains the units not declared as outliers, $p$
%                     is the total number of estimated parameters and $cor$
%                     is a correction factor to make the estimator
%                     consistent.
%         out.conflev = confidence level which is used to declare outliers.
%                       Remark: scalar out.conflev will be used to draw the
%                       horizontal lines (confidence bands) in the plots
%        out.outliers = vector containing the list of the units declared
%                       as outliers using confidence level specified in
%                       input scalar conflev.
%   out.outliersPval  =  p-value of the units declared as outliers.
%         out.singsub = Number of subsets wihtout full rank. Notice that if
%                       this number is greater than 0.1*(number of
%                       subsamples) a warning is produced on the screen
%            out.invXX = $cov(\beta)/\hat \sigma^2$. p-by-p, square matrix.
%                       If the model is linear out.invXX  is equal to
%                       $(X'X)^{-1}$, else out.invXX is equal to $(A'A)^{-1}$
%                       where $A$ is the matrix of partial derivatives. More
%                       precisely:
%                       \[
%                       a_{i,j}=\frac{\partial \eta_i(x_i, \hat \beta)}{\partial \hat \beta_j}
%                       \]
%                       where
%                       \begin{eqnarray}
%                       y_i & = & \eta(x_i,\beta)+ \epsilon_i  \\
%                           & = & \eta_i +\epsilon_i \\
%                           & = & \eta(x_i,\hat \beta)+ e_i  \\
%                           & = & \hat \eta_i + e_i
%                       \end{eqnarray}
% out.LastHarmonicPval = combined p value for the two coefficients of the
%                        last harmonic (this p value comes from an F test).
% out.LevelShiftPval  =  p value of the level shift which takes into
%                       account (this pvalue comes from Bonferronization to
%                       take it account that if you impose a level shift,
%                       this component is always found).
%            out.y    = response vector y.
%            out.X    = data matrix X containing trend, seasonal, expl
%                       (with autoregressive component) and
%                       lshift, if the model is linear or linearized
%                       version of $\eta(x_i, \beta)$ if the model is non
%                       linear containing in the columns partial
%                       derivatives evaluated in correspondence of
%                       out.B(:,1) with respect to each parameter. In other
%                       words, the $i,j$-th element of out.X is
%                       \[
%                       \frac{\partial \eta_i(x_i, \hat \beta)}{\partial \hat \beta_j}
%                       \]
%                       $j=1, 2, \ldots, p$, $i \in S_m$.
%                       The size of this matrix is:
%                       n-length(out.outliers)-by-p
%                       The field is present only if option
%                       yxsave is set to 1.
%           out.class = 'LTSts'.
%
%  Optional Output:
%
%            C        : cell  containing the indices of the subsamples
%                       extracted for computing the estimate (the so called
%                       elemental sets) for each tentative level shift
%                       position.
%                       C{1} is associated with the subsamples for
%                       first tentative level shift position;
%                       C{2} is associated with the subsamples for
%                       second tentative level shift position;
%                       ...
%                       C{end} is associated with the subsamples for
%                       last tentative level shift position;
%
% See also LXS, wedgeplot
%
% References:
%
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [RPRH]
%
%
% Copyright 2008-2021.
% Written by Marco Riani, Domenico Perrotta, Peter
% Rousseeuw and Mia Hubert
%
%
%<a href="matlab: docsearchFS('LTSts')">Link to the help function</a>
%
%$LastChangedDate:: 2019-12-15 21:09:21 #$: Date of the last commit

% Examples:


%{
    % Simulated data with linear trend and level shift.
    % No seasonal component.
    rng('default')
    T=45;
    a=1;
    b=0.8;
    sig=1;
    seq=(1:T)';
    y=a+b*seq+sig*randn(T,1);
    % Add a level shift in the simulated series
    y(round(T/2):end)=y(round(T/2):end)+10;
    % model with a linear trend, non seasonal and level shift
    model=struct;
    model.trend=1;
    model.seasonal=0;
    % Potential level shift position is investigated in positions:
    % t=11, t=12, ..., t=T-10.
    model.lshift=11:T-10;
    out=LTSts(y,'model',model,'plots',1);
    % Using the notation of the paper RPRH: A=1, B=1, G=0 and $\delta_1>0$.
    str=strcat('A=1, B=0, G=0, $\delta_2=',num2str(out.posLS),'$');
    title(findobj(gcf,'-regexp','Tag','LTSts:ts'),str,'interpreter','latex');
%}

%{
    % Airline data: linear trend + just one harmonic for seasonal component.
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl

    y=(y(:));
    yr = repmat((1949:1960),12,1);
    mo = repmat((1:12)',1,12);
    time = datestr(datenum(yr(:),mo(:),1));
    ts = timeseries(y(:),time,'name','AirlinePassengers');
    ts.TimeInfo.Format = 'dd-mmm-yyyy';
    tscol = tscollection(ts);
    % plot airline data
    plot(ts)
    % linear trend + just one harmonic for seasonal component
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=1;           % just one harmonic
    model.lshift=0;             % no level shift
    out=LTSts(y,'model',model,'dispresults',true);

    close all
    % Plot real and fitted values
    plot(y,'Linewidth',1.5);
    hold('on')
    plot(out.yhat,'r--','Linewidth',1.5)
    legend({'Real values','Fitted values'},'Location','SouthEast','interpreter','LaTeX','FontSize',14)
    numpar = {'model parameters:' , 'A=1, B=1, G=0, $\delta_1=0$'};
    title(gca,numpar,'interpreter','LaTeX','FontSize',16);
%}

%{
    % Model with linear trend and six harmonics for seasonal component.
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl

    y=(y(:));
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=6;           % six harmonics
    model.lshift=0;             % no level shift
    out=LTSts(y,'model',model);

    close all
    % Plot real and fitted values
    plot(y,'Linewidth',1.5);
    hold('on')
    plot(out.yhat,'r--','Linewidth',1.5)
    legend({'Real values','Fitted values'},'Location','SouthEast','interpreter','LaTeX','FontSize',14)
    numpar = {'model parameters:' , 'A=1, B=6, G=0, $\delta_1=0$'};
    title(gca,numpar,'interpreter','LaTeX','FontSize',16);

%}

%{
    % Model with linear trend, two harmonics for seasonal component and
    % varying amplitude using a linear trend (1).
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl

    y=(y(:));
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=102;         % two harmonics with time varying seasonality
    model.lshift=0;             % no level shift
    out=LTSts(y,'model',model);

    close all
    % Plot real and fitted values
    plot(y,'Linewidth',1.5);
    hold('on')
    plot(out.yhat,'r--','Linewidth',1.5)
    legend({'Real values','Fitted values'},'Location','SouthEast','interpreter','LaTeX','FontSize',14)
    numpar = {'model parameters:' , 'A=1, B=2, G=1, $\delta_1=0$'};
   title(gca,numpar,'interpreter','LaTeX','FontSize',16);
%}

%{
    % Model with linear trend, six harmonics for seasonal component and
    % varying amplitude using a linear trend (2).
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl

    y=(y(:));
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=106;         % six harmonics with linear time varying seasonality
    model.lshift=0;             % no level shift
    % out=fitTSLS(y,'model',model);
    out=LTSts(y,'model',model);

    close all
    % Plot real and fitted values
    plot(y,'Linewidth',1.5);
    hold('on')
    plot(out.yhat,'r--','Linewidth',1.5)
    legend({'Real values','Fitted values'},'Location','SouthEast','interpreter','LaTeX','FontSize',14)
    numpar = {'model parameters:' , 'A=1, B=6, G=1, $\delta_1=0$'};
   title(gca,numpar,'interpreter','LaTeX','FontSize',16);

%}

%{
    % Contaminated time series with upward level shift.
    % Model with linear trend, six harmonics for seasonal component and
    % varying amplitude using a linear trend).
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl

    y=(y(:));
    yLS=y;
    yLS(55:end)=yLS(55:end)+130;
    model=struct;
    model.trend=1;                          % linear trend
    model.s=12;                             % monthly time series
    model.seasonal=1;
    model.lshift=14:length(yLS)-13;         % impose level shift
    out=LTSts(yLS,'model',model);

    close all
    % Plot real and fitted values
    plot(yLS,'Linewidth',1.5);
    hold('on')
    plot(out.yhat,'r--','Linewidth',1.5)
    legend({'Real values','Fitted values'},'Location','SouthEast','interpreter','LaTeX','FontSize',14)
    % Using the notation of the paper RPRH: A=1, B=1, G=0 and $\delta_1>0$.
    str=strcat('A=1, B=1, G=0, $\delta_2=',num2str(out.posLS),'$');
    numpar = {'model parameters:' , str};
    title(gca,numpar,'interpreter','LaTeX','FontSize',16);

%}

%{
    % Contaminated time series with downward level shift.
    % Model with linear trend, three harmonics for seasonal component and
    % varying amplitude using a linear trend).
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl

    y=(y(:));
    yLS=y;
    yLS(35:end)=yLS(35:end)-300;
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=103;
    model.lshift=26:length(yLS)-25;
    out=LTSts(yLS,'model',model,'msg',0);

    close all
    % Plot real and fitted values
    plot(yLS,'Linewidth',1.5);
    hold('on')
    plot(out.yhat,'r--','Linewidth',1.5)
    legend({'Real values','Fitted values'},'Location','SouthEast','interpreter','LaTeX','FontSize',14)
    % Using the notation of the paper RPRH: A=1, B=3, G=1 and $\delta_1>0$.
    str=strcat('A=1, B=3, G=1, $\delta_2=',num2str(out.posLS),'$');
    numpar = {'model parameters:' , str};
    title(gca,numpar,'interpreter','LaTeX','FontSize',16);
%}

%{
    % Model with an explanatory variable using log-transformed series.
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl

    y=(y(:));
    y1=log(y);
    % Model with linear trend, two harmonics for seasonal component and
    % varying amplitude using a linear trend).
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=106;
    model.lshift=0;
    model.X=randn(length(y),1);
    out=LTSts(y1,'model',model);

    close all
    % Plot real and fitted values
    plot(y1,'Linewidth',1.5);
    hold('on')
    plot(out.yhat,'r--','Linewidth',1.5)
    legend({'Real values','Fitted values'},'Location','SouthEast','interpreter','LaTeX','FontSize',14)
    % Using the notation of the paper RPRH: A=1, B=6, G=1 and $\delta_1>0$.
    str=strcat('A=1, B=6, G=1, $\delta_1=0$');
    numpar = {'model parameters:' , str};
    title(gca,numpar,'interpreter','LaTeX','FontSize',16);
%}

%{
    %% Example 1 used in the paper RPRH.
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Two short level shifts in opposite directions and an isolated outlier.
    % Add a level shift contamination plus some outliers.
    y1=y(:);
    y1(50:55)=y1(50:55)-300;
    y1(70:75)=y1(70:75)+300;
    y1(90:90)=y1(90:90)+300;
    % Create structure specifying model
    model=struct;
    model.trend=2;                  % quadratic trend
    model.s=12;                     % monthly time series
    model.seasonal=204;             % number of harmonics
    model.lshift=41:length(y1)-40;  % position where monitoring level shift
    model.X='';
    % Create structure lts specifying lts options
    lshiftlocref=struct;
    % Set window length for local refinement.
    lshiftlocref.wlength=10;
    % Set tuning constant to use insde Huber rho function
    lshiftlocref.huberc=1.5;
    % Estimate the parameters
    [out]=LTSts(y1,'model',model,'nsamp',500,...
       'plots',1,'lshiftlocref',lshiftlocref,'msg',0);
    % Using the notation of the paper RPRH: A=2, B=4, G=2 and $\delta_1>0$.
    str=strcat('A=2, B=4, G=2, $\delta_2=',num2str(out.posLS),'$');
    numpar = {'model parameters:' , str};
    axeslast=findobj('-regexp','Tag','LTSts:ts');
    title(axeslast(end),numpar,'interpreter','LaTeX','FontSize',16);

    % generate the wedgeplot
    % wedgeplot(out,'transpose',true,'extradata',[y1 out.yhat]);
%}

%{
    %% Example 2 used in the paper RPRH.
    % A persisting level shift and three isolated outliers, two of which in
    % proximity of the level shift.
        % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    close all
    y1=y(:);
    y1(68:end)=y1(68:end)+1300;
    y1(67)=y1(67)-600;
    y1(45)=y1(45)-800;
    y1(68:69)=y1(68:69)+800;
    % Create structure specifying model
    model=struct;
    model.trend=2;                  % quadratic trend
    model.s=12;                     % monthly time series
    model.seasonal=204;             % number of harmonics
    model.lshift=41:length(y1)-40;  % position where monitoring level shift
    model.X='';
    % Create structure lts specifying lts options
    lshiftlocref=struct;
    % Set window length for local refinement.
    lshiftlocref.wlength=10;
    % Set tuning constant to use insde Huber rho function
    lshiftlocref.huberc=1.5;
    % Estimate the parameters
    [out, varargout]=LTSts(y1,'model',model,'nsamp',500,...
       'plots',1,'lshiftlocref',lshiftlocref,'msg',0);

    % Using the notation of the paper RPRH: A=2, B=4, G=2 and $\delta_1>0$.
    str=strcat('A=2, B=4, G=2, $\delta_2=',num2str(out.posLS),'$');
    numpar = {'model parameters:' , str};
    title(findobj('-regexp','Tag','LTSts:ts'),numpar,'interpreter','LaTeX','FontSize',16);

    % generate the wedgeplot
    % wedgeplot(out,'transpose',true,'extradata',[y1 out.yhat]);

%}

%{
    %% Example 3 used in the paper RPRH.
    % A persisting level shift preceded and followed in the proximity by
    % other two short level shifts, and an isolated outlier.
        % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    y1=y(:);
    y1(50:55)=y1(50:55)-300;
    y1(68:end)=y1(68:end)-700;
    y1(70:75)=y1(70:75)+300;
    y1(90:90)=y1(90:90)+300;
    % Create structure specifying model
    model=struct;
    model.trend=2;                  % quadratic trend
    model.s=12;                     % monthly time series
    model.seasonal=204;             % number of harmonics
    model.lshift=41:length(y1)-40;  % position where monitoring level shift
    model.X='';
    % Create structure lts specifying lts options
    lshiftlocref=struct;
    % Set window length for local refinement.
    lshiftlocref.wlength=10;
    % Set tuning constant to use insde Huber rho function
    lshiftlocref.huberc=1.5;
    close all
    % Estimate the parameters
    [out, varargout]=LTSts(y1,'model',model,'nsamp',500,...
       'plots',2,'lshiftlocref',lshiftlocref,'msg',0);
    % Using the notation of the paper RPRH: A=2, B=4, G=2 and $\delta_1>0$.
    str=strcat('A=2, B=4, G=2, $\delta_2=',num2str(out.posLS),'$');
    numpar = {'model parameters:' , str};
    title(findobj('-regexp','Tag','LTSts:ts'),numpar,'interpreter','LaTeX','FontSize',16);

    % generate the wedgeplot
    % wedgeplot(out,'transpose',true,'extradata',[y1 out.yhat]);

%}

%{

    %% Examples 4 and 5 used in the paper RPRH: trade data.
    close all; clear all;
    % the datasets
    load('P12119085');
    load('P17049075');
    Y4 = P12119085{:,1};
    Y5 = P17049075{:,1};

    % the model
    model           = struct;
    model.trend     = 1;
    model.seasonal  = 102;
    model.s         = 12;
    model.lshift    = 14:length(Y4)-13;

    % LTSts
    out4 = LTSts(Y4,'model',model,'plots',0,'dispresults',true,'msg',0);
    out5 = LTSts(Y5,'model',model,'plots',0,'dispresults',true,'msg',0);

    % the wedgeplot with the time series and the detected outliers and
    % level shift
    wedgeplot(out4,'extradata',Y4,'titl','P12119085, imports of plants from KN to UK');
    wedgeplot(out5,'extradata',Y5,'titl','P17049075, imports of sugars from UA to LT');

    % Forecasts with a 99.9 per cent confidence level
    nfore=10;
    outfore4 = forecastTS(out4,'model',model,'nfore',nfore,'conflev',0.999,'titl','LTSts forecast for P12119085, imports of plants from KN to UK');
    outfore5 = forecastTS(out5,'model',model,'nfore',nfore,'conflev',0.999,'titl','LTSts forecast for P17049075, imports of sugars from UA to LT');

    % Comparing with FS (needs conflev option)

    outLTS4 = LTSts(Y4,'model',model,'plots',1,'conflev',0.99,'msg',0);
    title(findobj(gcf,'Tag','LTSts:ts'),'P12119085, LTS with conflev=0.99');
    
    outFRS4 = FSRts(Y4,'model',model,'plots',1);
    title(findobj(gcf,'Tag','FSRts:ts'),'P12119085, FS with default conflev');

    outLTS5 = LTSts(Y5,'model',model,'plots',1,'conflev',0.99,'msg',0);
    title(findobj(gcf,'Tag','LTSts:ts'),'P17049075, LTS with conflev=0.99');

    outFRS5 = FSRts(Y5,'model',model,'plots',1);
    title(findobj(gcf,'Tag','FSRts:ts'),'P17049075, FS with default conflev');

%}

%% Beginning of code

% Input parameters checking

% setting global variable yin
yin = y;

% Extract size of the data
T = length(yin);

% seq is the vector which will contain linear time trend
seq   = (1:T)';
one   = ones(T,1);
zerT1 = false(T,1);

if nargin<1
    error('FSDA:LTSts:MissingInputs','Input time series is missing');
end

% Set up values for default model
modeldef         =struct;
modeldef.trend   =1;        % linear trend
modeldef.s       =12;       % monthly time series
modeldef.seasonal=1;        % just one harmonic
modeldef.X       =[];       % no extra explanatory variable
modeldef.lshift  =0;        % no level shift
modeldef.ARp     =0;        % no autoregressive component

% h to be implemented for LTS

% Set the default value for h (the default is 75 per cent of the data)
hdef    = round(0.75*T);
hmin    = floor(0.5*T);
bdpdef  = 1-hdef/T;
nsampdef= 1000;

% default value for ALS iterations
reftolALSdef   = 1e-03;
refstepsALSdef = 50;

% default values for structure which contains the parameters associated
% with local level shift refinement
lshiftlocrefdef         = struct;
lshiftlocrefdef.wlength = 15;
lshiftlocrefdef.typeres = 1;
lshiftlocrefdef.huberc  = 2;

% nbestindexesdef is a positive integer which specifies how many indices of
% the smallest values of the target functions we want to retain.
nbestindexesdef=3;

% dispresultsdef Boolean about display results.
dispresultsdef=false;



%% User options

% singsub= scalar which will contain the number of singular subsets which
% are extracted (that is the subsets of size p which are not full rank)
singsub=0;

% initialize brob which will be the vector of estimated robust regression
% coefficients
brob=-99*ones(T,1);
chktrim=1;

if coder.target('MATLAB')
    options=struct('intercept',true,'lts','','nsamp',nsampdef,'h',hdef,...
        'bdp',bdpdef,'plots',0,'model',modeldef,...
        'conflev',0.975,'msg',true,'yxsave',false,...
        'SmallSampleCor',2,'nocheck',false,...
        'reftolALS',reftolALSdef,'refstepsALS',refstepsALSdef,...
        'lshiftlocref',lshiftlocrefdef,'nbestindexes',nbestindexesdef,...
        'dispresults',dispresultsdef);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:LTSts3:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present in
        % structure options Remark: the nocheck option has already been dealt
        % by routine chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:LTSts:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
        
        % Extract the names of the optional arguments
        chklist=varargin(1:2:length(varargin));
        
        % Check whether the user has selected both h and bdp.
        chktrim=sum(strcmp(chklist,'h')+2*strcmp(chklist,'bdp'));
        if chktrim ==3
            error('FSDA:LTSts:TooManyArgs','Both input arguments bdp and h are provided. Only one is required.')
        end
    end
    
end
% Write in structure 'options' the options chosen by the user
for i=1:2:length(varargin)
    options.(varargin{i})=varargin{i+1};
end


% Default values for the optional parameters are set inside structure
% 'options'

if coder.target('MATLAB')
    if ~isequal(options.model,modeldef)
        fld=fieldnames(options.model);
        
        % Check if user options inside options.model are valid options
        chkoptions(modeldef,fld)
        for i=1:length(fld)
            modeldef.(fld{i})=options.model.(fld{i});
        end
    end
    
    model = modeldef;
else
    model=options.model;
end

% Get model parameters
s        = model.s;           % get periodicity of time series
trend    = model.trend;       % get kind of  trend
seasonal = model.seasonal;    % get number of harmonics
lshift   = model.lshift;      % get level shift

% nbestindexes = indexes of the best  nbestindexes solutions for each
% tentative position of level shift
nbestindexes=options.nbestindexes;

% Check if the optional user parameters are valid.
if s <=0
    error('FSDA:LTSts:WrongInput','s= %.0f is the periodicity of the time series (cannot be negative or 0)',s)
end

if isempty(intersect(trend,0:3))
    error('FSDA:LTSts:WrongInput','Trend must assume the following values: 0  1 or 2 or 3')
end

% Construct the matrices which are fixed in each step of the minimization
% procedure
Seq = [one seq seq.^2 seq.^3];

% Define matrix which contains linear quadratic of cubic trend
intercept=options.intercept;
if intercept ==true
    Xtrend = Seq(:,1:trend+1);
else
    Xtrend = Seq(:,2:trend+1);
end
ntrend = size(Xtrend,2);

% seasonal component
yhatseaso=0;
if seasonal >0
    sstring=sprintf('%.0f',seasonal);
    % sstring=num2str(seasonal); TODO
    if seasonal>100
        varampl=real(str2double(sstring(1)));
        seasonal=real(str2double(sstring(2:3)));
    else
        varampl=0;
    end
    
    if seasonal < 1 || seasonal >floor(s/2)
        stoprint=floor(s/2);
        error('FSDA:LTSts:WrongInput','Seasonal component must be an integer between 1 and %.0f', stoprint)
    end
    
    Xseaso=zeros(T,seasonal*2);
    for j=1:seasonal
        Xseaso(:,2*j-1:2*j)=[cos(j*2*pi*seq/s) sin(j*2*pi*seq/s)];
    end
    % Remark: when s is even the sine term disapperas for j=s/2 and so the
    % maximum number of trigonometric terms is s-1
    if seasonal==(s/2)
        Xseaso=Xseaso(:,1:end-1);
    end
    nseaso=size(Xseaso,2);
else
    nseaso=0;
    varampl=0;
    Xseaso=[];
end

X = model.X;

% Order of the autoregressive component

ARp=model.ARp;
ARp=ARp(1);
if ARp>6
    disp('Number of autoregressive component is too big and can create model instability: it is set to 6');
    ARp=6;
end
if ARp>0
    % Ylagged = matrix which contains lagged values of Y
    Ylagged=zeros(T,ARp);
    for j=1:ARp
        Ylagged(:,j)=[y(1:j); y(1:end-j)];
    end
    X=[Ylagged X];
end

% nexpl = number of potential explanatory variables
isemptyX=isempty(X);
if isemptyX
    nexpl=0;
else
    nexpl=size(X,2);
end

% pini = number of parameters in the linear model without level shifts nor
% varying amplitude
% ntrend = number of trend parameters,
% nseaso = number of parameters associated with the harmonics,
% nexpl = number of explanatory variables,
pini=ntrend+nseaso+nexpl;

% p = total number of parameters in the model
% nini +
% varampl = number of parameters involving time varying trend,
% + 2 additional parameters if there is a level shift component
lshiftYN=0;
if lshift(1)~=0
    lshiftYN=1;
end
p=pini+varampl+lshiftYN*2;


% lshift=-1 is not valid in MATLAB C coder
if coder.target('MATLAB')
    % if lshift=-1, then tentative level shifts are considered for positions
    % p+1, p+2, ..., T-p-1
    if length(lshift)==1 && lshift==-1
        lshift=(p+1):(T-p);
    end
end

% indexes of linear part of seasonal component
if seasonal <6
    indlinsc=(trend+2):(trend+1+seasonal*2);
else
    indlinsc=(trend+2):(trend+1+seasonal*2-1);
end

otherind=setdiff(1:p,indlinsc);
if lshiftYN==1
    otherind=otherind(1:end-1);
end

% If the number of all possible subsets is <10000 the default is to extract
% all subsets otherwise just 10000. Notice that we use bc, a fast version
% of nchoosek. One may also use the approximation
% floor(exp(gammaln(n+1)-gammaln(n-p+1)-gammaln(pini+1))+0.5)
ncomb=bc(T,pini);


% And check if the optional user parameters are reasonable.

% Check h and bdp The user has only specified h: no need to specify bdp.
if chktrim==1
    if options.h < hmin
        error('FSDA:LTSts:WrongInput',['The LTS must cover at least ' int2str(hmin) ' observations.'])
    elseif options.h >= T
        error('FSDA:LTSts:WrongInput','h is greater or equal to the number of non-missings and non-infinites.')
    end
    bdp=1-options.h/T;
    
    % the user has only specified bdp: h is defined accordingly
elseif chktrim==2
    bdp=options.bdp;
    if bdp <= 0
        error('FSDA:LTSts:WrongInput','Attention: bdp should be larger than 0');
    end
    
    nalpha=floor(T*(1-bdp));
    options.h=nalpha;
else
    bdp=-99;
end

% Check number of subsamples to extract
if options.nsamp>ncomb
    if options.msg==true
        disp('Number of subsets to extract greater than (n p). It is set to (n p)');
    end
    options.nsamp=0;
elseif  options.nsamp<0
    error('FSDA:LTSts:WrongInput','Number of subsets to extract must be 0 (all) or a positive number');
end


h=floor(options.h);         % Number of data points on which estimates are based
nsamp=options.nsamp;        % Number of subsets to extract
nsampsubsequentsteps=round(nsamp/2);

lts=options.lts;
SmallSampleCor=options.SmallSampleCor; % small sample correction factor
if varampl>0
    % Convergence criteria inside ALS loop
    reftolALS=options.reftolALS;
    refstepsALS=options.refstepsALS;
    if coder.target('MATLAB')
        verLess2016b=verLessThanFS(9.1);
    else
        verLess2016b=false;
    end
else
    verLess2016b=false;
    reftolALS=0;
    refstepsALS=0;
end

constr=0;

if ~isstruct(lts) && isempty(lts)
    refsteps=2;
    reftol=1e-6;
    bestr=20;
    refstepsbestr=50;
    reftolbestr=1e-8;
    
elseif isstruct(lts)
    if coder.target('MATLAB')
        ltsdef.refsteps=2;
        ltsdef.reftol=1e-6;
        ltsdef.bestr=20;
        ltsdef.refstepsbestr=50;
        ltsdef.reftolbestr=1e-8;
        
        % Control the appearance of the trajectories to be highlighted
        if ~isequal(lts,ltsdef)
            
            fld=fieldnames(lts);
            
            % Check if user options inside options.fground are valid options
            chkoptions(ltsdef,fld)
            for i=1:length(fld)
                ltsdef.(fld{i})=lts.(fld{i});
            end
        end
        
        % For the options not set by the user use their default value
        lts=ltsdef;
    end
    refsteps=lts.refsteps;
    reftol=lts.reftol;
    bestr=lts.bestr;
    refstepsbestr=lts.refstepsbestr;
    reftolbestr=lts.reftolbestr;
else
    error('FSDA:LTSts:WrongInput','Input option lts must be a structure or a empty value');
end


conflev=options.conflev;    % Confidence level which is used for outlier detection
msg=options.msg;            % Scalar which controls the messages displayed on the screen

if coder.target('MATLAB')
    % Get user values of warnings
    warnrank=warning('query','MATLAB:rankDeficientMatrix');
    warnsing=warning('query','MATLAB:singularMatrix');
    warnnear=warning('query','MATLAB:nearlySingularMatrix');
    % Set them to off inside this function at the end of the file they will be
    % restored to previous values
    warning('off','MATLAB:rankDeficientMatrix');
    warning('off','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix');
end

if lshiftYN==1
    % If a level shift is present, it is necessary to
    % reestimate a linear model each time with a different
    % level shift and, if  take the one which minimizes the target
    % function (residual sum of squares/2 = negative log
    % likelihood)
    LSH = lshift;
    % total number of subsets to pass to procedure subsets
    ncombLSH=bc(T-1,pini+1);
    
else
    LSH=0;
end


% ScaleLSH= estimate of the squared scale for each value of LSH which has been
% considered
numscale2LSH=[LSH' zeros(length(LSH),2)];

% yhatrobLSH = vector of fitted values for each value of LSH
yhatrobLSH=zeros(T,length(LSH));

% ilsh is a counter which is linked to the rows of LSH
ilsh=0;

% lLSH = length of tentative level shift positions
lLSH=length(LSH);


bestrdiv2=round(bestr/2);

% allnumscale2 will contain the best best estimates of the target function
% for a tentative value of level shift position
allnumscale2=zeros(bestr,1);

% Store all bestr target functions for each tentative level shift
% position (target function = truncated residual sum of squares)
ALLnumscale2=zeros(bestr,lLSH);

% Store the position of the indexes occupying nbestindexes best solutions of target
% function for each tentative level shift position
% 1-bestrdiv2       = solutions from fresh subsets.
% bestrdiv2+1-bestr = best solutions coming from previous tentative level
% shift position
NumScale2ind=zeros(nbestindexes,lLSH);


% Weights = units forming subset for the solution associated to the minimum
% scale for each value of LSH
Weights=false(T,lLSH);

brobLSH=zeros(p,lLSH);

% Construct matrix X (called Xsel) which contains the linear part of the model
if seasonal==0
    if isemptyX
        Xsel=Xtrend;
    else
        Xsel=[Xtrend X];
    end
else
    if isemptyX
        Xsel=[Xtrend Xseaso];
    else
        Xsel= [Xtrend Xseaso X];
    end
    % zero for varampl is automatically included because b0 is
    % initialized as a vector of zeroes b0=[b0;zeros(varampl,1)];
end

% WEIisum = matrix which will contain the number of times each units has
% been included into the best h-subset after two iterations
WEIisum=zeros(T,lLSH);

% WEIisumbest10 = matrix which will contain the number of times each units has
% been included into the best h-subsets among the bestr/2 best
WEIibest10sum=zeros(T,lLSH);
WEIibestrdiv2=zeros(T,bestrdiv2);

RES = nan(T,lLSH);

% Consistency factor based on the variance of the truncated normal
% distribution. 1-h/n=trimming percentage Compute variance of the truncated
% normal distribution.
a=norminv(0.5*(1+h/T));
%factor=1/sqrt(1-(2*a.*normpdf(a))./(2*normcdf(a)-1));
factor=1/sqrt(1-2*(T/h)*a.*normpdf(a));

% Initialize 2D or 3D array which stores indexes of extracted
% subsets for each tentative level shift position
if nargout>1
    Ccell=cell(lLSH,1);
    % Initialization of Ccell is necessary for MATLAB C coder
    zerlSH=zeros(T,pini);
    for i=1:lLSH
        Ccell{i}=zerlSH;
    end
end

if ~coder.target('MATLAB')
    bestyhattoadd=zeros(pini,pini);
    bestbetastoadd=bestyhattoadd;
    bestsubsettoadd=bestyhattoadd;
    ncombLSH=0;
    bsb=0;
    ibest=0;
    yhatrob=0;
    weightsst=false;
    posLS=0;
    Likloc=0;
end



for lsh=LSH
    ilsh=ilsh+1;
    
    sworst = Inf;
    
    % Xlshift = explanatory variable associated with
    % level shift Xlshift is 0 up to lsh-1 and 1 from
    % lsh to T
    Xlshift= [zeros(lsh-1,1);ones(T-lsh+1,1)];
    
    
    if ilsh>1
        nsamp=nsampsubsequentsteps;
        bestrLSH=bestrdiv2;
        bestnumscale2 = Inf * ones(bestrdiv2,1);
        bestbetas = zeros(bestrdiv2,p);
        bestyhat=zeros(T,bestrdiv2);
        bestsubset = zeros(bestrdiv2,pini+lshiftYN*2);
        
    else
        bestbetas = zeros(bestr,p);
        bestyhat=zeros(T,bestr);
        bestsubset = zeros(bestr,pini+lshiftYN*2);
        bestnumscale2 = Inf * ones(bestr,1);
        bestrLSH=bestr;
    end
    
    if lsh>0
        
        [Cini,nselected] = subsets(nsamp,T-1,pini+1,ncombLSH,msg);
        
        C=[lsh*ones(nselected,1) zeros(nselected,pini+1)];
        
        
        % Make sure that observation lsh is always included in the subset
        % and that the subset contains at least one unit smaller than lsh
        for r=1:nselected
            Cr=Cini(r,:);
            % Observations greater or equal than lsh will be increased by one
            boo=Cr>=lsh;
            Cr(boo)=Cr(boo)+1;
            % Make sure there is at least one observation smaller than lsh
            boo=Cr<lsh;
            % if sum(boo)==0 then in the subset there is no observation
            % which is smaller than lsh
            if sum(boo)<1
                Cr(1)=randsample(lsh-1,1);
            end
            C(r,2:end)=Cr;
        end
        
    else
        % If there is no level shift component
        [Cini,nselected] = subsets(nsamp,T,pini,ncomb,msg);
        C=Cini;
    end
    
    % Store indexes of extracted subsets if nargout is greater than 1
    if nargout>1
        Ccell{ilsh}=Cini;
    end
    
    % yhatall= matrix which will contain fitted values for each extracted
    % subset
    % yhatall=zeros(T,nselected);
    
    % WEIi = matrix which will contain indication of the units forming best
    % h subset. Each column refers to a subset
    WEIi=zeros(T,nselected);
    
    % ij is a scalar used to ensure that the best first bestr solutions are
    % stored in order to be brought to full convergence
    % subsets are stored
    ij=1;
    
    % Loop through all nselected subsamples
    for i=1:nselected
        % Initialize b0 as vector of zeroes for each subset
        % The order of the elements of b0 is as follows
        % 1) trend elements (if present). If the trend is order two r are
        % r+1 coefficients if the intercept is present otherwise there are
        % just r components (Xtrend)
        % 2) linear part of seasonal component 2, 4, 6, ..., s-2, s-1 coefficients
        % (if present)   (Xseaso)
        % 3) coefficients associated with the matrix of explanatory variables
        % which have a potential effect on the time series under study (X)
        % 4) non linear part of seasonal component, that is varying amplitude.
        % If varying amplitude is of order k there are k coefficients (if
        % present) (Seq)
        % 5) level shift component (if present). In this case there are two
        % coefficients, the second (which is also the last element of
        % vector b0) is an integer which specifies the time in which level
        % shift takes place and the first (which is also the penultime
        % element of vector b0) is a real number which identifies the
        % magnitude of the upward (downward) level shift (Xlshift)
        beta0=zeros(p,1);
        
        % extract a subset of size p
        index = C(i,:);
        
        if lsh==0
            Xlshift=[];
        end
        
        Xfinal=[Xsel Xlshift];
        % Preliminary OLS estimates (including tentative level shift) based
        % just on the units forming subset
        bsb=index(:);
        betaini=Xfinal(bsb,:)\yin(bsb);
        
        % Check if betaini contains NaN
        if ~any(isnan(betaini))
            
            % The first pini components are associated with
            % trend and seasonal (without varying
            % amplitude) and explanatory variables
            beta0(1:pini)=betaini(1:pini);
            
            if lsh>0
                % The last two components of beta0 are the associated with
                % level shift. More precisely penultimate position is for the
                % coefficient of level shift and, final position is the integer
                % which specifies the starting point of level shift
                beta0(end-1:end)=[betaini(end) lsh];
            end
            
            if varampl>0
                if verLess2016b==1
                    [betaout]=ALSbsxfun(beta0);
                else
                    [betaout]=ALS(beta0);
                end
                %  betaoutCHK=ALSbsxfun(beta0);
                % dd=1;
            else
                betaout=beta0;
                
                %disp(['lsh' num2str(lsh)])
                %disp(beta0)
                %disp('------')
            end
            
            % Compute  fitted values (for all units). Therefore recall function
            % lik but this time computed using all observations
            bsb=seq;
            % Procedure lik computes yhat (fitted values for all the
            % observations using parameter estimates based on bsb). vector yhat
            % will be used inside procedure IRWLSreg as starting value of the
            % iterations (concentration steps)
            lik(betaout);
            beta=betaout;
            
            % 1(a) ii. -  Now apply concentration steps
            tmp = IRWLSreg(yin,beta,refsteps,reftol,h);
            
            % Store weights
            WEIi(:,i)=tmp.weights;
            
            % Store fitted values for each subset
            % yhatall(:,i)=tmp.yhat;
            
            betarw = tmp.betarw;
            numscale2rw = tmp.numscale2rw;
            
            % 1(c) Consider only the 10 subsets that yield the lowest objective
            % function so far.
            if ij > bestrLSH
                
                if numscale2rw < sworst
                    
                    % Store numscale2rw, betarw and indexes of the units
                    % forming the best subset for the current iteration
                    
                    % Find position of the maximum value of previously
                    % stored best numerator of squared scaled
                    [~,ind] = max(bestnumscale2);
                    
                    bestnumscale2(ind) = numscale2rw;
                    bestbetas(ind,:)   = betarw';
                    bestsubset(ind,:)  = index;
                    bestyhat(:,ind)    = yhat;
                    % sworst = best scale among the bestr found up to now
                    sworst             = max(bestnumscale2);
                end
            else
                
                bestnumscale2(ij)  = numscale2rw;
                bestbetas(ij,:) = betarw';
                bestsubset(ij,:)= index;
                bestyhat(:,ij)=yhat;
                % sworst = best scale among the bestr found up to now
                sworst = max(bestnumscale2);
                ij = ij+1;
                brob = 1;
            end
        end
    end
    
    if brob(1)==-99
        error('FSDA:LTSts:NoFullRank','No subset had full rank. Please increase the number of subsets or check your design matrix X')
    else
    end
    
    % Store for each tentative level shift the number of times each unit
    % belonged to the best subset
    WEIisum(:,ilsh)=sum(WEIi,2);
    
    % 1 (b)
    % With the 0 subsets that yield the lowest objective function so far.
    % Apply C-steps to these until full convergence.
    
    % perform C-steps on best 'bestr' solutions, till convergence or for a
    % maximum of refstepsbestr steps using a convergence tolerance as
    % specified by scalar reftolbestr
    
    
    % If ilsh >1 it is necessary also to consider the 10 best solutions from
    % step j-1
    if ilsh==1
        bestyhatall=bestyhat;
        bestbetasall=bestbetas;
        bestsubsetall=bestsubset;
    else
        bestyhatall=[bestyhat bestyhattoadd];
        bestbetasall=[bestbetas; bestbetastoadd];
        bestsubsetall=[bestsubset; bestsubsettoadd];
    end
    
    % numsuperbestscale2 = numerator of estimate of super best squared
    % scale
    numsuperbestscale2 = Inf;
    
    % Just to have an idea about y and yhat for a particular lsh value
    % plot([y bestyhat(:,1)])
    
    
    for i=1:bestr
        yhat=bestyhatall(:,i);
        tmp = IRWLSreg(yin,bestbetasall(i,:)',refstepsbestr,reftolbestr,h);
        
        % Store information about the units forming best h subset among the
        % 10 best
        WEIibestrdiv2(:,i)=tmp.weights;
        
        allnumscale2(i,1)=tmp.numscale2rw;
        % allscales(i,2)=tmp.betarw(end);
        
        if tmp.numscale2rw < numsuperbestscale2
            % brob = superbestbeta
            brob = tmp.betarw;
            % bs = superbestsubset, units forming best subset according to
            % fastlts
            % bs = bestsubsetall(i,:);
            yhatrob=tmp.yhat;
            numsuperbestscale2=tmp.numscale2rw;
            ibest=i;
            weightsst=tmp.weights;
        end
    end
    
    % Store the bestrdiv2 best values of target function
    [~,numscale2ssorind]=sort(allnumscale2);
    bestyhattoadd=bestyhatall(:,numscale2ssorind(1:bestrdiv2));
    bestbetastoadd=bestbetasall(numscale2ssorind(1:bestrdiv2),:);
    % The last element of estimated beta coefficients is the point in
    % which level shift takes place. This has to be increased by one
    % unit. Please note that betas are stored in rows therefore we have
    % to change the last column
    bestbetastoadd(:,end)=bestbetastoadd(:,end)+1;
    
    bestsubsettoadd=bestsubsetall(numscale2ssorind(1:bestrdiv2),:);
    
    numscale2LSH(ilsh,2:3)=[numsuperbestscale2 ibest];
    yhatrobLSH(:,ilsh)=yhatrob;
    brobLSH(:,ilsh)=brob;
    
    % plot(seq,[y yhatrob])
    % title(['Level shift in step t=' num2str(LSH(ilsh))])
    ALLnumscale2(:,ilsh)=allnumscale2;
    
    scaledres = (yin-yhatrob)/sqrt(numsuperbestscale2/h);
    RES(:,ilsh) = scaledres;
    
    
    weightsst = (weightsst | abs(scaledres)<2.58*factor);
    % disp(sum(weightsst))
    Weights(:,ilsh) = weightsst;
    
    % Store the indexes among the bestr best, forming the bestrdiv2 best
    % estimates of the target function (target function = numerator of
    % squared scale)
    NumScale2ind(:,ilsh)=numscale2ssorind(1:nbestindexes);
    
    WEIibest10sum(:,ilsh)=sum(WEIibestrdiv2,2);
    if lsh>0 && msg ==true
        fprintf('Level shift for t=%.0f\n',lsh)
    end
end

% save RES to output structure (these residuals can be used for example to
% prouduce the double wedge plot, see function wedgeplot for more details)
out.RES = RES;

Weimod=double(Weights);
for j=1:size(Weimod,2)
    boo=Weimod(:,j)==1;
    Weimod(boo,j)=seq(boo);
    Weimod(~boo,j)=NaN;
end

% Store units forming best h subset
out.Hsubset=Weimod;

[~,minidx]=min(numscale2LSH(:,2));
brobbest=brobLSH(:,minidx);

% Pass from numerator of squared estimate of the scale to proper scale
% estimate
sh0=sqrt(numscale2LSH(minidx,2)/h);

% Consistency factor
s0=sh0*factor;

% Apply small sample correction factor of Pison et al.
s0=s0*sqrt(corfactorRAW(1,T,h/T));

if  lsh>0
    % Compute the residuals locally just changing the position of the level
    % shift
    bstar=brobbest;
    
    lshiftlocref=options.lshiftlocref;
    if isfield(lshiftlocref,'wlength')
        k=lshiftlocref.wlength;
    else
        k=15;
    end
    
    if isfield(lshiftlocref,'typeres')
        typeres=lshiftlocref.typeres;
    else
        typeres=1;
    end
    
    if isfield(lshiftlocref,'huberc')
        huberc=lshiftlocref.huberc;
    else
        huberc=2;
    end
    
    tloc=bstar(end)-k:bstar(end)+k;
    % Reduce width of tloc dinamically
    LSHmin=min(LSH);
    LSHmax=max(LSH);
    % make sure that tloc is in the range LSHmin and LSHmax
    while (max(tloc)>LSHmax  || min(tloc)<LSHmin )
        if k==0
            break
        end
        k=k-1;
        tloc=bstar(end)-k:bstar(end)+k;
    end
    
    
    bsb=tloc(:);
    Likloc=[tloc' zeros(length(tloc),3)];
    ij=0;
    
    for j=tloc(1):tloc(end)
        ij=ij+1;
        btmp=bstar;
        btmp(end)=j;
        
        Xlshift= [zeros(j-1,1);ones(T-j+1,1)];
        
        lik(btmp);
        
        resbsb=(yin(bsb)-yhat)/sh0;
        Likloc(ij,2)=sum((HUrho(resbsb,huberc)).^2);
        Likloc(ij,3)=sum((yin(bsb)-yhat).^2);
        
    end
    % Use Huberized residual sum of squares to find minimum
    [~,locmin]=min(Likloc(:,typeres+1));
    posLS=Likloc(locmin,1);
    Xlshift=[zeros(posLS-1,1);ones(T-posLS+1,1)];
    brobfinal=bstar;
    brobfinal(end)=posLS;
else
    brobfinal= brobbest;
end
bsb=seq;

% Compute fitted values using final estimate of beta for all the
% observations
lik(brobfinal);

% REWEIGHTING STEP

% residuals = Raw residuals using final estimate of beta
residuals=yin-yhat;

% Find the units with the smallest absolute p+1 residuals (before
% reweighting step)
[~,IndBestRes]=sort(abs(residuals));
nofullrank=true;
bs=IndBestRes(1:p+1);
ij=0;

while nofullrank
    bs=IndBestRes(1:p+ij);
    if rank(zscore(Xsel(bs,2:end)))<pini-1
        ij=ij+1;
    else
        nofullrank = false;
    end
end

%if the robust s0 is too small, compute it with a set of different methods:
%Qn, Sn, std and the interquantile difference for increasing percentages
%([0.25-0.75], [0.26-0.76], ...)
if abs(s0) < 1e-7
    if msg==true
        disp('Attention: there was an exact fit. Robust estimate of s^2 is <1e-7')
    end
    [~,~,s0]=zscoreFS(residuals,'median','Qn');
    if s0==0
        [~,~,s0]=zscoreFS(residuals,'median','Sn');
    end
    if s0==0
        [~,~,s0]=zscoreFS(residuals,'median','std');
    end
    if s0==0
        absr=abs(residuals);
        for j=75:1:99
            s0=prctile(absr,j)- prctile(absr,100-j);
            if s0>0
                break
            end
        end
    end
    
    %weights = abs(residuals)<=1e-7;
    % stdres = residuals/s0;
else
end

stdres = residuals/s0;
if SmallSampleCor==1
    plinear=pini+lshiftYN;
    robest='LTS';
    eff=[];
    rhofunc='';
    sizesim=0;
    Tallis=1;
    if T<50
        Ttouse=50;
    else
        Ttouse=T;
    end
    
    if bdp==-99
        bdp=1-options.h/T;
    end
    thresh=RobRegrSize(Ttouse,plinear,robest,rhofunc,bdp,eff,sizesim,Tallis);
    extracoeff=sqrt(thresh/chi2inv(0.99,1));
    weights = abs(stdres)<=sqrt(chi2inv(0.99,1))*extracoeff;
    
elseif  SmallSampleCor==2
    weights=GYfilt(stdres,'iterating',false,'alpha',0.99,'centering',true,'niter',10);
elseif  SmallSampleCor==3
    weights=GYfilt(stdres,'iterating',true,'alpha',0.99,'centering',true,'niter',10);
elseif SmallSampleCor==4
    weights = abs(stdres)<=sqrt(chi2inv(0.99,1));
else
    error('FSDA:ltsTS:WrongInputOpt','wrong small sample cor factor')
end


% else
%     % There is an approximate perfect fit for the first h observations.
%     % We consider as outliers all units with residual greater than 1e-7.
%     weights = abs(residuals)<=1e-7;
%
%     %     % Store the weights
%     %     out.weights=weights;
%
%
%     % s is set to 0
% %     s0=0;
% %
% %     % Standardized residuals are artificially set equal to raw residuals.
% %     stdres=residuals;
% end


% weights is a boolean vector.
bsb=seq(weights);

% Store bsb to use in order to find sum of squares of residuals for
% reduced model.
bsbModSel=bsb;

% Find new estimate of beta using only observations which have
% weight equal to 1. Notice that new brob overwrites old brob
% computed previously.


if varampl==0 && lshiftYN==0 % In this case the model is linear
    % Function lik constructs fitted values and residual sum of
    % squares
    betaout = Xsel(bsb,:) \ yin(bsb);
    % update fitted values
    yhat = Xsel * betaout;
    
    % find fitted values using all observations
    yhat =  Xsel * betaout;
    s2=sum((yin(bsb)-yhat(bsb)).^2)/(h-size(Xsel,2));
    invXX=inv(Xsel'*Xsel);
    covB=s2*invXX; %#ok<MINV>
    Xlin=Xsel;
    
elseif   varampl==0 && lshiftYN==1
    % In this case there is just level shift however we do not redo
    % the non linear estimation but a simple LS
    
    Xseldum=[Xsel Xlshift];
    betaout = Xseldum(bsb,:) \ yin(bsb);
    
    % find fitted values using all observations
    yhat =  Xseldum * betaout;
    s2=sum((yin(bsb)-yhat(bsb)).^2)/(h-size(Xseldum,2));
    invXX=inv(Xseldum(bsb,:)'*Xseldum(bsb,:));
    covB=s2*invXX; %#ok<MINV>
    Xlin=Xseldum;
else % model is non linear because there is time varying amplitude in seasonal component
    Xtrendf=Xtrend(bsb,:);
    Xseasof=Xseaso(bsb,:);
    if ~isempty(X)
        Xf=X(bsb,:);
    end
    Seqf=Seq(bsb,:);
    yf=yin(bsb);
    
    % Find new estimate of scale using only observations which have
    % weight equal to 1.
    weights=false(T,1);
    weights(bsb)=true;
    
    if coder.target('MATLAB')
        
        if lshiftYN==1
            Xlshiftf=Xlshift(bsb);
            [betaout,~,Xlin,covB,MSE,~]  = nlinfit(Xtrendf,yf,@likyhat,brobfinal(1:end-1));
        else
            Xlshiftf=0;
            [betaout,~,Xlin,covB,MSE,~]  = nlinfit(Xtrendf,yf,@likyhat,brobfinal);
            % [betaout,R,J,covB,MSE,ErrorModelInfo]  = nlinfit(Xtrendf,yf,@likyhat,brobfinal);
            % Note that MSE*inv(J'*J) = covB
        end
        
        % yfitFS = likyhat(betaout,Xtrendf);
        % nans=false(length(yfitFS),1);
        % sqweights=ones(length(yfitFS),1);
        % fdiffstep=1.0e-05*0.6655*ones(length(betaout),1);
        % J = getjacobianFS(betaout,fdiffstep,@likyhat,yfitFS,nans,sqweights);
    else  % MATLAB CCODER PART nlinfit replaced by lsqcurvefit
        optionsLSQcurvefit = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt');
        if lshiftYN==1
            Xlshiftf=Xlshift(bsb);
            ub=Inf(length(brobfinal)-1,1);
            lb=-ub;
            [betaout,~,~,~,~,~,Xlin] = lsqcurvefit(@likyhat,brobfinal(1:end-1),Xtrendf,yf,lb,ub,optionsLSQcurvefit);
        else
            Xlshiftf=0;
            % [betaoutCHK,resnorm,residual,exitflag,output,lambda,XlinCHK]= lsqcurvefit(@likyhat,brobfinal,Xtrendf,yf);
            ub=Inf(length(brobfinal),1);
            lb=-ub;
            [betaout,~,~,~,~,~,Xlin] = lsqcurvefit(@likyhat,brobfinal,Xtrendf,yf,lb,ub,optionsLSQcurvefit);
        end
        covB=eye/length(betaout);
        MSE=1;
        %         MSE=(residuals'*residuals)/length(betaout);
        %         covB=MSE*inv(XlinCHK'*XlinCHK)*MSE;
        
    end
    invXX=covB/MSE;
    
    % Now compute again vector yhat using final vector betaout
    bsb=seq;
    lik(betaout);
    
end

% Store beta standard error, t stat and p values
sebetaout=sqrt(diag(covB));
tout=betaout./sebetaout;
dfe=T-length(betaout);
pval=2*(tcdf(-abs(tout), dfe));
B=[betaout sebetaout tout pval];

out.B=B;

if lsh>0
    % Store position of level shift
    out.posLS=posLS;
else
    if ~coder.target('MATLAB')
        out.posLS=[];
    end
end

% Computation of reweighted residuals.
residuals=yin-yhat;

% s2full
s2full=residuals(bsbModSel)'*residuals(bsbModSel);

% s0 =sqrt(MSE)
s0=sqrt(sum(weights.*residuals.^2)/(sum(weights)-1));
% Compute new standardized residuals.

% Apply consistency factor to reweighted estimate of sigma
hrew=sum(weights);
if hrew<T
    % Make sure that hrew has at least T/2 observations
    if hrew<T/2
        hrew=T/2;
    end
    
    % factor=consistencyfactor(hrew,n,1);
    a=norminv(0.5*(1+hrew/T));
    %factor=1/sqrt(1-(2*a.*normpdf(a))./(2*normcdf(a)-1));
    factor=1/sqrt(1-2*(T/hrew)*a.*normpdf(a));
    % Apply small sample correction factor to reweighted estimate
    % of sigma
    factor=factor*sqrt(corfactorREW(1,T,hrew/T));
else
    factor=1;
end

s0=s0*factor;
if s0==0
    stdres=residuals;
else
    stdres=residuals/s0;
end

% Declare as outliers the observations which have a standardized
% residual greater than cutoff. REMARK: while the first threshold
% was based on the Student T (with modified degrees of freedom), in
% this second round the threshold is based on the Normal. Notice
% that: sqrt(chi2inv(0.975,1)) = tinv(0.9875,\infinity) =
% norminv(0.9875)

out.yhat=yhat;


%% Store quantities in the out structure

outliers     = abs(stdres)>norminv((conflev+1)/2);
out.outliers = seq(outliers);

%decomment the following two lines to get outlier pvalues
p_all = normcdf(-abs(stdres));
out.outliersPval = p_all(outliers);

% Store robust estimate of s
out.scale = s0;

% Store the 20 best estimates of the scale for each tentative level shift
% which is considered
out.numscale2 = ALLnumscale2;

% Store indices forming the bestrdiv2 best estimates of the target function
out.BestIndexes = NumScale2ind;

% Store scaled residuals
out.residuals=stdres;

% Store units forming best initial subset of p-1 observations
out.bs=bs;

% Store list of units declared as outliers
% out.outliers=seq(weights==0);

% Store confidence level which is used to draw the horizontal lines in the
% plot
out.conflev=options.conflev;

% Store the number of observations that have determined the LTS (LMS)
% estimator, i.e. the value of h.
out.h=h;

% Store vector of weights (values equal to 1 are associated with units
% parteciapting to the fit)
out.weights=weights;

% Store number of singular subsets
out.singsub=singsub;
if msg==true
    if singsub/nselected>0.1
        percexcl=100*singsub/nselected;
        disp('------------------------------')
        % disp(['Warning: Number of subsets without full rank equal to ' num2str(100*singsub/nselected) '%'])
        fprintf('Warning: Number of subsets without full rank equal to %.1f%%\n',percexcl)
        
    end
end
% Store information about the class of the object

out.class='LTSts';

if lsh>0
    % Store local improvement of the likelihood
    out.Likloc=Likloc;
else
    if ~coder.target('MATLAB')
        out.Likloc=0;
    end
end

% Store response
out.y=yin;

if options.yxsave == true
    if options.intercept==true
        % Store X (without the column of ones if there is an intercept)
        out.X=Xlin(:,2:end);
    else
        out.X=Xlin;
    end
else
    if ~coder.target('MATLAB')
        out.X=[];
    end
end

out.invXX=invXX;

dispresults=options.dispresults;

% b_trend = {'b_trend1'; 'b_trend2'; 'b_trend3'; 'b_trend4'};
% b_seaso = {'b_cos1'; 'b_sin1'; 'b_cos2'; 'b_sin2'; ...
%     'b_cos3'; 'b_sin3'; 'b_cos4'; 'b_sin4'; ...
%     'b_cos5'; 'b_sin5'; 'b_cos6'};
% b_AR =    {'b_AR1'; 'b_AR2'; 'b_AR3'; 'b_AR4'; 'b_AR5'; 'b_AR6'};
% b_X  =    {'b_X1'; 'b_X2'; 'b_X3'; 'b_X4'; 'b_X5'; 'b_X6'};
% b_varampl = {'b_varampl'; 'b_varamp2'; 'b_varamp3'};
% b_lshift  = {'b_lshift' ; 't_lshift'};
%
% if ARp>0
%      b_expl=[b_AR(1:ARp); b_X(1:nexpl-ARp)];
% else
%     b_expl=b_X;
% end
%
% if seasonal>0
%     if 2*seasonal==s
%         lab=[b_trend(1:trend+1); b_seaso];
%     else
%         lab=[b_trend(1:trend+1); b_seaso(1:2*seasonal)];
%     end
% else
%     lab=b_trend(1:trend+1);
% end
%
% if nexpl>0
%     lab=[lab;b_expl(1:nexpl)];
% end
% if varampl>0
%     lab=[lab;b_varampl(1:varampl)];
%     posvarampl=length(lab)-varampl+1:length(lab);
% else
%     posvarampl=[];
% end
% if lshiftYN==1
%     lab=[lab; b_lshift(1)];
% end


b_trend = ['b_trend1'; 'b_trend2'; 'b_trend3'; 'b_trend4'];
b_seaso =['b_cos1  '; 'b_sin1  '; 'b_cos2  '; 'b_sin2  '; ...
    'b_cos3  '; 'b_sin3  '; 'b_cos4  '; 'b_sin4  '; ...
    'b_cos5  '; 'b_sin5  '; 'b_cos6  '];
b_AR =    ['b_AutoR1'; 'b_AutoR2'; 'b_AutoR3'; 'b_AutoR4'; 'b_AutoR5'; 'b_AutoR6'];
b_X  =    ['b_explX1'; 'b_explX2'; 'b_explX3'; 'b_explX4'; 'b_explX5'; 'b_explX6'];
b_varampl = ['b_varaml'; 'b_varam2'; 'b_varam3'];
b_lshift  = ['b_lshift' ; 't_lshift'];

% Code generation does not support string arrays
% b_trend = ["b_trend1"; "b_trend2"; "b_trend3"; "b_trend4"];
% b_seaso = ["b_cos1"; "b_sin1"; "b_cos2"; "b_sin2"; ...
%     "b_cos3"; "b_sin3"; "b_cos4"; "b_sin4"; ...
%     "b_cos5"; "b_sin5"; "b_cos6"];
% b_AR =    ["b_AR1"; "b_AR2"; "b_AR3"; "b_AR4"; "b_AR5"; "b_AR6"];
% b_X  =    ["b_X1"; "b_X2"; "b_X3"; "b_X4"; "b_X5"; "b_X6"];
% b_varampl = ["b_varampl"; "b_varamp2"; "b_varamp3"];
% b_lshift  = ["b_lshift" ; "t_lshift"];

if ARp>0
    b_expl=[b_AR(1:ARp,:); b_X(1:nexpl-ARp,:)];
else
    b_expl=b_X;
end

if seasonal>0
    if 2*seasonal==s
        lab=[b_trend(1:trend+1,:); b_seaso];
    else
        lab=[b_trend(1:trend+1,:); b_seaso(1:2*seasonal,:)];
    end
else
    lab=b_trend(1:trend+1,:);
end

if nexpl>0
    lab=[lab;b_expl(1:nexpl,:)];
end
if varampl>0
    lab=[lab;b_varampl(1:varampl,:)];
    posvarampl=length(lab)-varampl+1:length(lab);
else
    posvarampl=[];
end
if lshiftYN==1
    lab=[lab; b_lshift(1,:)];
end

% if verLessThan ('matlab','8.2.0')
% else
% Store matrix B in table format (with labels for rows and columns)
out.Btable=array2table(B,'RowNames',string(lab),'VariableNames',{'Coeff','SE','t','pval'});
% end

if dispresults
    if coder.target('MATLAB')
        disp(out.Btable)
    else
        
        %         bhat=out.B(:,1);
        %         se=out.B(:,2);
        %         tstat=out.B(:,3);
        %         pval=out.B(:,4);
        %         %         disp('           Coeff.     SE         t-stat       p-values');
        %         fprintf('%s%.3f',lab,bhat)
        %     if verLessThan ('matlab','8.2.0')
        %         disp('           Coeff.     SE         t-stat       p-values');
        %         disp( [char(lab) num2str([bhat se tstat pval])]);
        %     else
        % disp([table(lab) table(bhat) table(se) table(tstat) table(pval)]);
        %     end
    end
    if lshiftYN==1
        fprintf('Level shift position t=%.0f\n',posLS)
    end
end

%% Create plots

if coder.target('MATLAB')
    plots=options.plots;        % Plot of residuals equal to 1
    
    % plots = 1 generates a figure with two panels: one with the time series
    % and another with the residuals against index number; plots =2 produces
    % also a number of other informative plots; else no plot is produced.
    if plots>=1
        
        % some general plot settings
        vlt15 = verLessThan('matlab', '7.15');
        clr = 'bkrgmcy';
        syb = {'-','--','-.',':','-','--','-.'};
        FontSize    = 14;
        SizeAxesNum = 14;
        
        % slightly increase the range of the time series axis values
        mine = min(yin(:));
        maxe = max(yin(:));
        delta = (maxe-mine)*0.1;
        yaxlim = [mine - delta ; maxe + delta];
        % the next check is introduced because if the two elements of Ylim are
        % the same (which happens if the series is constant), the set (gca)
        % some lines below produce errors.
        if yaxlim(1) == yaxlim(2)
            yaxlim(2)=yaxlim(2)+0.01*yaxlim(2);
        end
        % Time series + fitted values
        figure
        htmp = subplot(2,1,1);
        plot(yin, 'Color',clr(1),'LineStyle',syb{1},'LineWidth',1);
        hold('on');
        plot(yhat,'Color',clr(2),'LineStyle',syb{2},'LineWidth',1);
        
        set(htmp,'Tag','LTSts:ts');
        %xlabel('Time','FontSize',FontSize);
        ylabel('Real and fitted values','FontSize',FontSize,'interpreter','none');
        if ~vlt15
            set(gca,'FontSize',SizeAxesNum,'Ylim' , yaxlim,'Box','on','BoxStyle','full');
        else
            set(gca,'FontSize',SizeAxesNum,'Ylim' , yaxlim,'Box','on');
        end
        xtickval = get(htmp,'XTick');
        xticklab = get(htmp,'XTickLabel');
        set(htmp,'XTickMode','manual');
        drawnow;
        
        % mark outliers with their severity
        if ~isempty(residuals)
            seq = 1:T;
            quant = sqrt(chi2inv(conflev,1));
            resboo=out.residuals(out.outliers);
            th=8;resboo(abs(resboo)>th)=th;
            %Rescale residuals in the interval [0 3]
            sizeout=3*(abs(resboo)-quant)/(th-quant);
            for i=1:length(sizeout)
                plot(seq(out.outliers(i)),yin(out.outliers(i),1),'x','LineWidth',sizeout(i),'Color','r', 'MarkerFaceColor','k');
            end
        end
        
        % plot the vertical line of the level shift position and the associated
        % label on the X axis
        if isfield(out,'posLS') && ~isempty(out.posLS)
            line(out.posLS*ones(2,1) , yaxlim , 'LineStyle' , ':' , 'LineWidth' , 1.5 , 'Color' , 'k');
            text(out.posLS , yaxlim(1) , num2str(out.posLS) , 'HorizontalAlignment' , 'Center' , 'VerticalAlignment' ,  'Top');
        end
        
        % Index plot of robust residuals
        h2=subplot(2,1,2);
        laby='Robust lts residuals';
        labx='Index number';
        resindexplot(out.residuals,'conflev',conflev,'laby',laby,'labx',labx,'numlab',out.outliers,'h',h2,'title','');
        drawnow;
        set(get(gca,'Xlabel'),'interpreter','none');
        set(get(gca,'Ylabel'),'interpreter','none');
        if ~vlt15
            set(h2,'FontSize',SizeAxesNum,'Box','on','BoxStyle','full');
        else
            set(h2,'FontSize',SizeAxesNum,'Box','on');
        end
        set(h2,'XTick',xtickval,'XTickLabel',xticklab,'XTickMode','manual');
    end
    
    if plots==2 && lsh>0
        
        % Values of the target function for each tentative level shift position
        figure;
        boxplot(ALLnumscale2(:,1:end),LSH(1:end)','labelorientation','inline');
        % boxplot uses text to put the labels on the X axes labeling, therefore
        % we have to use findobj here to fix the font size
        txt = findobj(gca,'Type','text');
        %set(gca,'XTickLabel',{' '}); % this would delete all x labels
        nx = numel(txt);
        if nx > 20
            txt2 = txt(1:mod(nx,20):nx,:);
            delete(txt(setdiff(1:nx,1:mod(nx,20):nx),:));
        end
        set(txt2,'FontSize',SizeAxesNum,'VerticalAlignment', 'Middle');
        hold('on');
        plot(numscale2LSH(:,2));
        set(gca,'Fontsize',SizeAxesNum);
        xlabel('Position of level shift','FontSize',FontSize,'interpreter','none');
        title('Target function values','interpreter','none','FontSize',FontSize+2);
        ylim([min(ALLnumscale2(:)), prctile(ALLnumscale2(:),90)]);
        
        % Level Shift local refinement
        figure;
        sb1 = subplot(2,1,1);
        plot(Likloc(:,1),Likloc(:,2));
        %xlabel('Position of level shift','FontSize',FontSize);
        ylabel('Raw residuals','FontSize',FontSize,'interpreter','none');
        set(gca,'Fontsize',SizeAxesNum);
        subplot(2,1,2);
        plot(Likloc(:,1),Likloc(:,3));
        xlabel('Position of level shift','FontSize',FontSize,'interpreter','none');
        ylabel('Huber rho residuals','FontSize',FontSize,'interpreter','none');
        set(gca,'Fontsize',SizeAxesNum);
        title(sb1,'Level Shift local refinement','interpreter','none','FontSize',FontSize+2);
        
        %     plot(LSH,NumScale2ind','o')
        %     set(gca,'FontSize',1)
        %     ylabel(['Indexes of the best ' num2str(nbestindexes) ' solutions'])
        %     xlabel('Position of level shift')
        
        % Best solutions
        figure;
        one=ones(lLSH,1);
        for j=1:nbestindexes
            text(LSH,NumScale2ind(j,:)',num2str(j*one),'FontSize',12-j*1.5);
        end
        xlim([LSH(1) LSH(end)]);
        ylim([1 bestr]);
        ylabel(['Indexes of the best ' num2str(nbestindexes) ' solutions'],'FontSize',FontSize,'interpreter','none');
        xlabel('Position of level shift','FontSize',FontSize,'interpreter','none');
        set(gca,'Fontsize',SizeAxesNum);
        hold('on');
        plot([LSH(1) LSH(end)],bestrdiv2*ones(2,1)+0.5);
        title('Best solutions','interpreter','none','FontSize',FontSize+2);
        
        % units forming best h-subset
        figure;
        plot(LSH,Weimod','ko','MarkerSize',4);
        xlabel('Position of level shift','FontSize',FontSize,'interpreter','none');
        ylabel('Index number','FontSize',FontSize,'interpreter','none');
        title('o = units forming best h-subset','interpreter','none','FontSize',FontSize+2);
        set(gca,'Ytick',10:10:T,'Fontsize',SizeAxesNum);
        
        % Frequency of inclusion inside subset
        figure;
        subplot(2,1,1);
        bar(WEIisum(:,locmin)/nselected);
        title({'Frequency of inclusion in the h subset:' , [' after ' num2str(refsteps) ' iterations']},'interpreter','none','FontSize',FontSize+2);
        %ylabel('Frequency','FontSize',FontSize);
        %xlabel('Index number','FontSize',FontSize);
        set(gca,'Xtick',1:10:T,'Fontsize',SizeAxesNum);
        
        subplot(2,1,2);
        
        bar(WEIibest10sum(:,locmin)/size(bestyhatall,2));
        title(['among the ' num2str(size(bestyhatall,2)) ' best subsets'],'interpreter','none','FontSize',FontSize+2);
        %ylabel('Frequency','FontSize',FontSize);
        xlabel('Index number','FontSize',FontSize,'interpreter','none');
        set(gca,'Xtick',1:10:T,'Fontsize',SizeAxesNum);
    end
end

% Part of the code to find the F test for the final harmonic of the seasonal
% component part
bsb=bsbModSel;

pval=[];
if seasonal>0 && seasonal<6
    % selWithoutLastHarmonic = indexes of the linear part of the model after excluding the last harmonic
    selWithoutLastHarmonic=[1:ntrend+nseaso-2 ntrend+nseaso+1:size(Xsel,2)];
    
    if varampl==0 && lshiftYN==0 % In this case the model is linear
        % Function lik constructs fitted values and residual sum of
        % squares
        betaout = Xsel(bsb,selWithoutLastHarmonic) \ yin(bsb);
        % update fitted values
        yhat = Xsel(:,selWithoutLastHarmonic) * betaout;
        
        s2reduced=sum((yin(bsb)-yhat(bsb)).^2);
        
    elseif   varampl==0 && lshiftYN==1
        % In this case there is just level shift however we do not redo
        % the non linear estimation but a simple LS
        Xselreduced= Xsel(:,selWithoutLastHarmonic);
        Xseldum=[Xselreduced  Xlshift];
        betaout = Xseldum(bsb,:) \ yin(bsb);
        
        % find fitted values using all observations
        yhat =  Xseldum * betaout;
        
        s2reduced=sum((yin(bsb)-yhat(bsb)).^2);
        
    else % model is non linear because there is time varying amplitude in seasonal component
        Xtrendf=Xtrend(bsb,:);
        
        % Remove the last harmonic from Xseaso
        seasonal=seasonal-1;
        if seasonal==0
            Xseaso=[];
            Xseasof=[];
            yhatseaso=0;
        else
            Xseaso=Xseaso(:,1:end-2);
            Xseasof=Xseaso(bsb,:);
        end
        
        if ~isempty(X)
            Xf=X(bsb,:);
        end
        Seqf=Seq(bsb,:);
        yf=yin(bsb);
        
        lasind=length(brobfinal);
        
        selWithoutLastHarmonic=[1:ntrend+nseaso-2 ntrend+nseaso+1:lasind];
        
        % If there is no seasonality it is also necessary to
        % remove the non linear part of the seasonal component
        % that is, it is necessary to select the elements of vector selWithoutLastHarmonic
        % apart from those which are in posvarampl
        if seasonal==0
            selWithoutLastHarmonic=setdiff(selWithoutLastHarmonic,posvarampl);
            varampl=0;
        end
        
        if coder.target('MATLAB')
            if lshiftYN==1
                Xlshiftf=Xlshift(bsb);
                [betaout,~,~,~,~,~]  = nlinfit(Xtrendf,yf,@likyhat,brobfinal(selWithoutLastHarmonic(1:end-1)));
            else
                [betaout,~,~,~,~,~]  = nlinfit(Xtrendf,yf,@likyhat,brobfinal(selWithoutLastHarmonic));
            end
        else
            % TODO nlinfit not supported by MATLAB C Coder
            if lshiftYN==1
                betaout= brobfinal(selWithoutLastHarmonic(1:end-1));
            else
                betaout=brobfinal(selWithoutLastHarmonic);
            end
        end
        
        
        lik(betaout);
        % Computation of residuals.
        residuals=yin(bsb)-yhat;
        s2reduced=sum(residuals.^2);
    end
    v1=2;
    numFtest=(s2reduced-s2full)/v1;
    v2=(length(bsb)-p);
    denFtest=s2full/v2;
    pval=1-fcdf(numFtest/denFtest,v1,v2);
elseif seasonal>0
    % In presence of 6 harmonics, the last one is just made up of a single
    % variable, therefore the p value is just the p value of the assocaited
    % t-stat
    pval=B(ntrend+nseaso,4);
end
out.LastHarmonicPval=pval;

if lshiftYN==1
    lsdet=FSRinvmdr([length(LSH) abs(B(end,3))],p-1);
    LevelShiftPval=1-lsdet(1,2);
    out.LevelShiftPval=LevelShiftPval;
else
    if ~coder.target('MATLAB')
        out.LevelShiftPval=[];
    end
end

if coder.target('MATLAB')
    % Restore the previous state of the warnings
    warning(warnrank.state,'MATLAB:rankDeficientMatrix');
    warning(warnsing.state,'MATLAB:singularMatrix');
    warning(warnnear.state,'MATLAB:nearlySingularMatrix');
end

% check about the y global variable
if ~isequal(yin,yin)
    error('FSDA:LTSts:yDiscrepancy','y should not change in this code. Please check if the global variable has been misused.');
end

%% The part below contains subfunctions which are used only inside this file

% ALS computes Alternating Least Squares estimate of beta starting from
% vector beta0. The rows which are used are those specified in global
% variable bsb
    function [newbeta,exitflag]=ALS(beta0)
        iter        = 0;
        betadiff    = 9999;
        newbeta=beta0;
        oldbeta=beta0;
        % exitflag = flag which informs about convergence. exitflag =0
        % implies normal convergence, else no convergence has been obtained
        exitflag=0;
        
        % Define all the relevant matrices before the loop
        Seqbsb=Seq(bsb,1);
        Xseasobsb=Xseaso(bsb,:);
        Xtrendbsb=Xtrend(bsb,:);
        yinbsb=yin(bsb);
        indnlseaso=(trend+2+nexpl):(trend+2+nexpl+varampl-1);
        Seqbsbvarampl=Seq(bsb,2:varampl+1);
        
        if isemptyX
            if lshiftYN==1
                Xlshiftbsb=Xlshift(bsb);
                XtrendXbsbXseasonXlshift=[Xtrendbsb Seqbsbvarampl Xlshiftbsb];
                XtrendbsbXbsbXlshiftbsb=[Xtrendbsb Xlshiftbsb];
                indnlseasoc=[1:trend+1 trend+2+nexpl+varampl];
            else
                XtrendXbsbXseasonXlshift=[Xtrendbsb Seqbsbvarampl];
                XtrendbsbXbsbXlshiftbsb=Xtrendbsb;
                indnlseasoc=1:(trend+1);
            end
        else
            Xbsb= X(bsb,:);
            XtrendbsbXbsb=[Xtrendbsb Xbsb];
            if lshiftYN==1
                Xlshiftbsb=Xlshift(bsb);
                XtrendXbsbXseasonXlshift=[XtrendbsbXbsb Seqbsbvarampl Xlshiftbsb];
                XtrendbsbXbsbXlshiftbsb=[Xtrendbsb Xbsb Xlshiftbsb];
                indnlseasoc=[1:trend+1+nexpl trend+2+nexpl+varampl];
            else
                XtrendXbsbXseasonXlshift=[XtrendbsbXbsb Seqbsbvarampl];
                XtrendbsbXbsbXlshiftbsb=[Xtrendbsb Xbsb];
                indnlseasoc=1:trend+1+nexpl;
            end
        end
        
        while ( (betadiff > reftolALS) && (iter < refstepsALS) )
            iter = iter + 1;
            
            % b2378 estimate of linear part of seasonal component
            b2378=newbeta(indlinsc);
            % at= yhatseaso = fitted values for linear part of seasonal
            % component
            at=Xseasobsb*b2378;
            
            % OLS to estimate coefficients of trend + expl variables + non lin coeff of
            % seasonal + coefficient of fixed level shift
            % trlshift is the matrix of explanatory variables
            XtrendXbsbXseasonXlshift(:,indnlseaso)=at.*Seqbsbvarampl;
            
            % b0145 = coefficients of intercept trend + expl var + non
            % linear part of seasonal component + level shift
            b0145=XtrendXbsbXseasonXlshift\(yinbsb-at) ;
            
            % Now find new coefficients of linear part of seasonal
            % component in the regression of y-trend-expl-lsihft versus
            % vector which contains non linear part of seasonal component
            % which multiplies each column of matrix Xseaso (linear part of
            % seasonal component)
            yhatnlseaso = Seqbsb + Seqbsbvarampl * b0145(indnlseaso);
            b2378 = (yhatnlseaso.*Xseasobsb) \ (yinbsb - XtrendbsbXbsbXlshiftbsb * b0145(indnlseasoc));
            
            % Store new value of beta
            newbeta(indlinsc)=b2378;
            newbeta(otherind)=b0145;
            
            % betadiff is linked to the tolerance (specified in scalar
            % reftol)
            betadiff = norm(oldbeta - newbeta,1) / norm(newbeta,1);
            
            oldbeta=newbeta;
            
            % exit from the loop if the new beta has singular values. In
            % such a case, any intermediate estimate is not reliable and we
            % can just keep the initialbeta and initial scale.
            if (any(isnan(newbeta)))
                newbeta = beta0;
                exitflag=-1;
                break
            end
        end
    end

    function [newbeta,exitflag]=ALSbsxfun(beta0)
        iter        = 0;
        betadiff    = 9999;
        newbeta=beta0;
        oldbeta=beta0;
        % exitflag = flag which informs about convergence. exitflag =0
        % implies normal convergence, else no convergence has been obtained
        exitflag=0;
        
        while ( (betadiff > reftolALS) && (iter < refstepsALS) )
            iter = iter + 1;
            
            % b2378 estimate of linear part of seasonal component
            b2378=newbeta(indlinsc);
            % at= yhatseaso = fitted values for linear part of seasonal
            % component
            at=Xseaso(bsb,:)*b2378;
            
            % OLS to estimate coefficients of trend + expl variables + non lin coeff of
            % seasonal + coefficient of fixed level shift
            % trlshift is the matrix of explanatory variables
            if isemptyX
                if lshiftYN==1
                    tr_expl_nls_lshift=[Xtrend(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1)) Xlshift(bsb)];
                else
                    tr_expl_nls_lshift=[Xtrend(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1))];
                end
            else
                if lshiftYN==1
                    tr_expl_nls_lshift=[Xtrend(bsb,:) X(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1)) Xlshift(bsb)];
                else
                    tr_expl_nls_lshift=[Xtrend(bsb,:) X(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1))];
                end
            end
            % b0145 = coefficients of intercept trend + expl var + non
            % linear part of seasonal component + level shift
            b0145=tr_expl_nls_lshift\(yin(bsb)-at) ;
            
            % Now find new coefficients of linear part of seasonal
            % component in the regression of y-trend-expl-lsihft versus
            % vector which contains non linear part of seasonal component
            % which multiplies each column of matrix Xseaso (linear part of
            % seasonal component)
            yhatnlseaso=Seq(bsb,1)+ Seq(bsb,2:varampl+1)*b0145((trend+2+nexpl):(trend+2+nexpl+varampl-1));
            if isemptyX
                if lshiftYN==1
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(yin(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-Xlshift(bsb)*b0145(end));
                else
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(yin(bsb)-Xtrend(bsb,:)*b0145(1:trend+1));
                end
            else
                if lshiftYN==1
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(yin(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-X(bsb,:)*b0145((trend+2):(trend+1+nexpl)) - Xlshift(bsb)*b0145(end));
                else
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(yin(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-X(bsb,:)*b0145((trend+2):(trend+1+nexpl)));
                end
            end
            
            newbeta(indlinsc)=b2378;
            
            newbeta(otherind)=b0145;
            
            % betadiff is linked to the tolerance (specified in reftol)
            betadiff = norm(oldbeta - newbeta,1) / norm(newbeta,1);
            
            oldbeta=newbeta;
            
            % exit from the loop if the new beta has singular values. In
            % such a case, any intermediate estimate is not reliable and we
            % can just keep the initialbeta and initial scale.
            if (any(isnan(newbeta)))
                newbeta = beta0;
                exitflag=-1;
                break
            end
        end
    end

% lik computes the objective function (residual sum of squares/2 = negative
% log likelihood) which must be minimized for the units specified inside
% global variable bsb. Note that given that yhat is global it is possible
% to call this function to compute fitted values for the units specified in bsb
    function obj=lik(beta0)
        
        yhattrend=Xtrend(bsb,:)*beta0(1:trend+1);
        npar=trend+1;
        
        if seasonal >0
            if seasonal<s/2
                yhatseaso=Xseaso(bsb,:)*beta0(npar+1:npar+seasonal*2);
                npar=npar+seasonal*2;
            else
                yhatseaso=Xseaso(bsb,:)*beta0(npar+1:npar+seasonal*2-1);
                npar=npar+seasonal*2-1;
            end
            
            if varampl>0
                Xtre=1+Seq(bsb,2:varampl+1)*beta0((npar+1+nexpl):(npar+varampl+nexpl));
                yhatseaso=Xtre.*yhatseaso;
                npar=npar+varampl;
            end
        end
        
        if isemptyX
            yhatX=0;
        else
            % Note the order of coefficients is trend, linear part of
            % seasonal component, expl variables, non linear part of
            % seasonal component, level shift
            yhatX=X(bsb,:)*beta0(npar+1-varampl:npar+nexpl-varampl);
            npar=npar+nexpl;
        end
        
        if lshiftYN==1
            %  \beta_(npar+1)* I(t \geq \beta_(npar+2)) where beta_(npar+1)
            %  is a real number and \beta_(npar+2) is a integer which
            %  denotes the period in which level shift shows up
            yhatlshift=beta0(npar+1)*Xlshift(bsb);
        else
            yhatlshift=0;
        end
        
        % Fitted values from trend (yhattrend), (time varying) seasonal
        % (yhatseaso), explanatory variables (yhatX) and level shift
        % component (yhatlshift)
        yhat=yhattrend+yhatseaso+yhatX+yhatlshift;
        
        %         % Additional regression due to the presence of the autoregressive
        %         % component
        %         if ARp>0
        %             Yhatlagged=zeros(length(bsb),ARp);
        %             for jj=1:ARp
        %                 Yhatlagged(:,jj)=[NaN(jj,1); yhat(1:end-jj)];
        %             end
        %             Yhatlagged=Yhatlagged(ARp+1:end,:);
        %             yinbsb=yin(bsb);
        %             blagged=Yhatlagged\yinbsb(ARp+1:end);
        %             yhat(ARp+1:end)=Yhatlagged*blagged;
        %         end
        
        % obj = sum of squares of residuals/2 = negative log likelihood
        obj=sum((yin(bsb)-yhat).^2)/2;
        % format long
        % disp(obj)
    end


% likyhat computes fitted values using vector of regression coefficients
% beta0. Note that matrices Xtrendf, Xseasof, Seqf, Xf contain n-k rows.
% This function is called in the very last step of the procedure when
% routine nlinfit is invoked. Please, note the difference beween likyhat
% and lik
    function objyhat=likyhat(beta0,Xtrendf)
        
        yhattrend=Xtrendf*beta0(1:trend+1);
        
        npar=trend+1;
        
        if seasonal >0
            if seasonal<s/2
                yhatseaso=Xseasof*beta0(npar+1:npar+seasonal*2);
                npar=npar+seasonal*2;
            else
                yhatseaso=Xseasof*beta0(npar+1:npar+seasonal*2-1);
                npar=npar+seasonal*2-1;
            end
            
            if varampl>0
                Xtre=1+Seqf(:,2:varampl+1)*beta0((npar+1+nexpl):(npar+varampl+nexpl));
                yhatseaso=Xtre.*yhatseaso;
                npar=npar+varampl;
            end
        end
        
        if isemptyX
            yhatX=0;
        else
            % Note the order of coefficients is trend, linear part of
            % seasonal component, expl variables, non linear part of
            % seasonal component, level shift
            yhatX=Xf(:,:)*beta0(npar+1-varampl:npar+nexpl-varampl);
            npar=npar+nexpl;
        end
        
        if lshiftYN==1
            %  \beta_(npar+1)* I(t \geq \beta_(npar+2)) where beta_(npar+1)
            %  is a real number and \beta_(npar+2) is a integer which
            %  denotes the period in which level shift shows up
            
            yhatlshift=beta0(npar+1)*Xlshiftf;
        else
            yhatlshift=0;
        end
        
        % objhat = fitted values from trend (yhattrend), (time varying) seasonal
        % (yhatseaso), explanatory variables (yhatX) and level shift
        % component (yhatlshift)
        objyhat=yhattrend+yhatseaso+yhatX+yhatlshift;
    end

% -------------------------------------------------------------------
% subfunction IRWLSreg
% -------------------------------------------------------------------

    function outIRWLS = IRWLSreg(y,initialbeta,refsteps,reftol,h)
        %IRWLSreg (iterative reweighted least squares) does refsteps
        %refining steps from initialbeta
        %
        %  Required input arguments:
        %
        %    y:         A vector with n elements that contains the response
        %               variable. It can be both a row or column vector.
        %  initialbeta: vector containing initial estimate of beta
        %   refsteps  : scalar, number of refining (IRLS) steps
        %   reftol    : relative convergence tolerance
        %               Default value is 1e-7
        %      h      : scalar. number of observations with smallest
        %               residuals to consider
        %
        %           GLOBAL VARIABLES REQUIRED
        %    yhat :     A vector with T elements (fitted values for all the
        %               observations)
        %  Output:
        %
        %  The output consists of a structure 'outIRWLS' containing the
        %  following fields:
        %      betarw  : p x 1 vector. Estimate of beta after refsteps
        %                refining steps
        %  numscale2rw : scalar. Sum of the smallest h squared residuals
        %                from final iteration (after refsteps refining
        %                step).It is the numerator of the estimate of the
        %                squared scale.
        %     weights  : n x 1 vector. Weights assigned to each observation
        %               In this case weights are 0,1. 1 for the units
        %               associated with the smallest h squared residuals
        %               from final iteration 0 for the other units.
        %   exitflag   : scalar which informs about convergence. exitflag =
        %               0 implies normal convergence
        
        % For performance reasons, the output structure is created only at
        % the end
        % outIRWLS = struct('betarw',[],'yhat',[],'weights',[],'exiflag',[],'numscale2rw',[]);
        
        % Residuals for the initialbeta
        res = y - yhat;
        
        % Squared residuals for all the observations
        r2 = res.^2;
        
        % Ordering of squared residuals
        [r2s , i_r2s] = sort(r2);
        
        % ininumscale2 = initial value for trimmed sum of squares of
        % residuals
        ininumscale2  = sum(r2s(1:h));
        
        
        % Initialize parameters for the refining steps loop
        exitfl      =0;   newbeta=0; numscale2=0; % MATLAC C coder initialization
        iter        = 0;
        betadiff    = 9999;
        if lshiftYN==1
            beta=initialbeta(1:end);
        else
            beta        = initialbeta;
        end
        
        
        while ( (betadiff > reftol) && (iter < refsteps) )
            iter = iter + 1;
            
            if constr==1
                % Constrained sum of the smallest squared residuals
                % Constrained in the sense that initialbeta(end) is always
                % forced to be in the h subset
                
                % Check that unit initialbeta(end) belongs to subset in each
                % concentration step
                if sum(i_r2s(1:h)==initialbeta(end))==0
                    bsb=[i_r2s(1:h-1); initialbeta(end)];
                else
                    % i_r2s= units with smallest h squared residuals
                    bsb = i_r2s(1:h);
                    % new coefficients based on units with smallest h squared
                    % residuals
                end
            elseif constr ==2
                % Check that units initialbeta(end) and initialbeta(end)-1
                % belong to subset in each concentration step
                booLS=sum(i_r2s(1:h)==initialbeta(end));
                booLSprev=sum(i_r2s(1:h)==initialbeta(end)-1);
                
                if booLS ==0 && booLSprev ==0
                    bsb=[i_r2s(1:h-2); initialbeta(end)-1; initialbeta(end) ];
                elseif booLS ==0
                    bsb=[i_r2s(1:h-1); initialbeta(end)];
                elseif  booLSprev ==0
                    bsb=[i_r2s(1:h-1); initialbeta(end)-1];
                else
                    bsb=i_r2s(1:h);
                end
            else
                bsb=i_r2s(1:h);
            end
            
            if varampl==0 && lshiftYN==0 % In this case the model is linear
                % Function lik constructs fitted values and residual sum of
                % squares
                newbeta = Xsel(bsb,:) \ y(bsb);
                % update residuals
                yhat = Xsel * newbeta;
                exitfl=0;
                
            elseif   lshiftYN==1
                
                if varampl>0
                    % No minimization is used but just ALS
                    
                    if verLess2016b==1
                        [newbeta,exitfl]=ALSbsxfun(initialbeta);
                    else
                        [newbeta,exitfl]=ALS(initialbeta);
                    end
                    
                    
                    % Construct vector of fitted values for all the
                    % observations
                    bsb=seq;
                    lik(newbeta);
                else
                    % If there is just level shift
                    % we update estimate of beta using simple LS
                    Xseld=[Xsel Xlshift];
                    % newb = new estimate of beta just using units forming
                    % subset (newb does not contain the position of level
                    % shift in the last position)
                    newb = Xseld(bsb,:)\ y(bsb);
                    % yhat = vector of fitted values for all obs
                    yhat=Xseld*newb;
                    % newbeta = new estimate of beta  just using units
                    % forming subset (newb also contains as last element
                    % the position of level shift)
                    newbeta=[newb; initialbeta(end)];
                    exitfl=0;
                end
                
                
            else % model is non linear because there is just the time varying amplitude in seasonal component
                
                % Use Alternative least squares to update beta (just using
                % the units forming subset)
                
                if verLess2016b==1
                    [newbeta,exitfl]=ALSbsxfun(beta);
                else
                    [newbeta,exitfl]=ALS(beta);
                end
                
                
                % Call lik  with bsb=seq in order to create the vector
                % of fitted values (yhat) using all the observations
                bsb=seq;
                lik(newbeta);
            end
            % disp([beta newbeta])
            
            % betadiff is linked to the tolerance (specified in scalar
            % reftol)
            betadiff = norm(beta - newbeta,1) / norm(beta,1);
            
            % exit from the loop if new beta contains nan In
            % such a case, any intermediate estimate is not reliable and we
            % can just keep the initialbeta and initial scale.
            if (any(isnan(newbeta))) || exitfl ~=0
                newbeta   = initialbeta;
                numscale2 = ininumscale2;
                break
            end
            
            % update residuals
            res = y - yhat;
            r2= res.^2;
            % Ordering of all new squared residuals
            [~ , i_r2s] = sort(r2);
            % update beta
            beta = newbeta;
            
        end
        
        % newbeta = the final estimate of beta to be stored in outIRWLS.betarw
        %outIRWLS.betarw = newbeta;
        
        % yhat = the final fitted values for all the observations using
        % final estimate of beta, to be stored in outIRWLS.yhat
        %outIRWLS.yhat=yhat;
        
        if exitfl==0
            
            if constr==1
                if sum(i_r2s(1:h)==initialbeta(end))==0
                    bsb=[i_r2s(1:h-1); initialbeta(end)];
                else
                    % i_r2s= units with smallest h squared residuals
                    bsb = i_r2s(1:h);
                    % new coefficients based on units with smallest h squared
                    % residuals
                end
            elseif constr ==2
                
                % Force both initialbeta(end) and initialbeta(end)-1 to
                % belong to the subset
                booLS=sum(i_r2s(1:h)==initialbeta(end));
                booLSprev=sum(i_r2s(1:h)==initialbeta(end)-1);
                
                if booLS ==0 && booLSprev ==0
                    bsb=[i_r2s(1:h-2); initialbeta(end)-1; initialbeta(end) ];
                elseif booLS ==0
                    bsb=[i_r2s(1:h-1); initialbeta(end)];
                elseif  booLSprev ==0
                    bsb=[i_r2s(1:h-1); initialbeta(end)-1];
                else
                    bsb=i_r2s(1:h);
                end
            else
                bsb=i_r2s(1:h);
            end
            numscale2 = sum(r2(bsb));
            
            % numscale2 = the final estimate of trimmed sum of squares of
            % residuals, to be stored in outIRWLS.numscale2rw
            %outIRWLS.numscale2rw = numscale2;
        else
            %outIRWLS.numscale2rw = numscale2;
        end
        % weights = the final estimate of the weights for each observation,
        % to be stored in outIRWLS.weights. In this case weights are 0,1. 1
        % for the units associated with the units formig subset from  final
        % iteration 0 for the other units.
        weights=zerT1;
        weights(bsb)=true;
        %outIRWLS.weights=weights;
        
        % exitfl = the exit flag to be stored in outIRWLS.exiflag
        %outIRWLS.exiflag=exitfl;
        
        % Store all output variables
        outIRWLS = struct('betarw',newbeta,'yhat',yhat,'weights',weights,'exiflag',exitfl,'numscale2rw',numscale2);
        
    end

if nargout>1
    varargout=Ccell;
end
end

%% corfactorRAW function
function rawcorfac=corfactorRAW(p,n,alpha)

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
elseif 0.875 < alpha && alpha < 1
    fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
else
    fp_alpha_n=1;
end
rawcorfac=1/fp_alpha_n;
if rawcorfac <=0 || rawcorfac>50
    rawcorfac=1;
    % if msg==true
    disp('Warning: problem in subfunction corfactorRAW')
    disp(['Correction factor for covariance matrix based on simulations found =' num2str(rawcorfac)])
    disp('Given that this value is clearly wrong we put it equal to 1 (no correction)')
    disp('This may happen when n is very small and p is large')
    % end
end
end

%% corfactorREW function
function rewcorfac=corfactorREW(p,n,alpha)

if p > 2
    coeffrewqpkwad875=[-0.544482443573914,1.25994483222292,2;-0.343791072183285,1.25159004257133,3]';
    coeffrewqpkwad500=[-1.02842572724793,1.67659883081926,2;-0.26800273450853,1.35968562893582,3]';
    y1_500=1+(coeffrewqpkwad500(1,1)*1)/p^coeffrewqpkwad500(2,1);
    y2_500=1+(coeffrewqpkwad500(1,2)*1)/p^coeffrewqpkwad500(2,2);
    y1_875=1+(coeffrewqpkwad875(1,1)*1)/p^coeffrewqpkwad875(2,1);
    y2_875=1+(coeffrewqpkwad875(1,2)*1)/p^coeffrewqpkwad875(2,2);
    y1_500=log(1-y1_500);
    y2_500=log(1-y2_500);
    y_500=[y1_500;y2_500];
    A_500=[1,log(1/(coeffrewqpkwad500(3,1)*p^2));1,log(1/(coeffrewqpkwad500(3,2)*p^2))];
    coeffic_500=A_500\y_500;
    y1_875=log(1-y1_875);
    y2_875=log(1-y2_875);
    y_875=[y1_875;y2_875];
    A_875=[1,log(1/(coeffrewqpkwad875(3,1)*p^2));1,log(1/(coeffrewqpkwad875(3,2)*p^2))];
    coeffic_875=A_875\y_875;
    fp_500_n=1-(exp(coeffic_500(1))*1)/n^coeffic_500(2);
    fp_875_n=1-(exp(coeffic_875(1))*1)/n^coeffic_875(2);
else
    if p == 2
        fp_500_n=1-(exp(3.11101712909049)*1)/n^1.91401056721863;
        fp_875_n=1-(exp(0.79473550581058)*1)/n^1.10081930350091;
    end
    if p == 1
        fp_500_n=1-(exp(1.11098143415027)*1)/n^1.5182890270453;
        fp_875_n=1-(exp(-0.66046776772861)*1)/n^0.88939595831888;
    end
end
if 0.5 <= alpha && alpha <= 0.875
    fp_alpha_n=fp_500_n+(fp_875_n-fp_500_n)/0.375*(alpha-0.5);
elseif 0.875 < alpha && alpha < 1
    fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
else
    fp_alpha_n=1;
end
rewcorfac=1/fp_alpha_n;
if rewcorfac <=0 || rewcorfac>50
    rewcorfac=1;
    %  if msg==true
    disp('Warning: problem in subfunction corfactorREW');
    disp(['Correction factor for covariance matrix based on simulations found =' num2str(rewcorfac)]);
    disp('Given that this value is clearly wrong we put it equal to 1 (no correction)');
    disp('This may happen when n is very small and p is large');
    %  end
end
end
%FScategory:REG-Regression