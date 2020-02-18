function out = regressts(y,varargin)
%regressts computes estimates of regression parameters for a time series models
%
%<a href="matlab: docsearchFS('regressts')">Link to the help function</a>
%
% Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%
% Optional input arguments:
%
%      model :  model type. Structure. A structure which specifies the model
%               which will be used. The model structure contains the following
%               fields:
%               model.s = scalar (length of seasonal period). For monthly
%                         data s=12 (default), for quartely data s=4, ...
%               model.trend = scalar (order of the trend component).
%                       trend = 1 implies linear trend with intercept (default),
%                       trend = 2 implies quadratic trend ...
%                       Admissible values for trend are, 0, 1, 2 and 3.
%               model.seasonal = scalar (integer specifying number of
%                        frequencies, i.e. harmonics, in the seasonal
%                        component. Possible values for seasonal are
%                        $1, 2, ..., [s/2]$, where $[s/2]=floor(s/2)$.
%                        For example:
%                        if seasonal =1 (default) we have:
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
%                        For example, seasonal =101 implies a seasonal
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
%               model.X  =  matrix of size T-by-nexpl containing the
%                         values of nexpl extra covariates which are likely
%                         to affect y.
%               model.posLS = positive integer which specifies to position
%                         to include the level shift component.
%                         For example if model.posLS =13 then the
%                         explanatory variable $I(t \geq 13})$ is created.
%                         If this field is not present or if it is empty,
%                         the level shift component is not included.
%               model.B  = column vector or matrix containing the initial
%                         values of parameter estimates which have to be used in the
%                         maximization procedure. If model.B is a matrix,
%                         then initial estimates are extracted from the
%                         first colum of this matrix. If this field is
%                         empty or if this field is not present, the
%                         initial values to be used in the maximization
%                         procedure are referred to the OLS parameter
%                         estimates of the linear part of the model. The
%                         parameters associated to time varying amplitude
%                         are initially set to 0.
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
%                               model.posLS='';
%
%      bsb :   units forming subset. Vector.
%                m x 1 vector.
%               The default value of bsb is 1:length(y), that is all units are
%               used to compute parameter estimates.
%               Example - 'bsb',[3 5 20:30]
%               Data Types - double
%
%smallsamplecor: small sample correction factor. Boolean. Boolean which
%               defines whether to use or not small sample correction
%               factor to inflate the scale estimate just in case option
%               bsb is used and length(bsb)<length(y).  If it is true
%               the small sample correction factor is used. The default
%               value of smallsamplecor is 1, that is the correction is
%               used. See
%               http://users.ugent.be/~svaelst/publications/corrections.pdf
%               for further details about the correction factor. The
%               default value of smallsamplecor is false.
%               Example - 'smallsamplecor',false
%               Data Types - logical
%
%     asymptcor: asymptotic consistency factor. Boolean. Boolean which
%               defines whether to use or not consistency correction
%               factor to inflate the scale estimate just in case option
%               bsb is used and length(bsb)<length(y).  If it is true
%               the asmptotic consistency is used. The default
%               value of asymptcor is false, that is the asymptotic
%               consistency factor is not used.
%               Example - 'asymptcor',false
%               Data Types - logical
%
%  plots :      Plot on the screen. Scalar. If equal to one a two panel plot 
%               appears on the scree. The top panel contains real and
%               fitted value. The bottom panel contains scaled residuals
%               with a 99.9 per cent confidence band, else (default) no
%               plot is shown.
%               Example - 'plots',1
%               Data Types - double
%               Remark: the plot which is produced is very simple. In order
%               to control a series of options in this plot and in order to
%               connect it dynamically to the other forward plots it is
%               necessary to use function mdrplot.
%
%  nocheck:     Check input arguments. Boolean.
%               If nocheck is equal to true no check is performed on
%               supplied structure model
%               Example - 'nocheck',false
%               Data Types - logical
%
%  dispresults : Display results of final fit. Boolean. If dispresults is
%               true,  labels of coefficients, estimated coefficients,
%               standard errors, tstat and p-values are shown on the
%               screen in a fully formatted way. The default value of
%               dispresults is false.
%               Example - 'dispresults',true
%               Data Types - logical
%
%   StartDate : The time of the first observation.
%               Numeric vector of length 2. Vector with two integers, which
%               specify a natural time unit and a (1-based) number of
%               samples into the time unit. For example, if model.s=12
%               (that is the data are monthly) and the first observation
%               starts in March 2016, then StartDate=[2016,3]; Similarly,
%               if models.s=4 (that is the data are quarterly) and the first
%               observation starts in the second quarter or year 2014, then
%               StartData=[2014,2]. The information in option StartDate
%               will be used to create in the output the dates inside the
%               time series object.
%                 Example - 'StartDate',[2016,3]
%                 Data Types - double
%
% Output:
%
%  out :     A structure containing the following fields
%
%             out.B =   Matrix containing estimated beta coefficients,
%                       (including the intercept when options.intercept=1)
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
%                        on the time series under study (X);
%                        4) non linear part of seasonal component, that is
%                        varying amplitude. If varying amplitude is of order
%                        k there are k coefficients (if present);
%                        5) level shift component (if present). In this case
%                        there are two coefficients, the second (which is
%                        also the last element of vector beta) is an integer
%                        which specifies the time in which level shift takes
%                        place and the first (which is also the penultime
%                        element of vector beta) is a real number which
%                        identifies the magnitude of the upward (downward)
%                        level shift;
%                       2nd col = standard errors;
%                       3rd col = t-statistics;
%                       4th col = p values.
%            out.yhat = vector of fitted values after final (NLS=non linear
%                       least squares) step.
%                       $ (\hat \eta_1, \hat \eta_2, \ldots, \hat \eta_T)'$
%               out.e = Vector T-by-1 containing the raw residuals.
%       out.residuals = Vector T-by-1 containing the scaled residuals.
%           out.scale = scale estimate of the residuals
%                     \[
%                     \hat \sigma = cor \times \sum_{i \in S_m} [y_i- \eta(x_i,\hat \beta)]^2/(m-p)
%                     \]
%                     where $S_m$ is a set of cardinality $m$ which
%                     contains the units belonging to bsb and  $p$
%                     is the total number of estimated parameters and $cor$
%                     is a correction factor to make the estimator
%                     consistent (see input options smallsamplecor and
%                     asymptcor).
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
%            out.covB = $cov(\beta)$. p-by-p, square matrix containing
%                       variance covariance matrix of estimated coefficients.
%            out.y    = response vector y.
%            out.X    = data matrix X containing trend, seasonal, expl and
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
%                       n-length(out.outliers)-by-p.
%         out.Exflag  = Reason nlinfit stopped. Integer.
%                       If the model is non linear out.Exflag is equal to 1
%                       if the maximization procedure did not produce
%                       warnings or the warning was different from
%                       "ILL Conditiioned Jacobian". For any other warning
%                       which is produced (for example,
%                       "Overparameterized", "IterationLimitExceeded",
%                       'MATLAB:rankDeficientMatrix") out.Exflag is equal
%                       to -1;
%       out.ExflagALS  = Reason ALS routine stopped. Integer.
%                       If in the iterations missing or NaN are found
%                       out.ExflagALS=-1.
%                       Note that ALS routine is used just in case there
%                       was no convergence inside routine nlinfit in order
%                       to provide a better set of initial parameter
%                       estimates before retrying nlinfit.
%                       If there was immediate convergence in nlinfit
%                       out.ExflagALS is empty.
%       out.iterALS    = Number of iterations in the ALS routine. Intger
%                       between 1 and 10000 (maximum number of iterations).
%                       Note that ALS routine is used just in case there
%                       was no convergence inside routine nlinfit in order
%                       to provide a better set of initial parameter
%                       estimates before retrying nlinfit.
%                       If there was immediate convergence in nlinfit
%                       out.iterALS is empty.
%
%
% See also: regressB, regressH, LTSts, FSRts, FSRtsmdr, forecastTS
%
% References:
%
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [RPRH]
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('regressts')">Link to the help function</a>
%
%$LastChangedDate:: 2018-06-21 15:29:09 #$: Date of the last commit

% Examples:

%{
    %% regressTS with all the default options.
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
    % linear trend + just one harmonic for seasonal
    out=regressts(y);    
%}

%{
    %% Example of the use of input option model and plots.
    % Model with linear trend, two harmonics for seasonal component and
    % varying amplitude using a linear trend.
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
    out=regressts(y,'model',model,'plots',1);    
%}

%{
    % Example of the use of input option bsb and dispresults.
    % Model with linear trend, two harmonics for seasonal component and
    % varying amplitude using a linear trend.
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
    % Fit is based just on the first 40 observations
    out=regressts(y,'model',model,'bsb',[1:40],'dispresults', true, 'plots',1);    
%}

%{
    % Example of the use of input option  StartDate.
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
    StartDate=[1949,1]
    % Imposed level shift in position t=60 and 4 harmonics
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=104;         % four harmonics with time varying seasonality
    model.posLS=60;             % level shift in position t=60
    out=regressts(y,'model',model,'plots',1,'StartDate',StartDate);    
%}

%{
    % Compare scaled residuals using or not correction factors.
    % Example of the use of input option model and plots.
    % Model with linear trend, two harmonics for seasonal component and
    % varying amplitude using a linear trend.
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
    model.seasonal=104;         % two harmonics with time varying seasonality
    bsini=1:72;
    out=regressts(y,'model',model,'bsb',bsini);    
    outCR=regressts(y,'model',model,'bsb',bsini,'smallsamplecor',true,'asymptcor',true);   
    h1=subplot(2,1,1);
    title('Estimate of the scale without correction factors')
    resindexplot(out.residuals,'h',h1)
    h2=subplot(2,1,2);
    resindexplot(outCR.residuals,'h',h2)
    title('Estimate of the scale with correction factors')
%}  

%% Beginning of code 

% Input parameters checking


% Set up values for default model
modeldef          = struct;
modeldef.trend    = 1;
modeldef.s        = 12;       % monthly time series
modeldef.seasonal = 1;
modeldef.X        = [];      % no explanatory variables
modeldef.posLS   = [];       % no level shift
modeldef.B        = [];      % empty initial parameter values
nocheck           = false;
StartDate         = '';

T=length(y);
bsbini=1:T;


%% User options

dispresultsdef=false;

options=struct('model',modeldef,'nocheck',0,'dispresults',dispresultsdef,...
    'StartDate',StartDate,'bsb',bsbini,'plots',0,...
    'smallsamplecor',false,'asymptcor',false);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:regressts:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin<1
    error('FSDA:regressts:missingInputs','response y is missing');
end

if nargin >1
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    
    % And check if the optional user parameters are reasonable.
    if ~isequal(options.model,modeldef)
        fld=fieldnames(options.model);
        
        if nocheck == false
            % Check if user options inside options.model are valid options
            chkoptions(modeldef,fld)
        end
        for i=1:length(fld)
            modeldef.(fld{i})=options.model.(fld{i});
        end
    end
    
end

% Get options chosen by the user
dispresults=options.dispresults;
model = modeldef;
bsb=options.bsb;
plots=options.plots;
smallsamplecor=options.smallsamplecor;
asymptcor=options.asymptcor;
StartDate = options.StartDate;

% Get model parameters
trend    = model.trend;       % get kind of  trend
s        = model.s;           % get periodicity of time series
seasonal = model.seasonal;    % get number of harmonics

%% Set up the model
if isfield(model,'posLS') && ~isempty(model.posLS)
    lshift   = model.posLS;
    posLS =lshift;
else
    lshift=0;
end

T=length(y);
mm=length(bsb);

% seq is the vector which will contain linear time trend
seq   = (1:T)';
one   = ones(T,1);

% Construct the matrices which are fixed in each step of the minimization
% procedure
Seq = [one seq seq.^2 seq.^3];

% Define matrix which contains linear quadratic of cubic trend
Xtrend = Seq(:,1:trend+1);

ntrend = size(Xtrend,2);

% seasonal component
if seasonal >0
    sstring=num2str(seasonal);
    if seasonal>100
        varampl=str2double(sstring(1));
        seasonal=str2double(sstring(2:3));
    else
        varampl=0;
    end
    
    if seasonal < 1 || seasonal >floor(s/2)
        error('FSDA:regressts:WrongInput',['Seasonal component must be an integer between 1 and ' num2str(floor(s/2))])
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
    Xseaso=[];
    nseaso=0;
    varampl=0;
end

X=model.X;
isemptyX=isempty(X);
if isemptyX
    % nexpl = number of potential explanatory variables
    nexpl=0;
else
    nexpl=size(X,2);
end

% Define the explanatory variable associated to the level shift component
if lshift>0
    % Xlshift = explanatory variable associated with
    % level shift Xlshift is 0 up to lsh-1 and 1 from
    % lsh to T
    Xlshift= [zeros(posLS-1,1);ones(T-posLS+1,1)];
else
    Xlshift =[];
end

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

if lshift>0
    Xsel=[Xsel Xlshift];
end

% pini = number of parameters in the linear model without level shifts nor
% varying amplitude
% ntrend = number of trend parameters,
% nseaso = number of parameters associated with the harmonics,
% nexpl = number of explanatory variables,
% 1 parameter for fixed level shift position
pini=ntrend+nseaso+nexpl+(lshift>0);

% p = total number of parameters in the model
% nini +
% varampl = number of parameters involving time varying trend,
p=pini+varampl;

% indexes of linear part of seasonal component
if seasonal <6
    indlinsc=(trend+2):(trend+1+seasonal*2);
else
    indlinsc=(trend+2):(trend+1+seasonal*2-1);
end

otherind=setdiff(1:p,indlinsc);


%% Find estimate of beta and residuals
Exflag=1;
ExflagALS=[];
iterA=[];

if varampl==0  % In this case the model is linear
    
    % compute parameter estimates is based on subset
    betaout=Xsel(bsb,:)\y(bsb);
    
    
    % Compute fitted values for all n units
    yhat = Xsel * betaout;
    
    % Compute estimate of \sigma^2 just using units belonging to subset
    s2=sum((y(bsb)-yhat(bsb)).^2)/(mm-size(Xsel,2));
    invXX=inv(Xsel'*Xsel);
    % covB = covariance matrix of parameter estimates
    covB=s2*invXX; %#ok<MINV>
    
else % model is non linear because there is time varying amplitude in seasonal component
    
    Xseasof=Xseaso(bsb,:);
    if ~isempty(X)
        Xf=X(bsb,:);
    end
    Seqf=Seq(bsb,:);
    
    
    yf=y(bsb);
    Xtrendf=Xtrend(bsb,:);
    
    if lshift>0
        Xlshiftf=Xlshift(bsb);
    end
    
    % If model contains a field named B than use the first column of this field
    % as initial parameter value, else use OLS estimate based on linear part of
    % the model
    if isfield(model,'B') && ~isempty(model.B)
        b=model.B(:,1); % get initial estimate of parameter values
    else
        
        % initial value of parameter estimates is based on subset
        b=Xsel(bsb,:)\y(bsb);
        if varampl>0
            if lshift>0
                b=[b(1:end-1); 0.01*zeros(varampl,1); b(end)];
                
            else
                b=[b; 0.01*zeros(varampl,1)];
            end
        end
    end
    
    
    % MaxIter = Maximum number of iterations in the maximization procedure
    MaxIter=1000;
    % TolX = Convergence tolerance in the maximization procedure
    TolX=1e-7;
    % DisplayLevel
    DisplayLevel='';
    nlinfitOptions=statset('Display',DisplayLevel,'MaxIter',MaxIter,'TolX',TolX);
    
    warning('off','stats:nlinfit:Overparameterized');
    warning('off','stats:nlinfit:IterationLimitExceeded');
    warning('off','stats:nlinfit:IllConditionedJacobian')
    warning('off','MATLAB:rankDeficientMatrix')
    
    
    iterALS=0;
    while iterALS < 2
        % [betaout,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(....)
        [betaout,~,~,covB,s2,~]  = nlinfit(Xtrendf,yf,@likyhat,b,'options',nlinfitOptions);
        % Note that MSE*inv(J'*J) = covB where MSE is the fourth element
        % returned by nlinfit
        invXX=covB/s2;
        
        % Capture ID of last warning message
        [~,ID] = lastwarn;
        
        % If lastwarn is not empty and ID of last warning is different from
        % ILL Conditiioned Jacobian then try ALS to see if it is possible
        % to find a better set of starting values for parameter estimates.
        if iterALS == 0 && ~isempty(lastwarn) && ~strcmp(ID,'stats:nlinfit:IllConditionedJacobian')
            lastwarn('')
            % ID='';
            % Compute a maximum of 10000 iterations using a stopping tolerance for the iterations equal to 1e-7
            [b,ExflagALS,iterA]=ALS(y,b,10000,1e-7);
            iterALS=iterALS+1;
        else
            iterALS=2;
        end
    end
    
    % Note that MSE*inv(J'*J) = covB
    [~,ID] = lastwarn;
    
    % If lastwarn is empty it means that there was full convergence and
    % Exflag=1. Exflag is -1 if there was a warning message which was
    % different from Ill Conditioned jacobian.
    if ~isempty(lastwarn) && ~strcmp(ID,'stats:nlinfit:IllConditionedJacobian')
        Exflag=-1;
    end
    
    % Compute yhat=fitted values based on all the observations
    bsb=1:T;
    yhat=lik(betaout);
    
end


factor=1;
if mm<T
    if asymptcor == true
        a=norminv(0.5*(1+mm/T));
        %factor=1/sqrt(1-(2*a.*normpdf(a))./(2*normcdf(a)-1));
        factor=1/sqrt(1-2*(T/mm)*a.*normpdf(a));
    end
    if smallsamplecor == true
        factor=factor*sqrt(corfactorRAW(1,T,mm/T));
    end
else
    factor=1;
end


e=y-yhat;  % e = vector of raw residuals for all units using b estimated using subset
% scale = estimate of the scale
scale=factor*sqrt(s2);

% residual = vector of scaled residuals
residuals =e./scale;

% If the model is non linear, compute the linearized version of matrix X,
% that is the

if varampl>0
    fdiffstep=[];
    Xsel = getjacobianFS(betaout,fdiffstep,@lik,yhat);
end

% Store beta standard error, t stat and p values
sebetaout=sqrt(diag(covB));
tout=betaout./sebetaout;
dfe=T-length(betaout);
pval=2*(tcdf(-abs(tout), dfe));
B=[betaout sebetaout tout pval];

out=struct;
out.Exflag=Exflag;
out.ExflagALS=ExflagALS;
out.iterALS=iterA;
out.y=y;
out.yhat=yhat;
out.X=Xsel;
out.B=B;
out.covB=covB;
out.invXX=invXX;
out.scale=scale;
out.e=e;
out.residuals=residuals;

if dispresults
    
    b_trend={'b_trend1'; 'b_trend2'; 'b_trend3'; 'b_trend4'};
    b_seaso={'b_cos1'; 'b_sin1'; 'b_cos2'; 'b_sin2'; 'b_cos3'; 'b_sin3'; ...
        'b_cos4'; 'b_sin4'; 'b_cos5'; 'b_sin5'; 'b_cos6'};
    b_expl={'b_X1'; 'b_X2'; 'b_X3'; 'b_X4'; 'b_X5'; 'b_X6'};
    b_varampl={'b_varampl'; 'b_varamp2'; 'b_varamp3'};
    b_lshift={'b_lshift'; 't_lshift'};
    
    if seasonal>0
        if 2*seasonal==s
            lab=[b_trend(1:trend+1); b_seaso];
        else
            lab=[b_trend(1:trend+1); b_seaso(1:2*seasonal)];
        end
    else
        lab=b_trend(1:trend+1);
    end
    
    if nexpl>0
        lab=[lab;b_expl(1:nexpl)];
    end
    if varampl>0
        lab=[lab;b_varampl(1:varampl)];
    end
    if lshift>0
        lab=[lab; b_lshift(1)];
    end
    
    if verLessThan ('matlab','8.2.0')
        disp('           Coeff.     SE         t-stat       p-values');
        disp( [char(lab) num2str([betaout sebetaout tout pval])]);
    else
        disp([table(lab) table(betaout) table(sebetaout) table(tout) table(pval)]);
    end
    if lshift>0
        disp(['Level shift position t=' num2str(posLS)])
    end
end

%% Create plots
if plots==1
    % label the x axis using appropriate dates
    if ~isempty(StartDate)
        IniYear=StartDate(1);
        FinalYear=IniYear+ceil(T/s);
        [years, months] = meshgrid(IniYear:FinalYear, 1:12/s:12);
        years=years(1:T);
        months=months(1:T);
        % Convert date and time to serial date number
        datesnumeric=datenum(years(:), months(:), 1);
    else
        datesnumeric=(1:T)';
    end
    
    % Plot of y and yhat in top panel and scaled residuals in the bottom panel
    % Time series + fitted values
    figure
    subplot(2,1,1)
    % Plot original time series
    plot(datesnumeric(1:T),y,'k')
    hold('on')
    % plot fitted values
    plot(datesnumeric(1:T),yhat,'b-')
    
    
    title('Fit','interpreter','LaTex','FontSize',14)
    if ~isempty(StartDate)
        datetick('x','mmm-yy');
        if ~verLessThanFS(8.4)
            set(gca,'XTickLabelRotation',90);
        end
    end
    xlabel('Time')
    ylabel('Real and fitted values')
    
    % Index plot of robust residuals
    h2=subplot(2,1,2);
    laby='Scaled residuals';
    conflev=0.999;
    resindexplot(residuals,'conflev',conflev,'laby',laby,'h',h2,'title','');
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
        
        if lshift >0
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

    function yhat=lik(beta0)
        
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
        
        if lshift >0
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
        
    end




% ALS computes Alternating Least Squares estimate of beta starting from
% vector beta0. The rows which are used are those specified in global
% variable bsb
    function [newbeta,exitflag,iter]=ALS(y,beta0,maxiterALS,maxtolALS)
        iter        = 0;
        betadiff    = 9999;
        newbeta=beta0;
        oldbeta=beta0;
        % exitflag = flag which informs about convergence. exitflag =0
        % implies normal convergence, else no convergence has been obtained
        exitflag=0;
        
        while ( (betadiff > maxtolALS) && (iter < maxiterALS) )
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
                if lshift>0
                    tr_expl_nls_lshift=[Xtrend(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1)) Xlshift(bsb)];
                else
                    tr_expl_nls_lshift=[Xtrend(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1))];
                end
            else
                if lshift>0
                    tr_expl_nls_lshift=[Xtrend(bsb,:) X(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1)) Xlshift(bsb)];
                else
                    tr_expl_nls_lshift=[Xtrend(bsb,:) X(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1))];
                end
            end
            % b0145 = coefficients of intercept trend + expl var + non
            % linear part of seasonal component + level shift
            b0145=tr_expl_nls_lshift\(y(bsb)-at) ;
            
            % Now find new coefficients of linear part of seasonal
            % component in the regression of y-trend-expl-lsihft versus
            % vector which contains non linear part of seasonal component
            % which multiplies each column of matrix Xseaso (linear part of
            % seasonal component)
            yhatnlseaso=Seq(bsb,1)+ Seq(bsb,2:varampl+1)*b0145((trend+2+nexpl):(trend+2+nexpl+varampl-1));
            if isemptyX
                if lshift>0
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(y(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-Xlshift(bsb)*b0145(end));
                else
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(y(bsb)-Xtrend(bsb,:)*b0145(1:trend+1));
                end
            else
                if lshift>0
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(y(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-X(bsb,:)*b0145((trend+2):(trend+1+nexpl)) - Xlshift(bsb)*b0145(end));
                else
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(y(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-X(bsb,:)*b0145((trend+2):(trend+1+nexpl)));
                end
            end
            
            
            
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



    function J = getjacobianFS(beta,fdiffstep,modelFS,yfit)
        function yplus = call_model_nested(betaNew)
            yplus = modelFS(betaNew);
        end
        J = statjacobianFS(@call_model_nested, beta, fdiffstep, yfit);
    end % function getjacobian

    function J = statjacobianFS(func, theta, DerivStep, y0)
        %STATJACOBIAN Estimate the Jacobian of a function
        
        % J is a matrix with one row per observation and one column per model
        % parameter. J(i,j) is an estimate of the derivative of the i'th
        % observation with respect to the j'th parameter.
        
        % For performance reasons, very little error checking is done on the input
        % arguments. This function makes the following assumptions about inputs:
        %
        % * func is the model function and is a valid function handle that accepts
        %   a single input argument of the same size as theta.
        % * theta is vector or matrix of parameter values. If a matrix, each row
        %   represents a different group or observation (see "Grouping Note" below)
        %   and each column represents a different model parameter.
        % * DerivStep (optional) controls the finite differencing step size. It may
        %   be empty, scalar, or a vector of positive numbers with the number of
        %   elements equal to the number model parameters.
        % * y0 (optional) is the model function evaluated at theta. A value of []
        %   is equivalent to omitting the argument and results in the model being
        %   evaluated one additional time.
        %
        % Example 1: NLINFIT
        %   NLINFIT is used to estimate the parameters b(1) and b(2) for the model
        %   @(b,T) b(1)*sin(b(2)*T), given data at T=1:5. NLINFIT needs the
        %   Jacobian of the model function with respect to b(1) and b(2) at each T.
        %   To do this, it constructs a new function handle that is only a function
        %   of b and that "burns-in" the value of T (e.g. model2 = @(b) model1(b,T)).
        %   It then calls STATJACOBIAN with the new function handle to obtain a
        %   matrix J, where J(i,j) is an estimate of the derivative of the model
        %   with respect to the j'th parameter evaluated at T(i) and b.
        %
        % Example 2: NLMEFIT or NLMEFITSA with group-specific parameters
        %   NLMEFIT requires the Jacobian of the model function with respect to two
        %   parameters evaluated at group-specific values. (Group-specific
        %   parameters can arise, for example, from using the default FEConstDesign
        %   and REConstDesign options.) NLMEFIT calls STATJACOBIAN passing in a
        %   matrix of parameter values theta, with one row per group, where
        %   theta(i,j) represents a parameter value for i'th group and j'th
        %   parameter. STATJACOBIAN returns a matrix J, where J(i,j) is an estimate
        %   of the derivative of the model with respect to the j'th parameter,
        %   evaluated for observation i with parameter values theta(rowIdx(i),:),
        %   which are the parameter values for the observation's group.
        %
        % Example 3: NLMEFIT with observation-specific parameters
        %   NLMEFIT requires the Jacobian of the model function with respect to two
        %   parameters evaluated at observation-specific values. (Observation-
        %   specific parameters can arise, for example, from using the FEObsDesign
        %   or REObsDesign options.) NLMEFIT calls STATJACOBIAN passing in a matrix
        %   of parameter values theta, with one row per observation, where
        %   theta(i,j) represents a parameter value for the i'th observation and
        %   j'th parameter. In this case, rowIdx is 1:N, where N is the number of
        %   observations. STATJACOBIAN returns a matrix J, where J(i,j) is an
        %   estimate of the derivative of the model with respect to the j'th
        %   parameter, evaluated for observation i with parameter values
        %   theta(i,:), which are the parameter values for the observation.
        
        % Use the appropriate class for variables.
        classname = class(theta);
        
        % Handle optional arguments, starting with y0 since it will be needed to
        % determine the appropriate size for a default groups.
        if nargin < 4 || isempty(y0)
            y0 = func(theta);
        end
        
        % When there is only one group, ensure that theta is a row vector so
        % that vectoriation works properly. Also ensure that the underlying
        % function is called with an input with the original size of theta.
        thetaOriginalSize = size(theta);
        theta = reshape(theta, 1, []);
        
        func = @(theta) func(reshape(theta, thetaOriginalSize));
        
        % All observations belong to a single group; scalar expansion allows us
        % to vectorize using a scalar index.
        rowIdx = 1;
        
        [numThetaRows, numParams] = size(theta);
        
        if nargin < 3 || isempty(DerivStep)
            % Best practice for forward/backward differences:
            DerivStep = repmat(sqrt(eps(classname)), 1, numParams);
            % However, NLINFIT's default is eps^(1/3).
        elseif isscalar(DerivStep)
            DerivStep = repmat(DerivStep, 1, numParams);
        end
        
        delta = zeros(numThetaRows, numParams, classname);
        J = zeros(numel(y0), numParams, classname);
        for ii = 1:numParams
            % Calculate delta(:,ii), but remember to set it back to 0 at the end of the loop.
            delta(:,ii) = DerivStep(ii) * theta(:,ii);
            deltaZero = delta(:,ii) == 0;
            if any(deltaZero)
                % Use the norm as the "scale", or 1 if the norm is 0.
                nTheta = sqrt(sum(theta(deltaZero,:).^2, 2));
                delta(deltaZero,ii) = DerivStep(ii) * (nTheta + (nTheta==0));
            end
            thetaNew = theta + delta;
            yplus = func(thetaNew);
            dy = yplus(:) - y0(:);
            J(:,ii) = dy./delta(rowIdx,ii);
            delta(:,ii) = 0;
        end
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
end
if 0.875 < alpha && alpha < 1
    fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
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
%FScategory:REG-Regression