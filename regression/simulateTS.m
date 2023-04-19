function [out] = simulateTS(T,varargin)
%simulateTS simulates a time series with trend, time varying seasonal, level shift and irregular component
%
%<a href="matlab: docsearchFS('simulateTS')">Link to the help function</a>
%
% simulateTS simulates a time series with trend (up to third order),
% seasonality (constant or of varying amplitude) with a different number of
% harmonics and a level shift. Moreover, it is possible to add to the
% series the effect of explanatory variables.
%
%  Required input arguments:
%
%         T  :  time series length. Scalar. T is a positive integer
%               which defines the length of the simulated time series.
%
%  Optional input arguments:
%
%      model :  model type. Structure. A structure which specifies the model
%               used to simulate the time series. The structure contains
%               the following fields:
%               model.trend = scalar (order of the trend component).
%                       trend = 1 implies linear trend with intercept,
%                       trend = 2 implies quadratic trend, etc.
%                       If this field is empty the simulated time series
%                       will not contain a trend. The default value
%                       of model.trend is 1.
%               model.trendb = vector of doubles containining the beta
%                       coefficients for the trend. For example model.trend=1
%                       and model.trendb=[3.2 2] generate a linear trend
%                       of the kind 3.2+2*t.
%                       If this field is an empty double the simulated time
%                       series will not contain a trend. The default value
%                       of model.trendb is [0 1] that is a slope equal to 1
%                       and intercept equal to 0.
%               model.s = scalar greater than zero which specifies the
%                       length of the seasonal period. For monthly
%                       data (default) s=12, for quartely data s=4, ...
%                       The default value of model.s is 12, which is for
%                       monthly data.
%               model.seasonal = scalar. Integer specifying the number of
%                        frequencies, i.e. harmonics, in the seasonal
%                        component. Possible values are $1, 2, ..., [s/2]$,
%                        where $[s/2]=floor(s/2)$.
%                        For example:
%                        if seasonal = 1 (default) we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 \sin (2 \pi t/s)$;
%                        if seasonal = 2 we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 \sin (2 \pi t/s)
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
%                        For example, seasonal = 201 implies a seasonal
%                        component which just uses one frequency
%                        which grows in a quadratic way over time as
%                        follows:
%                        $(1+\beta_3 t + \beta_4  t^2)\times( \beta_1 \cos(
%                        2 \pi t/s) + \beta_2 \sin ( 2 \pi t/s))$.
%                       If this field is an empty double (default) the
%                       simulated time series will not contain a seasonal
%                       component.
%               model.seasonalb = vector of doubles containing the beta
%                       coefficients for the seasonal component.
%                       For example model.seasonal = 201 and model
%                       model.trendb = [1.2 2.3 3.4 4.5] generates a
%                       seasonal component of the kind:
%                        $(1+ 3.4 t + 4.5  t^2)\times( 1.2 \cos(
%                        2 \pi t/s) + 2.3 \sin ( 2 \pi t/s))$.
%                       If this field is an empty double (default) the
%                       simulated time series will not contain a seasonal
%                       component.
%               model.X  = scalar or matrix of size T-by-nexpl. If model.X
%                       is a matrix of size T-by-nexpl, it contains the
%                       values of nexpl extra covariates which
%                       affect y. If model.X is a scalar equal to k,
%                       where k=1, 2, ... k explanatory variables using
%                       random numbers from the normal distribution are
%                       generated. If this field is an empty double
%                       (default) the simulated time series will not
%                       contain explanatory variables.
%               model.Xb = vector of doubles containing the beta
%                       coefficients for the explanatory variables.
%                       For example model.X = 2 and model.Xb = [4,5]
%                       generate two additional explanatory variables
%                       of the kind: $ 4*randn(T,1) + 5*randn(T,1) $.
%                       If this field is an empty double (default) the
%                       simulated time series will not contain explanatory
%                       variables.
%               model.ARIMAX = boolean variable (true/false) that specifies
%                       the model for the AR component (when present). If
%                       true, then the autoregressive component involves
%                       the dependent variable, as in expression (1) of
%                       https://it.mathworks.com/help/econ/convert-between-armax-and-regarima-models.html#btw30mu
%                       with N(L)=1. If false, then the autoregressive
%                       component involves the residual term, as in
%                       expression (2) of the same page with A(L)=1. If
%                       this field is empty or not present, the default
%                       value is false.
%               model.ARp = vector with non negative integer numbers
%                       specifying the autoregressive
%                       components. For example:
%                        model.ARp=[1 2] means a AR(2) process;
%                        model.ARp=2 means just the lag 2 component;
%                        model.ARp=[1 2 5 8] means AR(2) + lag 5 + lag 8;
%                        model.ARp=0 (default) means no autoregressive component.
%               model.ARb = vector of doubles containing the beta
%                       coefficients for the autoregressive component. For
%                       example model.ARb = [0.5 -0.2] in combination with
%                       model.ARp=[1 2] and ARIMAX=true generates an AR(2)
%                       time series of the kind: $y_t = 0.5 y_{t-1} - 0.2
%                       y_{t-2}$ + seasonal + lshift + $\epsilon_t$.
%                       model.ARb must have the same number of elements of
%                       model.ARp. If this field is an empty double
%                       (default) the simulated time series will not have
%                       an autoregressive component.
%               model.lshift = scalar greater than 0 which specifies the
%                       position where to include a level shift component.
%                       If this field is an empty double (default) the
%                       simulated time series will not contain a level shift.
%               model.lshiftb = scalar double which specifies the magnitude
%                       of the level shift component.
%                       For example model.lshift = 26 and model.lshiftb = 3
%                       generates the following explanatory variable
%                        $ [zeros(25,1) + 3*ones(T-25+1,1)] $.
%                       If this field is an empty double (default) the
%                       simulated time series will not contain a level
%                       shift.
%               model.nsim = scalar wich defines the number of time series
%                       to simulate. If this field is empty or not present,
%                       the default value of 1 is used.
%               model.residuals = T x nsim matrix of doubles of length T
%                       containing the values of the residuals to use in the
%                       simulation. Default value is [].
%                model.sigmaeps = scalar wich defines the standard error of
%                       the residuals. If model.residuals is not empty, the
%                       value of model.sigmaeps will be ignored. Default
%                       value is [].
%                model.signal2noiseratio = scalar wich defines the ratio
%                       between the variance of the systematic part of the
%                       model (signal) and the variance of the noise
%                       (irregular model). The greater is this value, the
%                       smaller is the effect of the irregular component.
%                       If one between model.residuals and model.sigmaeps
%                       is not empty, the value of model.signal2noiseratio
%                       will be ignored. If both model.residuals and
%                       model.sigmaeps are empty and model.ARIMAX is true,
%                       then an algorithm will search the best value of the
%                       standard error of the residuals that guarantees the
%                       desired model.signal2noiseratio. If this field is
%                       empty or not present the default value of 1 is used.
%                 Example - 'model', model
%                 Data Types - struct
%               Remark: the default model is for monthly data with a linear
%               trend with slope 1 and intercept 0, no seasonal, no
%               level shift and a signal to noise ratio equal to 1, that is
%                               model=struct;
%                               model.s=[];
%                               model.trend=1;
%                               model.trendb=[0 1];
%                               model.X=[];
%                               model.lshift=[];
%                               model.signal2noiseratio=1;
%
%       plots : Plots on the screen. Scalar.
%               If plots == 1 a six panel plot appears on the screen.
%               Top left panel contains the simulated time series y.
%               y=TR+SE+X+LS+I of the first simulated time series.
%               Top central panel contains the signal component (that is
%               trend + seasonal + explanatory variables + level shift = TR
%               + SE + LS + X).
%               Top right panel contains the trend component (TR).
%               Bottom left panel contains the (time varying) seasonal
%               component (SE).
%               Bottom central panel contains the level shift component
%               (LS).
%               Bottom right panel contains the explanatory variable component
%               (X) if it is present, otherwise, it contains the irregular
%               (I) component of the first simulation.
%               The default value of plot is 0, that is no plot is shown on
%               the screen.
%                 Example - 'plots',1
%                 Data Types - double
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
%      nocheck: Check input arguments. Boolean. If nocheck is true no
%               check is performed on supplied input. Otherwise (default)
%               every input of the structure model is checked.
%               Example - 'nocheck',false
%               Data Types - double
%
% FileNameOutput : save simulated time series to txt file. Character.
%               If FileNameOutput is empty (default) nothing is saved on
%               the disk, else FileNameOutput will contain the path where
%               to save the file on the disk.
%               Example - 'FileNameOutput',['C:' filesep 'myoutput' fielsep 'savesimdata.txt']
%               Data Types - Character
%
% samescale :   same ylim  in the output plot.  Logical.
%               If true (default), all underlying components of the time
%               series are shown in the plot with the same scale.
%               Example - 'samescale',false
%               Data Types - logical
%
%  Output:
%
%         out:   structure which contains the following fields:
%
%                out.y = the simulated time series.
%                   Matrix of dimension T x nsim, which is sum of trend +
%                   (time varying) seasonal + explanatory variables + level
%                   shift + irregular = TR+SE+X+LS+I.
%                out.signal = signal (TR+SE+X+LS).
%                   Column vector of length T, which is sum of trend +
%                   (time varying) seasonal + explanatory variables + level
%                   shift. Signal = out.y - out.irregular.
%                out.trend = trend (TR).
%                   Column vector of length T which contains the trend
%                   component.
%                out.seasonal = (time varying) seasonal (SE).
%                   Column vector of length T which contains the seasonal
%                   component. If there is no seasonal component
%                   outyhatseaso=0.
%                out.X = explanatory variables (X).
%                   Column vector of length T which contains the component
%                   associated to the explanatory variables.
%                   If there is no explanatory variable, out.X=0.
%                out.lshift = level shift (LS).
%                   Column vector of length T which contains the level
%                   shift component.
%                   If there is no level shift component out.lshift=0.
%                out.irregular = irregular component (I).
%                   Matrix of dimension T x nsim which contains the irregular
%                   component.
%                   When the signal to noise ratio tends to infinity the
%                   irregular component tends to 0.
%                out.model = structure. The model used to simulate the time
%                   series.
%
% See also LTSts, wedgeplot
%
% References:
%
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [RPRH]
%
%
% Copyright 2008-2023.
% Written by Marco Riani, Domenico Perrotta, Peter Rousseeuw and Mia Hubert
%
%
%<a href="matlab: docsearchFS('simulateTS')">Link to the help function</a>
%
%$LastChangedDate:: 2018-02-19 17:38:15 #$: Date of the last commit

% Examples:

%{
    %% Simulated time series with linear trend.
    % A time series of 100 observations is simulated from a model which
    % contains a linear trend (with slope 1 and intercept 0), no seasonal
    % component, no explanatory variables and a signal to noise ratio egual
    % to 1 (the default).
    out=simulateTS(100,'plots',1);
%}

%{
    % Same as above, but without homogenizing the y-scale.
    close;
    out=simulateTS(100,'plots',1,'samescale',false);
%}


%{
    %% Simulated time series with a linear time varying seasonal component.
    % A time series of 100 observations is simulated from a model which
    % contains no trend, a linear time varying seasonal component with
    % three harmonics, no explanatory variables and a signal to noise ratio
    % egual to 20
    rng('default')
    rng(1)
    model=struct;
    model.trend=[];
    model.trendb=[];
    model.seasonal=103;
    model.seasonalb=40*[0.1 -0.5 0.2 -0.3 0.3 -0.1 0.222];
    model.signal2noiseratio=20;
    T=100;
    out=simulateTS(T,'model',model,'plots',1);
%}

%{
    %% Simulated time series with a quadratic time varying seasonal component.
    % A time series of 100 observations is simulated from a model which
    % contains no trend, a quadratic time varying seasonal component with
    % one harmonic, no explanatory variables and a signal to noise ratio
    % egual to 20
    rng(1)
    model=struct;
    model.trend=[];
    model.trendb=[];
    model.seasonal=201;
    model.seasonalb=40*[0.1 -0.5 10.222 -10];
    model.signal2noiseratio=20;
    T=100;
    out=simulateTS(T,'model',model,'plots',1);
%}

%{
    %% Simulated time series with quadratic trend, fixed seasonal and level shift.
    % A time series of 100 observations is simulated from a model which
    % contains a quadratic trend, a seasonal component with two harmonics
    % no explanatory variables and a level shift in position 30 with size
    % 5000 and a signal to noise ratio egual to 20
    rng(1)
    model=struct;
    model.trend=2;
    model.trendb=[5,10,-3];
    model.seasonal=2;
    model.seasonalb=100*[2 4 0.1 8];
    model.signal2noiseratio=20;
    model.lshift=30;
    model.lshiftb=5000;
    T=100;
    out=simulateTS(T,'model',model,'plots',1);
%}

%{
    %% Simulated time series with quadratic trend, fixed seasonal and LS.
    % A time series of 100 observations is simulated from a model
    % which contains a quadratic trend, a linear time varying seasonal
    % component with two harmonics no explanatory variables and a level
    % shift in position 30 with size -10000 and a signal to noise ratio
    % egual to 20
    rng(1)
    model=struct;
    model.trend=2;
    model.trendb=[5,10,-3];
    model.seasonal=102;
    model.seasonalb=100*[2 4 0.1 8 0.001];
    model.signal2noiseratio=20;
    model.lshift=30;
    model.lshiftb=-10000;
    T=100;
    out=simulateTS(T,'model',model,'plots',1);
%}

%{
    %% Simulated time series with quadratic trend, fixed seasonal, LS and
    % two explanatory variables.
    % A time series of 100 observations is simulated from a model
    % which contains a quadratic trend, a linear time varying seasonal
    % component with two harmonics, two explanatory variables and a level
    % shift in position 30 with size -40000 and a signal to noise ratio
    % egual to 10
    rng(1)
    model=struct;
    model.trend=2;
    model.trendb=[5,10,-3];
    model.seasonal=102;
    model.seasonalb=100*[2 4 0.1 8 0.001];
    model.signal2noiseratio=10;
    model.lshift=30;
    model.lshiftb=-40000;
    model.X=2;
    model.Xb=[10000 20000];
    T=100;
    out=simulateTS(T,'model',model,'plots',1);
%}

%{
    % Example of the use of option FileNameOutput.
    % In this example the simulated time series is saved into a file named
    % ysimout.txt in the current folder
    FileNameOutput=[pwd filesep 'ysimout.txt'];
    T=100;
    out=simulateTS(T,'FileNameOutput',FileNameOutput);
%}

%{
    %% Example of the use of option StartDate.
    % Suppose that the inital observation refers to February 2016.
    StartDate=[2016 2];
    % The x axis of the plots contains the dates using format mmm-yyyy
    rng(1)
    model=struct;
    model.trend=2;
    model.trendb=[5,10,-3];
    model.seasonal=102;
    model.seasonalb=100*[2 4 0.1 8 0.001];
    model.signal2noiseratio=10;
    model.lshift=30;
    model.lshiftb=-40000;
    model.X=2;
    model.Xb=[10000 20000];
    T=100;
    out=simulateTS(T,'model',model,'plots',1,'StartDate',StartDate);
%}

%{
    % Example of the use of option samescale.
    % Use a different scale for each panel in the output plot.
    rng(1)
    model=struct;
    model.trend=2;
    model.trendb=[5,10,-3];
    model.seasonal=102;
    model.seasonalb=100*[2 4 0.1 8 0.001];
    model.signal2noiseratio=10;
    model.lshift=30;
    model.lshiftb=-40000;
    model.X=2;
    model.Xb=[10000 20000];
    T=100;
    out=simulateTS(T,'model',model,'plots',1,'samescale',false);
%}

%{
    % Simulated data with linear trend and errors with AR(2) component.
    % No seasonal component.
    rng('default')
    rng(100)
    model=struct;
    model.trend=1;
    model.trendb=[5,1000];
    model.seasonal='';
    model.signal2noiseratio=10;
    model.ARp=[1 2];
    model.ARb=[0.2 0.7];
    T=100;
    out=simulateTS(T,'model',model,'plots',1);
    y=out.y;

    % The lines below just work if the econometric toolbox is present
    % The autocorrelation function of y shows just two peaks in
    % correspondence of the first two lags
    if exist('parcorr','file')>0
        figure
        parcorr(out.y)
    end
%}

%{
    % Example of simulation and fitting.
    % Simulated data with linear trend, errors with AR(2) component and 1
    % explanatory variable.
    rng(100)
    model=struct;
    model.trend=1;
    model.trendb=[5,1000];
    model.seasonal='';
    model.signal2noiseratio=10;
    model.ARp=[1 2];
    model.ARb=[0.2 0.7];
    T=100;
    X=1e+2*randn(T,1);
    model.X=X;
    model.Xb=100;
    out=simulateTS(T,'model',model,'plots',1);
    y=out.y;
    % Fit a model with linear trend, AR(3) and the true expl. variable
    model=struct;
    model.trend=1;
    model.seasonal=0;
    model.lshift=0;
    model.X=X;
    model.ARp=[1 2 3];
    out=LTSts(y,'model',model,'plots',1,'dispresults',true);
%}

%% Beginning of code

% seq is the vector which will contain linear time trend
seq   = (1:T)';
one   = ones(T,1);

if nargin<1
    error('FSDA:simulateTS:MissingInputs','Input time series is missing');
end

% Set up values for default model
modeldef          = struct;
modeldef.trend    = 1;        % linear trend with intercept
modeldef.trendb   = [0;1];    % coefficients for the trend
modeldef.s        = 12;       % monthly time series
modeldef.seasonal = [];       % no seasonal component
modeldef.seasonalb= [];       %
modeldef.X        = [];       % no explanatory variables
modeldef.Xb       = [];       %
modeldef.lshift   = [];       % no level shift
modeldef.lshiftb  = [];       %
modeldef.ARIMAX   = false;    % autoregressive component on the residuals
modeldef.ARp      = [];       % no autoregressive component
modeldef.ARb      = [];       %
modeldef.nsim     = 1;        % simulate only 1 time series
modeldef.residuals= [];       % T x nsim matrix of residuals
modeldef.sigmaeps = [];       % std err of residuals
modeldef.signal2noiseratio=[]; % ratio between variances of signal and irregular parts
nocheck           = false;
plots             = 0;
FileNameOutput    = '';
StartDate         = '';

options=struct('model',modeldef,...
    'nocheck',nocheck,'plots',plots,...
    'FileNameOutput',FileNameOutput,...
    'StartDate',StartDate,'samescale',true);


%% User options

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)

    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:simulateTS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end

    % Check if all the specified optional arguments were present in
    % structure options Remark: the nocheck option has already been dealt
    % by routine chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:simulateTS:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    nocheck         = options.nocheck;
    plots           = options.plots;
    FileNameOutput  = options.FileNameOutput;
    StartDate       = options.StartDate;
    samescale       = options.samescale;
end

% Default values for the optional parameters are set in structure 'options'
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

%% Get model parameters

model    = modeldef;
trend    = model.trend;       % get kind of  trend
trendb   = model.trendb;      % coefficients for trend
s        = model.s;           % get periodicity of time series
seasonal = model.seasonal;    % get number of harmonics
seasonalb= model.seasonalb;   % coefficients for seasonal component
X        = model.X;
Xb       = model.Xb;
lshift   = model.lshift;      % get level shift
lshiftb  = model.lshiftb;     % get coefficient of level shift
ARIMAX   = model.ARIMAX;      % model for AR component
ARp      = model.ARp;         % lags of the autoregressive component
ARb      = model.ARb;         % coefficients of the autoregressive component
residuals= model.residuals;   % T x nsim matrix of residuals
sigmaeps = model.sigmaeps;    % std err of the residuals
signal2noiseratio=model.signal2noiseratio;
nsim     = model.nsim;

Seq = [one seq seq.^2 seq.^3];

if nocheck == false
    % Check if the optional user parameters are valid.
    if s <=0
        error('FSDA:simulateTS:WrongInput',['s=' num2str(s) 'is the periodicity of the time series (cannot be negative)'])
    end
end

%% checks on the trend component
if ~isempty(trend)
    if nocheck == false

        if isempty(intersect(trend,0:3))
            error('FSDA:LTSts:WrongInput','Trend must assume the following values: 0  1 or 2 or 3')
        end

        if isempty(trendb) && ~isempty(trend)
            disp('Warning: option trend has not been specified but the beta coefficients for the trend have been specified')
            error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
        end
    end
    if length(trendb)~=trend+1
        disp('Warning: trend order does not match with the length of the beta coefficients supplied for trend')
        disp('For example: trend order equal to 0 must have just one beta coefficient.')
        disp('For example: trend order equal to 1 must have two beta coefficients.')
        disp('...')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    else
        Xtrend = Seq(:,1:trend+1);
    end
else
    % Controls on the trend component
    if ~isempty(trendb)
        disp('Warning: option trend has been specified but the beta coefficients for the trend have not been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    Xtrend=Seq(:,1);
    trendb=0;
end

%% checks on the explanatory variables
% Define nexpl = number of explanatory variables
isemptyX=isempty(X);
if isemptyX
    % nexpl = number of potential explanatory variables
    nexpl=0;
elseif isscalar(X)
    nexpl=X;
    X=randn(T,nexpl);
else
    nexpl=size(X,2);
end

if nocheck == false
    % Controls on the explanatory variables component
    if isempty(X) && ~isempty(Xb)
        disp('Warning: option X has not been specified but the beta coefficients for the explanatory variables have been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    if isempty(Xb) && ~isempty(X)
        disp('Warning: option X has been specified but the beta coefficients for the explanatory variables have not been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    if ~isempty(X) && length(Xb)~=nexpl
        % Define matrix which contains explanatory variables
        disp(['Warning: option X has been specified but the length of ' ...
            'beta coefficients for X is not in agreeement'])
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    % Controls on level shift component
    if isempty(lshift) && ~isempty(lshiftb)
        disp('Warning: option lshift has not been specified but the beta coefficient for the level shift has been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    if isempty(lshiftb) && ~isempty(lshift)
        disp('Warning: option lshift has been specified but the beta coefficient for the level shift has not been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    if ~isempty(lshift) && lshift<=1 && lshift> T
        error('FSDA:simulateTS:WrongInput','Level shift must be an interger in the range 1 T')
    end

end


%% checks on the seasonal component
if seasonal==0; seasonal=[]; end
if ~isempty(seasonal) && seasonal >0
    sstring=num2str(seasonal);
    if seasonal>100
        varampl=str2double(sstring(1));
        seasonal=str2double(sstring(2:3));
    else
        varampl=0;
    end

    if seasonal < 1 || seasonal >floor(s/2)
        error('FSDA:simulateTS:WrongInput',['Seasonal component must be an integer between 1 and ' num2str(floor(s/2))])
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
end

if nocheck == false
    % Controls on the seasonal  component
    if isempty(seasonal) && ~isempty(seasonalb)
        disp('Warning: option seasonal has been specified but the beta coefficients for the seasonal have not been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    if isempty(seasonalb) && ~isempty(seasonal)
        disp('Warning: option seasonal has not been specified but the beta coefficients for the seasonal have been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    if ~isempty(seasonal) && length(seasonalb)~=nseaso+varampl
        % Define matrix which contains linear
        disp(['Warning: option seasonal has been specified but the length of ' ...
            'beta coefficients for the seasonal component is not in agreeement'])
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end
end


%% checks on the AR component
if ARp==0; ARp=[]; end
if nocheck == false
    % Controls on the AR component
    if isempty(ARp) && ~isempty(ARb)
        disp('Warning: option ARp has been specified but the beta coefficients for the autoregressive terms have not been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    if isempty(ARb) && ~isempty(ARp)
        disp('Warning: option seasonal has not been specified but the beta coefficients for the seasonal have been specified')
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end

    if ~isempty(ARp) && length(ARb)~=length(ARp)
        % Define matrix which contains linear
        disp(['Warning: option ARp has been specified but the length of ' ...
            'beta coefficients for the autoregressive component is not in agreeement'])
        error('FSDA:simulateTS:WrongInput','Specify the requested component together with the requested coefficients')
    end
end


%% checks on residuals definition
if isempty(residuals)
    if isempty(sigmaeps) && isempty(signal2noiseratio)
        signal2noiseratio=1; % No input for residual -> default setting
        disp(['None of the inputs for simulating residuals has been specified. ' ...
            'The default value signal2noiseratio=1 will be considered.'])
    elseif ~isempty(sigmaeps) && ~isempty(signal2noiseratio)
        disp(['Warning: both sigmaeps and signal2noiseratio have been specified. ' ...
            'The value of signal2noiseratio will be ignored.'])
    end
else
    if nocheck == false
        if size(residuals,1)~=T || size(residuals,2)~=nsim
            disp('Warning: The residual matrix must be T x nsim')
            error('FSDA:simulateTS:WrongInput','Specify a residual matrix with the correct dimension')
        end
    end
    if ~isempty(sigmaeps) || ~isempty(signal2noiseratio)
        disp(['Warning: even though one between sigmaeps and signal2noiseratio has been specified, ' ...
            'only the input residuals will be considered in the simulation.'])
    end
end


%%  Define the explanatory variable associated to the level shift component
if lshift>0
    if lshift>T-1
        disp('Warning: podition of level shift cannot be greater than T-1')
        error('FSDA:simulateTS:WrongInput','Wrong position of level shift component')
    end
    % Xlshift = explanatory variable associated with level shift Xlshift is
    % 0 up to lsh-1 and 1 from lsh to T
    Xlshift= [zeros(lshift-1,1);ones(T-lshift+1,1)];
end

if varampl>0
    seasonalb_linpart=seasonalb(1:end-varampl);
    seasonalb_nonlinpart=seasonalb(end-varampl+1:end);
else
    seasonalb_linpart=seasonalb;
    seasonalb_nonlinpart=[];
end

%% Compute the underlying components of the signal of the time series
% beta0 = vector which will contain all the coefficients of the model.
% The order of the coefficients is:
% - trend component,
% - linear part of the seasonal component,
% - explanatory variables,
% - non linear part of the seasonal component,
% - level shift
% - beta=[trendb(:); seasonalb_linpart(:); Xb(:); seasonalb_nonlinpart(:); lshiftb(:)];

yhattrend = Xtrend * trendb(:);

if seasonal > 0
    yhatseaso = Xseaso * seasonalb_linpart(:);

    if varampl > 0
        Xtre      = 1 + Seq(:,2:varampl+1)*seasonalb_nonlinpart(:);
        yhatseaso = Xtre .* yhatseaso;
    end
else
    yhatseaso=0;
end

if isemptyX
    yhatX=0;
else
    % Note the order of coefficients is trend, linear part of seasonal
    % component, expl variables, non linear part of seasonal component,
    % level shift
    yhatX = X*Xb(:);
end

if lshift >0
    yhatlshift = lshiftb * Xlshift;
else
    yhatlshift=0;
end

% Simulated values, formed by: trend (yhattrend), (time varying) seasonal
% (yhatseaso), explanatory variables (yhatX) and level shift (yhatlshift)
signal = yhattrend + yhatseaso + yhatX + yhatlshift;

%% Final simulated time series y, with or without autoregressive component
if isempty(ARb)
    % No autoregressive part
    if isempty(residuals)
        % Calculate the std err of the residuals (if not provided)
        if isempty(sigmaeps)
            varsignal = var(signal);
            if varsignal>0
                sigmaeps= sqrt(varsignal/signal2noiseratio);
            else
                sigmaeps= sqrt(1/signal2noiseratio);
            end
        end
        % Simulate the irregular component
        residuals = sigmaeps * randn(T,nsim);
    end
    % y is the final simulated time series (signal + irregular component)
    irregular = residuals;
    y = signal*ones(1,nsim) + irregular;
else
    % Add autoregressive part to the irregular.
    % Be reasonable with the number of autoregressive components ...
    if length(ARb)>6
        disp('Number of autoregressive components is too big and can create model instability: it is set to 6');
        ARp=ARp(1:6);
        ARb=ARb(1:6);
    end
    ARb_all=zeros(1,max(ARp));
    ARb_all(ARp)=ARb;
    if ARIMAX==false
        % AR components on residuals
        % Calculate the std err of the residuals (if not provided)
        if isempty(residuals)
            if isempty(sigmaeps)
                varsignal = var(signal);
                if varsignal>0
                    sigmaeps= sqrt(varsignal/signal2noiseratio);
                else
                    sigmaeps= sqrt(1/signal2noiseratio);
                end
            end
        end
        if exist('regARIMA','file') && isempty(residuals)>0
            % Use the Econometric toolbox if present: regARIMA
            Mdl1 = regARIMA('Intercept',0,'AR',num2cell(ARb_all),'Beta',1,'Variance',sigmaeps^2);
            % arima generates ARIMAX models
            % Mdl1 = arima('AR',num2cell(ARb), 'Beta',1,'Variance',sigmaeps^2);
            [y , irregular] = simulate(Mdl1,T,'X',signal,'NumPaths',nsim);
        else
            % Simulate the data y(t) without the Econometrics toolbox.
            % For more details see files
            % R2019b\toolbox\econ\econ\@regARIMA\simulate.m
            % R2019b\toolbox\econ\econ\@ARIMA\simulate.m
            % R2019b\toolbox\econ\econ\+internal\+econ\simulateStandardizedVariates.m
            LagsAR = 1:length(ARb_all);
            LagsMA = 0;
            coefficients = [0  ARb_all  1];  % ARIMA coefficient vector
            I      = ones(1,nsim);
            maxPQ  = length(ARb_all);
            Y      = zeros(T+maxPQ,nsim);
            E      = Y;
            if isempty(residuals)
                Z      = randn(T,nsim);
                E((maxPQ + 1:end),:)=Z * sigmaeps;
            else
                E((maxPQ + 1:end),:)=residuals;
            end
            for t = (maxPQ + 1):T+maxPQ
                data   = [I;  Y(t - LagsAR,:);  E(t - LagsMA,:)];
                Y(t,:) = coefficients * data;
            end
            Y = Y(maxPQ + 1:T+maxPQ,:);
            irregular = E(maxPQ+1:T+maxPQ,:);
            y = signal*ones(1,nsim) + Y;
        end
    else
        % AR components on dependent variable
        if ~isempty(residuals) % residuals provided
            [y,irregular,~]=simulateARIMAX(T,nsim,ARb_all,residuals,signal);
        elseif ~isempty(sigmaeps) % sigmaeps fixed
            [y,irregular,~]=simulateARIMAX(T,nsim,ARb_all,sigmaeps,signal);
        else % search sigmaeps depending on signal2noiseratio
            varSignal=var(signal);
            if varSignal>eps
                disp('Iterative search of sigmaeps depending on the desired signal2noise ratio.')
                S2NR=zeros(1000,1);
                se=zeros(1001,1);
                se(1)=varSignal;
                niter=0;
                conv=0;
                while conv==0 && niter<1000
                    niter=niter+1;
                    [y,irregular,S2NR(niter)]=simulateARIMAX(T,nsim,ARb_all,se(niter),signal);
                    if abs(S2NR(niter)-signal2noiseratio)/signal2noiseratio<1e-5
                        conv=1;
                        sigmaeps=se(niter);
                    else
                        c=S2NR(niter)/signal2noiseratio;
                        se(niter+1)=se(niter)*c;
                        if se(niter+1)>10e50
                            conv=-1;
                        end
                    end
                end
                if conv==1
                    disp(['After ' num2str(niter) ' iterations, the desired value of signal2noise ratio has been reached with sigmaeps = ' num2str(sigmaeps)])
                elseif conv==-1
                    rows=(S2NR>eps);
                    [~,pos]=min(abs(signal2noiseratio-S2NR(rows)));
                    disp('Warning: the desired signal2noiseratio is too low, impossible to obtain without sigmaeps going to infinity.')
                    disp(['The value of signal2noise ratio closest to the desired one is ' num2str(S2NR(pos)) ', obtained with sigmaeps = ' num2str(se(pos))])
                else
                    [~,pos]=min(abs(signal2noiseratio-S2NR));
                    sigmaeps=se(pos);
                    disp(['After ' num2str(niter) ' iterations, the value of signal2noise ratio closest to the desired one is ' num2str(S2NR(pos)) ','])
                    disp(['obtained with sigmaeps = ' num2str(sigmaeps)])
                end
            else
                disp('The signal has no variance. Simulations are generated with sigmaeps = 1.')
                [y,irregular,~]=simulateARIMAX(T,nsim,ARb_all,1,signal);
            end
        end

    end
end


%% Output structure
out=struct;
out.y=y;
out.signal=signal;
out.trend=yhattrend;
out.seasonal=yhatseaso;
out.X=yhatX;
out.lshift=yhatlshift;
out.irregular=irregular;
out.model = model;

% Write to file the simulated data
if ~isempty(FileNameOutput)
    writematrix(y,FileNameOutput,'delimiter',';')
end

% label the x axis using appropriate dates
if ~isempty(StartDate)
    IniYear=StartDate(1);
    FinalYear=IniYear+ceil(T/s);
    [years, months] = meshgrid(IniYear:FinalYear, 1:12/s:12);
    years=years(1:T);
    months=months(1:T);
    % Convert date and time to serial date number
    % datesnumeric=datenum(years(:), months(:), 1);
    datesnumeric=datetime(years(:), months(:), 1);
    % datetime(StartDate(1),StartDate(2),1)+calmonths(0:(T-1))
else
    datesnumeric=1:T;
end

%% Create plots
if plots==1
    % if nsim>1, we plot only the first simulated series
    % some general plot settings
    vlt15 = verLessThan('matlab', '7.15');
    %clr = 'bkrgmcy';
    %syb = {'-','--','-.',':','-','--','-.'};
    FontSize    = 12;
    SizeAxesNum = 12;

    if samescale
        % yscale to keep uniform across the plots
        [minV,maxV]=minmax(y(:,1),signal,yhattrend,yhatseaso,yhatlshift,yhatX,y(:,1)-signal);
    end

    % Minimum value for xlim
    if isempty(StartDate)
        minxlim=1;
        maxxlim=T;
    else
        minxlim=min(datesnumeric);
        maxxlim=max(datesnumeric);
    end

    % Time series + fitted values
    sb1 = subplot(2,3,1);
    plot(datesnumeric,y(:,1));
    if samescale
        ylim([minV,maxV]);
    end
    title({'Final simulated data',''},'interpreter','none','FontSize',FontSize+2);
    xlim([minxlim,maxxlim]);
    if ~isempty(StartDate)
        if ~verLessThanFS('8.4')
            set(gca,'XTickLabelRotation',90);
        end
        datetick('x','mmm-yy','keepticks');
    end

    sb2 = subplot(2,3,2);
    plot(datesnumeric,signal);
    if samescale
        ylim([minV,maxV]);
    end

    title({'TR+SE+LS+X',''},'interpreter','none','FontSize',FontSize+2);
    xlim([minxlim,maxxlim]);
    if ~isempty(StartDate)
        if ~verLessThanFS('8.4')
            set(gca,'XTickLabelRotation',90);
        end
        datetick('x','mmm-yy','keepticks');
    end

    sb3 = subplot(2,3,3);
    plot(datesnumeric,yhattrend);
    if samescale
        ylim([minV,maxV]);
    end

    title({'Trend (TR)',''},'interpreter','none','FontSize',FontSize+2);
    xlim([minxlim,maxxlim]);
    if ~isempty(StartDate)
        if ~verLessThanFS('8.4')
            set(gca,'XTickLabelRotation',90);
        end
        datetick('x','mmm-yy','keepticks');
    end

    sb4 = subplot(2,3,4);
    plot(datesnumeric,yhatseaso);
    if samescale
        ylim([minV,maxV]);
    end
    xlim([minxlim,maxxlim]);
    if isdatetime(datesnumeric)
        datetick('x','mmm-yy','keepticks');
        if ~verLessThanFS('8.4')
            set(gca,'XTickLabelRotation',90);
        end
    end

    title({'Seasonal (SE)',''},'interpreter','none','FontSize',FontSize+2);

    sb5 = subplot(2,3,5);
    plot(datesnumeric,yhatlshift);
    if samescale
        ylim([minV,maxV]);
    end

    title({'Level shift (LS)',''},'interpreter','none','FontSize',FontSize+2);
    xlim([minxlim,maxxlim]);
    if ~isempty(StartDate)
        if ~verLessThanFS('8.4')
            set(gca,'XTickLabelRotation',90);
        end
        datetick('x','mmm-yy','keepticks');
    end

    sb6 = subplot(2,3,6);
    if yhatX~=0
        plot(datesnumeric,yhatX);
        if samescale
            ylim([minV,maxV]);
        end
        xlim([minxlim,maxxlim]);

        %  xlim([minxlim,T]);
        title({'Explanatory variables (X)',''},'interpreter','none','FontSize',FontSize+2);
    else
        plot(datesnumeric,y(:,1)-signal);
        if samescale
            ylim([minV,maxV]);
        end
        xlim([minxlim,maxxlim]);

        title({'Irregular (I)',''},'interpreter','none','FontSize',FontSize+2);
    end
    if ~isempty(StartDate)
        if ~verLessThanFS('8.4')
            set(gca,'XTickLabelRotation',90);
        end
        datetick('x','mmm-yy','keepticks');
    end

    if ~vlt15
        set([sb1 , sb2, sb3, sb4, sb5, sb6] ,'FontSize',SizeAxesNum,'Box','on','BoxStyle','full');
    else
        set([sb1 , sb2, sb3, sb4, sb5, sb6] ,'FontSize',SizeAxesNum,'Box','on');
    end

end
    function [minV,maxV]=minmax(varargin)
        % returns the minimum and the maximum of all vectors in varargin,
        % omitting vectors with only one value (scalars).
        N = nargin;
        minmaxout = nan(N,2);
        for i1 = 1:N
            p1 = varargin{i1};
            if length(p1)>1
                minmaxout(i1,:) = [min(p1) max(p1)];
            end
        end
        minV = min(minmaxout(:,1),[],'omitnan');
        maxV = max(minmaxout(:,2),[],'omitnan');
    end
end


function [Y,irregular,s2nr]=simulateARIMAX(T,nsim,ARcoeff,errorDef,signal)
ARterms=zeros(length(ARcoeff),nsim);
coefficients = [1 ARcoeff];  % ARIMA coefficient vector
Y=zeros(T,nsim);
if max(size(errorDef))==1
    Z=randn(T,nsim);
    irregular =Z*errorDef;
    varNoise=errorDef^2;
else
    irregular=errorDef;
    varNoise=mean(var(irregular));
end
for t = 1:T
    data   = [signal(t)*ones(1,nsim); ARterms];
    Y(t,:) = coefficients*data +irregular(t,:);
    ARterms=[Y(t,:); ARterms(1:end-1,:)];
end
varSignal=mean(var(Y-irregular));
s2nr=varSignal/varNoise;
end



%FScategory:REG-Regression