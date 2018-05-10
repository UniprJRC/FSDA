function [outFORE] = forecastTS(outEST,varargin)
%Forecast for a time series with trend, time varying seasonal, level shift and irregular component
%
%<a href="matlab: docsearchFS('forecastTS')">Link to the help function</a>
%
% forecastTS produces forecasts for a time series with trend (up to third order),
% seasonality (constant or of varying amplitude) with a different number of
% harmonics and a level shifta and explanatory variables.
%
%  Required input arguments:
%
%  outEST :     A structure containing the output of routine LTSts.
%               Structure.
%               Structure containing the following fields.
%          outEST.B =   Matrix containing estimated beta coefficients,
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
%          outEST.posLS = scalar associated with best tentative level shift
%                       position. If this field does not exist, forecasts
%                       are done assuming no level shift.
%         outEST.invXX = $cov(\beta)/\hat \sigma^2$. p-by-p, square matrix.
%                       If the model is linear out.invXX  is equal to
%                       $(X'X)^{-1}$, else out.invXX is equal to $(A'A)^{-1}$
%                       where $A$ is the matrix of partial derivatives. More
%                       precisely:
%                       \[
%                       a_{i,j}=\frac{\partial \eta_i(\hat \beta)}{\partial \hat \beta_j}
%                       \]
%                       where
%                       \begin{eqnarray}
%                       y_i & = & \eta(x_i,\beta)+ \epsilon_i  \\
%                           & = & \eta_i +\epsilon_i \\
%                           & = & \eta(x_i,\hat \beta)+ e_i  \\
%                           & = & \hat \eta_i + e_i
%                       \end{eqnarray}
%         outEST.yhat = vector of fitted values after final (NLS=non linear
%                       least squares) step:
%                       $ (\hat \eta_1, \hat \eta_2, \ldots, \hat \eta_T)'$
%        outEST.scale = Final scale estimate of the residuals
%                     \[
%                     \hat \sigma = cor \times \sum_{i \in S_m} [y_i- \eta(x_i,\hat \beta)]^2/(m-p)
%                     \]
%                     where $S_m$ is a set of cardinality $m$ which
%                     contains the units not declared as outliers and $p$
%                     is the total number of estimated parameters and cor
%                     is a correction factor to make the estimator
%                     consistent.
%                     REMARK: structure outEST can be conveniently created
%                     by function LTSts.
%                 Data Types - struct
%
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
%               model.s = scalar greater than zero which specifies the
%                       length of the seasonal period. For monthly
%                       data (default) s=12, for quartely data s=4, ...
%                       The default value of model.s is 12 (that is monthly
%                       data are assumed)
%               model.seasonal = scalar (integer specifying number of
%                        frequencies, i.e. harmonics, in the seasonal
%                        component. Possible values for seasonal are
%                        $1, 2, ..., [s/2]$, where $[s/2]=floor(s/2)$.
%                        For example:
%                        if seasonal = 1 (default) we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 sin ( 2 \pi t/s)$;
%                        if seasonal = 2 we have:
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
%                        For example, seasonal = 201 implies a seasonal
%                        component which just uses one frequency
%                        which grows in a quadratic way over time as
%                        follows:
%                        $(1+\beta_3 t + \beta_4  t^2)\times( \beta_1 \cos(
%                        2 \pi t/s) + \beta_2 \sin ( 2 \pi t/s))$.
%                       If this field is an empty double (default) the
%                       simulated time series will not contain a seasonal
%                       component.
%               model.X  = explanatory variabels. Matrix of size T-by-nexpl. If model.X
%                       is a matrix of size T-by-nexpl, it contains the
%                       values of nexpl extra covariates which
%                       affect y.
%                       If this field is an empty double (default) there is
%                       no effect of explanatory variables.
%                 Example - 'model', model
%                 Data Types - struct
%
%       nfore  : number of forecasts. Scalar.
%               Positive integer which defines the number of forecasts. The
%               default value of nfore is 24.
%               Example - 'nfore',12
%               Data Types - double
%
%      conflev : confidence level for the confidence bands. Scalar.
%                A number between 0 and 1 which defines the confidence
%                level which is used to produce the bands. The default
%                value of conflev is 0.99.
%               Example - 'conflev',0.999
%               Data Types - double
%
%       plots : Plots on the screen. Scalar.
%               If plots == 1 a plot with the real time series  with fitted
%               values and forecasts (with confidence bands) will appear on
%               the screen.
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
%
% FileNameOutput : save simulated time series to txt file. Character.
%               If FileNameOutput is empty (default) nothing is saved on
%               the disk, else FileNameOutput will contain the path where
%               to save the file on the disk.
%               Example - 'FileNameOutput',['C:' filesep 'myoutput' fielsep 'savesimdata.txt']
%               Data Types - Character
%
%  dispresults : Display results of final fit. Boolean. If dispresults is
%               true,  labels of coefficients, estimated coefficients,
%               standard errors, tstat and p-values are shown on the
%               screen in a fully formatted way. The default value of
%               dispresults is false.
%               Example - 'dispresults',true
%               Data Types - logical
%
%  Output:
%
%     outFORE:   structure which contains the following fields:
%
%                outFORE.signal = vector of length (length(y)+nfore) containing
%                   predictive values in sample and out of sample.
%                   Predictive values = TR+SE+X+LS.
%                outFORE.trend = vector of length (length(y)+nfore) containing
%                   estimated trend (TR) in sample and out of sample.
%                   If this component is not present, it is equal to 0.
%                outFORE.seasonal = vector of length (length(y)+nfore) containing
%                   estimated seasonal component (SE) in sample and out of sample.
%                   If this component is not present, it is equal to 0.
%                outFORE.lshift = vector of length (length(y)+nfore) containing
%                   level shift (LS) in sample and out of sample.
%                   If this component is not present, it is equal to 0.
%                outFORE.X = vector of length (length(y)+nfore)
%                   containing the effecf of the explanatory variables.
%                   If this component is not present, it is equal to 0.
%                outFORE.confband  = matrix of size (length(y)+nfore)-by-2
%                   containing lower and upper confidence bands of the
%                   forecasts. The confidence level of the bands is
%                   splecified is input parameter conflev. Note that the
%                   first length(y) rows of this matrix are equal to NaN.
%               outFORE.datesnumeric = vector of length (length(y)+nfore)
%                   containing the dates in numeric format.
%
%
%
% See also LTSts, wedgeplot, simulateTS
%
% References:
%
% Rousseeuw, P.J., Perrotta D., Riani M., Hubert M. (2018), Robust
% Monitoring of Time Series with Application to Fraud Detection,
% Submitted.
%
%
% Copyright 2008-2017.
% Written by Marco Riani, Domenico Perrotta, Peter Rousseeuw and Mia Hubert
%
%
%<a href="matlab: docsearchFS('forecastTS')">Link to the help function</a>
%
%$LastChangedDate:: 2018-02-19 17:38:15 #$: Date of the last commit

% Examples:

%{
    %% Linear time varying seasonal component.
    close all
    rng(1)
    model=struct;
    model.trend=1;
    model.seasonal=103;
    modelSIM=model;
    modelSIM.trendb=[0 0];
    modelSIM.seasonalb=40*[0.1 -0.5 0.2 -0.3 0.3 -0.1 0.222];
    modelSIM.signal2noiseratio=20;
    T=100;
    % Simulate
    outSIM=simulateTS(T,'model',modelSIM,'plots',1);
    ySIM=outSIM.y;
    % Estimate
    outEST=LTSts(ySIM,'model',model,'plots',1);
    % Forecast
    outFORE=forecastTS(outEST,'model',model,'plots',1);
%}

%{
    %% Quadratic trend and constant seasonal.
    close all
    rng(1)
    model=struct;
    model.trend=2;
    model.seasonal=3;
    modelSIM=model;
    modelSIM.trendb=[100 10 -0.05];
    modelSIM.seasonalb=400*[0.1 -0.5 0.2 -0.3 0.3 -0.1];
    modelSIM.signal2noiseratio=1;
    T=100;
    % Simulate
    outSIM=simulateTS(T,'model',modelSIM,'plots',1);
    ySIM=outSIM.y;
    % Estimate
    outEST=LTSts(ySIM,'model',model,'plots',1);
    % Forecast
    outFORE=forecastTS(outEST,'model',model,'plots',1);
%}

%{
    %% Simulated time series with quadratic trend, fixed seasonal and level shift.
    % A time series of 100 observations is simulated from a model which
    % contains a quadratic trend, a seasonal component with two harmonics
    % no explanatory variables and a level shift in position 30 with size
    % 5000 and a signal to noise ratio egual to 20
    close all
    rng(1)
    model=struct;
    model.trend=2;
    model.seasonal=2;
    model.lshift=30;
    modelSIM=model;

    modelSIM.trendb=[5,10,-3];
    modelSIM.seasonalb=100*[2 4 0.1 8];
    modelSIM.signal2noiseratio=20;
    modelSIM.lshiftb=10000;

    T=100;
    % Simulate
    outSIM=simulateTS(T,'model',modelSIM,'plots',1);
    ySIM=outSIM.y;
    % Estimate
    %  model.lshift=5 implies that LS is investigated from position 5
    model.lshift=5;
    outEST=LTSts(ySIM,'model',model,'plots',1,'msg',0);
    % Forecast
    outFORE=forecastTS(outEST,'model',model,'plots',1);
%}

%{
    %% Contaminated airline data (1).
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417   % Jan
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
    y=y(:);
    % Contaminate the first 20 observations
    y(1:20)=y(1:20)+200;
    close all
    % Model with linear trend, three harmonics for seasonal component and
    % varying amplitude using a linear trend. Search for a level shift
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=103;         % three harmonics with linear time varying seasonality
    model.lshift=10;            % search for level shift
    out=LTSts(y,'model',model,'plots',1,'dispresults',true,'msg',0);

    % 3 years forecasts
    nfore=36;
    StartDate=[1949 1];
    conflev=0.999; % Wide confidence level for the forecast
    outFORE=forecastTS(out,'model',model,'nfore',nfore,'StartDate',StartDate,'conflev',conflev);
%}

%{
    %% Contaminated airline data (2).
    close all
    % In this example we estimate a model without the seasonal component
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417   % Jan
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
    y=y(:);
    % Contaminate the first 20 observations
    y(1:20)=y(1:20)+200;

    % Model with linear trend and no seasonal component. Search for a level shift
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=[];          % no seasonal component
    model.lshift=10;            % search for level shift
    out=LTSts(y,'model',model,'plots',1,'dispresults',true,'msg',0);

    % 3 years forecasts
    nfore=36;
    StartDate=[1949 1];
    conflev=0.999; % Wide confidence level for the forecast
    outFORE=forecastTS(out,'model',model,'nfore',nfore,'StartDate',StartDate,'conflev',conflev);
%}

%{
    % Example with simulated data.
    % Simulate data with linear trend, time varying seasonal component, and
    % level shift
    rng(1)
    model=struct;
    model.trend=1;
    model.trendb=[1 1];
    model.seasonal=103;
    model.seasonalb=40*[0.5 -0.5 0.3 -0.3 0.1 -0.1 0.222];
    model.lshift=40;
    model.lshiftb=13000;
    model.signal2noiseratio=20;
    T=150;
    FileNameOutput=[pwd filesep 'ysimout.txt'];
    outSIM=simulateTS(T,'model',model,'FileNameOutput',FileNameOutput,...
        'plots',true,'samescale',true);
    y=outSIM.y;
    % Data contamination
    y(131:140)=y(131:140)-29000;

    % Estimation
    modelEST=struct;
    modelEST.trend=1;
    modelEST.seasonal=103;
    modelEST.lshift=30;
    outEST=LTSts(y,'model',modelEST,'dispresults',true,'plots',0);

    % Forecasting
    % nfore= number of forecasts;
    nfore=36;
    % Forecasts with a 99.9 per cent confidence level
    OUTfore=forecastTS(outEST,'model',modelEST,'nfore',nfore,'conflev',0.999);
%}


%% Beginning of code
if nargin<1
    error('FSDA:forecastTS:MissingInputs','Input structure is missing');
end

% Set up values for default model
modeldef          = struct;
modeldef.trend    = 1;
modeldef.s        = 12;       % monthly time series
modeldef.seasonal = [];
modeldef.X        = [];       % no explanatory variables
modeldef.lshift   = [];       % no level shift
nocheck           = false;
plots             = 1;
FileNameOutput    = '';
StartDate         = '';
dispresults     =false;
nfore           =24;    % number of forecasts
conflev         = 0.99; % default confidence level for the forecasts

options=struct('model',modeldef,...
    'dispresults',dispresults,'nfore',nfore,'plots',plots,...
    'FileNameOutput',FileNameOutput,'StartDate',StartDate,'conflev',conflev);


%% User options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:forecastTS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present in
    % structure options Remark: the nocheck option has already been dealt
    % by routine chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:forecastTS:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    plots=options.plots;
    dispresults=options.dispresults;
    nfore=options.nfore;
    FileNameOutput=options.FileNameOutput;
    StartDate=options.StartDate;
    conflev=options.conflev;
end

% Default values for the optional parameters are set inside structure
% 'options'
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

model = modeldef;

% Get model parameters
trend    = model.trend;       % get kind of  trend
s        = model.s;           % get periodicity of time series
seasonal = model.seasonal;    % get number of harmonics

if isfield(outEST,'posLS')
    lshift   = outEST.posLS;
else
    lshift=0;
end

y=outEST.y;
n=length(y);
T = n+nfore;

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
        error('FSDA:forecastTS:WrongInput',['Seasonal component must be an integer between 1 and ' num2str(floor(s/2))])
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
    Xlshift= [zeros(outEST.posLS-1,1);ones(T-outEST.posLS+1,1)];
else
    Xlshift =[];
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
% + 2 additional parameters is there is a level shift component
p=pini+varampl+(lshift>0);

% Now compute again vector yhat using final vector betaout
bsb=seq;
betaout=outEST.B(:,1);

if length(betaout)~=p
    disp('Warning: number of supplied regression parameters is not in agreement with those of input structure model')
    disp(['Number of supplied regression parameters=' num2str(length(betaout))])
    disp(['Number of parameters in input strcuture model=' num2str(p)])
    disp('...')
    error('FSDA:forecast:WrongInput','Wrong input model')
end

[yhat,yhattrend,yhatseaso,yhatX,yhatlshift]=lik(betaout);

% yfitFS = likyhat(betaout,Xtrendf);
yfitFS = yhat;

if varampl==0
    J=[Xtrend Xseaso X Xlshift];
else
    fdiffstep=[];
    J = getjacobianFS(betaout,fdiffstep,@lik,yfitFS);
end

confband=NaN(length(yhat),2);
if ~isempty(conflev)
    invXX=outEST.invXX;
    quant=tinv(1-(1-conflev)/2,n-length(betaout));
    for i=n+1:n+nfore
        se=quant*outEST.scale*sqrt(J(i,:)*invXX*(J(i,:)'));
        confband(i,1)=yhat(i)-se;
        confband(i,2)=yhat(i)+se;
    end
end

outFORE=struct;
outFORE.signal=yhat;
outFORE.trend=yhattrend;
outFORE.seasonal=yhatseaso;
outFORE.X=yhatX;
outFORE.lshift=yhatlshift;
outFORE.confband=confband;

yhatwithbands=[yhat confband];

% Write to file the simulated data
if ~isempty(FileNameOutput)
    dlmwrite(FileNameOutput,yhatwithbands,'delimiter','\t','precision','%.12f')
end

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
    datesnumeric=1:T;
end

% Raw numbers associated with the dates on the x axis
outFORE.datesnumeric=datesnumeric;


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
    
    bhat=outEST.B(:,1);
    se=outEST.B(:,2);
    tstat=outEST.B(:,3);
    pval=outEST.B(:,4);
    if verLessThan ('matlab','8.2.0')
        disp('           Coeff.     SE         t-stat       p-values');
        disp( [char(lab) num2str([bhat se tstat pval])]);
    else
        disp([table(lab) table(bhat) table(se) table(tstat) table(pval)]);
    end
    if lshift>0
        disp(['Level shift position t=' num2str(outEST.posLS)])
    end
end

if plots==1
    figure;
    yfore=outFORE.signal;
    % Plot original time series
    plot(datesnumeric(1:n),y,'k')
    hold('on')
    % plot the original series
    plot(datesnumeric(1:n),yfore(1:n),'b-')
    % plot the forecasts
    plot(datesnumeric(n+1:n+nfore),yfore(n+1:end),'r')
    
    % plot the signal (TR+LS+X)
    % plot(datesnumeric,outFORE.trend+outFORE.lshift+outFORE.X,'color','m')
    
    title('Fit and forecasts from LTS','interpreter','LaTex','FontSize',14)
    plot(datesnumeric(n+1:n+nfore),confband(n+1:end,:),'r--')
    if ~isempty(StartDate)
        datetick('x','mmm-yy');
        if ~verLessThanFS(8.4)
            set(gca,'XTickLabelRotation',90);
        end
    end
    ax=axis;
    ylimits=ax(3:4);
    line(0.5*sum(datesnumeric(n:n+1))*ones(2,1),ylimits,'color','g')
    
end

    function [yhat,yhattrend,yhatseaso,yhatX,yhatlshift]=lik(beta0)
        
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
        else
            yhatseaso=0;
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
%FScategory:REG-Regression