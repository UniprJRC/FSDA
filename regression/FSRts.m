function [out]=FSRts(y,varargin)
%FSRts is an automatic adaptive procedure to detect outliers in time series
%
%<a href="matlab: docsearchFS('FSRts')">Link to the help function</a>
%
% Required input arguments:
%
%           y:  Time series to analyze. Vector. A row or a column vector
%               with T elements, which contains the time series.
%
% Optional input arguments:
%
%      model :  Model type. Structure. A structure which specifies the model
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
%                         values of parameter estimates which have to be
%                         used in the maximization procedure. If model.B is
%                         a matrix, then initial estimates are extracted
%                         from the first colum of this matrix. If this
%                         field is empty or if this field is not present,
%                         the initial values to be used in the maximization
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
%       nsamp   : Number of subsamples which will be extracted to find the
%                 robust estimator. Scalar. If nsamp=0 all subsets will be extracted.
%                 They will be (n choose p).
%                 If the number of all possible subset is <1000 the
%                 default is to extract all subsets otherwise just 1000.
%                 Example - 'nsamp',1000
%                 Data Types - double
%
%       lms     : Criterion to use to find the initial subset to initialize
%                 the search. Scalar,  vector or structure. lms specifies
%                 the criterion to use to find the initial subset to
%                 initialize the search (LTS with concentration steps or
%                 subset supplied directly by the user).
%                 The default value is 1 (Least trimmed squares
%                 is computed to initialize the search).
%                 If lms is a struct it is possible to control a
%                 series of options for concentration steps (for more
%                 details see function LTSts.m)
%                 If, on the other hand, the user wants to initialize the
%                 search with a prespecified set of units there are two
%                 possibilities:
%                 1) lms can be a vector with length greater than 1 which
%                 contains the list of units forming the initial subset. For
%                 example, if the user wants to initialize the search
%                 with units 4, 6 and 10 then lms=[4 6 10];
%                 2) lms is a struct which contains a field named bsb which
%                 contains the list of units to initialize the search. For
%                 example, in the case of simple regression through the
%                 origin with just one explanatory variable, if the user
%                 wants to initialize the search with unit 3 then
%                 lms=struct; lms.bsb=3;
%                 Example - 'lms',1
%                 Data Types - double
%
%           h   : The number of observations that have determined the least
%                 trimmed squares estimator. Scalar. h is an integer
%                 greater or equal than p but smaller then n. Generally if
%                 the purpose is outlier detection h=[0.5*(n+p+1)] (default
%                 value). h can be smaller than this threshold if the
%                 purpose is to find subgroups of homogeneous observations.
%                 In this function the LTSts estimator is used just to
%                 initialize the search. The default value of h which is
%                 used is round(0.75*T).
%                 Example - 'h',round(n*0.55)
%                 Data Types - double
%
%       plots   : Plot on the screen. Scalar.
%                 If plots=1 (default) the plot of minimum deletion
%                 residual with envelopes based on T observations and the
%                 scatterplot matrix with the outliers highlighted is
%                 produced together with a two panel plot.
%                 The upper panel contains the orginal time series with
%                 fitted values. The bottom panel will contain the plot
%                 of robust residuals against index number. The confidence
%                 level which is used to draw the horizontal lines associated
%                 with the bands for the residuals is 0.999.
%                 If plots=2 the user can also monitor the intermediate
%                 plots based on envelope superimposition.
%                 Else no plot is produced.
%                 Example - 'plots',1
%                 Data Types - double
%
%      init :   Start of monitoring point. Scalar.
%               It specifies the point where to initialize the search and
%               start monitoring required diagnostics. If it is not
%               specified it is set equal floor(0.5*(T+1))
%               Example - 'init',100 starts monitoring from step m=100
%               Data Types - double
%
%       nocheck : Check input arguments inside input option model.
%               As default nocheck=false.
%               Example - 'nocheck',true
%               Data Types - boolean
%
%    bivarfit : Superimpose bivariate least square lines. Character. This option adds
%                 one or more least squares lines, based on
%                 SIMPLE REGRESSION of y on Xi, to the plots of y|Xi.
%                 bivarfit = ''
%                   is the default: no line is fitted.
%                 bivarfit = '1'
%                   fits a single ols line to all points of each bivariate
%                   plot in the scatter matrix y|X.
%                 bivarfit = '2'
%                   fits two ols lines: one to all points and another to
%                   the group of the genuine observations. The group of the
%                   potential outliers is not fitted.
%                 bivarfit = '0'
%                   fits one ols line to each group. This is useful for the
%                   purpose of fitting mixtures of regression lines.
%                 bivarfit = 'i1' or 'i2' or 'i3' etc. fits
%                   an ols line to a specific group, the one with
%                   index 'i' equal to 1, 2, 3 etc. Again, useful in case
%                   of mixtures.
%               Example - 'bivarfit','2'
%               Data Types - char
%
%       multivarfit : Superimpose multivariate least square lines. Character.
%                 This option adds one or more least square lines, based on
%                 MULTIVARIATE REGRESSION of y on X, to the plots of y|Xi.
%                 multivarfit = ''
%                   is the default: no line is fitted.
%                 multivarfit = '1'
%                   fits a single ols line to all points of each bivariate
%                   plot in the scatter matrix y|X. The line added to the
%                   scatter plot y|Xi is avconst + Ci*Xi, where Ci is the
%                   coefficient of Xi in the multivariate regression and
%                   avconst is the effect of all the other explanatory
%                   variables different from Xi evaluated at their centroid
%                   (that is overline{y}'C))
%                 multivarfit = '2'
%                   equal to multivarfit ='1' but this time we also add the
%                   line based on the group of unselected observations
%                   (i.e. the normal units).
%               Example - 'multivarfit','1'
%               Data Types - char
%
%      labeladd : Add outlier labels in plot. Character. If this option is
%                 '1',  we label the outliers with the
%                 unit row index in matrices X and y. The default value is
%                 labeladd='', i.e. no label is added.
%               Example - 'labeladd','1'
%               Data Types - char
%
%       nameX  : Add variable labels in plot. Cell array of strings. Cell
%                 array of strings of length p containing the labels of
%                 the variables of the regression dataset. If it is empty
%                 (default) the sequence X1, ..., Xp will be created
%                 automatically
%               Example - 'nameX',{'NameVar1','NameVar2'}
%               Data Types - cell
%
%       namey  :  Add response label. Character. String containing the
%                 label of the response
%               Example - 'namey','NameOfResponse'
%               Data Types - char
%
%       ylim   :  Control y scale in plot. Vector. Vector with two elements
%                 controlling minimum and maximum on the y axis. 
%                 Default value is '' (automatic scale)
%               Example - 'ylim',[0,10] sets the minimum value to 0 and the
%               max to 10 on the y axis
%               Data Types - double
%
%       xlim   : Control x scale in plot. Vector. Vector with two elements
%               minimum and maximum on the x axis. Default value is ''
%               (automatic scale)
%               Example - 'xlim',[0,10] sets the minimum value to 0 and the
%               max to 10 on the x axis
%               Data Types - double
%
%      bonflev  : Signal to use to identify outliers. Scalar. Option to be
%                used if the distribution of the data is
%                 strongly non normal and, thus, the general signal
%                 detection rule based on consecutive exceedances cannot be
%                 used. In this case bonflev can be:
%                 - a scalar smaller than 1 which specifies the confidence
%                   level for a signal and a stopping rule based on the
%                   comparison of the minimum MD with a
%                   Bonferroni bound. For example if bonflev=0.99 the
%                   procedure stops when the trajectory exceeds for the
%                   first time the 99% bonferroni bound.
%                 - A scalar value greater than 1. In this case the
%                   procedure stops when the residual trajectory exceeds
%                   for the first time this value.
%                 Default value is '', which means to rely on general rules
%                 based on consecutive exceedances.
%               Example - 'bonflev',0.99
%               Data Types - double
%
%       msg    :  Level of output to display. Scalar. It controls whether
%                 to display or not messages on the screen
%                 If msg==1 (default) messages are displayed on the screen about
%                   step in which signal took place
%                 else no message is displayed on the screen.
%               Example - 'msg',1
%               Data Types - double
%
% bsbmfullrank : Dealing with singluar X matrix. Scalar. This option tells
%                 how to behave in case subset at step m
%                 (say bsbm) produces a singular X. In other words,
%                 this options controls what to do when rank(X(bsbm,:)) is
%                 smaller then number of explanatory variables. If
%                 bsbmfullrank =1 (default) these units (whose number is
%                 say mnofullrank) are constrained to enter the search in
%                 the final n-mnofullrank steps else the search continues
%                 using as estimate of beta at step m the estimate of beta
%                 found in the previous step.
%               Example - 'bsbmfullrank',1
%               Data Types - double
%
%        tag    : tags to the plots which are created. 
%                 character or cell array of characters.
%                 This option enables to add a tag to the plots which are
%                 created. The default tag names are:
%                 fsr_mdrplot for the plot of mdr based on all the
%                 observations;
%                 fsr_yXplot for the plot of y against each column of X
%                 with the outliers highlighted;
%                 fsr_resuperplot for the plot of resuperimposed envelopes. The
%                 first plot with 4 panel of resuperimposed envelopes has
%                 tag fsr_resuperplot1, the second  fsr_resuperplot2 ...
%                 If tag is character or a cell of characters of length 1,
%                 it is possible to specify the tag for the plot of mdr
%                 based on all the observations;
%                 If tag is a cell of length 2 it is possible to control
%                 both the tag for the plot of mdr based on all the
%                 observations and the tag for the yXplot with outliers
%                 highlighted.
%                 If tag is a cell of length 3 the third element specifies
%                 the names of the plots of resuperimposed envelopes.
%                 Example - 'tag',{'plmdr' 'plyXplot'};
%                 Data Types - char or cell
%
% Output:
%
%         out:   structure which contains the following fields
%
% out.ListOut  = k x 1 vector containing the list of the units declared as
%                outliers or NaN if the sample is homogeneous
% out.outliers = out.ListOut. This field is added for homogeneity with the
%                other robust estimators.
% out.beta   = Matrix containing estimated beta coefficients,
%                       standard errors, t-stat and p-values
%                       The content of matrix out.beta is as follows:
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
% out.scale  =  scalar containing the estimate of the scale (sigma).
% out.residuals= n x 1 vector containing the estimates of the robust
%                scaled residuals.
% out.fittedvalues= n x 1 vector containing the fitted values.
% out.mdr    =  (n-init) x 2 matrix
%               1st col = fwd search index
%               2nd col = value of minimum deletion residual in each step
%               of the fwd search
%   Exflag  :   Reason nlinfit stopped. Integer matrix.
%               (n-init+1) x 2 matrix containing information about the
%               result  of the maximization procedure.
%               If the model is non linear out.Exflag(i,2) is equal to 1
%               if at step out.Exflag(i,1) the maximization procedure did not produce
%               warnings or the warning was different from
%               "ILL Conditiioned Jacobian". For any other warning
%               which is produced (for example,
%               "Overparameterized", "IterationLimitExceeded",
%               'MATLAB:rankDeficientMatrix") out.Exflag(i,2) is equal
%               to -1;
% out.Un     =  (n-init) x 11 matrix which contains the unit(s) included
%               in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one.
%               Un(1,2) for example contains the unit included in step
%               init+1.
%               Un(end,2) contains the units included in the final step
%               of the search.
% out.nout    = 2 x 5 matrix containing the number of times mdr went out
%               of particular quantiles.
%               First row contains quantiles 1 99 99.9 99.99 99.999.
%               Second row contains the frequency distribution.
% out.constr  = This output is produced only if the search found at a
%               certain step X is a singular matrix. In this case the
%               search runs in a constrained mode, that is including the
%               units which produced a non singular matrix in the last n-constr
%               steps. out.constr is a vector which contains the list of
%               units which produced a singular X matrix
% out.Exflag  = Reason nlinfit stopped. Integer matrix.
%               (n-init+1) x 2 matrix containing information about the
%               result  of the maximization procedure.
%               If the model is non linear out.Exflag(i,2) is equal to 1
%               if at step out.Exflag(i,1) the maximization procedure did not produce
%               warnings or the warning was different from
%               "ILL Conditiioned Jacobian". For any other warning
%               which is produced (for example,
%               "Overparameterized", "IterationLimitExceeded",
%               'MATLAB:rankDeficientMatrix") out.Exflag(i,2) is equal
%               to -1;
% out.class  =  'FSRts'.
%
% See also: FSR, LTSts, FSRtsmdr
%
% References:
%
% Riani, M., Atkinson, A.C. and Cerioli, A. (2009), Finding an unknown
% number of multivariate outliers, "Journal of the Royal Statistical
% Society Series B", Vol. 71, pp. 201-221.
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [RPRH]
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('FSRts')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-05-31 10:53:11 #$: Date of the last commit

% Examples:

%{
    % FSRts with all default options.
    % Reset the random generator
    rng('default')
    n=200;
    rng(123456);
    outSIM=simulateTS(n,'plots',1);
    % Uncontaminated data
    y=outSIM.y;

    % FSRts and LTSts on Uncontaminated data (conflev defaults differ)
    [outFSRu] = FSRts(y,'plots',1);
    [outLTSu] = LTSts(y,'plots',1,'conflev',0.99);
    outFSRu.outliers
    outLTSu.outliers'

    % Contaminated data
    close all
    ycont=y;
    ycont(10:15) = ycont(10:15)+2*mean(ycont)*sign(ycont(10:15));

    % FSRts and LTSts on contaminated data (conflev defaults differ)
    [outFSR] = FSRts(ycont,'plots',1);
    [outLTS] = LTSts(ycont,'plots',1,'conflev',0.99);
    outFSR.outliers
    outLTS.outliers'
%}

%{
    %% FSRts with optional arguments.
    rng(1)
    model=struct;
    model.trend=[];
    model.trendb=[];
    model.seasonal=103;
    model.seasonalb=40*[0.1 -0.5 0.2 -0.3 0.3 -0.1 0.222];
    model.signal2noiseratio=20;
    T=100;
    outSIM=simulateTS(T,'model',model,'plots',1);
    y=outSIM.y;
    model1=struct;
    model1.trend=1;              % linear trend
    model1.s=12;                 % monthly time series
    model1.lshift=0;             % No level shift
    model1.seasonal=104;         % Four harmonics
    out=FSRts(y,'model',model1);
%}

%{
    % FSRts with optional arguments in time series with outliers.
    % A time series of 100 observations is simulated from a model which
    % contains no trend, a linear time varying seasonal component with
    % three harmonics, no explanatory variables and a signal to noise ratio
    % egual to 20
    rng(1)
    model=struct;
    model.trend=[];
    model.trendb=[];
    model.seasonal=103;
    model.seasonalb=40*[0.1 -0.5 0.2 -0.3 0.3 -0.1 0.222];
    model.signal2noiseratio=20;
    T=100;
    outSIM=simulateTS(T,'model',model,'plots',1);
    y=outSIM.y;
    y(80:90)=y(80:90)+15000;
    model1=struct;
    model1.trend=1;              % linear trend
    model1.s=12;                 % monthly time series
    model1.lshift=0;
    model1.seasonal=104;
    out=FSRts(y,'model',model1);
%}

%{
    %% Example of the use of option lms as a vector.
    % A time series of 100 observations is simulated from a model which
    % contains no trend, a linear time varying seasonal component with
    % three harmonics, no explanatory variables and a signal to noise ratio
    % egual to 20
    rng(1)
    model=struct;
    model.trend=[];
    model.trendb=[];
    model.seasonal=103;
    model.seasonalb=40*[0.1 -0.5 0.2 -0.3 0.3 -0.1 0.222];
    model.signal2noiseratio=20;
    T=100;
    outSIM=simulateTS(T,'model',model);
    y=outSIM.y;
    % Contaminate the series.
    y(80:90)=y(80:90)+15000;
    model1=struct;
    model1.trend=1;              % linear trend
    model1.s=12;                 % monthly time series
    model1.seasonal=104;
    % Initialize the search with the first 20 units.
    out=FSRts(y,'model',model1,'lms',1:20);
%}

%{
    %% Example of the use of option lms as struct.
    % A time series of 100 observations is simulated from a model which
    % contains no trend, a linear time varying seasonal component with
    % three harmonics, no explanatory variables and a signal to noise ratio
    % egual to 20
    rng(1)
    model=struct;
    model.trend=[];
    model.trendb=[];
    model.seasonal=102;
    model.seasonalb=40*[0.1 -0.5 0.2 0.3 0.01];
    model.signal2noiseratio=20;
    model.lshift=30;
    model.lshiftb=2000;
    T=100;
    outSIM=simulateTS(T,'model',model,'plots',1);
    y=outSIM.y;
    % Contaminate the series.
    y(80:90)=y(80:90)+2000;
    model1=struct;
    model1.trend=1;              % linear trend
    model1.s=12;                 % monthly time series
    model1.seasonal=104;
    lms=struct;
    lms.bsb=[1:20 80:85];
    lms.posLS=30;
    % Initialize the search with the units inside lms.bsb.
    out=FSRts(y,'model',model1,'lms',lms);
%}

%{
    %% Automatic outlier and level shift detection.
    % A time series of 100 observations is simulated from a model which
    % contains no trend, a linear time varying seasonal component with
    % three harmonics, no explanatory variables and a signal to noise ratio
    % egual to 20
    rng(1)
    model=struct;
    model.trend=[];
    model.trendb=[];
    model.seasonal=102;
    model.seasonalb=40*[0.1 -0.5 0.2 0.3 0.01];
    model.signal2noiseratio=20;
    model.lshift=30;
    model.lshiftb=2000;
    T=100;
    outSIM=simulateTS(T,'model',model,'plots',1);
    y=outSIM.y;
    % Contaminate the series.
    y(80:90)=y(80:90)+2000;
    model1=struct;
    model1.trend=1;              % linear trend
    model1.s=12;                 % monthly time series
    model1.seasonal=104;
    model1.lshift=-1;
    % Automatically serch for outliers and level shift
    out=FSRts(y,'model',model1,'msg',0);
%}


%% Beginning of code

T=length(y);
nsampdef= 1000;

init=floor(0.5*(T+1));

% Set up values for default model
% Set up values for default model
modeldef         =struct;
modeldef.trend   =1;        % linear trend
modeldef.s       =12;       % monthly time series
modeldef.seasonal=1;        % just one harmonic
modeldef.X       ='';       % no extra explanatory variable
modeldef.lshift  =0;        % no level shift
hdef    = round(0.75*T);
% tag 
tagdef='pl_fsr';

options=struct('nsamp',nsampdef,'model',modeldef,'lms',1,'plots',1,...
    'init',init,'h',hdef,...
    'labeladd','','bivarfit','','multivarfit','',...
    'xlim','','ylim','','nameX','','namey','',...
    'msg',1,'nocheck',false,'bonflev','',...
    'bsbmfullrank',1,'tag',tagdef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRts:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

init=options.init;
h=options.h;
lms=options.lms;
bsbmfullrank=options.bsbmfullrank;
nsamp=options.nsamp;
msg=options.msg;
model=options.model;

%% Start of the forward search

seq=1:T;

iter=0;

% Use as initial subset the one supplied by the user or the best according
% to LMS or LTS

if length(lms)>1 || (isstruct(lms) && isfield(lms,'bsb'))
    if length(lms)>1
        bs=lms;
    else
        bs=lms.bsb;
    end
    if init<length(bs)
        init=length(bs);
    end
    
    modelmdr=model;
    
    if isstruct(lms) && isfield(lms,'posLS')
        modelmdr.posLS=lms.posLS;
    end
    
    
    if isfield(modelmdr,'lshift')
        modelmdr=rmfield(modelmdr,'lshift');
    end
    
    % Compute Minimum Deletion Residual for each step of the search
    [mdr,Un,bb,Bols,S2,Exflag] = FSRtsmdr(y,bs,'model',modelmdr,'init',init,'plots',0,'nocheck',true,'msg',msg);
    
    if size(mdr,2)<2
        if length(mdr)>=T/2
            disp('More than half of the observations produce a singular X matrix')
            disp('X is badly defined')
            disp('If you wish to run the procedure using for updating the values of beta of the last step in which there was full rank use option bsbmfullrank=0')
            out.ListOut  = setdiff(seq,mdr);
        else
            disp('Bad starting point which produced a singular matrix, please restart the search from a different starting point or use option bsbmfullrank=0 ')
            
        end
        
        out.ListOut=NaN;
        out.outliers=NaN;
        out.mdr = NaN;
        out.Un  = NaN;
        out.nout= NaN;
        out.beta=NaN;
        out.scale=NaN;
        out.fittedvalues=NaN;
        out.residuals=NaN;
        out.Exflag=NaN;
        out.class='FSRts';
        return
    end
    
else % Initial subset is not supplied by the user
    
    % Find initial subset to initialize the search
    [out]=LTSts(y,'model',model,'h',h,'nsamp',nsamp,'msg',msg,'yxsave',1);
    
    if out.scale==0
        disp('More than half of the observations produce a linear model with a perfect fit')
        % Just return the outliers found by LXS
        %out.ListOut=out.outliers;
        %return
    end
    bs=out.bs;
    
    % bs=out.bs;
    mdr=0;
    constr='';
    
    % Prepare model for minimum deletion residual.
    %
    modelmdr=model;
    modelmdr.B=out.B;
    
    if isfield(modelmdr,'lshift')
        if model.lshift~=0
            modelmdr.posLS=out.posLS;
        end
        modelmdr=rmfield(modelmdr,'lshift');
    end
    
    while size(mdr,2)<2 && iter <6
        % Compute Minimum Deletion Residual for each step of the search
        % The instruction below is surely executed once.
        [mdr,Un,bb,Bols,S2,Exflag]=FSRtsmdr(y,bs,'model',modelmdr,'init',init,'plots',0,'nocheck',true,'msg',msg,'constr',constr,'bsbmfullrank',bsbmfullrank);
        
        % If FSRtsmdr runs without problems mdr has two columns. In the second
        % column it contains the value of the minimum deletion residual
        % monitored in each step of the search
        
        % If mdr has just one columns then one of the following two cases took place:
        % isnan(mdr)=1 ==> in this case initial subset was not full rank
        % mdr has just one column ==> in this case, even if the initial
        %    subset was full rank, the search has found at a certain step
        %    m<n/2 a list of units which produce a singular matrix. In this
        %    case LXS is rerun excluding these units which gave rise to a
        %    singular matrix
        
        if size(mdr,2)<2
            if length(mdr)>=T/2
                disp('More than half of the observations produce a singular X matrix')
                disp('If you wish to run the procedure using for updating the values of beta of the last step in which there was full rank use option bsbmfullrank=0')
                
                out.ListOut = setdiff(seq,mdr);
                
                return
            elseif isnan(mdr(1,1))
                % INITIAL SUBSET WAS NOT FULL RANK
                % restart LXS without the units forming
                % initial subset
                bsb=setdiff(seq,out.bs);
                [out]=LTSts(y(bsb),'model',model,'h',h,'nsamp',nsamp,'msg',msg,'yxsave',1);
                
                bs=bsb(out.bs);
                
                
            else
                % INITIAL SUBSET WAS FULL RANK BUT THE SEARCH HAS FOUND A
                % SET OF OBSERVATIONS CONSTR <n/2  WHICH PRODUCED A SINGULAR
                % MATRIX. IN THIS CASE NEW LXS IS BASED ON  n-constr OBSERVATIONS
                iter=iter+1;
                bsb=setdiff(seq,mdr);
                constr=mdr;
                % [out]=LXS(y(bsb),X(bsb,:),'lms',lms,'nsamp',nsamp,'nocheck',true,'msg',msg);
                [out]=LTSts(y(bsb),'model',model,'h',h,'nsamp',nsamp,'msg',msg,'yxsave',1);
                
                bs=bsb(out.bs);
            end
        end
    end
    
    
end


if iter >=5
    %     out.mdr = NaN;
    %     out.Un  = NaN;
    %     out.nout= NaN;
    error('FSDA:FSRts:NoConv','No convergence')
end

INP=struct;
INP.y=y;

if length(lms)>1 || (isstruct(lms) && isfield(lms,'bsb'))
    p=size(Bols,2)-1;
else
    INP.X=out.X;
    p=size(out.B(:,1),1);
end
INP.model=modelmdr;
INP.n=T;
INP.p=p;
INP.mdr=mdr;
INP.init=init;
INP.Un=Un;
INP.bb=bb;
INP.Bcoeff=Bols;
INP.S2=S2(:,1:2);

%% Call core function which computes exceedances to thresholds of mdr
[out]=FSRcore(INP,'ts',options);

% compute and store in output structure the robust scaled residuals
out.fittedvalues = out.yhat;
out.residuals    = (y-out.fittedvalues)/out.scale;
out=rmfield(out,'yhat');

out.class  =  'FSRts';
out.Exflag=Exflag;

end

%FScategory:REG-Regression