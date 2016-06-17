function out = FSMtra(Y,varargin)
%FSMtra computes MLE of transformation parameters.
%It uses the multivariate version of the parametric family of power
%transformations.
%
%
%<a href="matlab: docsearchFS('FSMtra')">Link to the help function</a>
%
% Required input arguments:
%
% Y :           Input data. Matrix.
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values are automatically excluded from the computations.
%               Data Types - single|double
%
% Optional input arguments:
%
%    family :   string which identifies the family of transformations which
%               must be used. Character. Possible values are 'BoxCox'
%               (default) or 'YJ'.
%               The Box-Cox family of power transformations equals
%               $(y^{\lambda}-1)/\lambda$ for $\lambda$ not equal to zero,
%               and $\log(y)$ if $\lambda = 0$.
%               The Yeo-Johnson (YJ) transformation is the Box-Cox
%               transformation of $y+1$ for nonnegative values, and of
%               $|y|+1$ with parameter 2-lambda for y negative.
%               The basic power transformation returns $y^{\lambda}$ if
%               $\lambda$ is not zero, and $\log(\lambda)$  otherwise.
%               Remark. BoxCox and the basic power family can be used just
%               if input y is positive. Yeo-Johnson family of
%               transformations does not have this limitation.
%               Example - 'family','YJ'
%               Data Types - char
%   init    :   Point where to start monitoring required diagnostics. Scalar.
%               Note that if bsb is suppliedinit>=length(bsb). If init is
%               not specified it will be set equal to floor(n*0.6).
%               Example - 'init',50
%               Data Types - double
%       bsb :   It contains the units forming initial subset. Vector. The
%               default value of bsb is '' that is the initial subset is
%               found through the intersection of robust bivariate ellipses
%               This option is useful if a forced start is required.
%               Example - 'bsb',[4 6 9]
%               Data Types - double
%        rf :   confidence level for bivariate ellipses. Scalar. Default is
%               0.9. If bsb is not empty this argument is ignored.
%               Example - 'rf',0.99
%               Data Types - double
%   ColToTra:   It specifies the variables which must be
%               transformed. Vector. It is a k x 1 integer vector.
%               Example - 'ColToTra',[1 3]
%               Data Types - double
%        la0:   It contains set of transformation
%               parameters for the k ColtoTra. Vector. It is a k x 1
%               vector.  The ordering of Mahalanobis distances at each step
%               of the forward search uses variables transformed with la0.
%               la0 empty is equivalent to its default value
%               la0=ones(length(ColToTra),1).
%               Example - 'la0',[-1 0]
%               Data Types - double
%  onelambda:   If onelambda=1, a common value lambda is estimated
%               for all variables specified in ColToTra. Scalar.
%               Example - 'onelambda',0
%               Data Types - double
%   optmin  :  It contains the options dealing with the
%               maximization algorithm. Structure. Use optimset to set
%               these options. Notice that the maximization algorithm which
%               is used is fminunc is the optimization toolbox is present
%               else is fminsearch.
%               Example -'optmin.Display','off'
%               Data Types - double
%     speed :   If speed=1 the initial value at step m of
%               the maximization procedure is the
%               final value at step m-1 else it is la0. Scalar. Default
%               value 1. The maximization procedure is fminunc or fminsearch.
%               Example -'speed',0
%               Data Types - double
%   colnames:   It contains the names of
%               the variables of the dataset. Cell array of strings. Cell
%               array of strings of length v. If colnames is empty then the
%               sequence 1:v is created to label the variables.
%               Example -'colnames', {'1' '2' '3' '4' '5' '10' '11' '12' '13'};
%               Data Types - cell array of strings
%   prolik  :   It specifies whether it is necessary to
%               monitor the profile log likelihood of the transformation
%               parameters at selected steps of the search. Scalar or
%               structure.
%               If prolik is a scalar, the plot of the profile loglikelihoods
%               is produced at step m=n with all default parameters.
%               If prolik is a structure it may contain the following fields:
%                   prolik.steps = vector containing the steps of the fwd
%                                  search for which profile logliks have to
%                                  be plotted. The default value of steps
%                                  is n;
%                   prolik.clev = scalar between 0 and 1 determining
%                                 confidence level for each element of
%                                 lambda based on the asymptotic chi1^2 of
%                                 twice the loglikelihood ratio.
%                                 The default confidence level is 0.95;
%                   prolik.xlim = vector with two elements determining
%                                 minimum and maximum values of lambda in
%                                 the plots of profile loglikelihoods. The
%                                 default value of xlim is [-2 2];
%                   prolik.LineWidth = line width of the vertical lines
%                                 defining confidence levels of the
%                                 transformation parameters.
%                 Example -'prolik',7
%                 Data Types - double
%   plotsmle:   It specifies whether it is necessary to
%               plot the maximum likelihood estimates of the transformation
%               parameters. Scalar or structure. Three horizontal lines
%               associated respectively with values -1, 0 and 1 are added
%               to the plot.
%               If prolik is a scalar, the plot of the monitoring of
%               maximum likelihood estimates of transformation parameters
%               is produced on the screen with all the default options.
%               If plotsmle is a structure, it may contain the following fields:
%                 plotsmle.xlim = minimum and maximum on the x axis;
%                 plotsmle.ylim = minimum and maximum on the y axis;
%                 plotsmle.LineWidth = Line width of the trajectories of
%                                      mle of transformation parameters;
%                 plotsmle.LineStyle = cell containing Line styles of
%                                      the trajectories of mle of
%                                      transformation parameters;
%                 plotsmle.LineWidthEnv=Line width of the horizontal
%                                      lines;
%                 plotsmle.Tag       = tag of the plot (default is pl_mle);
%                 plotsmle.FontSize  = font size of the text labels which
%                                      identify the trajectories.
%                 Example -'plotsmle',1
%                 Data Types - double
%   plotslrt:   It specifies whether it is necessary to
%               plot the likelihood ratio test. Scalar or structure.
%               If plotslrt is a scalar, the plot of the monitoring of
%               likelihood ratio test is produced on the screen with all
%               default options.
%               If plotslrt is a strucure, it may contain the following fields:
%                   plotslrt.xlim     = minimum and maximum on the x axis;
%                   plotslrt.ylim     = minimum and maximum on the y axis;
%                   plotslrt.LineWidth= Line width of the trajectory of lrt of
%                                       transformation parameters;
%                   plotslrt.conflev  = vector which defines the confidence
%                                       levels of the horizontal line for
%                                       the likelihood ratio test (default
%                                       is conflev=[0.95 0.99]);
%                   plotslrt.LineWidthEnv= Line width of the horizontal lines;
%                   plotslrt.Tag      = tag of the plot (default is pl_lrt).
%                 Example -'plotslrt',1
%                 Data Types - double
%  msg  :      It controls whether to display or not messages
%               about great interchange on the screen.  Scalar.
%               If msg==1 (default) messages are displayed on the screen
%               else no message is displayed on the screen.
%                 Example -'msg',1
%                 Data Types - double
%
% Remark:       The user should only give the input arguments that have to
%               change their default value.
%               The name of the input arguments needs to be followed by
%               their value. The order of the input arguments is of no
%               importance.
%
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations. y can be both a row of column vector.
%
% Output:
%
%         out:   structure which contains the following fields
%
%   out.MLEtra= n-init+1 x v matrix which contains the monitoring of
%               MLE of transformation parameters:
%               1st col = fwd search index (from init to n);
%               2nd col = MLE of variable 1;
%               3rd col = MLE of variable 2;
%               ...;
%               (v+1)th col = MLE of variable v.
%   out.LIKrat= n-init+1 x 2 = matrix which contains the monitoring of
%               likelihood ratio for testing H0:\lambda=la0:
%               1st col = fwd search index (from init to n);
%               2nd col = value of the likelihood ratio.
%   out.Exflag= n-init+1 x 2 = matrix which contains the monitoring of
%               the integer identifying the reason why the maximization
%               algorithm terminated. See help page fminunc of the
%               optimization toolbox  for the list of values of exitflag
%               and the corresponding reasons the algorithm terminated:
%               1st col = fwd search index (from init to n);
%               2nd col = the value that describes the exit condition
%   out.Un =    (n-init) x 11 Matrix which contains the unit(s)
%               included in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Un(1,2) for example contains
%               the unit included in step init+1 Un(end,2) contains the
%               units included in the final step of the search
%
% See also FSMmmd.m, FSM.m
%
% References:
%
%   Atkinson Riani and Cerioli (2004), Exploring multivariate data with the
%   forward search Springer Verlag, New York.
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSMtra')">Link to the help function</a>
% Last modified 31-05-2016

% Examples:
%{
    % FSMtra with all default options.
    % Baby food data.
    load('baby.mat');
    Y=baby.data(:,6:end);
    % FS based on untrasnformed data H_0:\lambda=1 for all variables
    % Plot of mle of transformation parameters with all default options
    % Compare the output with Figure 4.7 p. 167 of ARC (2004)
    [out]=FSMtra(Y);
%}

%{
    %% FSMtra with optional arguments.
    % Plot the maximum likelihood estimates of the transformation parameters.
    % Baby food data.
    load('baby.mat');
    Y=baby.data(:,6:end);
    % FS based on untrasnformed data H_0:\lambda=1 for all variables
    % Plot of mle of transformation parameters with all default options
    % Compare the output with Figure 4.7 p. 167 of ARC (2004)
    [out]=FSMtra(Y,'plotsmle',1);
%}


%{
    %% Personalized options for plotsmle.
    % Baby food data.
    load('baby.mat');
    Y=baby.data(:,6:end);
    % Personalized options for plotsmle
    plotsmle=struct;
    plotsmle.LineWidth=3;
    plotsmle.LineWidthEnv=3;
    plotsmle.FontSize=14;
    plotsmle.ylim=[-0.4 0.4];
    [out]=FSMtra(Y,'plotsmle',plotsmle);
%}

%{
    % FSMtra based on log trasnformed data.
    % Baby food data.
    load('baby.mat');
    Y=baby.data(:,6:end);
    % FS based on log trasnformed data H_0:\lambda=0 for all variables
    % Plot of mle of transformation parameters with all default options
    v=size(Y,2);
    plotsmle=struct;
    plotsmle.ylim=[-0.4 1];
    [out]=FSMtra(Y,'la0',zeros(v,1),'init',11,'plotsmle',plotsmle);
%}

%{
    %% Monitoring the profile log likelihood of the transformation parameters.
    % Baby food data ignoring the regression structure.
    load('baby.mat');
    Y=baby.data(:,6:end);
    % FS based on log trasnformed data H_0:\lambda=0 for all variables
    % Plot of mle of transformation parameters with all default options
    v=size(Y,2);
    plotsmle=struct;
    plotsmle.ylim=[-0.4 1];
    prolik=struct;
    prolik.steps=26;
    prolik.xlim=[-1 1];
    [out]=FSMtra(Y,'la0',zeros(v,1),'init',11,'prolik',prolik);
%}

%{
    % Swiss heads data, Example 1.
    % FSMtra based on untransformed data $H_0:\lambda=1$ for all variables
    load('head.mat');
    Y=head.data;
    [out]=FSMtra(Y,'plotsmle',1);
%}

%{
    % Swiss heads data, Example 2.
    % FSMtra based on untransformed data $H_0:\lambda=1$ for all variables
    load('head.mat');
    Y=head.data;
    % Analysis of profile loglikelihoods at step m=198
    prolik=struct;
    prolik.steps=198;
    prolik.xlim=[-3 5];
   [out]=FSMtra(Y,'prolik',prolik);
%}

%{
    % Swiss heads data, Example 3.
    % FSMtra based on untransformed data $H_0:\lambda=1$ for all variables
    load('head.mat');
    Y=head.data;
    % Monitoring of likelihood ratio test
    % Compare the output with Figure 4.13 p. 172 of ARC (2004)
    [out]=FSMtra(Y,'plotslrt',1);
%}

%{
    % Swiss heads data, Example 4.
    load('head.mat');
    Y=head.data;
    % FS based on untransformed data H_0:\lambda=1 for variable 4
    % Monitoring of likelihood ratio test
    % Compare the output with Figure 4.14 p. 173 of ARC (2004)
    [out]=FSMtra(Y,'ColToTra',4,'plotslrt',1);
%}

%{
    % Mussels data, Example 1.
    % FSMtra based on untransformed data $H_0:\lambda=1$ for all variables
    load('mussels.mat');
    Y=mussels.data;
    % Compare plot of mle with Figure 4.19 p. 178 of ARC (2004)
    % Compare plot of lrt with Figure 4.18 p. 178 of ARC (2004)
    [out]=FSMtra(Y,'plotsmle',1,'plotslrt',1);
%}

%{
    % Mussels data, Example 2.
    load('mussels.mat');
    Y=mussels.data;
    % FSMtra based on with H_0:\lambda=[1 0.5 1 0 1/3]
    % Compare plot of mle with Figure 4.21 p. 178 of ARC (2004)
    % Compare plot of lrt with Figure 4.20 p. 178 of ARC (2004)
    [out]=FSMtra(Y,'la0',[1 0.5 1 0 1/3],'plotsmle',1,'plotslrt',1);
%}


%{
    % Swiss bank notes data, Example 1.
    load('swiss_banknotes')
    Y=swiss_banknotes.data;
    n=size(Y,1);
    Y1=repmat(max(Y),n,1);
    Y=Y./Y1;
    % FS using just one value of lambda for all the variables
    % Compare plot of lrt with left panel of Figure 4.69 p. 225 of ARC (2004)
    [out]=FSMtra(Y,'init',40,'onelambda',1,'plotslrt',1);
%}

%{
    % Swiss bank notes data, Example 2.
    load('swiss_banknotes')
    Y=swiss_banknotes.data;
    n=size(Y,1);
    Y1=repmat(max(Y),n,1);
    Y=Y./Y1;
    % FS using just one value of lambda for all the variables
    % Search starts with the first 20 genuine notes
    % Compare plot of lrt with central panel of Figure 4.69 p. 225 of ARC (2004)
    [out]=FSMtra(Y,'init',20,'onelambda',1,'bsb',1:20,'plotslrt',1);
%}

%{
    % Swiss bank notes data, Example 3.
    load('swiss_banknotes')
    Y=swiss_banknotes.data;
    Y=Y(1:100,:);
    % FS using just one value of lambda for all the variables
    % Monitoring of mle of lambda (Figure 4.66 p. 223 of ARC (2004))
    % Monitoring of lrt (Figure 4.67 p. 223 of ARC (2004))
    plotsmle=struct;
    plotsmle.ylim=[-1.5 2.5];
    % Profile loglikelihoods at steps m=90 and m=100
    % (Figure 4.68 p. 224 of ARC (2004))
    prolik=struct;
    prolik.steps=[90 100];
    prolik.xlim=[-3.2 3.2];
    [out]=FSMtra(Y,'onelambda',1,'plotsmle',plotsmle,'plotslrt',1,'prolik',prolik);
%}

%{
    % Swiss bank notes data, Example 4.
    load('swiss_banknotes')
    Y=swiss_banknotes.data;
    n=size(Y,1);
    Y1=repmat(max(Y),n,1);
    Y=Y./Y1;
    % FS using just one value of lambda for all the variables
    % Search starts with the first 20 forged notes
    % Compare plot of lrt with right panel of Figure 4.69 p. 225 of ARC (2004)
    [out]=FSMtra(Y,'init',20,'onelambda',1,'bsb',101:120,'plotslrt',1);
%}

%{

    %% Track records data.
    load('recordfg');
    Y=recordfg.data;
    n=size(Y,1);
    Y1=repmat(max(Y),n,1);
    Y=Y./Y1;
    la0=[-1 -2 -3 -4];
    tags={'lrt-1' 'lrt-2' 'lrt-3' 'lrt-4'};
    plotslrt=struct;
    plotslrt.ylim=[0 21];
    ii=1;
    for la=la0;
        plotslrt.Tag=tags{ii};
        [out]=FSMtra(Y,'plotslrt',plotslrt,'onelambda',1,'la0',la);
        ii=ii+1;
    end
    disp(out)
% Compare these 4 plots with Figure 4.50 p. 207 of ARC (2004)
%}

%{
    % Emilia Romagna data, Example 1.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Extract demographic variables
    Y1=Y(:,[1 2 3 4 5 10 11 12 13]);
    colnames={'1' '2' '3' '4' '5' '10' '11' '12' '13'};
    [out]=FSMtra(Y1,'plotsmle',1,'colnames',colnames);
    % Compare the plot with Figure 4.31 p. 188 of ARC (2004)

%}

%{
    % Emilia Romagna data, Example 2.
    % Yeo and Johnson family of transformations is
    % used. In this case there is no need to correct 0 values
    load('emilia2001')
    Y=emilia2001.data;

    % Extract demographic variables
    Y1=Y(:,[1 2 3 4 5 10 11 12 13]);
    colnames={'1' '2' '3' '4' '5' '10' '11' '12' '13'};
    plotslrt=struct;
    plotslrt.ylim=[-8.2 8.2];
    ColToComp=[1 3 5 9];
    la0=[0 0.25 0 0.5 0.5 0 0 0.5 0.25];
    [out]=FSMfan(Y1,la0,'ColToComp',ColToComp,'plotslrt',plotslrt,'colnames',colnames,'family','YJ');
%}

%{
    % Emilia Romagna data, Example 3.
    % Demographic data
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Extract demographic variables
    Y1=Y(:,[1 2 3 4 5 10 11 12 13]);
    colnames={'1' '2' '3' '4' '5' '10' '11' '12' '13'};
    la0=[0 0.25 0 0.5 0.5 0 0 0.5 0.25];
    prolik=struct;
    prolik.steps=[331];
    prolik.xlim=[-1 1];

    plotslrt=struct;
    plotslrt.ylim=[4 21];
 
   [out]=FSMtra(Y1,'plotsmle',1,'plotslrt',plotslrt,'la0',la0,'colnames',colnames,'prolik',prolik);
    % Compare the plots with Figures 4.32, 4.33 and 4.34 p. 189-191 of ARC (2004)

%}

%{
    % Emilia Romagna data, Example 4.
    % Modified wealth variables.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify wealth variables
    Y(:,16)=100-Y(:,16);
    Y(:,23)=100-Y(:,23);
    % Extract wealth variables
    Y1=Y(:,[14:23]);
    colnames={'14' '15' '16' '17' '18' '19' '20' '21' '22' '23'};
    la0=[0 1 0.25 1 1 0.5 -0.5 0.25 0.25 -1];
    [out]=FSMtra(Y1,'plotslrt',1,'la0',la0,'colnames',colnames);
    % Compare the plot with left panel of Figure 4.38 p. 188 of ARC (2004)

%}

%{
    % Emilia Romagna data, Example 5.
    % Modified wealth variables.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify wealth variables
    Y(:,16)=100-Y(:,16);
    Y(:,23)=100-Y(:,23);
    % Extract wealth variables
    Y1=Y(:,[14:23]);
    colnames={'14' '15' '16' '17' '18' '19' '20' '21' '22' '23'};
    la0=[0.5 1 0.25 1 1 0.5 -0.5 0.25 0.25 -1];
    [out]=FSMtra(Y1,'plotslrt',1,'la0',la0,'colnames',colnames);
    % Compare the plot with Figure 4.40 p. 196 of ARC (2004)

%}

%{
    % Emilia Romagna data, Example 6.
    % Work variables.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21 25 26];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Extract work variables
    Y1=Y(:,[6:9 24:28]);
    colnames={'6' '7' '8' '9' '24' '25' '26' '27' '28'};
    la0=[0.25,0,2,-1,0,1.5,0.5,1,1];
    [out]=FSMtra(Y1,'plotsmle',1,'plotslrt',1,'la0',la0,'colnames',colnames);
    % Compare the plots with Figures 4.41 p. 197 and left panel of Figure
    % 4.42 of ARC (2004)
%}

%{
    % Emilia Romagna data, Example 7.
    % Modified work variables.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify variables 25 and 26
    Y(:,25)=100-Y(:,25);
    Y(:,26)=100-Y(:,26);
    % Extract work variables
    Y1=Y(:,[6:9 24:28]);
    colnames={'6' '7' '8' '9' '24' '25' '26' '27' '28'};
    la0=[0.25,0,2,-1,0,0,1.5,1,1];
    [out]=FSMtra(Y1,'plotsmle',1,'plotslrt',1,'la0',la0,'colnames',colnames);
%}


%% Input parameters checking
% Extract size of the data
[n,v]=size(Y);
% Initialize matrix which will contain Mahalanobis distances in each step
seq=(1:n)';
one=ones(n,1);

if nargin<1
    error('FSDA:FSMtra:missingInputs','Initial data matrix is missing');
end
family='BoxCox';
hdef=floor(n*0.6);
options=struct('init',hdef,'bsb','','ColToTra','','la0','','onelambda',0,'rf',0.9,...
    'speed',1,'optmin',optimset,'msg',1,'colnames','',...
    'prolik','','plotsmle','','plotslrt','','family',family);


UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSMtra:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

%init1=options.init;
if nargin > 1
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
    family=options.family;
end


Xnormfamily=strcat('norm',family);
familytra=str2func(Xnormfamily);

% And check if the optional user parameters are reasonable.
% check init
init1=options.init;
msg=options.msg;

if  init1 <v+1;
    mess=sprintf(['Attention : init1 should be larger than v. \n',...
        'It is set to v+1.']);
    disp(mess);
    init1=v+1;
elseif init1>n;
    mess=sprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n.']);
    disp(mess);
    init1=n;
end

% ColToTra = vector which contains the numbers associated to the variables
% which have to be transformed
ColToTra=options.ColToTra;
if isempty(ColToTra)
    ColToTra=(1:v);
elseif max(ColToTra)>v || min(ColToTra)<1
    error('FSDA:FSMtra:WorngColtToTra',['The columns to transform are not in the range. 1-' num2str(v)]);
end

ctra=length(ColToTra);

% la0= values of transformation parameter(s) which must be tested
la0=options.la0;

% Specify whether a common value of lambda is required
onelambda=options.onelambda;

% Initialize matrix which will contain the MLE of the transformation parameters
% monitored in each step
if onelambda==1
    MLEtra=cat(2,(init1:n)',NaN(n-init1+1,1));
    if strcmp(family,'BoxCox')
        loglik=@likonelambdaBoxCox;
    elseif strcmp(family,'YJ')
        loglik=@likonelambdaYJ;
    else
        warning('FSDA:FSMtra:WrongFamily','Transformation family which has been chosen is not supported')
        error('FSDA:FSMtra:WrongFamily','Supported values are BoxCox or YaoJohnson or basicpower')
    end
    
    if isempty(la0)
        la0=1;
    end
else
    if isempty(la0)
        la0=ones(1,ctra);
    end
    MLEtra=cat(2,(init1:n)',NaN(n-init1+1,v));
    
    if strcmp(family,'BoxCox')
        loglik=@likBoxCox;
    elseif strcmp(family,'YJ')
        loglik=@likYJ;
    else
        warning('FSDA:FSMtra:WrongFamily','Transformation family which has been chosen is not supported')
        error('FSDA:FSMtra:WrongFamily','Supported values are BoxCox or YaoJohnson')
    end
end

% lainit starting value for minimization at step m=init
lainit=la0;
laout=la0;

% if speed=1 starting value for minimization at step m>init is final value at
% step m-1 else vector la0 is used as starting value in each step
speed=options.speed;

% Options to use in the minimization;
optmin=options.optmin;
optmin.Display='off';
% optmin=optimset('Display','off');
% optmin=optimset('Display','on');
% optmin=optimset(defopt,newopts);
warning('off','optim:fminunc:SwitchingMethod');
%warning('on','optim:fminunc:SwitchingMethod');

% Check if minimization toolbox is installed in current computer
typemin=exist('fminunc','file');

% Plot of profile loglikelihood
prolik=options.prolik;


if length(la0) ~= length(ColToTra) && onelambda~=1
    error('FSDA:FSMtra:WrongLambda','Length of vector lambda must be equal the length of the vector which contains the columns to transform');
elseif onelambda==1 && length(la0)>1
    error('FSDA:FSMtra:WrongLambda','onelambda=1 therefore lambda must a scalar');
end

if onelambda==1
    % Ytr contains the matrix of transformed values.
    Ytr=feval(familytra,Y,ColToTra,la0*ones(length(ColToTra),1));
else
    Ytr=feval(familytra,Y,ColToTra,la0);
end

if ~isempty(options.colnames)
    colnames=options.colnames;
else
    colnames=cellstr(num2str((1:v)'));
end


if ~isempty(prolik)
    if ctra<=4;
        nr=2;
        nc=2;
    elseif ctra<=6
        nr=2;
        nc=3;
    elseif ctra<=9
        nr=3;
        nc=3;
    elseif ctra<=16
        nr=4;
        nc=4;
    elseif ctra<=25
        nr=5;
        nc=5;
    else
        error('FSDA:FSMtra:TooManyVars','Profile loglikelihood cannot be displayed for more than 25 variables')
    end
    
    if isstruct(prolik)
        % Check if confidence level has been specified by the
        % user
        
        fprolik=fieldnames(prolik);
        
        % specify steps where it is necessary to plot the
        % profile loglikelihoods
        d=find(strcmp('steps',fprolik));
        if d>0
            steps=prolik.steps;
        else
            steps=n;
        end
        
        d=find(strcmp('clev',fprolik));
        if d>0
            clev=chi2inv(prolik.clev,1);
        else
            clev=chi2inv(0.95,1);
        end
        
        d=find(strcmp('xlim',fprolik));
        if d>0
            xlimla=prolik.xlim;
        else
            xlimla=[-2 2];
        end
        
    else
        steps=n;
        clev=chi2inv(0.95,1);
        xlimla=[-2 2];
        
    end
    
    % seqla = vector which contains the xlimits for lambda in the profile
    % loglikelihood plots
    seqla=xlimla(1):0.01:xlimla(2);
    
end



% Confidence level for robust bivariate ellipses
rf=options.rf;

% Find initial subset to initialize the search
if isempty(options.bsb)
    [fre]=unibiv(Ytr,'rf',rf);
    fre=sortrows(fre,[3 4 2]);
    bsb=fre(1:v+1,1);
    ini0=v+1;
else
    bsb=options.bsb;
    ini0=length(bsb);
end

% Initialize matrix which stores in each step the integer identifying the
% reason why the algorithm terminated
Exflag=MLEtra(:,1:2);

% Initialize matrix which will contain the likelihood ratio for
% \lambda=\lambda0 monitored in each step
LIKrat=MLEtra(:,1:2);

%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2,(init1+1:n)',NaN(n-init1,10));


% Initialize the matrix which contains the values of MD based on
% transformed obervations Ytr
mala=[seq zeros(n,1)];


if (rank(Ytr(bsb,:))<v)
    warning('FSDA:FSMtra:NoFullRank','The supplied initial subset is not full rank matrix');
    disp('FS loop will not be performed')
    out=struct;
else
    
    for mm = ini0:n
        
        % Ytrb =subset of transformed values
        Ytrb=Ytr(bsb,:);
        
        % Yb = subset of original values
        Yb=Y(bsb,:);
        
        % Find vector of means of subset of transformed values
        % Note that ym is a row vector
        ym=mean(Ytrb);
        
        
        % Squared Mahalanobis distances computed using QR decomposition
        % based on transformed values
        Ym=Ytr-one*ym;
        [~,R]=qr(Ym(bsb,:),0);
        
        
        % Remark: u=(Ym/R)' should be much faster than u=inv(R')*Ym';
        u=(Ym/R)';
        % Compute square Mahalanobis distances on transformed data
        mala(:,2)=sqrt(((mm-1)*sum(u.^2)))';
        
        
        if (mm>=init1);
            
            if speed==1;
                lainit=laout+1e-4*randn(size(lainit));
            end
            
            
            try
                fval='';
                
                % In order to find minimum it is possible to use MATLAB
                % function fminsearch or function fminunc from the optimization
                % toolbox
                if typemin==2;
                    [laout,fval,exitflag]  = fminunc(loglik,lainit,optmin);
                else
                    [laout,fval,exitflag]  = fminsearch(loglik,lainit,optmin);
                end
                %laout=5*sin(laout);
            catch
                disp(['Warning: non convergence at step mm=' num2str(mm)])
                if isempty(fval)
                    exitflag=nan;
                    laout=nan;
                    fval=nan;
                end
            end
            % Store maximum likelihood estimates of transformation parameters
            Exflag(mm-init1+1,2)=exitflag;
            
            % Store value of the likelihood ratio for H0
            if onelambda==1;
                % Store common maximum likelihood estimate of tranformation parameter
                MLEtra(mm-init1+1,2)=laout;
                Ytrb0=feval(familytra,Yb,ColToTra,la0*ones(length(ColToTra),1));
            else
                % Store maximum likelihood estimates of tranformation parameters
                MLEtra(mm-init1+1,1+ColToTra)=laout;
                Ytrb0=feval(familytra,Yb,ColToTra,la0);
            end
            
            LIKrat(mm-init1+1,2)=mm*(log(det(cov(Ytrb0)))-fval);
        end
        
        
        % Create plots of profile loglikelihood at selected steps of the fwd search
        if ~isempty(prolik)
            
            if ~isempty(intersect(mm,steps))
                % -fval*mm = value of the profile loglik in correspondence
                % of the maximum
                figure;
                
                pmax=-fval*mm-0.5*clev;
                
                if onelambda==1
                    
                    Loglikj=[seqla', zeros(length(seqla),1)];
                    % Loop over lambda (values defined in seqla)
                    il=1;
                    for lai=seqla
                        Ybpro=feval(familytra,Yb,ColToTra,lai*ones(ctra,1));
                        Loglikj(il,2)=-mm*log(det(cov(Ybpro)));
                        il=il+1;
                    end
                    plot(Loglikj(:,1),Loglikj(:,2))
                    xlabel('\lambda')
                    xlim([seqla(1) seqla(end)])
                    
                    title(['Profile loglikelihood: unique value of transformation parameter  at step m='  num2str(mm)],'FontSize',12);
                    % Construct (1-alpha)% confidence level. The default is to
                    % construct 95% confidence interval, that is find the
                    % points where the profile loglikelihood curve has
                    % decreased by 3.84/2
                    
                    v=axis;
                    
                    selinf=Loglikj(Loglikj(:,1)<laout,:);
                    admiss=selinf(selinf(:,2)<pmax,:);
                    if ~isempty(admiss)
                        cflow=admiss(end,:);
                        % line in correspondence of lower confidence interval
                        line([cflow(1,1);cflow(1,1)],[v(3);cflow(1,2)],'color','r','LineWidth',0.1,'Tag','env');
                    end
                    
                    selsup=Loglikj(Loglikj(:,1)>laout,:);
                    admiss=selsup(selsup(:,2)<pmax,:);
                    if ~isempty(admiss)
                        cfsup=admiss(1,:);
                        % line in correspondence of upper confidence interval
                        line([cfsup(1,1);cfsup(1,1)],[v(3);cfsup(1,2)],'color','r','LineWidth',0.1,'Tag','env');
                    end
                    
                    % line in correspondence of MLE of transformation
                    % parameter
                    line([laout;laout],[v(3);-fval*mm],'color','r','LineWidth',0.1,'Tag','env');
                    
                else
                    
                    % Ymax contains normalized Box Cox transformed values of Y in
                    % correspondence of maximum likelihood estimate of
                    % transformation parameters
                    Ybmax=feval(familytra,Yb,ColToTra,laout);
                    
                    % loop over the variables which have to be transformed
                    ij=1;
                    for jj=1:ctra
                        Loglikj=[seqla', zeros(length(seqla),1)];
                        cj=ColToTra(jj);
                        
                        % Loop over lambda
                        il=1;
                        for lai=seqla
                            Ybpro=Ybmax;
                            Ybpro(:,cj)=Yb(:,cj);
                            Ybpro=feval(familytra,Ybpro,cj,lai);
                            Loglikj(il,2)=-mm*log(det(cov(Ybpro)));
                            il=il+1;
                        end
                        subplot(nr,nc,ij)
                        plot(Loglikj(:,1),Loglikj(:,2))
                        ylabel(['Var. ' colnames{cj}])
                        xlabel('\lambda')
                        xlim([seqla(1) seqla(end)])
                        % Add the (1-alpha)% confidence interval
                        v=axis;
                        
                        
                        if ij==1;
                            %curax=gca;
                            yl=get(gca,'Ylim');
                        else
                            % curax=gca;
                            set(gca,'Ylim',yl);
                        end
                        ylim([yl(1) yl(2)])
                        
                        ij=ij+1;
                        
                        selinf=Loglikj(Loglikj(:,1)<laout(jj),:);
                        admiss=selinf(selinf(:,2)<pmax,:);
                        if ~isempty(admiss)
                            cflow=admiss(end,:);
                            % line in correspondence of lower confidence interval
                            % line([cflow(1,1);cflow(1,1)],[v(3);cflow(1,2)],'color','r','LineWidth',0.1,'Tag','env');
                            line([cflow(1);cflow(1,1)],[yl(1);cflow(1,2)],'color','r','LineWidth',0.1,'Tag','env');
                        end
                        
                        selsup=Loglikj(Loglikj(:,1)>laout(jj),:);
                        admiss=selsup(selsup(:,2)<pmax,:);
                        if ~isempty(admiss)
                            cfsup=admiss(1,:);
                            % line in correspondence of upper confidence interval
                            % line([cfsup(1,1);cfsup(1,1)],[v(3);cfsup(1,2)],'color','r','LineWidth',0.1,'Tag','env');
                            line([cfsup(1,1);cfsup(1,1)],[yl(1);cfsup(1,2)],'color','r','LineWidth',0.1,'Tag','env');
                        end
                        
                        % line in correspondence of MLE of transformation
                        % parameter. Draw it only if it is inside the
                        % xlimits
                        if laout(jj)>v(1) && laout(jj)<v(2)
                            % line([laout(jj);laout(jj)],[v(3);-fval*mm],'color','r','LineWidth',0.1,'Tag','env');
                            line([laout(jj);laout(jj)],[yl(1);-fval*mm],'color','r','LineWidth',0.1,'Tag','env');
                        end
                    end
                    % Now add the general title to the set of subplots
                    ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
                    set(get(ax,'Title'),'Visible','on')
                    title(['Profile loglikelihoods for the ' num2str(ctra) ' transformation parameters  at step m='  num2str(mm)],'FontSize',12);
                end
                set(gcf,'Name',['Profile loglikelihoods of transformation parameters at step m=' num2str(mm)]);
                
            end
        end
        
        
        zs=sortrows(mala,2);
        if mm<n
            % eval('mm');
            
            % store units forming old subset in vector oldbsb
            oldbsb=bsb;
            
            % the dimension of subset increases by one unit.
            % vector bsb contains the indexes corresponding to the units of
            % the new subset
            bsb=zs(1:mm+1,1);
            
            if (mm>=init1);
                unit=setdiff(bsb,oldbsb);
                if (length(unit)<=10)
                    Un(mm-init1+1,2:(length(unit)+1))=unit;
                else
                    if msg==1
                        disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        disp(['Number of units which entered=' int2str(length(unit))]);
                    end
                    Un(mm-init1+1,2:end)=unit(1:10);
                end
            end
        end
    end % close FS loop
    
    lla=length(laout);
    
    %% Plot of monitoring MLE of transformation parameters
    if ~isempty(options.plotsmle)
        
        plotsmle=options.plotsmle;
        
        if isstruct(plotsmle)
            
            fplotsmle=fieldnames(plotsmle);
            
            d=find(strcmp('xlim',fplotsmle));
            if d>0
                xlimx=plotsmle.xlim;
            else
                xlimx='';
            end
            
            d=find(strcmp('ylim',fplotsmle));
            if d>0
                ylimy=plotsmle.ylim;
            else
                ylimy='';
            end
            
            d=find(strcmp('LineWidth',fplotsmle));
            if d>0
                LineWidth=plotsmle.LineWidth;
            else
                LineWidth=2;
            end
            
            % LineWidthEnv = line width of the trajectories associated with -1 0 1 in
            % the plot of monitoring of MLE of transformation parameters or the
            % horizontal lines associated with the asymptotic confidence levels
            % for the likelihood ratio test
            
            d=find(strcmp('LineWidthEnv',fplotsmle));
            if d>0
                LineWidthEnv=plotsmle.LineWidthEnv;
            else
                LineWidthEnv=1;
            end
            
            % Specify the line type for the trajectories of MLE of
            % transformation parameters
            d=find(strcmp('LineStyle',fplotsmle));
            if d>0
                LineStyle=plotsmle.LineStyle;
            else
                slin=repmat({'-';'--';':';'-.'},ceil(lla/4),1);
                LineStyle=slin(1:lla);
            end
            
            
            
            d=find(strcmp('Tag',fplotsmle));
            if d>0
                tag=plotsmle.Tag;
            else
                tag='pl_mle';
            end
            
            % Font size for text messages associated with the trajectories of
            % MLE of trasnformation parameters
            
            d=find(strcmp('FontSize',fplotsmle));
            if d>0
                FontSize=plotsmle.FontSize;
            else
                FontSize=12;
            end
        else
            
            xlimx='';
            ylimy='';
            LineWidth=2;
            LineWidthEnv=1;
            tag='pl_mle';
            FontSize=12;
            slin=repmat({'-';'--';':';'-.'},ceil(lla/4),1);
            LineStyle=slin(1:lla);
        end
        
        % Specify where to send the output of the monitoring of MLE of
        % transformation parameters
        hmle=findobj('-depth',1,'tag',tag);
        if (~isempty(hmle))
            clf(hmle);
            figure(hmle)
            axes;
        else
            figure;
            set(gcf,'Name','MLE of transformation parameters');
        end
        
        
        
        if onelambda==1
            plot(MLEtra(:,1),MLEtra(:,2),'LineWidth',LineWidth);
        else
            plot1=plot(MLEtra(:,1),MLEtra(:,ColToTra+1),'LineWidth',LineWidth);
            % Add labels at the end of the search
            text(n*ones(1,lla),MLEtra(end,ColToTra+1),colnames(ColToTra),'FontSize',FontSize,'HorizontalAlignment','Left');
            % Add labels at the beginning of the search
            text(MLEtra(1,1)*ones(1,lla),MLEtra(1,ColToTra+1),colnames(ColToTra),'FontSize',FontSize,'HorizontalAlignment','Right');
            
            % Specify the line type for the trajectories of MLE of
            % transformation parameters
            set(plot1,{'LineStyle'},LineStyle);
        end
        
        set(gcf,'Tag',tag)
        
        if ~isempty(xlimx)
            xlim(xlimx);
        end
        
        if ~isempty(ylimy)
            ylim(ylimy);
        end
        
        v=axis;
        if max(max(MLEtra(:,2:end)))>1;
            line([v(1),v(2)],[1,1],'color','r','LineWidth',LineWidthEnv,'Tag','env');
        end
        minMLE=min(min(MLEtra(:,2:end)));
        if minMLE<-1;
            line([v(1),v(2)],[-1,-1],'color','r','LineWidth',LineWidthEnv,'Tag','env');
        elseif minMLE<0;
            line([v(1),v(2)],[0,0],'color','r','LineWidth',LineWidthEnv,'Tag','env');
        end
        
        xlabel('Subset size m');
        ylabel(['MLE of transformation parameters H_0: \lambda=' mat2str(la0) ]);
    end
    
    
    %% Plot of likelihood ratio test
    plotslrt=options.plotslrt;
    if ~isempty(plotslrt)
        
        if isstruct(plotslrt)
            
            fplotslrt=fieldnames(plotslrt);
            
            d=find(strcmp('xlim',fplotslrt));
            if d>0
                xlimx=plotslrt.xlim;
            else
                xlimx='';
            end
            
            d=find(strcmp('ylim',fplotslrt));
            if d>0
                ylimy=plotslrt.ylim;
            else
                ylimy='';
            end
            
            d=find(strcmp('LineWidth',fplotslrt));
            if d>0
                LineWidth=plotslrt.LineWidth;
            else
                LineWidth=2;
            end
            
            % LineWidthEnv = line width of the
            % horizontal lines associated with the asymptotic confidence levels
            % for the likelihood ratio test
            
            d=find(strcmp('LineWidthEnv',fplotslrt));
            if d>0
                LineWidthEnv=plotslrt.LineWidthEnv;
            else
                LineWidthEnv=1;
            end
            
            d=find(strcmp('conflev',fplotslrt));
            if d>0
                conflev=plotslrt.conflev;
            else
                conflev=[0.95 0.99];
            end
            
            d=find(strcmp('Tag',fplotslrt));
            if d>0
                tag=plotslrt.Tag;
            else
                tag='pl_lrt';
            end
        else
            
            xlimx='';
            ylimy='';
            LineWidth=2;
            LineWidthEnv=1;
            tag='pl_lrt';
            conflev=[0.95 0.99];
        end
        
        % Plot of the likelihood ratio
        % Specify where to send the output of the monitoring of lrt
        hlrt=findobj('-depth',1,'tag',tag);
        if (~isempty(hlrt))
            clf(hlrt);
            figure(hlrt)
            axes;
        else
            figure;
        end
        set(gcf,'Name',['Likelihood ratio $\lambda$=' mat2str(la0)]);
        
        
        plot(LIKrat(:,1),LIKrat(:,2),'LineWidth',LineWidth);
        
        set(gcf,'Tag',tag)
        
        if ~isempty(xlimx)
            xlim(xlimx);
        end
        
        if ~isempty(ylimy)
            ylim(ylimy);
        end
        
        
        
        xlabel('Subset size m');
        ylabel(['Lik ratio test for \lambda=' mat2str(la0) ]);
        
        if onelambda==1
            quant=chi2inv(conflev,1);
        else
            quant=chi2inv(conflev,length(ColToTra));
        end
        v=axis;
        line([v(1),v(2)],quant(1)*[1,1],'color','r','LineWidth',LineWidthEnv,'Tag','env');
        line([v(1),v(2)],quant(2)*[1,1],'color','r','LineWidth',LineWidthEnv,'Tag','env');
        
    end
    
    
    %% Store quantities in structure out
    % MLEtra=MLE of tramsformation parameters
    out.MLEtra=MLEtra;
    % Exflag = reason why maximization algorithm stopped
    out.Exflag=Exflag;
    % LIKrat = Likelihood ratio test
    out.LIKrat=LIKrat;
    % Un = Units entering the subset
    out.Un=Un;
    
end


% lik computes the likelihood when different lambdas are possible for
% different variables
    function dZ=likBoxCox(la)
        Z=Yb;
        
        for j=1:ctra;
            laj=la(j);
            cj=ColToTra(j);
            Gj=exp(mean(log(Yb(:,cj))));
            
            if laj~=0
                Z(:,cj)=(Yb(:,cj).^laj-1)/(laj*(Gj^(laj-1)));
            else
                Z(:,cj)=Gj*log(Yb(:,cj));
            end
        end
        
        
        dZ=log(det(cov(Z)));
        % disp(dZ);
    end


% lik computes the likelihood when just one value of lambda is allowed
    function dZ=likonelambdaBoxCox(la)
        Z=Yb;
        
        for j=1:ctra;
            cj=ColToTra(j);
            Gj=exp(mean(log(Yb(:,cj))));
            
            if la~=0
                Z(:,cj)=(Yb(:,cj).^la-1)/(la*(Gj^(la-1)));
            else
                Z(:,cj)=Gj*log(Yb(:,cj));
            end
        end
        
        dZ=log(det(cov(Z)));
    end



% lik computes the likelihood when different lambdas are possible for
% different variables
    function dZ=likYJ(la)
        
        Z=Yb;
        
        for j=1:ctra;
            laj=la(j);
            cj=ColToTra(j);
            
            nonnegs = Yb(:,cj) > 0;
            negs = ~nonnegs;
            % YJ transformation is the Box-Cox transformation of
            % y+1 for nonnegative values of y
            if laj ~=0
                Z(nonnegs,cj)= ((Yb(nonnegs,cj)+1).^laj-1)/laj;
            else
                Z(nonnegs,cj)= log(Yb(nonnegs,cj)+1);
            end
            
            % YJ transformation is the Box-Cox transformation of
            %  |y|+1 with parameter 2-lambda for y negative.
            if 2-laj~=0
                Z(negs,cj) = - ((-Yb(negs,cj)+1).^(2-laj)-1)/(2-laj);
            else
                Z(negs,cj) = -log(-Yb(negs,cj)+1);
            end
            
            % transformation is normalized so that its
            % Jacobian will be 1
            Z(:,cj)=Z(:,cj) * (exp(mean(log(   (1 + abs(Yb(:,cj))).^(2 * nonnegs - 1)) )))^(1 - laj);
            
        end
        
        dZ=log(det(cov(Z)));
        % disp(dZ);
    end


% lik computes the likelihood when just one value of lambda is allowed
    function dZ=likonelambdaYJ(la)
        Z=Yb;
        
        for j=1:ctra;
            cj=ColToTra(j);
            
            
            nonnegs = Yb(:,cj) >= 0;
            negs = ~nonnegs;
            % YJ transformation is the Box-Cox transformation of
            % y+1 for nonnegative values of y
            if la ~=0
                Z(nonnegs,cj)= ((Yb(nonnegs,cj)+1).^la-1)/la;
            else
                Z(nonnegs,cj)= log(Yb(nonnegs,cj)+1);
            end
            
            % YJ transformation is the Box-Cox transformation of
            %  |y|+1 with parameter 2-lambda for y negative.
            if 2-la~=0
                Z(negs,cj) = - ((-Yb(negs,cj)+1).^(2-la)-1)/(2-la);
            else
                Z(negs,cj) = -log(-Yb(negs,cj)+1);
            end
            
            % If Jacobian ==true the transformation is normalized so that its
            % Jacobian will be 1
            Z(:,cj)=Z(:,cj) * (exp(mean(log(   (1 + abs(Yb(:,cj))).^(2 * nonnegs - 1)) )))^(1 - la);
        end
        
        dZ=log(det(cov(Z)));
    end
end
%FScategory:MULT-Transformations
