function [out]=FSRfan(y,X,varargin)
%FSRfan monitors the values of the score test statistic for each lambda
%
%<a href="matlab: docsearchFS('FSRfan')">Link to the help function</a>
%
% The transformations for negative and positive responses were determined
% by Yeo and Johnson (2000) by imposing the smoothness condition that the
% second derivative of zYJ(λ) with respect to y be smooth at y = 0. However
% some authors, for example Weisberg (2005), query the physical
% interpretability of this constraint which is oftern violated in data
% analysis. Accordingly, Atkinson et al (2019) and (2020) extend the
% Yeo-Johnson transformation to allow two values of the transformations
% parameter: $\lambda_N$ for negative observations and $\lambda_P$ for
% non-negative ones.
% FSRfan monitors:
% 1) the t test associated with the constructed variable computed assuming
% the same transformation parameter for positive and negative observations
% fixed. In short we call this test, "global score test for positive
% observations".
% 2) the t test associated with the constructed variable computed assuming
% a different transformation for positive observations keeping the value of
% the transformation parameter for negative observations fixed. In short we
% call this test, "test for positive observations".
% 3) the t test associated with the constructed variable computed assuming
% a different transformation for negative observations keeping the value of
% the transformation parameter for positive observations fixed. In short we
% call this test, "test for negative observations".
% 4) the F test for the joint presence of the two constructed variables
% described in points 2) and 3.
% 4) the F likelihood ratio test based on the MLE of $\lambda_P$ and
% $\lambda_N$. In this case the residual sum of squares of the null model
% bsaed on a single trasnformation parameter $\lambda$ is compared with the
% residual sum of squares of the model based on data transformed data using
% MLE of $\lambda_P$ and $\lambda_N$.
%
%
%  Required input arguments:
%
%    y:         Response variable. Vector. A vector with n elements that
%               contains the response
%               variable.  It can be either a row or a column vector.
%    X :        Predictor variables. Matrix. Data matrix of explanatory
%               variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables.
%               Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing or
%               infinite values will automatically be excluded from the
%               computations. NOTICE THAT THE INTERCEPT MUST ALWAYS BE INCLUDED.
%
%  Optional input arguments:
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%       nocheck :   Check input arguments. Boolean.
%                   If nocheck is equal to true no check is performed
%                   on matrix y and matrix X. Notice that y and X are left
%                   unchanged. In other words the additional column of ones
%                   for the intercept is not added. As default nocheck=false.
%                   Example - 'nocheck',true
%                   Data Types - boolean
%
%           la  :   values of the transformation parameter for which it is
%                   necessary to compute the score test. Vector.
%                   Default value of lambda is la=[-1 -0.5 0 0.5 1]; that
%                   is the five most common values of lambda
%                   Example - 'la',[-1 -0.5]
%                   Data Types - double
%
%           h   :   The number of observations that have determined the
%                   least trimmed (median of) squares estimator. Integer.
%                   Generally h is an integer greater or equal than
%                   [(n+size(X,2)+1)/2] but smaller then n
%                   Example - 'h',5
%                   Data Types - double
%
%       nsamp   :   Number of subsamples which will be extracted to find
%                   the robust estimator. Scalar.
%                   If nsamp=0 all subsets will be extracted. They will be
%                   (n choose p). Remark: if the number of all possible
%                   subset is <1000 the default is to extract all subsets
%                   otherwise just 1000. If nsamp is a matrix of size
%                   r-by-p, it contains in the rows the subsets which sill
%                   have to be extracted. For example, if p=3 and nsamp=[ 2
%                   4 9; 23 45 49; 90 34 1]; the first subset is made up of
%                   units [2 4 9], the second subset of units [23 45 49]
%                   and the third subset of units [90 34 1];
%                   Example - 'nsamp',1000
%                   Data Types - double
%
%       lms     :   Criterion to use to find the initlal
%                 subset to initialize the search. Scalar or Vector.
%                 If lms=1 (default) Least Median of Squares is
%                 computed, else Least trimmed of Squares is computed.
%                 If, lms is matrix with size
%                 p-1+intercept-by-length(la) it contains in column
%                 j=1,..., lenght(la) the list of units forming the initial
%                 subset for the search associated with la(j). In this last
%                 case previous input option nsamp is ignored.
%                 Example - 'lms',1
%                 Data Types - double
%
%        family :   string which identifies the family of transformations which
%                   must be used. Character. Possible values are 'BoxCox'
%                   (default), 'YJ', 'YJpn' or 'YJall'.
%                   The Box-Cox family of power transformations equals
%                   $(y^{\lambda}-1)/\lambda$ for $\lambda$ not equal to zero,
%                   and $\log(y)$ if $\lambda = 0$.
%                   The Yeo-Johnson (YJ) transformation is the Box-Cox
%                   transformation of $y+1$ for nonnegative values, and of
%                   $|y|+1$ with parameter $2-\lambda$ for $y$ negative.
%                   Remember that BoxCox can be used just
%                   if input y is positive. Yeo-Johnson family of
%                   transformations does not have this limitation.
%                   If family is 'YJpn' Yeo-Johnson family is applied but in
%                   this case it is also possible to monitor (in the output
%                   arguments out.Scorep and out.Scoren) the score test
%                   respectively for positive and negative observations.
%                   If family is 'YJall', it is also possible to monitor
%                   the joint F test for the presence of the two
%                   constructed variables for positive and negative
%                   observations.
%                   Example - 'family','YJ'
%                   Data Types - char
%
%      scoremle: likelihood ratio test for the two different transformation
%                parameters $\lambda_P$ and $\lambda_N$. Boolean.
%                If scoremle is true it is possible to compute the
%                likelihood ratio test. In this case the residual sum of
%                squares of the null model bsaed on a single trasnformation
%                parameter $\lambda$ is compared with the residual sum of
%                squares of the model based on data transformed data using
%                MLE of $\lambda_P$ and $\lambda_N$. If scoremle is true it
%                is possible through following option usefmin, to control
%                the parameters of the optmization routine.
%               Example - 'scoremle',true
%               Data Types - logical
%
%    usefmin :  use solver to find MLE of lambda. Boolean or struct.
%               if usefmin is true or usefmin is a struct it is
%               possible to use MATLAB solvers fminsearch or fminunc to
%               find the maximum likelihood estimates of $\lambda_P$ and
%               $\lambda_N$. The default value of usefmin is false that is
%               solver is not used and the likelihood is evaluated at the
%               grid of points with steps 0.01.
%               If usefmin is a structure it may contain the following
%               fields:
%               usefmin.MaxIter = Maximum number of iterations (default is 1000).
%               usefmin.TolX   = Termination tolerance for the parameters
%                   (default is 1e-7).
%               usefmin.solver = name of the solver. Possible values are
%                   'fminsearch' (default) and 'fminunc'. fminunc needs the
%                   optimization toolbox.
%               usefmin.displayLevel = amount of information displayed by
%                   the algorithm. possible values are 'off' (displays no
%                   information, this is the default), 'final' (displays
%                   just the final output) and 'iter' (displays iterative
%                   output to the command window).
%               Example - 'usefmin',true
%               Data Types - boolean or struct
%
%       init    :   Search initialization. Scalar.
%                   It specifies the initial subset size to start
%                   monitoring the value of the score test, if init is not
%                   specified it will be set equal to:
%                    p+1, if the sample size is smaller than 40;
%                    min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%                    Example - 'init',100 starts monitoring from step m=100
%                    Data Types - double
%
%       plots   :  Plot on the screen. Scalar.
%                   If plots=1 the fan plot is produced
%                   else (default) no plot is produced.
%                   Example - 'plots',1
%                   Data Types - double
%                   REMARK: all the following options work only if plots=1
%
%       conflev :   Confidence level. Scalar or vector. Confidence level
%                   for the bands (default is 0.99, that is we plot two
%                   horizontal lines in correspondence of value -2.58 and
%                   2.58).
%                   Example - 'conflev',[0.9 0.95 0.99]
%                   Data Types - double
%
%       titl    :   a label for the title. Character.
%                   default: 'Fan plot'
%                   Example - 'titl','my title'
%                   Data Types - char
%
%       labx    :   a label for the x-axis. Character.
%                   default: 'Subset size m'
%                   Example - 'labx','my labx'
%                   Data Types - char
%
%       laby    :   a label for the y-axis. Character.
%                   default:'Score test statistic'
%                   Example - 'laby','my laby'
%                   Data Types - char
%
%       xlimx   :   Minimum and maximum of the x axis. Vector.
%                   Default value is [init n]
%                   Example - 'xlimx',[0 1]
%                   Data Types - double
%
%       ylimy   :  Minimum and maximum of the y axis. Vector.
%                   Default value for ylimy(1)=max(min(score_test),-20).
%                   Default value for ylimy(2)=min(max(score_test),20).
%                   Example - 'ylimx',[0 1]
%                   Data Types - double
%
%       lwd     :   linewidth of the curves which
%                   contain the score test. Scalar.
%                   Default line width=2.
%                   Example - 'lwd',5
%                   Data Types - double
%
%       lwdenv  :   width of the lines associated
%                   with the envelopes. Scalar.
%                   Default is lwdenv=1.
%                   Example - 'lwdenv',5
%                   Data Types - double
%
%       FontSize:   font size of the labels of  the axes. Scalar.
%                   Default value is 12.
%                   Example - 'FontSize',20
%                   Data Types - double
%
%    SizeAxesNum:   Scalar which controls the size of the numbers of the
%                   axes. Scalar.
%                   Default value is 10.
%                  Example - 'SizeAxesNum',12
%                  Data Types - double
%
%         msg   : Level of output to display. Scalar.
%                   scalar which controls whether to display or not
%                   messages on the screen. Scalar.
%                   If msg==1 (default) messages are
%                   displayed on the screen about estimated time to compute
%                   the LMS (LTS) for each value of lamabda else no message
%                   is displayed on the screen
%                  Example - 'msg',1
%                  Data Types - double
%
%       tag     :   handle of the plot which is about to be created.
%                   Character. The default is to use tag 'pl_fan'. Notice
%                   that if the program finds a plot which has a tag equal
%                   to the one specified by the user, then the output of
%                   the new plot overwrites the existing one in the same
%                   window else a new window is created Example -
%                   'tag','mytag' Data Types - char
%  Output:
%
%         out:   structure which contains the following fields
%
%  out.Score  = (n-init+1) x length(la)+1 matrix containing the values of the
%               score test for each value of the transformation parameter:
%               1st col = fwd search index;
%               2nd col = value of the score test in each step of the
%               fwd search for la(1);
%               ...........
%               end col = value of the score test in each step of the fwd
%               search for la(end).
%  out.Scorep = (n-init+1) x length(la)+1 matrix containing the values of the
%               score test for positive observations for each value of the
%               transformation parameter.
%               1st col = fwd search index;
%               2nd col = value of the (positive) score test in each step
%               of the fwd search for la(1);
%               ...........
%               end col = value of the (positive) score test in each step
%               of the fwd search for la(end).
%               Note that this output is present only if input option
%               family is 'YJpn' or 'YJall'.
% out.Scoren  = (n-init+1) x length(la)+1 matrix containing the values of the
%               score test for positive observations for each value of the
%               transformation parameter:
%               1st col = fwd search index;
%               2nd col = value of the (negative) score test in each step
%               of the fwd search for la(1);
%               ...........
%               end col = value of the (negative) score test in each step
%               of the fwd search for la(end).
%               Note that this output is present only if input option
%               family is 'YJpn' or 'YJall'.
% out.Scoreb  = (n-init+1) x length(la)+1 matrix containing the values of the
%               score test for the joint presence of both constructed
%               variables (associated with positive and negative
%               observations) for each value of the transformation
%               parameter.  In this case the reference distribution is the
%               $F$ with 2 and subset_size-p degrees of freedom.
%               1st col = fwd search index (subset_size);
%               2nd col = value of the score test in each step
%               of the fwd search for la(1);
%               ...........
%               end col = value of the score test in each step
%               of the fwd search for la(end).
%               Note that this output is present only if input option
%               family is 'YJall'
% out.Scoremle  = (n-init+1) x length(la)+1 matrix containing the values of the
%               (score) likelihood ratio test for the joint presence of both constructed
%               variables (associated with positive and negative
%               observations) for each value of the transformation
%               parameter.  In this case the reference distribution is the
%               $F$ with 2 and subset_size-p degrees of freedom.
%               1st col = fwd search index (subset_size);
%               2nd col = value of the score test in each step
%               of the fwd search for la(1);
%               ...........
%               end col = value of the score test in each step
%               of the fwd search for la(end).
%               Note that this output is present only if input option
%               scoremle is true
%    out.laMLE = (n-init+1) x 2*length(la)+1 matrix containing the values of the
%               maximum ikelihood estimate of laP and laN.
%               Columns 2:3 are associated with  the search which has
%               ordered the data using to la(1);
%               .........
%               Columns 2*length(la):2*length(la)+1 are associated with
%               the search which has ordered the data using to
%               la(length(la)).
%               Note that out.laMLE(end,2)=out.laMLE(end,2)=...=out.laMLE(end,2*length(la))
%               because all these variables contain the MLE of laP based on
%               all the observations. Similarly notice that
%               out.laMLE(end,3)=out.laMLE(end,5)=...=out.laMLE(end,2*length(la)+1)
%               because all these variables contain the MLE of laN based on
%               all the observations. This output is present only if input
%               option scoremle is true.
%  out.la     = vector containing the values of lambda for which fan plot
%               is constructed
%  out.bs     = matrix of size p x length(la) containing the units forming
%               the initial subset for each value of lambda
%  out.Un     = cell of size length(la).
%               out.Un{i} is a n-init) x 11 matrix which contains the
%               unit(s) included in the subset at each step in the search
%               associated with la(i).
%               REMARK: in every step the new subset is compared with the old subset. Un
%               contains the unit(s) present in the new subset but not in
%               the old one Un(1,:) for example contains the unit included
%               in step init+1 ... Un(end,2) contains the units included in the
%               final step of the search
%  out.y      = A vector with n elements that contains the response
%               variable which has been used
%  out.X      = Data matrix of explanatory variables
%               which has been used (it also contains the column of ones if
%               input option intercept was missing or equal to 1)
%
%
% See also: Score, ScoreYJ, ScoreYJpn
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
% Atkinson, A.C. and Riani, M. (2002a), Tests in the fan plot for robust,
% diagnostic transformations in regression, "Chemometrics and Intelligent
% Laboratory Systems", Vol. 60, pp. 87-100.
% Atkinson, A.C. Riani, M., Corbellini A. (2019), The analysis of
% transformations for profit-and-loss data, Journal of the Royal
% Statistical Society, Series C, "Applied Statistics",
% https://doi.org/10.1111/rssc.12389
% Atkinson, A.C. Riani, M. and Corbellini A. (2020), The Box-Cox
% Transformation: Review and Extensions, "Statistical Science", in press.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRfan')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% FSRfan with all default options.
    % Store values of the score test statistic
    % for the five most common values of $\lambda$.
    % Produce also a fan plot and display it on the screen.
    % Common part to all examples: load wool dataset.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % Function FSRfan stores the score test statistic.
    % In this case we use the five most common values of lambda are considered
    [out]=FSRfan(y,X);
    fanplot(out);
    %The fan plot shows the log transformation is diffused throughout the data and does not depend on the presence of particular observations.
%}

%{
    % FSRfan with optional arguments.
    % Produce a personalized fan plot with required font sizes
    % for labels and axes.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    [out]=FSRfan(y,X,'plots',1,'FontSize',16,'SizeAxesNum',16);
%}

%{
    % Example specifying $\lambda$.
    % Produce a fan plot for each value of $\lambda$ inside vector la.
    % Extract in matrix Un the units which entered the search in each step
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    la=[-1 -0.5 0 0.5];
    [out]=FSRfan(y,X,'la',la,'plots',1);
    Unsel=cell2mat(out.Un);
    lla=length(la);
    nr=size(Unsel,1)/lla;
    Un=[Unsel(1:nr,1) reshape(Unsel(:,2),nr,lla)];
%}

%{
    % Example specifying the confidence level and the initial
    % starting point for monitoring.
    % Construct fan plot specifying the confidence level and the initial
    % starting point for monitoring.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    [out]=FSRfan(y,X,'init',size(X,2)+2,'nsamp',0,'conflev',0.95,'plots',1);
%}

%{
    % Example with starting point based on LTS.
    % Extraction of all subsamples, construct
    % fan plot specifying the confidence level and the initial starting
    % point for monitoring based on p+2 observations strong line width for
    % lines associated with the confidence bands.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    [out]=FSRfan(y,X,'init',size(X,2)+2,'nsamp',0,'lms',0,'lwdenv',3,'plots',1);
%}

%{
    %% Fan plot using fidelity cards data.
    % In the example, la is the vector contanining the most common values
    % of the transformation parameter.
    % Store the score test statistics for the specified values of lambda
    % and automatically produce the fan plot
    XX=load('loyalty.txt');
    namey='Sales'
    y=XX(:,end);
    nameX={'Number of visits', 'Age', 'Number of persons in the family'};
    X=XX(:,1:end-1);
    % la = vector contanining the most common values of the transformation
    % parameter
    la=[0 1/3 0.4 0.5];
    % Store the score test statistics for the specified values of lambda
    % and automatically produce the fan plot
    [out]=FSRfan(y,X,'la',la,'init',size(X,2)+2,'plots',1,'lwd',3);
    % The fan plot shows the even if the third root is the best value of the
    % transformation parameter at the end of the search in earlier steps it
    % lies very close to the upper rejection region. The best value of the
    % transformation parameter seems to be the one associated with l=0.4
    % which is always the confidence bands but at the end of search, due to
    % the presence of particular observations it goes below the lower
    % rejection line.
%}

%{
    %% Compare BoxCox with Yeo and Johnson transformation.
    % Store values of the score test statistic
    % for the five most common values of $\lambda$.
    % Produce also a fan plot and display it on the screen.
    % Common part to all examples: load wool dataset.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % Store the score test statistic using Box Cox transformation.
    [outBC]=FSRfan(y,X,'nsamp',0);
    % Store the score test statistic using Yeo and Johnson transformation.
    [outYJ]=FSRfan(y,X,'family','YJ','nsamp',0);
    fanplot(outBC,'titl','Box Cox');
    fanplot(outYJ,'titl','Yeo and Johnson','tag','YJ');
    disp('Maximum difference in absolute value')
    disp(max(max(abs(outYJ.Score-outBC.Score))))
%}


%{
    %% Example of monitoring of score test for positive and negative obseravations.
    rng('default')
    rng(10)
    close all
    n=200;

    X=randn(n,3);
    beta=[ 1; 1; 1];
    sig=0.5;
    y=X*beta+sig*randn(n,1);

    disp('Fit in the true scale')
    disp('R2 of the model in the true scale')
    if verLessThanFS(8.1)
        out=regstats(y,X,'linear',{'rsquare'});
        disp(out.rsquare)
    else
        outlmori=fitlm(X,y);
        disp(outlmori.Rsquared.Ordinary)
    end
    [~,~,BigAx]=yXplot(y,X,'tag','ori');
    title(BigAx,'Data in the original scale')


    % Find the data to transform
    la=-0.25;
    ytra=normYJ(y,[],la,'inverse',true);
    if any(isnan(ytra))
        disp('response with missing values')
    end

    disp('Fit in the transformed scale')
    disp('R2 of the model in the wrong (inverse) scale')
    if verLessThanFS(8.1)
        out=regstats(ytra,X,'linear',{'rsquare'});
        disp(out.rsquare)
    else
        outlmtra=fitlm(X,ytra);
        disp(outlmtra.Rsquared.Ordinary)
    end
    [~,~,BigAx]=yXplot(ytra,X,'tag','tra','namey','Data to transform (zoom of y [0 500])','ylimy',[0 500]);
    title(BigAx,'Data in the inverse scale')

    la=[ -0.5 -0.25 0];
    out=FSRfan(ytra,X,'la',la,'family','YJpn','plots',1,'init',round(n/2),'msg',0);
    title('Extended fan plot')
%}

%{
    %% Example of monitoring all score tests (also the F test).
    rng('default')
    rng(1000)
    close all
    n=200;

    X=randn(n,3);
    beta=[ 1; 1; 1];
    sig=0.5;
    y=X*beta+sig*randn(n,1);

    % Find the data to transform
    la=0;
    ytra=normYJ(y,[],la,'inverse',true);
    if any(isnan(ytra))
        disp('response with missing values')
    end

    la=[ -0.1 0 0.1];
    % In this case  family is YJall
    out=FSRfan(ytra,X,'la',la,'family','YJall','plots',1,'init',round(n/2),'msg',0);
    
%}

%{
   %% Comparison of  F test based on constructed variables with F test based on MLE.
    rng('default')
    rng(100)
    close all
    n=200;

    X=randn(n,3);
    beta=[ 1; 1; 1];
    sig=0.5;
    y=X*beta+sig*randn(n,1);

    % Find the data to transform
    la=0;
    ytra=normYJ(y,[],la,'inverse',true);
    if any(isnan(ytra))
        disp('response with missing values')
    end

    la=[ -0.1 0 0.1];
    % Monitor test based on MLE using option scoremle
    scoremle= true;
    out=FSRfan(ytra,X,'la',la,'family','YJall','plots',1,'init',round(n/2),'msg',false,'scoremle',true);
%}

%% Beginning of code

% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

% If the number of all possible subsets is <1000 the default is to extract
% all subsets, otherwise just 1000.
ncomb=bc(n,p);
nsamp=min(1000,ncomb);

% REMARK: a fast approximation of the bc computed above is:
% ncomb=floor(exp( gammaln(n+1) - gammaln(n-p+1) - gammaln(p+1) ) + .5);

h=floor(0.5*(n+p+1));
if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));
end
family='BoxCox';
scoremle=false;
usefmin=true;
plo=1;
lms=1;
conflev=0.99;
msg=1;
tag='pl_fan';
la=[-1 -0.5 0 0.5 1];
nocheck=false;
lwd=2;
lwdenv=1;
FontSize=12;
xlimx=[];
ylimy=[];
SizeAxesNum=10;
labx='Subset size m';
laby='Score test statistic';
intercept=true;
titl='Fan plot';

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    options=struct('la',la,'h',h,...
        'nsamp',nsamp,'lms',lms,'plots',plo,'init',init,'conflev',conflev,...
        'titl',titl,'labx',labx,...
        'laby',laby,'xlimx',xlimx,'ylimy',ylimy,'lwd',lwd,...
        'lwdenv',lwdenv,'FontSize',FontSize,'SizeAxesNum',SizeAxesNum,...
        'tag',tag,'intercept',intercept,'msg',msg,'nocheck',nocheck,'family',family,...
        'scoremle',scoremle,'usefmin',usefmin);
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRfan:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    
    h=options.h;
    lms=options.lms;
    plo=options.plots;
    nsamp=options.nsamp;
    msg=options.msg;
    nocheck=options.nocheck;
    family=options.family;
    scoremle=options.scoremle;
    usefmin=options.usefmin;
    tag=options.tag;
    la=options.la;
    init=options.init;
    lwd=options.lwd;
    lwdenv=options.lwdenv;
    conflev=options.conflev;
    labx=options.labx;
    laby=options.laby;
    titl=options.titl;
    
end

if strcmp(family,'BoxCox')
    BoxCox=1;
elseif strcmp(family,'YJ')
    BoxCox=0;
elseif strcmp(family,'YJpn')
    BoxCox=-1;
elseif strcmp(family,'YJall')
    BoxCox=-2;
else
    warning('FSDA:FSRfan:WrongFamily','Transformation family which has been chosen is not supported')
    error('FSDA:FSRfan:WrongFamily','Supported values are BoxCox or YJ or YJpn or YJall')
end

% Specify where to send the output of the current procedure if options plot
% =1
if plo==1
    h1=findobj('-depth',1,'tag',tag);
    if (~isempty(h1))
        clf(h1);
        figure(h1)
        axes;
    else
        figure;
        % include specified tag in the current plot
        set(gcf,'tag',tag);
    end
end

if  init <p+1
    fprintf(['Attention : init should be larger than p+1. \n',...
        'It is set to p+2.']);
    init=p+2;
end



%% Start of the forward search
zeron1=false(n,1);

% Initialization of the n x 1 Boolean vector which contains a true in
% correspondence of the units belonging to subset in each step
bsbT=zeron1;

seq=(1:n)';

%  Unlai is a Matrix whose 2:11th col contains the unit(s) just included.
Unlai = cat(2 , (init+1:n)' , NaN(n-init,10));

% Un = cell which will contain the matrices Unlai for each value of lambda
Un=cell(length(la),1);

lla=length(la);

% Initialize matrix which will contain the score test
Sco=[((init):n)'  NaN(n-init+1,lla)];

if BoxCox==-1
    Scop=Sco;
    Scon=Sco;
elseif BoxCox==-2
    Scop=Sco;
    Scon=Sco;
    Scob=Sco;
end

if scoremle == true
    Scomle=Sco;
    laMLE=[((init):n)'  NaN(n-init+1,2*lla)];
    % SSElaMLE=laMLE(:,1:2);
end

% The second column of matrix r will contain the OLS residuals at each step
% of the forward search
r=[seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100 = 1000*(1:1:ceil(n/1000));
seq100boo=false(n,1);
seq100boo(seq100)=true;


binit=zeros(p,lla);

% Preextract subsample once and for all for all values of lambda;
if lla>1
    [nsamp] = subsets(nsamp,n,p);
end

% loop over the values of \lambda
for i=1:lla
    
    if BoxCox==1
        % Construct transformed z according to power tansformation
        if abs(la(i))<1e-8
            z=log(y);
        else
            z=y.^la(i);
        end
    else
        z=normYJ(y,1,la(i),'Jacobian',false);
    end
    
    % Find initial subset to initialize the search using as y transformed
    % vector z
    if size(lms,1)==1
        [out]=LXS(z,X,'lms',lms,'h',h,'nsamp',nsamp,'nocheck',true,'msg',msg);
        bsb=out.bs;
        % Store information about the units forming subset for each value of
        % lambda
        binit(:,i)=bsb';
    else
        bsb=lms(:,i);
        % Store information about the units forming subset for each value of
        % lambda
        binit(:,i)=bsb;
    end
    
    bsbT(bsb)=true;
    
    % bsb=[1 8 12 15];
    %ini0 = initial value for forward search loop
    ini0=length(bsb);
    
    % FS loop for a particular value of vector la
    zb=z(bsb);
    yb=y(bsb);
    Xb=X(bsb,:);
    
    % last correctly computed beta oefficients
    blast=NaN(p,1);
    
    if nocheck==false && (rank(Xb)~=p)
        warning('FSRfan:message','The provided initial subset does not form full rank matrix');
        % FS loop will not be performed
    else
        for mm=ini0:n
            % if n>1000 show every 100 steps the fwd search index
            if  msg==1 && seq100boo(mm) == true
                % OLD CODE if length(intersect(mm,seq100))==1
                disp(['m=' int2str(mm)]);
            end
            
            
            if (mm>=init)
                if BoxCox==1
                    % Compute and store the value of the score test
                    [outSC]=Score(yb,Xb,'la',la(i),'nocheck',true);
                    % Store score test for the units belonging to subset
                    Sco(mm-init+1,i+1)=outSC.Score;
                    
                elseif BoxCox==0
                    % Compute and store the value of the score test using Yeo
                    % and Johnson transformation (just the global test)
                    [outSC]=ScoreYJ(yb,Xb,'la',la(i),'nocheck',true);
                    % Store score test for the units belonging to subset
                    Sco(mm-init+1,i+1)=outSC.Score;
                    
                else % in this case BoxCox==-1 || BoxCox==-2
                    %[outSC]=ScoreYJ(yb,Xb,'la',la(i),'nocheck',true);
                    % [outSCpn]=ScoreYJpn(yb,Xb,'la',la(i),'nocheck',true);
                    % [outSCpn]=ScoreYJpn(yb,Xb,'la',la(i),'nocheck',true);
                    %                     if i==1
                    if mm==init
                        clear cachedlahatPreviousStep
                    end
                    
                    [outSCpn]=ScoreYJall(yb,Xb,'la',la(i),'scoremle',scoremle,'nocheck',true,'usefmin',usefmin);
                    if scoremle == true
                        laMLE(mm-init+1,i*2:i*2+1)=outSCpn.laMLE;
                    end
                end
                
                if BoxCox<=-1
                    Sco(mm-init+1,i+1)=outSCpn.Score(1,1);
                    Scop(mm-init+1,i+1)=outSCpn.Score(1,2);
                    Scon(mm-init+1,i+1)=outSCpn.Score(1,3);
                end
                if BoxCox==-2
                    Scob(mm-init+1,i+1)=outSCpn.Score(1,4);
                end
                
                if scoremle == true
                    Scomle(mm-init+1,i+1)=outSCpn.Score(1,5);
                end
                
            end
            
            if nocheck==true
                NoRankProblem=true;
            else
                % Compute b using transformed vector zb
                NoRankProblem=(rank(Xb) == p);
            end
            
            if NoRankProblem  % rank is ok
                b=Xb\zb;
                blast=b;   % Store correctly computed b for the case of rank problem
            else   % number of independent columns is smaller than number of parameters
                warning('FSR:FSRfan','Rank problem in step %d: Beta coefficients are used from the most recent correctly computed step',mm);
                b=blast;
            end
            
            
            % e= (n x 1) vector of residuals for all units using b estimated
            % using subset and transformed response
            e=z-X*b;
            
            % r_i =e_i^2
            r(:,2)=e.^2;
            
            
            if mm<n
                
                % store units forming old subset in vector oldbsb
                oldbsbT=bsbT;
                
                % order the r_i and include the smallest among the units
                %  forming the group of potential outliers
                % ord=sortrows(r,2);
                [~,ord]=sort(r(:,2));
                
                % bsb= units forming the new  subset
                bsb=ord(1:(mm+1),1);
                
                bsbT=zeron1;
                bsbT(bsb)=true;
                
                
                Xb=X(bsb,:);  % subset of X
                yb=y(bsb);    % subset of y
                zb=z(bsb);    % subset of z
                
                if mm>=init
                    
                    % unit = vector containing units which just entered subset;
                    % unit=setdiff(bsb,oldbsb);
                    % new instruction to find unit
                    unit=seq(bsbT & ~oldbsbT);
                    
                    if length(unit)<=10
                        Unlai(mm-init+1,2:(length(unit)+1))=unit;
                    else
                        if msg==1
                            disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        end
                        Unlai(mm-init+1,2:end)=unit(1:10);
                    end
                end
            end
            
        end
    end  % rank check
    % Store in cell Un matrix Unlai
    Un{i}=Unlai;
end

%% Structure returned by function FSRfan
out=struct;

out.Score=Sco;
out.la=la;
out.bs=binit;
out.Un=Un;
out.y=y;
out.X=X;
if BoxCox==-1
    out.Scorep=Scop;
    out.Scoren=Scon;
elseif BoxCox==-2
    out.Scorep=Scop;
    out.Scoren=Scon;
    out.Scoreb=Scob;
end

if scoremle == true
    out.Scoremle=Scomle;
    out.laMLE=laMLE;
end

if plo==1
    if BoxCox == -2 && scoremle == false
        subplot(2,1,1);
    end
    
    % Specify the line type for the units inside vector units
    slin={'-';'--';':';'-.'};
    slin=repmat(slin,ceil(lla/4),1);
    
    % Specify the color for the trajectories
    ColorOrd=[{[0 0 1]}; {[1 0 0]}; {[1 0 1]}; {[1 1 0]}; {[0 1 0]}; {[0 1 1]}];
    ColorOrd=repmat(ColorOrd,4,1);
    
    
    % plot the lines associated with the score test lwd = line width of the
    % trajectories which contain the score test
    plot1=plot(Sco(:,1),Sco(:,2:end),'LineWidth',lwd);
    
    set(plot1,{'Color'}, ColorOrd(1:lla,:));
    
    set(plot1,{'LineStyle'},slin(1:lla));
    
    if BoxCox <= -1
        hold('on')
        plotp=plot(Scop(:,1),Scop(:,2:end),'LineWidth',lwd);
        set(plotp,{'LineStyle'},{'--'});
        % set(plotp,{'LineStyle'},slin(1:lla));
        %  set(plotp,{'Color'}, ColorOrd(1:lla,:));
        % Green color for positive values of y
        set(plotp,{'Color'},{'g'})
        
        plotn=plot(Scon(:,1),Scon(:,2:end),'LineWidth',lwd);
        set(plotn,{'LineStyle'},slin(1:lla));
        set(plotn,{'LineStyle'},{'--'});
        % set(plotn,{'Color'}, ColorOrd(1:lla,:));
        set(plotn,{'Color'},{'k'})
        
        set(plot1,{'LineStyle'},{'-'});
        
    end
    
    % Confidence bands lwdenv = line width of the curves associated with
    % the envelopes
%     v=axis;
%     quant=sqrt(chi2inv(conflev,1));
%     line([v(1),v(2)],[quant,quant],'color','r','LineWidth',lwdenv);
%     line([v(1),v(2)],[-quant,-quant],'color','r','LineWidth',lwdenv);
 
    rangeaxis=axis;
quant = sqrt(chi2inv(conflev,1));
numconflev=length(conflev);
V=repmat([rangeaxis(1);rangeaxis(2)],1,2*numconflev);
QUANT=[[quant;quant],[ -quant;-quant]];
line(V, QUANT,'LineWidth',lwdenv,'color','r','LineWidth',lwdenv);
 
    
    
    if size(la,2)>1
        la=la';
    end
    text(n*ones(lla,1),Sco(end,2:end)',num2str(la));
    
    % set the x and y axis
    
    if ~isempty(xlimx)
        xlim(xlimx);
    end
    
    if isempty(ylimy)
        
        if BoxCox <= -1
            Scog=[Sco(:,2:end);Scop(:,2:end);Scon(:,2:end)];
            maxSco=max([max(Scog) quant]);
            minSco=min([min(Scog) -quant]);
            ylim1=max(-20,minSco);
            ylim2=min(20,maxSco);
        else
            % Use default limits for y axis
            ylim1=max(-20,min([min(Sco(:,2:end)) -quant]));
            ylim2=min(20,max([max(Sco(:,2:end)) quant]));
        end
        ylim([ylim1 ylim2]);
    else
        
        % Use limits specified by the user
        ylim(ylimy);
    end
    
    
    % Main title of the plot and labels for the axes
    
    title(titl);
    
    % FontSize = font size of the axes labels
    
    % Add to the plot the labels for values of la Add the horizontal lines
    % representing asymptotic confidence bands
    xlabel(labx,'Fontsize',FontSize);
    ylabel(laby,'Fontsize',FontSize);
    
    % FontSizeAxesNum = font size for the axes numbers
    set(gca,'FontSize',SizeAxesNum)
    box on
    
    
    
    
    if BoxCox==-2
        if scoremle == false
            subplot(2,1,2)
        else
            figure
            subplot(2,1,1)
        end
        
        plot1=plot(Scob(:,1),Scob(:,2:end),'LineWidth',lwd);
        set(plot1,{'Color'}, ColorOrd(1:lla,:));
        
        set(plot1,{'LineStyle'},slin(1:lla));
        text(n*ones(lla,1),Scob(end,2:end)',num2str(la));
        
        % Confidence bands lwdenv = line width of the curves associated with
        % the envelopes
        lwdenv=options.lwdenv;
        quant=conflev;
        
        % Compute and superimpose envelopes based on the F distribution
        EnvF=[Scob(:,1) zeros(size(Scob,1),length(quant))];
        
        for i=1:size(Scob,1)
            ast=finv(quant,2,Scob(i,1)-p);
            EnvF(i,2:end)=ast;
        end
        
        line(EnvF(:,1),EnvF(:,2:end),'LineStyle','-','Color','r','LineWidth',lwdenv);
        
        % Add text associated with the envelopes which have been used
        % text(repmat(n,length(quant),1),EnvF(end,2:end)',num2str(quant'));
        
        % Main title of the plot and labels for the axes
        labx=options.labx;
        laby=options.laby;
        
        
        title('F test for both constructed variables','Interpreter','latex','FontSize',14);
        
        % FontSize = font size of the axes labels
        FontSize =options.FontSize;
        
        % Add to the plot the labels for values of la Add the horizontal lines
        % representing asymptotic confidence bands
        xlabel(labx,'Fontsize',FontSize);
        ylabel(laby,'Fontsize',FontSize);
        
        % FontSizeAxesNum = font size for the axes numbers
        SizeAxesNum=options.SizeAxesNum;
        set(gca,'FontSize',SizeAxesNum)
        box on
        
        if isempty(ylimy)
            
            % Use default limits for y axis
            ylim1=0;
            ylim2=min(100,max(max([Scob(:,2:end) EnvF(:,2:end)] )));
            ylim([ylim1 ylim2]);
        else
            
            % Use limits specified by the user
            ylim(ylimy);
        end
        
        if scoremle == true && BoxCox == -2
            subplot(2,1,2)
        elseif scoremle ==true
            figure
        end
        
        if scoremle ==true
            plot2=plot(Scomle(:,1),Scomle(:,2:end),'LineWidth',lwd);
            set(plot2,{'Color'}, ColorOrd(1:lla,:));
            
            set(plot2,{'LineStyle'},slin(1:lla));
            text(n*ones(lla,1),Scomle(end,2:end)',num2str(la));
            
            % Confidence bands lwdenv = line width of the curves associated with
            % the envelopes
            lwdenv=options.lwdenv;
            quant=conflev;
            
            % Compute and superimpose envelopes based on the F distribution
            EnvF=[Scomle(:,1) zeros(size(Scomle,1),length(quant))];
            
            for i=1:size(Scomle,1)
                ast=finv(quant,2,Scomle(i,1)-p);
                EnvF(i,2:end)=ast;
            end
            
            line(EnvF(:,1),EnvF(:,2:end),'LineStyle','-','Color','r','LineWidth',lwdenv);
            
            % Add text associated with the envelopes which have been used
            % text(repmat(n,length(quant),1),EnvF(end,2:end)',num2str(quant'));
            
            % Main title of the plot and labels for the axes
            labx=options.labx;
            laby=options.laby;
            
            
            title('F test based on MLE of $\lambda_P$ and $\lambda_N$','interpreter','latex','FontSize',14);
            
            % FontSize = font size of the axes labels
            FontSize =options.FontSize;
            
            % Add to the plot the labels for values of la Add the horizontal lines
            % representing asymptotic confidence bands
            xlabel(labx,'Fontsize',FontSize);
            ylabel(laby,'Fontsize',FontSize);
            
            % FontSizeAxesNum = font size for the axes numbers
            SizeAxesNum=options.SizeAxesNum;
            set(gca,'FontSize',SizeAxesNum)
            box on
            
            if isempty(ylimy)
                
                % Use default limits for y axis
                ylim1=0;
                ylim2=min(100,max(max([Scomle(:,2:end) EnvF(:,2:end)] )));
                ylim([ylim1 ylim2]);
            else
                
                % Use limits specified by the user
                ylim(ylimy);
            end
        end
    end
end


end
%FScategory:REG-Transformations
