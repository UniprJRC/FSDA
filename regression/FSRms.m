function outms = FSRms(y,X,varargin)
%FSRms performs robust model selection using flexible trimming in linear regression
%
%<a href="matlab: docsearchFS('FSRms')">Link to the help function</a>
%
% Required input arguments:
%
%       y:      Response variable. Vector. A vector with n elements that contains the response variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%       X :     Predictor variables. Matrix. Data matrix of explanatory variables (also called
%               'regressors') of dimension (n x (bigP-1)).
%               The intercept will be added in automatic way, so that the
%               dimension of the full model is bigP
%               Rows of X represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
%
% Optional input arguments:
%
%    intercept :  Indicator for constant term. true (default) | false. 
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%      init   : Search initialization. Scalar. 
%               It specifies the initial subset size to start
%               monitoring the required quantities, if
%               init is not specified it set equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%               Example - 'init',100 starts monitoring from step m=100 
%               Data Types - double
%
%         h   : The number of observations that have determined the least
%               trimmed squares estimator. Scalar.
%               h is an integer greater or
%               equal than [(n+p+1)/2] but smaller then n
%                 Example - 'h',round(n*0,75) 
%                 Data Types - double
%
%     nsamp   : Number of subsamples which will be extracted to find the
%                 robust estimator. Scalar.
%                   Number of subsamples which will be extracted to find the
%               robust estimator. If nsamp=0 all subsets will be extracted.
%               They will be (n choose smallp).
%               Remark. if the number of all possible subset is <1000 the
%               default is to extract all subsets otherwise just 1000.
%                 Example - 'nsamp',1000 
%                 Data Types - double
%
%         lms : Criterion to use to find the initlal
%                 subset to initialize the search. Scalar.
%                   If lms=1 (default) Least Median of Squares is
%               computed, else Least Trimmed of Squares is computed.
%                 Example - 'lms',1 
%                 Data Types - double
%
%     nocheck : Check input arguments. Boolean.
%               If nocheck is equal to true no check is performed on
%               matrix y and matrix X. Note that y and X are left
%               unchanged. In other words the additional column of ones
%               for the intercept is not added. As default nocheck=false.
%               Example - 'nocheck',true 
%               Data Types - boolean
%
%    smallpint: submodels to consider. Vector. It specifies which submodels 
%               (number of variables) must be considered. 
%               The default is to consider all models
%               from size 2 to size bigP-1. In other words, as default,
%               smallpint=(bigP-1):-1:2.
%               When smallpint=2 all submodels including one explanatory
%               variable and the constant will be considered.
%               When smallpint=3 all submodels including two explanatory
%               variables and a constant will be considered. ....
%               Example - 'smallpint',3 
%               Data Types - double
%
%      labels : names of the explanatory variables. Cell array of strings.
%               Cell array of strings of length bigP-1 containing the
%               names of the explanatory variables.
%                If labels is a missing
%               value the following sequence of strings will be
%               automatically created for labelling the column of matrix X
%               (1,2,3,4,5,6,7,8,9,A,B,C,D,E,E,G,H,I,J,K,...,Z)
%               Example - 'labels',{'1','2'} 
%               Data Types - cell
%
%     fin_step: Initial and final step. Vector with two elements.
%               Initial and final step of the search which has to be
%               monitored to choose the best models as specified in scalar
%               first_k.
%               The first element of the vector specifies the initial step of the search
%               which has to be monitored to choose the best models as
%               specified in scalar first_k below. The second element
%               specifies the ending point of the central part of the
%               search. This information will be used to create the
%               candlestick Cp plot.
%               If the elements of fin_step are integers greater or equal 1
%               they refer to the number of steps. For example if
%               fin_step=[10 3] the program considers the last 10 steps to
%               choose the best models and the central part of the search
%               is defined up to step n-3.
%               If the elements of fin_step are real numbers
%               alpha (0<alpha<0.5) in the interval (0 0.5] then the
%               program considers the last round(n*alpha) steps.
%               As default fin_step(1)=round(n*0.2) that is the last 20%
%               of the steps are considered.
%               As default fin_step(2)=round(n*0.05) that
%               is the central part of the search extends up to 95% of the
%               observations
%               Example - 'fin_step',[1 50] 
%               Data Types - double
%
%      first_k: Number of models to consider. Scalar. Number of best models to
%               consider in each of the last fin_step. 
%               For example if
%               first_k=5 in each of the fin_step the models which had
%               the 5 smallest values of Cp are considered. As default
%               first_k=3
%               Example - 'first_k',5
%               Data Types - double
%
%       ignore: submodels to consider. Scalar. 
%               If ignore=1, when dealing with p explanatory
%               variables, the submodels of the models with p+1
%               explanatory variables which were considered irrelevant
%               according to option ExclThresh, are not considered. As
%               default ignore=1, because this saves computational time.
%               If ignore is different from 1, for each p all submodels of
%               size p which contain a constant are considered
%               Example - 'ignore',1
%               Data Types - double
%
%   ExclThresh:  Exclusion threshold. Scalar.
%               It has effect only if ignore=1.
%               Exclusion threshold associated to the uppper
%               percentage point of the F distribution of Cp which defines
%               the threshold for each p declaring models as irrelevant.
%               The default value of ExclThresh is 0.99999 that is the
%               models whose minimum value of Cp in the part of the
%               search defined by fin_step is above ExclThresh are
%               stored for each p. If option ignore=1, the submodels with
%               p-1 explanatory variables which are contained inside the
%               models considered irrelevant are not considered
%               Example - 'ExclThresh',0.9
%               Data Types - double
%
%     meanmed : Boxes of tha candles. Scalar. It specifies how to construct
%               the boxes of the candles.
%               If meanmed=1 boxes are constructed using mean and median
%               else using the first and third quartile.
%               Example - 'meanmed',1
%               Data Types - double
%
%       plots :  Plot on the screen. Scalar.
%               If plot==1 a candlestick Cp plot is created on the screen
%               else (default) no plot is shown on the screen
%               The options below only work when plots=1
%               Example - 'plots',1
%               Data Types - double
%
%          rl : spread of the candles around
%               each integer value defining the size of the submodels.
%               Scalar.
%               For example if rl=0.4 for each smallp candles are spread in the
%               interval [smallp-rl smallp+rl]. The default value of rl
%               is 0.4. rl does not have to be greater than 0.45 otherwise
%               the candles overlap
%               Example - 'rl',0.3
%               Data Types - double
%
%     quant   :  quantiles for the horizontal lines
%               associated with the confidence bands of Cp. Vector.
%               The default is to plot 2.5% and
%               97.5% envelopes. In other words the default is
%               quant=[0.025;0.975];
%               Example - 'quant',[0.01;0.99]
%               Data Types - double
%
%  CandleWidth:  width of the boxes associated with
%               the central part of the search. Scalar.
%               The default width is 0.05;
%               Example - 'CandleWidth',0.01
%               Data Types - double
%
%   LineWidth : Line Width (in points) for the vertical lines outside the 
%               boxes of the candles. Scalar.
%               The default LineWidth is 0.5 points.
%               Example - 'LineWidth',0.01
%               Data Types - double
%
%     ylimy   : minimum and maximum
%               on the y axis. Vector.
%               Default value is [-2 50] (automatic scale)
%               Example - 'ylimy',[0 10]
%               Data Types - double
%
%     xlimx   :  minimum and maximum
%               on the x axis. Vector.
%               Default value is '' (automatic scale)
%               Example - 'xlimy',[0 10]
%               Data Types - double
%
%
%
% Output:
%
%         outms:   structure which contains the following fields
%
%        outms.stor = k x 9 matrix containing statistics which can be used to create the candles
%               1st col: max Cp values; 
%               2nd col: min Cp values; 
%               3rd col: average Cp values; 
%               4nd col: median Cp values; 
%                   Remark: the information in the first 4 columns is
%                   referred to the central part of the search. 
%               5th col: x coordinates (or size of submodel); 
%               6th col: number of explanatory variables of the submodel
%               7th col: y coordinate of final Cp; 
%               8th col: units entering the final step of the search; 
%               9th col: maximum Cp value during the (central and final
%               part of the) search (This information is used to print the
%               labels on top of each model). 
%        outms.outl = r x 4 matrix containing information about 'influential
%               units' or empty matrix. 
%               Influential units in this context are defined as the units
%               which enter the subset in the final part of the search and
%               bring the value of Cp below the minimum or above the
%               maximum value of the central part of the search. 
%               1st col: x coordinates; 
%               2nd col: y coordinates; 
%               3rd col: step of entry into subset; 
%               4nd col: unit number. 
%               If matrix outl contains more columns they are ignored. 
%        outms.siz  = vector of length 2 containing information about n (number of
%               units of the sample and bigP, number of explanatory
%               variables, including the constant, in the full model). This
%               information is necessary to compute the envelopes. 
%         outms.MAL = (n-init+1) x (k+1) matrix
%                 Mallows Cp monitored along the search for the selected
%                 models. 
%                   1st col is fwd search index; 
%                   2nd col is associated with first selected model; 
%                   3rd col is associated with second selected model; 
%                   ............; 
%                   (k+1)th col is associated with k-th selected model
%                   Notice that k<=(n choose smallp) and that all
%                   models contain the constant. 
%       outms.LAB   =    cell array of strings of length k containing the labels
%                   of the models which have been extracted. First element
%                   of LAB is associated with second column of matrix
%                   MAL ... last element of LAB is associated with last
%                   column of matrix MAL
%
%       Remark: the loop through all values of smallp works backwards in
%       the sense that first all possible submodels of size bigP-1 are
%       considered, then the models with size bigP-2 ...
%
%
%
% See also FSRcp.m, cdsplot.m
%
% References:
%
%   Riani M. and Atkinson A.C. (2010), Robust Model Selection with Flexible Trimming,
%   "Computational Statistics and Data Analysis", Vol. 54, p. 3300-3312.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRms')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{  
    % FSRms with all default options.
    % Common part to all examples: load Ozone dataset, tranform the 
    % response using logs and add a time trend.

    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    [Cpms]=FSRms(y,X);
%}

%{
    %% FSRms with optional arguments.
    % Perform robust model selection and show the generalized candlestick
    % plot.

    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpms]=FSRms(y,X,'labels',labels,'plots',1);
%}

%{
    %% Reproduce the candlestick plot given in Figure 5 of Riani and Atkinson (2010).
    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    n=length(y);
    fin_step=floor([n*0.1 n*0.02]);
    outms=FSRms(y,X,'fin_step',fin_step,'plots',1,'labels',labels,'smallpint',[4:7]);

    %The figure has slightly changed and certainly there can be some random
    %fluctuations due to the number of subset which have been used to initialize
    %the search for each model. However, The indication of the previous two
    %Figures does not change at all: the values of smallp of 4 or 5 should yield
    %a satisfactory model. For smallp = 4 the best model has the trend, x3 and
    %x7, although the plot shows the values of Cp(m) increasing towards the end
    %of the search. By far the most stable model for smallp= 5 adds x2 to these
    %variables.
%}
%{
    % Considering all submodels.
    % Perform robust model selection and show the generalized candlestick
    % plot considering all submodels for each smallp from 2 to size(X).

    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpms]=FSRms(y,X,'labels',labels,'ignore',0,'plots',1);
%}

%{
    %% Comparing results of different settings.
    % Perform robust model selection and show the generalized candlestick
    % plot. Restric attention to the models with size in the interval 4:6
    % Compare the results using ignore=1 with those with ignore=0
    % default option ignore=1.

    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpms]=FSRms(y,X,'smallpint',4:6,'labels',labels,'plots',1);
    % ignore=0
    [Cpms]=FSRms(y,X,'ignore',0,'smallpint',4:6,'labels',labels,'plots',1);
    % with ignore=0 but changing the threshold for excluding models
    [Cpms]=FSRms(y,X,'smallpint',4:6,'labels',labels,'plots',1,'ExclThresh',0.99999999999999);
%}

%{
    % Changing confidence bands.
    % Same options as before but using different confidence bands.

    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    qu=[0.01 0.5 0.99];
    [Cpms]=FSRms(y,X,'smallpint',4:6,'labels',labels,'plots',1,'quant',qu);
%}

%{
    % Personalized LineWidth and CandleWidth.

    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    LineWidth=2;
    CandleWidth=0.03;
    [Cpms]=FSRms(y,X,'smallpint',4:6,'labels',labels,'plots',1,'LineWidth',LineWidth,'CandleWidth',CandleWidth);
%}

%{
    % Input fin_step supplied as fraction (1).
    % For example when fin_step=[0.3 0.1] the central part of the search
    % goes from m=round(n*0.7)=56 to m=round(n*0.9)=72 and the final part
    % of the search goes from m=73 to m=80.
    
    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpms]=FSRms(y,X,'smallpint',4:6,'labels',labels,'plots',1,'fin_step',[0.3 0.1]);
%}

%{
    % Input fin_step supplied as fraction (2).
    % For example when fin_step=[0.36 0.06] the central part of the search
    % goes from m=round(n*0.64)=51 to m=round(n*0.94)=75 and the final part of the search goes from
    % m=76 to m=80.

    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpms]=FSRms(y,X,'smallpint',4:6,'labels',labels,'plots',1,'fin_step',[0.36 0.06]);
%}

%{
    % Input fin_step supplied as integers.
    % For example when fin_step=[20 5] the central part of the search
    % goes from m=n-20=61 to m=n-5=75 and the final part of the search goes from
    % m=76 to m=80.
    % It is worthwhile to notice that independently on how fin_step is
    % chosen, the message of the generalized candlestick plot remains the
    % same. In other words, the best two models with 5 variables are always
    % (Time,4,5,6) and (Time,2,4,5)
    % while two reasonable models with 6 variables are (Time,2,4,5,6) and
    % (Time,2,3,4,5).

    X=load('ozone.txt');
    X(:,end)=log(X(:,end));
    X=[(-40:39)' X];
    y=X(:,end);
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpms]=FSRms(y,X,'smallpint',4:6,'labels',labels,'plots',1,'fin_step',[25 5], 'CandleWidth',0.01);
%}

%% Beginning of code

% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

smallpdef=p-1:-1:2;
nsampdef=1000;
hdef=floor(0.5*(n+p+1));
if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));  
end
options=struct('nocheck',0,'intercept',true,'h',hdef,...
    'nsamp',nsampdef,'lms',1,'init',init,'meanmed',1,...
    'smallpint',smallpdef,'labels','','fin_step',round(n*[0.2 0.05]),'first_k',3,...
    'ignore',1,'ExclThresh',0.99999,...
    'plots',0,'rl',0.4,'quant',[0.025 0.975],'LineWidth',0.5,'CandleWidth',0.05,'xlimx','','ylimy','');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRms:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
nsamp=options.nsamp;
% Information about plotting or not the trajectories
plots=options.plots;

% Information about the labels of the variables
labels=options.labels;

if isempty(labels)
    % Default labels for the explanatory variables
    labels={'1','2','3','4','5','6','7','8','9','A','B','C',...
        'D','E','F','G','H','I','J','K','L','M','N','O',...
        'P','Q','R','S','T','U','V','X','Y','Z'};
    labels=labels(1:p-1);
end

% Information about model exclusion
ignore=options.ignore;
% Threshold for excluding model (only if ignore=1)
ExclThresh=options.ExclThresh;

% outl=matrix which will contain information about the 'influential units'
outl=zeros(100,6);
ii=1;

% ini = scalar which initializes the rows of matrces stor, storLAB and
% the columns of matrix storFWD
% ini is a guess about the number of searches which will be selected
ini=500;

% matrix which will contain information about the selected searches
% 1st col = max values
% 2nd col = min
% 3rd col = mean
% 4th col = median
% 5th col = xcoordinates
% 6th col = small p
% 7th col = final value of Cp
% 8th col = unit which entered the last step of the search
% 9th col = overall maximum (central +final part +final step)
stor=zeros(ini,9);
jj=1;

% storLAB = vector which contains the labels associated to the rows of
% matrix stor
% stor has as many rows as matrix storLAB
storLAB=cell(ini,1);

% storLAB matrix which will contain in the columns the selected searches
% Number of columns of matrix storFWD is equal to the number of rows of
% matrices stor and storLAB
storFWD=zeros(n-3*p,ini);


% Tolerance for candles
% For example if rl=0.4 for each small p candles are spread in the
% interval [p-rl p+rl]
rl=options.rl;

% specfy how to construct the boxes of the candles
% if meanmed=1 boxes are constructed using mean and median
% else we use first and third quartile
meanmed=options.meanmed;


% xlimp Set of values of small p for which to calculate Cp
xlimp=options.smallpint;
xlimp=sort(xlimp,'descend');

% quant contains the quantiles which will be used in the candlestick Cpplot
quant=options.quant;


% Information about the part of the search which has to be analyzed
fin_step=options.fin_step;

% Check if fin_step has been supplied in terms of fractions or of steps
if min(fin_step)<1
    fin_step=round(n*fin_step);
end

% Information about the number of best searches which have to be analyzed
first_k=options.first_k;

% Initalize matrix Excl (the matrix which contains the irrelevant models)
Excl='';

%% Beginning of procedure (loop through all values of smallp)

% Loop through all values of smallp
for smallp=xlimp
    
    % Note discuss 
    outCp=FSRcp(y,X,smallp,'h',h,'nsamp',nsamp,'lms',lms,'init',init,'nocheck',1,...
        'labels',labels,'fin_step',fin_step(1),'first_k',first_k,'Excl',Excl,'ExclThresh',ExclThresh);
    
    % If outCp is an empty structure there no submodels of interest
    if isempty(outCp)
        break;
    end
    
    if ignore
        Excl=outCp.Ajout;
    else
        Excl='';
    end
    
    
    % Set of matrices Un for selected searches
    Unsel=outCp.UnAll;
    
    % Transform cell Unsel in a big unique matrix
    Unall=cell2mat(Unsel);
    
    % Find last unit which entered each search
    lastentry=Unall(Unall(:,1)==n,2);
    
    
    MAL=outCp.MAL;
    lsto=size(MAL(:,2:end),2);
    % lsto = number of selected searches
    
    % Now store required values for robust Cp plot
    % Extract the central part of selected searches in matrix MALcent
    MALcent=MAL((end-fin_step(1)):(end-fin_step(2)),2:end);
    
    % Store max,min,mean and median for selected part of the search
    % Each row refers to a different search
    if meanmed==1
        stat=vertcat(max(MALcent),min(MALcent),mean(MALcent),median(MALcent))';
    else
        stat=vertcat(max(MALcent),min(MALcent),quantile(MALcent,0.75),quantile(MALcent,0.25))';
    end
    
    % Create sequence to spread values
    if lsto==1
        seqr=smallp;
    elseif lsto==2 && rl>0.2
        seqr=vertcat(smallp-0.2,smallp+0.2);
    else
        seqr=(smallp+(-rl:2*rl/(lsto-1):rl))';
    end
    
    % Store: statistics, spread xcoordinates, smallp,  final values of Cp,
    % the unit which entered the last step of the search
    % and in the last column overall maximum (from step n-fin_step(1)+1 to
    % n)
    stor(jj:jj+lsto-1,:)=[stat seqr smallp*ones(lsto,1) MAL(end,2:end)' lastentry max(MAL((end-fin_step(1)):end,2:end))'];
    %disp(stor(jj:jj+lsto-1,:))
    %pause;
    % Store labels associated to each search which is stored

   storLAB(jj:jj+lsto-1)=outCp.LAB;
    
    % Store selected searches in big matrix storFWD
    storFWD(:,jj:jj+lsto-1)=MAL(:,2:end);
    
    
    jj=jj+lsto;
    
    % Now for each selected search compute the 'influential units'
    % i.e. the values in the final part of the search (excluding last
    % value) which are below min or above max
    % Now store required values for robust Cp plot
    MALfin=MAL((end-fin_step(2)+1):(end-1),:);
    
    for j=2:(lsto+1)
        % Extract values above max (>stat(j-1,1)) or below min (stat(j-1,2))
        selj=MALfin(MALfin(:,j)>stat(j-1,1) | MALfin(:,j)<stat(j-1,2),[1 j]);
        
        % If selj is not empty it means that in the final part of the search
        % there have been values smaller than min or grater than max
        % In this case store for each outlier for the Cp plot
        % the x coordinates of the outliers,
        % the y coordinates of the outliers,
        % the step the search in which the outliers took place
        % the associated unit which entered the subset at that particular step
        % the value of small p
        if ~isempty(selj)
            rselj=size(selj,1);
            Unselj=Unsel{j-1};
            
            % Find the required row(s) of matrix Unselj
            % [c, ia, ib] = intersect(Unselj(:,1),selj(:,1));
            [~, ia] = intersect(Unselj(:,1),selj(:,1));
            
            
            % xcoord                             y coordinates  steps units             smallp
            outl(ii:ii+rselj-1,:)=[repmat(seqr(j-1),size(selj,1),1) selj(:,[2 1]) Unselj(ia,2:3) repmat(smallp,size(selj,1),1)];
            ii=ii+rselj;
            
        end
        
    end
    disp(['small p=' int2str(smallp)]);
    
    
end

if ii>1
    outl=outl(1:ii-1,:);
end

if jj>1
    stor=stor(1:jj-1,:);
    storLAB=storLAB(1:jj-1,:);
    storFWD=storFWD(:,1:jj-1);
end

%% Draw the figures
outms=struct;


% Store information about the number of units and the number of bigP
% variables
outms.siz=[n p];
% Store the searches of the best models
outms.MAL=storFWD;
% Store the labels of the best models
outms.LAB=storLAB;
% Store the influential units of the best models
outms.outl=outl;
% Store all the quantities necessary to draw the generalized candlestick
% plot of the best models
outms.stor=stor;

if plots
    % Width of the candles
    width=options.CandleWidth;
    % Width of the vertical lines
    LineWidth=options.LineWidth;
    
    %% Draw the figures
    figure;
    hold on
    hi=stor(:,1);
    lo=stor(:,2);
    av=stor(:,3);
    me=stor(:,4);
    
    index=stor(:,5)';
    
    avs = get(gca, 'colororder');
    color = avs(1, :);
    back = get(gca, 'color');
    
    m = jj-1; % length(hi(:));
    
    % Need to pad all inputs with NaN's to leave spaces between day data
    tmp = nan;
    nanpad = tmp(1, ones(1, m));
    hilo = [hi'; lo'; nanpad];
    
    % Plot lines between high and low value of Cp for each model
    indhilo = index(ones(3, 1), :);
    plot(indhilo(:), hilo(:), 'color', color,'LineWidth',LineWidth)
    
    
    avpad = [av(:)';nanpad];
    avpad = avpad(:)';
    mepad = [me(:)'; nanpad];
    mepad = mepad(:)';
    
    % Create boundaries for filled regions
    xbottom = index - width;
    xbotpad = [xbottom(:)'; nanpad];
    xbotpad = xbotpad(:)';
    xtop = index + width;
    xtoppad = [xtop(:)'; nanpad];
    xtoppad = xtoppad(:)';
    ybottom = min(avpad, mepad);
    ytop = max(avpad, mepad);
    
    % If the median is less than the average, box is empty
    i = find(mepad(:) <= avpad(:));
    boxes(i) = patch([xbotpad(i); xbotpad(i); xtoppad(i); xtoppad(i)],...
        [ytop(i); ybottom(i); ybottom(i); ytop(i)],...
        back, 'edgecolor', color);
    
    % If the median price is greater than the average, box is filled
    i = find(mepad(:) > avpad(:));
    boxes(i) = patch([xbotpad(i); xbotpad(i); xtoppad(i); xtoppad(i)],...
        [ytop(i); ybottom(i); ybottom(i); ytop(i)],...
        color, 'edgecolor', color); %#ok
    
    
    % Add the outliers
    plot(outl(:,1),outl(:,2),'*r');
    
    % Plot final value as a blue circle
    plot(stor(:,5),stor(:,7),'ob');
    
    % Superimpose bands based on the F distribution
    smallp=unique(stor(:,6));
    Env=zeros(length(smallp)*3-1,length(quant)+1);
    Env(:)=NaN;
    % ienv = index which is linked to the rows of matrix Env
    ienv=1;
    
    for i=1:length(smallp)
        ast=finv(quant,p-smallp(i),n-p);
        ast=ast*(p-smallp(i))+2*smallp(i)-p;
        
        Env(ienv,:)=horzcat(smallp(i)-0.5,ast);
        Env(ienv+1,:)=horzcat(smallp(i)+0.5,ast);
        ienv=ienv+3;
    end
    
    plot(Env(:,1),Env(:,2:end),'b','tag','bar');
    
    
    % Add vertical lines in order to separate the different values of p using black color
    xl1=min(stor(:,6))-0.5;
    xl2=max(stor(:,6))+0.5;
    seql=xl1:xl2;
    ons=ones(1,length(seql));
    ycoord=vertcat(-100*ons,100*ons);
    xcoord=repmat(seql,2,1);
    % plot the vertical lines
    line(xcoord,ycoord,'Color','k','tag','bar');
    
    ylimy=options.ylimy;
    
    % Define limits for y axis
    if isempty(ylimy)
        ylim([0 min(3*p,max(stor(:,9)))]);
    else
        ylim(ylimy);
    end
    
    % Define limits for x axis
    xlimx=options.xlimx;
    if isempty(xlimx)
        xlim([xl1 xl2]);
    else
        xlim(xlimx);
    end
    
    % Write quantiles which are used to compute the envelopes
    %  text(repmat(xlim2+0.05,1,2),Env(end,2:end),num2str(quant'));
    
    % Add labels associated with the different selected models
    % 9th col of matrix stor contains overall maximum of Cp values
    text(stor(:,5),stor(:,9)+0.5,storLAB,'Rotation',90,'FontSize',14);
end

end
%FScategory:REG-ModelSelection