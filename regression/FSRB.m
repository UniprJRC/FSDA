function [out]=FSRB(y,X,varargin)
%FSRB gives an automatic outlier detection procedure in Bayesian linear regression
%
%<a href="matlab: docsearch('frsb')">Link to the help function</a>
%
% Required input arguments:
%
%    y: A vector with n elements that contains the response variable. y can
%       be both a row of column vector.
%    X: Data matrix of explanatory variables (also called 'regressors') of
%       dimension (n x p-1). Rows of X represent observations, and columns
%       represent variables.
%       Missing values (NaN's) and infinite values (Inf's) are allowed,
%       since observations (rows) with missing or infinite values will
%       automatically be excluded from the computations.
%
% Optional input arguments:
%
%   intercept   : If 1, a model with constant term will be fitted (default),
%                 if 0, no constant term will be included.
%    bayes      : a structure which specifies prior information
%               Strucure bayes contains the following fields
%               beta0:  p-times-1 vector containing prior mean of \beta
%               R    :  p-times-p positive definite matrix which can be
%                       interepreted as X0'X0 where X0 is a n0 x p matrix
%                       coming from previous experiments (assuming that the
%                       intercept is included in the model
%
%               The prior distribution of tau0 is a gamma distribution with
%               parameters a and b, that is
%                     p(tau0) \propto \tau^{a0-1} \exp (-b0 \tau)
%                         E(tau0)= a0/b0
%
%               tau0 : scalar. Prior estimate of tau=1/ \sigma^2 =a0/b0
%               n0   : scalar. Sometimes it helps to think of the prior
%                      information as coming from n0 previous experiments.
%                      Therefore we assume that matrix X0 (which defines
%                      R), was made up of n0 observations.
%              REMARK: if structure bayes is not supplied the default
%                      values which are used are
%                      beta0= zeros(p,1)  % vector of zeros
%                      R=eye(p);          % Identity matrix
%                      tau0=1/1e+6;       % Very large value for the
%                                         % prior variance, that is a very
%                                         % small value for tau0
%                      n0=1;              % just one prior observation
%
%
%       plots   : Scalar.
%                 If plots=1 (default) the plot of minimum deletion
%                 residual with envelopes based on n observations and the
%                 scatterplot matrix with the outliers highlighted is
%                 produced.
%                 If plots=2 the user can also monitor the intermediate
%                 plots based on envelope superimposition.
%                 else no plot is produced.
%       init    : scalar which specifies the initial subset size to start
%                 monitoring exceedances of minimum deletion residual, if
%                 init is not specified it set equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%       nocheck : Scalar. If nocheck is equal to 1 no check is performed on
%                 matrix y and matrix X. Notice that y and X are left
%                 unchanged. In other words the additional column of ones
%                 for the intercept is not added. As default nocheck=0.
%    bivarfit : This option adds one or more least square lines, based on
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
%                 bivarfit = 'i1' or 'i2' or 'i3' etc.
%                   fits an ols line to a specific group, the one with
%                   index 'i' equal to 1, 2, 3 etc. Again, useful in case
%                   of mixtures.
%       multivarfit : This option adds one or more least square lines, based on
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
%      labeladd : If this option is '1',  we label the outliers with the
%                 unit row index in matrices X and y. The default value is
%                 labeladd='', i.e. no label is added.
%       nameX  :  cell array of strings of length p containing the labels of
%                 the variables of the regression dataset. If it is empty
%                 (default) the sequence X1, ..., Xp will be created
%                 automatically
%       namey  :  character containing the label of the response
%       ylim   :  vector with two elements controlling minimum and maximum
%                 on the y axis. Default value is '' (automatic scale)
%       xlim   :  vector with two elements controlling minimum and maximum
%                 on the x axis. Default value is '' (automatic scale)
%      bonflev  : option to be used if the distribution of the data is
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
%       msg    :  scalar which controls whether to display or not messages
%                 on the screen
%                 If msg==1 (default) messages are displayed on the screen about
%                   step in which signal took place and ....
%                 else no message is displayed on the screen
%
%
% Output:
%
%  The output consists of a structure 'out' containing the following fields:
% out.ListOut=  k x 1 vector containing the list of the units declared as
%               outliers or NaN if the sample is homogeneous
% out.mdrB    =  (n-init) x 2 matrix
%               1st col = fwd search index
%               2nd col = value of Bayesian minimum deletion residual in
%               each step of the fwd search
% out.Un     =  (n-init) x 11 Matrix which contains the unit(s) included
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
%               certain step a non singular matrix X. In this case the
%               search run in a constrained mode, that is including the
%               units which produced a singular matrix in the last n-constr
%               steps. out.constr is a vector which contains the list of
%               units which produced a singular X matrix
%
% See also: FSRBmdr, FSR.m
%
% References:
%
% Chaloner and Brant (1988). A Bayesian Approach to Outlier Detection and
% Residual Analysis, Biometrika, Vol 75 pp. 651-659.
% Riani M., Corbellini A., Atkinson A.C. (2015), Very Robust Bayesian
% Regression for Fraud Detection, submitted
% Atkinson A.C., Corbellini A., Riani M., (2015) Robust Bayesian
% Regression, submitted
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearch('fsrb')">Link to the help page for this function</a>
%
% Last modified 06-Feb-2015

% Examples:

%{
    % Example of Houses Price
    % load dataset
    load hprice.txt;
    
    % setup parameters
    n=size(hprice,1);
    y=hprice(:,1);
    X=hprice(:,2:5);
    n0=5;

    % set \beta components
    beta0=0*ones(5,1);
    beta0(2,1)=10;
    beta0(3,1)=5000;
    beta0(4,1)=10000;
    beta0(5,1)=10000;

    % \tau
    s02=1/4.0e-8;
    tau0=1/s02;

    % R prior settings
    R=2.4*eye(5);
    R(2,2)=6e-7;
    R(3,3)=.15;
    R(4,4)=.6;
    R(5,5)=.6;
    R=inv(R);

    % define a Bayes structure with previous data
    bayes=struct;
    bayes.R=R;
    bayes.n0=n0;
    bayes.beta0=beta0;
    bayes.tau0=tau0;
    intercept=1;

    % function call
    outBA=FSRB(y,X,'bayes',bayes,'msg',0,'plots',1,'init',round(n/2),'intercept', intercept)
%}

%{
    % Fishery Example with Empirical prior
    close all
    % nsamp = number of subsamples to use in the frequentist analysis of first
    % year, in order to find initial subset using LMS
    nsamp=3000;
    % threshold to be used to increase susbet of good units
    threshold=300;
    bonflev=0.99; % Bonferroni confidence level to be used for first year
    bonflevB=0.99; % Bonferroni confidence level to be used for subsequent years

    % Load 2002 Fishery dataset
    Fishery2002=load('Fishery2002.txt');
    y02=Fishery2002(:,3);
    X02=Fishery2002(:,2);

    n02=length(y02);
    seq02=1:n02;

    % frequentist Forward Search, 1st year
    [out02]=FSR(y02,X02,'nsamp',nsamp,'plots',1,'msg',0,'init',round(n02*3/4),'bonflev',bonflev);

    % In what follows
    % g stands for good units
    % i stand for intermediate units (i.e. units whose raw residual is smaller
    % than threshold)
    % o stands for outliers
    % gi stands for good +intermediate units

    % u02g = good units
    % n02g = number of good units
    u02g=setdiff(seq02,out02.ListOut);
    n02g=length(u02g);

    X02g=[ones(length(u02g),1) X02(u02g,:)];
    y02g=y02(u02g);
    % b02g = regression coefficients just using g units
    b02g=X02g\y02g;
    % res02 = squared raw residuals for all units using b02g
    res02=(y02-[ones(length(X02),1) X02]*b02g).^2;
    res02o=res02(out02.ListOut);
    % sel= boolean vector which is true for the intermediate units 
    % (units whose squared residual is below the threshold)
    sel=res02o<threshold^2;
    % u02i = vector containing intermediate units (that is outliers whose
    % residual is smaller than threshold)
    u02i=out02.ListOut(sel);
    % u02o = vector containing outliers whose residual is out of the threshold
    u02o=out02.ListOut(~sel);
    % u02gi = g + i units
    if ~isempty(u02i)
        u02gi=[u02g u02i];
    else
        u02gi=u02g;
    end
    % n02gi = number of good + intermediate units
    n02gi=length(u02gi);

    % plotting section
    hold('off')
    % good units, plotted as (+)
    plot(X02(u02g)',y02(u02g)','Marker','+','LineStyle','none','Color','b')
    hold('on')
    % intermediate units plotted as (X)
    plot(X02(u02i)',y02(u02i)','Marker','X','MarkerSize',9,'LineWidth',2,'LineStyle','none','Color','m')

    % outliers, plotted as (O)
    plot(X02(u02o)',y02(u02o)','Marker','o','LineStyle','none','Color','r')
    xlabel('Quantity');
    ylabel('Price');
    title('Frequentist - 2002');

    % S202gi = estimated of sigma^2 using g+i units
    S202gi=sum(res02(u02gi))/(n02gi-2);


    % X02gi = X matrix referred to good + intermediate units
      X02gi=[ones(n02gi,1) X02(u02gi,:)];
    % y02gi = y vector referred to good + intermediate units
      y02gi=y02(u02gi);

    % bayes = structure which contains prior information to be used in year
    % 2003
    bayes=struct;
    bayes.beta0=b02g; % beta prior is beta based on g units
    tau0=1/S202gi; % tau0 is based on g + i units
    bayes.tau0=tau0;  
    R=X02g'*X02g; % R is based on g units
    bayes.n0=n02gi; % n0 is based on g + i units
    bayes.R=R;

    %% 2003 
    % Load 2003 Fishery dataset
    Fishery2003=load('Fishery2003.txt');
    y03=Fishery2003(:,3);
    X03=Fishery2003(:,2);
    n03=length(y03);
    seq03=1:n03;

    % Bayesian Forward Search, 2nd year
    out03=FSRB(y03,X03,'bayes',bayes,'msg',0,'plots',1,'init',round(n03/2),'bonflev',bonflevB);

    u03g=setdiff(seq03,out03.ListOut);
    n03g=length(u03g);

    % compute beta coefficient for year 2003 just using good units
    X03g=[ones(n03g,1) X03(u03g,:)];
    y03g=y03(u03g);
    b03g=X03g\y03g;

    res03=(y03-[ones(length(X03),1) X03]*b03g).^2;
    res03o=res03(out03.ListOut);
    sel=res03o<threshold^2;
    % u03i = units to add to the good units subset (intermediate units)
    u03i=out03.ListOut(sel);
    % u03o =  outliers out of the threshold
    u03o=out03.ListOut(~sel);
    if ~isempty(u03i)
        u03gi=[u03g u03i];
    else
        u03gi=u03g;
    end
    n03gi=length(u03gi);
    X03gi=[ones(n03gi,1) X03(u03gi,:)];
    y03gi=y03(u03gi,:);


    % plotting section
    hold('off')
    % good units, plotted as (+)
    plot(X03(u03g)',y03(u03g)','Marker','+','LineStyle','none','Color','b')
    hold('on')
    % units below the threshold, plotted as (X)
    plot(X03(u03i)',y03(u03i)','Marker','X','MarkerSize',9,'LineWidth',2,'LineStyle','none','Color','m')

    % outliers, plotted as (O)
    plot(X03(u03o)',y03(u03o)','Marker','o','LineStyle','none','Color','r')
    set(gca,'FontSize',14)
    xlabel('QUANTITY in tons')
    ylabel('VALUE in 1000 euro')
    title('2003');

    % Definition of bayes structure (based on 2002 and 2003)
    bayes=struct;
    X02gX03g=[X02g; X03g];
    y02gy03g=[y02g; y03g];
    n02gn03g=n02g+n03g;
    % b0203g prior estimate of beta for year 2004 is computed using good units
    % for years 2002 and 2003
    b0203g=X02gX03g\y02gy03g;

    bayes.beta0=b0203g;
    % R is just referred to good units for years 2002 and 2003
    R=X02gX03g'*X02gX03g;
    bayes.R=R;
    % n0 is referred to g + i units in 2002 and 2003
    bayes.n0=n02gi+n03gi;

    X02giX03gi=[X02gi; X03gi];
    y02giy03gi=[y02gi; y03gi];
    % n02gin03gi = number of g+i units in 2002 and 2003 
    n02gin03gi=n02gi+n03gi;

    % res = residuals for g+i units using b0203g
    res=y02giy03gi-X02giX03gi*b0203g;
    S203gi=sum(res.^2)/(n02gin03gi-2);
    % estimate of tau is based on g + i units
    tau0=1/S203gi;
    bayes.tau0=tau0;


    %% Load 2004 Fishery dataset
    Fishery2004=load('Fishery2004.txt');
    y04=Fishery2004(:,3);
    X04=Fishery2004(:,2);
    n04=length(y04);
    seq04=1:n04;

    % Bayesian Forward Search, 3rd year
    out04=FSRB(y04,X04,'bayes',bayes,'msg',0,'plots',1,'init',round(n04/2),'bonflev',bonflevB);

    u04g=setdiff(seq04,out04.ListOut);
    n04g=length(u04g);

    X04g=[ones(n04g,1) X04(u04g,:)];
    y04g=y04(u04g);
    % b04g = beta based on good units for year 2004
    b04g=X04g\y04g;

    res04=(y04-[ones(length(X04),1) X04]*b04g).^2;
    % res04o squared residuals for the tentative outliers
    res04o=res04(out04.ListOut);
    % we keep statistical units below the treshold
    sel=res04o<threshold^2;
    % u04i = units to add to the good units subset (intermediate units)
    u04i=out04.ListOut(sel);
    % u04o = units outliers out of the threshold
    u04o=out04.ListOut(~sel);
    if ~isempty(u04i)
        u04gi=[u04g u04i];
    else
        u04gi=u04g;
    end
    n04gi=length(u04gi);

    % plotting section
    hold('off')
    % good units, plotted as (+)
    plot(X04(u04g)',y04(u04g)','Marker','+','LineStyle','none','Color','b')
    hold('on')
    % units below the treshold, plotted as (X)
    plot(X04(u04i)',y04(u04i)','Marker','X','MarkerSize',9,'LineWidth',2,'LineStyle','none','Color','m')

    % outliers, plotted as (O)
    plot(X04(u04o)',y04(u04o)','Marker','o','LineStyle','none','Color','r')
    set(gca,'FontSize',14)
    xlabel('QUANTITY in tons')
    ylabel('VALUE in 1000 euro')


    %% frequentist Forward Search, 3rd year
    out04=FSR(y04,X04,'nsamp',nsamp,'plots',1,'msg',0,'init',round(n04/2),'bonflev',bonflev);
    xlabel('QUANTITY in tons')
    ylabel('VALUE in 1000 euro')
    title('Frequentist - 2004');

%}


%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));
end

% ini0=init;
bayesdef=struct;
bayesdef.beta0=zeros(p,1);
bayesdef.R=eye(p);
bayesdef.tau0=1/1e+6;
bayesdef.n0=1;


options=struct('plots',1,'init',init,...
    'labeladd','','bivarfit','','multivarfit','',...
    'xlim','','ylim','','nameX','','namey','','msg',1, ...
    'nocheck',0,'intercept',1,'bonflev','', 'bayes',bayesdef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRB:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


% Ovewrite (if provided) in structure 'options' the options chosen by the user
% if user provides the variables the rewrite the strucure with such
% variables

if nargin > 2
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end


init=options.init;
bayes=options.bayes;

beta0=bayes.beta0;
R=bayes.R;
tau0=bayes.tau0;
n0=bayes.n0;


%% Start of the forward search
[mdrB,Un,bb,~,~] = FSRBmdr(y, X, beta0, R, tau0, n0,'nocheck',1);

%% Call core function which computes exceedances to thresholds of mdr
[out]=FSRcore(y,X,n,p,mdrB,init,Un,bb,options);

end

