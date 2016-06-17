function [out , varargout] = FSRBr(y, X, varargin)
%Bayesian forward search in linear regression reweighted
%
%
%<a href="matlab: docsearchFS('FSRBr')">Link to the help function</a>
%
%   FSRBr uses the units not declared as outliers by  FSRB to produce a robust fit.
%   The units whose residuals exceeds the threshold determined by option
%   alpha are declared as outliers. Moreover, it is possible in option
%   R2th to modify the estimate of sigma2 which is used to declare
%   the outliers. This is useful when there is almost a perfect fit in the
%   data, the estimate of the error variance is very small and therefore
%   there is the risk of declaring as outliers very small deviations from
%   the robust fit. In this case the estimate of sigma2 is corrected in
%   order to achieve a value of R2 equal to R2th.
%
%
%  Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :          Predictor variables. Matrix. Matrix of explanatory
%               variables (also called 'regressors') of dimension n x (p-1)
%               where p denotes the number of explanatory variables
%               including the intercept.
%               Rows of X represent observations, and columns represent
%               variables. By default, there is a constant term in the
%               model, unless you explicitly remove it using input option
%               intercept, so do not include a column of 1s in X. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations.
%
%Optional input arguments:
%
%
%       alpha: test size. Scalar. Number between 0 and 1 which
%              defines test size to declare the outliers
%                 Example - 'alpha',0.01
%                 Data Types - double
%       R2th : R2 threshold. Scalar. Scalar which defines the value R2 does
%              have to exceed. For example if R2 based on good observations
%              is 0.92 and R2th is 0.90 the estimate of the variance of the
%              residuals which is used to declare the outliers is adjusted
%              in order to have a value of R2 which is equal to 0.90.
%                 Example - 'R2th',0.95
%                 Data Types - double
%fullreweight: Option to declare outliers. Boolean. If fullreweight is true
%              (default option), the list of outliers refers to all the
%              units whose residuals is above the threshold else if it is
%              false the outliers are the observations which by procedure
%              FSR had been declared outliers and have a residual greater
%              than threshold
%                 Example - 'fullreweight',true
%                 Data Types - boolean
%    plotsPI  : Plot of prediction intervals. Scalar. If plotsPI =1 and
%               the number of regressors (excluding the constant term) is
%               equal 1, it is possible to see on the screen the yX scatter
%               with superimposed the prediction intervals using a
%               confidence level 1-alpha, else no plot is shown on the
%               screen
%                 Example - 'plotsPI',1
%                 Data Types - double
%   intercept   :  Indicator for constant term. Scalar.
%                       If 1, a model with constant term will be fitted
%                       (default), if 0, no constant term will be included.
%                        Example - 'intercept',1
%                       Data Types - double
%    bayes      : Prior information. Structure.
%
%                       It contains the following fields
%               bayes.beta0=  p-times-1 vector containing prior mean of \beta
%               bayes.R    =  p-times-p positive definite matrix which can be
%                       interpreted as X0'X0 where X0 is a n0 x p matrix
%                       coming from previous experiments (assuming that the
%                       intercept is included in the model.
%
%               The prior distribution of $\tau_0$ is a gamma distribution with
%               parameters $a_0$ and $b_0$, that is
%               \[
%                     p(\tau_0) \propto \tau^{a_0-1} \exp (-b_0 \tau)
%                       \qquad E(\tau_0) = a_0/b_0
%               \]
%               bayes.tau0 = scalar. Prior estimate of
%                       \[ \tau=1/ \sigma^2 =a_0/b_0 \]
%               bayes.n0   = scalar. Sometimes it helps to think of the prior
%                      information as coming from n0 previous experiments.
%                      Therefore we assume that matrix X0 (which defines
%                      R), was made up of n0 observations.
%              REMARK if structure bayes is not supplied the default
%                      values which are used are.
%                      beta0= zeros(p,1):  Vector of zeros.
%                      R=eye(p):           Identity matrix.
%                      tau0=1/1e+6:        Very large value for the
%                                          prior variance, that is a very
%                                          small value for tau0.
%                      n0=1:               just one prior observation.
%
%               $\beta$ is assumed to have a normal distribution with
%               mean $\beta_0$ and (conditional on $\tau_0$) covariance
%               $(1/\tau_0) (X_0'X_0)^{-1}$.
%               $\beta \sim N(    \beta_0, (1/\tau_0) (X_0'X_0)^{-1}    )$
%
%                     Example - bayes=struct;bayes.R=R;bayes.n0=n0;bayes.beta0=beta0;bayes.tau0=tau0;
%                     Data Types - double
% plots   :    Plot on the screen. Scalar.
%                 If plots=1 (default) the plot of minimum deletion
%                 residual with envelopes based on n observations and the
%                 scatterplot matrix with the outliers highlighted is
%                 produced.
%                 If plots=2 the user can also monitor the intermediate
%                 plots based on envelope superimposition.
%                 Else no plot is produced.
%                 Example - 'plots',1
%                 Data Types - double
%       init    :  Search initialization. Scalar.
%                   scalar which specifies the initial subset size to start
%                   monitoring exceedances of minimum deletion residual, if
%                   init is not specified it set equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%                   Example - 'init',100 starts monitoring from step m=100
%                   Data Types - double
%   nocheck : Check input arguments. Scalar.
%                    If nocheck is equal to 1 no check is performed on
%                    matrix y and matrix X. Notice that y and X are left
%                    unchanged. In other words the additional column of ones
%                     for the intercept is not added. As default nocheck=0.
%                   Example - 'nocheck',1
%                   Data Types - double
%    bivarfit :  Superimpose bivariate least square lines. Character.
%                   This option adds one or more least square lines, based on
%                   SIMPLE REGRESSION of y on Xi, to the plots of y|Xi.
%                  bivarfit = ''
%                   is the default: no line is fitted.
%                  bivarfit = '1'
%                   fits a single ols line to all points of each bivariate
%                   plot in the scatter matrix y|X.
%                  bivarfit = '2'
%                   fits two ols lines: one to all points and another to
%                   the group of the genuine observations. The group of the
%                   potential outliers is not fitted.
%                  bivarfit = '0'
%                   fits one ols line to each group. This is useful for the
%                   purpose of fitting mixtures of regression lines.
%                  bivarfit = 'i1' or 'i2' or 'i3' etc.
%                   fits an ols line to a specific group, the one with
%                   index 'i' equal to 1, 2, 3 etc. Again, useful in case
%                   of mixtures.
%                 Example - 'bivarfit',2
%                 Data Types - char
%       multivarfit : Superimpose multivariate least square lines. Character.
%                   This option adds one or more least square lines, based on
%                   MULTIVARIATE REGRESSION of y on X, to the plots of y|Xi.
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
%                 Example - 'multivarfit','1'
%                 Data Types - char
%      labeladd : Add outlier labels in plot. Character.
%                 If this option is '1',  we label the outliers with the
%                 unit row index in matrices X and y. The default value is
%                 labeladd='', i.e. no label is added.
%                 Example - 'labeladd','1'
%                 Data Types - char
%       nameX  :  Add variable labels in plot. Cell array of strings.
%                 cell array of strings of length p containing the labels of
%                 the variables of the regression dataset. If it is empty
%                 (default) the sequence X1, ..., Xp will be created
%                 automatically
%                 Example - 'nameX',{'NameVar1','NameVar2'}
%                 Data Types - cell
%       namey  :  Add response label. Character.
%               character containing the label of the response
%               Example - 'namey','NameOfResponse'
%               Data Types - char
%       ylim   :   Control y scale in plot. Vector.
%                   vector with two elements controlling minimum and maximum
%                 on the y axis. Default value is '' (automatic scale)
%               Example - 'ylim','[0,10]' sets the minim value to 0 and the
%               max to 10 on the y axis
%               Data Types - double
%       xlim   :   Control x scale in plot. Vector.
%                  vector with two elements controlling minimum and maximum
%                 on the x axis. Default value is '' (automatic scale)
%               Example - 'xlim','[0,10]' sets the minim value to 0 and the
%               max to 10 on the x axis
%               Data Types - double
%      bonflev  : Signal to use to identify outliers. Scalar.
%                   option to be used if the distribution of the data is
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
%       msg    :  Level of output to display. Scalar.
%               scalar which controls whether to display or not messages
%                 on the screen
%                 If msg==1 (default) messages are displayed on the screen about
%                   step in which signal took place and ....
%                 else no message is displayed on the screen
%               Example - 'msg',1
%               Data Types - double
%
%
% Output:
%
%         out:   structure which contains the following fields
%
% out.outliers=  k x 1 vector containing the list of the units declared
%               outliers by procedure FSR or NaN if the sample is
%               homogeneous
% out.beta   =  p-by-1 vector containing the estimated regression parameter
%               by procedure FSR
% out.outliersr=  k1 x 1 vector containing the list of the units declared
%               outliers after the reweighting step or NaN if the sample is
%               homogeneous
% out.betar  =  p-by-1 vector containing the estimated regression parameter
%               after the reweighting step
% out.rstud =  n-by-2 matrix.
%              First column = studentized residuals
%              Second column = p-values (computed using as reference
%              distribution the Student t)
%
%Optional Output:
%
%           xnew : new points. Vector. Vector with a number of new points where to evaluate the
%                  prediction interval. xnew is a vector.
%          ypred : values predicted by the fitted model on xnew. Vector of
%                  length(xnew)
%           yci  : Confidence intervals. A two-column matrix with each row providing
%                  one interval. 
%
% See also: FSRB, FSR
%
% References:
%
% Atkinson A.C., Corbellini A., Riani M., (2016) Robust Bayesian
% Regression, submitted
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('FSRBr')">Link to the help page for this function</a>
% Last modified 14-06-2016
%
% Examples:
%
%{
        %% Example of FSRB for international trade data.
        % Bayesian FS to fit the group of undervalued flows.
        load('fishery');
        X = fishery.data(:,1);
        y = fishery.data(:,2);
        [n,p] = size(X);
        X = X + 0.000001*randn(n,1);

        % id = undervalued flows
        id = (y./X < 9.5);
        
        % my prior on beta
        mybeta = median(y(id)./X(id));
        bmad   = mad(y(id)./X(id));
        ListIn = find(y./X <= mybeta + bmad);
        % gscatter(X,y,id)


        rr          = y(ListIn)-(X(ListIn,:) * mybeta);
        numS2cl     = rr'*rr;
        dfe         = length(ListIn)-p;
        S2cl        = numS2cl/dfe;
        bayes       = struct;
        bayes.R     = X(ListIn,:)'*X(ListIn,:);
        bayes.tau0  = 1/S2cl;
        bayes.n0    = length(ListIn);
        bayes.beta0 = mybeta;

        % fit based on Bayesian FS with prior on the underdeclared flows
        [out_B, xnew1 , ypred1, yci1]   = FSRBr(y,X,'bayes',bayes,'intercept',0,'alpha',0.01,'bonflev',0.999,'fullreweight',false,'plotsPI',1,'plots',0);
         
        h1 = allchild(gca); a1 = gca; f1 = gcf;

        % fit based on traditional FS
        [out, xnew2 , ypred2, yci2]   = FSRr(y,X,'intercept',0,'alpha',0.01,'bonflev',0.999,'fullreweight',false,'plotsPI',1,'plots',0);

        h2 = allchild(gca); a2 = gca; f2 = gcf;

        % move the figure above into a single one with two panels
        hh = figure; ax1 = subplot(2,1,1); ax2 = subplot(2,1,2);
        copyobj(h1,ax1); title(ax1,get(get(a1,'title'),'string'));
        copyobj(h2,ax2); title(ax2,get(get(a2,'title'),'string'));
        figsize = get(hh, 'Position'); 
        set(hh,'Position',figsize);
        close(f1); close(f2);

        disp(['Bayesian FS fit    = ' num2str(out_B.betar) ' using a prior based on undervalued flows']);
        disp(['Traditional FS fit = ' num2str(out.betar)]);

%}


%{
        % Example of FSRB for international trade data (explore options).
        % Bayesian FS to fit the group of undervalued flows.
        load('fishery');
        X = fishery.data(:,1);
        y = fishery.data(:,2);
        [n,p] = size(X);
        X = X + 0.000001*randn(n,1);

        % id = undervalued flows
        id = (y./X < 9.5);
        
        % my prior on beta
        mybeta = median(y(id)./X(id));
        bmad   = mad(y(id)./X(id));
        ListIn = find(y./X <= mybeta + bmad);
        % gscatter(X,y,id)


        rr          = y(ListIn)-(X(ListIn,:) * mybeta);
        numS2cl     = rr'*rr;
        dfe         = length(ListIn)-p;
        S2cl        = numS2cl/dfe;
        bayes       = struct;
        bayes.R     = X(ListIn,:)'*X(ListIn,:);
        bayes.tau0  = 1/S2cl;
        bayes.n0    = length(ListIn);
        bayes.beta0 = mybeta;
        alpha=0.0001;
        bonflev=0.99999;
        % fit based on Bayesian FS with prior on the underdeclared flows
        [out_B, xnew1 , ypred1, yci1]   = FSRBr(y,X,'bayes',bayes,'intercept',0,'alpha',alpha,'bonflev',bonflev,'fullreweight',false,'plotsPI',1,'plots',0);
         
        h1 = allchild(gca); a1 = gca; f1 = gcf;

        % fit based on traditional FS
        [out, xnew2 , ypred2, yci2]   = FSRr(y,X,'intercept',0,'alpha',alpha,'bonflev',bonflev,'fullreweight',false,'plotsPI',1,'plots',0);

        h2 = allchild(gca); a2 = gca; f2 = gcf;

        % move the figure above into a single one with two panels
        hh = figure; ax1 = subplot(2,1,1); ax2 = subplot(2,1,2);
        copyobj(h1,ax1); title(ax1,get(get(a1,'title'),'string'));
        copyobj(h2,ax2); title(ax2,get(get(a2,'title'),'string'));
        figsize = get(hh, 'Position'); 
        set(hh,'Position',figsize);
        close(f1); close(f2);

        disp(['Bayesian FS fit    = ' num2str(out_B.betar) ' using a prior based on undervalued flows']);
        disp(['Traditional FS fit = ' num2str(out.betar)]);

%}

%% Beginning of code

% The first four options below are specific for this function, all the others
% refer to routine FSRB
n=length(y);
init       = round(n*0.5);

options     = struct('plotsPI',0,'alpha',0.05,'fullreweight',true,'R2th',1,...
    'plots',1,'init',init,...
    'labeladd','','bivarfit','','multivarfit','',...
    'xlim','','ylim','','nameX','','namey','','msg',0, ...
    'nocheck',0,'intercept',1,'bonflev','', 'bayes','');


UserOptions = varargin(1:2:length(varargin));

if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('Error:: number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if all the specified optional arguments were present
    % in structure options
    inpchk       = isfield(options,UserOptions);
    WrongOptions = UserOptions(inpchk==0);
    
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('Error:: in total %d non-existent user options found.', length(WrongOptions));
    end
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for j=1:2:length(varargin);
        options.(varargin{j}) = varargin{j+1};
    end
end

init=options.init;
labeladd=options.labeladd;
bivarfit=options.bivarfit;
multivarfit=options.multivarfit;
xlim=options.xlim;
ylim=options.ylim;
nameX=options.nameX;
namey=options.namey;
msg=options.msg;
nocheck=options.nocheck;
intercept=options.intercept;
bonflev=options.bonflev;
bayes=options.bayes;
plots=options.plots;

[outFSRB]=FSRB(y,X,...
    'plots',plots,'init',init,...
    'labeladd',labeladd,'bivarfit',bivarfit,'multivarfit',multivarfit,...
    'xlim',xlim,'ylim',ylim,'nameX',nameX,'namey',namey,'msg',msg, ...
    'nocheck',nocheck,'intercept',intercept,'bonflev',bonflev, 'bayes',bayes);

% Initialize structure out
out.outliers=outFSRB.ListOut;
out.beta=outFSRB.beta;
beta=outFSRB.beta';
seq=1:n;

% Options specific for this function
alpha=options.alpha;
fullreweight=options.fullreweight;
R2th=options.R2th;
plotsPI=options.plotsPI;

% ListOut vector containing the outliers
if ~isnan(outFSRB.ListOut)
    ListOut=outFSRB.ListOut;
    ListIn=setdiff(seq,ListOut);
else
    ListIn=seq;
    ListOut='';
end
nlistIn=length(ListIn);


% Find S2 using units not declared as outliers using FSR
%in case intercept=1:
if intercept == 1
    X = [ones(n,1) X];
end
% p= number of explanatory variables
p=size(X,2);
% Xb = subset of X referred to good units
Xb=X(ListIn,:);
% res=raw residuals for all observartions
res=y-X*beta;
% resb= raw residuals for good observations
resb=res(ListIn);
% numS2b = numerator of the estimate of the error variance (referred to
% subset)
numS2b=(resb'*resb);
%ytildeb = deviation from the mean (if intercept is present) for subset
if intercept==1
    ytildeb=y(ListIn)-mean(y(ListIn));
else
    ytildeb=y(ListIn);
end

% devtotb = total deviance referred to subset
devtotb=ytildeb'*ytildeb;
% compute R2b = R squared referred to susbet;
R2b=1-numS2b/devtotb;

% Correct the value of the deviance of residuals (numerator of S2) if R2
% is greater than R2th
if R2b >R2th
    numS2b=devtotb*(1-R2th);
end
dfe=nlistIn-p;
S2b=numS2b/dfe;

% studres= vector which will contain squared (appropriately studentized) residuals for all n units.
% For the units non declared as outliers by FS they will be squared studentized
% residuals (that is at the denominator we have (1-h)), while for the units declared
% as outliers by FS, they are deletion residuals (that is at the denominator
% we have (1+h)).
studres2=zeros(n,1);

mAm=Xb'*Xb;

if ~isempty(ListOut)
    % Take units not belonging to bsb
    Xncl = X(ListOut,:);
    % Find leverage for units not belonging to good observations
    % mmX=inv(mAm);
    % hi = sum((Xncl*mmX).*Xncl,2);
    hi=sum((Xncl/mAm).*Xncl,2);
    studres2(ListOut)= ((res(ListOut).^2)./(1+hi));
end
hi=sum((Xb/mAm).*Xb,2);
studres2(ListIn)=((resb.^2)./(1-hi));
studres2=studres2/S2b;

% The final outliers are the units declared as outiers by FSR for which
% observations r(ncl) is greater than the confidence threshold
if fullreweight
    % rncl boolean vector which contains true for the unit whose
    % squared stud residual exceeds the F threshold
    rncl=studres2>finv(1-alpha, 1, dfe);
else
    %rncl=boolean vector which contains true for the units which had
    %been declared as outliers by FSR and whose squared stud residual
    %exceeds the F threshold
    rncl=false(n,1);
    rncl(ListOut)=studres2(ListOut)>finv(1-alpha, 1, dfe);
end
% outliersr = list of units declared as outliers after reweighting step
outliersr=seq(rncl);
out.outliersr=outliersr;

if isequal(out.outliers,outliersr)
    out.betar=outFSRB.beta;
else
    weights=~rncl;
    betar = X(weights,:) \ y(weights);
    out.betar=betar';
end

% Find p-values of studentized residuals

% Store studentized residuals
% and the associated p-values
if verLessThan('matlab','8.3.0')
    rstud=[sign(res).*sqrt(studres2) 1 - fcdf(studres2,1,dfe)];  
else
    rstud=[sign(res).*sqrt(studres2) fcdf(studres2,1,dfe,'upper')];
end
out.rstud=rstud;

if nargout > 0 || plotsPI==1
    
    minX=min(X(:,end));
    maxX=max(X(:,end));
    
    xnew=(minX:((maxX-minX)/1000):maxX)';
    if intercept==1
        xnew=[ones(length(xnew),1) xnew];
        hasintercept=true;
    else
        hasintercept=false;
    end
    % Var cov matrix of regression coefficients
    Sigma = (inv(mAm))*S2b;
    
    sim   = false;
    pred  = true;
    [ypred, yci] = predci(xnew,beta,Sigma,S2b,dfe,alpha,sim,pred,hasintercept);
    
    varargout = {xnew , ypred, yci};
    
    if plotsPI==1
        
        %     PI_LOWER_BOUND = yci(:,1);
        %     PI_UPPER_BOUND = yci(:,2);
        figure('name','Bayesian FS: Prediction Interval');
        hold('on');
        plot(X(:,end),y,'o')
        plot(xnew(:,end),ypred)
        
        plot(xnew(:,end),yci(:,1))
        plot(xnew(:,end),yci(:,2))
        
        if R2th < 1
            if R2b > R2th
                tit2 = ['with variance of residuals adjusted to correct R2 from ' num2str(R2b,4) ' to ' num2str(R2th,4)];
            else
                tit2 = ['variance of residuals not adjusted, as R2=' num2str(R2b,4) ' < R2th=' num2str(R2th,4)];
            end
        else
            tit2 = '';
        end
        title({[num2str((1-alpha)*100,4) '% Prediction Interval of the Bayesian FS'] , tit2});
        
    end
    
end

end

%FScategory:REG-Bayes
