function [out , varargout] = FSRBr(y, X, varargin)
%Bayesian forward search in linear regression reweighted
%
%
%<a href="matlab: docsearchFS('FSRBr')">Link to the help function</a>
%
%   FSRBr uses the units not declared by outliers by  FSRB to produce a robust fit.
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
%           xnew = vector with a number of new points where to evaluate the
%                  prediction interval. xnew is a vector.
%          ypred = values predicted by the fitted model on xnew. Vector of
%                  length(xnew)
%           yci  = Confidence intervals. A two-column matrix with each row providing
%                  one interval. 
%
% See also: FSR
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('FSRBr')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:
%
%{
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
rstud=[sign(res).*sqrt(studres2) fcdf(studres2,1,dfe,'upper')];
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


