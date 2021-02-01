function [out] = mdpdReda(y, X, varargin)
%mdpdReda allows to monitor  Minimum Density Power Divergence criterion to parametric regression problems.
%
%<a href="matlab: docsearchFS('mdpdReda')">Link to the help function</a>
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
%                 Data Types - double
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
%                 Data Types - double
%
%  Optional input arguments:
%
%  tuningpar    : tuning parameter. Scalar or Vector.
%               tuningpar may refer to $\alpha$ (default) or to breakdown
%               point (depending on input option alphaORbdp.
%               As the tuning parameter $\alpha$ (bdp) decreases the robustness
%               of the Minimum Density Power Divergence estimator decreases
%               while its efficiency increases (Basu et al., 1998). For
%               $\alpha=0$ the MDPDE becomes the Maximum Likelihood
%               estimator, while for $\alpha=1$ the divergence yields the
%               $L_2$ metric and the estimator minimizes the $L_2$ distance
%               between the densities, e.g., Scott (2001), Durio and Isaia
%               (2003).  The sequence is forced to be monotonically
%               decreasing, e.g. alpha=[1 0.9 0.5 0.01]. The default for
%               tuningpar is a sequence from 1 to 0 with step -0.01.
%                 Example - 'tuningpar',[1 0.8 0.5 0.4 0.3 0.2 0.1]
%                 Data Types - double
%
%   alphaORbdp  : ctuning refers to $\alpha$ or to breakdown point.
%                 Character. Character which specifies what are the values
%                 in input option tuningpar. If this option is not
%                 specified or it is equal to 'alpha' then program assumes
%                 that the values of tuningpar refer to 'alpha', elseif
%                 this option is equal to 'bdp', program assumes that the
%                 values of tuninpar refer to breakdownpoint.
%                 Example - 'alphaORbdp','bdp'
%                 Data Types - char
%
%   modelfun   : non linear function to use.
%                function_handle or empty value (default). If
%                modelfun is empty the link between $X$ and $\beta$ is assumed
%                to be linear else it is necessary to specify a function
%                (using @) that accepts two arguments, a coefficient vector
%                and the array X and returns the vector of fitted values
%                from the non linear model y. For example, to specify the
%                hougen (Hougen-Watson) nonlinear regression function, use
%                the function handle @hougen.
%                 Example - 'modelfun', modelfun where modelfun = @(beta,X) X*beta(1).*exp(-beta(2)*X);
%                 Data Types - function_handle or empty value
%
%  theta0       :  empty value or vector containing initial values for the
%                 coefficients (beta0 and sigma0) just in case modelfun is
%                 non empty. If modelfun is empty this argument is ignored
%                 and LMS solution will be used as initial solution for the
%                 minimization.
%                 Example - 'beta0',[0.5 0.2 0.1]
%                 Data Types - double
%
%  intercept :  Indicator for constant term. true (default) | false.
%               If true, and modelfun is empty (that is if the link between
%               X and beta is linear) a model with constant term will be
%               fitted (default), else no constant term will be included.
%               This argument is ignored if modelfun is not empty.
%               Example - 'intercept',true
%               Data Types - boolean
%
%     conflev :  Confidence level. Scalar.
%               Confidence level which is used to declare units as outliers.
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975.
%                 Example - 'conflev',0.99
%                 Data Types - double
%
%       plots : Plot on the screen. Scalar.
%               If plots = 1, generates a plot of the monitoring of residuals
%               against alpha.
%                 Example - 'plots',0
%                 Data Types - single | double
%
%  Output:
%
%  out :     A structure containing the following fields
%
%            out.Beta = matrix containing the mpdp estimator of regression
%                       coefficients for each value of alpha
%            out.Scale= vector containing the estimate of the scale
%                       (sigma) for each value of alpha.
%              out.RES= n x length(alpha) matrix containing the robust
%                       scaled residuals for each value of bdp
%        out.Outliers = Boolean matrix containing the list of
%                       the units declared as outliers for each value of
%                       alpha using confidence level specified in input
%                       scalar conflev
%         out.conflev = confidence level which is used to declare outliers.
%         out.alpha   = vector which contains the values of alpha which
%                       have been used. To each value of alpha corresponds
%                       a value of bdp (see out.bdp).
%         out.bdp     = vector which contains the values of bdp which
%                       have been used. To each value of bdp corresponds
%                       a value of alpha (see out.alpha).
%            out.y    = response vector y. The field is present if option
%                       yxsave is set to 1.
%            out.X    = data matrix X. The field is present if option
%                       yxsave is set to 1.
%           out.class = 'MDPDReda'
%             out.Fval = Value of the objective function and reason fminunc
%                       or fminsearch stopped and v. Matrix.
%                       length(alpha)-by-3 matrix.
%                       The first column contains the values of alpha which
%                       have been considered.
%                       The second column contains the values of the
%                       objective function at the solution.
%                       The third column contins the details about
%                       convergence. A value greater then 0 denotes normal
%                       convergence. See help of functions fminunc.m or
%                       fminsearch.m for further details.
%
%
% More About:
%
%
% We assume that the random variables $Y|x$ are distributed as normal $N(
% \eta(x,\beta), \sigma_0)$ random variable with density function $\phi$.
% Note that if the model is linear $\eta(x,\beta)= x^T \beta$. The estimate
% of the vector $\theta_\alpha=(\beta_1, \ldots, \beta_p)^T$  (Minimum
% Density Power Divergence Estimate) is given by:
%
% \[
%  \mbox{argmin}_{\beta, \sigma} \left[ \frac{1}{\sigma^\alpha \sqrt{ (2
%  \pi)^\alpha(1+\alpha)}}-\frac{\alpha+1}{\alpha} \frac{1}{n} \sum_{i=1}^n
%  \phi^\alpha (y_i| \eta (x_i, \beta), \sigma)
%  \right]
%  \]
%
% As the tuning paramter $\alpha$ increases, the robustness of the Minimum
% Density Power Divergence Estimator (MDPDE) increases while its efficieny
% decreases (Basu et al. 1998). For $\alpha=0$ the MDPDE becomes the Maximum
% Likelihood Estimator, while for $\alpha=1$ the estimator minimizes the
% $L_2$ distance between the densities (Durio and Isaia, 2003),
%
%
%
% See also Sregeda, mdpdR
%
%  References:
%
%   Basu, A., Harris, I.R., Hjort, N.L. and Jones, M.C., (1998), Robust and
%   efficient estimation by minimizing a density power divergence,
%   "Biometrika", Vol. 85, pp. 549-559.
%   Durio,A.,Isaia,E.D.(2003), A parametric regression model by minimum L2
%   criterion: a study on hydrocarbon pollution of electrical transformers,
%   "Developments in Applied Statistics, Metodoloski Zvezki", Vol. 19, pp.
%   69-83.
%   Durio A., Isaia E.D. (2011), The Minimum Density Power Divergence
%   Approach in Building Robust Regression Models, "Informatica", Vol. 22,
%   pp. 43-56.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdpdReda')">Link to the help function</a>
%
%$LastChangedDate:: 2019-05-14 16:04:25 #$: Date of the last commit
%
%
% Examples:
%
%

%{
    % Call of mdpdReda with all default options.
    % Simulate a regression model.
    n=100;
    p=3;
    sig=0.01;
    eps=randn(n,1);
    X=randn(n,p);
    bet=3*ones(p,1);
    y=X*bet+sig*eps;
    % Contaminate the first 10 observations.
    y(1:10)=y(1:10)+0.05;
    [out] = mdpdReda(y, X,'plots',1);
%}

%{
    % Example of use of option alpha.
    n=100;
    p=3;
    sig=0.01;
    eps=randn(n,1);
    X=randn(n,p);
    bet=3*ones(p,1);
    y=X*bet+sig*eps;
    % Contaminate the first 10 observations.
    y(1:10)=y(1:10)+0.05;
    [out] = mdpdReda(y, X,'plots',1,'alphaORbdp','alpha','tuningpar',[1 0.8 0.2 0]);
%}

%{
    % mdpdReda applied to example 1 of Durio and Isaia (2011).
    % 600 points generated according to the model
    % Y=0.5*X1+0.5*X2+eps
    % and n2 = 120 points (outliers), drawn from the model
    % X1,X2~U(0,1) eps~N(0,0.1^2)
    n=600;
    p=2;
    sig=0.1;
    eps=randn(n,1);
    X=rand(n,p);
    bet=0.5*ones(p,1);
    y=X*bet+sig*eps;
    [out] = mdpdReda(y,X ,'plots',1);
%}

%{
    %% mdpdReda applied to Forbes data.
    load('forbes.txt');
    y=forbes(:,2);
    X=forbes(:,1);
    [outalpha0] = mdpdReda(y, X, 'plots',1);
%}

%{
    %% mdpdReda applied to multiple regression data.
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);
    [out] = mdpdReda(y, X, 'plots',1);
%}



%% Beginning of code

if nargin<2
    error('FSDA:mdpdReda:missingInputs','y or X is missing')
end

modelfun='';
theta0='';
intercept=1;
conflev=0.975;
plots=0;
alphaORbdp='alpha';
tuningpar=1:-0.01:0;

if nargin>2
    options=struct('alphaORbdp',alphaORbdp,'tuningpar',tuningpar,'intercept',intercept,'modelfun',modelfun,...
        'theta0',theta0,'conflev',conflev,...
        'plots',plots);
    
    UserOptions=varargin(1:2:length(varargin));
    
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mdpdReda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present in
    % structure options Remark: the nocheck option has already been dealt
    % by routine chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:mdpdReda:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    tuningpar=options.tuningpar;
    tuningpar=sort(tuningpar(:),'descend');
    alphaORbdp=options.alphaORbdp;
    
    if strcmp(alphaORbdp,'alpha')
        if min(tuningpar)<0
            error('FSDA:mdpdReda:WrongInputOpt','minimum value of alpha must be zero')
        end
        alphavec=tuningpar;
        % Store corresponding values of bdp,
        bdpvec= PDc(alphavec);
        
    elseif strcmp(alphaORbdp,'bdp')
        % In this case tuning paramter is breakdown point
        if min(tuningpar)<0 || max(tuningpar)>0.5
            error('FSDA:mdpdReda:WrongInputOpt','bdp must be in the interval [0 0.5]')
        end
        bdpvec = tuningpar;
        % convert the values of bdp into values of alpha
        alphavec=PDbdp(tuningpar);
    else
        error('FSDA:mdpdReda:WrongInputOpt','alphaORbdp must be ''bdp'' or ''alpha''')
    end
    
    intercept=options.intercept;
    modelfun=options.modelfun;
    theta0=options.theta0;
    conflev=options.conflev;
    plots=options.plots;
end

n=length(y);
conflev = (conflev+1)/2;
seq = 1:n;

% If the link is linear and there is the intercept add the column of ones
if isempty(modelfun) && intercept==1
    X=[ones(n,1) X];
end
p=size(X,2);

MaxIter=1000;
DisplayLevel='';
% nlinfitOptions=statset('Display',DisplayLevel,'MaxIter',MaxIter,'TolX',1e-7);
nlinfitOptions=statset('Display',DisplayLevel,...
    'MaxIter',MaxIter,'TolX',1e-8,'TolTypeX','rel',...
    'TolFun',1e-8,'TolTypeFun','rel');


% Use LMS solution as starting values of the parameters
if isempty(theta0)
    LS=false;
    if LS==true
        beta0  = X\y;
        yhat  = X*beta0;
        resMLE= y-yhat;
        sigma0 = sqrt(resMLE'*resMLE/(n-p));
    else
        outini=LXS(y,X,'nocheck',1,'msg',0,'nsamp',1000);
        beta0=outini.beta;
        sigma0=outini.scale;
    end
    
    theta0=[beta0;sigma0];
end


% Define matrices which will store relevant quantities
lalphavec=length(alphavec);
% Beta= matrix which will contain beta coefficients
Beta=zeros(p,lalphavec);
% Scale = vector which will contain the estimate of the scale
Scale=zeros(lalphavec,1);
Residuals=zeros(n,lalphavec);
Outliers=false(n,lalphavec);
Fval=[alphavec zeros(lalphavec,2)];

% Given that likfmin only accepts objective functions that depend only
% on a single variable (in this case betsigma)
% Vector of regression coefficients and scale

for jj=1:length(alphavec)
    alphajj=alphavec(jj);
    likfminOneParameter = @(betsigma)likfmin(betsigma, modelfun, y, X, alphajj);
    
    if isempty(modelfun) && alphajj == 0 % MLE of beta and sigma
        bhat=X\y;
        yhat=X*bhat;
        resMLE=y-yhat;
        scale=sqrt(resMLE'*resMLE/(n-p));
        residuals=resMLE/scale;
        exitflag=1;
    else
        [betaout,fval,exitflag]  = fminsearch(likfminOneParameter,theta0,nlinfitOptions);
        
        theta0=betaout;
        bhat=betaout(1:end-1);
        scale=betaout(end);
        if isempty(modelfun)
            yhat=X*bhat;
        else
            yhat=modelfun(bhat,X);
        end
        residuals=(y-yhat)/scale;
    end
    outliers=seq( abs(residuals)>norminv(conflev) );
    
    Residuals(:,jj)=residuals;
    Beta(:,jj)=bhat;
    Scale(jj)=scale;
    Outliers(outliers,jj)=true;
    Fval(jj,2:3)=[fval exitflag];
end

out.Beta = Beta;
out.Scale = Scale;
out.RES=Residuals;
out.Fval=Fval;

% Store in output structure the outliers
out.Outliers = Outliers;
% Store values of alphavec which have been used
% to each value of alphavec corresponds a value of bdpvec
out.alpha=alphavec;
out.bdp=bdpvec;

out.class='MDPDReda';

if intercept==1
    % Store X (without the column of ones if there is an intercept)
    out.X=X(:,2:end);
else
    out.X=X;
end
% Store response
out.y=y;
out.conflev=conflev;

% Plot residuals as function of the break down point
if plots==1
    laby='Scaled MPDP residuals';
    resfwdplot(out);
    ylabel(laby)
end


% likfmin = Objective function to call with fminunc or fminsearch
    function objyhat=likfmin(betsigma, modelfun, y, X, alpha)
        
        bet=betsigma(1:end-1);
        if isempty(modelfun)
            eta=X*bet;
        else
            eta=modelfun(bet,X);
        end
        
        res=y-eta;
        %   sigma=sqrt(res'*res/(n-p));
        sigma=betsigma(end);
        objyhat=1/(sigma^alpha *sqrt((2*pi)^alpha *(1+alpha))) ...
            -((alpha+1)/(alpha*n))*sum((normpdf(res,0,sigma)).^alpha);
    end

end

%FScategory:REG-Regression