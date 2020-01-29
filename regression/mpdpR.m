function [out] = mpdpR(y, X, alpha, varargin)
%mpdpR allows to apply Minimum Density Power Divergence criterion to parametric regression problems.
%
%<a href="matlab: docsearchFS('mpdpR')">Link to the help function</a>
%
%
% Required input arguments:
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
%  alpha :      tuning parameter. Non negative scalar in the interval (0 1].
%               The robustness of the Minimum Density Power Divergence
%               estimator increases as the tuning parameter $\alpha$
%               increases, but at the same time its efficiency decreases
%               (Basu et al., 1998). For $\alpha=0$ the MDPDE becomes the
%               Maximum Likelihood estimator, while for $\alpha=1$ the
%               divergence yields the $L_2$ metric and the estimator
%               minimizes the $L_2$ distance between the densities, e.g.,
%               Scott (2001), Durio and Isaia (2003).
%                 Data Types - double
%
%
%  Optional input arguments:
%
%   modelfun   : non linear function to use. It is a function_handle or
%                an empty value (default). If modelfun is empty the link
%                between $X$ and $\beta$ is assumed to be linear, otherwise
%                it is necessary to specify a function (using @) that
%                accepts two arguments, a coefficient vector and the array
%                X and returns the vector of fitted values from the non
%                linear model y. For example, to specify the hougen
%                (Hougen-Watson) nonlinear regression function, use the
%                function handle @hougen.
%                 Example - 'modelfun', modelfun
%                 where modelfun = @(beta,X) X*beta(1).*exp(-beta(2)*X);
%                 Data Types - function_handle or empty value.
%  theta0       : Initial point. Vector or empty valu.
%                 Empty value or vector containing initial values for the
%                 coefficients (beta0 and sigma0) just in case modelfun is
%                 non empty. If modelfun is empty this argument is ignored.
%                 Example - 'beta0',[0.5 0.2 0.1]
%                 Data Types - double
%  intercept   :  Indicator for constant term. Scalar. If 1, and modelfun is
%                 empty (that is if the link between X and beta is linear)
%                 a model with constant term will be fitted (default), else
%                 no constant term will be included. This argument is
%                 ignored if modelfun is not empty.
%                 Example - 'intercept',1
%                 Data Types - double
%  dispresults :  Display results on the screen. Boolean. If dispresults is
%                 true (default) it is possible to see on the screen table
%                 Btable.
%                 Example - 'dispresults',false
%                 Data Types - Boolean
%      conflev :  Confidence level which is used to declare units as outliers.
%                 Scalar. Usually conflev=0.95, 0.975 0.99 (individual
%                 alpha) or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous
%                 alpha). Default value is 0.975
%                 Example - 'conflev',0.99
%                 Data Types - double
%       plots  :  Plot on the screen. Scalar.
%                 If plots = 1, generates a plot with the residuals
%                 against index number. The confidence level used to draw
%                 the confidence bands for the residuals is given by the
%                 input option conflev. If conflev is not specified a
%                 nominal 0.975 confidence interval is used.
%                 Example - 'plots',0
%                 Data Types - single | double
%       yxsave :  store X and y. Scalar. Scalar that is set to 1 to request
%                 that the response vector y and data matrix X are saved
%                 into the output structure out. Default is 0, i.e. no
%                 saving is done.
%                 Example - 'yxsave',1
%                 Data Types - double
%      MaxIter : maximum number of iterations allowed. Positive integer.
%                The default value is 1000-
%                Example - 'MaxIter',100
%                Data Types - double
%        TolX  : Tolerance for declaring convergence. Scalar.
%                The default value of TolX is 1e-7.
%                Example - 'TolX',1e-8
%                Data Types - double
%
%  Output:
%
%  out :     A structure containing the following fields
%
%            out.beta  = vector containing the MPDP estimator of regression
%                        coefficients.
%            out.scale = scalar containing the estimate of the scale
%                        (sigma).
%        out.residuals = n x 1 vector containing the estimates of the
%                        scaled residuals. The residuals are robust or not
%                        depending on the input value alpha.
%    out.fittedvalues = n x 1 vector containing the fitted values.
%        out.outliers = this output is present only if conflev has been
%                       specified. It is a vector containing the list of
%                       the units declared as outliers using confidence
%                       level specified in input scalar conflev.
%         out.conflev = confidence level which is used to declare outliers.
%                       Remark: scalar out.conflev will be used to draw the
%                       horizontal line (confidence band) in the plot.
%            out.y    = response vector Y. The field is present if option
%                       yxsave is set to 1.
%            out.X    = data matrix X. The field is present if option
%                       yxsave is set to 1.
%           out.class = 'Sreg'
%          out.Btable = table containing estimated beta coefficients,
%                       standard errors, t-stat and p-values
%                       The content of matrix B is as follows:
%                       1st col = beta coefficients and sigma (in the last
%                       element). This output is present just if input
%                       option dispresults is true.
%                       2nd col = standard errors;
%                       3rd col = t-statistics;
%                       4th col = p-values.
%       out.exitflag = Reason fminunc or fminsearch stopped. Integer.
%                       A value greater then 0 denotes normal convergence.
%                       See help of functions fminunc.m or fminsearch.m for
%                       further details.
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
% See also Sreg, mpdp
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
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mpdpR')">Link to the help function</a>
%
%$LastChangedDate:: 2019-05-14 16:04:25 #$: Date of the last commit
%
%
% Examples:
%
%

%{
    % Call of mpdpR with all default options.
    % Simulate a regression model.
    rng('default')
    rng(1000)
    n=1000;
    p=3;
    sig=0.01;
    eps=randn(n,1);
    X=randn(n,p);
    bet=3*ones(p,1);
    y=X*bet+sig*eps;
    [out] = mpdpR(y, X, 0.01);
%}

%{
    % Example of use of option yxsave.
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
    [out] = mpdpR(y, X, 1,'yxsave',1);
    resindexplot(out)
%}

%{
    %% Compare MLE estimator with MDPD estimator (uncontamindated data).
    % Scenario as in example 1 of Durio and Isaia (2011).
    % 600 points generated according to the model
    % Y=0.5*X1+0.5*X2+eps
    % and n2 = 120 points (outliers), drawn from the model
    % X1,X2~U(0,1) eps~N(0,0.1^2)
    close all;
    n   = 600;
    p   = 2;
    sig = 0.1;
    eps = randn(n,1);
    X   = rand(n,p);
    bet = 0.5*ones(p,1);
    y   = X*bet+sig*eps;
    [outalpha1] = mpdpR(y, X, 1);
    h1 = subplot(2,1,1);
    resindexplot(outalpha1,'h',h1);
    title('alpha=1','FontSize',15);

    h2 = subplot(2,1,2);
    [outalpha0] = mpdpR(y, X, 0);
    resindexplot(outalpha0,'h',h2);
    title('alpha=0','FontSize',15);

%}

%{
    % Compare MLE estimator with MDPD estimator (EX2).
    % Scenario as in example 2 of Durio and Isaia (2011).
    % 480 points generated according to the model
    % Y=0.5*X1+0.5*X2+eps
    % and n2 = 120 points (outliers), drawn from the model
    % Y =0.7X1 +0.7X2 + eps
    % X1,X2~U(0,1) eps~N(0,0.1^2)
    close all;
    sig  = 0.1;
    p    = 2;
    n1   = 480;
    eps1 = randn(n1,1);
    X1   = rand(n1,p);
    bet1 = 0.5*ones(p,1);
    y1   = X1*bet1+sig*eps1;
    n2   = 120;
    eps2 = randn(n2,1);
    X2   = rand(n2,p);
    bet2 = 0.7*ones(p,1);
    y2   = X2*bet2+sig*eps2;
    y    = [y1;y2];
    X    = [X1;X2];
    group=2*ones(n,1);
    group(1:n1)=1;
    yXplot(y,X,group)
    [out] = mpdpR(y, X, 1);
    h1 = subplot(2,1,1);
    resindexplot(out,'h',h1);
    title('alpha=1','FontSize',15);
    n=n1+n2;
    h2 = subplot(2,1,2);
    [outalpha0] = mpdpR(y, X, 0);
    resindexplot(outalpha0,'h',h2);
    title('alpha=0','FontSize',15)
    % Compare robust and MLE estimate
    disp(table(outalpha0.beta,out.beta,'VariableNames',{'MLE alpha=0' 'MD alpha=1'}))
%}

%{
    % Compare MLE estimator with MDPD estimator (EX3).
    % Scenario as in example 3 of Durio and Isaia (2011)
    % 480 points generated according to the model
    % Y=0.25*X1+0.25*X2+0.25*X3+0.25*X4+eps
    % and n2 = 120 points (outliers), drawn from the model
    % Y=0.35*X1+0.35*X2+0.35*X3+0.35*X4+eps
    % X1,X2,X3,X4~U(0,1) eps~N(0,0.1^2)
    close all;
    sig  = 0.1;
    p    = 4;
    n1   = 480;
    eps1 = randn(n1,1);
    X1   = rand(n1,p);
    bet1 = 0.25*ones(p,1);
    y1   = X1*bet1+sig*eps1;
    n2   = 120;
    eps2 = randn(n2,1);
    X2   = rand(n2,p);
    bet2 = 0.35*ones(p,1);
    y2   = X2*bet2+sig*eps2;
    y    = [y1;y2];
    X    = [X1;X2];
    group= 2*ones(n,1);
    group(1:n1)=1;
    yXplot(y,X,group)
    [out] = mpdpR(y, X, 1);
    h1 = subplot(2,1,1);
    resindexplot(out,'h',h1);
    title('alpha=1','FontSize',15);
    n = n1+n2;
    h2=subplot(2,1,2);
    % MLE estimate
    [outalpha0] = mpdpR(y, X, 0);
    resindexplot(outalpha0,'h',h2);
    title('alpha=0','FontSize',15);
    % Compare robust and MLE estimate
    disp(table(outalpha0.beta,out.beta,'VariableNames',{'MLE alpha=0' 'MD alpha=1'}))
%}

%{

    % Compare MLE estimator with MDPD estimator (EX4).
    % Scenario as in example 4 of Durio and Isaia (2011)
    % 180 points generated according to the model
    % Y  = 0.25*X1+eps
    % X1~U(0,0.5) eps~N(0,0.1^2)
    % and n2 = 20 points (outliers), drawn from the model
    % Y  = 0.25*X2+eps
    % X2~U(0.5,1) eps~N(0,0.1^2)
    % and m points (m=5, 10, 20, 30 40, 50)
    % Y  = 0.7*X3+eps3
    % X3~U(0.7,1) eps3~N(0,0.05^2)
    close all;
    sig  = 0.1;
    p    = 1;
    n1   = 180;
    eps1 = randn(n1,1);
    X1   = rand(n1,p)*0.5;
    y1   = X1+sig*eps1;
    n2   = 20;
    eps2 = randn(n2,1);
    X2   = rand(n2,p)*0.5+0.5;
    y2   = X2+sig*eps2;
    % Additional m points
    m    = 5;
    X3   = rand(m,p)*0.3+0.7;
    eps3 = randn(m,1);
    y3   = X3+0.05*eps3;
    y    = [y1;y2;y3];
    X    = [X1;X2;X3];
    group= 3*ones(n1+n2+m,1);
    group(1:n1)=1;
    group(n1+1:n1+n2)=2;

    yXplot(y,X,group)
    [out] = mpdpR(y, X, 1);
    h1=subplot(2,1,1);
    resindexplot(out,'h',h1);
    title('alpha=1','FontSize',15);
    n=n1+n2;
    h2=subplot(2,1,2);
    % MLE estimate
    [outalpha0] = mpdpR(y, X, 0);
    resindexplot(outalpha0,'h',h2);
    title('alpha=0','FontSize',15);
    % Compare robust and MLE estimate
    disp(table(outalpha0.beta,out.beta,'VariableNames',{'MLE alpha=0' 'MD alpha=1'}))
%}

%{
    %% MPDP applied to Forbes data.
    % Interactive_example
    clearvars;close all;
    % scatterplot of data: one point looks outlying
    load('forbes.txt');
    y=forbes(:,2);
    X=forbes(:,1);
    h1=subplot(2,1,1);
    [out] = mpdpR(y, X, 1);
    resindexplot(out,'h',h1);
    title('alpha=1','FontSize',15);

    h2=subplot(2,1,2);
    % MLE estimate
    [outalpha0] = mpdpR(y, X, 0);
    resindexplot(outalpha0,'h',h2);
    title('alpha=0','FontSize',15);
%}

%% Beginning of code

if nargin<3
    error('FSDA:mdpdR:missingInputs','y or alpha missing')
end

if ~isscalar(alpha)
    error('FSDA:mdpdR:WrongInputOpt','alpha should be a non negative scalar')
else
    if alpha<0
        error('FSDA:mdpdR:WrongInputOpt','alpha should be a non negative scalar')
    end
end


modelfun = '';
theta0   = '';
dispresults = false;
intercept   = 1;
conflev     = 0.975;
plots       = 0;
yxsave      = 0;
MaxIter     = 1000;
TolX        = 1e-7;

if nargin>3
    options=struct('intercept',intercept,'modelfun',modelfun,...
        'theta0',theta0,'dispresults',dispresults,'conflev',conflev,...
        'plots',plots,'yxsave',yxsave,'MaxIter',MaxIter,'TolX',TolX);
    
    UserOptions=varargin(1:2:length(varargin));   
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mpdpR:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present in
    % structure options Remark: the nocheck option has already been dealt
    % by routine chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:mpdpR:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    intercept   = options.intercept;
    modelfun    = options.modelfun;
    theta0      = options.theta0;
    dispresults = options.dispresults;
    conflev     = options.conflev;
    plots       = options.plots;
    yxsave      = options.yxsave;
end

n=length(y);

% If the link is linear and there is the intercept add the column of ones
if isempty(modelfun) && intercept==1
    X=[ones(n,1) X];
end
p=size(X,2);

if isempty(modelfun) && alpha == 0 
    % MLE of beta and sigma
    bhat  = X\y;
    yhat  = X*bhat;
    resMLE= y-yhat;
    scale = sqrt(resMLE'*resMLE/(n-p));
    residuals = resMLE/scale;
    exitflag  = 1;
else
    % Use linear squares as starting values of the parameters
    if isempty(theta0)
        beta0  = X\y;
        yhat0  = y-(X*beta0);
        sigma0 = sqrt(yhat0'*yhat0/(n-p));
        theta0 = [beta0;sigma0];
    else
        if isempty(modelfun)
            if length(theta0)~=p+1
                error('FSDA:mpdpR:WrongDim',...
                    ['Wrong dimension for theta0, it must be a vector with length = ', num2str(p+1)]);
            else
                % Just in case input theta0 is a rwo vector
                theta0=theta0(:);
            end
        end
    end
    
    DisplayLevel='';
    nlinfitOptions=statset('Display',DisplayLevel,'MaxIter',MaxIter,'TolX',TolX);
    
    % Given that likfmin only accepts objective functions that depend only
    % on a single variable (in this case betsigma)
    % Vector of regression coefficients and scale
    likfminOneParameter = @(betsigma)likfmin(betsigma, modelfun, y, X, alpha);
    if dispresults == true
        [betaout,~,exitflag,~,~,covB] = ...
            fminunc(likfminOneParameter,theta0,nlinfitOptions);
    else
        [betaout,~,exitflag]  = ...
            fminsearch(likfminOneParameter,theta0,nlinfitOptions);
    end
    
    bhat = betaout(1:end-1);
    scale = betaout(end);
    if isempty(modelfun)
        yhat=X*bhat;
    else
        yhat=modelfun(bhat,X);
    end
    residuals=(y-yhat)/scale;
end

% Compute scaled residuals
out=struct;
out.beta=bhat;
out.residuals=residuals;
out.fittedvalues=yhat;
out.exitflag=exitflag;
out.scale=scale;

if dispresults==true
    if isempty(modelfun) && alpha == 0
        se=sqrt(diag(scale^2*inv(X'*X))); %#ok<MINV>
    else
        se=sqrt(diag(inv(covB(1:end-1,1:end-1))));
    end
    % Show the estimated results
    tstat  = bhat./se;
    Btable = table(bhat,se,tstat);
    bnames = cellstr(num2str((1:length(bhat))','b%d'));
    
    Btable.Properties.RowNames=bnames;
    out.Btable=Btable;
    disp(Btable)
end

% Store in output structure the outliers found with confidence level conflev
out.conflev = conflev;

conflev = (conflev+1)/2;
seq = 1:n;
out.outliers = seq( abs(out.residuals)>norminv(conflev) );


if yxsave
    if options.intercept==1
        % Store X (without the column of ones if there is an intercept)
        out.X=X(:,2:end);
    else
        out.X=X;
    end
    % Store response
    out.y=y;
end

out.class = 'Sreg';

% Plot resindexplot with outliers highlighted
if plots==1
    laby='Scaled MPDP residuals';
    resindexplot(out.residuals,'conflev',out.conflev,'laby',laby,'numlab',out.outliers);
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