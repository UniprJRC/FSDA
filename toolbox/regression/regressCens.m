function [out] = regressCens(y,X, varargin)
%regressCens computes estimates of regression parameters under the censored (Tobit) model
%
%
%<a href="matlab: docsearchFS('regressCens')">Link to the help function</a>
%
%
% The model is estimated by Maximum Likelihood (ML) assuming a Gaussian
% (normal) distribution of the error term. The maximization of the
% likelihood function is done by function fminunc of the optimization toolbox.
%
%
% Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Data Types - array or table
%
%  X :          Predictor variables in the regression equation. Matrix. Matrix of
%               explanatory variables (also called 'regressors') of
%               dimension n x (p-1) where p denotes the number of
%               explanatory variables including the intercept. Rows of X
%               represent observations, and columns represent variables. By
%               default, there is a constant term in the model, unless you
%               explicitly remove it using input option intercept, so do
%               not include a column of 1s in X. Missing values (NaN's) and
%               infinite values (Inf's) are allowed, since observations
%               (rows) with missing or infinite values will automatically
%               be excluded from the computations.
%               Data Types - array or table
%
% Optional input arguments:
%
%      bsb :   units forming subset. Vector.
%                m x 1 vector.
%               The default value of bsb is 1:length(y), that is all units are
%               used to compute parameter estimates.
%               Example - 'bsb',[3 5 20:30]
%               Data Types - double
%
%
%  dispresults : Display results of final fit. Boolean. If dispresults is
%               true,  labels of coefficients, estimated coefficients,
%               standard errors, tstat and p-values are shown on the
%               screen in a fully formatted way. The default value of
%               dispresults is false.
%               Example - 'dispresults',true
%               Data Types - logical
%
%    left :     left limit for the censored dependent variable. Scalar.
%               If set to -Inf, the dependent variable is assumed to be not
%               left-censored; default value of left is zero (classical
%               Tobit model).
%               Example - 'left',1
%               Data Types - double
%
% initialbeta : initial estimate of beta. Vector.
%               p x 1 vector. If initialbeta is not supplied (default) standard least
%               squares is used to find initial estimate of beta
%               Example - 'initialbeta',[3 8]
%               Data Types - double
%
% intercept :   Indicator for constant term. true (default) | false.
%               Indicator for the constant term (intercept) in the fit,
%               specified as the comma-separated pair consisting of
%               'Intercept' and either true to include or false to remove
%               the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%
%  nocheck:     Check input arguments. Boolean.
%               If nocheck is equal to true no check is performed on
%               supplied structure model
%               Example - 'nocheck',false
%               Data Types - logical
%
%
%    right :    right limit for the censored dependent variable. Scalar.
%               If set to Inf, the dependent variable is assumed to be not
%               right-censored; default value of right is Inf (classical
%               Tobit model).
%               Example - 'right',800
%               Data Types - double
%
%   optmin  :   It contains the options dealing with the
%               maximization algorithm. Structure.
%               Use optimset to set these options.
%               Notice that the maximization algorithm which is used is
%               fminunc if the optimization toolbox is present else is
%               fminsearch.
%                 Example -'optmin.Display','off'
%                 Data Types - double
%
%  Output:
%
%         out:   structure which contains the following fields
%
%           out.Beta  = p-by-3 matrix containing:
%                       1st col = Estimates of regression coefficients;
%                       2nd col = Standard errors of the estimates of
%                           regression coefficients;
%                       3rd col = t-tests of the estimates of  regression
%                           coefficients;
%            out.LogL = scalar. Value of the maximized log likelihood
%                   (ignoring constants)
%         out.Exflag  = Reason fminunc stopped. Integer.
%                      out.Exflag is equal to 1
%                       if the maximization procedure did not produce
%                       warnings or the warning was different from
%                       "ILL Conditiioned Jacobian". For any other warning
%                       which is produced (for example,
%                       "Overparameterized", "IterationLimitExceeded",
%                       'MATLAB:rankDeficientMatrix") out.Exflag is equal
%                       to -1;
%          out.X      = design matrix which has been used. This output is
%                       present just in X is a table containing categorical
%                       variables.
%
% More About:
%
%
% The issue is one where data is censored such that while we observe the
% value, it is not the true value, which would extend beyond the range of
% the observed data. This is very commonly seen in cases where the
% dependent variable has been given some arbitrary cutoff at the lower or
% upper end of the range, often resulting in floor or ceiling effects
% respectively. The conceptual idea is that we are interested in modeling
% the underlying latent variable that would not have such restriction if it
% was actually observed.
%
% In the standard Tobit model (Tobin 1958), we have a dependent variable $y$ that is left-censored at zero:
% \begin{eqnarray}
% y_i^* & = & x_i^{\prime} \beta+\varepsilon_i \\
% y_i   =  &  0            & \text { if } y_i^* \leq 0 \\
% y_i   =  &  y_i^*        &  \text { if } y_i^*>0
% \end{eqnarray}
%
% Here the subscript $i=1, \ldots, n$ indicates the observation,
% $y_i^*$ is an unobserved ("latent") variable,
% $x_i$ is a vector of explanatory variables,
% $\beta$ is a vector of unknown parameters, and
% $\varepsilon_i$ is a disturbance term.
%
% The censored regression model is a generalisation of the standard Tobit model.
% The dependent variable can be either left-censored, right-censored, or
% both left-censored and right-censored, where the lower and/or upper limit
% of the dependent variable can be any number:
% \begin{eqnarray}
%  y_i^*=x_i^{\prime} \beta+\varepsilon_i     \\
%  y_i & = &     a       \qquad  \text { if } y_i^* \leq a    \\
%      & = &   y_i^*     \qquad  \text { if } a<y_i^*<b    \\
%      & = &    b        \qquad  \text { if } y_i^* \geq b
%\end{eqnarray}
%
% Here $a$ is the lower limit and $b$ is the upper limit of the dependent variable. If $a=-\infty$ or $b=\infty$, the dependent variable is not left-censored or right-censored, respectively.
%
% Censored regression models (including the standard Tobit model) are usually estimated by the Maximum Likelihood (ML) method. Assuming that the disturbance term $\varepsilon$ follows a normal distribution with mean 0 and variance $\sigma^2$, the log-likelihood function is
% \begin{aligned}
% \log L=\sum_{i=1}^N  & \left[  I_i^a \log \Phi\left(\frac{a-x_i^{\prime} \beta}{\sigma}\right)+I_i^b \log \Phi\left(\frac{x_i^{\prime} \beta-b}{\sigma}\right)  \right. \\
% & \left.+\left(1-I_i^a-I_i^b\right)\left(\log \phi\left(\frac{y_i-x_i^{\prime} \beta}{\sigma}\right)-\log \sigma\right)\right],
% \end{aligned}
%
% where $\phi( \cdot)$ and $\Phi( \cdot )$ denote the probability density
% function and the cumulative distribution function, respectively, of
% the standard normal distribution, and $I_i^a$ and $I_i^b$ are indicator
% functions with
%
% \begin{eqnarray}
% I_i^a & = &  1  \text { if } y_i=a   \\
%       & = &  0   \text { if } y_i>a  \\
% I_i^b & = &  1  \text { if } y_i=b   \\
%       & = &  0  \text { if } y_i<b   \\
% \end{eqnarray}
%
% In this file the censored log likelihood above is maximized with respect
% to the parameter vector $(\beta', \sigma)$ using routine fminunc of the
% optimization toolbox
%
% See also regress, regstats, regressB, regressH, regressts
%
% References:
%
% Greene, W.H. (2008), "Econometric Analysis, Sixth Edition", Prentice Hall, pp. 871-875.
%
% Henningsen, A. (2012), Estimating Censored Regression Models in R using
% the censReg Package,
% [https://cran.r-project.org/web/packages/censReg/vignettes/censReg.pdf]
%
% Kleiber C., Zeileis A. (2008), "Applied Econometrics with R", Springer, New York.
%
% Tobin, J. (1958), Estimation of Relationships for Limited Dependent
% Variables, "Econometrica", 26, pp. 24-36.
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('regressCens')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Tobit regression using the affairs dataset.
    % In the example of Kleiber and Zeileis (2008, p. 142), the number of a
    % person's extramarital sexual inter-courses ("affairs") in the past year
    % is regressed on the person's age, number of years married, religiousness,
    % occupation, and won rating of the marriage. The dependent variable is
    % left-censored at zero and not right-censored. Hence this is a standard
    % Tobit model which can be estimated by the following lines
    load affairs.mat
    % Show the description of this table
    wraptextFS(affairs.Properties.Description,'width',80)
    % Define X and y 
    X=affairs(:,["age" "yearsmarried" "religiousness" "occupation" "rating"]);
    y=affairs(:,"affairs");
    out=regressCens(y,X,"dispresults",true);
%}

%{
    %% Example of right censoring.
    % When left=-Inf and right=0 indicates that there is no left-censoring but
    % there is a right censoring at 0. The same model as above but with the
    % negative number of affairs as the dependent variable can be estimated  by
    load affairs.mat
    X=affairs(:,["age" "yearsmarried" "religiousness" "occupation" "rating"]);
    y=affairs(:,"affairs").*(-1);
    out=regressCens(y,X,"dispresults",true,"left",-Inf,"right",0);
    % This estimation returns beta parameters that have the opposite sign of
    % the beta parameters in the original model, but the estimate of sigma does
    % not change.
%}

%{
    %% Example of right censoring with X a table with categorical variables.
    % When left=-Inf and right=0 indicates that there is no left-censoring but
    % there is a right censoring at 0. The same model as above but with the
    % negative number of affairs as the dependent variable can be estimated  by
    load affairs.mat
    X=affairs(:,["age" "yearsmarried" "religiousness" "occupation" "rating"]);
    y=affairs(:,"affairs").*(-1);
    out=regressCens(y,X,"dispresults",true,"left",-Inf,"right",0);
    % This estimation returns beta parameters that have the opposite sign of
    % the beta parameters in the original model, but the estimate of sigma does
    % not change.
%}

%{
    %% Another example with right censoring.
    % This example is taken from https://stats.oarc.ucla.edu/r/dae/tobit-models/
    % Consider the situation in which we have a measure of academic aptitude
    % (scaled 200-800) which we want to model using reading and math test
    % scores, as well as, the type of program the student is enrolled in
    % (academic, general, or vocational). The problem here is that students who
    % answer all questions on the academic aptitude test correctly receive a
    % score of 800, even though it is likely that these students are not
    % “truly” equal in aptitude. The same is true of students who answer all of
    % the questions incorrectly. All such students would have a score of 200,
    % although they may not all be of equal aptitude.
    XX=readtable("https://stats.idre.ucla.edu/stat/data/tobit.csv","ReadRowNames",true);
    XX.prog=categorical(XX.prog);
    % The dataset contains 200 observations. The academic aptitude variable is
    % "apt", the reading and math test scores are read and math respectively. The
    % variable prog is the type of program the student is in, it is a
    % categorical (nominal) variable that takes on three values, academic
    % general, and vocational (prog = 3). 
    % The scatterplot matrix is shown below
    spmplot(XX(:,[1 2 4]));
    % Now let’s look at the data descriptively. Note that in this dataset, the
    % lowest value of apt is 352. That is, no students received a score of 200
    % (the lowest score possible), meaning that even though censoring from
    % below was possible, it does not occur in the dataset.
    summary(XX)
    % Define y and X
    y=XX(:,"apt");
    X=XX(:,["read", "math" "prog"]);
    % Call regressCens
    out=regressCens(y,X,'right',800,'left',-Inf,'dispresults',true);
    % Show the plot of fitted vs residuals
    fitted=out.X*out.Beta(1:end-1,1);
    residuals=(y{:,1}-fitted)/out.Beta(end,1);
    figure
    scatter(fitted,residuals)
    xlabel('Fitted values (yhat)')
    ylabel('Residuals')
%}

%% Beginning of code
if istable(y)
    % namey=y.Properties.VariableNames;
    y=y{:,1};
end

if istable(X)

    tableX=true;

    catColumns = varfun(@iscategorical, X, 'OutputFormat', 'uniform');
    if any(catColumns)
        saveX=true;
        Xtmp=[];
        lab=[];

        for j=1:size(X,2)
            vnamej=X.Properties.VariableNames(j);
            if catColumns(j)==true
                Xcatj=X{:,j};
                Xj = dummyvar(Xcatj);
                dummynames=categories(Xcatj);
                colnamesj = strcat(vnamej,'_',dummynames);
                Xtmp=[Xtmp,Xj(:,2:end)];  %#ok<AGROW>
                lab=[lab; colnamesj(2:end)];  %#ok<AGROW>
            else
                Xtmp=[Xtmp,X{:,j}]; %#ok<AGROW>
                lab=[lab; vnamej];  %#ok<AGROW>
            end
        end
        X=Xtmp;
    else
        lab=X.Properties.VariableNames';
        X=X{:,:};
        saveX=false;
    end
else
    tableX=false;
    saveX=false;
end

% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = aux.chkinputR(y,X,nnargin,vvarargin);

%% User options
bsb=1:n;
dispresults=false;
left=0;
right=10^20;


initialbeta=[];
intercept=true;

options=struct('nocheck',false,'dispresults',dispresults,...
    'bsb',bsb,'left',left,'right',right, ...
    'initialbeta',initialbeta,'intercept',intercept,'optmin',optimset);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:regressCens:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    aux.chkoptions(options,UserOptions)
end

if nargin<2
    error('FSDA:regressCens:missingInputs','response y or X is missing');
end

if nargin >2

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    % Get options chosen by the user
    dispresults=options.dispresults;
    bsb=options.bsb;
    left=options.left;
    right=options.right;
    initialbeta=options.initialbeta;
    intercept=options.intercept;
    optmin=options.optmin;
else
    optmin=optimset;
end

optmin.Display='off';
optmin.MaxFunEvals=10000;

if isempty(bsb)
    Xin=X;
    yin=y;
else
    Xin=X(bsb,:);
    yin=y(bsb);
    n=length(yin);
end

if any(yin>right)
   error('FSDA:regressCens:WrgInpt','Observations with y greater than right truncation point')
end
if any(yin<left)
    error('FSDA:regressCens:WrgInpt','Observations with y smaller than left truncation point')
end

obsBelow = yin == left;
obsAbove = yin == right;
obsBetween = ~obsBelow & ~obsAbove;


if sum(obsBelow) + sum(obsAbove) == 0
    if length(y)>50 & (left~=-Inf && isfinite(right))
        warning('FSDA:regressCens:WrongInputOpt',"there are no censored observations")
    end
end

if sum(obsBetween) == 0
    warning('FSDA:regressCens:WrongInputOpt',"there are no uncensored observations")
end

if range(yin)>0

    if isempty(initialbeta)
        initialbeta=Xin\yin;
    end


    e=(yin-Xin*initialbeta);
    sigmaini=e'*e/(n-p);

    % lainit = column vector containing initial estimate of beta and sigma
    thetainit=[initialbeta;sigmaini];

    iter=0;
    conv=false;
    % optmin.TolX=1e-05;
    while conv ==false && iter <=5
        [betaout,fval,Exflag,~,~,hessian]  = fminunc(@loglik,thetainit,optmin);
        if Exflag >0
            conv=true;
        else
            thetainit=thetainit.*unifrnd(0.9,1.1,p+1,1);
            iter=iter+1;
        end
    end

    sebetaout=sqrt(diag(inv(hessian)));
    tout=betaout./sebetaout;

else
    betaini=Xin\yin;

    e=(yin-Xin*betaini);

    sigmaini=e'*e/(n-p);
    betaout=[betaini; sigmaini];

    sebetaout=[sqrt(sigmaini*diag(inv(X'*X))); Inf];
    tout=Inf(p+1,1);
    Exflag=0;
    fval=NaN;
end

dfe=n-p-1;
pval=2*(tcdf(-abs(tout), dfe));
Beta=[betaout sebetaout tout pval];


out=struct;
out.Beta=Beta;
out.Exflag=Exflag;
out.LogL=-fval;

if dispresults == true
    if tableX==true
        if intercept==true
            lab=["(Intercept)"; string(lab); "sigma"];
        else
            lab=[string(lab) "sigma"];
        end
    else
        if intercept==true
            lab=["(Intercept)"; "x"+(1:(p-1))'; "sigma"];
        else
            lab=["x"+(1:p)'; "sigma"];
        end
    end
    disp('Observations:')
    Total=n; LeftCensored= sum(obsBelow);
    Uncensored=sum(obsBetween);
    RightCensored=sum(obsAbove);
    disp([table(Total) table(LeftCensored) table(Uncensored) table(RightCensored)])
    disp('Coefficients')
    Coeff=[table(betaout) table(sebetaout) table(tout) table(pval)];
    Coeff.Properties.RowNames=lab;
    disp(Coeff)
    disp(['Number of observations: ' num2str(n) , ', Error degrees of freedom:' num2str(dfe)]);
    disp(['Log-likelihood: ' num2str(-fval)]);

    R2=corr(y,X*betaout(1:end-1));
    disp(['R-squared: ' num2str(R2)])
    saveX=true;
end

if saveX== true
    out.X=X;
end

    function dZ=loglik(theta)
        beta=theta(1:end-1);
        sigma=theta(end);
        yhat=Xin*beta;
        dZ=-sum(obsBelow.*log(normcdf((left-yhat)/sigma)+1e-12)+...
            obsAbove.*log(normcdf((yhat-right)/sigma)+1e-12)+...
            +obsBetween.*(log( normpdf((yin-yhat)/sigma)+1e-12 )-log(sigma)));
    end
end