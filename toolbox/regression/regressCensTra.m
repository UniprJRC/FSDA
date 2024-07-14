function [out] = regressCensTra(y,X, varargin)
%regressCensTra computes signed sqrt LR test for lambda in the censored (Tobit) model
%
%
%<a href="matlab: docsearchFS('regressCensTra')">Link to the help function</a>
%
% This routines estimates the regression coefficients, sigma and the
% trasformation parameter lambda. The model is estimated by Maximum
% Likelihood (ML) assuming a Gaussian (normal) distribution of the error
% term in the transformed scale. The maximization of the likelihood
% function is done by function fmincon of the optimization toolbox.
% It also computes the signed sqrt likelihood ratio test of
% $H_0:\lambda=\lambda_0$.
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
%               m x 1 vector of integers or logical vector of length n.
%               The default value of bsb is 1:n, that is all units are
%               used to compute parameter estimates.
%               Example - 'bsb',[3 5 20:30]
%               Data Types - double or logical vector
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
%
% initialbeta : initial estimate of beta. Vector.
%               (p-1) x 1 vector. If initialbeta is not supplied (default) standard least
%               squares is used to find initial estimate of beta
%               Example - 'initialbeta',[3 8]
%               Data Types - double
%
%  initialla : initial estimate of lambda. scalar.
%               Initial value of the transformation parameter which has to
%               be used in the optimization procedure
%               Example - 'initialbeta',[3 8]
%               Data Types - double
%
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
%      la0:     Transformation parameter to test. Scalar.
%               Value of the transformation parameter which has to be
%               tested using the signed sqrt root lik ratio test. The
%               default value is la0 is 1. For example to test log
%               transformation use la0=0. The family of transformations
%               which is used is Yeo and Johnson.
%                 Example - 'la0',0.5
%                 Data Types - double
%
%    left :     left limit for the censored dependent variable. Scalar.
%               If set to -Inf, the dependent variable is assumed to be not
%               left-censored; default value of left is zero (classical
%               Tobit model).
%               Example - 'left',1
%               Data Types - double
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
%           out.Beta  = (p+2)-by-3 matrix containing:
%                       1st col = Estimates of regression coefficients,
%                       sigma and lambda;
%                       2nd col = Standard errors of the estimates
%                       3rd col = t-tests of the estimates
%          out.beta0  = vector containing the estimates of the regression
%                       coefficients under z(lambda_0)
%           out.jacla0 = value of the Jacobian which has been used in order
%                       to normalize the observations under H_0. The
%                       Jacobian just uses uncensored observations.
%            out.LogL = scalar. Value of the maximized log likelihood using
%                       lambda
%                       (ignoring constants).
%         out.Exflag  = Reason fminunc stopped in the maximization of the
%                       unconstrained likelihood. Integer.
%                      out.Exflag is equal to 1
%                       if the maximization procedure did not produce
%                       warnings or the warning was different from
%                       "ILL Conditioned Jacobian". For any other warning
%                       which is produced (for example,
%                       "Overparameterized", "IterationLimitExceeded",
%                       'MATLAB:rankDeficientMatrix") out.Exflag is equal
%                       to -1;
%         out.signLR = Signed sqrt of the lik ratio test of lambda=la0.
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
% Here $a$ is the lower limit and $b$ is the upper limit of the dependent
% variable trasformed. If $a=-\infty$ or $b=\infty$, the dependent variable is not
% left-censored or right-censored, respectively.
%
% Censored regression models (including the standard Tobit model) are usually estimated by the Maximum Likelihood (ML) method. Assuming that the disturbance term $\varepsilon$ follows a normal distribution with mean 0 and variance $\sigma^2$, the log-likelihood function is
% \begin{aligned}
% \log L=\sum_{i=1}^N  & \left[  I_i^a \log \Phi\left(\frac{a-z_i(\lambda)^{\prime} \beta}{\sigma}\right)+I_i^b \log \Phi\left(\frac{x_i^{\prime} \beta-b}{\sigma}\right)  \right. \\
% & \left.+\left(1-I_i^a-I_i^b\right)\left(\log \phi\left(\frac{y_i-z_i(\lambda)^{\prime} \beta}{\sigma}\right)-\log \sigma\right)\right],
% \end{aligned}
%
% where $\phi( \cdot)$ and $\Phi( \cdot )$ denote the probability density
% function and the cumulative distribution function, respectively, of
% the standard normal distribution, and $I_i^a$ and $I_i^b$ are indicator
% functions with. $z_i(\lambda)$ are the transformed observations using Yeo and
% Johnson transformation normalized using the Jacabian based on uncensored
% observations.
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
% See also regressCens, regress, regstats, regressB, regressH, regressts
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
%<a href="matlab: docsearchFS('regressCensTra')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Tobit regression using the affairs dataset.
    % In the example of Kleiber and Zeileis (2008, p. 142), the number of a
    % person's extramarital sexual inter-courses ("affairs") in the past year
    % is regressed on the person's age, number of years married, religiousness,
    % occupation, and won rating of the marriage. The dependent variable is
    % left-censored at zero and not right-censored.
    % We show below the estimates of lambda, beta, sigma and transformation
    % parameter lambda
    load affairs.mat
    X=affairs(:,["age" "yearsmarried" "religiousness" "occupation" "rating"]);
    y=affairs(:,"affairs");
    out=regressCensTra(y,X,"dispresults",true);
%}

%{
    % Example of right censoring.
    % When left=-Inf and right=0 indicates that there is no left-censoring but
    % there is a right censoring at 0. The same model as above but with the
    % negative number of affairs as the dependent variable can be estimated  by
    load affairs.mat
    X=affairs(:,["age" "yearsmarried" "religiousness" "occupation" "rating"]);
    y=affairs(:,"affairs").*(-1);
    out=regressCensTra(y,X,"dispresults",true,"left",-Inf,"right",0);
    % This estimation does not return anymore beta parameters that have the opposite sign of
    % the beta parameters in the original model, because of lambda
%}

%{
    % Example of right censoring with X a table with categorical variables.
    % When left=-Inf and right=0 indicates that there is no left-censoring but
    % there is a right censoring at 0. The same model as above but with the
    % negative number of affairs as the dependent variable can be estimated  by
    load affairs.mat
    X=affairs(:,["age" "yearsmarried" "religiousness" "occupation" "rating"]);
    y=affairs(:,"affairs").*(-1);
    out=regressCensTra(y,X,"dispresults",true,"left",-Inf,"right",0);
%}

%{
    % Another example with right censoring.
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
    out=regressCensTra(y,X,'right',800,'left',-Inf,'dispresults',true);
    disp('Value of the signed sqrt lik ratio test')
    disp(out.signLR)
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
la0=1;
initialla=[];

options=struct('nocheck',false,'dispresults',dispresults,...
    'bsb',bsb,'left',left,'right',right, 'la0',la0, ...
    'initialbeta',initialbeta,'initialla',initialla,'intercept',intercept,'optmin',optimset);

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
    la0=options.la0;
    initialla=options.initialla;
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
    error('FSDA:regrCensTra:WrgInpt','Observations with y greater than right truncation point')
end
if any(yin<left)
    error('FSDA:regrCensTra:WrgInpt','Observations with y smaller than left truncation point')
end

yeqleft = yin == left;
yeqright = yin== right;
obsBetween = ~yeqleft & ~yeqright;

% bsbT if a logical vector of length n which contains true in correspondence of the units 
% which are not censored
% bsbT=true(n,1);
% bsbT(yeqleft)=false;
% bsbT(yeqright)=false;


if sum(yeqleft) + sum(yeqright) == 0
    if length(y)>50 & (left~=-Inf && isfinite(right))
        warning('FSDA:regressCens:WrongInputOpt',"there are no censored observations")
    end
end

if sum(obsBetween) == 0
    warning('FSDA:regressCens:WrongInputOpt',"there are no uncensored observations")
end


if isempty(initialbeta)
    initialbeta=Xin\yin;
end


e=(yin-Xin*initialbeta);
sigmaini=e'*e/(n-p);

% yeqleft=yin==left;
% yeqright=yin==right;
% obsBetween=(yin>left & yin<right);

% lainit = column vector containing initial estimate of beta and sigma
if isempty(initialla)
    lainit=[initialbeta;sigmaini;la0];
else
    lainit=[initialbeta;sigmaini;initialla];
end

iter=0;
conv=false;

lb=-Inf(length(lainit),1);
ub=+Inf(length(lainit),1);
% Lower bound for sigma
lb(end-1)=1e-12;

lb(end,1)=-3;
ub(end,1)=3;

% optmin.TolX=1e-05;
while conv ==false && iter <=5

    [betaout,fval,Exflag,~,~,~,hessian]  = fmincon(@loglik,lainit,[],[],[],[],lb,ub,[],optmin);

    if Exflag >0
        conv=true;
    else
        lainit=lainit.*unifrnd(0.9,1.1,p+2,1);
        iter=iter+1;
    end
end


sebetaout=sqrt(diag(inv(hessian)));
tout=betaout./sebetaout;

% Likelihood maximization under H_0
lb(end,1)=la0;
ub(end,1)=la0;
lainit(end)=la0;
[beta0,fval0]=fmincon(@loglik,lainit,[],[],[],[],lb,ub,[],optmin);
lambdahat=betaout(end,1);

% Compute the signed sqrt likelihood ratio test
signLR=sign(lambdahat-la0)*sqrt(fval0-fval);

% jacla0 = Jacobian based on la0 using uncensored observations
nonnegs = y >= 0;
jacla0=(exp(mean(log(   (1 + abs(y(obsBetween))).^(2 * nonnegs(obsBetween) - 1)) )))^(1 - la0);


dfe=n-p-1;
pval=2*(tcdf(-abs(tout), dfe));
Beta=[betaout sebetaout tout pval];


out=struct;
out.Beta=Beta;
out.Exflag=Exflag;
out.LogL=-fval;
out.signLR=signLR;
out.jacla0=jacla0;
out.beta0=beta0;


if dispresults == true
    if tableX==true
        if intercept==true
            lab=["(Intercept)"; string(lab); "sigma"; "lambda"];
        else
            lab=[string(lab) "sigma" "lambda"];
        end
    else
        if intercept==true
            lab=["(Intercept)"; "x"+(1:(p-1))'; "sigma"; "lambda"];
        else
            lab=["x"+(1:p)'; "sigma" "lambda"];
        end
    end
    disp('Observations:')
    Total=n; LeftCensored= sum(yeqleft);
    Uncensored=sum(obsBetween);
    RightCensored=sum(yeqright);
    disp([table(Total) table(LeftCensored) table(Uncensored) table(RightCensored)])
    disp('Coefficients')
    Coeff=[table(betaout) table(sebetaout) table(tout) table(pval)];
    Coeff.Properties.RowNames=lab;
    disp(Coeff)
    disp(['Number of observations: ' num2str(n) , ', Error degrees of freedom:' num2str(dfe)]);
    disp(['Log-likelihood: ' num2str(-fval)]);

    R2=corr(y,X*betaout(1:end-2));
    disp(['R-squared: ' num2str(R2)])
    saveX=true;
end

if saveX== true
    out.X=X;
end

    function dZ=loglik(theta)
        % Extract from theta, beta sigma and lambda
        beta=theta(1:end-2);
        sigma=theta(end-1);
        lambda=theta(end);

        yhat=Xin*beta;
        yout=normYJ([yin;left;right],1,lambda,'Jacobian',true,'bsb',obsBetween);
        ytra=yout(1:end-2);
        % leftf and rightf are the truncation points in the transformed scale
        leftf=yout(end-1);
        rightf=yout(end);
        dZ=-sum(yeqleft.*log(normcdf((leftf-yhat)/sigma)+1e-12)+...
            yeqright.*log(normcdf((yhat-rightf)/sigma)+1e-12)+...
            +obsBetween.*(log( normpdf((ytra-yhat)/sigma)+1e-12 )-log(sigma)));
    end
end