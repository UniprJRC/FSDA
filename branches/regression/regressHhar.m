function [out]=regressHhar(y,X,Z,varargin)
%Fits a multiple linear regression model with Harvey heteroskedasticity
%
%<a href="matlab: docsearchFS('regressHhar')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    y:         Response variable. Vector. A vector with n elements that contains the response variable.
%               It can be either a row or column vector.
%    X :        Predictor variables in the regression equation. Matrix. Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables.
%               By default, there is a constant term in the model, unless
%               you explicitly remove it using option intercept, so do not
%               include a column of 1s in X.
%     Z :       Predictor variables in the skedastic equation. Matrix. n x r matrix or vector of length r.
%               If Z is a n x r matrix it contains the r variables which
%               form the scedastic function as follows: 
%               \[
%               \sigma^2_i = exp(\gamma_0 + \gamma_1 Z_{i1} + ...+
%               \gamma_{r} Z_{ir}). 
%               \]
%               If Z is a vector of length r it contains the indexes of the
%               columns of matrix X which form the scedastic function as
%               follows
%               \[
%               \sigma^2_i = exp(\gamma_0 + \gamma_1 X(i,Z(1)) + ...+
%               \gamma_{r} X(i,Z(r)))
%               \]
%               Therefore, if for example the explanatory variables
%               responsible for heteroscedisticity are columns 3 and 5
%               of matrix X, it is possible to use both the sintax
%                    regressH(y,X,X(:,[3 5]))
%               or the sintax
%                    regressH(y,X,[3 5]). 
%
%               Remark: Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%
%  Optional input arguments:
%
%   intercept : Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               else no constant term will be included.
%               Example - 'intercept',1
%               Data Types - double
% initialbeta : initial estimate of beta. Vector.
%               p x 1 vector. If initialbeta is not supplied (default) standard least
%               squares is used to find initial estimate of beta
%               Example - 'initialbeta',[3.6 8.1]
%               Data Types - double
% initialgamma: initial estimate of gamma. Vector.
%                vector of length (r+1). If initialgamma is not supplied (default)  initial estimate
%               of gamma is nothing but the OLS estimate in a regression
%               where the response is given by squared residuals and the
%               regressors are specified in input object Z (this regression
%               also contains a constant term).
%               Example - 'initialgamma',[0.6 2.8]
%               Data Types - double
%     maxiter : Maximum number of iterations to find model paramters. Scalar.
%               If not defined, maxiter is fixed to 200. Remark: in order
%               to obtain the FGLS estimator (two step estimator) it is
%               enough to put maxiter=1.
%               Example - 'maxiter',8
%               Data Types - double
%     tol     : The tolerance for controlling convergence. Scalar.
%               If not defined, tol is fixed to 1e-8. Convergence is
%               obtained if $||d_{old}-d_{new}||/||d_{new}||<1e-8$ where d is the
%               vector of length p+r+1 which contains regression and scedastic
%               coefficients $d=(\beta' \gamma')'$; while $d_{old}$ and $d_{new}$ are
%               the values of d in iterations t and t+1 t=1,2, ..., maxiter
%               Example - 'tol',0.0001
%               Data Types - double
%    msgiter : Level of output to display. Scalar.
%               If msgiter=1 it is possible to see the estimates of
%               the regression and scedastic parameters together with their
%               standard errors and the values of Wald, LM and
%               Likelihood ratio test, and the value of the maximized
%               loglikelihood. If msgiter>1 it is also possible to see
%               monitor the estimates of the coefficients in each step of
%               the iteration. If msgiter<1 nothing is displayed on the
%               screen
%               Example - 'msgiter',0
%               Data Types - double
%  nocheck:   Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
%  Output:
%
%  The output consists of a structure 'out' containing the following fields
%
%           out.Beta  = p-by-3 matrix containing: 
%                       1st col = Estimates of regression coefficients; 
%                       2nd col = Standard errors of the estimates of regr
%                       coeff; 
%                       3rd col = t-tests of the estimates of regr coeff.
%           out.Gamma = (r+1)-by-3 matrix containing: 
%                       1st col = Estimates of scedastic coefficients; 
%                       2nd col = Standard errors of the estimates of
%                       scedastic coeff; 
%                       3rd col = t tests of the estimates of scedastic
%                       coeff. 
%                       Remark: the first row of matrix out.Gamma is
%                       referred to the estimate of \( \sigma^2 \). In
%                       other words 
%                       \[ 
%                       \hat \sigma^2= \exp(\gamma_1) 
%                       \]
%              out.WA = scalar. Wald test
%              out.LR = scalar. Likelihood ratio test
%              out.LM = scalar. Lagrange multiplier test
%            out.LogL = scalar. Complete maximized log likelihood
%
%  More About:
%
%   This routine implements Harvey’s (1976) model of
%   multiplicative heteroscedasticity which is a very flexible, general
%   model that includes most of the useful formulations as special cases.
%   The general formulation is: 
%   
%       $\sigma^2_i =\sigma^2 \exp(z_i \alpha)$
%   
%   Let $z_i$ include a constant term so that \( z_i'=(1 \; q_i) \) where \( q_i \) is the
%   original set of variables which are supposed to explain
%   heteroscedasticity. This routine automatically adds a column of 1 to
%   input matrix Z (therefore Z does not have to include a constant term).
%   Now let 
%   \[  
%   \gamma'=[\log \sigma^2 \alpha'] = [ \gamma_0, \ldots, \gamma_r]. 
%   \] 
%   Then the model is simply
%   \[
%   \sigma^2_i = \exp(\gamma' z_i)
%   \]
%   Once the full parameter vector is estimated \( \exp( \gamma_0)\) provides the
%   estimator for \( \sigma^2 \). 
%
%  The model is:  
%               \[
%                 y=X \times\beta+ \epsilon,  \qquad 
%                 \epsilon \sim N(0, \;  \Sigma)
%               \]
%               \[
%                   \Sigma=diag(\sigma_1^2, ..., \sigma_n^2)    \qquad
%                   \sigma_i^2=\exp(z_i^T \times \gamma)        \qquad
%                   var(\epsilon_i)=\sigma_i^2                  \qquad
%                   i=1, ..., n    
%               \]                                                      
%               $\beta$ = p-by-1 vector which contains regression
%               parameters.                                                 
%               $\gamma$ = (r+1)-times-1 vector $\gamma_0, \ldots,
%               \gamma_r$ (or written in MATLAB language $(\gamma(1),
%               \ldots, \gamma(r+1))$) which contains skedastic parameters.
%               $X$ = n-by-p matrix containing explanatory variables in the
%               mean equation (including the constant term if present).
%               $Z$ = n-by-(r+1) matrix containing the explanatory
%               variables in the skedastic equation. $Z= (z_1^T, \ldots,
%               z_n^T)^T$ 
%               $z_i^T=(1, z_{i,1}, \ldots, z_{i,r})$ (or written in MATLAB
%               language  $z_i=(z(1), \ldots, z(r+1)$).
%               REMARK1: given that the first element of $z_i$ is equal to 1
%               $\sigma_i^2$ can be written as
%               \[
%               \sigma_i^2 = \sigma^2 \times \exp(z_i(2:r+1)*\gamma(2:r+1))
%                          = \exp(\gamma(1))*\exp(z_i(2:r+1)*\gamma(2:r+1))
%               \]
%               that is, once the full parameter vector $\gamma$ containing
%               the skedastic parameters is estimated $\exp( \gamma(1))$
%               provides the estimator for $\sigma^2$. 
%               REMARK2: if $Z=log(X)$ then
%               \[
%                            \sigma^2_i= \exp(z_i^T \times \gamma) =
%                           \prod_{j=1}^p x_{ij}^{\gamma_j}     \qquad  
%                            j=1, ..., p
%               \]
%               REMARK3: if there is just one explanatory variable (say $x
%               =(x_1 \ldots, x_n)$) which is responsible for
%               heteroskedasticity and the model is
%               \[
%               \sigma^2_i=( \sigma^2 \times x_i^\alpha)
%               \]
%               then it is necessary to supply $Z$ as $Z=log(x)$. In this
%               case, given that the program automatically adds a column of
%               ones to $Z$
%               \[
%                  \exp(Z(i,1) \times \gamma(1) +Z(i,2) \times \gamma(2))=
%                  \exp(\gamma(1))*x_i^{\gamma(2)}
%               \]
%               therefore the $\exp$ of the first element of vector
%               $\gamma$ (namely exp(gamma(1))) is the estimate of
%               $\sigma^2$ while the second element of vector $\gamma$
%               (namely gamma(2)) is the estimate of $\alpha$
%
%
% See also regressHart, regressH
%
% References:
%
%   Greene W.H. (1987), Econometric Analysis, Prentice Hall. [5th edition,
%   section 11.7.1 pp. 232-235, 7th edition, section  9.7.1 pp. 280-282].
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('regressHhar')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%
%{
    %% regressHhar with all default options.
    % Monthly credit card expenditure for 100 individuals.
    % Results in structure "OUT" coincides with "Maximum Likelihood
    % Estimates" of table 11.3, page 235, 5th edition of Greene (1987).
    % Results in structure "OLS" coincide with "Ordinary Leat Squares
    % Estimates" of table 11.3, page 235, 5th edition of Greene (1987).
 
    load('TableF91_Greene');
    data=TableF91_Greene.data;
    n=size(data,1);

    % Linear regression of monthly expenditure on a constant, age, income
    % its square and a dummy variable for home ownership using the 72
    % observations for which expenditure was nonzero produces the residuals
    % plotted below

    X=zeros(n,4);
    X(:,1)=data(:,3);%age
    X(:,2)=data(:,6);% Own rent (dummy variable)
    X(:,3)=data(:,4);% Income
    X(:,4)=(data(:,4)).^2; %Income  square
    y=data(:,5); % Monthly expenditure

    % Select the 72 observations for which expenditure was nonzero
    sel=y>0;
    X=X(sel,:);
    y=y(sel);
    whichstats={'r','tstat'};
    OLS=regstats(y,X,'linear',whichstats);
    r=OLS.r;

    disp('Ordinary Least Squares Estimates')
    LSest=[OLS.tstat.beta OLS.tstat.se OLS.tstat.t OLS.tstat.pval];
    disp(LSest)
    disp('Multiplicative Heteroskedasticity Model')
    % The variables which enter the scedastic function are Income and
    % Income square (that is columns 3 and 4 of matrix X)
    out=regressHhar(y,X,[3 4]);

    % Plot OLS residuals against Income (This is nothing but Figure 11.1 of
    % Green (5th edition) p. 216)
    plot(X(:,4),r,'o')
    xlabel('Income')
    ylabel('OLS residuals')
    grid on
%}

%{
    % regressHhar with optional arguments.
    % The data in Appendix Table F6.1 were used in a study of efficiency in
    % production of airline services in Greene (2007a).
    % See p. 557 of Green (7th edition).
    % Results in structure "out.Beta" coincide with those of
    % table 14.3 page 557, 7th edition of Greene (2007).
    % (line of the table which starts with MLE)

    load('TableF61_Greene');
    Y=TableF61_Greene.data;

    Q=log(Y(:,4));
    Pfuel=log(Y(:,5));
    Loadfactor=Y(:,6);
    n=size(Y,1);
    X=[Q Q.^2 Pfuel];
    y=log(Y(:,3));

    whichstats={'beta', 'r','tstat'};
    OLS=regstats(y,X,'linear',whichstats);

    disp('Ordinary Least Squares Estimates')
    LSest=[OLS.tstat.beta OLS.tstat.se OLS.tstat.t OLS.tstat.pval];
    disp(LSest)

    % Estimate a multiplicative heteroscedastic model and print the
    % estimates of regression and scedastic parameters together with LM, LR
    % and Wald test
    out=regressHhar(y,X,Loadfactor,'msgiter',1,'test',1);
%}

%{
    %% FGLS estimator.
    % Estimate a multiplicative heteroscedastic model using just one iteration
    % that is find FGLS estimator (two step estimator).
    % Data are monthly credit card expenditure for 100 individuals.
    % Results in structure "out" coincide with estimates of row
    % "$\sigma^2_i=\sigma^2 \exp(z' \alpha)$" in table 11.2, page 231, 5th edition of
    % Greene (1987).
    
    load('TableF91_Greene');
    data=TableF91_Greene.data;
  
    n=size(data,1);

    % Linear regression of monthly expenditure on a constant, age, income and
    % its square and a dummy variable for home ownership using the 72
    % observations for which expenditure was nonzero produces the residuals
    % plotted plotted below

    X=zeros(n,4);
    X(:,1)=data(:,3);%age
    X(:,2)=data(:,6);% Own rent (dummy variable)
    X(:,3)=data(:,4);% Income
    X(:,4)=(data(:,4)).^2; %Income  square
    y=data(:,5); % Monthly expenditure

    % Select the 72 observations for which expenditure was nonzero
    sel=y>0;
    X=X(sel,:);
    y=y(sel);
    out=regressHhar(y,X,[3 4],'msgiter',1,'maxiter',1);
%}

%% Beginning of code
nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% maxiter = maximum number of iterations in the iterative procedure
maxiterdef=200;
% Scalar defining convergence tolerance of iterative procedure
toldef=1e-08;
test=0;

options=struct('intercept',1,'maxiter',maxiterdef,...
    'initialbeta','','initialgamma','','tol',toldef,'msgiter',0,'test',test,'type','har','nocheck',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:regressHhar:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:regessHhar:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
end
% if test =1 procedure returns the results of Likelihood ratio, Lagrange
% multiplier and Wald test, together with the maximed log likelihood
test=options.test;

initialbeta=options.initialbeta;
if isempty(initialbeta)
    %Initialization of beta using standard least squares
    b0=X\y;
else
    b0=initialbeta;
end

r=y-X*b0;
oldbeta=b0;

% logL_R = -2*restricted loglikelihood  (without ln 2 \pi)
% -2* ( (-n/2)*[1+ln(r'r/n)   +ln 2 \pi]
logL_R= n*(1+log(r'*r/n));

%Initialization of gamma
% Z = n-by-r matrix which contains the explanatory variables for
% heteroskedasticity or vector containing the indexes of the columns of the
% explanatory variables responsible for heteroskedasticity
if size(Z,1)==n
    Z=[ones(n,1) Z];
else
    % Check if intercept is true
    intercept=options.intercept;
    if intercept==1
        Z=[ones(n,1) X(:,Z+1)];
    else
        Z=[ones(n,1) X(:,Z)];
    end
end
% Estimate of gamma is nothing but the OLS estimate in a regression where
% the response is given by squared residuals and the regressors are a set
% of variables Z supplied by the user (Z matrix may contain a subset of X
% or another set of explanatory variables)
% Note that the response is log(r.^2) because the scedastic function is
% \sigma_i^2 = exp( \gamma_0 + \gamma_1*Z(i,2) + ... + \gamma_{r}*Z(i,r))
initialgamma=options.initialgamma;
if isempty(initialgamma)
    gamma0=Z\log(r.^2);
else
    gamma0= initialgamma;
end

oldgamma=gamma0;

tol=options.tol;
%Iterative procedure
cont=1;
% iter = scalar which counts the number of iteration needed to achieve
% convergence
iter=0;

maxiter=options.maxiter;

% d is the column vector which contains the full set of
% parameters (beta and gamma)
dold=[oldbeta;oldgamma];

th=20;

while cont==1 && iter<maxiter
    iter=iter+1;
    
    % 1. Estimate the disturbance variance
    % \sigma^2_i with \exp(\gamma_i' z_i)
    % sigma2hat = vector of length n which contains the n estimates of the
    % disturbance variance
    Zoldgamma=Z*oldgamma;
    
    Zoldgamma(Zoldgamma>th)=th;
    Zoldgamma(Zoldgamma<-th)=-th;
    
    sigma2hat=exp(Zoldgamma);
    % sigma2hati(sigma2hati>1e+8)=1e+8;
    
    % 2. Compute \beta_{t+1} (new estimate of \beta using FGLS)
    sqweights = sigma2hat.^(-0.5);
    
    % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
    Xw = bsxfun(@times, X, sqweights);
    yw = y .* sqweights;
    
    % Estimate of beta from (re)weighted regression (RWLS)
    newbeta = Xw\yw;
    
    % 3. Update estimate of \gamma adding the vector of slopes in the least
    % squares regression of [\epsilon_i^2/exp(z_i^' \gamma)-1] on z_i
    newres2=(y-X*newbeta).^2;
    % newgamma=oldgamma+ Z\(newres2./exp(Zoldgamma)-1);
    newgamma=oldgamma+ Z\(newres2./sigma2hat-1);
    
    
    dnew=[newbeta;newgamma];
    
    % Show estimate of beta in each step of the iteration
    if options.msgiter>1
        disp(['Values of hat beta iteration' num2str(iter)])
        disp(newbeta)
        disp(['Values of hat gamma iteration' num2str(iter)])
        disp(newgamma)
    end
    
    
    % Check if convergence has been achieved
    if sum((dnew-dold).^2)/sum(dold.^2)>tol
        cont=1;
        oldgamma=newgamma;
        dold=dnew;
    else
        cont=0;
    end
end

% Display a warning if convergence has not been achieved
if iter==maxiter && options.maxiter >1
    warning('FSDA:regressHhar:NoConvergence','Maximum number of iteration has been reached without convergence')
end

% Store results

% Store beta coefficients, standard errors and corresponding tstats
Beta=zeros(p,3);
Beta(:,1)=newbeta;
% Compute standard errors of beta coefficients
Beta(:,2)=sqrt(diag(inv(Xw'*Xw)));
% Compute t-stat of beta coefficients
Beta(:,3)=Beta(:,1)./Beta(:,2);
out.Beta=Beta;

% Store parameters of scedastic function with associated standard errors
Gamma=zeros(length(newgamma),3);
Gamma(:,1)=newgamma;
% Find standard errors of elements of \gamma
% ZZ=Z'*Z;
% ZZ=0.5*ZZ(2:end,2:end);
covZ=2*inv(Z'*Z); %#ok<MINV>
Gamma(:,2)=diag(sqrt(covZ));
% % Estimate of sigma^2 Remark: \sigma^2 = exp( \hat gamma_1)
% Gamma(1,1)=exp(newgamma(1));
% % Find standard error of estimate of \sigma^2
% % The asymptotic variance of \sigma^2 is [exp(\gamma_1)]^2* Asy. Var
% % (\gamma_1)
% Gamma(1,2)=Gamma(1,1)*Gamma(1,2);
Gamma(:,3)=Gamma(:,1)./Gamma(:,2);

% Store inside out structure standard error of regression and heteroskedastic parameters
out.Gamma=Gamma;
out.alpha=out.Gamma(end,1);
out.sigma2=exp(out.Gamma(1,1));

if test==1
    % Wald test
    % This test is computed extracting from the full parameter vector \gamma
    % and its estimated asymptotic covariance matrix, the subvector \hat alpha
    % and its asymptotic covariance matrix
    % In the case of multiplicative heteroscedasticity
    % \gamma=[ln \sigma^2 alpha'] so alpha are all elements of vector gamma but
    % the first
    % Wald = \hat \alpha' \left{[0 I] [2 (Z'Z)]^{-1} [0 I] \right}^{-1} \hat \alpha
    % See p. 556 of Greene 7th edition. Note that on p. 556 of Greene there is
    % a missing inverse after the right curly bracket.
    alpha=Gamma(2:end,1);
    Varalpha=covZ(2:end,2:end);
    % Wald=alpha'*inv(Varalpha)*alpha;
    WA=alpha'*(Varalpha\alpha);
    
    % logL_U = -2*unrestricted loglikelihood
    Zgamma=Z*newgamma;
    logL_U= sum(Zgamma)+sum(newres2./exp(Zgamma));
    
    % LR = Likelihood ratio test
    LR=logL_R-logL_U;
    
    % Complete maximized log likelihood
    LogL= -(logL_U+n*log(2*pi))/2;
    
    % Lagrange multiplier test
    % Take residuals from OLS model and form reponse variable h (nx1) where
    % the ith element of vector h is given by
    % h_i= e_i^2/(e'e/n) -1  and z_i is the i-th row of matrix Z
    % i=1, ..., n
    h=r.^2/(sum(r.^2)/n)-1;
    % Regress h on Z and find the explained sum of squares
    % bh = vector of regression coefficients from regression of h on Z
    bh=Z\h;
    % Zbh = fitted values from the regression of h on Z
    Zbh=Z*bh;
    % LM is nothing but one-half times the explained sum of squares in the
    % linear regression of the variable h on Z
    LM=Zbh'*Zbh/2;
    
    % Below it is possible to find two alternative (inefficient) ways of
    % computing the LM test
    %
    %     % The row below is an inefficient way of computing the LM test
    %     % See equation 9.28 p. 276 of Greene 7th edition
    %     LM=h'*Z*inv(Z'*Z)*Z'*h/2;
    %
    %     % An additional alternative way of computing LM is as follows
    %     % vg is row vector of length columns of Z
    %     % vg = \sum v_i*z_i where v_i is a scalar equal to
    %     % h_i= e_i^2/(e'e/n) -1  and z_i is the i-th row of matrix Z
    %     vg = bsxfun(@times, Z(:,2:end), h);
    %     LM=sum(vg,1)*inv((n-1)*cov(Z(:,2:end)))*(sum(vg,1)')/2;
    
    % Store value of Likelihood ratio test
    out.LR = LR;
    % Store value of Lagrange multiplier test
    out.LM = LM;
    % Store value of Wald test
    out.WA = WA;
    % Store value of complete maximized log likelihood
    out.LogL=LogL;
    
end

msgiter=options.msgiter;
if msgiter ==1
    if maxiter>1
        disp('Regression parameters beta')
        disp('Coeff.   SE     t-stat')
        disp(Beta)
        disp('Scedastic parameters gamma')
        disp('Coeff.   SE ')
        disp(Gamma)
        if test==1
            disp('Tests')
            disp(['Likelihood ratio test    (LR)=' num2str(LR)])
            disp(['Lagrange multiplier test (LM)=' num2str(LM)])
            disp(['Wald test                (Wa)=' num2str(WA)])
            disp(['Complete maximized log likelihood=' num2str(LogL)])
        end
    else
        disp('Regression parameters beta')
        disp('Coeff.   SE     t-stat')
        disp(Beta)
        disp('Scedastic parameters gamma from first iteration')
        disp('Coeff.')
        disp(gamma0)
    end
end
end
%FScategory:REG-Hetero
