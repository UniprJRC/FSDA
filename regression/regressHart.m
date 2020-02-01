function [out]=regressHart(y,X,Z,varargin)
%regressHart fits a multiple linear regression model using ART heteroskedasticity
%
%<a href="matlab: docsearchFS('regressHart')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    y:         Response variable. Vector. A vector with n elements that contains the response variable.
%               It can be either a row or column vector.
%    X :        Predictor variables in the regression equation. Matrix.
%               Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables.
%               By default, there is a constant term in the model, unless
%               you explicitly remove it using option intercept, so do not
%               include a column of 1s in X.
%     Z :       Predictor variables in the skedastic equation. Matrix. n x r matrix or vector of length r.
%               If Z is a n x r matrix it contains the r variables which
%               form the scedastic function as follows: 
%
%               $\omega_i = 1 + exp(\gamma_0 + \gamma_1 Z(i,1) + ...+
%               \gamma_{r} Z(i,r))$. 
%
%               If Z is a vector of length r it contains the indexes of the
%               columns of matrix X which form the scedastic function as
%               follows: 
%
%               $\omega_i = 1 +  exp(\gamma_0 + \gamma_1 X(i,Z(1)) + ...+
%               \gamma_{r} X(i,Z(r)))$. 
%
%               Therefore, if for example the explanatory variables
%               responsible for heteroscedisticity are columns 3 and 5
%               of matrix X, it is possible to use both the sintax 
%                    regressHart(y,X,X(:,[3 5]))
%               or the sintax 
%                    regressHart(y,X,[3 5]). 
%
%               Remark: Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%
%  Optional input arguments:
%
%   intercept :  Indicator for constant term. true (default) | false. 
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
% initialbeta : initial estimate of beta. Vector.
%               p x 1 vector. If initialbeta is not supplied (default) standard least
%               squares is used to find initial estimate of beta
%               Example - 'initialbeta',[3.6 8.1]
%               Data Types - double
%
% initialgamma: initial estimate of gamma. Vector.
%               vector of length (r+1). If initialgamma is not supplied
%               (default)  initial estimate of gamma is nothing but the OLS
%               estimate in a regression where the response is given by
%               squared residuals and the regressors are specified in input
%               object Z (this regression also contains a constant term).
%               Example - 'initialgamma',[0.6 2.8]
%               Data Types - double
%
%     maxiter : Maximum number of iterations to find model paramters. Scalar.
%               If not defined, maxiter is fixed to 200. Remark: in order
%               to obtain the FGLS estimator (two step estimator) it is
%               enough to put maxiter=1.
%               Example - 'maxiter',8
%               Data Types - double
%
%     tol     : The tolerance for controlling convergence. Scalar.
%               If not defined, tol is fixed to 1e-8. Convergence is
%               obtained if $||d_{old}-d_{new}||/||d_{new}||<1e-8$ where d is the
%               vector of length p+r+1 which contains regression and scedastic
%               coefficients $d=(\beta' \gamma')'$; while $d_{old}$ and $d_{new}$ are
%               the values of d in iterations t and t+1 t=1,2,...,maxiter.
%               Example - 'tol',0.0001
%               Data Types - double
%
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
%
%
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
%          out.sigma2 = scalar. Estimate of $\sigma^2$ (sum of squares of
%                       residuals divided by n in the transformed scale)
%              out.WA = scalar. Wald test
%              out.LR = scalar. Likelihood ratio test
%              out.LM = scalar. Lagrange multiplier test
%            out.LogL = scalar. Complete maximized log likelihood
%
%
%   DETAILS. This routine implements art heteroscedasticity
%
%
%  More About:
%   The model is: 
%               $y=X*\beta+ \epsilon,
%               \epsilon ~ N(0, \Sigma) = N(0, \sigma^2*\Omega)$; 
%
%                   $\Omega=diag(\omega_1, ..., \omega_n)$; 
%
%                   $\omega_i=1+exp(z_i^T*\gamma)$; 
%
%                   $\Sigma=diag(\sigma_1^2, ...,
%                   \sigma_n^2)=diag(\sigma^2*\omega_1, ...,
%                   \sigma^2*\omega_n)$; 
%
%                   $var(\epsilon_i)=\sigma^2_i = \sigma^2 \omega_i \;\;\; i=1,
%                   ..., n$. 
%
%               $\beta$ = vector which contains regression parameters; 
%               $\gamma$= vector which contains skedastic parameters. 
%               REMARK 1: if $Z=log(X)$ then $1+exp(z_i^T*\gamma) =
%                         1+exp(\gamma_1)* \prod x_{ij}^{\gamma_j} \;\; j=1,
%                         ..., p-1$. 
%               REMARK2: if there is just one explanatory variable (say x)
%               which is responsible for heteroskedasticity and the model is
%               $\sigma_i=\sigma_2(1+ \theta*x_i^\alpha)$
%               then it is necessary to to supply Z as $Z=log(x)$. In this
%               case, given that the program automatically adds a column of
%               ones to Z: 
%                  $exp(Z_{1i}*\gamma_1 +Z_{2i}*\gamma_2)=
%                  exp(\gamma_1)*x_{1i}^{\gamma_2}$
%               therefore $exp(\gamma_1)$ is the estimate of $\theta$ while
%               $\gamma_2$ is the estimate of $\alpha$
%
% See also regress, regstats
%
%
% References:
%
% Greene, W.H. (1987), "Econometric Analysis", Prentice Hall. [5th edition,
% section 11.7.1 pp. 232-235, 7th edition, section  9.7.1 pp. 280-282]
% Atkinson, A.C., Riani, M. and Torti, F. (2016), Robust methods for
% heteroskedastic regression, "Computational Statistics and Data Analysis",
% Vol. 104, pp. 209-222, http://dx.doi.org/10.1016/j.csda.2016.07.002 [ART]
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('regressHart')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% regressHart with all default options.
    % The data in Appendix Table F6.1 were used in a study of efficiency in
    % production of airline services in Greene (2007a). See p. 557 of Green (7th edition).
    % Common part to all examples: load TableF61_Greene dataset.

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
    out=regressHart(y,X,Loadfactor);
%}

%{
    % regressHart with optional arguments.
    % Estimate a multiplicative heteroscedastic model and print the
    % estimates of regression and scedastic parameters together with LM, LR
    % and Wald test
    load('TableF61_Greene');
    Y=TableF61_Greene.data;

    Q=log(Y(:,4));
    Pfuel=log(Y(:,5));
    Loadfactor=Y(:,6);
    n=size(Y,1);
    X=[Q Q.^2 Pfuel];
    y=log(Y(:,3));
    out=regressHart(y,X,Loadfactor,'msgiter',1,'test',1);
%}



%% Beginning of code
<<<<<<< HEAD

=======
>>>>>>> master
nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% maxiter = maximum number of iterations in the iterative procedure
maxiterdef=10000;
% Scalar defining convergence tolerance of iterative procedure
toldef=1e-08;
%toldef=1e-13;
test=0;

options=struct('intercept',1,'maxiter',maxiterdef,'type','art',...
    'initialbeta','','initialgamma','','tol',toldef,'nocheck',0,'msgiter',0,'test',test);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:regressHart:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:regressHart:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
end
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
sigma2=r'*r/n;

% logL_R = -2*restricted loglikelihood  (without ln 2 \pi)
% -2* ( (-n/2)*[1+ln(r'r/n)   +ln 2 \pi]
logL_R= n*(1+log(sigma2));



%Initialization of gamma
% Z = n-by-r matrix which contains the explanatory variables for
% heteroskedasticity
if size(Z,1)==n
    Z=[ones(n,1) Z];
else
    % Check if interecept was true
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
    %    gamma0=Z\(r.^2/mean(r.^2)-1);
    gamma0=Z\(n*r.^2/sum(r.^2)-1);
    %    gamma0=Z\log(r.^2);
    %  gamma0=[0; 0];
    % gamma0=[log(50000); 1.9];
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

delt=1;

% Given that exp(z'\gamma) it is necessary in the estimatio procedure to
% put a threshold to z'\gamma. scalar |th| is the maximum value allowed during the
% estimation procedure for each element of gamma
th=8;

% d (dold and d new) are the column vectors which contains the full set of
% parameters (beta and gamma)
dold=[oldbeta;oldgamma];
while cont==1 && iter<maxiter
    iter=iter+1;
    
    % 1. Estimate the disturbance variance
    % \sigma^2_i with \exp(\gamma_i' z_i)
    % sigma2hat1 = vector of length n which contains the n estimates of the
    % disturbance variance
    Zoldgamma=Z*oldgamma;
    
    greatth=Zoldgamma>th;
    Zoldgamma(greatth)=th;
    
    smallth=Zoldgamma<-th;
    Zoldgamma(smallth)=-th;
    
    expZgamma=exp(Zoldgamma);
    
    omegahat=(1+expZgamma);
    
    % 2. Compute \beta_{t+1} (new estimate of \beta using FGLS)
    % Remark: as usual exp log in MATLAB is much more efficient than ^
    % sqweights = omegahat.^(-0.5);
    sqweights = exp(-0.5*log(omegahat));
    
    
    % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
    Xw = bsxfun(@times, X, sqweights);
    yw = y .* sqweights;
    
    % Estimate of beta from (re)weighted regression (RWLS)
    newbeta = Xw\yw;
    
    % 3. Update estimate of \gamma adding the vector of slopes in the least
    % squares regression of [\epsilon_i^2/exp(z_i^' \gamma)-1] on z_i
    newres2=(yw-Xw*newbeta).^2;
    newsigma2=sum(newres2)/n;
    
    Qweights=(expZgamma./(1+expZgamma));
    Zq = bsxfun(@times, Z, Qweights);
    %yq = (newres2./(newsigma2*(1+Zgamma))-1).* Qweights;
    % yq = (newres2./(newsigma2*(1+Zgamma))-1) ./ Qweights;
    newres2ori=(y-X*newbeta).^2;
    yq = (newres2ori./(newsigma2*(1+expZgamma))-1);
    
    newgamma=oldgamma+ delt*Zq\yq;
    
    
    dnew=[newbeta;newgamma];
    
    % Show estimate of beta in each step of the iteration
    if options.msgiter>1
        disp(['Values of hat beta iteration' num2str(iter)])
        disp(newbeta)
        disp(['Values of hat gamma iteration' num2str(iter)])
        disp(newgamma)
        if test==1
            disp('Tests')
%TODO
%             disp(['Likelihood ratio test    (LR)=' num2str(LR)])
%             disp(['Lagrange multiplier test (LM)=' num2str(LM)])
%             disp(['Wald test                (Wa)=' num2str(WA)])
%             disp(['Complete maximized log likelihood=' num2str(LogL)])
        end
    end
    
    % lik=-0.5*(sum(log(sigma2hati))+sum((newres2./sigma2hati)));
    %disp(lik)
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
    warning('FSDA:regressHart:NoConvergence','Maximum number of iterations has been reached without convergence')
end

if max(abs(newgamma))==th
    disp(['Max. value of \gamma during iterative procedure has been reached n.obs=' num2str(n)])
end

% Store results

% Store beta coefficients, standard errors and corresponding tstats
Beta=zeros(p,3);
Beta(:,1)=newbeta;
% Compute standard errors of beta coefficients
Beta(:,2)=sqrt(newsigma2*diag(inv(Xw'*Xw)));
% Compute t-stat of beta coefficients
Beta(:,3)=Beta(:,1)./Beta(:,2);
out.Beta=Beta;
out.BetaOLD=newbeta;
out.sigma2=newsigma2;

% Store parameters of scedastic function with associated standard errors
Gamma=zeros(length(newgamma),3);
Gamma(:,1)=newgamma;
% Find standard errors of elements of \gamma
% ZZ=Z'*Z;
% ZZ=0.5*ZZ(2:end,2:end);
covZ=2*inv(Zq'*Zq); %#ok<MINV>
Gamma(:,2)=diag(sqrt(covZ));
Gamma(:,3)=Gamma(:,1)./Gamma(:,2);


if test==1
    % Wald test
    % This test is computed extracting from the full parameter vector \gamma
    % and its estimated asymptotic covariance matrix, the subvector \hat alpha
    % and its asymptotic covariance matrix
    % In the case of art heteroscedasticity
    %               \omega_i = 1 +  exp(\gamma_0 + \gamma_1 X(i,Z(1)) + ...+
    %               \gamma_{r} X(i,Z(r)))
    % so alpha are all elements of vector gamma but
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
    logL_U= n*sum(log(newsigma2*(1+exp(Zgamma))))+sum(newres2ori./(newsigma2*(1+exp(Zgamma))));
    
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

% Store inside out structure standard error of regression and heteroskedastic parameters
out.Gamma=Gamma;

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

% The two lines below are temporary just to have the connection with power
% model
% out.alphaOLD=out.Gamma(end,1);
% out.GammaOLD=exp(out.Gamma(1,1));
end
%FScategory:REG-Hetero
