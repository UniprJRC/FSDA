function out=regressHart(y,X,sel,varargin)
%regressHart fits a multiple linear regression model using art heteroskedasticity
%
%<a href="matlab: docsearchFS('regresshart')">Link to the help function</a>
%
%  Required input arguments:
%
%    y:         A vector with n elements that contains the response variable.
%               It can be either a row or column vector.
%    X :        Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables.
%               By default, there is a constant term in the model, unless
%               you explicitly remove it using option intercept, so do not
%               include a column of 1s in X.
%     Z :       n x r matrix or vector of length r.
%               If Z is a n x r matrix it contains the r variables which
%               form the scedastic function as follows
%
%               \sigma^2_i = 1 + exp(\gamma_0 + \gamma_1 Z(i,1) + ...+ \gamma_{r} Z(i,r))
%
%               If Z is a vector of length r it contains the indexes of the
%               columns of matrix X which form the scedastic function as
%               follows
%
%               \sigma^2_i = 1 +  exp(\gamma_0 + \gamma_1 X(i,Z(1)) + ...+
%               \gamma_{r} X(i,Z(r)))
%
%               Therefore, if for example the explanatory variables
%               responsible for heteroscedisticity are columns 3 and 5
%               of matrix X, it is possible to use both the sintax
%                    regressHart(y,X,X(:,[3 5]))
%               or the sintax
%                    regressHart(y,X,[3 5])
%
%               Remark: Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%
%  Optional input arguments:
%
%   intercept : If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
% initialbeta : p x 1 vector containing initial estimate of beta. If
%               initialbeta is not supplied (default) standard least
%               squares is used to find initial estimate of beta
% initialgamma: vector of length (r+1) containing initial estimate of gamma.
%               If initialgamma is not supplied (default)  initial estimate
%               of gamma is nothing but the OLS estimate in a regression
%               where the response is given by squared residuals and the
%               regressors are specified in input object Z (this regression
%               also contains a constant term).
%     maxiter : scalar. Maximum number of iterations to find model paramters.
%               If not defined, maxiter is fixed to 200. Remark: in order
%               to obtain the FGLS estimator (two step estimator) it is
%               enough to put maxiter=1.
%     tol     : scalar. The tolerance for controlling convergence.
%               If not defined, tol is fixed to 1e-8. Convergence is
%               obtained if ||d_old-d_new||/||d_new||<1e-8 where d is the
%               vector of length p+r+1 which contains regression and scedastic
%               coefficients d=(\beta' \gamma')'; while d_old and d_new are
%               the values of d in iterations t and t+1 t=1,2, ..., maxiter
%    msgiter : scalar. If msgiter=1 it is possible to see the estimates of
%               the regression and scedastic parameters together with their
%               standard errors and the values of Wald, LM and
%               Likelihood ratio test, and the value of the maximized
%               loglikelihood. If msgiter>1 it is also possible to see
%               monitor the estimates of the coefficients in each step of
%               the iteration. If msgiter<1 nothing is displayed on the
%               screen
%
%  Output:
%
%  The output consists of a structure 'out' containing the following fields:
%
%           out.Beta  : p-by-3 matrix containing
%                       1st col = Estimates of regression coefficients
%                       2nd col = Standard errors of the estimates of regr coeff
%                       3rd col = t-tests of the estimates of regr coeff
%           out.Gamma : (r+1)-by-3 matrix containing
%                       1st col = Estimates of scedastic coefficients
%                       2nd col = Standard errors of the estimates of scedastic coeff
%                       3rd col = t tests of the estimates of scedastic coeff
%                       Remark: the first row of matrix out.Gamma is
%                       referred to the estimate of \sigma
%              out.WA : scalar. Wald test
%              out.LR : scalar. Likelihood ratio test
%              out.LM : scalar. Lagrange multiplier test
%            out.LogL : scalar. Complete maximized log likelihood
%
%
%   DETAILS. This routine implements art heteroscedasticity 
%
% See also regress, regstats
%
% References:
%
%   Atkinson A.C., Riani M. and Torti F. (2015), Robust methods for
%   heteroskedastic regression, submitted (ART)
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('regresshart')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    %The data in Appendix Table F6.1 were used in a study of efficiency in
    %production of airline services in Greene (2007a).
    % See p. 557 of Green (7th edition)

    % Monthly credit card expenditure for 100 individuals.

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
    out=regressHart(y,X,Loadfactor,'msgiter',1);
%}

%% Beginning of code
nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% maxiter = maximum number of iterations in the iterative procedure
maxiterdef=10000;
% Scalar defining convergence tolerance of iterative procedure
toldef=1e-08;
%toldef=1e-13;

options=struct('intercept',1,'maxiter',maxiterdef,...
    'initialbeta','','initialgamma','','tol',toldef,'msgiter',0,'type','art');

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
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
    
end
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
if size(sel,1)==n
    Z=[ones(n,1) sel];
else
    % Check if interecept was true
    intercept=options.intercept;
    if intercept==1
        Z=[ones(n,1) X(:,sel+1)];
    else
        Z=[ones(n,1) X(:,sel)];
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

% d is the column vector which contains the full set of
% parameters (beta and gamma)
dold=[oldbeta;oldgamma];
while cont==1 && iter<maxiter
    iter=iter+1;
    
    % 1. Estimate the disturbance variance
    % \sigma^2_i with \exp(\gamma_i' z_i)
    % sigma2hat1 = vector of length n which contains the n estimates of the
    % disturbance variance
    Zgamma=exp(Z*oldgamma);
    sigma2hati=(1+Zgamma);
    % sigma2hati(sigma2hati>1e+8)=1e+8;
    
    % 2. Compute \beta_{t+1} (new estimate of \beta using FGLS)
    sqweights = sigma2hati.^(-0.5);
    
    % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
    Xw = bsxfun(@times, X, sqweights);
    yw = y .* sqweights;
    
    % Estimate of beta from (re)weighted regression (RWLS)
    newbeta = Xw\yw;
    
    % 3. Update estimate of \gamma adding the vector of slopes in the least
    % squares regression of [\epsilon_i^2/exp(z_i^' \gamma)-1] on z_i
    newres2=(yw-Xw*newbeta).^2;
    newsigma2=sum(newres2)/n;
    
    Qweights=(Zgamma./(1+Zgamma));
    Zq = bsxfun(@times, Z, Qweights);
    %yq = (newres2./(newsigma2*(1+Zgamma))-1).* Qweights;
    % yq = (newres2./(newsigma2*(1+Zgamma))-1) ./ Qweights;
    newres2ori=(y-X*newbeta).^2;
    yq = (newres2ori./(newsigma2*(1+Zgamma))-1);
    
    newgamma=oldgamma+ delt*Zq\yq;
    % newgamma(abs(newgamma)>10)=10;
    th=8;
    newgamma(newgamma>th)=th;
    newgamma(newgamma<-th)=-th;
    
    dnew=[newbeta;newgamma];
    
    % Show estimate of beta in each step of the iteration
    if options.msgiter>1
        % disp(['Values of hat beta iteration' num2str(iter)])
        % disp(newbeta)
        disp(['Values of hat gamma iteration' num2str(iter)])
        disp(newgamma)
    end
    
    % lik=-0.5*(sum(log(sigma2hati))+sum((newres2./sigma2hati)));
    %disp(lik)
    % Check if convergence has been achieved
    if sum((dnew-dold).^2)/sum(dold.^2)>tol;
        cont=1;
        oldgamma=newgamma;
        dold=dnew;
    else
        cont=0;
    end
end

% Display a warning if convergence has not been achieved
if iter==maxiter && options.maxiter >1
    warning('FSDA:regressHart:NoConvergence','Maximum number of iteration has been reached without convergence')
end

% Store results

% Store beta coefficients, standard errors and corresponding tstats
BetaM=zeros(p,3);
BetaM(:,1)=newbeta;
% Compute standard errors of beta coefficients
BetaM(:,2)=sqrt(diag(inv(Xw'*Xw)));
% Compute t-stat of beta coefficients
BetaM(:,3)=BetaM(:,1)./BetaM(:,2);
out.BetaM=BetaM;
out.Beta=newbeta';
out.sigma2=newsigma2;

% Store parameters of scedastic function with associated standard errors
GammaM=zeros(length(newgamma),3);
GammaM(:,1)=newgamma;
% Find standard errors of elements of \gamma
% ZZ=Z'*Z;
% ZZ=0.5*ZZ(2:end,2:end);
covZ=2*inv(Zq'*Zq); %#ok<MINV>
GammaM(:,2)=diag(sqrt(covZ));
GammaM(:,3)=GammaM(:,1)./GammaM(:,2);


% Store inside out structure standard error of regression and heteroskedastic parameters
out.GammaM=GammaM;
% The two lines below are temporary just to have the connection with power
% model
out.alpha=out.GammaM(end,1);
out.Gamma=exp(out.GammaM(1,1));
end

