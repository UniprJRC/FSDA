function out=regressB(y, X, beta0, R, tau0, n0, varargin)
%regressB computes Bayesian estimates of regression parameters
%
%<a href="matlab: docsearch('regressb')">Link to the help function</a>
%
% Required input arguments:
%
%               SAMPLE INFORMATION
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n1, where n1 is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :          Predictor variables. Matrix. Matrix of explanatory variables (also called
%               'regressors') of dimension n1 x (p-1) where p denotes the
%               number of explanatory variables including the intercept.
%               Rows of X represent observations, and columns represent
%               variables. By default, there is a constant term in the
%               model, unless you explicitly remove it using input option
%               intercept, so do not include a column of 1s in X. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations.
%               Remark: note that here we use symbol n1 instead of
%               traditional symbol n because we want to better separate
%               sample information coming from n1 values to prior
%               information coming from n0 previous experiments.
%
%
%               PRIOR INFORMATION
%               $\beta$ is assumed to have a normal distribution with
%               mean $\beta_0$ and (conditional on $\tau_0$) covariance
%               $(1/\tau_0) (X_0'X_0)^{-1}$.
%               $\beta \sim N(    \beta_0, (1/\tau_0) (X_0'X_0)^{-1}    )$
%
%   beta0 :     Prior mean of $\beta$. p-times-1 vector.
%    R    :     Matrix associated with covariance matrix of $\beta$. p-times-p
%               positive definite matrix.
%               It can be interpreted as X0'X0 where X0 is a n0 x p
%               matrix coming from previous experiments (assuming that the
%               intercept is included in the model)
%
%               The prior distribution of $\tau_0$ is a gamma distribution with
%               parameters $a_0$ and $b_0$, that is
%
%                \[     p(\tau_0) \propto \tau^{a_0-1} \exp (-b_0 \tau)
%                       \qquad   E(\tau_0)= a_0/b_0               \]
%
%    tau0 :     Prior estimate of tau. Scalar. Prior estimate of $\tau=1/ \sigma^2 =a_0/b_0$.
%      n0 :     Number of previous experiments. Scalar. Sometimes it helps
%               to think of the prior information as coming from n0
%               previous experiments. Therefore we assume that matrix X0
%               (which defines R), was made up of n0 observations.
%
%  Optional input arguments:
%
%   intercept : Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%               Example - 'intercept',1
%               Data Types - double
%      bsb :   units forming subset. Vector.
%                m x 1 vector.
%               The default value of bsb is 1:n1, that is all n1 units are
%               used to compute beta1
%               Example - 'bsb',[3 5]
%               Data Types - double
%               REMARK: if bsb='' (empty value) just prior information is
%               used
%      c  :    it can be used to control the prior information
%               about beta. Scalar.
%               Scalar between 0 (excluded) and 1 (included).
%               The covariance matrix of the prior distribution
%               of $\beta$ is
%               \[
%               (1/tau0)* (c X0'X0)^{-1} = (1/tau0)* (c*R)^{-1}
%               \]
%               therefore multiplication of $R$ by $c$ (with $c<1$) increases the
%               variance of $\beta$.
%               The default value of $c$ is 1. Note also that if we adjust
%               the information for $\beta$ multiplying its covariance matrix
%               by $c$ we need to make a similar adjustment to the variance
%               of $\tau$. Simple mathematical calculations show that the new
%               variance of $\tau$ is equal to the old variance of $\tau$
%               divided by $c$. So, once again, given that $c < 1$, the new
%               variance of $\tau$ increases.
%               The value of $c$ will usually be the default, i,e., 1.
%               It was included to clarify the difference from Zellner's
%               procedure and has the effect of changing n0. For we can write
%               n0' = cn0 and use n0' in the calculations. However, it may
%               be useful in scaling X0 if the prior information is chosen
%               according to a design that is not of the appropriate size
%               to represent the amount of prior knowledge.
%               Example - 'c',1.2
%               Data Types - double
%   stats:   additional statistics. Scalar.
%               If stats =1 the following additional statistics are
%               computed:
%               1) Bayesian p-values.
%               2) highest posterior density intervals (HPDI) for each value
%               of input option conflev.
%               3) posterior odds for beta_j=0.
%               4) posterior model probability of the model which excludes
%               variable j
%               Example - 'stats',1
%               Data Types - double
%  conflev:     confidence levels to be used to
%               compute HPDI. Vector.
%               This input option is used just if input
%               stats=1. The default value of conflev is [0.95 0.99] that
%               is 95 per cent and 99 per cent HPDI confidence intervals
%               are computed
%               Example - 'conflev',[0.99 0.999]
%               Data Types - double
%
% Output:
%
%         out:   structure which contains the following fields
%
%   out.beta1=    p x 1 vector containing posterior mean (conditional on
%               $\tau_0$) of $\beta$ (regression coefficents)
%               $ \beta_1 = (c \times R + X'X)^{-1} (c \times R \times \beta_0 + X'y)$
%   out.a1    = shape parameter of the posterior gamma distribution of
%               $\tau$ ($\sigma^2$). Scalar.
%               $a1 = 0.5 (c*n_0 + n_1)$
%               The posterior distribution of $\tau$ is a gamma distribution
%               with parameters $a_1$ and $b_1$.
%               The posterior distribution of $\sigma^2$ is an inverse-gamma distribution
%               with parameters $a_1$ and $b_1$.
%   out.b1    = scale parameter of the posterior gamma distribution of
%               $\tau$ ($\sigma^2$). Scalar.
%               \[
%                   b1 = 0.5 * ( n_0 / \tau_0 + (y-X*\beta_1)'y +(\beta_0-\beta_1)'*c*R*\beta_0 )
%               \]
%               The posterior distribution of $\tau$ is a gamma distribution
%               with parameters $a_1$ and $b_1$.
%               The posterior distribution of $\sigma^2$ is an inverse-gamma distribution
%               with parameters $a_1$ and $b_1$.
%               Remark: note that if bsb is supplied X'X and X'y and n1 are
%               transformed respectively into Xm'Xm and Xm'ym and m where
%               Xm=X(bsb,:), ym=y(bsb) and m=length(bsb), therefore all the
%               previous quantities are estimated just using the units
%               forming subset
%   out.tau1  = posterior estimate of $\tau$. Scalar.
%               $\tau_1$ is obtained as $a_1/b_1$
%  out.covbeta1=    p x p matrix containing covariance matrix
%               (conditional on $\tau_1$) of $\beta_1$
%               $covbeta1 = cov(\beta_1)= (1/\tau_1) * (c*R + X'X)^{-1}$
%               where $\tau_1$ is obtained as $a_1/b_1$ (that is through the gamma
%               parameters of the posterior distribution of $\tau$)
% out.sigma21  = posterior estimate of $\sigma^2$. Scalar. $\sigma^2_1$ is
%               obtained as $b_1/(a_1-1)$ (that is through the inverse gamma
%               parameters of the posterior distribution of $\sigma^2$).
%               Remember that if $X\sim IG(a,b)$, then  $E(X)=b/(a-1)$.
%    out.res =   n1-times-2 matrix.
%               1st column = raw residuals
%               res(i,1) is computed as y(i) - X(i,:)*beta1 where beta1 is
%               computed using the units forming subset.
%               In the Bayesian approach they are the posterior means of
%               the $\epsilon_i$ and can be interpreted as point estimates of
%               the $\epsilon_i$.
%               2nd col = deletion residuals (just computed for the units
%               which do not form the subset).
%               res(i,2) with i \not \in  subset
%               is computed as
%               (y(i)-X(i,:)*beta1) * sqrt ((a1/b1)/(1+hii))
%               where
%               hii = X(i,:)* (c*R + Xm'*Xm)^{-1} * X(i,:)'.
%
%       The additional output which follows is produced just if input
%       scalar stats is equal 1
%
%  out.beta1HPD =   p-by-2*length(conflev) matrix HPDI of \hat \beta.
%               (HPDI stands for Highest Posterior Density Interval)
%               1st column =lower bound of HPDI associated with conflev(1).
%               2st column =upper bound of HPDI associated with conflev(1).
%               ...
%               2*length(conflev)-1 column = lower bound of HPDI associated
%               with conflev(end).
%               2*length(conflev) column (last column) = upper bound of
%               HPDI associated with conflev(end).
%  out.tau1HPD  =   1-by-2*length(conflev) matrix HPDI of $\hat \tau_1$.
%               1st element =lower bound of HPDI associated with conflev(1).
%               2st element =upper bound of HPDI associated with conflev(1).
%               ...
%               2*length(conflev)-1 element = lower bound of HPDI associated
%               with conflev(end).
%               2*length(conflev) element (last element) = upper bound of
%               HPDI associated with conflev(end).
%               Remark: confidence levels are based on the Gamma distribution
% out.sigma21HPD =   1-by-2*length(conflev) matrix HPDI of $\hat \sigma^_1$.
%               1st element =lower bound of HPDI associated with conflev(1).
%               2st element =upper bound of HPDI associated with conflev(1).
%               ...
%               2*length(conflev)-1 element = lower bound of HPDI associated
%               with conflev(end).
%               2*length(conflev) element (last element) = upper bound of
%               HPDI associated with conflev(end).
%               Remark: confidence levels are based on the inverse-gamma distribution
%    out.Bpval =   p-by-1 vector containing Bayesian p-values.
%               p-value = P(|t| > | \hat \beta se(beta) |)
%               = prob. of beta different from 0
%  out.postodds = p-by-1 vector which contains posterior odds for
%               betaj=0.
%               For example the posterior odd of betaj=0 is p(y| model
%               which contains all expl variables except the one associated
%               with betaj) divided by p(y| model which contains all expl
%               variables).
%               Warning: postodds can be computed just if n1>=p
% out.modelprob =   p-by-1 vector which contains  posterior model probability
%               of the model which excludes variable j.
%               For example if modelprob(j)= 0.28, that is if the
%               probability of the model which does not contain variable j
%               is equal to 0.28, it means that there is a 28 per cent
%               chance that beta_j=0 and a 72 per cent chance that it is
%               not.
%               Warning: modelprob can be computed just if n1>=p
%
%
%  More About:
%  The density of the gamma distribution (which we denote with $G(a,b)$) with parameters $a$ and $b$ is
%     \[
%     f_G(x,a,b)=  \frac{b^a}{\Gamma(a)}
%      x^{a-1} \exp(-b x)
%     \]
%     With this parametrization $G(a,b)$ has mean $a/b$ and variance $a/b^2$.
%
%     The density of the inverse gamma distribution defined over the support $x>0$ with shape parameter $a$ and scale parameter $b$, which we denote with $IG(a, b)$ is
%     \[
%     f_{IG}(x, a, b) = \frac{b^a}{\Gamma(a)}  (1/x)^{a +1} \exp (-b/x)
%     \]
%     The mean  (for $a>1$) is $b/(a-1)$ and the variance (for $a>2$) is $\frac{b^2}{(a-1)^2(a-2)}$.
%
%     With the above parameterizations if $X\sim G(a,b)$, then $Y=1/X$  has
%     an $IG(a, b)$.
%
%     The posterior distribution of
%     $\beta$ conditional on $\tau$ is
%     $N(\hat \beta_1, (1/\tau) (R+ X^TX)^{-1})$
%     where
%     \begin{eqnarray*}
%       \hat \beta_1 &=& (R+ X^TX)^{-1} (R \beta_0+X^T y) \\
%          &=&   (R+ X^TX)^{-1} (R \beta_0 + X^TX \hat \beta) \\
%          &=&   (I-A) \beta_0 + A \hat \beta,
%     \end{eqnarray*}
%     with $A= (R+ X^TX)^{-1} X^TX$. This last expression shows that the posterior estimate $\hat \beta_1$ is a matrix weighted average of the prior mean $\beta_0$ and the classical OLS estimate  $\hat \beta$, with weights $I-A$ and $A$. If prior information is strong elements of $R$ will be large, and $A$ will be small, so that the posterior mean gives most weight to the prior mean. In the classical approach these weights are fixed, while with the forward search technique as the subset size grows, the weight assigned to A becomes stronger and stronger, so we can dynamically see how the estimate changes as the effect of the prior becomes less important.
%
%     The posterior distribution of $\tau$ is $G(a_1, b_1)$ where
%     \begin{equation}\label{a1}
%     a_1=a+\frac{n}{2}= \frac{1}{2} (n_0 + n -p)
%     \end{equation}
%      and
%
%     \begin{equation}\label{b1}
%                        b_1 = 0.5 \left( \frac{n_0-p}{\tau_0} + (y-X \beta_1)'y +(\beta_0-\beta_1)'R\beta_0 \right)
%     \end{equation}
%     and, analogously, the posterior distribution of $\sigma^2$ is $IG(a_1, b_1)$.
%     The posterior estimates of $\tau$ and $\sigma^2$ are respectively:
%     \[
%     \hat \tau_1 = a_1/b_1, \qquad    \hat \sigma^2_1 = b_1 /(a_1-1)
%     \]
%
%     The posterior marginal distribution of $\beta$ is multivariate $T$ with parameters
%     \[
%     \hat \beta_1, (1/\tau_1)
%      \{R + X^T X \}^{-1}, n_0 +n - p
%     \]
%
% See also regress.m,
%
% References:
%
%   Chaloner and Brant (1988) Biometrika, Vol 75 pp. 651-659.
%
%   Koop G., Bayesian Econometrics (2003) Chapt. 3, WIley, NJ
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearch('regressB')">Link to the help function</a>
%
% Last modified 06-Feb-2015
%
% Examples:

%{
    %% regressB with all default options.
    % Common part to all examples: definition of the prior information.
    rng('default')
    rng(100) % set seed for replicability
    p=3;
    n0=30;
    X0=[ones(n0,1) randn(n0,p-1)];
    R=X0'*X0;
    beta0=zeros(p,1);
    tau0=1;

    % SAMPLE INFORMATION
    n1=100;
    X=randn(n1,p-1);
    y=randn(n1,1);

    % Run regressB using all n1 data and use the default value of c=1
    out=regressB(y, X, beta0, R, tau0, n0);
%}

%{
    % regressB with optional arguments.
    % Run regressB and compute new estimate of beta using just the first 20
    % observations and use a value of c equal to 1.2.
    bsb=1:20;
    c=1.2;
    out=regressB(y, X, beta0, R, tau0, n0,'bsb',bsb,'c',c);
%}

%{
    % Example of the use of input option stats.
    bsb=1:20;
    stats=true;
    out=regressB(y, X, beta0, R, tau0, n0,'bsb',bsb,'stats',stats);
%}

%{
    %% Example of Houses Price.
    % Compare the output with Table 3.3 of Koop (2004) p. 52.
    
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

    out=regressB(y, X, beta0, R, tau0, n0,'stats',1);
    disp(out)
%}

%% Beginning of code

nnargin=nargin;
vvarargin=varargin;
[y,X,n1,p] = chkinputRB(y,X,nnargin,vvarargin);

% default arguments values
bsbini=1:n1;
cini=1;
stats=0;

options=struct('intercept',1,'bsb',bsbini,'c',cini,'stats',stats,...
    'conflev',[0.95 0.99],'nocheck',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:regressB:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
    
end


c=options.c;
bsb=options.bsb;
stats=options.stats;

nbsb=numel(bsb);
Xbsb=X(bsb,:);
ybsb=y(bsb,:);

if nbsb>0 && nbsb<n1;
    aco=norminv(0.5*(1+nbsb/n1));
    corr=1-2*(n1./nbsb).*aco.*normpdf(aco);
    corr=sqrt(corr);
else
    corr=1;
end

% if stats=1 compute the ingredients to compute Bayesian confidence
% (credible) intervals
if stats==1
    cRXX1=(c*R+Xbsb'*Xbsb);
    cRXX1inv=inv(cRXX1);
end

Xbsb = Xbsb/corr;
ybsb=ybsb/corr;


% beta1=inv(c*R+Xbsb'*Xbsb)*(c*R*beta0+Xbsb'*ybsb);
XXbsb=Xbsb'*Xbsb;
cRXX=(c*R+XXbsb);
cRXXinv=inv(cRXX);
beta1=cRXXinv*(c*R*beta0+Xbsb'*ybsb); %#ok<MINV>

% The posterior distribution of \tau  is gamma distribution with parameters
% a1 and b1
a1 = 0.5 *(c*n0 + nbsb -p);

% Remark: notice that if bsb is empty 
% beta1=beta0, that is posterior beta is equal to prior beta
% and -beta1'*R*beta0 +beta0'*R*beta0=0 
b1 = 0.5 * ( c*(n0-p) / tau0 + ((ybsb-Xbsb*beta1)'*ybsb -beta1'*R*beta0) +beta0'*R*beta0 );

% Notation is as follows:
% tau1 = posterior estimate of $\tau$ (mean of the posterior distribution of
%        $\tau$). Post distrib of $\tau$ is  $\tau_1 \sim G(a_1, b_1)$
% sigma21 = posterior estimate of $\sigma^2$ (mean of the posterior distribution of
%        $\sigma^2$) Post distrib of $\tau$ is  $\tau_1 \sim IG(a_1, b_1)$

% tau1 = E[G(a_1, b_1)]
tau1=a1/b1;  % posterior estimate of $\tau$
% sigma21 = E[IG(a_1, b_1)]
sigma21=b1/(a1-1);   % posterior estimate of $\sigma^2$

% covbeta1 = (1/tau1) * inv(c*R + Xbsb'*Xbsb);
covbeta1 = (1/tau1)*cRXXinv; %#ok<MINV>



res=NaN(n1,2);

res(:,1)=y-X*beta1;


if nbsb<n1
    seq=1:n1;
    ncl=setdiff(seq,bsb);
    Xncl = X(ncl,:);
    
    % hi= element i,i of matrix H = Xncl * (c*R + Xm'*Xm)^{-1} * Xncl'
    % hi = X(i,:)* (c*R + Xbsb'*Xbsb)^{-1} * X(i,:)'
    % with i \not in bsb
    hi = sum((Xncl*cRXXinv).*Xncl,2);   %#ok<MINV>
    
    res(ncl,2) = sqrt(tau1)*(res(ncl,1)./sqrt(1+hi))/corr;
end



if stats==1 && nbsb>0
    %posterior probability that each element of beta is positive
    %as well as 95 and 99 HPDIs for each element of beta
    k=length(beta1);
    n0nbsb=n0+nbsb-p;
    
    % Bpval = Bayesian p-values
    % p-value = P(|t| > | \hat \beta se(beta) |)
    % = prob. of beta different from 0
    ci=sqrt((1/abs(tau1))*diag(cRXX1inv));
    tstat = -abs(beta1)./ci;
    Bpval = 2*tcdf(tstat,n0nbsb);
    
    
    % Compute highest posterior density interval for each value of
    % vector conflev
    conflev=options.conflev;
    conflev=1-(1-conflev)/2;
    % Tinvcdf = required quantiles of T distribution
    % consider just upper quantiles due to simmetry
    Tinvcdf=tinv(conflev,n0nbsb);
    % IGinvcdf = required quantiles of Inverse Gamma distribution
    IGinvcdf=inversegaminv([1-conflev conflev],a1,b1);
    % Ginvcdf = required quantiles of Gamma distribution
    Ginvcdf=gaminv([1-conflev conflev],a1,1/b1);
    
    beta1HPD=zeros(k,2*length(conflev));
    tau1HPD=zeros(1,2*length(conflev));
    sigma2HPD=tau1HPD;
    
    % The first two columns of matrix beta1HPD refer to conflev(1)
    % Columns three and four of matrix beta1HPD refer to conflev(2) ...
    lconflev=length(conflev);
    for j=1:lconflev
        beta1HPD(:,j*2-1:j*2)=[ beta1-Tinvcdf(j)*ci  beta1+Tinvcdf(j)*ci];
        tau1HPD(:,j*2-1:j*2)=[Ginvcdf(j) Ginvcdf(j+lconflev)];
        sigma2HPD(:,j*2-1:j*2)=[IGinvcdf(j) IGinvcdf(j+lconflev)];
    end
    
    
    % Computation of posterior odds for betaj=0
    
    %log of marginal likelihood for the model if prior is informative
    s02=1/tau0;
    Rinv=inv(R);
    % bols = ols beta coefficients
    bols=Xbsb\ybsb;
    % df = degrees of freedom
    df=nbsb-k;
    % resols = ols residuals of units forming subset
    resols=(ybsb-Xbsb*bols);
    % estimate of sigma^2 using units of the subsets
    s2 = resols'*resols/df;
    % deltab = difference between beta ols and beta of the prior
    deltab=(bols-beta0);
    % v1s12 = ingredient to compute the marginal likelihood (see lmarglik
    % below). v1s12 is nothing but equation (3.33), p. 41 of Koop (2004)
    v1s12 = n0*s02 + df*s2 + deltab'* ((Rinv + inv(XXbsb))\deltab);
    % logcj = log of (cj) see equation (3.35), p. 41 of Koop (2004)
    logcj=gammaln(.5*n0nbsb) + .5*n0*log(n0*s02)- gammaln(.5*n0) -.5*nbsb*log(pi);
    % lmarglik = log of the marginal likelihood for the full modell
    % lmarglik = log (y|full model which contains all expl variables)
    lmarglik=logcj + .5*log(det(cRXX1inv)/det(Rinv)) - .5*n0nbsb*log(v1s12);
    
    if n1>=p
        % postodds = vector which contains posterior odds for betaj=0
        % For example the posteriorodd of beta0=0 is p(y| model which contains
        % all expl variables except the one associated with beta0) divided by
        % p(y| model which contains all expl variables)
        postodds=zeros(k,1);
        seq1k=1:k;
        for j=seq1k
            selj=setdiff(seq1k,j);
            Rj=R(selj,selj);
            Xbsbj=Xbsb(:,selj);
            Rinvj=inv(Rj);
            cRXX1inv=inv(c*Rj+Xbsbj'*Xbsbj);
            beta0j=beta0(selj);
            bolsj=Xbsbj\ybsb;
            resolsj=ybsb-Xbsbj*bolsj;
            s2j = resolsj'*resolsj/(nbsb-(k-1));
            v1s12 = n0*s02 + (df+1)*s2j + (bolsj-beta0j)'* ((Rinvj + inv(Xbsbj'*Xbsbj))\(bolsj-beta0j));
            % v1s12 = v0*s02 + v*s2 + (bols-beta0)'*inv(capv0 + inv(X'*X))*(bols-beta0);
            logcj=gammaln(.5*n0nbsb) + .5*n0*log(n0*s02)- gammaln(.5*n0) -.5*nbsb*log(pi);
            % lmarglikj = log of the marginal likelihood for the model
            % which does not contain variable j
            % lmarglik = log (y|model which does not contain variable j)
            lmarglikj=logcj + .5*log(det(cRXX1inv)/det(Rinvj)) - .5*n0nbsb*log(v1s12);
            postodds(j)=exp(lmarglikj-lmarglik);
            % Remark: The noninformative choice, p(M_1)=p(M_2)=0.5 is commonly made
            % where  p(M_1) is the model without variable j and p(M_2) is the
            % full model with all the explanatory variables
        end
        
        % Remark: in order to obtain the posterior model probabilities for each
        % restricted model it is enough to do the following calculation
        % For example if modelprob(j)= 0.28, that is if the probability of the
        % model which does not contain variable j is equal to 0.28, it means
        % that there is a 28% chance that beta_j=0 and a 72% chance that it is
        % not.
        modelprob = postodds./(1+postodds);
    else
        warning('FSDA:regressB:Wrongn1','n1 must not be smaller than p to compute modelprob and postodds.');
        postodds=NaN(k,1);
        modelprob=postodds;
    end
    
    out=struct;
    out.beta1=beta1;
    out.a1=a1;
    out.b1=b1;
    out.tau1=tau1; % Posterior estimate of tau
    out.covbeta1=covbeta1; % cov matrix of beta1 |tau1
    out.sigma21=sigma21;   % posterior estimate of $\sigma^2$
    out.res=res;
    out.beta1HPD=beta1HPD;  % HPDI for beta
    out.tau1HPD=tau1HPD;   % HPDI for tau
    out.sigma21HPD=sigma2HPD;  % HPDI for Sigma2
    
    out.Bpval=Bpval;
    out.postodds=postodds;
    out.modelprob=modelprob;
else
    % ordinary output structure
    out=struct;
    out.beta1=beta1; % posterior estimates of beta
    out.a1=a1;       % shape parameter of Gamma (InverseGamma) distribution
    out.b1=b1;      % scale parameter of Gamma (InverseGamma) distribution
    out.tau1=tau1;   % posterior estimate of $\tau$
    out.covbeta1=covbeta1; % cov matrix of $\beta_1 |\tau_1$
    out.sigma21=sigma21;   % posterior estimate of $\sigma^2$
    out.res=res;
end
end