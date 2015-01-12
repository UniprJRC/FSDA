function out=regressB(y, X, beta0, R, tau0, n0, varargin)
%regressB computes Bayesian estimates of regression parameters
%
%<a href="matlab: docsearch('regressB')">Link to the help function</a>
%
% Required input arguments:
%
%               SAMPLE INFORMATION
%    y:         A vector with n1 elements that contains the response variable.
%               It can be either a row or a column vector.
%    X :        Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n1 x p-1). Rows of X represent observations, and
%               columns represent variables.
%               Remark: note that here we use symbol n1 instead of
%               traditional symbol n because we want to better separate
%               sample information coming from n1 values to prior
%               information coming from n0 previous experiments.
%
%
%               PRIOR INFORMATION
%               \beta is assumed to have a normal distribution with
%               mean \beta0 and (conditional on tau0) covariance
%               (1/tau0) (X0'X0)^{-1}
%               \beta~N(    beta0, (1/tau0) (X0'X0)^{-1}    )
%
%   beta0 :     p-times-1 vector containing prior mean of \beta
%    R    :     p-times-p positive definite matrix
%               which can be interepreted as X0'X0 where X0 is a n0 x p
%               matrix coming from previous experiments (assuming that the
%               intercept is included in the model
%
%               The prior distribution of tau0 is a gamma distribution with
%               paramters a and b, that is
%                     p(tau0) \propto \tau^{a0-1} \exp (-b0 \tau)
%                         E(tau0)= a0/b0
%
%
%    tau0 :     scalar. Prior estimate of tau=1/ \sigma^2 =a0/b0
%      n0 :     scalar. Sometimes it helps to think of the prior
%               information as coming from n0 previous experiments.
%               Therefore we assume that matrix X0 (which defines R), was
%               made up of n0 observations.
%
% Optional input arguments:
%
%
%  Optional input arguments:
%
%   intercept : If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%      bsb :    m x 1 vector containing the units forming subset.
%               The default value of bsb is 1:n1, that is all n1 units are
%               used to compute beta1
%               If bsb=NaN just prior information is used
%      c  :     scalar between 0 (excluded) and 1 (included) which can be
%               used to control the prior information
%               about beta. The covariance matrix of the prior distribution
%               of \beta is (1/tau0)* (c X0'X0)^{-1} = (1/tau0)* (c*R)^{-1}
%               therefore multiplication of R by c (with c<1) increases the
%               variance of \beta.
%               The default value of c is 1. Note also that if we adjust
%               the information for \beta multiplying its covariance matrix
%               by c we need to make a similar adjustment to the variance
%               of \tau. Simple mathematical calculations show that the new
%               variance of \tau is equal to the old variance of \tau
%               divided by c. So, once again, given that c < 1, the new
%               variance of \tau increases.
%               The value of c will usually be the default, i,e., 1.
%               It was included to clarify the difference from Zellner's
%               procedure and has the effect of changing n0. For we can write
%               n0' = cn0 and use n0' in the calculations. However, it may
%               be useful in scaling X0 if the prior information is chosen
%               according to a design that is not of the appropriate size
%               to represent the amount of prior knowledge.
%
%   stats:      Boolean, default=0, Additional statistics taken from the book:
%               'Bayesian Econometrics' by Gary Koop, Chapt. 3
%
%
% Output:
%
%  The output consists of a structure 'out' containing the following fields:
%     beta1:    p x 1 vector containing posterior mean (conditional on
%               tau0) of \beta (regression coefficents)
%               beta1 = (c*R + X'X)^{-1} (c*R*beta0 + X'y)
%  covbeta1:    p x p matrix containing posterior covariance matrix
%               (conditional on tau1) of \beta
%               covbeta1 = (1/tau1) * (c*R + X'X)^{-1}
%               where tau1 is defined as a1/b1 (that is through the gamma
%               parameters of the posterior distribution of \tau)
%
%               The posterior distribution of \tau is a gamma distribution
%               with parameters a1 and b1
%     a1    :   scalar parameter of the posterior gamma distribution of tau
%               a1 = 0.5 (c*n0 + n1)
%     b1    :   scalar parameter of the posterior gamma distribution of tau
%               b1 = 0.5 * ( n0 / tau0 + (y-X*beta1)'y +(beta0-beta1)'*c*R*beta0 )
%
%               Remark: note that if bsb is supplied X'X and X'y and n1 are
%               transformed respectively into Xm'Xm and Xm'ym and m where
%               Xm=X(bsb,:), ym=y(bsb) and m=length(bsb), therefore all the
%               previous quantities are estimated just using the units
%               forming subset
%       res :   n1-times-2 matrix
%               1st column = raw residuals
%               res(i,1) is computed as y(i) - X(i,:)*beta1 where beta1 is
%               computed using the units forming subset
%               In the Bayesian approach are they are the posterior means of the \epsilon_i
%               and can be interpreted as point estimates of the \epsilon_i
%               2nd col = deletion residuals (just computed for the units
%               which do not form the subset).
%               res(i,2) with i \not \in  subset
%               is computed as
%               (y(i)-X(i,:)*beta1) * sqrt ((a1/b1)/(1+hii))
%               where
%               hii = X(i,:)* (c*R + Xm'*Xm)^{-1} * X(i,:)'
%
%   probpos:    Posterior Odds probability for \beta_j=0
%   bhpdi95:    Bayesian Highest Posterior Density Intervals, 95% quantiles 
%   bhpdi99:    Bayesian Highest Posterior Density Intervals, 99% quantiles 

% See also FSRmmd.m,
%
% References:
%
%   Chaloner and Brant (1988) Biometrika, Vol 75 pp. 651-659.
%
%   Gary Koop, Bayesian Econometrics (2003) Chapt. 3
%
% Copyright 2008-2015.
% FSDA toolbox.
%
%<a href="matlab: docsearch('regressB')">Link to the help function</a>
%
% Last modified 01-Jan-2015
%
% Examples:

%{
%Common part to all examples:
% PRIOR INFORMATION
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
% Run regressB and compute new estimate of beta using just the first 20
% observations and use a value of c equal to 1.2
bsb=1:20;
c=1.2;
out=regressB(y, X, beta0, R, tau0, n0,'bsb',bsb,'c',c);
%}

%% Beginning of code
% nargin number of input arguments
if nargin<5
    error('Some input elements are missing');
end


nnargin=nargin;
vvarargin=varargin;
[y,X,n] = chkinputR(y,X,nnargin,vvarargin);

% default arguments values
bsbini=1:n;
cini=1;
stats=0;

options=struct('intercept',1,'bsb',bsbini,'c',cini,'stats',stats);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('Error:: number of supplied options is invalid. Probably values for some parameters are missing.');
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

nbsb=numel(bsb);
Xbsb=X(bsb,:);
ybsb=y(bsb,:);

if nbsb>0 && nbsb<n;
    aco=norminv(0.5*(1+nbsb/n));
    corr=1-2*(n./nbsb).*aco.*normpdf(aco);
    corr=sqrt(corr);
else
    corr=1;
end

cRX1=(c*R+Xbsb'*Xbsb);
cRX1inv=inv(cRX1);
Xbsb = Xbsb/corr;
ybsb=ybsb/corr;


% beta1=inv(c*R+Xbsb'*Xbsb)*(c*R*beta0+Xbsb'*ybsb);
cRX=(c*R+Xbsb'*Xbsb);
cRXinv=inv(cRX);
beta1=cRXinv*(c*R*beta0+Xbsb'*ybsb); %#ok<MINV>

% The posterior distribution of \tau  is gamma distribution with parameters
% a1 and b1
a1 = 0.5 *(c*n0 + nbsb);

b1 = 0.5 * ( n0 / tau0 + ((ybsb-Xbsb*beta1)'*ybsb -beta1'*R*beta0) +beta0'*R*beta0 );

% tau1 = posterior mean of \tau
tau1=a1/b1;

% covbeta1 = (1/tau1) * inv(c*R + Xbsb'*Xbsb);
covbeta1 = (1/tau1)*cRXinv; %#ok<MINV>

res=nan(n,2);

res(:,1)=y-X*beta1;


if nbsb<n
    seq=1:n;
    ncl=setdiff(seq,bsb);
    Xncl = X(ncl,:);
    
    % hi= element i,i of matrix H = Xncl * (c*R + Xm'*Xm)^{-1} * Xncl'
    % hi = X(i,:)* (c*R + Xbsb'*Xbsb)^{-1} * X(i,:)'
    % with i \not in bsb
    hi = sum((Xncl*cRXinv).*Xncl,2);   %#ok<MINV>
    
    res(ncl,2) = sqrt(tau1)*(res(ncl,1)./(1+hi))/corr;
    % res(ncl,3)= res(ncl,2)/corr;
end

if stats==1
    %posterior probability that each element of beta is positive
    %as well as 95 and 99 HPDIs for each element of beta
    k=length(beta1);
    probpos=zeros(k,1);
    bhpdi95=zeros(k,2);
    bhpdi99=zeros(k,2);
    s12=1/abs(tau1);
    v1=n0+n;
    %
    % get quantiles of t for calculating HPDIs
    invcdf95=tdis_inv(.975,v1);
    invcdf99=tdis_inv(.995,v1);
    %
    for i = 1:k
        tnorm = -beta1(i,1)/sqrt(s12*cRX1inv(i,i));
        probpos(i,1) = 1 - tdis_cdf(tnorm,v1);
        % column 1
        bhpdi95(i,1) = beta1(i,1)-invcdf95*sqrt(s12*cRX1inv(i,i));
        % column 2
        bhpdi95(i,2) = beta1(i,1)+invcdf95*sqrt(s12*cRX1inv(i,i));
        bhpdi99(i,1) = beta1(i,1)-invcdf99*sqrt(s12*cRX1inv(i,i));
        bhpdi99(i,2) = beta1(i,1)+invcdf99*sqrt(s12*cRX1inv(i,i));
    end
    
    % output structure with additional statistics
    % see Gary Koop Book
    out=struct;
    out.beta1=beta1;
    out.tau1=tau1;
    out.covbeta1=covbeta1;
    out.a1=a1;
    out.b1=b1;
    out.res=res;
    out.probpos=probpos;
    out.bhpdi95=bhpdi95;
    out.bhpdi99=bhpdi99;
else
    % ordinary output structure
    out=struct;
    out.beta1=beta1;
    out.tau1=tau1;
    out.covbeta1=covbeta1;
    out.a1=a1;
    out.b1=b1;
    out.res=res;
    
end
