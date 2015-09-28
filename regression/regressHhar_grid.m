function [out] = regressHhar_grid(y,X,Z,varargin)
%grid search to find minimum for ART model
%
%
%<a href="matlab: docsearch('regressHhar_grid')">Link to the help function</a>
%
% Required input arguments:
%
%
%    y:         Response variable. Vector. A vector with n elements that contains the response variable.
%               It can be either a row or column vector.
%    X :        Predictor variables in the regression equation. Matrix. Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables.
%               By default, there is a constant term in the model, unless
%               you explicitly remove it using option intercept, so do not
%               include a column of 1s in X.
%     Z :       Predictor variable in the skedastic equation. Vector. 
%               n x 1 vector of length containing the unique predictor in the skedastic equation
%
%               $\omega_i = 1 + exp(\gamma_0 + \gamma_1 Z(i,1) + ...+
%               \gamma_{r} Z(i,r))$. 
%
% Optional input arguments:
%
%  intercept :   Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%               Example - 'intercept',1
%               Data Types - double%  plots :    Plot on the screen. Scalar.
%               If equal to one a plot of profile loglikelihood
%               appears  on the screen 
%               else (default) no plot is shown.
%                 Example - 'plots',1
%                 Data Types - double
%  nocheck:   Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%  nocheck:   Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
%
% Copyright 2008-2015.
%
% Written by FSDA team
%
%
% Last modified 06-Feb-2015


nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

alpha=0.1:0.1:4;
options=struct('intercept',1,'alpha',alpha,'nocheck',0,'plots',0);

if nargin > 3
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end
alpha=options.alpha;

lalpha=length(alpha);

LogL=[alpha' zeros(lalpha,1)];

% Z = n-by-1 vector which contains the explanatory variables for
% heteroskedasticity

sigma2all=zeros(lalpha,1);

betaall=zeros(p,lalpha);

ij=1;
for i_alpha=1:length(alpha)
    omega2hat=real(Z.^alpha(i_alpha));
        sqweights = omega2hat.^(-0.5);
        
        % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
        Xw = bsxfun(@times, X, sqweights);
        yw = y .* sqweights;
        
        % Estimate of beta from (re)weighted regression (RWLS)
        newbeta = Xw\yw;
        
        newres2=(yw-Xw*newbeta).^2;
        % s2 = MLE of sigma2 on transformed data
        sigma2=sum(newres2)/n;
        sigma2all(ij)=sigma2;
        betaall(:,ij)=newbeta;
        % Construct the value of the negative loglikelihood
        
        loglik=sum(log(omega2hat))+n*log(sigma2);
        LogL(ij,2)= loglik;
    ij=ij+1;
end
[LogLmin,indmin]=min(LogL(:,2));

plots=options.plots;
if plots==1
 plot(LogL(:,1),LogL(:,2))
end
out=struct;

out.Beta=betaall(:,indmin);
out.Gamma=[log(sigma2all(indmin)); alpha(indmin)];
out.LogLmin=LogLmin;
out.LogL=LogL;
out.sigma2=sigma2all(indmin);
end