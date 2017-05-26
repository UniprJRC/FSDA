function [out] = regressHart_grid(y,X,Z,varargin)
%grid search to find minimum for ART model
%
%<a href="matlab: docsearch('regressHart_grid')">Link to the help function</a>
%
% Required input arguments:
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
%               $\omega_i = 1 + Z^\alpha \theta $. 
%
% Optional input arguments:
%
%  intercept :   Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               else no constant term will be included.
%               Example - 'intercept',1
%               Data Types - double
%  plots :    Plot on the screen. Scalar.
%               If equal to one a plot of profile loglikelihood
%               appears  on the screen 
%               else (default) no plot is shown.
%                 Example - 'plots',1
%                 Data Types - double
%   alpha  :   coefficient in the skedastic equation. Vector. Vector of
%               length r containing the values of coefficient alpha to
%               consider in the grid search. The default value for alpha is
%               alpha=0.1:0.1:4;
%               Example - 'alpha',0.1:0.1:5
%               Data Types - double
%   theta  :   coefficient in the skedastic equation. Vector. Vector of
%               length r containing the values of coefficient theta to
%               consider in the grid search. The default value for theta is
%               theta=[0.001 0.01 0.1 1 1.71 10:120 500 1000 5000 10000 50000];
%               Example - 'theta',0.1:0.1:5
%               Data Types - double
%  nocheck:   Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
%         out:   structure which contains the following fields
%
%out.Beta =  estimate of coefficients of regression equation. 
%            px1 vector. 
%out.Gamma= vector of length 2 containing estimates of skedastic
%           coefficients
%           1st elements = log (sigma2)
%           2nd element = estimate of alpha
%           The skedastic equation is 
%           omegahat=1+exp(out.Gamma(1,1))*exp(Z*out.Gamma(2,1))
%out.LogLmin= value of minimized  negative log lik. Scalar.
%out.sigma2= estimate of sigma2
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('regressHart_grid')">Link to the help function</a>
% Last modified 31-05-2016

% Examples:

%
%{
    % Random grid search (using simulated homoskedastic data).
    n=50;
    p=3;
    y=randn(n,1);
    X=randn(n,p);
    Z=exp(randn(n,1));
    HET_grid = regressHart_grid(y,X,exp(Z))

    %compare the estimates obtained with the grid with the estimates obtained
    %without the grid.
    HET_nogrid = regressHart(y,X,exp(Z));
    diff_Beta = HET_grid.Beta - HET_nogrid.Beta(:,1);
%}

%{
    %% regressHart_grid with all default options.
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
    out_grid=regressHart_grid(y,X,Loadfactor);

    %compare the estimates obtained with the grid with the estimates obtained
    %without the grid.
    out_nogrid = regressHart(y,X,Loadfactor);
    diff_Beta = out_grid.Beta - out_nogrid.Beta(:,1);
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
    out_grid=regressHart_grid(y,X,Loadfactor,'msgiter',1,'test',1);

    %compare the estimates obtained with the grid with the estimates obtained
    %without the grid.
    out_nogrid = regressHart(y,X,Loadfactor,'msgiter',1,'test',1);
    diff_Beta = out_grid.Beta - out_nogrid.Beta(:,1);
%}

%% Beginning of code

nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

alpha=0.1:0.1:4;
theta=[0.001 0.01 0.1 1 1.71 10:120 500 1000 5000 10000 50000];

options=struct('intercept',1,...
    'alpha',alpha,'theta',theta,'plots',0,'nocheck',0);


if nargin > 3
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

% alpha and theta if supplied by the user provide the grid of values of alpha
% and gamma for which the likelihood must be evaluated
alpha=options.alpha;
theta=options.theta;

% Z = n-by-1 vector which contains the explanatory variables for

lalpha=length(alpha);
ltheta=length(theta);
lthetalalpha=ltheta*lalpha;

sigma2all=zeros(1,lthetalalpha);
betaall=zeros(p,lthetalalpha);
    
% Initialization of Xw;
Xw=X;

ij=1;
omegahatall=zeros(n,lthetalalpha);

% logLtmp=zeros(ltheta*lalpha,3);

for i_alpha=1:lalpha
    Zialpha=real(Z.^alpha(i_alpha));
    omegahat=1+Zialpha*theta;
    omegahatall(:,(i_alpha-1)*ltheta+1:i_alpha*ltheta)=omegahat;
    % sqweights=omegahat.^(-0.5);
    sqweights=exp(-0.5*log(omegahat));
    
    for j_theta=1:ltheta
        
        
        % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
        if p==2
            Xw(:,1)=X(:,1).*sqweights(:,j_theta);
            Xw(:,2)=X(:,2).*sqweights(:,j_theta);
        elseif p==1
            Xw=X.*sqweights(:,j_theta);
        else
            Xw = bsxfun(@times, X, sqweights(:,j_theta));
        end
        yw = y .* sqweights(:,j_theta);
        
        % Estimate of beta from weighted regression (transformed space)
        beta = Xw\yw;
        
        res2=(yw-Xw*beta).^2;
        % sigma2 = MLE of sigma2 using transformed data
        sigma2=sum(res2)/n;
        sigma2all(ij)=sigma2;
        betaall(:,ij)=beta;
                
               %  logliktmp=sum(log(omegahat(:,j_theta)))+n*log(sigma2);
               %  logLtmp(ij,:)=[alpha(i_alpha) theta(j_theta) logliktmp];
        
        ij=ij+1;
    end
end
loglik=sum(log(omegahatall),1)+n*log(sigma2all);
[~,indmin]=min(loglik);

% [~,indmintmp]=min(logLtmp(:,3));


out=struct;
% Store in a column vector of estimate regression coefficients 
out.Beta=betaall(:,indmin);

% Store values of gamma and alpha which maximized the likelihood
aa=repmat(alpha,ltheta,1);
gg=repmat(theta',lalpha,1);
out.Gamma= [log(gg(indmin)); aa(indmin)];

% Store estimate of sigma2 (sum of squares of residuals/n in the
% transformed scale)
out.sigma2=sigma2all(indmin);

out.GammaOLD=gg(indmin);
out.alphaOLD=aa(indmin);

%Store value of maximized log likelihood
out.LogLmin=loglik(indmin);

plots=options.plots;

if plots==1
    aa=reshape(loglik,length(theta),length(alpha));
    mina=min(min(aa));
    maxa=max(max(aa));
    for jj=1:40
        subplot(5,8,jj)
        
        plot(aa(:,jj))
        xlim([1 ltheta])
        
        %   plot(gam,aa(:,jj))
        % xlim([-10000 gam(end)])
        ylim([mina maxa])
        title(['\alpha=' num2str(alpha(jj))])
    end
end

end
