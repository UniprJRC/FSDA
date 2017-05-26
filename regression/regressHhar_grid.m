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
%               else no constant term will be included.
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
%
% Output:
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
%            omegahat=exp(Z*HET.Gamma(2:end,1))
%out.LogLmin= value of minimized negative log lik. Scalar.
% out. LogL= Value of Log lik for each value of alpha.
%            1st col = values of alpha
%            2nd col = values of log lik
%out.sigma2= estimate of sigma2
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('regressHhar_grid')">Link to the help function</a>
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
    HETgrid = regressHhar_grid(y,X,exp(Z))

    %compare the estimates obtained with the grid with the estimates obtained
    %without the grid.
    HETnogrid = regressHhar(y,X,exp(Z))
    diff_Beta = HETgrid.Beta - HETnogrid.Beta(:,1);
%}

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
    logX = log(X(:,1)+0.0001)
    out_grid=regressHhar_grid(y,X,logX);

    % Plot OLS residuals against Income (This is nothing but Figure 11.1 of
    % Green (5th edition) p. 216)
    plot(X(:,4),r,'o')
    xlabel('Income')
    ylabel('OLS residuals')
    grid on

    %compare the estimates obtained with the grid with the estimates obtained
    %without the grid.
    out_nogrid = regressHhar(y,X,logX);
    diff_Beta = out_grid.Beta - out_nogrid.Beta(:,1);

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
    out_grid=regressHhar_grid(y,X,Loadfactor,'msgiter',1,'test',1);

    %compare the estimates obtained with the grid with the estimates obtained
    %without the grid.
    out_nogrid = regressHhar(y,X,Loadfactor,'msgiter',1,'test',1);
    diff_Beta = out_grid.Beta - out_nogrid.Beta(:,1);
%}

%{
    %% FGLS estimator.
    % Estimate a multiplicative heteroscedastic model using just one iteration
    % that is find FGLS estimator (two step estimator).
    % Data are monthly credit card expenditure for 100 individuals.
    % Results in structure "out" coincide with estimates of row
    % "\sigma^2*exp(z'*\alpha)" in table 11.2, page 231, 5th edition of
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
    out_grid=regressHhar_grid(y,X,3,'msgiter',1,'maxiter',1);

    %compare the estimates obtained with the grid with the estimates obtained
    %without the grid.
    out_nogrid = regressHhar(y,X,3,'msgiter',1,'maxiter',1);
    diff_Beta = out_grid.Beta - out_nogrid.Beta(:,1);
%}

%% Beginning of code
nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

alpha=0.1:0.1:4;
options=struct('intercept',1,'alpha',alpha,'nocheck',0,'plots',0);

if nargin > 3
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
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

