function out=regressH(y,X,Z,varargin)
%regressH fits a multiple linear regression model with heteroskedasticity
%
%<a href="matlab: docsearchFS('regressH')">Link to the help function</a>
%
%  Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
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
%     Z :       Predictor variables in the skedastic equation. Matrix. n x
%               r matrix or vector of length r. If Z is a n x r matrix it
%               contains the r variables which form the scedastic function.
%               If Z is a vector of length r it contains the indexes of
%               the columns of matrix X which form the scedastic function.
%               Therefore, if for example the explanatory variables
%               responsible for heteroscedisticity are columns 3 and 5
%               of matrix X, it is possible to use both the sintax
%                    regressH(y,X,X(:,[3 5]))
%               or the sintax
%                    regressH(y,X,[3 5])
%
%  Optional input arguments:
%
%   type:       Parametric function to be used in the skedastic equation.
%               String.
%               If type is 'arc' (default) than the skedastic function is
%               modelled as follows
%               \[
%               \sigma^2_i = \sigma^2 (1 + \exp(\gamma_0 + \gamma_1 Z(i,1) +
%                           \cdots + \gamma_{r} Z(i,r)))
%               \]
%               on the other hand, if type is 'har' then traditional
%               formulation due to Harvey is used as follows
%               \[
%               \sigma^2_i = \exp(\gamma_0 + \gamma_1 Z(i,1) + \cdots +
%                           \gamma_{r} Z(i,r)) =\sigma^2 (\exp(\gamma_1
%                           Z(i,1) + \cdots + \gamma_{r} Z(i,r))
%               \]
%               Remark. Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%               Example - 'type','har' 
%               Data Types - string
%
% intercept :   Indicator for constant term. true (default) | false. 
%               Indicator for the constant term (intercept) in the fit,
%               specified as the comma-separated pair consisting of
%               'Intercept' and either true to include or false to remove
%               the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
% initialbeta : initial estimate of beta. Vector.
%               p x 1 vector. If initialbeta is not supplied (default) standard least
%               squares is used to find initial estimate of beta
%               Example - 'initialbeta',[3 8] 
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
%
%     maxiter : Maximum number of iterations to find model paramters. Scalar. 
%               If not defined, maxiter is fixed to 200. Remark: in order
%               to obtain the FGLS estimator (two step estimator) it is
%               enough to put maxiter=1.
%               Example - 'maxiter',8 
%               Data Types - double
%
%     tol     : The tolerance for controlling convergence. Scalar
%               If not defined, tol is fixed to 1e-8. Convergence is
%               obtained if \( ||d_{old}-d_{new}||/||d_{new}||<1e-8 \)
%               where \( d \) is the vector of length p+r+1 which contains
%               regression and scedastic coefficients \( d=(\beta' \;
%               \gamma')' \); while \( d_{old} \) and \(d_{new} \) are the
%               values of d in iterations t and t+1 t=1,2, ..., maxiter
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
%  Output:
%
%         out:   structure which contains the following fields
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
%              out.WA = scalar. Wald test
%              out.LR = scalar. Likelihood ratio test
%              out.LM = scalar. Lagrange multiplier test
%            out.LogL = scalar. Complete maximized log likelihood using
%                       'har' or 'art' heteroskedasticity 
%
% See also regress, regstats
%
% References:
%
% Greene, W.H. (1987), "Econometric Analysis", Prentice Hall. [5th edition,
% section 11.7.1 pp. 232-235, 7th edition, section  9.7.1 pp. 280-282]
% Atkinson, A.C., Riani, M. and Torti, F. (2016), Robust methods for
% heteroskedastic regression, "Computational Statistics and Data Analysis",
% Vol. 104, pp. 209-222, http://dx.doi.org/10.1016/j.csda.2016.07.002 [ART]
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('regressH')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% regressH with all default options.
    % The data in Appendix Table F6.1 were used in a study of efficiency in
    % production of airline services in Greene (2007a).
    % See p. 557 of Green (7th edition). 
    
    load('TableF61_Greene');
    Y=TableF61_Greene{:,:};

    Q=log(Y(:,4));
    Pfuel=log(Y(:,5));
    Loadfactor=Y(:,6);
    n=size(Y,1);
    X=[Q Q.^2 Pfuel];
    y=log(Y(:,3));
    out=regressH(y,X,Loadfactor);
%}

%{
    % regressH with optional arguments.
    % Using the same data of the previous example and the traditional Harvey 
    % formulation for the skedastic function, we replicate in 
    % structure "out.Beta" the same results contained in table 9.2, 
    % page 282, 7th edition of Greene (2007) (lines "Iterated").
    load('TableF61_Greene');
    Y=TableF61_Greene{:,:};

    Q=log(Y(:,4));
    Pfuel=log(Y(:,5));
    Loadfactor=Y(:,6);
    n=size(Y,1);
    X=[Q Q.^2 Pfuel];
    y=log(Y(:,3));
    out=regressH(y,X,Loadfactor,'type','har');
%}

%{
    % Monthly credit card expenditure for 100 individuals.
    % Results in structure "OUT" coincides with "Maximum Likelihood
    % Estimates" of table 11.3, page 235, 5th edition of Greene (1987).
    % Results in structure "OLS" coincide with "Ordinary Leat Squares
    % Estimates" of table 11.3, page 235, 5th edition of Greene (1987).
 
    load('TableF91_Greene');
    data=TableF91_Greene{:,:};
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
    % The variables which enter the skedastic function are Income and
    % Income square (that is columns 3 and 4 of matrix X)
    out=regressH(y,X,[3 4],'msgiter',0,'type','har');

    % Plot OLS residuals againt Income (This is nothing but Figure 11.1 of
    % Green (5th edition) p. 216)
    plot(X(:,4),r,'o')
    xlabel('Income')
    ylabel('OLS residuals')
    grid on
%}

%{
    %% Comparing Harvey's and ART models.
    % Data are monthly credit card expenditure for 100 individuals.
    % Results in structure "out" coincides with estimates of row 
    % "\sigma^2*exp(z'*\alpha)" in table 11.2, page 231, 5th edition of
    % Greene (1987).
    
    load('TableF91_Greene');
    data=TableF61_Greene{:,:};
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

    % Compare output from Harvey's model with the one of ART
    outHAR=regressH(y,X,[3 4],'msgiter',1,'type','har');
    outART=regressH(y,X,[3 4],'msgiter',1,'type','art');

%}

%% Beginning of code

if nargin>3
    options=struct('type','art','intercept',1,'maxiter',100,...
        'initialbeta','','initialgamma','','tol',1e-7); %#ok<NASGU>

% check if input option type exists
    chklist=varargin(1:2:length(varargin));
    
    chktype = find(strcmpi('type',chklist)); 
    if ~isempty(chktype) && strcmp(varargin{2*chktype},'har') ==1
        out=regressHhar(y,X,Z,varargin{:});
    else
        out=regressHart(y,X,Z,varargin{:});
    end
else
     out=regressHart(y,X,Z);
end
%FScategory:REG-Hetero
