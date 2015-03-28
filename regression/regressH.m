function out=regressH(y,X,Z,varargin)
%regressH fits a multiple linear regression model with heteroskedasticity
%
%<a href="matlab: docsearchFS('regressh')">Link to the help function</a>
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
%               form the scedastic function 
%               If Z is a vector of length r it contains the indexes of the
%               columns of matrix X which form the scedastic function 
%               Therefore, if for example the explanatory variables
%               responsible for heteroscedisticity are columns 3 and 5
%               of matrix X, it is possible to use both the sintax
%                    regressH(y,X,X(:,[3 5]))
%               or the sintax
%                    regressH(y,X,[3 5])
%   type:       string specifying the parametric function to be used in the skedastic equation
%               If type is 'arc' (default) than the skedastic function is
%               modelled as follows
%
%               \sigma^2_i = \sigma^2 (1 + exp(\gamma_0 + \gamma_1 Z(i,1) +
%                           ...+ \gamma_{r} Z(i,r)))
%
%               on the other hand, if type is 'har' then traditional
%               formulation due to Harvey is used as follows
%
%               \sigma^2_i = exp(\gamma_0 + \gamma_1 Z(i,1) + ...+
%                           \gamma_{r} Z(i,r)) =\sigma^2 (exp(\gamma_1
%                           Z(i,1) + ...+ \gamma_{r} Z(i,r))
%               
%
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
%   DETAILS. This routine implements Harvey’s (1976) model of
%   multiplicative heteroscedasticity which is a very flexible, general
%   model that includes most of the useful formulations as special cases.
%   The general formulation is
%       \sigma^2_i =\sigma^2 exp(z_i \alpha)
%   Let z_i include a constant term so that z_i'=(1 q_i) where q_i is the
%   original set of variables which are supposed to explain
%   heteroscedasticity. This routine automatically adds a column of 1 to
%   input matrix Z (therefore Z does not have to include a constant term).
%   Now let \gamma'=[ln \sigma^2 \alpha']. Then the model is simply
%       \sigma^2_i = exp(\gamma_' z_i)
%   Once the full parameter vector is estimated \exp( \gamma_1) provides the
%   estimator for \sigma^2
%
% See also regress, regstats
%
% References:
%
%   Greene W.H.(1987): Econometric Analysis (5th edition, section 11.7.1
%   pp. 232-235), (7th edition, section  9.7.1 pp. 280-.
%   Prentice Hall,.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('regressh')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%
%{
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

%
%{
    %The data in Appendix Table F6.1 were used in a study of efficiency in
    %production of airline services in Greene (2007a).
    % See p. 557 of Green (7th edition)

    % Monthly credit card expenditure for 100 individuals.

    % Results in structure "out.Beta" coincide with "Iterated" lines in
    % table 9.2, page 282, 7th edition of Greene (2007).


    load('TableF61_Greene');
    Y=TableF61_Greene.data;

    Q=log(Y(:,4));
    Pfuel=log(Y(:,5));
    Loadfactor=Y(:,6);
    n=size(Y,1);
    X=[Q Q.^2 Pfuel];
    y=log(Y(:,3));


    % Estimate a multiplicative heteroscedastic model and print the
    % estimates of regression and scedastic parameters together with LM, LR
    % and Wald test
    out=regressH(y,X,Loadfactor,'type','har');
%}

%{
    % Estimate a multiplicative heteroscedastic model using just one iteration
    % that is find FGLS estimator (two step estimator)

    % Monthly credit card expenditure for 100 individuals.

    % Results in structure "out" coincides with estimates of row 
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

    % Compare output from Harvey's model with the one of ART
    outHAR=regressH(y,X,[3 4],'msgiter',1,'type','har');
    outART=regressH(y,X,[3 4],'msgiter',1,'type','art');

%}

%% Beginning of code

if nargin>3
% check if input option type exists
    chklist=varargin(1:2:length(varargin));
    
    chktype = find(strcmpi('type',chklist)); 
    if ~isempty(chktype) && strcmp(varargin{2*chktype},'har') ==1
        out=regressHhar(y,X,Z,varargin{:});
    else
        out=regressHart(y,X,Z,varargin{:});
    end
else
     out=regressHart(y,X,Z,varargin);
end

