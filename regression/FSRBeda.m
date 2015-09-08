function [out] = FSRBeda(y, X, varargin)
%FSRBeda enables to monitor several quantities in each step of the Bayesian search
%
%<a href="matlab: docsearchFS('fsrbeda')">Link to the help function</a>
%
% Required input arguments:
%
%   y:          A vector with n elements that contains the response variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%   X :         Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables. Missing values (NaN's) and
%               infinite values (Inf's) are allowed, since observations
%               (rows) with missing or infinite values will automatically
%               be excluded from the computations.
%               PRIOR INFORMATION
%               \beta is assumed to have a normal distribution with
%               mean \beta0 and (conditional on tau0) covariance
%               (1/tau0) (X0'X0)^{-1}
%               \beta~N(    beta0, (1/tau0) (X0'X0)^{-1}    )
%
%
% Optional input arguments:
%
%   intercept :  Indicator for constant term. Scalar.
%                     If 1, a model with constant term will be fitted (default),
%                     if 0, no constant term will be included.
%               Example - 'intercept',1 
%               Data Types - double
%    bayes    : It specifies prior information. Structure.
%               A structure which specifies prior information
%               Strucure bayes contains the following fields
%               bayes.beta0= ( p-times-1 vector containing prior mean of \beta)
%               bayes.R  =  (p-times-p positive definite matrix which can be
%                       interepreted as X0'X0 where X0 is a n0 x p matrix
%                       coming from previous experiments (assuming that the
%                       intercept is included in the model)
%
%               The prior distribution of tau0 is a gamma distribution with
%               parameters a and b, that is
%                     p(tau0) \propto \tau^{a0-1} \exp (-b0 \tau)
%                         E(tau0)= a0/b0
%
%               tau0 (scalar. Prior estimate of tau=1/ \sigma^2 =a0/b0)
%               n0   ( scalar. Sometimes it helps to think of the prior
%                      information as coming from n0 previous experiments.
%                      Therefore we assume that matrix X0 (which defines
%                      R), was made up of n0 observations)
%                  Example - bayes=struct;bayes.R=R;bayes.n0=n0;bayes.beta0=beta0;bayes.tau0=tau0;
%                  Data Types - double
%              REMARK: if structure bayes is not supplied the default
%                      values which are used are
%                      beta0= zeros(p,1)  % vector of zeros
%                      R=eye(p);          % Identity matrix
%                      tau0=1/1e+6;       % Very large value for the
%                                         % prior variance, that is a very
%                                         % small value for tau0
%                      n0=1;              % just one prior observation
%       bsb   : list of units forming the initial subset. Vector.
%                if bsb=0 then the procedure starts with p
%               units randomly chosen else if bsb is not 0 the search will
%               start with m0=length(bsb). The default value of bsb is ''
%               that is in the first step just prior information is used.
%               Example - bsb=[2 5 1];
%               Data Types - double
%        init : Search initialization. Scalar. 
%               scalar, specifies the point where to start monitoring
%               required diagnostics. if init is not specified it will be
%               set equal to :
%                 p+1, if the sample size is smaller than 40;
%                 min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%                   Example - 'init',100 starts monitoring from step m=100 
%                   Data Types - double
%      nocheck: Check input arguments. Scalar.
%               Scalar. If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0. See
%               routine chkinputRB.m for the details of the operations.
%               Example - 'nocheck',1 
%               Data Types - double
%  conflev:   confidence levels to be used to compute HPDI. Vector.
%               This input option is used just if input
%               stats=1. The default value of conflev is [0.95 0.99] that
%               is 95% and 99% HPDI confidence intervals are computed.
%               Example - 'conflev',[0.90 0.93] 
%               Data Types - double
% Remark:       The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%
%               Missing values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations. y can
%               be both a row of column vector.
%
% Output:
%
%   The output consists of a structure 'out' containing the following
%   fields:
%   out.RES=    Scaled residuals. Matrix.
%               n x (n-init+1) matrix containing the monitoring of
%               scaled residuals
%               1st row = residual for first unit ......
%               nth row = residual for nth unit.
%   out.LEV=    Leverage. Matrix.
%               (n+1) x (n-init+1) = matrix containing the monitoring of
%               leverage
%               1st row = leverage for first unit ......
%               nth row = leverage for nth unit.
%    out.BB=    n x (n-init+1) matrix containing the information about the units belonging
%               to the subset at each step of the forward search.
%               1st col = indexes of the units forming subset in the initial step
%               ...
%               last column = units forming subset in the final step (all units)
%   out.mdrB=     n-init x 3 matrix which contains the monitoring of Bayesian
%               minimum deletion residual or (m+1)ordered residual  at each
%               step of the forward search.
%               1st col = fwd search index (from init to n-1)
%               2nd col = minimum deletion residual
%               3rd col = (m+1)-ordered residual
%               Remark: these quantities are stored with sign, that is the
%               min deletion residual is stored with negative sign if
%               it corresponds to a negative residual
%   out.msrB=       n-init+1 x 3 = matrix which contains the monitoring of
%               maximum studentized residual or m-th ordered residual
%               1st col = fwd search index (from init to n)
%               2nd col = maximum studentized residual
%               3rd col = (m)-ordered studentized residual
%   out.Coo=        (n-init+1) x 2 matrix containing the monitoring of Cook or
%               modified Cook distance in each step of the forward search
%               1st col = fwd search index (from init to n)
%               2nd col = monitoring of Cook distance
%  out.beta1 =  (n-init+1) x (p+1) matrix containing the monitoring of
%               posterior mean of \beta (regression coefficents)
%               beta1 = (c*R + X'X)^{-1} (c*R*beta0 + X'y)
%    out.Gam    =  (n-init+1) x 3 matrix containing
%               1st col = fwd search index (from init to n)
%               2nd col = shape parameter a1 of the posterior gamma distribution of tau
%               3rd col = scale parameter b1 of the posterior gamma distribution of tau
%               Remark: a1 = 0.5 (c*n0 + m) where m is subset size
%                       b1 = 0.5 * ( n0 / tau0 + (y-X*beta1)'y +(beta0-beta1)'*c*R*beta0 )%  out.covbeta1=    p x p x (n-init+1) 3D array containing posterior covariance matrix
%               (conditional on tau1) of \beta
%               covbeta1 = (1/tau1) * (c*R + X'X)^{-1}
%               where tau1 is defined as a1/b1 (that is through the gamma
%               parameters of the posterior distribution of \tau)
%               The posterior distribution of \tau is a gamma distribution
%               with parameters a1 and b1
%    out.S21=   (n-init+1) x 3 matrix containing the monitoring of
%               posterior estimate of $\sigma^2$ and $\tau$  
%               in each step of the forward search
%               1st col = fwd search index (from init to n)
%               2nd col = monitoring of $\sigma^2_1$ (posterior estimate of
%               $\sigma^2$)
%               3rd col = monitoring $\tau_1$ (posterior estimate of $\tau$)
%     out.Bpval =   (n-init+1) x (p+1) containing Bayesian p-values.
%               p-value = P(|t| > | \hat \beta se(beta) |)
%               = prob. of beta different from 0
%               1st col = fwd search index (from init to n)
%               2nd col = p-value for first variable
%               ...
%               (p+1) col = p-value for p-th variable
%out.beta1HPD   =  (n-init+1)-by-2-by-p 3D array.
%               Bhpd(:,:,1) = lower and upper HPDI of first element of
%               $\beta_1$ (posterior estimate of $\beta$)
%               ...
%               Bhpd(:,:,p) = lower and upper HPDI of  last element
%               of $\beta_1$  (posterior estimate of $\beta$)
%   out.tau1HPD =  (n-init+1) x 3 containing HPDI for $\tau_1$
%               1st col = fwd search index (from init to n)
%               2nd col = lower value of HPDI 
%               3rd col = upper value of HPDI 
% out.sigma21HPD =  (n-init+1) x 3 containing HPDI for $\sigma^2_1$
%               1st col = fwd search index (from init to n)
%               2nd col = lower value of HPDI 
%               3rd col = upper value of HPDI 
%  out.postodds =   (n-init+1)-by-(p+1) matrix which contains posterior odds for betaj=0
%               For example the posterior odd of beta0=0 is p(y| model which contains
%               all expl variables except the one associated with beta0) divided by
%               p(y| model which contains all expl variables)
%               1st col = fwd search index (from init to n)
%               2nd col = posterior odd for beta1
%               ...
%               (p+1) col = posterior odd for betap
% out.modelprob =  (n-init+1)-by-(p+1) matrix which contains
%               posterior model probability of the model which excludes
%               variable j. For example if modelprob(j)= 0.28, that is if
%               the probability of the model which does not contain
%               variable j is equal to 0.28, it means that there is a 28%
%               chance that beta_j=0 and a 72% chance that it is not.
%               1st col = fwd search index (from init to n)
%               2nd col = posterior model prob of the model which excludes beta1
%               ...
%               (p+1) col = posterior model prob of the model which
%               excludes betap
%    out.Un=        (n-init) x 11 matrix which contains the unit(s)
%               included in the subset at each step of the fwd search
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Un(1,2) for example contains
%               the unit included in step init+1 Un(end,2) contains the
%               units included in the final step of the search
%     out.y=        A vector with n elements that contains the response
%               variable which has been used
%     out.X=       Data matrix of explanatory variables
%               which has been used (it also contains the column of ones if
%               input option intercept was missing or equal to 1)
%out.class =    string FSRBeda.
%
%
% See also FSRB, regressB, FSRBmdr
%
% References:
%
%   Atkinson A.C., Corbellini A. and Riani M. (2015), Robust Bayesian
%   Regression, submitted
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('fsrbeda')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % FSRBeda with all default options.
    % Common part to all examples: load Houses Price Dataset.
    load hprice.txt;
    
    % setup parameters
    n=size(hprice,1);
    y=hprice(:,1);
    X=hprice(:,2:5);
    [out]=FSRBeda(y,X)
%}

%{
    % FSRBeda with optional arguments.
    load hprice.txt;
    
    % setup parameters
    n=size(hprice,1);
    y=hprice(:,1);
    X=hprice(:,2:5);
    
    bayes=struct;
    n0=5;
    bayes.n0=n0;

    % set \beta components
    beta0=0*ones(5,1);
    beta0(2,1)=10;
    beta0(3,1)=5000;
    beta0(4,1)=10000;
    beta0(5,1)=10000;
    bayes.beta0=beta0;

    % \tau
    s02=1/4.0e-8;
    tau0=1/s02;
    bayes.tau0=tau0;

    % R prior settings
    R=2.4*eye(5);
    R(2,2)=6e-7;
    R(3,3)=.15;
    R(4,4)=.6;
    R(5,5)=.6;
    R=inv(R);
    bayes.R=R;
    [out]=FSRBeda(y,X,'bayes',bayes)    
%}

%{
    %% Plot posterior estimates of beta and sigma2.
    % Plot posterior estimates of beta and sigma2 in the interval (subset size) [20 125].
    % In this example for the house price data we monitor the forward plots
    % of the parameters of the linear model adding 95% and 99% HPD
    % regions. The first 8 panels refer to the elements of $\beta_1$
    % and the bottom right-hand panel refer to the estimate of sigma2.
    % The horizontal lines correspond to prior values    
    % The vertical lines refer to the step prior to the introduction of the first outlier
    % init = point to start monitoring diagnostics along the FS
    % run routine FSRB in order to find the outliers automatically
    load hprice.txt;
    
    % setup parameters
    n=size(hprice,1);
    y=hprice(:,1);
    X=hprice(:,2:5);
    
    bayes=struct;
    n0=5;
    bayes.n0=n0;

    % set \beta components
    beta0=0*ones(5,1);
    beta0(2,1)=10;
    beta0(3,1)=5000;
    beta0(4,1)=10000;
    beta0(5,1)=10000;
    bayes.beta0=beta0;

    % \tau
    s02=1/4.0e-8;
    tau0=1/s02;
    bayes.tau0=tau0;

    % R prior settings
    R=2.4*eye(5);
    R(2,2)=6e-7;
    R(3,3)=.15;
    R(4,4)=.6;
    R(5,5)=.6;
    R=inv(R);
    bayes.R=R;       
    outBA=FSRB(y,X,'bayes',bayes', 'plots',0);
    dout=n-length(outBA.ListOut);
    
    % init = initial point to start monitoring
    init=20;
    xlimL=init; % lower value of xlim
    xlimU=125;  % upper value of xlim 

    outBAeda=FSRBeda(y,X,'bayes',bayes,'init',init, 'conflev', [0.95 0.99]);

    % Set font size, line width and line style
    figure;
    lwd=2.5;
    FontSize=14;
    linst={'-','--',':','-.','--',':'};

    for j=1:5
        my_subplot=subplot(3,2,j);
        hold('on')
        % plot 95% and 99% HPD  trajectories
        plot(outBAeda.beta1(:,1),outBAeda.beta1HPD(:,1:2,j),'LineStyle',linst{4},'LineWidth',lwd,'Color','b')
        plot(outBAeda.beta1(:,1),outBAeda.beta1HPD(:,3:4,j),'LineStyle',linst{4},'LineWidth',lwd,'Color','r')

        % plot posterior estimate of beta1_j
        plot(outBAeda.beta1(:,1),outBAeda.beta1(:,j+1)','LineStyle',linst{1},'LineWidth',lwd,'Color','k')

        % Add the horizontal line which corresponds to prior values
        xL = get(my_subplot,'XLim');
        line(xL,[beta0(j) beta0(j)],'Color','r','LineWidth',lwd);

        % Set ylim
        ylimU=max([outBAeda.beta1HPD(:,4,j); beta0(j)]);
        ylimL=min([outBAeda.beta1HPD(:,3,j); beta0(j)]);
        ylim([ylimL ylimU])

        % Add vertical line in correspondence of the step prior to the
        % entry of the first outlier
        line([dout; dout],[ylimL; ylimU],'Color','r','LineWidth',lwd);

        % Set xlim
        xlim([xlimL xlimU]);


        ylabel(['$\hat{\beta_' num2str(j-1) '}$'],'Interpreter','LaTeX','FontSize',20,'rot',-360);
        set(gca,'FontSize',FontSize);
        if j>4
            xlabel('Subset size m','FontSize',FontSize);
        end
    end

    % Subplot associated with the monitoring of sigma^2
    subplot(3,2,6);
    %figure()
    hold('on')
    % 99%
    plot(outBAeda.sigma21HPD(:,1),outBAeda.sigma21HPD(:,4:5),'LineStyle',linst{4},'LineWidth',lwd,'Color','r')
    % 95%
    plot(outBAeda.sigma21HPD(:,1),outBAeda.sigma21HPD(:,2:3),'LineStyle',linst{2},'LineWidth',lwd,'Color','b')
    % Plot 1\/tau1
    plot(outBAeda.S21(:,1),1./outBAeda.S21(:,3),'LineWidth',lwd,'Color','k')
    ylabel('$\hat{\sigma}^2$','Interpreter','LaTeX','FontSize',20,'rot',-360);
    set(gca,'FontSize',FontSize);

    % Set ylim
    ylimU=max([outBAeda.sigma21HPD(:,5); s02]);
    ylimL=min([outBAeda.sigma21HPD(:,4); s02]);
    ylim([ylimL ylimU])

    % Set xlim
    xlim([xlimL xlimU]);
    xL = get(my_subplot,'XLim');

    % Add the horizontal line which corresponds to prior value of $\sigma^2$
    line(xL,[s02 s02],'Color','r','LineWidth',lwd);

    % Add vertical line in correspondence of the step prior to the
    % entry of the first outlier
    line([dout; dout],[ylimL; ylimU],'Color','r','LineWidth',lwd);
    xlabel('Subset size m','FontSize',FontSize);

    % Add multiple title
    suplabel(['Housing data; forward plots in the interval ['...
        num2str(xlimL) ',' num2str(xlimU) ...
        ']  of HPD regions for \beta and \sigma^2'],'t')
%}

%{
    %% Plot posterior estimates of beta and sigma2 in the last steps.
    % Plot posterior estimates of beta and sigma2 in the interval (subset size) [250 n+1].
    % In this example for the house price data we monitor the forward plots
    % of the parameters of the linear model adding 95% and 99% HPD
    % regions. The first 8 panels refer to the elements of $\beta_1$
    % and the bottom right-hand panel refer to the estimate of sigma2.
    % The horizontal lines correspond to prior values    
    % The vertical lines refer to the step prior to the introduction of the first outlier
    % init = point to start monitoring diagnostics along the FS
    % run routine FSRB in order to find the outliers automatically
    load hprice.txt;
    
    % setup parameters
    n=size(hprice,1);
    y=hprice(:,1);
    X=hprice(:,2:5);
    
    bayes=struct;
    n0=5;
    bayes.n0=n0;

    % set \beta components
    beta0=0*ones(5,1);
    beta0(2,1)=10;
    beta0(3,1)=5000;
    beta0(4,1)=10000;
    beta0(5,1)=10000;
    bayes.beta0=beta0;

    % \tau
    s02=1/4.0e-8;
    tau0=1/s02;
    bayes.tau0=tau0;

    % R prior settings
    R=2.4*eye(5);
    R(2,2)=6e-7;
    R(3,3)=.15;
    R(4,4)=.6;
    R(5,5)=.6;
    R=inv(R);
    bayes.R=R;
    outBA=FSRB(y,X,'bayes',bayes', 'plots',0);
    dout=n-length(outBA.ListOut);
    
    % init = initial point to start monitoring
    init=250;
    xlimL=init; % lower value fo xlim
    xlimU=n+1;  % upper value of xlim 

    outBAeda=FSRBeda(y,X,'bayes',bayes,'init',init, 'conflev', [0.95 0.99]);

    % Set font size, line width and line style
    figure;
    lwd=2.5;
    FontSize=14;
    linst={'-','--',':','-.','--',':'};

    for j=1:5
        my_subplot=subplot(3,2,j);
        hold('on')
        % plot 95% and 99% HPD  trajectories
        plot(outBAeda.beta1(:,1),outBAeda.beta1HPD(:,1:2,j),'LineStyle',linst{4},'LineWidth',lwd,'Color','b')
        plot(outBAeda.beta1(:,1),outBAeda.beta1HPD(:,3:4,j),'LineStyle',linst{4},'LineWidth',lwd,'Color','r')

        % plot posterior estimate of beta1_j
        plot(outBAeda.beta1(:,1),outBAeda.beta1(:,j+1)','LineStyle',linst{1},'LineWidth',lwd,'Color','k')

        % Add the horizontal line which corresponds to prior values
        xL = get(my_subplot,'XLim');
        line(xL,[beta0(j) beta0(j)],'Color','r','LineWidth',lwd);

        % Add vertical line in correspondence of the step prior to the
        % entry of the first outlier
        line([dout; dout],[ylimL; ylimU],'Color','r','LineWidth',lwd);

        % Set ylim
        ylimU=max([outBAeda.beta1HPD(:,4,j); beta0(j)]);
        ylimL=min([outBAeda.beta1HPD(:,3,j); beta0(j)]);
        ylim([ylimL ylimU])

        % Set xlim
        xlim([xlimL xlimU]);


        ylabel(['$\hat{\beta_' num2str(j-1) '}$'],'Interpreter','LaTeX','FontSize',20,'rot',-360);
        set(gca,'FontSize',FontSize);
        if j>4
            xlabel('Subset size m','FontSize',FontSize);
        end
    end

    % Subplot associated with the monitoring of sigma^2
    subplot(3,2,6);
    %figure()
    hold('on')
    % 99%
    plot(outBAeda.sigma21HPD(:,1),outBAeda.sigma21HPD(:,4:5),'LineStyle',linst{4},'LineWidth',lwd,'Color','r')
    % 95%
    plot(outBAeda.sigma21HPD(:,1),outBAeda.sigma21HPD(:,2:3),'LineStyle',linst{2},'LineWidth',lwd,'Color','b')
    % Plot 1/tau1
    plot(outBAeda.S21(:,1),1./outBAeda.S21(:,3),'LineWidth',lwd,'Color','k')
    ylabel('$\hat{\sigma}^2$','Interpreter','LaTeX','FontSize',20,'rot',-360);
    set(gca,'FontSize',FontSize);

    % Set ylim
    ylimU=max([outBAeda.sigma21HPD(:,5); s02]);
    ylimL=min([outBAeda.sigma21HPD(:,4); s02]);
    ylim([ylimL ylimU])

    % Set xlim
    xlim([xlimL xlimU]);
    xL = get(my_subplot,'XLim');

    % Add the horizontal line which corresponds to prior value of $\sigma^2$
    line(xL,[s02 s02],'Color','r','LineWidth',lwd);

    % Add vertical line in correspondence of the step prior to the
    % entry of the first outlier
    line([dout; dout],[ylimL; ylimU],'Color','r','LineWidth',lwd);
    xlabel('Subset size m','FontSize',FontSize);

    % Add multiple title
    suplabel(['Housing data; forward plots in the interval ['...
        num2str(xlimL) ',' num2str(xlimU) ...
        ']  of HPD regions for \beta and \sigma^2'],'t')
%}

%{
    % Plot of HPDI for BankProfit data.
    XX=load('BankProfit.txt');
    R=load('BankProfitR.txt');

    X=XX(:,1:end-1);
    y=XX(:,end);
    % Load prior information
    beta0=zeros(10,1);
    beta0(1,1)=-0.5;        % 
    beta0(2,1)=9.1;         % Number of products (NUMPRO)
    beta0(3,1)=0.001;       % direct revenues  (DIRREV)
    beta0(4,1)=0.0002;      % indirect revenues   (INDREV)
    beta0(5,1)=0.002;       %  savings accounts   SAVACC
    beta0(6,1)=0.12;        %  number of operations   NUMOPE
    beta0(7,1)=0.0004;      %  total amount of operations  TOTOPE
    beta0(8,1)=-0.0004;     %  Bancomat POS
    beta0(9,1)=1.3;         %  Number of cards   NUMCAR
    beta0(10,1)=0.00004;    %  Amount in cards   TOTCAR

    % \tau
    s02=10000;
    tau0=1/s02;


    % number of obs in which prior was based
    n0=1500;

    bayes=struct;
    bayes.R=R;
    bayes.n0=n0;
    bayes.beta0=beta0;
    bayes.tau0=tau0;
    n=length(y);

    % run routine FSRB in order to find the outliers automatically
    outBA=FSRB(y,X,'bayes',bayes', 'plots',0);
    dout=n-length(outBA.ListOut);

    % init = point to start monitoring diagnostics along the FS
    init=1700;
    xlimL=init;
    xlimU=n+1;

    outBAeda=FSRBeda(y,X,'bayes',bayes,'init',init, 'conflev', [0.95 0.99]);
    % Set font size, line width and line style
    figure;
    lwd=2.5;
    FontSize=14;
    linst={'-','--',':','-.','--',':'};
    nr=4;
    nc=3;
    for j=1:length(beta0)
        my_subplot=subplot(nr,nc,j);
        hold('on')
        % plot 95% and 99% HPD  trajectories
        plot(outBAeda.beta1(:,1),outBAeda.beta1HPD(:,1:2,j),'LineStyle',linst{4},'LineWidth',lwd,'Color','b')
        plot(outBAeda.beta1(:,1),outBAeda.beta1HPD(:,3:4,j),'LineStyle',linst{4},'LineWidth',lwd,'Color','r')

        % plot posterior estimate of beta1_j
        plot(outBAeda.beta1(:,1),outBAeda.beta1(:,j+1)','LineStyle',linst{1},'LineWidth',lwd,'Color','k')

        % Add the horizontal line which corresponds to prior values
        xL = get(my_subplot,'XLim');
        line(xL,[beta0(j) beta0(j)],'Color','r','LineWidth',lwd);


        % Set ylim
        ylimU=max([outBAeda.beta1HPD(:,4,j); beta0(j)]);
        ylimL=min([outBAeda.beta1HPD(:,3,j); beta0(j)]);
        ylim([ylimL ylimU])

        % Set xlim
        xlim([xlimL xlimU]);

            % Add vertical line in correspondence of the step prior to the
        % entry of the first outlier
        line([dout; dout],[ylimL; ylimU],'Color','r','LineWidth',lwd);

        ylabel(['$\hat{\beta_' num2str(j-1) '}$'],'Interpreter','LaTeX','FontSize',20,'rot',-360);
        set(gca,'FontSize',FontSize);
        if j>nr*(nc-1)
            xlabel('Subset size m','FontSize',FontSize);
        end
    end

    % Subplot associated with the monitoring of sigma^2
    subplot(4,3,11);
    %figure()
    hold('on')
    % 99%
    plot(outBAeda.sigma21HPD(:,1),outBAeda.sigma21HPD(:,4:5),'LineStyle',linst{4},'LineWidth',lwd,'Color','r')
    % 95%
    plot(outBAeda.sigma21HPD(:,1),outBAeda.sigma21HPD(:,2:3),'LineStyle',linst{2},'LineWidth',lwd,'Color','b')
    % Plot 1\/tau1
    plot(outBAeda.S21(:,1),1./outBAeda.S21(:,3),'LineWidth',lwd,'Color','k')
    ylabel('$\hat{\sigma}^2$','Interpreter','LaTeX','FontSize',20,'rot',-360);
    set(gca,'FontSize',FontSize);

    % Set ylim
    ylimU=max([outBAeda.sigma21HPD(:,5); s02]);
    ylimL=min([outBAeda.sigma21HPD(:,4); s02]);
    ylim([ylimL ylimU])

    ylimL=8000;
    ylimU=14000;
    ylim([ylimL ylimU])

    % Set xlim
    xlim([xlimL xlimU]);
    xL = get(my_subplot,'XLim');

    % Add the horizontal line which corresponds to prior value of $\sigma^2$
    line(xL,[s02 s02],'Color','r','LineWidth',lwd);

    % Add vertical line in correspondence of the step prior to the
    % entry of the first outlier
    line([dout; dout],[ylimL; ylimU],'Color','r','LineWidth',lwd);
    xlabel('Subset size m','FontSize',FontSize);


    % Add multiple title
    suplabel(['Bank profit data; forward plots in the interval ['...
        num2str(xlimL) ',' num2str(xlimU) ...
        ']  of HPD regions for \beta and \sigma^2'],'t')
%}

%% Input parameters checking

if nargin>6
    
    % chklist = vector containing the names of optional arguments
    chklist=varargin(1:2:length(varargin));
    
    % chkint gives the position of the option nocheck in vector chklist
    % It is empty if the option nocheck is not specified by the user
    chkint=find(strcmp('nocheck',chklist));
else
    chkint='';
end

% If nocheck=1 skip checks on y and X
if ~isempty(chkint) && cell2mat(varargin(2*chkint))==1;
    [n,p]=size(X);
else
    
    
    
    % Now n is the new number of non missing observations
    n=length(y);
    
    
    % Now check if the intercept has to be added
    if nargin > 6
        
        % chklist = vector containing the names of optional arguments
        chklist=varargin(1:2:length(varargin));
        
        % chkint is non empty if the user has specified the option intercept
        % chkint gives the position of the option intercept in vector chklist
        % It is empty if the option intercept is not specified by the user
        chkint=find(strcmp('intercept',chklist));
        
        
        % if a value for the intercept has been specified
        % and this value is equal to 1
        % then add the column of ones to matrix X
        if isempty(chkint) || cell2mat(varargin(2*chkint))==1
            X=cat(2,ones(n,1),X); % add column of ones
        end
    else
        % If the user has not specified a value for the intercept than the
        % column of ones is automatically attached
        X=cat(2,ones(n,1),X); % add column of ones
    end
    
    % p is the number of parameters to be estimated
    p=size(X,2);
    
end

%% User options

% init = scalar which specifies where to start monitoring required statistics
if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));
end

% Initialize bsb as empty. That is in the initial step, the estimate is
% just based on prior values
bsb='';

% ini0=init;
bayesdef=struct;
bayesdef.beta0=zeros(p,1);
bayesdef.R=eye(p);
bayesdef.tau0=1/1e+6;
bayesdef.n0=1;


conflevdef=[0.95 0.99];
options=struct('intercept',1,'init',init,'bayes',bayesdef,'bsb',bsb,...
    'nocheck',0,'conflev',conflevdef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRBeda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

if nargin > 3
    
    % We now overwrite inside structure options the default values with
    % those chosen by the user
    % Notice that in order to do this we use dynamic field names
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end
conflev=options.conflev;
bayes=options.bayes;

beta0=bayes.beta0;
R=bayes.R;
tau0=bayes.tau0;
n0=bayes.n0;

bsb=options.bsb;
if bsb==0
    bsb=randsample(n,p);
    Xb=X(bsb,:);
else
    Xb=X(bsb,:);
end

% is bsb is empty ini0 =0
ini0=length(bsb);

% check init
init=options.init;
if init<ini0;
    mess=sprintf(['Attention : init should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    disp(mess);
    init=ini0;
elseif init>=n;
    mess=sprintf(['Attention : init should be smaller than n. \n',...
        'It is set to n-1.']);
    disp(mess);
    init=n-1;
end

%% Declare matrices to store quantities

% sequence from 1 to n
seq=(1:n)';

% complementary of bsb
ncl=setdiff(seq,bsb);

% The second column of matrix R will contain the OLS residuals
% at each step of the forward search
r=[seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100=100*(1:1:ceil(n/100));

zer=NaN(n-init,2);
zer1=NaN(n-init+1,2);

% Matrix beta1 will contain posterior beta coeff. in each step of the fwd search
% The first column will contain the fwd search index
beta1=[(init:n)' NaN(n-init+1,p)];

% covbeta1 3D array will contain covariance matrix of beta1
covbeta1=NaN(p,p,n-init+1);

%    Gam    :   (n-init+1) x 3 matrix containing
%               1st col = fwd search index (from init to n)
%               2nd col = parameter a1 of the posterior gamma distribution of tau
%               3rd col = parameter b1 of the posterior gamma distribution of tau
%               Remark: a1 = 0.5 (c*n0 + m) where m is subset size
%                       b1 = 0.5 * ( n0 / tau0 + (y-X*beta1)'y +(beta0-beta1)'*c*R*beta0 )
Gam=beta1(:,1:3);


% S21=(n-init+1) x 3  matrix which will contain
% 1st col = fwd search index
% 2nd col = posterior estimate of $\sigma^2$
% 3rd col = posterior estimate of $\tau$
S21=[(init:n)' zer1];

% mdr= (n-init) x 3 matrix
% 1st column = fwd search index
% 2nd col min deletion residual among observations non belonging to the
% subset
mdr=[(init:n-1)'  zer(:,1)];

% mdr= (n-init+1) x 3 matrix which will contain max studentized residual
%  among bsb and m-th studentized residual
msr=[(init:n)'  zer1(:,1)];

% coo= (n-init) x 2 matrix which will contain Cook distances
%  (2nd col)
coo=[((init+1):n)'  NaN(n-init,1)];


% Matrix RES will contain the resisuals for each unit in each step of the forward search
% The first row refers to the residuals of the first unit
RES=NaN(n,n-init+1);
RES(:)=NaN;

% Matrix BB will contain the units forming subset in each step of the forward search
% The first column contains the units forming subset at step init
% The first row is associated with the first unit
BB=RES;

% Matrix LEV will contain the leverage of the units forming subset in each step of the forward search
% The first column contains the leverage associated with the units forming subset at step init
% The first row is associated with the first unit
LEV=RES;

%  Un= Matrix whose 2nd column:11th col contain the unit(s) just included
Un=NaN(n-init,10);
Un=[(init+1:n)' Un];

% Bpval will contain (n-init+1) x (p+1) containing Bayesian p-values.
Bpval=beta1;
% Bhpd will contain (n-init+1)-by-2-by-p 3D array.
beta1HPD=NaN(n-init+1,2*length(conflev),p);
tau1HPD=[(init:n)' zeros(n-init+1,2*length(conflev))];
sigma21HPD=tau1HPD;

% postodds will contain (n-init+1)-by-(p+1) matrix which contains posterior
% odds for betaj=0
postodds=beta1;
% modelprob will contain (n-init+1)-by-(p+1) matrix which contains which
% contains posterior model probability of the model which excludes variable j.
modelprob=beta1;

%% Start of the forward search
for mm=ini0:n;
    
    % if n>200 show every 100 steps the fwd search index
    if n>200;
        if length(intersect(mm,seq100))==1;
            disp(['m=' int2str(mm)]);
        end
    end
    
    
    %b=Xb\yb;
    % call bayesian procedure
    if mm>=init
        [bayes]=regressB(y, X, beta0, R, tau0, n0, 'bsb', bsb,'stats',1,'conflev',conflev, 'nocheck',1);
    else
        [bayes]=regressB(y, X, beta0, R, tau0, n0, 'bsb', bsb,'stats',0, 'nocheck',1);
    end
    
    % why not bayes.res?
    b=bayes.beta1;
    
    if (mm>=init);
        
        % Store Units belonging to the subset
        BB(bsb,mm-init+1)=bsb;
        
        
        % Store posterior beta coefficients
        beta1(mm-init+1,2:p+1)=b';
        
        
        % Store leverage for the units belonging to subset
        % hi contains leverage for all units
        % It is a proper leverage for the units belonging to susbet
        % It is a pseudo leverage for the unit not belonging to the subset
        mAm=bayes.covbeta1;
        
        mmX=bayes.tau1*mAm;
        % Leverage for the units belonging to subset
        hi=sum((Xb*mmX).*Xb,2);
        
        LEV(bsb,mm-init+1)=hi;
        
        % Store cov matrix of beta
        covbeta1(:,:,mm-init+1)=bayes.covbeta1;
        
        % Store paramters a and b associated with the posterior estimate of tau
        Gam(mm-init+1,2:3)=[bayes.a1 bayes.b1];
        
        
        % Store p-values
        Bpval(mm-init+1,:)=[mm bayes.Bpval'];
        
        % Store HPD for $\beta_1$, $\tau_1$ and $\sigma^2_1$ 
        beta1HPD(mm-init+1,:,:)=bayes.beta1HPD';
        tau1HPD(mm-init+1,2:end)=bayes.tau1HPD;
        sigma21HPD(mm-init+1,2:end)=bayes.sigma21HPD;
        
        %
        postodds(mm-init+1,:)=[mm bayes.postodds'];
        modelprob(mm-init+1,:)=[mm bayes.modelprob'];
    end
    
    
    % store res. sum of squares/(mm-k)
    % Store estimate of \sigma^2 using units forming subset
    Sb=1/bayes.tau1;
    
    
    % e= vector of residuals for all units using b estimated using subset
    % e= y-X*b;
    e=bayes.res(:,1);
    
    if (mm>=init)
        % Store all residuals
        RES(:,mm-init+1)=e;
        
        % Store posterior estimates of sigma2 an tau just using the units
        % belonging to subset
        S21(mm-init+1,2:3)=[bayes.sigma21 bayes.tau1];
        resBSB=bayes.res(bsb,1);
        
        % Store maximum studentized residual
        % among the units belonging to the subset
        msrsel=sort(abs(resBSB)./sqrt(Sb*(1-hi)));
        msr(mm-init+1,2)=msrsel(mm);
        
    end
    
    r(:,2)=e.^2;
    
    if mm>init;
        % Store in the second column of matrix coo the Cook
        % distance
        bib=beta1(mm-init+1,2:p+1)-beta1(mm-init,2:p+1);
        
        coo(mm-init,2)=bib*(mAm\(bib'))/p;
        % bib*inv(bayes.covbeta1)*(bib')/(p)
        
        if length(unit)>5;
            unit=unit(1:5);
        end
    end
    
    if mm<n;
        if mm>=init;
            % ord = matrix whose first col (divided by S2(i)) contains the deletion residuals
            % for all units. For the units belonging to the subset these are proper deletion residuals
            % ord = [(r(:,2)./(1+hi)) e];
            ord = abs(bayes.res(ncl,2));
            
            % Store minimum deletion residual in 2nd col of matrix mdr
            selmdr=sort(ord,1);
            
            mdr(mm-init+1,2)=selmdr(1);
        end
        
        % store units forming old subset in vector oldbsb
        oldbsb=bsb;
        
        % order the r_i and include the smallest among the units
        %  forming the group of potential outliers
        ord=sortrows(r,2);
        
        % bsb= units forming the new  subset
        bsb=ord(1:(mm+1),1);
        
        Xb=X(bsb,:);  % subset of X
        
        if mm>=init;
            unit=setdiff(bsb,oldbsb);
            if length(unit)<=10
                Un(mm-init+1,2:(length(unit)+1))=unit;
            else
                disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                Un(mm-init+1,2:end)=unit(1:10);
            end
        end
        
        
        if mm < n-1;
            % ncl= units forming the new noclean
            ncl=ord(mm+2:n,1);
        end
    end
    
end

%% Structure returned by function FSReda
out=struct;
% Scale residuals multiplying by the final posterior estimate of $\tau$
out.RES=RES*sqrt(S21(end,3));
out.LEV=LEV;
out.BB=BB;
out.mdr=mdr;
out.msr=msr;
out.beta1=beta1;
out.covbeta1=covbeta1;
out.Gam=Gam;
out.S21=S21;
out.coo=coo;
out.Bpval=Bpval;
out.beta1HPD=beta1HPD;
out.tau1HPD=tau1HPD;
out.sigma21HPD=sigma21HPD;
out.postodds=postodds;
out.modelprob=modelprob;
out.Un=Un;
out.y=y;
out.X=X;
out.class='FSRBeda';
end

