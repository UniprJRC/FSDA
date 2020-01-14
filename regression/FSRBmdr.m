function [mdrB,Un,BB,BBayes,S2Bayes] = FSRBmdr(y, X, beta0, R, tau0, n0, varargin)
%FSRBmdr computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search.
%
%
%<a href="matlab: docsearchFS('FSRBmdr')">Link to the help function</a>
%
% Required input arguments:
%
%  y:           Response variable. Vector.  A vector with n elements that contains the response variable.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :          Predictor variables. Matrix.  Data matrix of explanatory
%               variables (also called
%               'regressors') of dimension (n x p-1).
%               Rows of X represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
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
%               It can be interpreted as $X_0'X_0$ where $X_0$ is a n0 x p
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
%
% Optional input arguments:
%
%      bsb :   units forming initial subset. Vector.
%                 m x 1 vector containing the units forming initial subset. The
%               default value of bsb is '' (empty value), that is we
%               initialize the search just using prior information.
%                 Example - 'bsb',[3,6,9]
%                 Data Types - double
%  init :       Search initialization. Scalar.
%               It specifies the point where to start monitoring
%               required diagnostics. If it is not specified it is set
%               equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%               The minimum value of init is 0. In this case in the first
%               step we start monitoring at step m=0 (step just based on
%               prior information)
%               Example - 'init',100 starts monitoring from step m=100
%               Data Types - double
% intercept :  Indicator for constant term. true (default) | false. 
%               Indicator for the constant term (intercept) in the fit,
%               specified as the comma-separated pair consisting of
%               'intercept' and either true to include or false to remove
%               the constant term from the model.
%               Example - 'intercept',false
%               Data Types - boolean
%  plots :    Plot on the screen. Scalar.
%               If equal to one a plot of Bayesian minimum deletion residual
%               appears  on the screen with 1%, 50% and 99% confidence
%               bands else (default) no plot is shown.
%                 Example - 'plots',1
%                 Data Types - double
%               Remark: the plot which is produced is very simple. In order
%               to control a series of options in this plot and in order to
%               connect it dynamically to the other forward plots it is necessary to use
%               function mdrplot
%  nocheck:   Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%  msg  :    Level of output to display. Scalar.
%               It controls whether to display or not messages
%               about great interchange on the screen
%               If msg==1 (default) messages are displyed on the screen
%               else no message is displayed on the screen
%               Example - 'msg',1
%               Data Types - double
%   bsbsteps :  steps of the fwd search where to save the units forming subset. Vector.
%               If bsbsteps is 0 we store the units forming
%               subset in all steps. The default is store the units forming
%               subset in all steps if n<5000, else to store the units
%               forming subset at steps init and steps which are multiple
%               of 100. For example, if n=753 and init=6, units forming
%               subset are stored for m=init, 100, 200, 300, 400, 500 and
%               600.
%               Example - 'bsbsteps',[10,20,30]
%               Data Types - double
%
%  Remark:      The user should only give the input arguments that have to
%               change their default value.
%               The name of the input arguments needs to be followed by
%               their value. The order of the input arguments is of no
%               importance.
%
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations. y can be both a row of column vector.
%
% Output:
%
%   mdrB:       n x 2 matrix which contains the monitoring of minimum
%               deletion residual at each step of the forward search.
%               1st col = fwd search index (from 0 to n-1).
%               2nd col = minimum deletion residual.
%  Un:           (n-init) x 11 Matrix which contains the unit(s) included
%               in the subset at each step of the search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one.
%               Un(1,2) for example contains the unit included in step
%               init+1.
%               Un(end,2) contains the units included in the final step
%               of the search.
%  BB:           n x (n-init+1) matrix which contains the units belonging to the
%               subset at each step of the forward search.
%               1st col = index forming subset in the initial step;
%               ...;
%               last column = units forming subset in the final step
%               (i.e. all units).
%  BBayes:       posterior estimates of $\beta$. Matrix.
%               (n-init+1) x (p+1) matrix containing the monitoring o
%               posterior mean of $\beta$ (regression coefficents)
%               $\beta_1 = (c*R + X'X)^{-1} (c*R*\beta_0 + X'y)$.
%  S2Bayes :    posterior estimate of $\sigma^2$ and $\tau$. Matrix.
%               (n-init+1) x 3 matrix containing the monitoring of
%               posterior estimate of $\sigma^2$ and $\tau$
%               in each step of the forward search.
%               1st col = fwd search index (from init to n).
%               2nd col = monitoring of $\sigma^2_1$ (posterior estimate of
%               $\sigma^2$). 
%               3rd col = monitoring $\tau_1$ (posterior estimate of
%               $\tau$). Note that these posterior estimates during the fwd
%               search are corrected using the usual truncated variance
%               consistency factor.
%
% See also: FSRB, FSRmdr, regressB
%
% References:
% Chaloner, K. and Brant, R. (1988), A Bayesian Approach to Outlier
% Detection and Residual Analysis, "Biometrika", Vol. 75, pp. 651-659.
% Riani, M., Corbellini, A. and Atkinson, A.C. (2018), Very Robust Bayesian
% Regression for Fraud Detection, "International Statistical Review",
% http://dx.doi.org/10.1111/insr.12247
% Atkinson, A.C., Corbellini, A. and Riani, M., (2017), Robust Bayesian
% Regression with the Forward Search: Theory and Data Analysis, "Test",
% Vol. 26, pp. 869-886, DOI 10.1007/s11749-017-0542-6
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FSRBmdr')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % FSRBmdr with all default options.
    % Common part to all examples: load Houses Price Dataset.
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
    mdrB=FSRBmdr(y,X,beta0,R,tau0,n0);
%}

%{
    %% FSRBmdr with optional arguments.
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
    mdrB=FSRBmdr(y,X,beta0, R, tau0, n0,'plots',1);
%}

%{
    % Analyze units entering the search in the final steps.
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
    [mdr,Un]=FSRBmdr(y,X,beta0, R, tau0, n0,'plots',1);
%}

%{
    % Units forming subset in each step.
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
    [mdr,Un,BB]=FSRBmdr(y,X,beta0, R, tau0, n0,'plots',1);
%}

%{
    % Monitor $\hat  \beta$.
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
    [mdr,Un,BB,BBayes]=FSRBmdr(y,X,beta0, R, tau0, n0,'plots',1);
%}

%{
    %% Monitor $s^2$.
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
    [mdr,Un,BB,BBayes,S2Bayes]=FSRBmdr(y,X,beta0, R, tau0, n0,'plots',1);
%}

%{
    % Additional example.
    % We change n0 from 5 to 250 and impose FSRBmdr monitoring from step 300.
    load hprice.txt;
    
    % setup parameters
    n=size(hprice,1);
    y=hprice(:,1);
    X=hprice(:,2:5);
    n0=250;

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
    mdrB=FSRBmdr(y,X,beta0, R, tau0, n0,'init',300,'plots',1);
%}



%% Beginning of code

% Input parameters checking

if nargin < 6
    error('FSDA:FSRBmdr:missingInputs','Some Bayesian input (beta0, R, tau0, n0) is missing');
end

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputRB(y,X,nnargin,vvarargin);

%% User options

if n<40
    initdef=0;
else
    initdef=min([3*p+1 floor(0.5*(n+p+1)) n]);
end

% Default for vector bsbsteps which indicates for which steps of the fwd
% search units forming subset have to be saved
bsbstepdef='';

% If n is very large (>500), the step of the search is printed every 500 step
% seq500 is linked to printing
seq500=500*(1:1:ceil(n/500));


bsb='';
options=struct('bsb',bsb,'init',initdef,'intercept',1,...
    'plots',0,'nocheck',0,'msg',1,'bsbsteps',bsbstepdef);



if nargin > 7
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSRBmdr:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end
% And check if the optional user parameters are reasonable.
bsb=options.bsb;

if bsb==0
    bsb=randsample(n,p);
    Xb=X(bsb,:);
    yb=y(bsb);
else
    Xb=X(bsb,:);
    yb=y(bsb,:);
end


% check init
init=options.init;
if  init <0
    fprintf(['Attention : init should be greater or equal than 0. \n',...
        'It is set to 0.']);
    init=0;
elseif init<length(bsb)
    fprintf(['Attention : init1 should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    init=length(bsb);
elseif init>=n
    fprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n-1.']);
    init=n-1;
end


ini0=numel(bsb);

msg=options.msg;
bsbsteps=options.bsbsteps;


%% Initialise key matrices

% sequence from 1 to n.
seq=(1:n)';

% the set complementary to bsb.
ncl=setdiff(seq,bsb);

% The second column of matrix R will contain the Bayesian raw residuals at
% each step of the search
r=[seq zeros(n,1)];


% Matrix BBayes will contain the (Bayesian) beta coefficients in each step of
% the fwd search. The first column of BBayes contains the fwd search index
BBayes=[(init:n)' NaN(n-init+1,p)];     %initial value of beta coefficients is set to NaN

% S2 = (n-init1+1) x 3 matrix which will contain:
% 1st col = fwd search index
% 2nd col = S2= posterior estimate of sigma2
% 3rd col = R^2
S2Bayes=[(init:n)' NaN(n-init+1,2)];        %initial value of S2 (R2) is set to NaN

% mdr = (n-init1-1) x 2 matrix which will contain min deletion residual
% among nobsb r_i^*
mdrB=[(init:n-1)'  NaN(n-init,1)];      %initial value of mdr is set to NaN


% Matrix BB will contain the units forming subset in each step (or in
% selected steps) of the forward search. The first column contains
% information about units forming subset at step init1.
if isempty(bsbsteps)
    % Default for vector bsbsteps which indicates for which steps of the fwd
    % search units forming subset have to be saved
    if n<=5000
        bsbsteps = init:1:n;
    else
        bsbsteps = [init init+100-mod(init,100):100:100*floor(n/100)];
    end
    BB = NaN(n,length(bsbsteps),'single');
elseif bsbsteps==0
    bsbsteps=init:n;
    BB = NaN(n,n-init+1,'single');
else
    if min(bsbsteps)<init
        warning('FSDA:FSMbsb:WrongInit','It is impossible to monitor the subset for values smaller than init');
    end
    bsbsteps=bsbsteps(bsbsteps>=init);
    
    BB = NaN(n,length(bsbsteps),'single');
end

% ij = index which is linked with the columns of matrix BB. During the
% search every time a subset is stored inside matrix BB ij icreases by one
ij=1;

%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init+1:n)' , NaN(n-init,10));

%% Start of the forward search

%mj=1;
for mm=ini0:n
    
    % if n>5000 show every 500 steps the fwd search index
    if msg==1 && n>5000
        if length(intersect(mm,seq500))==1
            disp(['m=' int2str(mm)]);
        end
    end
    
    % call bayesian procedure
    [bayes]=regressB(y, X, beta0, R, tau0, n0, 'bsb', bsb, 'nocheck',1);
    
    % bayesian beta
    b=bayes.beta1;
    
    
    % please note:
    % FS uses standard residuals to advance but
    % at each iteration saves the bayesian residuals
    
    resBSB=yb-Xb*b;
    
    e=y-X*b;    % e = vector of residual for all units using b estimated using subset
    
    r(:,2)=e.^2;
    
    if (mm>=init)
        % Store units belonging to the subset
        if intersect(mm,bsbsteps)==mm
            BB(bsb,ij)=bsb;
            ij=ij+1;
        end
        
        % Stores beta coefficients if there is no rank problem
        BBayes(mm-init+1,2:p+1)=b';
        
        % Sb = posterior estimate of sigma2
        Sb=1/bayes.tau1;
        S2Bayes(mm-init+1,2)=Sb;
        
        % Store R2
        S2Bayes(mm-init+1,3)=1-var(resBSB)/var(yb);
        
        if mm<n
            % Store Bayesian mdr
            mdrB(mm-init+1,2)= min(abs(bayes.res(ncl,2)));
        end
    end
    
    if mm<n
        
        % store units forming old subset in vector oldbsb
        oldbsb=bsb;
        
        % order the r_i
        ord=sortrows(r,2);
        
        % bsb= units forming the new  subset
        bsb=ord(1:(mm+1),1);
        
        Xb=X(bsb,:);  % subset of X
        yb=y(bsb);    % subset of y
        
        if mm>=init
            unit=setdiff(bsb,oldbsb);
            
            % If the interchange involves more than 10 units, store only the
            % first 10.
            if length(unit)<=10
                Un(mm-init+1,2:(length(unit)+1))=unit;
            else
                if msg==1
                    disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                    disp(['Number of units which entered=' int2str(length(unit))]);
                    Un(mm-init+1,2:end)=unit(1:10);
                end
            end
        end
        
        if mm < n-1
            ncl=ord(mm+2:n,1);    % ncl= units forming the new noclean
        end
    end   % if mm<n
    %mj=mj+n0;
end  % for mm=ini0:n loop


% Plot minimum deletion residual with 1%, 50% and 99% envelopes
if options.plots==1
    quant=[0.01;0.5;0.99];
    % Compute theoretical envelops for minimum deletion residual based on all
    % the observations for the above quantiles.
    [gmin] = FSRenvmdr(n,p,'prob',quant);
    plot(mdrB(:,1),mdrB(:,2));
    
    % Superimpose 1%, 99%, 99.9% envelopes based on all the observations
    lwdenv=2;
    % Superimpose 50% envelope
    line(gmin(:,1),gmin(:,3),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
    % Superimpose 1% and 99% envelope
    line(gmin(:,1),gmin(:,2),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
    line(gmin(:,1),gmin(:,4),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
    
    xlabel('Subset size m');
    ylabel('Monitoring of minimum deletion residual');
end
end
%FScategory:REG-Bayes