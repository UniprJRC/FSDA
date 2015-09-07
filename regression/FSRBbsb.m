function [Un,BB] = FSRBbsb(y, X, beta0, R, tau0, n0, varargin)
%FSRBbsb returns the units belonging to the subset in each step of the Bayesian forward search
%
%
%<a href="matlab: docsearch('FSRBbsb')">Link to the help function</a>
%
% Required input arguments:
%
%  y:            A vector with n elements that contains the response variable.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :           Data matrix of explanatory variables (also called
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
%  intercept :   Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%               Example - 'intercept',1
%               Data Types - double
%   plots   :   Plot on the screen. Scalar. 
%               If plots=1 the monitoring units plot is displayed on the
%               screen. The default value of plots is 0 (that is no plot
%               is produced on the screen).
%               Example - 'plots',1
%               Data Types - double
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
%               If msg==1 (default) messages are displayed on the screen
%               else no message is displayed on the screen
%               Example - 'msg',1
%               Data Types - double
%
% Output:
%
%  Un:          Units included in each step. Matrix.
%               (n-init) x 11 Matrix which contains the unit(s) included
%               in the subset at each step of the search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one.
%               Un(1,2) for example contains the unit included in step
%               init+1.
%               Un(end,2) contains the units included in the final step
%               of the search.
%  BB:          Units belonging to subset in each step. Matrix.
%               n x (n-init+1) matrix which contains the units belonging to the
%               subset at each step of the forward search.
%               1st col = index forming subset in the initial step
%               ...
%               last column = units forming subset in the final step 
%               (i.e. all units).
%
% See also FSRbsb, FSRHbsb
%
% References:
% Chaloner and Brant (1988). A Bayesian Approach to Outlier Detection and
% Residual Analysis, Biometrika, Vol 75 pp. 651-659.
% Riani M., Corbellini A., Atkinson A.C. (2015), Very Robust Bayesian
% Regression for Fraud Detection, submitted
% Atkinson A.C., Corbellini A., Riani M., (2015) Robust Bayesian
% Regression, submitted
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearch('FSRBbsb')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    %% FSRBbsb with all default options.
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

    % $\tau$
    s02=1/4.0e-8;
    tau0=1/s02;

    % R prior settings
    R=2.4*eye(5);
    R(2,2)=6e-7;
    R(3,3)=.15;
    R(4,4)=.6;
    R(5,5)=.6;
    R=inv(R);
    [~,Un,BB]=FSRBmdr(y,X,beta0, R, tau0, n0);
    [Unchk,BBchk]=FSRBbsb(y,X,beta0, R, tau0, n0);
    % Test for equality BB and BBchk
    disp(isequaln(BB,BBchk))
    % Test for equality Un and Unchk
    disp(isequaln(Un,Unchk))
%}

%{
    %% Display the monitoring units plot.
    % Suppress all messages about interchange with option msg 
    [Unchk,BBchk]=FSRBbsb(y,X,beta0,R,tau0,n0,'plots',1,'msg',0);
%}

%% Input parameters checking

if nargin < 6
    error('FSDA:FSRBmdr:missingInputs','Some Bayesian input (beta0, R, tau0, n0) is missing');
end

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputRB(y,X,nnargin,vvarargin);

%% User options

if n<40
    initdef=p+1;
else
    initdef=min(3*p+1,floor(0.5*(n+p+1)));
end


% If n is very large (>500), the step of the search is printed every 500 step
% seq500 is linked to printing
seq500=500*(1:1:ceil(n/500));


bsb='';
options=struct('bsb',bsb,'init',initdef,'intercept',1,...
    'plots',0,'nocheck',0,'msg',1);



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
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end
% And check if the optional user parameters are reasonable.
bsb=options.bsb;

if bsb==0;
    bsb=randsample(n,p);
end


% check init
init=options.init;
if  init <0;
    mess=sprintf(['Attention : init1 should be greater or equal than 0. \n',...
        'It is set to 0.']);
    disp(mess);
    init=0;
elseif init<length(bsb);
    mess=sprintf(['Attention : init1 should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    disp(mess);
    init=length(bsb);
elseif init>=n;
    mess=sprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n-1.']);
    disp(mess);
    init=n-1;
end


ini0=numel(bsb);

msg=options.msg;
intercept=options.intercept;

%% Initialise key matrices

% sequence from 1 to n.
seq=(1:n)';

% The second column of matrix R will contain the Bayesian raw residuals at
% each step of the search
r=[seq zeros(n,1)];


% Matrix BB will contain the units forming subset in each step of the
% forward search. The first column contains information about units forming subset at
% step init1.
    BB = NaN(n,n-init+1);


%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init+1:n)' , NaN(n-init,10));

%% Start of the forward search

for mm=ini0:n;
    
    % if n>5000 show every 500 steps the fwd search index
    if msg==1 && n>5000;
        if length(intersect(mm,seq500))==1;
            disp(['m=' int2str(mm)]);
        end
    end
    
    % call bayesian procedure
    [bayes]=regressB(y, X(:,2:end), beta0, R, tau0, n0, 'bsb', bsb,'intercept',intercept);
    
    % bayesian beta
    b=bayes.beta1;
    
    e=y-X*b;    % e = vector of residual for all units using b estimated using subset
    
    r(:,2)=e.^2;
    
    if (mm>=init);
        % Store units belonging to the subset
                          BB(bsb,mm-init+1)=bsb;
    end
    
    if mm<n;
        
        % store units forming old subset in vector oldbsb
        oldbsb=bsb;
        
        % order the r_i
        ord=sortrows(r,2);
        
        % bsb= units forming the new  subset
        bsb=ord(1:(mm+1),1);
        
        if mm>=init;
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
        
    end   % if mm<n
    %mj=mj+n0;
end  % for mm=ini0:n loop

plots=options.plots;

if plots==1
    % Create the 'monitoring units plot'
    figure;
    seqr=[Un(1,1)-1; Un(:,1)];
    plot(seqr,BB','bx');
    xlabel('Subset size m');
    ylabel('Monitoring units plot');
end
end