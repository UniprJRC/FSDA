function out = regressts(y,varargin)
%regressts computes estimates of regression parameters for a time series models
%
%<a href="matlab: docsearchFS('regressts')">Link to the help function</a>
%
% Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :          Predictor variables. Matrix. Matrix of explanatory
%               variables (also called 'regressors') of dimension n x (p-1)
%               where p denotes the number of explanatory variables
%               including the intercept.
%               Rows of X represent observations, and columns represent
%               variables. By default, there is a constant term in the
%               model, unless you explicitly remove it using input option
%               intercept, so do not include a column of 1s in X. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations.
%  bsb     :    list of units forming the initial subset. Vector. If bsb=0
%               (default) then the procedure starts with p units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb)
%
% Optional input arguments:
%
%  init :       Search initialization. Scalar.
%               It specifies the point where to initialize the search and
%               start monitoring required diagnostics. If it is not
%               specified it is set equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%               Example - 'init',100 starts monitoring from step m=100
%               Data Types - double
%  intercept :  Indicator for constant term. Scalar. If 1, a model with
%               constant term will be fitted (default), else no constant
%               term will be included.
%               Example - 'intercept',1
%               Data Types - double
%  plots :      Plot on the screen. Scalar. If equal to one a plot of
%               minimum deletion residual appears  on the screen with 1%,
%               50% and 99% confidence bands else (default) no plot is
%               shown.
%               Example - 'plots',1
%               Data Types - double
%               Remark: the plot which is produced is very simple. In order
%               to control a series of options in this plot and in order to
%               connect it dynamically to the other forward plots it is
%               necessary to use function mdrplot.
%  nocheck:     Check input arguments. Scalar. If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additioanl column of ones for
%               the intercept is not added. As default nocheck=0. The
%               controls on h, alpha and nsamp still remain
%               Example - 'nocheck',1
%               Data Types - double
%  msg  :       Level of output to display. Scalar. It controls whether to
%               display or not messages about great interchange on the
%               screen If msg==1 (default)
%               messages are displayed on the screen
%               else no message is displayed on the screen
%               Example - 'msg',1
%               Data Types - double
%  constr :     Constrained search. Vector. r x 1 vector which contains the list of units which are
%               forced to join the search in the last r steps. The default
%               is constr=''.  No constraint is imposed
%               Example - 'constr',[1:10] forces the first 10 units to join
%               the subset in the last 10 steps
%               Data Types - double
% bsbmfullrank :What to do in case subset at step m (say bsbm) produces a
%               non singular X. Scalar.
%               This options controls what to do when rank(X(bsbm,:)) is
%               smaller then number of explanatory variables.
%               If bsbmfullrank = 1 (default is 1) these units (whose number
%               is say mnofullrank) are constrained to enter the search in
%               the final n-mnofullrank steps else the search continues
%               using as estimate of beta at step m the estimate of beta
%               found in the previous step.
%               Example - 'bsbmfullrank',1
%               Data Types - double
%      bsb :   units forming subset. Vector.
%                m x 1 vector.
%               The default value of bsb is 1:n1, that is all n1 units are
%               used to compute beta1
%               Example - 'bsb',[3 5]
%               Data Types - double
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
%  mdr:          n -init x 2 matrix which contains the monitoring of minimum
%               deletion residual at each step of the forward search.
%               1st col = fwd search index (from init to n-1).
%               2nd col = minimum deletion residual.
%               REMARK: if in a certain step of the search matrix is
%               singular, this procedure checks ohw many observations
%               produce a singular matrix. In this case mdr is a column
%               vector which contains the list of units for which matrix X
%               is non singular.
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
%               n x (n-init+1) or n-by-length(bsbsteps) matrix (depending on input
%               option bsbsteps) which contains information about the units
%               belonging to the subset at each step of the forward search.
%               BB has the following structure
%               1-st row has number 1 in correspondence of the steps in
%                   which unit 1 is included inside subset and a missing
%                   value for the other steps
%               ......
%               (n-1)-th row has number n-1 in correspondence of the steps in
%                   which unit n-1 is included inside subset and a missing
%                   value for the other steps
%               n-th row has number n in correspondence of the steps in
%                   which unit n is included inside subset and a missing
%                   value for the other steps
%               The size of matrix BB is n x (n-init+1) if option input
%               bsbsteps is 0 else the size is n-by-length(bsbsteps).
%  Bols:        OLS coefficents. Matrix.
%               (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward
%               search.
%  S2:          S2 and R2. Matrix.
%               (n-init+1) x 3 matrix containing the monitoring of S2 (2nd
%               column)and R2 (third column) in each step of the forward
%               search.
%
% See also: FSR, FSReda
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
% Atkinson, A.C. and Riani, M. (2006), Distribution theory and
% simulations for tests of outliers in regression, "Journal of
% Computational and Graphical Statistics", Vol. 15, pp. 460-476.
% Riani, M. and Atkinson, A.C. (2007), Fast calibrations of the forward
% search for testing multiple outliers in regression, "Advances in Data
% Analysis and Classification", Vol. 1, pp. 123-141.
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRmdr')">Link to the help function</a>
%
%$LastChangedDate:: 2018-06-21 15:29:09 #$: Date of the last commit

% Examples:

%{
    % FSRmdr with all default options.
    % Compute minimum deletion residual.
    % Monitor minimum deletion residual in each step of the forward search.
    % Common part to all examples: load fishery dataset.
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
     [mdr] = FSRmdr(y,X,out.bs);
     plot(mdr(:,1),mdr(:,2))
     title('Monitoring of minimum deletion residual')
%}

%{
    % FSRmdr with optional arguments.
    % Choose step to start monitoring.
    % Compute minimum deletion residual and start monitoring it from step
    % 60.
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr] = FSRmdr(y,X,out.bs,'init',60);
%}

%{
    % Analyze units entering the search in the final steps.
    % Compute minimum deletion residual and analyze the units entering
    % subset in each step of the fwd search (matrix Un).  As is well known,
    % the FS provides an ordering of the data from those most in agreement
    % with a suggested model (which enter the first steps) to those least in
    % agreement with it (which are included in the final steps).
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr,Un] = FSRmdr(y,X,out.bs);
%}





%% Input parameters checking

% Set up values for default model
modeldef          = struct;
modeldef.trend    = 1;
modeldef.s        = 12;       % monthly time series
modeldef.seasonal = [];
modeldef.X        = [];       % no explanatory variables
modeldef.posLS   = [];       % no level shift
modeldef.B        = [];        % empty initial parameter values
nocheck           = false;
StartDate         = '';

n=length(y);
bsbini=1:n;


%% User options

dispresultsdef=false;
options=struct('model',modeldef,'nocheck',0,'dispresults',dispresultsdef,...
    'StartDate',StartDate,'bsb',bsbini,'betaini',[],'plots',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:regressts:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin<1
    error('FSDA:regressts:missingInputs','Initial subset is missing');
end
%init1=options.init;
if nargin >1
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

% And check if the optional user parameters are reasonable.

% Default values for the optional parameters are set inside structure
% 'options'
if ~isequal(options.model,modeldef)
    fld=fieldnames(options.model);
    
    if nocheck == false
        % Check if user options inside options.model are valid options
        chkoptions(modeldef,fld)
    end
    for i=1:length(fld)
        modeldef.(fld{i})=options.model.(fld{i});
    end
    
end
dispresults=options.dispresults;
model = modeldef;
bsb=options.bsb;
plots=options.plots;

% Get model parameters
trend    = model.trend;       % get kind of  trend
s        = model.s;           % get periodicity of time series
seasonal = model.seasonal;    % get number of harmonics

B=model.B;


if isfield(model,'posLS')
    lshift   = model.posLS;
    posLS =lshift;
else
    lshift=0;
end

n=length(y);
T = n;

% seq is the vector which will contain linear time trend
seq   = (1:T)';
one   = ones(T,1);

% Construct the matrices which are fixed in each step of the minimization
% procedure
Seq = [one seq seq.^2 seq.^3];

% Define matrix which contains linear quadratic of cubic trend
Xtrend = Seq(:,1:trend+1);

ntrend = size(Xtrend,2);

% seasonal component
if seasonal >0
    sstring=num2str(seasonal);
    if seasonal>100
        varampl=str2double(sstring(1));
        seasonal=str2double(sstring(2:3));
    else
        varampl=0;
    end
    
    if seasonal < 1 || seasonal >floor(s/2)
        error('FSDA:FSRtsmdr:WrongInput',['Seasonal component must be an integer between 1 and ' num2str(floor(s/2))])
    end
    
    Xseaso=zeros(T,seasonal*2);
    for j=1:seasonal
        Xseaso(:,2*j-1:2*j)=[cos(j*2*pi*seq/s) sin(j*2*pi*seq/s)];
    end
    % Remark: when s is even the sine term disapperas for j=s/2 and so the
    % maximum number of trigonometric terms is s-1
    if seasonal==(s/2)
        Xseaso=Xseaso(:,1:end-1);
    end
    nseaso=size(Xseaso,2);
else
    Xseaso=[];
    nseaso=0;
    varampl=0;
end

X=model.X;
isemptyX=isempty(X);
if isemptyX
    % nexpl = number of potential explanatory variables
    nexpl=0;
else
    nexpl=size(X,2);
end

% Define the explanatory variable associated to the level shift component
if lshift>0
    % Xlshift = explanatory variable associated with
    % level shift Xlshift is 0 up to lsh-1 and 1 from
    % lsh to T
    Xlshift= [zeros(posLS-1,1);ones(T-posLS+1,1)];
else
    Xlshift =[];
end


% Construct matrix X (called Xsel) which contains the linear part of the model
if seasonal==0
    if isemptyX
        Xsel=Xtrend;
    else
        Xsel=[Xtrend X];
    end
else
    
    if isemptyX
        Xsel=[Xtrend Xseaso];
    else
        Xsel= [Xtrend Xseaso X];
    end
    % zero for varampl is automatically included because b0 is
    % initialized as a vector of zeroes b0=[b0;zeros(varampl,1)];
end

if lshift>0
    Xsel=[Xsel Xlshift];
end


% pini = number of parameters in the linear model without level shifts nor
% varying amplitude
% ntrend = number of trend parameters,
% nseaso = number of parameters associated with the harmonics,
% nexpl = number of explanatory variables,
% 1 parameter for fixed level shift position
pini=ntrend+nseaso+nexpl+(lshift>0);

% p = total number of parameters in the model
% nini +
% varampl = number of parameters involving time varying trend,
p=pini+varampl;

% indexes of linear part of seasonal component
if seasonal <6
    indlinsc=(trend+2):(trend+1+seasonal*2);
else
    indlinsc=(trend+2):(trend+1+seasonal*2-1);
end

otherind=setdiff(1:p,indlinsc);
% if lshift>0
%     otherind=otherind(1:end-1);
% end


% If model contains a field named B than use the first column of this field
% as initial parameter value, else use OLS estimate based on linear part of
% the model
if isfield('model',B) && ~isempty(model.B)
    b=model.B(:,1); % get initial estimate of parameter values
else
    
    % initial value of parameter estimates is based on subset
    bsel=Xsel(bsb,:)\y(bsb);
    if varampl>0
        if lshift>0
            b=[bsel(1:end-1); 0.01*zeros(varampl,1); bsel(end)];
            
        else
            b=[bsel; 0.01*zeros(varampl,1)];
        end
    end
end


% MaxIter=[];
MaxIter=1000;

DisplayLevel='';
nlinfitOptions=statset('Display',DisplayLevel,'MaxIter',MaxIter,'TolX',1e-7);

warning('off','stats:nlinfit:Overparameterized');
warning('off','stats:nlinfit:IterationLimitExceeded');
warning('off','stats:nlinfit:IllConditionedJacobian')
warning('off','MATLAB:rankDeficientMatrix')

%% Set up the model

% find estimate of beta and residuals
if varampl==0  % In this case the model is linear
    % Function lik constructs fitted values and residual sum of
    % squares
    betaout = Xsel(bsb,:) \ y(bsb);
    % update fitted values
    yhat = Xsel * betaout;
    
    s2=sum((y(bsb)-yhat(bsb)).^2)/(mm-size(Xsel,2));
    invXX=inv(Xsel'*Xsel);
    
else % model is non linear because there is time varying amplitude in seasonal component
    Xtrendf=Xtrend(bsb,:);
    Xseasof=Xseaso(bsb,:);
    if ~isempty(X)
        Xf=X(bsb,:);
    end
    Seqf=Seq(bsb,:);
    yf=y(bsb);
end

if lshift>0
    Xlshiftf=Xlshift(bsb);
end

Exflag=1;
exitflagALS=[];
iterA=[];
iterALS=0;
while iterALS < 2
    %   [b,exitflag,iter]=ALS(y,b,10000,1e-7);
    [betaout,~,Xsel,covB,s2,~]  = nlinfit(Xtrendf,yf,@likyhat,b,'options',nlinfitOptions);
    % Note that MSE*inv(J'*J) = covB
            invXX=covB/s2;

    
    [~,ID] = lastwarn;
    
    if iterALS == 0 && ~isempty(lastwarn) && ~strcmp(ID,'stats:nlinfit:IllConditionedJacobian')
        lastwarn('')
        % ID='';
        % [b,exitflag,iter]=ALS(y,b,10000,1e-7);
        [b,exitflagALS,iterA]=ALS(y,b,10000,1e-7);
        iterALS=iterALS+1;
    else
        iterALS=2;
    end
end

%                 [b,exitflag,iter]=ALS(y,b,10000,1e-7);
%                 [betaout,~,~,covB,s2,~]  = nlinfit(Xtrendf,yf,@likyhat,b,'options',nlinfitOptions);


% Note that MSE*inv(J'*J) = covB
[~,ID] = lastwarn;


if ~isempty(lastwarn) && ~strcmp(ID,'stats:nlinfit:IllConditionedJacobian')
    Exflag=0;
end



hrew=length(bsb);
T=n;
    if hrew<T
        % factor=consistencyfactor(hrew,n,1);
        a=norminv(0.5*(1+hrew/T));
        %factor=1/sqrt(1-(2*a.*normpdf(a))./(2*normcdf(a)-1));
        factor=1/sqrt(1-2*(T/hrew)*a.*normpdf(a));
        % Apply small sample correction factor to reweighted estimate
        % of sigma
        factor=factor*sqrt(corfactorREW(1,T,hrew/T));
    else
        factor=1;
    end
    
   
    bsb=1:n;
yhat=lik(betaout);


    
e=y-yhat;  % e = vector of residuals for all units using b estimated using subset
residuals =e./(factor*sqrt(s2));

out=struct;
out.s2=s2;
out.invXX=invXX;
out.e=e;
out.Exflag=Exflag;
out.Xsel=Xsel;
out.covB=covB;
out.exitflagALS=exitflagALS;
out.iterA=iterA;
out.beta=betaout;

if dispresults
    
    b_trend={'b_trend1'; 'b_trend2'; 'b_trend3'; 'b_trend4'};
    b_seaso={'b_cos1'; 'b_sin1'; 'b_cos2'; 'b_sin2'; 'b_cos3'; 'b_sin3'; ...
        'b_cos4'; 'b_sin4'; 'b_cos5'; 'b_sin5'; 'b_cos6'};
    b_expl={'b_X1'; 'b_X2'; 'b_X3'; 'b_X4'; 'b_X5'; 'b_X6'};
    b_varampl={'b_varampl'; 'b_varamp2'; 'b_varamp3'};
    b_lshift={'b_lshift'; 't_lshift'};
    
    if seasonal>0
        if 2*seasonal==s
            lab=[b_trend(1:trend+1); b_seaso];
        else
            lab=[b_trend(1:trend+1); b_seaso(1:2*seasonal)];
        end
    else
        lab=b_trend(1:trend+1);
    end
    
    if nexpl>0
        lab=[lab;b_expl(1:nexpl)];
    end
    if varampl>0
        lab=[lab;b_varampl(1:varampl)];
    end
    if lshift>0
        lab=[lab; b_lshift(1)];
    end
    
    bhat=outEST.B(:,1);
    se=outEST.B(:,2);
    tstat=outEST.B(:,3);
    pval=outEST.B(:,4);
    if verLessThan ('matlab','8.2.0')
        disp('           Coeff.     SE         t-stat       p-values');
        disp( [char(lab) num2str([bhat se tstat pval])]);
    else
        disp([table(lab) table(bhat) table(se) table(tstat) table(pval)]);
    end
    if lshift>0
        disp(['Level shift position t=' num2str(outEST.posLS)])
    end
end

if plots==1
    figure;
    % Plot original time series
    plot(datesnumeric(1:n),y,'k')
    hold('on')
    % plot fitted values
    plot(datesnumeric(1:n),yhat,'b-')
    
    
    title('Fit','interpreter','LaTex','FontSize',14)
    if ~isempty(StartDate)
        datetick('x','mmm-yy');
        if ~verLessThanFS(8.4)
            set(gca,'XTickLabelRotation',90);
        end
    end
    ax=axis;
    ylimits=ax(3:4);
    line(0.5*sum(datesnumeric(n:n+1))*ones(2,1),ylimits,'color','g')
   
    
    
    %% Create plots
% If plots is a structure, plot directly those chosen by the user; elseif
% plots is 1 a plot or residuals against index number appears else no plot
% is produced.
if plots>=1
    % Time series + fitted values
    figure
    subplot(2,1,1)
    plot([y yhat])
    xlabel('Time')
    ylabel('Real and fitted values')
    
    % Index plot of robust residuals
    h2=subplot(2,1,2);
    laby='Robust lts residuals';
    resindexplot(residuals,'conflev',conflev,'laby',laby,'h',h2,'title','');
end

end


% likyhat computes fitted values using vector of regression coefficients
% beta0. Note that matrices Xtrendf, Xseasof, Seqf, Xf contain n-k rows.
% This function is called in the very last step of the procedure when
% routine nlinfit is invoked. Please, note the difference beween likyhat
% and lik
    function objyhat=likyhat(beta0,Xtrendf)
        
        yhattrend=Xtrendf*beta0(1:trend+1);
        
        npar=trend+1;
        
        if seasonal >0
            if seasonal<s/2
                yhatseaso=Xseasof*beta0(npar+1:npar+seasonal*2);
                npar=npar+seasonal*2;
            else
                yhatseaso=Xseasof*beta0(npar+1:npar+seasonal*2-1);
                npar=npar+seasonal*2-1;
            end
            
            if varampl>0
                Xtre=1+Seqf(:,2:varampl+1)*beta0((npar+1+nexpl):(npar+varampl+nexpl));
                yhatseaso=Xtre.*yhatseaso;
                npar=npar+varampl;
            end
        end
        
        if isemptyX
            yhatX=0;
        else
            % Note the order of coefficients is trend, linear part of
            % seasonal component, expl variables, non linear part of
            % seasonal component, level shift
            yhatX=Xf(:,:)*beta0(npar+1-varampl:npar+nexpl-varampl);
            npar=npar+nexpl;
        end
        
        if lshift >0
            %  \beta_(npar+1)* I(t \geq \beta_(npar+2)) where beta_(npar+1)
            %  is a real number and \beta_(npar+2) is a integer which
            %  denotes the period in which level shift shows up
            
            yhatlshift=beta0(npar+1)*Xlshiftf;
        else
            yhatlshift=0;
        end
        
        % objhat = fitted values from trend (yhattrend), (time varying) seasonal
        % (yhatseaso), explanatory variables (yhatX) and level shift
        % component (yhatlshift)
        objyhat=yhattrend+yhatseaso+yhatX+yhatlshift;
    end

    function yhat=lik(beta0)
        
        yhattrend=Xtrend(bsb,:)*beta0(1:trend+1);
        npar=trend+1;
        
        if seasonal >0
            if seasonal<s/2
                yhatseaso=Xseaso(bsb,:)*beta0(npar+1:npar+seasonal*2);
                npar=npar+seasonal*2;
            else
                yhatseaso=Xseaso(bsb,:)*beta0(npar+1:npar+seasonal*2-1);
                npar=npar+seasonal*2-1;
            end
            
            if varampl>0
                Xtre=1+Seq(bsb,2:varampl+1)*beta0((npar+1+nexpl):(npar+varampl+nexpl));
                yhatseaso=Xtre.*yhatseaso;
                npar=npar+varampl;
            end
        end
        
        if isemptyX
            yhatX=0;
        else
            % Note the order of coefficients is trend, linear part of
            % seasonal component, expl variables, non linear part of
            % seasonal component, level shift
            yhatX=X(bsb,:)*beta0(npar+1-varampl:npar+nexpl-varampl);
            npar=npar+nexpl;
        end
        
        if lshift >0
            %  \beta_(npar+1)* I(t \geq \beta_(npar+2)) where beta_(npar+1)
            %  is a real number and \beta_(npar+2) is a integer which
            %  denotes the period in which level shift shows up
            yhatlshift=beta0(npar+1)*Xlshift(bsb);
        else
            yhatlshift=0;
        end
        
        % Fitted values from trend (yhattrend), (time varying) seasonal
        % (yhatseaso), explanatory variables (yhatX) and level shift
        % component (yhatlshift)
        yhat=yhattrend+yhatseaso+yhatX+yhatlshift;
        
    end

  


% ALS computes Alternating Least Squares estimate of beta starting from
% vector beta0. The rows which are used are those specified in global
% variable bsb
    function [newbeta,exitflag,iter]=ALS(y,beta0,maxiterALS,maxtolALS)
        
        
        iter        = 0;
        betadiff    = 9999;
        newbeta=beta0;
        oldbeta=beta0;
        % exitflag = flag which informs about convergence. exitflag =0
        % implies normal convergence, else no convergence has been obtained
        exitflag=0;
        
        while ( (betadiff > maxtolALS) && (iter < maxiterALS) )
            iter = iter + 1;
            
            % b2378 estimate of linear part of seasonal component
            b2378=newbeta(indlinsc);
            % at= yhatseaso = fitted values for linear part of seasonal
            % component
            at=Xseaso(bsb,:)*b2378;
            
            % OLS to estimate coefficients of trend + expl variables + non lin coeff of
            % seasonal + coefficient of fixed level shift
            % trlshift is the matrix of explanatory variables
            if isemptyX
                if lshift>0
                    tr_expl_nls_lshift=[Xtrend(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1)) Xlshift(bsb)];
                else
                    tr_expl_nls_lshift=[Xtrend(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1))];
                end
            else
                if lshift>0
                    tr_expl_nls_lshift=[Xtrend(bsb,:) X(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1)) Xlshift(bsb)];
                else
                    tr_expl_nls_lshift=[Xtrend(bsb,:) X(bsb,:) bsxfun(@times,at,Seq(bsb,2:varampl+1))];
                end
            end
            % b0145 = coefficients of intercept trend + expl var + non
            % linear part of seasonal component + level shift
            b0145=tr_expl_nls_lshift\(y(bsb)-at) ;
            
            % Now find new coefficients of linear part of seasonal
            % component in the regression of y-trend-expl-lsihft versus
            % vector which contains non linear part of seasonal component
            % which multiplies each column of matrix Xseaso (linear part of
            % seasonal component)
            yhatnlseaso=Seq(bsb,1)+ Seq(bsb,2:varampl+1)*b0145((trend+2+nexpl):(trend+2+nexpl+varampl-1));
            if isemptyX
                if lshift>0
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(y(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-Xlshift(bsb)*b0145(end));
                else
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(y(bsb)-Xtrend(bsb,:)*b0145(1:trend+1));
                end
            else
                if lshift>0
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(y(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-X(bsb,:)*b0145((trend+2):(trend+1+nexpl)) - Xlshift(bsb)*b0145(end));
                else
                    b2378=bsxfun(@times,yhatnlseaso,Xseaso(bsb,:))...
                        \(y(bsb)-Xtrend(bsb,:)*b0145(1:trend+1)-X(bsb,:)*b0145((trend+2):(trend+1+nexpl)));
                end
            end
            
            
            
            newbeta(indlinsc)=b2378;
            
            newbeta(otherind)=b0145;
            
            % betadiff is linked to the tolerance (specified in scalar
            % reftol)
            betadiff = norm(oldbeta - newbeta,1) / norm(newbeta,1);
            
            oldbeta=newbeta;
            
            % exit from the loop if the new beta has singular values. In
            % such a case, any intermediate estimate is not reliable and we
            % can just keep the initialbeta and initial scale.
            if (any(isnan(newbeta)))
                newbeta = beta0;
                exitflag=1;
                break
            end
        end
    end
end
%FScategory:REG-Regression