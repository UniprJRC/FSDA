function [Un,BB] = FSRtsbsb(y,bsb,varargin)
%FSRtsbsb returns the units belonging to the subset in each step of the forward search 
%
%<a href="matlab: docsearchFS('FSRtsbsb')">Link to the help function</a>
%
% Required input arguments:
%
%    y:         Time series to analyze. Vector. A row or a column vector
%               with T elements, which contains the time series.
%  bsb :        list of units forming the initial subset. Vector | 0. If
%               bsb=0 then the procedure starts with p units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb)
%
% Optional input arguments:
%
%  init :       Start of monitoring point. Scalar.
%               It specifies the point where to initialize the search and
%               start monitoring required diagnostics. If it is not
%               specified it is set equal floor(0.5*(T+1))
%               Example - 'init',100 starts monitoring from step m=100
%               Data Types - double
%
%      model :  model type. Structure. A structure which specifies the model
%               which will be used. The model structure contains the following
%               fields:
%               model.s = scalar (length of seasonal period). For monthly
%                         data s=12 (default), for quartely data s=4, ...
%               model.trend = scalar (order of the trend component).
%                       trend = 1 implies linear trend with intercept (default),
%                       trend = 2 implies quadratic trend ...
%                       Admissible values for trend are, 0, 1, 2 and 3.
%               model.seasonal = scalar (integer specifying number of
%                        frequencies, i.e. harmonics, in the seasonal
%                        component. Possible values for seasonal are
%                        $1, 2, ..., [s/2]$, where $[s/2]=floor(s/2)$.
%                        For example:
%                        if seasonal =1 (default) we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 sin ( 2 \pi t/s)$;
%                        if seasonal =2 we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 \sin ( 2 \pi t/s)
%                        + \beta_3 \cos(4 \pi t/s) + \beta_4 \sin (4 \pi t/s)$.
%                        Note that when $s$ is even the sine term disappears
%                        for $j=s/2$ and so the maximum number of
%                        trigonometric parameters is $s-1$.
%                        If seasonal is a number greater than 100 then it
%                        is possible to specify how the seasonal component
%                        grows over time.
%                        For example, seasonal =101 implies a seasonal
%                        component which just uses one frequency
%                        which grows linearly over time as follows:
%                        $(1+\beta_3 t)\times ( \beta_1 cos( 2 \pi t/s) +
%                        \beta_2 \sin ( 2 \pi t/s))$.
%                        For example, seasonal =201 implies a seasonal
%                        component which just uses one frequency
%                        which grows in a quadratic way over time as
%                        follows:
%                        $(1+\beta_3 t + \beta_4  t^2)\times( \beta_1 \cos(
%                        2 \pi t/s) + \beta_2 \sin ( 2 \pi t/s))$.
%                        seasonal =0 implies a non seasonal model.
%               model.X  =  matrix of size T-by-nexpl containing the
%                         values of nexpl extra covariates which are likely
%                         to affect y.
%               model.posLS = positive integer which specifies to position
%                         to include the level shift component.
%                         For example if model.posLS =13 then the
%                         explanatory variable $I(t \geq 13})$ is created.
%                         If this field is not present or if it is empty,
%                         the level shift component is not included.
%               model.B  = column vector or matrix containing the initial
%                         values of parameter estimates which have to be used in the
%                         maximization procedure. If model.B is a matrix,
%                         then initial estimates are extracted from the
%                         first colum of this matrix. If this field is
%                         empty or if this field is not present, the
%                         initial values to be used in the maximization
%                         procedure are referred to the OLS parameter
%                         estimates of the linear part of the model. The
%                         parameters associated to time varying amplitude
%                         are initially set to 0.
%                 Example - 'model', model
%                 Data Types - struct
%               Remark: the default model is for monthly data with a linear
%               trend (2 parameters) + seasonal component with just one
%               harmonic (2 parameters), no additional explanatory
%               variables and no level shift that is
%                               model=struct;
%                               model.s=12;
%                               model.trend=1;
%                               model.seasonal=1;
%                               model.X='';
%                               model.posLS='';
%
%  nocheck:     Check input arguments. Boolean.
%               If nocheck is equal to true no check is performed on
%               supplied structure model
%               Example - 'nocheck',false
%               Data Types - logical
%
%   bsbsteps :  Save the units forming subsets in selected steps. Vector.
%               It specifies for which steps of the fwd search it is
%               necessary to save the units forming subset. If bsbsteps is
%               0 we store the units forming subset in all steps. The
%               default is store the units forming subset in all steps if
%               n<=5000, else to store the units forming subset at steps
%               init and steps which are multiple of 100. For example, as
%               default, if n=7530 and init=6, units forming subset are
%               stored for
%               m=init, 100, 200, ..., 7500.
%               Example - 'bsbsteps',[100 200] stores the unis forming
%               subset in steps 100 and 200.
%               Data Types - double
%
%       plots   : Plot on the screen. Scalar.
%                 If plots=1 the monitoring units plot is displayed on the
%                 screen. The default value of plots is 0 (that is no plot
%                 is produced on the screen).
%                 Example - 'plots',1
%                 Data Types - double
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
%  BB:          Units belonging to subset in each step or selected steps. Matrix.
%               n-by-(n-init+1) or n-by-length(bsbsteps) matrix which
%               contains the units belonging to the subset at each step (or
%               in selected steps as specified by optional vector bsbsteps)
%               of the forward search.
%               More precisely:
%               BB(:,1) contains the units forming subset in step bsbsteps(1);
%               ....;
%               BB(:,end) contains the units forming subset in step  bsbsteps(end);
%               Row 1 of matrix BB is referred to unit 1;
%               ......;
%               Row n of matrix BB is referred to unit n;
%               Units not belonging to subset are denoted with NaN.
%
% See also FSRbsb, FSRBbsb, FSRHbsb
%
% See also: FSRts, LTSts, regressts
%
% References:
%
% Atkinson, A.C. and Riani, M. (2006), Distribution theory and
% simulations for tests of outliers in regression, "Journal of
% Computational and Graphical Statistics", Vol. 15, pp. 460-476.
% Riani, M. and Atkinson, A.C. (2007), Fast calibrations of the forward
% search for testing multiple outliers in regression, "Advances in Data
% Analysis and Classification", Vol. 1, pp. 123-141.
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [RPRH]
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRtsbsb')">Link to the help function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit

% Examples:

%{
    % FSRtsbsb with all default options.
    load('fishery');
    y=fishery.data(:,1);
    bsbini=[97    77    12     2    26    95    10    60    94   135     7    61   114];
    [Un,BB]=FSRtsbsb(y,bsbini);
%}

%{
    %% FSRtsbsb with optional arguments.
    % Load airline data
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960.
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl
    y=(y(:));
    % Define the model  and show the monitoring units plots.
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=104;         % four harmonics with time varying seasonality
    bsbini=[97    77    12     2    26    95    10    60    94   135     7    61   114];
    [Un,BB]=FSRtsbsb(y,bsbini,'model',model,'plots',1);
%}

%{
    %% Monitoring the units belonging to subset.
    % Load airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
         118  126  150  180  196  188  233  277  301  318  342  391    % Feb
         132  141  178  193  236  235  267  317  356  362  406  419    % Mar
         129  135  163  181  235  227  269  313  348  348  396  461    % Apr
         121  125  172  183  229  234  270  318  355  363  420  472    % May
         135  149  178  218  243  264  315  374  422  435  472  535    % Jun
         148  170  199  230  264  302  364  413  465  491  548  622    % Jul
         148  170  199  242  272  293  347  405  467  505  559  606    % Aug
         136  158  184  209  237  259  312  355  404  404  463  508    % Sep
         119  133  162  191  211  229  274  306  347  359  407  461    % Oct
         104  114  146  172  180  203  237  271  305  310  362  390    % Nov
         118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl
    y=(y(:));
    % Contaminates units 31:40
    y(31:40)=y(31:40)+200;
    % Define the model  and show the monitoring units plots.
    model=struct;
    model.trend=1;              % linear trend
    model.s=12;                 % monthly time series
    model.seasonal=104;         % four harmonics with time varying seasonality
    bsbini=[97    77    12     2    26    95    10    60    94   135     7    61   114];
    [Un,BB]=FSRtsbsb(y,0,'model',model,'plots',1);
    % Create the 'monitoring units plot'
    figure;
    seqr=[Un(1,1)-1; Un(:,1)];
    plot(seqr,BB','bx');
    xlabel('Subset size m');
    ylabel('Monitoring units plot');
    % The plot, which monitors the units belonging to subset in each step of
    % the forward search shows that independently of the initial starting
    % point the contaminated units (31:40) are always the last to enter the
    % forward search.
%}

%% Beginning of code 

% Input parameters checking

% Set up values for default model
modeldef          = struct;
modeldef.trend    = 1;
modeldef.s        = 12;       % monthly time series
modeldef.seasonal = [];
modeldef.X        = [];       % no explanatory variables
modeldef.posLS   = [];       % no level shift
modeldef.B        = [];        % empty initial parameter values


% User options
n=length(y);
init=floor(0.5*(n+1));

if init<length(bsb)
    init=length(bsb);
end

bsbstepdef='';

options=struct('init',init,'nocheck',0,'plots',0,...
    'bsbsteps',bsbstepdef,'model',modeldef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRtsbsb:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin<2
    error('FSDA:FSRtsbsb:missingInputs','Initial subset is missing');
end

if nargin >2
    % We now overwrite inside structure options the default values with
    % those chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end
nocheck=options.nocheck;
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

model = modeldef;

% Get model parameters
trend    = model.trend;       % get kind of  trend
s        = model.s;           % get periodicity of time series
seasonal = model.seasonal;    % get number of harmonics

if isfield(model,'posLS') && ~isempty(model.posLS)
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

if bsb==0
    Ra=1; nwhile=1;
    sizRandomSubsets=max([p+1 round(n/4)]);
    while and(Ra,nwhile<100)
        bsb=randsample(n,sizRandomSubsets);
        % bsbini=bsb;
        Xb=Xsel(bsb,:);
        Ra=(rank(Xb)<size(Xb,2));
        nwhile=nwhile+1;
    end
    if nwhile==100
        warning('FSDA:FSRtsmdr:NoFullRank','Unable to randomly sample full rank matrix');
    end
else
end

ini0=length(bsb);

% check init
init=options.init;
if  init <p+1
    fprintf(['Attention : init1 should be larger than p. \n',...
        'It is set to p+1.']);
    init=p+1;
elseif init<ini0
    fprintf(['Attention : init1 should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    init=ini0;
elseif init>=n
    fprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n-1.']);
    init=n-1;
end

% If model contains a field named B than use the first column of this field
% as initial parameter value, else use OLS estimate based on linear part of
% the model
if ~isempty(model.B)
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
% posvarampl = position of non linear term of time varying seasonal
% component inside b
posvarampl=p-varampl+1:p;
posvarampl=posvarampl-(lshift>0);

if varampl>0
    bprevious=b;
end

%% Initialise key matrices

% sequence from 1 to n.
seq = (1:n)';

% The second column of matrix R will contain the OLS residuals at each step
% of the forward search
r = [seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100 = 100*(1:1:ceil(n/100));

bsbsteps=options.bsbsteps;
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

%  UN is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init+1:n)' , NaN(n-init,10));

% The last correctly computed beta oefficients
blast=NaN(p,1);

Xb=Xsel(bsb,:);

% MaxIter=[];
MaxIter=1000;

DisplayLevel='';
nlinfitOptions=statset('Display',DisplayLevel,'MaxIter',MaxIter,'TolX',1e-7);

warning('off','stats:nlinfit:Overparameterized');
warning('off','stats:nlinfit:IterationLimitExceeded');
warning('off','stats:nlinfit:IllConditionedJacobian')
warning('off','MATLAB:rankDeficientMatrix')
%% Start of the forward search
if rank(zscore(Xb(:,2:end)))<pini-1
    warning('FSDA:FSRtsmdr:NoFullRank','Supplied initial subset does not produce full rank matrix');
    warning('FSDA:FSRtsmdr:NoFullRank','FS loop will not be performed');
    % FS loop will not be performed
else
    % ij = index which is linked with the columns of matrix BB. During the
    % search every time a subset is stored inside matrix BB ij icreases by one
    ij=1;
    
    for mm = ini0:n
         oldbsb=bsb;
         
        % if n>200 show every 100 steps the fwd search index
        if n>200
            if length(intersect(mm,seq100))==1
                disp(['m=' int2str(mm)]);
            end
        end
        
        % Store units belonging to the subset
        if (mm>=init)
            if intersect(mm,bsbsteps)==mm
                BB(oldbsb,ij)=oldbsb;
                ij=ij+1;
            end
        end
        
        % Compute beta coefficients using subset
        
          % Note that Xsel is the X matrix of the linearized version if the
        % model is non linear (that is it contains time varying amplitude)
        NoRankProblem=( rank(zscore(Xsel(bsb,2:end))) == size(Xsel,2)-1 );
        
        if NoRankProblem  % rank is ok
            
            % find estimate of beta and residuals
            if varampl==0  % In this case the model is linear
                % Function lik constructs fitted values and residual sum of
                % squares
                betaout = Xsel(bsb,:) \ y(bsb);
                % update fitted values
                yhat = Xsel * betaout;
                
                %s2=sum((y(bsb)-yhat(bsb)).^2)/(mm-size(Xsel,2));
                %invXX=inv(Xsel'*Xsel);
                
            else % model is non linear because there is time varying amplitude in seasonal component
                Xtrendf=Xtrend(bsb,:);
                Xseasof=Xseaso(bsb,:);
                if ~isempty(X)
                    Xf=X(bsb,:);
                end
                Seqf=Seq(bsb,:);
                yf=y(bsb);
                
                if lshift>0
                    Xlshiftf=Xlshift(bsb);
                end
                
                iterALS=0;
                while iterALS < 2
                    %   [b,exitflag,iter]=ALS(y,b,10000,1e-7);
                    [betaout,~,~,~,~,~]  = nlinfit(Xtrendf,yf,@likyhat,b,'options',nlinfitOptions);
                    % Note that MSE*inv(J'*J) = covB
                    [~,ID] = lastwarn;
                    
                    if iterALS == 0 && ~isempty(lastwarn) && ~strcmp(ID,'stats:nlinfit:IllConditionedJacobian')
                        lastwarn('')
                        % ID='';
                        % [b,exitflag,iter]=ALS(y,b,10000,1e-7);
                        [b]=ALS(y,b,10000,1e-7);
                        iterALS=iterALS+1;
                    else
                        iterALS=2;
                    end
                end
                
                % Note that MSE*inv(J'*J) = covB
                [~,ID] = lastwarn;
                
                if ~isempty(lastwarn) && strcmp(ID,'stats:nlinfit:ModelConstantWRTParam')
                    betaout=bprevious;
                end
                
                % clear "last warning"
                lastwarn('')
                
                
                % Now compute vector yhat for all the observations
                % using input  vector betaout
                bsb=seq;
                yhat=lik(betaout);
                
                % Xsel always has n rows and p columns and is referred to
                % all the units (this is the linearized version of matrix
                % X)
                fdiffstep=[];
                Xsel = getjacobianFS(betaout,fdiffstep,@lik,yhat);
            end
  
            % Check whether the estimate of b which has come out is
            % reasonable. An estimate of b is called unreasonable if
            % max(yhat)>2*max(y)  and min(yhat)<0.5*min(y)
            % Make sure that the coefficients of posvarampl are set to 0 if
            % they are greater than a certain threshold
            if max(abs(betaout(posvarampl)))>10
                betaout(posvarampl)=0;
            end
            
            if max(yhat)>2*max(y) && min(yhat)<0.5*min(y)
                b=bprevious;
            else
                b=betaout;
            end
            
            % Store correctly computed b for the case of rank problem
            bprevious=b;
        else   % number of independent columns is smaller than number of parameters
            if bsbmfullrank
                Xb=Xsel(bsb,:);
                Xbx=Xb;
                nclx=ncl;
                bsbx=zeros(n,1);
                bsbx(1:mm)=bsb;
                norank=1;
                while norank ==1
                    
                    norank=1;
                    % Increase the size of the subset by one unit iteratively until you
                    % obtain a full rank matrix
                    for i=1:length(nclx)
                        Xbb=[Xbx;Xsel(nclx(i),:)];
                        if length(nclx)==1
                            norank=0;
                            break
                        end
                        if rank(zscore(Xbb(:,2:end))) == size(Xbb,2)-1
                            norank=0;
                            break
                        else
                            bsbx(1:size(Xbb,1))=[bsbx(1:size(Xbb,1)-1);nclx(i)];
                            Xbx=Xsel(bsbx(1:size(Xbb,1)),:);
                            nclx=setdiff(seq,bsbx(1:size(Xbb,1)));
                            norank=1;
                            break
                        end
                    end
                end
                % check how many observations produce a singular X matrix
                bsbsing=bsbx(1:size(Xbb,1)-1);
                
                if msg==1
                    warning('FSDA:FSRtsmdr','Rank problem in step %d:',mm);
                    disp('Observations')
                    disp(bsbsing')
                    disp('produce a singular matrix')
                end
                Un=NaN;
                BB=NaN;
                return
                
            else
                disp(['Matrix without full rank at step m=' num2str(mm)])
                disp('Estimate of \beta which is used is based on previous step with full rank')
                b=blast;
                % disp([mm b'])
            end
        end
        
        
           e=y-yhat;  % e = vector of residuals for all units using b estimated using subset
        r(:,2)=e.^2;

        
        if mm<n
            
            % order the r_i and include the smallest among the units forming
            % the group of potential outliers
            ord=sortrows(r,2);
            
            % bsb= units forming the new subset
            bsb=ord(1:(mm+1),1);
            
            % yb=y(bsb);    % subset of y
            
            if mm>=init
                unit=setdiff(bsb,oldbsb);
                % If the interchange involves more than 10 units, store only the
                % first 10.
                if (size(unit,2)<=10)
                    Un(mm-init+1,2:(size(unit,1)+1)) = unit;
                else
                    disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                    Un(mm-init+1,2:end) = unit(1:10)';
                end
            end
        end
    end
end  % no rank
plots=options.plots;
if plots==1
    % Create the 'monitoring units plot'
    figure;
    plot(bsbsteps,BB','bx')
    xlabel('Subset size m');
    ylabel('Monitoring units plot');
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

    function J = getjacobianFS(beta,fdiffstep,modelFS,yfit)
        function yplus = call_model_nested(betaNew)
            yplus = modelFS(betaNew);
        end
        J = statjacobianFS(@call_model_nested, beta, fdiffstep, yfit);
    end % function getjacobian

    function J = statjacobianFS(func, theta, DerivStep, y0)
        %STATJACOBIAN Estimate the Jacobian of a function
        
        % J is a matrix with one row per observation and one column per model
        % parameter. J(i,j) is an estimate of the derivative of the i'th
        % observation with respect to the j'th parameter.
        
        % For performance reasons, very little error checking is done on the input
        % arguments. This function makes the following assumptions about inputs:
        %
        % * func is the model function and is a valid function handle that accepts
        %   a single input argument of the same size as theta.
        % * theta is vector or matrix of parameter values. If a matrix, each row
        %   represents a different group or observation (see "Grouping Note" below)
        %   and each column represents a different model parameter.
        % * DerivStep (optional) controls the finite differencing step size. It may
        %   be empty, scalar, or a vector of positive numbers with the number of
        %   elements equal to the number model parameters.
        % * y0 (optional) is the model function evaluated at theta. A value of []
        %   is equivalent to omitting the argument and results in the model being
        %   evaluated one additional time.
        %
        % Example 1: NLINFIT
        %   NLINFIT is used to estimate the parameters b(1) and b(2) for the model
        %   @(b,T) b(1)*sin(b(2)*T), given data at T=1:5. NLINFIT needs the
        %   Jacobian of the model function with respect to b(1) and b(2) at each T.
        %   To do this, it constructs a new function handle that is only a function
        %   of b and that "burns-in" the value of T (e.g. model2 = @(b) model1(b,T)).
        %   It then calls STATJACOBIAN with the new function handle to obtain a
        %   matrix J, where J(i,j) is an estimate of the derivative of the model
        %   with respect to the j'th parameter evaluated at T(i) and b.
        %
        % Example 2: NLMEFIT or NLMEFITSA with group-specific parameters
        %   NLMEFIT requires the Jacobian of the model function with respect to two
        %   parameters evaluated at group-specific values. (Group-specific
        %   parameters can arise, for example, from using the default FEConstDesign
        %   and REConstDesign options.) NLMEFIT calls STATJACOBIAN passing in a
        %   matrix of parameter values theta, with one row per group, where
        %   theta(i,j) represents a parameter value for i'th group and j'th
        %   parameter. STATJACOBIAN returns a matrix J, where J(i,j) is an estimate
        %   of the derivative of the model with respect to the j'th parameter,
        %   evaluated for observation i with parameter values theta(rowIdx(i),:),
        %   which are the parameter values for the observation's group.
        %
        % Example 3: NLMEFIT with observation-specific parameters
        %   NLMEFIT requires the Jacobian of the model function with respect to two
        %   parameters evaluated at observation-specific values. (Observation-
        %   specific parameters can arise, for example, from using the FEObsDesign
        %   or REObsDesign options.) NLMEFIT calls STATJACOBIAN passing in a matrix
        %   of parameter values theta, with one row per observation, where
        %   theta(i,j) represents a parameter value for the i'th observation and
        %   j'th parameter. In this case, rowIdx is 1:N, where N is the number of
        %   observations. STATJACOBIAN returns a matrix J, where J(i,j) is an
        %   estimate of the derivative of the model with respect to the j'th
        %   parameter, evaluated for observation i with parameter values
        %   theta(i,:), which are the parameter values for the observation.
        
        % Use the appropriate class for variables.
        classname = class(theta);
        
        % Handle optional arguments, starting with y0 since it will be needed to
        % determine the appropriate size for a default groups.
        if nargin < 4 || isempty(y0)
            y0 = func(theta);
        end
        
        % When there is only one group, ensure that theta is a row vector so
        % that vectoriation works properly. Also ensure that the underlying
        % function is called with an input with the original size of theta.
        thetaOriginalSize = size(theta);
        theta = reshape(theta, 1, []);
        
        func = @(theta) func(reshape(theta, thetaOriginalSize));
        
        % All observations belong to a single group; scalar expansion allows us
        % to vectorize using a scalar index.
        rowIdx = 1;
        
        [numThetaRows, numParams] = size(theta);
        
        if nargin < 3 || isempty(DerivStep)
            % Best practice for forward/backward differences:
            DerivStep = repmat(sqrt(eps(classname)), 1, numParams);
            % However, NLINFIT's default is eps^(1/3).
        elseif isscalar(DerivStep)
            DerivStep = repmat(DerivStep, 1, numParams);
        end
        
        delta = zeros(numThetaRows, numParams, classname);
        J = zeros(numel(y0), numParams, classname);
        for ii = 1:numParams
            % Calculate delta(:,ii), but remember to set it back to 0 at the end of the loop.
            delta(:,ii) = DerivStep(ii) * theta(:,ii);
            deltaZero = delta(:,ii) == 0;
            if any(deltaZero)
                % Use the norm as the "scale", or 1 if the norm is 0.
                nTheta = sqrt(sum(theta(deltaZero,:).^2, 2));
                delta(deltaZero,ii) = DerivStep(ii) * (nTheta + (nTheta==0));
            end
            thetaNew = theta + delta;
            yplus = func(thetaNew);
            dy = yplus(:) - y0(:);
            J(:,ii) = dy./delta(rowIdx,ii);
            delta(:,ii) = 0;
        end
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
                exitflag=-1;
                break
            end
        end
    end
end
%FScategory:REG-Regression