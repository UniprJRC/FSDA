function [out , varargout] = LXS(y,X,varargin)
%LXS computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators
%
%<a href="matlab: docsearchFS('LXS')">Link to the help function</a>
%
%  Required input arguments:
%
%    y:         Response variable. Vector. A vector with n elements that
%               contains the response
%               variable.  It can be either a row or a column vector.
%    X :        Predictor variables. Matrix. Data matrix of explanatory
%               variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables.
%
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%
%  Optional input arguments:
%
%   intercept :  Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               else no constant term will be included.
%               Example - 'intercept',1
%               Data Types - double
%           h : The number of observations that have determined the least
%                 trimmed squares estimator. Scalar.
%               The number of observations that have determined the least
%               trimmed squares estimator. h is an integer greater than p
%               (number of columns of matrix X including the intercept but
%               smaller then n. If the purpose is outlier detection than h
%               does not have to be smaller than [0.5*(n+p+1)] (default
%               value). On the other hand if the purpose is to find
%               subgroups of homogeneous observations h can be smaller than
%               [0.5*(n+p+1)]. If h <p+1 an error will be given.
%                 Example - 'h',round(n*0,75)
%                 Data Types - double
%         bdp :  breakdown point. Scalar.
%               It measures the fraction of outliers
%               the algorithm should
%               resist. In this case any value greater than 0 but smaller
%               or equal than 0.5 will do fine. If on the other hand the
%               purpose is subgroups detection then bdp can be greater than
%               0.5. In any case however n*(1-bdp) must be greater than
%               p. If this condition is not fulfilled an error will be
%               given. Please specify h or bdp not both.
%                 Example - 'bdp',0.4
%                 Data Types - double
%       nsamp : Number of subsamples which will be extracted to find the
%               robust estimator. Scalar.
%               If nsamp=0 all subsets will be extracted. They will be (n choose p).
%                 Example - 'nsamp',0
%                 Data Types - double
%               Remark: if the number of all possible subset is <1000 the
%               default is to extract all subsets, otherwise just 1000 if
%               fastLTS is used (lms=2 or lms is a structure) or 3000 for
%               standard LTS or LMS.
%       lms   : Criterion to use to find the initlal
%                 subset to initialize the search. Scalar, vector or structure.
%               If lms is a scalar = 1 (default) Least Median of Squares is
%                       computed,
%               else if lms is a scalar = 2 fast lts with the all default options is used
%               else if lms is a scalar different from 1 and 2 standard lts
%                       is used (without concentration steps)
%               else if lms is a struct fast lts (with concentration steps) is used.
%                  In this case the user can control the following options:
%                  lms.refsteps : scalar defining number of refining iterations in each
%                               subsample (default = 3). refsteps = 0 means
%                               "raw-subsampling" without iterations.
%                   lms.reftol  : scalar. Default value of tolerance for the refining steps
%                               The default value is 1e-6.
%                   lms.bestr   : scalar defining number of "best betas" to remember from the
%                               subsamples. These will be later iterated until convergence
%                               (default=5).
%             lms.refstepsbestr : scalar defining number of refining iterations for each
%                               best subset (default = 50).
%              lms.reftolbestr  : scalar. Default value of tolerance for the refining steps
%                               for each of the best subsets
%                               The default value is 1e-8.
%                 Example - 'lms',1
%                 Data Types - double
%       rew   : LXS reweighted. Scalar.
%                If rew=1 the reweighted version of LTS (LMS) is
%               used and the output quantities refer to the reweighted
%               version
%               else no reweighting is performed (default).
%                 Example - 'rew',1
%                 Data Types - double
%     conflev :  Confidence level which is
%               used to declare units as outliers. Scalar
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%                 Example - 'conflev',0.99
%                 Data Types - double
%       plots : Plot on the screen. Scalar or structure.
%               If plots = 1, a plot which shows the
%               robust residuals against index number is shown on the
%               screen. The confidence level which is used to draw the
%               horizontal lines associated with the bands for the
%               residuals is as specified in input option conflev. If
%               conflev is missing a nominal 0.975 confidence interval will
%               be used.
%                 Example - 'plots',1
%                 Data Types - double
%        msg  : It controls whether to display or not messages on the screen. Scalar.
%                If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the estimator
%               and the warnings about
%               'MATLAB:rankDeficientMatrix', 'MATLAB:singularMatrix' and
%               'MATLAB:nearlySingularMatrix' are set to off
%               else no message is displayed on the screen
%               Example - 'msg',1
%               Data Types - double
%      nocheck: Check input arguments. Scalar. If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additioanl column of ones for
%               the intercept is not added. As default nocheck=0. The
%               controls on h, bdp and nsamp still remain.
%               Example - 'nocheck',1
%               Data Types - double
%        nomes:  It controls whether to display or not on the screen
%               messages about estimated  time to compute LMS (LTS) . Scalar.
%               If nomes is equal to 1 no message about estimated
%               time to compute LMS (LTS) is displayed, else if nomes is
%               equal to 0 (default), a message about estimated time is
%               displayed.
%               Example - 'nomes',1
%               Data Types - double
%       yxsave : the response vector y and data matrix X are saved into the output
%                structure out. Scalar.
%               Default is 0, i.e. no saving is done.
%               Example - 'yxsave',1
%               Data Types - double
%
%
%       Remark: The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%
%
%  Output:
%
%  out :     A structure containing the following fields
%
%            out.rew  = Scalar if out.rew=1 all subsequent output refers to
%                       reweighted else no reweighting is done.
%            out.beta = Vector of beta LTS (LMS) coefficient estimates,
%                       including the intercept when options.intercept=1.
%                       out.beta=[intercept slopes].
%              out.bs = p x 1 vector containing the units forming subset
%                       associated with bLMS (bLTS).
%       out.residuals = Vector containing the standardized residuals from
%                       the regression.
%           out.scale = Scale estimate of the residuals.
%         out.weights = Vector like y containing weights. The elements of
%                       this vector are 0 or 1.
%                       These weights identify the h observations which are
%                       used to compute the final LTS (LMS) estimate.
%                       sum(out.weights)=h if there is not a perfect fit
%                       otherwise sum(out.weights) can be greater than h
%               out.h = The number of observations that have determined the
%                       LTS (LMS) estimator, i.e. the value of h.
%        out.outliers = vector containing the list of the units declared
%                       as outliers using confidence level specified in
%                       input scalar conflev
%         out.conflev = confidence level which is used to declare outliers.
%                       Remark: scalar out.conflev will be used
%                       to draw the horizontal lines (confidence bands) in the plots
%         out.singsub = Number of subsets wihtout full rank. Notice that if
%                       this number is greater than 0.1*(number of
%                       subsamples) a warning is produced on the screen
%           out.class = 'LTS' or 'LMS'.
%            out.y    = response vector Y. The field is present if option
%                       yxsave is set to 1.
%            out.X    = data matrix X. The field is present if option
%                       yxsave is set to 1.
%
%  Optional Output:
%
%            C        : Indexes of the extracted subsamples. Matrix.
%                       Matrix containing the indexes of the subsamples
%                       extracted for computing the estimate (the so called
%                       elemental sets). For example, if C(3,:)=[2 5 20],
%                       implies that the third extracted subsample is
%                       formed by units 2, 5 and 20.
%
%
% See also FSReda, Sreg, MMreg
%
% References:
%
%   Rousseeuw PJ, Leroy AM (1987), Robust regression and outlier detection,
%   Wiley.
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('LXS')">Link to the help function</a>
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%
%{
    % LXS with default input and output.
    % Compute LMS estimator without reweighting, add intercept to matrix X
    % and do not produce plots.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)+6;
    [out]=LXS(y,X);
%}

%{
    % LXS with optional arguments.
    % Compute LMS estimator, reweight and plot the residuals.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)+6;
    [out]=LXS(y,X,'rew',1,'plots',1);
%}

%{
    % LXS with optional output.
    % Generating the C matrix containing the indices of the subsamples
    % extracted for computing the estimate (the so called elemental sets).
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)+6;
    [out,C]=LXS(y,X);
%}

%{
    % Reweighted LTS estimator.
    % Compute reweighted LTS estimator and produce the plot of
    % residuals.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)+6;
    [out]=LXS(y,X,'rew',1,'lms',0,'plots',1);
%}

%{
    % Specifying the number of subsamples.
    % Compute LMS estimator, without plots using 20000 subsamples.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)+6;
    [out]=LXS(y,X,'nsamp',20000);
%}

%{
    % Specifying a confidence level.
    % Compute reweighted LMS and use a confidence level for outlier
    % detection equal to 0.999.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)+6;
    [out]=LXS(y,X,'rew',1,'conflev',0.999);
%}

%{
    % Using fast options.
    % Compute LTS using fast options
    % detection equal to 0.999.
    lms=struct;
    % Do 5 refining steps for each elemental subset
    lms.refsteps=5;
    % Store the best 10 subsets
    lms.bestr=10;
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)+6;
    [out]=LXS(y,X,'lms',lms,'plots',1);
%}

%{
    %% We compare the output of different procedures:  only the reweighted
    % LTS seems to detect half of the outlier with a Bonferroni
    %significance level.
    close all;
    state=100;
    randn('state', state);
    n=100;
    X=randn(n,3);
    bet=[3;4;5];
    y=3*randn(n,1)+X*bet;
    y(1:20)=y(1:20)+13;

    % Define nominal confidence level
    conflev=[0.99,1-0.01/length(y)];
    % Define number of subsets
    nsamp=3000;
    % Define the main title of the plots
    titl='';

    % LMS with no reweighting
    [outLMS]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1));
    h=subplot(2,2,1)
    laby='Scaled LMS residuals';
    resindexplot(outLMS.residuals,'h',h,'title',titl,'laby',laby,'numlab','','conflev',conflev)

    % LTS with no reweighting
    h=subplot(2,2,2);
    [outLTS]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'lms',0);
    laby='Scaled LTS residuals';
    resindexplot(outLTS.residuals,'h',h,'title',titl,'laby',laby,'numlab','','conflev',conflev);

    % LMS with reweighting
    [outLMSr]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'rew',1);
    h=subplot(2,2,3);
    laby='Scaled reweighted LMS residuals';
    resindexplot(outLMSr.residuals,'h',h,'title',titl,'laby',laby,'numlab','','conflev',conflev)

    % LTS with reweighting
    [outLTSr]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'rew',1,'lms',0);
    h=subplot(2,2,4);
    laby='Scaled reweighted LTS residuals';
    resindexplot(outLTSr.residuals,'h',h,'title',titl,'laby',laby,'numlab','','conflev',conflev);
    % By simply changing the seed to 543 (state=543), using a Bonferroni size of 1%, no unit is declared as outlier.
%}

%% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);
seq=(1:n)';

%% User options

% singsub= scalar which will contain the number of singular subsets which
% are extracted (that is the subsets of size p which are not full rank)
singsub=0;

% If the number of all possible subsets is <10000 the default is to extract
% all subsets otherwise just 10000.
% Notice that we use bc, a fast version of nchoosek. One may also use the
% approximation floor(exp(gammaln(n+1)-gammaln(n-p+1)-gammaln(p+1))+0.5)
ncomb=bc(n,p);
nsampdef=min(3000,ncomb);

% Set the "half" of the data points.
hdef=floor(0.5*(n+p+1));
bdpdef=1-hdef/n;

hmin=p+1;

% initialize brob which will be the vector of estimated robust regression
% coefficients
brob=-99;

options=struct('intercept',1,'nsamp',nsampdef,'h',hdef,'bdp',...
    bdpdef,'lms',1,'rew',0,'plots',0,'nocheck',0,'nomes',0,...
    'conflev',0.975,'msg',1,'yxsave',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % If 'lms' is part of UserOptions and it is 2 or it is a structure (fastLTS),
    % then decrease the default for nsamp (nsampdef and options.nsamp) from
    % 3000 to 1000.
    checklms2 = strcmp(UserOptions,'lms');
    if sum(checklms2)
        lmsval = vvarargin(2*find(checklms2));
        if isstruct(lmsval{1}) || lmsval{1}==2
            nsampdef = min(1000,ncomb);
            options.nsamp = nsampdef;
        end
    end
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:LXS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:LXS:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin > 2
    
    % Extract the names of the optional arguments
    chklist=varargin(1:2:length(varargin));
    
    % Check whether the user has selected both h and bdp.
    chktrim=sum(strcmp(chklist,'h')+2*strcmp(chklist,'bdp'));
    if chktrim ==3
        error('FSDA:LXS:TooManyInputArgs','Both input arguments bdp and h are provided. Only one is required.')
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    % And check if the optional user parameters are reasonable.
    
    % Check h and bdp
    % The user has only specified h: no need to specify bdp.
    if chktrim==1
        if options.h < hmin
            error('FSDA:LXS:Wrongh',['The LTS (LMS) must cover at least ' int2str(hmin) ' observations.'])
        elseif options.h >= n
            error('FSDA:LXS:Wrongh','h is greater or equal to the number of non-missings and non-infinites.')
        end
        
        % the user has only specified bdp: h is defined accordingly
    elseif chktrim==2
        if options.bdp <= 0
            error('FSDA:LXS:WrongBdp','Attention: bdp should be larger than 0');
        end
        
        % nalpha must be greater or equal than 50% of the observations
        nalpha=ceil(n*(1-options.bdp));
        
        if nalpha <p+1
            error('FSDA:LXS:WrongAlpha',['Attention: the specified trimming proportion is too high.\n',...
                'It is necessary to specify bdp in such way that n*(1-options.bdp)>=p+1.']);
        end
        options.h=nalpha;
    end
    
    % Check number of subsamples to extract
    if options.nsamp>ncomb
        if options.msg==1
            disp('Number of subsets to extract greater than (n p). It is set to (n p)');
        end
        options.nsamp=0;
    elseif  options.nsamp<0
        error('FSDA:LXS:WrongNsamp','Number of subsets to extract must be 0 (all) or a positive number');
    end
end

% Default values for the optional
% parameters are set inside structure 'options'

h=options.h;                % Number of data points on which estimates are based
if isempty(h)
    h=hdef;
end

plots=options.plots;        % Plot of residuals equal to 1
nsamp=options.nsamp;        % Number of subsets to extract
if isempty(nsamp)
    nsamp=nsampdef;
end

lms=options.lms;            % if options.lms==1 then LMS, else LTS

if ~isstruct(lms) && lms==2
    lms=struct;
    refsteps=3;
    reftol=1e-6;
    bestr=5;
    refstepsbestr=50;
    reftolbestr=1e-8;
    % ij is a scalar used to ensure that the best first bestr non singular
    % subsets are stored
    ij=1;
    
    % lmsopt is associated with the message about total computing time
    lmsopt=2;
    
elseif isstruct(lms)
    lmsdef.refsteps=3;
    lmsdef.reftol=1e-6;
    lmsdef.bestr=5;
    lmsdef.refstepsbestr=50;
    lmsdef.reftolbestr=1e-8;
    
    % Control the appearance of the trajectories to be highlighted
    if ~isequal(lms,lmsdef)
        
        fld=fieldnames(lms);
        
        % Check if user options inside options.fground are valid options
        chkoptions(lmsdef,fld)
        for i=1:length(fld)
            lmsdef.(fld{i})=lms.(fld{i});
        end
    end
    
    % For the options not set by the user use their default value
    lms=lmsdef;
    
    refsteps=lms.refsteps;
    reftol=lms.reftol;
    bestr=lms.bestr;
    refstepsbestr=lms.refstepsbestr;
    reftolbestr=lms.reftolbestr;
    
    bestbetas = zeros(bestr,p);
    bestsubset = bestbetas;
    bestscales = Inf * ones(bestr,1);
    sworst = Inf;
    
    % ij is a scalar used to ensure that the best first bestr non singular
    % subsets are stored
    ij=1;
    % lmsopt is associated with the message about total computing time
    lmsopt=2;
else
    % lmsopt is associated with the message about total computing time
    if lms==1
        lmsopt=1;
    else
        lmsopt=0;
    end
end


rew=options.rew;            % if options.rew==1 use reweighted version of LMS/LTS,
conflev=options.conflev;    % Confidence level which is used for outlier detection
conflev=(conflev+1)/2;
msg=options.msg;            % Scalar which controls the messages displayed on the screen

nomes=options.nomes;        % if options.nomes==1 no message about estimated time to compute LMS is displayed

% Get user values of warnings
warnrank=warning('query','MATLAB:rankDeficientMatrix');
warnsing=warning('query','MATLAB:singularMatrix');
warnnear=warning('query','MATLAB:nearlySingularMatrix');
% Set them to off inside this function
% At the end of the file they will be restored to previous values
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');


%% Extract in the rows of matrix C the indexes of all required subsets
[C,nselected] = subsets(nsamp,n,p,ncomb,msg);
% Store the indices in varargout
if nargout==2
    varargout={C};
end


% rmin will contain the minimum value of LMS (LTS)
rmin=Inf;

% initialise and start timer.
tsampling = ceil(min(nselected/100 , 1000));
time=zeros(tsampling,1);

%% Computation of LMS (LTS)
for i=1:nselected
    if i <= tsampling, ttic = tic; end
    
    % extract a subset of size p
    index = C(i,:);
    
    
    Xb=X(index,:);
    yb=y(index);
    
    % if rank(Xb)==p Warning: this instruction has been replaced by a
    % posteriori control on vector b
    % Compute the vector of coefficients using matrice Xb and yb
    b=Xb\yb;
    
    if ~isnan(b(1)) && ~isinf(b(1))
        
        if ~isstruct(lms)
            % Residuals for all observations using b based on subset
            r=y-X*b;
            
            % Squared residuals for all the observations
            r2=r.^2;
            
            % Ordering of squared residuals
            r2s  = sort(r2);
            
            if lms==1
                % LMS
                rrob=r2s(h);
            else
                % STANDARD LTS without concentration steps
                rrob=sum(r2s(1:h));
            end
            
            if rrob<rmin
                % rmin = smallest ordered quantile or smallest truncated sum.
                rmin=rrob;
                
                % brob = \beta_lms or \beta_lts
                brob=b;
                
                % bs = units forming best subset according to lms or lts
                bs=index;
            end
            
            
        else % in this case the user has chosen the FAST LTS (with concentration steps)
            
            tmp = IRWLSreg(y,X,b,refsteps,reftol,h);
            
            betarw = tmp.betarw;
            numscale2rw = tmp.numscale2rw;
            
            if ij > bestr
                
                if numscale2rw < sworst
                    % Find position of the maximum value of previously stored
                    % best scales
                    
                    [~,ind] = max(bestscales);
                    
                    % Store numscale2rw, betarw and indexes of the units forming the
                    % best subset for the current iteration
                    bestscales(ind)     = numscale2rw;
                    bestbetas(ind,:)    = betarw';
                    bestsubset(ind,:)   = index;
                    % sworst = the best scale among the bestr found up to now
                    sworst              = max(bestscales);
                end
            else
                bestscales(ij)  = numscale2rw;
                bestbetas(ij,:) = betarw';
                bestsubset(ij,:)= index;
                % sworst = the best scale among the bestr found up to now
                sworst = max(bestscales);
                ij = ij+1;
                brob = 1;
            end
            
        end
        
        
    else
        singsub=singsub+1;
    end
    
    if ~nomes
        if i <= tsampling
            
            % sampling time until step tsampling
            time(i) = toc(ttic);
        elseif i==tsampling+1
            % stop sampling and print the estimated time
            if msg==1
                switch lmsopt
                    case 1
                        fprintf('Total estimated time to complete LMS: %5.2f seconds \n', nselected*median(time));
                    case 2
                        fprintf('Total estimated time to complete FASTLTS: %5.2f seconds \n', nselected*median(time));
                    otherwise
                        fprintf('Total estimated time to complete LTS: %5.2f seconds \n', nselected*median(time));
                end
            end
        end
    end
end

if brob==-99
    error('FSDA:LXS:NoFullRank','No subset had full rank. Please increase the number of subsets or check your design matrix X')
else
end

if isstruct(lms)
    % perform C-steps on best 'bestr' solutions, till convergence or for a
    % maximum of refstepsbestr steps using a convergence tolerance as specified
    % by scalar reftolbestr
    
    superbestscale = Inf;
    
    for i=1:bestr
        tmp = IRWLSreg(y,X,bestbetas(i,:)',refstepsbestr,reftolbestr,h);
        
        if tmp.numscale2rw < superbestscale
            % sh0 = superbestscale
            sh0 = tmp.numscale2rw;
            % brob = superbestbeta
            brob = tmp.betarw;
            % bs = superbestsubset, units forming best subset according to fastlts
            bs = bestsubset(i,:);
            superbestscale=sh0;
        end
    end
    
    % Pass from numerator of squared estimate of the scale to proper scale
    % estimate
    sh0=sqrt(sh0/h);
else
    
    if lms==1
        
        % Estimate of scale based on h-quantile of all squared residuals
        sh0=sqrt(rmin);
    else
        
        % Estimate of scale based on the first h squared smallest residuals
        sh0=sqrt(rmin/h);
    end
    
end


% residuals = Raw residuals using robust estimate of beta
residuals=y-X*brob;


% Consistency factor
if lmsopt==1
    % Consistency factor based on the fact that med |z_i|/ \Phi^-1(0.75) is a
    % consistent estimator of \sigma when z_i ~ N(0, \sigma^2)
    % The additional factor 1+5/(n-p) was found by simulation by Rousseeuw and
    % Leroy (1987), see p. 202
    factor=1.4826*(1+5/(n-p));
    
    % Apply the consistency factor to the preliminary scale estimate
    s0=sh0*factor;
    
else
    % Consistency factor based on the variance of the truncated normal distribution.
    % 1-h/n=trimming percentage
    % Compute variance of the truncated normal distribution.
    a=norminv(0.5*(1+h/n));
    %factor=1/sqrt(1-(2*a.*normpdf(a))./(2*normcdf(a)-1));
    factor=1/sqrt(1-2*(n/h)*a.*normpdf(a));
    
    % Note that factor=sqrt(factor1)
    %     v=1;
    %     a=chi2inv(h/n,1);
    %     factor1=(h/n)/(chi2cdf(a,1+2));
    
    % Apply the asymptotic consistency factor to the preliminary scale estimate
    s0=sh0*factor;
    
    % Apply small sample correction factor of Pison et al.
    if h<n
        s0=s0*sqrt(corfactorRAW(1,n,h/n));
    end
    
    %         % Analysis of the small sample correction factor of Pison et al.
    %         rangen=20:100;
    %         corf=zeros(length(rangen),1);
    %         for i=1:length(rangen)
    %             corf(i)=sqrt(corfactorRAW(1,rangen(i),0.7));
    %         end
    %         plot(rangen',corf)
    %         disp('s0 after')
    %         disp(s0)
    
    
end

if abs(s0) > 1e-7
    
    % Assign weight=1 to the h units which show the smallest h squared
    % residuals
    [~ , indsorres2] = sort(residuals.^2);
    weights = zeros(n,1);
    weights(indsorres2(1:h)) = 1;
    
    % Initialize structure out
    out=struct;
    
    % Store inside structure out, the vector of the weights
    out.weights=weights;
    
    % Compute the Student T quantile threshold . If options.conflev=0.975,
    % 1.25% on the right and
    % 1.25% on the left, globally 2.5%.
    % m = factor which modifies the degrees of freedom used for computing
    % the quantile.
    m = 2*n / asvar(h,n);
    quantile=tinv(conflev,m);
    
    % Observations with a standardized residual smaller than the quantile
    % threshold have a weight equal to 1, else the weight is equal to 0.
    % REMARK: using this threshold, even if the sample is homogeneous,
    % you are willing to declare at least 2.5% units as outliers.
    % Remark: sqrt(chi2inv(0.975,1)) = tinv(0.9875,\infinity) = quantile
    stdres = residuals/s0;
    weights = abs(stdres)<=quantile;
    % weights is a boolean vector.
    
    
    %% Reweighting part
    if rew==1
        
        % Find new estimate of beta using only observations which have
        % weight equal to 1. Notice that new brob overwrites old brob
        % computed previously.
        
        brob = X(weights==1,:) \ y(weights==1);
        % The QR decomposition is equivalent to the above but less efficient:
        % [Q,R]=qr(X(weights==1,:),0);
        % brob = R\(Q'*y(weights==1));
        
        % Computation of reweighted residuals.
        residuals=y-X*brob;
        % Find new estimate of scale using only observations which have
        % weight equal to 1.
        
        s0=sqrt(sum(weights.*residuals.^2)/(sum(weights)-1));
        % Compute new standardized residuals.
        
        % Apply consistency factor to reweighted estimate of sigma
        hrew=sum(weights);
        if hrew<n
            % factor=consistencyfactor(hrew,n,1);
            a=norminv(0.5*(1+hrew/n));
            %factor=1/sqrt(1-(2*a.*normpdf(a))./(2*normcdf(a)-1));
            factor=1/sqrt(1-2*(n/hrew)*a.*normpdf(a));
            % Apply small sample correction factor to reweighted estimate of sigma
            factor=factor*sqrt(corfactorREW(1,n,hrew/n));
        else
            factor=1;
        end
        
        s0=s0*factor;
        stdres=residuals/s0;
        
        % Declare as outliers the observations which have a standardized
        % residual greater than cutoff.
        % REMARK: while the first threshold was based on the Student T
        % (with modified degrees of freedom), in this second round the
        % threshold is based on the Normal. Notice that:
        % sqrt(chi2inv(0.975,1)) = tinv(0.9875,\infinity) = norminv(0.9875)
        
        weights = abs(stdres)<=norminv(conflev);
        % The new vector of weights overwrites previous vector of weigths
        % before reweighting.
        
        % Store information about reweighting
        out.rew=1;
        
    else
        % The default is no reweighting
        out.rew=0;
        
    end
    
else % Perfect fit
    if msg==1
        disp('Attention: there was an exact fit. Robust estimate of s^2 is <1e-7')
    end
    % There is an approximate perfect fit for the first h observations.
    % We consider as outliers all units with residual greater than 1e-7.
    weights = abs(residuals)<=1e-7;
    
    % Store the weights
    out.weights=weights;
    
    % s is set to 0
    s0=0;
    
    % Standardized residuals are artificially set equal to raw residuals.
    stdres=residuals;
end

%% Store quantities in the out structure


% Store robust estimate of beta coefficients
out.beta=brob;

% Store robust estimate of s
out.scale=s0;

% Store standardized residuals
out.residuals=stdres;

% Store units forming initial subset
out.bs=bs;

% Store list of units declared as outliers
out.outliers=seq(weights==0);

% Store confidence level which is used to draw the horizontal lines in the
% plot
out.conflev=options.conflev;

% Store the number of observations that have determined the LTS (LMS)
% estimator, i.e. the value of h.
out.h=h;

% Store number of singular subsets
out.singsub=singsub;
if msg==1
    if singsub/nselected>0.1
        disp('------------------------------')
        disp(['Warning: Number of subsets without full rank equal to ' num2str(100*singsub/nselected) '%'])
    end
end
% Store information about the class of the object
if lmsopt==1
    out.class='LMS';
else
    out.class='LTS';
end

if options.yxsave
    if options.intercept==1
        % Store X (without the column of ones if there is an intercept)
        out.X=X(:,2:end);
    else
        out.X=X;
    end
    % Store response
    out.y=y;
end


%% Create plots
% If plots is a structure, plot directly those chosen by the user;
% elseif plots is 1 a plot or residuals against index number appears
% else no plot is produced.
if plots==1
    if lmsopt==1
        laby='Robust lms residuals';
    else
        laby='Robust lts residuals';
    end
    resindexplot(out.residuals,'conflev',options.conflev,'laby',laby,'numlab',out.outliers);
else
end

% Restore the previous state of the warnings
warning(warnrank.state,'MATLAB:rankDeficientMatrix');
warning(warnsing.state,'MATLAB:singularMatrix');
warning(warnnear.state,'MATLAB:nearlySingularMatrix');



%% The part below contains subfunctions which are used only inside this file

    function asymptvar=asvar(h,n)
        %asvar computes the new degrees of freedom for the Student T
        hn=h/n;
        qalfa=chi2inv(hn,1);
        c2=gamcdf(qalfa/2,3/2);
        c1=1/c2;
        c3=3*gamcdf(qalfa/2,5/2);
        asymptvar=(qalfa*hn-c2)^2;
        asymptvar=(c3-2*qalfa*c2+hn*(qalfa^2))-asymptvar;
        asymptvar=c1^2*asymptvar;
    end



end



% -------------------------------------------------------------------
% subfunction IRWLSreg
% -------------------------------------------------------------------

function outIRWLS = IRWLSreg(y,X,initialbeta,refsteps,reftol,h)
%IRWLSreg (iterative reweighted least squares) does refsteps refining steps from initialbeta
%
%  Required input arguments:
%
%    y:         A vector with n elements that contains the response variable.
%               It can be both a row or column vector.
%    X :        Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p). Rows of X represent observations, and
%               columns represent variables.
% initialbeta : p x 1 vector containing initial estimate of beta
%   refsteps  : scalar, number of refining (IRLS) steps
%   reftol    : relative convergence tolerance
%               Default value is 1e-7
%      h      : scalar. number of observations with smallest residuals to consider
%
%  Output:
%
%  The output consists of a structure 'outIRWLS' containing the following fields:
%      betarw  : p x 1 vector. Estimate of beta after refsteps refining steps
%  numscale2rw : scalar. Sum of the smallest h squared residuals from
%                final iteration (after refsteps refining step).It is the
%                numerator of the estimate of the squared scale.
%     weights  : n x 1 vector. Weights assigned to each observation
%               In this case weights are 0,1.
%               1 for the units associated with the smallest h squared residuals from
%               final iteration
%               0 for the other units.
%

n=size(y,1);

% Residuals for the initialbeta
res = y - X * initialbeta;

% Squared residuals for all the observations
r2 = res.^2;

% Ordering of squared residuals
[r2s , i_r2s] = sort(r2);
initialscale  = sum(r2s(1:h));

% Initialize parameters for the refining steps loop
iter        = 0;
betadiff    = 9999;
beta        = initialbeta;
scale       = Inf;

% update of weights moved at the end of the function
% weights=zeros(n,1);
% weights(i_r2s(1:h))=1;

while ( (betadiff > reftol) && (iter < refsteps) )
    iter = iter + 1;
    
    % i_r2s= units with smallest h squared residuals
    i_r2s = i_r2s(1:h);
    % new coefficients based on units with smallest h squared
    % residuals
    newbeta = X(i_r2s,:) \ y(i_r2s);
    
    % exit from the loop if the new beta has singular values. In such a
    % case, any intermediate estimate is not reliable and we can just
    % keep the initialbeta and initial scale.
    if (any(isnan(newbeta)))
        newbeta = initialbeta;
        scale = initialscale;
        break
    end
    
    % betadiff is linked to the tolerance (specified in scalar reftol)
    betadiff = norm(beta - newbeta,1) / norm(beta,1);
    
    % update residuals
    res = y - X * newbeta;
    % Ordering of all new squared residuals
    [r2s , i_r2s] = sort(res.^2);
    % sum of the smallest new squared residuals
    scale = sum(r2s(1:h));
    % update beta
    beta = newbeta;
    
end

% store final estimate of beta
outIRWLS.betarw = newbeta;

% store final estimate of scale
outIRWLS.numscale2rw = scale;

% store final estimate of the weights for each observation
% In this case weights are 0,1.
% 1 for the units associated with the smallest h squared residuals from
% final iteration
% 0 for the other units.
weights=zeros(n,1);
weights(i_r2s(1:h))=1;
outIRWLS.weights=weights;

end



%% corfactorRAW function
function rawcorfac=corfactorRAW(p,n,alpha)

if p > 2
    coeffqpkwad875=[-0.455179464070565,1.11192541278794,2;-0.294241208320834,1.09649329149811,3]';
    coeffqpkwad500=[-1.42764571687802,1.26263336932151,2;-1.06141115981725,1.28907991440387,3]';
    y1_500=1+(coeffqpkwad500(1,1)*1)/p^coeffqpkwad500(2,1);
    y2_500=1+(coeffqpkwad500(1,2)*1)/p^coeffqpkwad500(2,2);
    y1_875=1+(coeffqpkwad875(1,1)*1)/p^coeffqpkwad875(2,1);
    y2_875=1+(coeffqpkwad875(1,2)*1)/p^coeffqpkwad875(2,2);
    y1_500=log(1-y1_500);
    y2_500=log(1-y2_500);
    y_500=[y1_500;y2_500];
    A_500=[1,log(1/(coeffqpkwad500(3,1)*p^2));1,log(1/(coeffqpkwad500(3,2)*p^2))];
    coeffic_500=A_500\y_500;
    y1_875=log(1-y1_875);
    y2_875=log(1-y2_875);
    y_875=[y1_875;y2_875];
    A_875=[1,log(1/(coeffqpkwad875(3,1)*p^2));1,log(1/(coeffqpkwad875(3,2)*p^2))];
    coeffic_875=A_875\y_875;
    fp_500_n=1-(exp(coeffic_500(1))*1)/n^coeffic_500(2);
    fp_875_n=1-(exp(coeffic_875(1))*1)/n^coeffic_875(2);
else
    if p == 2
        fp_500_n=1-(exp(0.673292623522027)*1)/n^0.691365864961895;
        fp_875_n=1-(exp(0.446537815635445)*1)/n^1.06690782995919;
    end
    if p == 1
        fp_500_n=1-(exp(0.262024211897096)*1)/n^0.604756680630497;
        fp_875_n=1-(exp(-0.351584646688712)*1)/n^1.01646567502486;
    end
end
if 0.5 <= alpha && alpha <= 0.875
    fp_alpha_n=fp_500_n+(fp_875_n-fp_500_n)/0.375*(alpha-0.5);
elseif 0.875 < alpha && alpha < 1
    fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
else
    error('FSDA:LXS:WrongBdp','Condition 1-alpha>=0.5 not respected')
end

rawcorfac=1/fp_alpha_n;
if rawcorfac <=0 || rawcorfac>50
    rawcorfac=1;
    if msg==1
        disp('Warning: problem in subfunction corfactorRAW')
        disp(['Correction factor for covariance matrix based on simulations found =' num2str(rawcorfac)])
        disp('Given that this value is clearly wrong we put it equal to 1 (no correction)')
        disp('This may happen when n is very small and p is large')
    end
end
end

%% corfactorREW function
function rewcorfac=corfactorREW(p,n,alpha)

if p > 2
    coeffrewqpkwad875=[-0.544482443573914,1.25994483222292,2;-0.343791072183285,1.25159004257133,3]';
    coeffrewqpkwad500=[-1.02842572724793,1.67659883081926,2;-0.26800273450853,1.35968562893582,3]';
    y1_500=1+(coeffrewqpkwad500(1,1)*1)/p^coeffrewqpkwad500(2,1);
    y2_500=1+(coeffrewqpkwad500(1,2)*1)/p^coeffrewqpkwad500(2,2);
    y1_875=1+(coeffrewqpkwad875(1,1)*1)/p^coeffrewqpkwad875(2,1);
    y2_875=1+(coeffrewqpkwad875(1,2)*1)/p^coeffrewqpkwad875(2,2);
    y1_500=log(1-y1_500);
    y2_500=log(1-y2_500);
    y_500=[y1_500;y2_500];
    A_500=[1,log(1/(coeffrewqpkwad500(3,1)*p^2));1,log(1/(coeffrewqpkwad500(3,2)*p^2))];
    coeffic_500=A_500\y_500;
    y1_875=log(1-y1_875);
    y2_875=log(1-y2_875);
    y_875=[y1_875;y2_875];
    A_875=[1,log(1/(coeffrewqpkwad875(3,1)*p^2));1,log(1/(coeffrewqpkwad875(3,2)*p^2))];
    coeffic_875=A_875\y_875;
    fp_500_n=1-(exp(coeffic_500(1))*1)/n^coeffic_500(2);
    fp_875_n=1-(exp(coeffic_875(1))*1)/n^coeffic_875(2);
else
    if p == 2
        fp_500_n=1-(exp(3.11101712909049)*1)/n^1.91401056721863;
        fp_875_n=1-(exp(0.79473550581058)*1)/n^1.10081930350091;
    end
    if p == 1
        fp_500_n=1-(exp(1.11098143415027)*1)/n^1.5182890270453;
        fp_875_n=1-(exp(-0.66046776772861)*1)/n^0.88939595831888;
    end
end
if 0.5 <= alpha && alpha <= 0.875
    fp_alpha_n=fp_500_n+(fp_875_n-fp_500_n)/0.375*(alpha-0.5);
end
if 0.875 < alpha && alpha < 1
    fp_alpha_n=fp_875_n+(1-fp_875_n)/0.125*(alpha-0.875);
end
rewcorfac=1/fp_alpha_n;
if rewcorfac <=0 || rewcorfac>50
    rewcorfac=1;
    if msg==1
        disp('Warning: problem in subfunction corfactorREW');
        disp(['Correction factor for covariance matrix based on simulations found =' num2str(rewcorfac)]);
        disp('Given that this value is clearly wrong we put it equal to 1 (no correction)');
        disp('This may happen when n is very small and p is large');
    end
end
end
%FScategory:REG-Regression