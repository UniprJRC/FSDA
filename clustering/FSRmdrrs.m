function [mdrrs,BBrs]=FSRmdrrs(y,X,varargin)
%FSRmdrrs performs random start monitoring of minimum deletion residual
%
% The trajectories originate from many different random initial subsets and
% provide information on the presence of groups in the data. Linear
% regression structures are investigated by monitoring the minimum deletion
% residual outside the FS subset.
%
%<a href="matlab: docsearchFS('FSRmdrrs')">Link to the help function</a>
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
%               will automatically be excluded from the computations.%
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
%   bsbsteps :  Save the units forming subsets. Vector. It specifies for
%               which steps of the fwd search it
%               is necessary to save the units forming subsets. If bsbsteps
%               is 0 we store the units forming subset in all steps. The
%               default is store the units forming subset in all steps if
%               n<=500, else to store the units forming subset at steps
%               init and steps which are multiple of 100. For example, as
%               default, if n=753 and init=6,
%               units forming subset are stored for
%               m=init, 100, 200, 300, 400, 500 and 600.
%               Example - 'bsbsteps',[100 200] stores the unis forming
%               subset in steps 100 and 200.
%               Data Types - double
%               REMARK: vector bsbsteps must contain numbers from init to
%               n. if min(bsbsteps)<init a warning message will appear on
%               the screen.
%     nsimul :  number of random starts. Scalar. The default value is200.
%               Example - 'nsimul',300
%               Data Types - double
%  nocheck:     Check input arguments. Scalar. If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additioanl column of ones for
%               the intercept is not added. As default nocheck=0. The
%               controls on h, alpha and nsamp still remain
%               Example - 'nocheck',1
%               Data Types - double
%  constr :     Constrained search. Vector. r x 1 vector which contains the list of units which are
%               forced to join the search in the last r steps. The default
%               is constr=''.  No constraint is imposed
%               Example - 'constr',[1:10] forces the first 10 units to join
%               the subset in the last 10 steps
%               Data Types - double
%  plots :      Plot on the screen. Scalar. If equal to one a plot of
%               random starts minimum deletion residual appears  on the
%               screen with 1%, 50% and 99% confidence bands else (default)
%               no plot is shown.
%               Example - 'plots',1
%               Data Types - double
%               Remark: the plot which is produced is very simple. In order
%               to control a series of options in this plot and in order to
%               connect it dynamically to the other forward plots it is
%               necessary to use function mdrrsplot.
%   numpool :  use parallel computing and parfor. Scalar.
%               If numpool > 1, the routine automatically checks
%               if the Parallel Computing Toolbox is installed and
%               distributes the random starts over numpool parallel
%               processes. If numpool <= 1, the random starts are run
%               sequentially. By default, numpool is set equal to the
%               number of physical cores available in the CPU (this choice
%               may be inconvenient if other applications are running
%               concurrently). The same happens if the numpool value
%               chosen by the user exceeds the available number of cores.
%               REMARK 1: up to R2013b, there was a limitation on the
%               maximum number of cores that could be addressed by the
%               parallel processing toolbox (8 and, more recently, 12).
%               From R2014a, it is possible to run a local cluster of more
%               than 12 workers.
%               REMARK 2: Unless you adjust the cluster profile, the
%               default maximum number of workers is the same as the
%               number of computational (physical) cores on the machine.
%               REMARK 3: In modern computers the number of logical cores
%               is larger than the number of physical cores. By default,
%               MATLAB is not using all logical cores because, normally,
%               hyper-threading is enabled and some cores are reserved to
%               this feature.
%               REMARK 4: It is because of Remarks 3 that we have chosen as
%               default value for numpool the number of physical cores
%               rather than the number of logical ones. The user can
%               increase the number of parallel pool workers allocated to
%               the multiple start monitoring by:
%               - setting the NumWorkers option in the local cluster profile
%                 settings to the number of logical cores (Remark 2). To do
%                 so go on the menu "Home|Parallel|Manage Cluster Profile"
%                 and set the desired "Number of workers to start on your
%                 local machine".
%               - setting numpool to the desired number of workers;
%               Therefore, *if a parallel pool is not already open*,
%               UserOption numpool (if set) overwrites the number of
%               workers set in the local/current profile. Similarly, the
%               number of workers in the local/current profile overwrites
%               default value of 'numpool' obtained as feature('numCores')
%               (i.e. the number of physical cores).
%               Example - 'numpool',8
%               Data Types - double
%  cleanpool :  clean pool after execution. Scalar.
%               cleanpool is 1 if the parallel pool has to be cleaned after
%               the execution of the random starts. Otherwise it is 0.
%               The default value of cleanpool is 1.
%               Clearly this option has an effect just if previous option
%               numpool is > 1.
%               Example - 'clarnpool',false
%               Data Types - boolean
%  msg  :       Level of output to display. Scalar.
%               Scalar which controls whether to display or not messages
%               about random start progress. More precisely, if previous
%               option numpool>1, then a progress bar is displayed, on
%               the other hand a message will be displayed on the screen
%               when 10%, 25%, 50%, 75% and 90% of the random starts have
%               been accomplished
%               In order to create the progress bar when nparpool>1
%               the program writes on a temporary .txt file in the folder
%               where the user is working. Therefore it is necessary to
%               work in a folder where the user has write permission. If this
%               is not the case and the user (say) is working without write
%               permission in folder C:\Program Files\MATLAB the following
%               message will appear on the screen:
%                   Error using ProgressBar (line 57)
%                   Do you have write permissions for C:\Program Files\MATLAB?
%               Example - 'msg',1
%               Data Types - double
%
%  Remark:      The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%
%               The dataset can include missing values (NaN's) and infinite
%               values (Inf's), since observations (rows) with missing or
%               infinite values will be automatically excluded from the
%               computations. y can be both a row of column vector.
%
% Output:
%
%       mdrrs:  random start minimum deletion residual. Matrix.
%               (n-init)-by-(nsimul+1) matrix which contains the monitoring
%               of minimum deletion residual at each step of the forward
%               search for each random start.
%               1st col = fwd search index (from init to n-1).
%               2nd col = minimum deletion residual for random start 1.
%               ...
%               nsimul+1 col = minimum deletion residual for random start nsimul.
%       BBrs :  units belonging to subset. 3D array.
%               3D array which contains the units belonging to the subset
%               at the steps specified by input option bsbsteps.
%               If bsbsteps=0 BBrs has size n-by-(n-init+1)-by-nsimul.
%               In this case BBrs(:,:,j) with j=1, 2, ..., nsimul
%               has the following structure:
%               1-st row has number 1 in correspondence of the steps in
%                   which unit 1 is included inside subset and a missing
%                   value for the other steps;
%               ......
%               (n-1)-th row has number n-1 in correspondence of the steps
%                   in which unit n-1 is included inside subset and a
%                   missing value for the other steps;
%               n-th row has the number n in correspondence of the steps in
%                   which unit n is included inside subset and a missing
%                   value for the other steps.
%               If, on the other hand, bsbsteps is a vector which specifies
%               the steps of the search in which it is necessary to store
%               subset, BBrs has size n-by-length(bsbsteps)-by-nsimul.
%               In other words, BBrs(:,:,j) with j=1, 2, ..., nsimul has
%               the same structure as before, but now contains just
%               length(bsbsteps) columns.
%
%
% See also:     FSRmdr, FSMmmdrs, FSMmmd
%
% References:
%
%   Atkinson A.C., Riani M., and Cerioli A. (2006). Random Start Forward
%   Searches with Envelopes for Detecting Clusters in Multivariate Data,
%   in: ZANI S., CERIOLI A., RIANI M., VICHI M. EDS. Data Analysis,
%   Classification and the Forward Search, (pp. 163-172), 
%   ISBN:3-540-35977-x, BERLIN: Springer Verlag (GERMANY).
%
%   Atkinson A.C., Riani M., (2007), Exploratory Tools for Clustering
%   Multivariate Data. COMPUTATIONAL STATISTICS & DATA ANALYSIS, vol. 52,
%   pp. 272-285 ISSN: 0167-9473, doi:10.1016/j.csda.2006.12.034
%
%   Riani M., Cerioli A., Atkinson A.C., Perrotta D., Torti F. (2008),
%   Fitting Mixtures of Regression Lines with the Forward Search, 
%   in: Mining Massive Data Sets for Security F. Fogelman-Soulie et al. EDS,
%   (pp. 271-286), IOS Press, Amsterdam (The Netherlands).
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRmdrrs')">Link to the help function</a>
% Last modified 31-05-2016
%
% Examples:
%

%{
    % An example with simulated data with regression lines
    % with roughly the same number of observations.
    close all
    rng('default')
    rng(2);
    b1=[1 1];
    b2=[1 2.6];
    n1=40;
    n2=50;
    s1=0.1;
    s2=0.1;
    X1=rand(n1,1);
    X2=rand(n2,1);
    y1=randn(n1,1)*s1+b1(1)+b1(2)*X1;
    y2=randn(n2,1)*s2+b2(1)+b2(2)*X2;
    hold('on')
    plot(X1,y1,'o');
    plot(X2,y2,'o');
    title('Two simulated regression lines')
    y=[y1;y2];
    X=[X1;X2];
    figure
    % parfor of Parallel Computing Toolbox is used (if present in current computer)
    % and pool is not cleaned after % the execution of the random starts
    [mdrrs,BBrs]=FSRmdrrs(y,X,'constr','','nsimul',100,'init',10,'plots',1,'cleanpool',false);
    disp('The two peaks in the trajectories of minimum deletion residual (mdr).')
    disp('clearly show the presence of two groups.')
    disp('The decrease after the peak in the trajectories of mdr is due to the masking effect.')

%}

%{
    % Same example as before but now the values of n1 and n2 (size of the
    % two groups) have been increased. In this case it is possible to see
    % that there are two trajectories of minimum deletion residual which go
    % outside the envelopes in the central part of the search.
    % In this case the two groups have roughly the same size (n1=140 and n2=150)
    close all
    rng(2);
    b1=[1 1];
    b2=[1 2.6];
    n1=140;
    n2=150;
    s1=0.1;
    s2=0.1;
    X1=rand(n1,1);
    X2=rand(n2,1);
    y1=randn(n1,1)*s1+b1(1)+b1(2)*X1;
    y2=randn(n2,1)*s2+b2(1)+b2(2)*X2;
    hold('on')
    plot(X1,y1,'o');
    plot(X2,y2,'o');
    title('Two simulated regression lines')
    y=[y1;y2];
    X=[X1;X2];
    figure
    % parfor of Parallel Computing Toolbox is used (if present in current
    % computer) and pool is not cleaned after
    % the execution of the random starts
    % The number of workers which is used is the one specified
    % in the local/current profile
    [out]=FSRmdrrs(y,X,'constr','','nsimul',100,'init',10,'plots',1,'cleanpool',false);
    disp('The two peaks in the trajectories of minimum deletion residual (mdr).')
    disp('clearly show the presence of two groups.')
    disp('The decrease after the peak in the trajectories of mdr is due to the masking effect.')
%}

%{
    % Same example as before but now there is one group which has a size
    % much greater than the other (n1=60 and n2=150). In this case it is
    % possible to see that there is a trajectory of minimum deletion
    % residual which goes outside the envelope in steps 60-110. This
    % corresponds to the searches initialized using the units coming from
    % the smaller group. Note that due to the partial overlapping after the
    % peak in steps 60-110 there is a gradual decrease. When m is around
    % 160, most of the units from this group tend to get out of the subset.
    % Therefore the value of mdr becomes much smaller than it should be.
    % Please note the dip around step m=165, which is due to entrance of the
    % units of the second larger group. This trajectory just after the dip
    % collapses into the trajectory which starts from the second group.
    % Around steps 90-110 it is also possible to see two trajectories
    % inside the bands which collaps into one around m=120. Please use
    % mdrrsplot with option databrush in order to explore the units
    % belonging to subset. Here we limit ourselves to notice that around m
    % =180 all the units from second group are included into subset (plus
    % some of group 1 given that the two groups partially overlap). Also
    % notice once again the decrease in the unique trajectory of minimum
    % deletion residual after m around 180 which is due to the entry of the
    % units of the first smaller group.
    close all
    rng(2);
    b1=[1 1];
    b2=[1 2.6];
    n1=60;
    n2=150;
    s1=0.1;
    s2=0.1;
    X1=rand(n1,1);
    X2=rand(n2,1);
    y1=randn(n1,1)*s1+b1(1)+b1(2)*X1;
    y2=randn(n2,1)*s2+b2(1)+b2(2)*X2;
    hold('on')
    plot(X1,y1,'o');
    plot(X2,y2,'o');
    title('Two simulated regression lines')
    y=[y1;y2];
    X=[X1;X2];
    figure
    % parfor of Parallel Computing Toolbox is used (if present in current
    % computer). Parallel pool is closed after the execution of the random starts
    [out]=FSRmdrrs(y,X,'constr','','nsimul',100,'init',10,'plots',1);
%}

%{
    % Random start for fishery dataset: two regression structures,
    % difficult to identify because of a dense area.
    load('fishery.txt');
    y=fishery(:,2);
    X=fishery(:,1);
    % parfor of Parallel Computing Toolbox is used (if installed)
    figure
    []=FSRmdrrs(y,X,'nsimul',100,'plots',1);
%}

%{
    % Random start for fishery dataset: just store information about the
    % units forming subset for each random start at specified steps
    load('fishery.txt');
    y=fishery(:,2);
    X=fishery(:,1);
    % parfor of Parallel Computing Toolbox is used (if present in current
    % computer)
    figure
    [mdrrs,BBrs]=FSRmdrrs(y,X,'nsimul',100,'plots',1,'bsbsteps',[10 300 600]);
    % sum(~isnan(BBrs(:,1,1)))
    %
    % ans =
    %
    %     10
    %
    % sum(~isnan(BBrs(:,2,1)))
    %
    % ans =
    %
    %    300
    %
    % sum(~isnan(BBrs(:,3,1)))
    %
    % ans =
    %
    %    600
%}

%{
    % Random start for fishery dataset: two regression structures,
    % difficult to identify because of a dense area.
    load('fishery.txt');
    y=fishery(:,2);
    X=fishery(:,1);
    % traditional for loop is used
    [mdrrs,BBrs]=FSRmdrrs(y,X,'nsimul',100,'plots',1,'numpool',0);
%}
%% Input parameters checking

nnargin   = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

% check how many physical cores are available in the computer (warning:
% function 'feature' is undocumented; however, FSDA is automatically
% monitored for errors and other inconsistencies at each new MATLAB
% release).
numpool = feature('numCores');

% Default for vector bsbsteps which indicates for which steps of the fwd
% search units forming subset have to be saved
initdef   = p+1;
if n<=500
    bsbstepdef = initdef:n;
else
    bsbstepdef = [initdef 100:100:100*floor(n/100)];
end

nsimuldef = 200; % nsimuldef = default number of random starts
options   = struct('intercept',1,'init',initdef,'plots',0,'nocheck',0,'msg',1,...
    'constr','','nsimul',nsimuldef,'numpool',numpool, 'cleanpool', true, ...
    'bsbsteps',bsbstepdef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRmdrrs:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    chkbsbsteps = strcmp(UserOptions,'bsbsteps');
else
    UserOptions=0;
end

init        = options.init;
intercept   = options.intercept;
msg         = options.msg;
constr      = options.constr;
plots       = options.plots;
nsimul      = options.nsimul;
cleanpool   = options.cleanpool;
numpool     = options.numpool;
bsbsteps    = options.bsbsteps;

if numpool < 1
    numpool = 1;
end

% Initialize structures to store statistics
if bsbsteps == 0
    BBrs  = zeros(n,n-init+1,nsimul);
else
    
    if ~isempty(bsbsteps(bsbsteps<init))
        % The following warning is shown just if the user has supplied vector
        % bsbsteps
        if sum(chkbsbsteps)
            warning('FSDA:FSRmdrrs:Wronginit','It is not possible to store subset for values of m smaller than init')
        end
        bsbsteps=bsbsteps(bsbsteps>=init);
    end
    
    BBrs=zeros(n,length(bsbsteps),nsimul);
end

mdrrs = [(init:n-1)' zeros(n-init,nsimul)];

%% Check MATLAB environment

% Use the Parallel Computing Toolbox if it is installed and if the parallel
% pool available consists of more than one worker.
vPCT = ver('distcomp');
if numpool > 1 && ~isempty(vPCT)
    usePCT = 1;
else
    usePCT = 0;
end

% Check for the MATLAB release in use. From R2013b, 'parpool' has taken the
% place of 'matlabpool'.
if verLessThan('matlab', '8.2.0')
    usematlabpool = 1;
else
    usematlabpool = 0;
end

%% Prepare the parallel pool

if usePCT==1 % In this case Parallel Computing Toolbox Exists
    
    % First check if there is a parallel pool open. If this is the case,
    % then the pool will be used. To keep it open for later reuse is
    % useful, as opening a pool takes some time. Variable 'pworkers' is 0
    % if there is no parallel pool open; otherwise it contains the number
    % of workers allocated for the parallel pool.
    if usematlabpool
        pworkers = matlabpool('size'); %#ok<DPOOL>
    else
        ppool = gcp('nocreate');
        if isempty(ppool)
            pworkers = 0;
            % If the user has not specified numpool, then the number of
            % workers which will be used is the one set in the current
            % profile
            if max(strcmp(UserOptions,'numpool')) ~= 1
                pworkersLocProfile=parcluster();
                numpool=pworkersLocProfile.NumWorkers;
            end
            
            % Therefore if a parallel pool is not open,  UserOption numpool
            % (if set) overwrites the number of workers set in the
            % local/current profile. Similarly, the number of workers in
            % the local/current profile overwrites default value of 'numpool' obtained as
            % feature('numCores') (i.e. the number of phisical cores)
        else
            pworkers = ppool.Cluster.NumWorkers;
        end
    end
    
    if pworkers > 0
        % If a parallel pool is already open, ensure that numpool is not
        % larger than the number of workers allocated to the parallel pool.
        numpool = min(numpool,pworkers);
    else
        % If there is no parallel pool open, create one with numpool workers.
        if usematlabpool
            matlabpool('open',numpool); %#ok<DPOOL>
        else
            parpool(numpool);
        end
    end
    
end

%% Monitoring minimum deletion residual with random starts

% monitor execution time, without counting the opening/close of the parpool
tstart = tic;

if usePCT == 1 && msg == 1
    progbar = ProgressBar(nsimul);
else
    % In the parfor, 'progbar' will not be instanciated if usePCT is 0. In
    % this case, as a measure of precaution, the MATLAB interpreter
    % generates an error, to force the user to treat the case. This
    % assignment is a workaround to avoid this type of error.
    progbar = 9999;
end

if numpool == 1
    % the following re-assignement of numpool from 1 to 0 is necessary,
    % because the 'parfor' statement with numpool = 1 opens a parallel
    % pool of 1 worker instead of reducing the iteration to a simple and
    % faster 'for' statement.
    numpool = 0;
end

parfor (j = 1:nsimul , numpool)
    [mdr,~,BB] = FSRmdr(y,X,0,'init',init,'intercept',intercept,...
        'nocheck',1,'msg',0,'constr',constr,'bsbsteps',bsbsteps);
    if size(mdr,2)>1
        % Store units forming subset at each step
        BBrs(:,:,j) = BB;
        % Store minimum deletion residual
        mdrrs(:,j+1) = mdr(:,2);
    else
        BBrs(:,:,j) = NaN;
        % Store minimum deletion residual
        mdrrs(:,j+1) = NaN;
    end
    if msg==1
        if usePCT == 1
            progbar.progress;  %#ok<PFBNS>
        else
            if j==nsimul/2 || j==nsimul/4  || j==nsimul*0.75 || j==nsimul*0.9
                disp(['Simul nr. ' num2str(nsimul-j) ' n=' num2str(n)]);
                % note that a parfor, used in sequential mode, iterates in
                % the inverse order (this explains why we display
                % the simulation number nsimul-j rather than j).
            end
        end
    end
    
end

tend = toc(tstart);

if msg==1
    if usePCT == 1
        progbar.stop;
    end
    disp(['Total time required by the multiple start monitoring: ' num2str(tend) ' seconds']);
end

% close parallel jobs if necessary
if usePCT == 1 && cleanpool == true
    if usematlabpool
        matlabpool('close'); %#ok<DPOOL>
    else
        delete(gcp);
    end
end

%% Plot statistic with random starts

if plots==1
    
    tagfig  = 'pl_mdrrs';
    ylab    = 'Minimum Deletion Residual';
    xlab    = 'Subset size m';
    p       = size(X,2);
    
    hold('on');
    
    % Plot lines of empirical quantiles
    plot(mdrrs(:,1), mdrrs(:,2:end), 'tag',tagfig);
    
    % Compute teoretical quantiles for minimum deletion residual using
    % order statistics
    quantilesT = FSRenvmdr(n,p,'exact',1,'init',init);
    
    % Plots lines of theoretical quantiles
    line(quantilesT(:,1),quantilesT(:,2:4), ...
        'LineStyle','-','Color','r','LineWidth',2,'tag','env');
    
    ax = get(gca,'YLim');
    if ax(2)>20
        mdrrstmp = mdrrs(:,2:end);
        mdrrstmp(mdrrstmp>20) = NaN;
        maxylim  = max(max(mdrrstmp));
        minylim  = min(min(mdrrstmp));
        set(gca,'YLim',[minylim maxylim]);
    end
    
    % axes labels
    xlabel(xlab);
    ylabel(ylab);
    
end

end
%FScategory:CLUS-RobClaREG