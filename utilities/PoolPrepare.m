function [numpool, tstart, progbar, usePCT, usematlabpool] = PoolPrepare(numpool,pariter,UserOptions)
%PoolPrepare prepares a pool of MATLAB instances for executing code in parallel
%
%<a href="matlab: docsearchFS('PoolPrepare')">Link to the help function</a>
%
% PoolPrepare and PoolClose are used respectively to open and close a
% prespecified number of parallel MATLAB sessions, which need to be
% distributed over the physical cores where MATLAB is running.
%
%
%  Required input arguments:
%
%      numpool:     The number of parallel sessions to open. Integer. If
%                   numpool is not defined, then it is set equal to the
%                   number of physical cores in the computer.
%                   Data Types - scalar
%
%      pariter:     Number of parfor loops that need to be monitored.
%                   Integer. If pariter > 0, then the 'pariter' parallel
%                   instancies executed in a parfor statement will be
%                   monitored with a progress bar. If pariter is 0 or is
%                   not defined, then the progression of the parallel
%                   execution is not monitored.
%                   Data Types - scalar
%
%  UserOptions:     Structure containing the user options of the calling
%                   function. Cell array of strings. It is used, for
%                   example, to check if the user has specified numpool or
%                   not, and proceed accordingly (i.e. use the number of
%                   workers set in the current MATLAB profile, rather then
%                   allocate the numpool MATLAB instances requested by the
%                   user).
%                   Data Types - cell array of strings
%
% Optional input arguments:
%
% Output:
%
%       numpool:    The number of parallel sessions actually opened. Integer. 
%                   They may differ from the request of the user, depending 
%                   on the computer configuration.
%                   Data Types - double
%
%        tstart:    Time stamp to be given as input to PoolClose. Double.  
%                   Records the internal computer time at the end of the
%                   execution of the PoolPrepare function, so that to
%                   monitor the overall execution time of the statements
%                   embedded between PoolPrepare and PoolClose.
%                   Data Types - double
%
%       progbar:    To be given as input to PoolClose. Structure or integer. 
%                   Contains the status of the progress bar used to monitor
%                   the progression of the parallel execution.
%                   Data Types - struct | double
%
% Optional Output:
%
%        usePCT:    Boolean indicating if the parallel computing toolbox is
%                   installed. Scalar {0,1}. Parpool checks for the
%                   existence of the parallel computing toolbox. 'usePCT'
%                   returns the result of the check to PoolClose, to avoid
%                   additional unnecessary checks.
%                   Data Types - integer | logical
%
% usematlabpool:    Boolean indicating the use of 'usematlabpool' or 'parpool'.
%                   Scalar {0,1}. Boolean indicating if the pool of MATLAB
%                   instances is created using 'matlabpool' or 'parpool',
%                   depending on the MATLAB version installed. From R2013b
%                   'parpool' is used. Earlier releases use 'usematlabpool'.
%                   Data Types - integer | logical
%
%
% See also: PoolClose, parfor
%
% References:
%
% Copyright 2008-2016.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('PoolPrepare')">Link to the help page for this function</a>
%
% Last modified 31-05-2016

% Examples:

%{
    % Sequential vs parallel run.
    n = 50000;
    x = randn(1,n) ;
    y = zeros(1,n);

    % sequential run
    tic
    for i = 1 : n
        y(i) = std(x(1:i));
    end
    fprintf('\n\n\n  Normal for: %f secs \n \n ',toc);

    % parallel run
    numpool = 4;
    pariter = n;
    UserOptions = {};
    [numpool, tstart, progbar, usePCT, usematlabpool] = ...
            PoolPrepare(numpool,pariter,UserOptions);

    parfor i = 1 : n
        y(i) = std(x(1:i));
    end

    cleanpool = 1; % this closes the pool of MATLAB sessions
    tend = PoolClose(cleanpool, tstart, progbar, usePCT,  usematlabpool);

    fprintf('\n\n\n      parFor: %f secs\n\n',tend);
%}

%{
    % Sequential vs parallel run (show tstart).
    n = 50000;
    x = randn(1,n) ;
    y = zeros(1,n);

    % sequential run
    tic
    for i = 1 : n
        y(i) = std(x(1:i));
    end
    fprintf('\n\n\n  Normal for: %f secs \n \n ',toc);

    % parallel run
    numpool = 4;
    pariter = n;
    UserOptions = {};
    [numpool, tstart, progbar, usePCT, usematlabpool] = ...
            PoolPrepare(numpool,pariter,UserOptions);
    disp(tstart)

    parfor i = 1 : n
        y(i) = std(x(1:i));
    end

    cleanpool = 1; % this closes the pool of MATLAB sessions
    tend = PoolClose(cleanpool, tstart, progbar, usePCT,  usematlabpool);

    fprintf('\n\n\n      parFor: %f secs\n\n',tend);
%}


%{
    % Sequential vs parallel run (show progbar).
    n = 50000;
    x = randn(1,n) ;
    y = zeros(1,n);

    % sequential run
    tic
    for i = 1 : n
        y(i) = std(x(1:i));
    end
    fprintf('\n\n\n  Normal for: %f secs \n \n ',toc);

    % parallel run
    numpool = 4;
    pariter = n;
    UserOptions = {};
    [numpool, tstart, progbar, usePCT, usematlabpool] = ...
            PoolPrepare(numpool,pariter,UserOptions);
    % show progrbar
    disp(progbar)

    parfor i = 1 : n
        y(i) = std(x(1:i));
    end

    cleanpool = 1; % this closes the pool of MATLAB sessions
    tend = PoolClose(cleanpool, tstart, progbar, usePCT,  usematlabpool);

    fprintf('\n\n\n      parFor: %f secs\n\n',tend);
%}

%{
    % Sequential vs parallel run (show usePCT).
    n = 50000;
    x = randn(1,n) ;
    y = zeros(1,n);

    % sequential run
    tic
    for i = 1 : n
        y(i) = std(x(1:i));
    end
    fprintf('\n\n\n  Normal for: %f secs \n \n ',toc);

    % parallel run
    numpool = 4;
    pariter = n;
    UserOptions = {};
    [numpool, tstart, progbar, usePCT, usematlabpool] = ...
            PoolPrepare(numpool,pariter,UserOptions);
    disp(usePCT)

    parfor i = 1 : n
        y(i) = std(x(1:i));
    end

    cleanpool = 1; % this closes the pool of MATLAB sessions
    tend = PoolClose(cleanpool, tstart, progbar, usePCT,  usematlabpool);

    fprintf('\n\n\n      parFor: %f secs\n\n',tend);
%}

%{
    % Sequential vs parallel run (show usematlabpool).
    n = 50000;
    x = randn(1,n) ;
    y = zeros(1,n);

    % sequential run
    tic
    for i = 1 : n
        y(i) = std(x(1:i));
    end
    fprintf('\n\n\n  Normal for: %f secs \n \n ',toc);

    % parallel run
    numpool = 4;
    pariter = n;
    UserOptions = {};
    [numpool, tstart, progbar, usePCT, usematlabpool] = ...
            PoolPrepare(numpool,pariter,UserOptions);
    disp(usematlabpool)

    parfor i = 1 : n
        y(i) = std(x(1:i));
    end

    cleanpool = 1; % this closes the pool of MATLAB sessions
    tend = PoolClose(cleanpool, tstart, progbar, usePCT,  usematlabpool);

    fprintf('\n\n\n      parFor: %f secs\n\n',tend);
%}

%% Beginning of code

% 'pariter' is the number of upcoming parallel calculations to monitor with
% a progress bar. If it is 0 or is not defined, then the progression of the
% parallel execution is not monitored.
if ~isempty(pariter) && isnumeric(pariter)
    pariter = round(pariter);
else
    pariter = 0;
end

% if numpool is not defined, check the number of physical cores in the
% computer. This is to make this function independent of previous checks.
if numpool==0 || isempty(numpool)
    numpool = feature('numCores');
end

% If the Parallel Computing Toolbox is installed, then either 'matlabpool'
% or (from R2013b) 'parpool' must exist, and 'usematlabpool' is set to 1 or
% 0 accordingly. Otherwise (i.e. if the Parallel Computing Toolbox is not
% installed) 'usematlabpool' is set to NaN.
if exist('parpool','file')==2
    usematlabpool = 0;
elseif exist('matlabpool','file')==2
    usematlabpool = 1;
else
    usematlabpool = nan;
end

if numpool > 1 && ~isnan(usematlabpool)
    usePCT = 1;
else
    usePCT = 0;
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
                pworkersLocProfile = parcluster();
                numpool = pworkersLocProfile.NumWorkers;
            end
            
            % Therefore if a parallel pool is not open,  UserOption numpool
            % (if set) overwrites the number of workers set in the
            % local/current profile. Similarly, the number of workers in
            % the local/current profile overwrites default value of
            % 'numpool' obtained as feature('numCores') (i.e. the number of
            % phisical cores)
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

if usePCT == 1 && pariter > 0
    progbar = ProgressBar(pariter);
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

% PoolPrepare and PoolClose monitor the overall execution time of the code
% (parallel and not) between the two instances, without counting the
% opening/close of the parpool
tstart = tic;

end

%FScategory:UTIGEN
