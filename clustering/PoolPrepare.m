function [numpool, tstart, progbar, usePCT, usematlabpool] = PoolPrepare(numpool,pariter,UserOptions)
%PoolPrepare prepares a pool of MATLAB instances for executing code in parallel 
%
%<a href="matlab: docsearchFS('PoolPrepare')">Link to the help function</a>
%
% Additional detailed comments here .......
%
%  Required input arguments:
%
%      numpool:    Function name. String. The function to be checked.
%
%      pariter:    
%
%  UserOptions:
%            
% Output:
%
%       numpool:   
%
%        tstart:    
%                
%       progbar:
%
% Optional Output:
%
%              usePCT:
%
%       usematlabpool:
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('PoolPrepare')">Link to the help page for this function</a>
%
% Last modified 15-Feb-2016

% Examples:
%{

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
% if isfunction('parpool')
%     usematlabpool = 0;
% elseif isfunction('matlabpool')
%     usematlabpool = 1;
% else
%     usematlabpool = nan;
% end

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

% % Use the Parallel Computing Toolbox if it is installed and if the parallel
% % pool available consists of more than one worker.
% vPCT = ver('distcomp');
% if numpool > 1 && ~isempty(vPCT)
%     usePCT = 1;
% else
%     usePCT = 0;
% end
%
% % Check for the MATLAB release in use. From R2013b, 'parpool' has taken the
% % place of 'matlabpool'.
% if verLessThan('matlab', '8.2.0')
%     usematlabpool = 1;
% else
%     usematlabpool = 0;
% end


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


