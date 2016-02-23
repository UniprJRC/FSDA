function [tend] = PoolClose(cleanpool, tstart, progbar, usePCT,  usematlabpool)
%PoolClose closes the pool of MATLAB instances opened with PoolPrepare to execute code in parallel 
%
%<a href="matlab: docsearchFS('PoolClose')">Link to the help function</a>
%
% Additional detailed comments here .......
%
%  Required input arguments:
%
%      cleanpool:    Function name. String. The function to be checked.
%
%         tstart:    
%
%        progbar:
%
%         usePCT:    
%
%
% Optional Input arguments:
%
%         usePCT:    
%
%  usematlabpool:    
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('PoolClose')">Link to the help page for this function</a>
%
% Last modified 15-Feb-2016

% Examples:
%{

%}
%

%% Beginning of code
% This is to make 'PoolClose' independent from the specification of usePCT
% and usematlabpool returned by 'PoolPrepare', i.e. if not provided in input.
if nargin < 4 || (isempty(usePCT) && isempty(usematlabpool))
    % If the Parallel Computing Toolbox is installed, then either 'matlabpool'
    % or (from R2013b) 'parpool' must exist, and 'usematlabpool' is set to 1 or
    % 0 accordingly. Otherwise (i.e. if the Parallel Computing Toolbox is not
    % installed) 'usematlabpool' is set to NaN.
    if isfunction('parpool')
        usematlabpool = 0;
    elseif isfunction('matlabpool')
        usematlabpool = 1;
    else
        usematlabpool = nan;
    end
    if numpool > 1 && ~isnan(usematlabpool)
        usePCT = 1;
    else
        usePCT = 0;
    end
end

% PoolPrepare and PoolClose monitor the overall execution time of the code
% (parallel and not) between the two instances, without counting the
% opening/close of the parpool
tend = toc(tstart);

if progbar ~= 9999
    if usePCT == 1
        progbar.stop;
    end
    disp(['Total time required: ' num2str(tend) ' seconds']);
end

% close parallel jobs if necessary
if usePCT == 1 && cleanpool == true
    if usematlabpool
        matlabpool('close'); %#ok<DPOOL>
    else
        delete(gcp);
    end
end
end