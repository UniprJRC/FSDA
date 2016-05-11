function [check , location] = isfunction(funstr)
%isfunction checks if a function exists. 
%
%<a href="matlab: docsearchFS('isfunction')">Link to the help function</a>
%
% Returns 1 if the function 'funstr' exists, 0 otherwise. Also returns the
% location of the function, which is empty if funstr does not exist. Works
% for both MATLAB and third party functions, e.g. those in FSDA. REMARK:
% does not work for built-in functions such as sin, cos, etc.
%
%  Required input arguments:
%
%      funstr:    Function name. String. The function to be checked.
%                 Data Types - String
%           
% Optional input arguments:
%
% Output:
%
%       check:   Flag indicating the file existance. Boolean {1,0}.
%                1 if function 'funstr' exists, 0 otherwise.
%                Data Types - Logical
%
%    location:   Function path. String.  The location of the function, 
%                which is empty if 'funstr'  does not exist. 
%                Data Types - String
%
% See also: exist, which
%
%
% References:
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('isfunction')">Link to the help page for this function</a>
%
% Last modified 15-Feb-2016

% Examples:
%{
    % check if a function of the statistical toolbox exists
    [check , location] = isfunction('regress')
    
    % check if a function of the FSDA toolbox exists
    [check , location] = isfunction('FSR')
%}

%% 

fhandle = str2func(funstr);
finfo = functions(fhandle);
if (~isempty(finfo.file))
    check = 1;
else
    check = 0;
end
location = finfo.file;
end
%FScategory:UTIGEN