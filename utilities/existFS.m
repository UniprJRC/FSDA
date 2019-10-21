function result = existFS(FileName)
%Check if file exists and puts answer in a cached persistent variable
%
%
%<a href="matlab: docsearchFS('existFS')">Link to the help function</a>
%
%
% The answer (true/number) is cached for better
% performance. In order words, the first time existFS is called the
% answer is stored in persistent variable
% named cachedexistFS
%
%
% Required input arguments:
%
% FileName :  FileName whose existence has to be verified. double.
%               Example - 'myFile.mex'
%               Data Types - char
%
%
%
% Optional input arguments:
%
% Output:
%
%    result : True or false. Boolean. result is true if file exists in the
%               MATLAB path.
%
%
%  See also verLessThan.m, verLessThanFS.m
%
% References:
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('existFS')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit
%
%
% Examples
%
%{
    % Test whether a mex file exists.
    out=existFS('DfM.mexw64');
    if out == true
        disp('Mex file exists');
    else
        disp('Mex file does not exist');
    end
%}

%% Beginning of code


if nargin<1
    error('FSDA:existFS:MissingFile','FileName to test must be specified');
end

% We cache the presence of the mex file for better performance. That is we
% avoid calling function exist for each iteration
persistent cachedexistFS;

if ~isempty(cachedexistFS)
    result=cachedexistFS;
else
    exFile=exist(FileName,'file');
    if exFile==3 || exFile==2
        result =true;
        cachedexistFS = result;
    else
        result = false;
        cachedexistFS = false;
    end
    
end
%FScategory:UTIGEN