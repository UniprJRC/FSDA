function result = existFS(FileName)
%Check if file exist and put answer in a cached persistent variable
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
%    result : True or false. Boolean. result is true if file exists in the MATLAB path.
%
%
%  See also verLessThan.m
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

% We cache the MATLAB version number (double format) for better performance.
persistent cachedexistFS;

if ~isempty(cachedexistFS)
    result=cachedexistFS;
else
    exFile=exist(FileName,'file');
    if exFile==3 || exFile==2
        result =true;
        cachedexistFS = result;
    end
    
end
%FScategory:UTIGEN