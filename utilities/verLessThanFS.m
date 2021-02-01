function result = verLessThanFS(vernumber)
%verLessThanFS compares version of MATLAB to specified version number
%
%
%<a href="matlab: docsearchFS('verLessThanFS')">Link to the help page for this function</a>
%
% function verLessThanFS is much faster than verLessThan because it just
% checks the version of MATLAB and calls directly the relevant built in
% functions. The number containing the MATLAB version is cached for better
% performance. In order words, the first time verLessThanFS is called the
% number associated to the MATLAB version is stored in persistent variable
% named cachedMatlabVerFS
%
%
% Required input arguments:
%
% vernumber :  version of MATLAB to test. double.
%               Number containing the version of MATLAB to test again
%               current version.
%               Example - 8.3
%               Data Types - double
%
%
%
% Optional input arguments:
%
% Output:
%
%    result : True or false. Boolean. result is true if the current version of
%               version of MATLAB is older than the
%               version specified by vernumber, and false otherwise. 
%
%
%  See also verLessThan.m
%
% References:
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('verLessThanFS')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit
%
%
% Examples
%
%{
    % Test whether the current version is older than MATLAB version 8.4.
    numbertotest = 8.4;
    out=verLessThanFS(numbertotest);
    if out == true 
        disp(['current version is older than ' num2str(numbertotest)]);
    else
        disp(['current version is not older than ' num2str(numbertotest)]);
    end
%}

%% Beginning of code


if nargin<1
    error('FSDA:verLessThanFS:Missingnumber','MATLAB number to test must be specified');
end

% We cache the MATLAB version number (double format) for better performance.
persistent cachedMatlabVerFS;

if ~isempty(cachedMatlabVerFS)
    doubleMatlabversion=cachedMatlabVerFS;
else
    % locate the position of the Contents.m file
    currentFName =[matlabroot filesep 'toolbox' filesep 'matlab' filesep 'general' filesep 'Contents.m'];
    fid = fopen(currentFName,'r');
    fgets(fid);
    % get the line containing the line '% MATLAB Version X.X (RXXXX) dd-MMM-yyyy 
    charMatlabversion = fgets(fid);
    fclose(fid);
    % MATLAB version is in positions 18:21
    charMatlabversion = charMatlabversion(18:21);
    doubleMatlabversion=str2double(charMatlabversion);
    cachedMatlabVerFS = doubleMatlabversion;
end

if  doubleMatlabversion<vernumber
    result =true;
else
    result = false;
end
end
%FScategory:UTIGEN