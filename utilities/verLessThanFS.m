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
% vector with two numbers associated to the MATLAB version is stored in
% persistent variable named cachedMatlabVerFS.
%
%
% Required input arguments:
%
% vernumber :  version of MATLAB to test. scalar double or character or vector of doubles.
%               Number containing the version of MATLAB to test again
%               current version. If version is a character, it must be in
%               the format major.minor.revision or major.minor.
%               If version is a scalar double what is before the full stop is
%               the major revision and the FIRST number after the full stop
%               is the minor revision. If version is a vector of doubles
%               the first element is the major revision and the second
%               element is the minor revision.
%               REMARK: Note that this function only considers major
%               and minor version and not the revision version.
%               Example - 8.3 or [8 3] or '8.3'
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
    % In this case the input argument is a scalar double
    numbertotest = 8.4;
    out=verLessThanFS(numbertotest);
    if out == true
        disp(['current version is older than ' num2str(numbertotest)]);
    else
        disp(['current version is not older than ' num2str(numbertotest)]);
    end
%}

%{
    % In this case the input argument is a vector with two elements 
    % Test whether the current version is older than MATLAB version 9.10.
    % Note that when the minor revision to test is greater than 2 it is necessary 
    % to use a vector with two elements.
    numbertotest = [9 10];
    out=verLessThanFS(numbertotest);
    if out == true
        disp(['current version is older than ' num2str(numbertotest)]);
    else
        disp(['current version is not older than ' num2str(numbertotest)]);
    end
%}

%{
    % Test whether the current version is older than MATLAB version 9.12
    % In this case the input argument is a vector with two elements 
    numbertotest = '9.11';
    out=verLessThanFS(numbertotest);
    if out == true
        disp(['current version is older than ' num2str(numbertotest)]);
    else
        disp(['current version is not oder than ' num2str(numbertotest)]);
    end
%}

%% Beginning of code


if nargin<1
    error('FSDA:verLessThanFS:Missingnumber','MATLAB number to test must be specified');
end

% We cache the MATLAB version number (double format) for better performance.
persistent cachedMatlabVerFS;

if isempty(cachedMatlabVerFS)
    % locate the position of the Contents.m file
    currentFName =[matlabroot filesep 'toolbox' filesep 'matlab' filesep 'general' filesep 'Contents.m'];
    fid = fopen(currentFName,'r');
    fgets(fid);
    % get the line containing the line '% MATLAB Version X.X (RXXXX) dd-MMM-yyyy
    charMatlabversion = fgets(fid);
    fclose(fid);
    % MATLAB version is in positions 18:21
    charMatlabversion = charMatlabversion(18:21);
    doubleMatlabversion = getParts(charMatlabversion);
    cachedMatlabVerFS = doubleMatlabversion;
else
    doubleMatlabversion=cachedMatlabVerFS;
end

if isnumeric(vernumber)
    if isscalar(vernumber)
        numberToTest(1)=floor(vernumber);
        str = num2str(vernumber,'%0.1f'); 
        str = str(find(str=='.')+1:end); %// keep only digits after decimal point
        numberToTest(2)=str2double(str);
    else
        numberToTest=vernumber;
    end
else
    numberToTest=getParts(vernumber);
end

if numberToTest(1) ~= doubleMatlabversion(1)     % major version
    result = doubleMatlabversion(1) <numberToTest(1);
else                                  % minor version
    result = doubleMatlabversion(2) < numberToTest(2);
end


function parts = getParts(V)
parts = sscanf(V, '%d.%d.%d')';
if length(parts) < 3
    parts(3) = 0; % zero-fills to 3 elements
end


%FScategory:UTIGEN