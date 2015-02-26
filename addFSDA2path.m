function addFSDA2path(FSDApath)
%Add FSDA toolbox to path
%
%  Required input arguments:
%
%    FSDA:      A string containing the path which contains the root folder of FSDA toolbox
%
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
% Last modified 06-Feb-2015

% Examples:

%
%{
%      If FSDA has been installed in D:\matlab\FSDA
%      in order to include the required subfolders of FSDA use the following sintax
 
       addFSDA2path('D:\matlab\FSDA')
%}

%{
        % The expression fileparts(which('docsearchFS.m')) locates the main folder
        % where FSDA is installed
        addFSDA2path(fileparts(which('docsearchFS.m')))
%}

%% Beginning of code

if ispc ==1
    addp=[FSDApath ';' FSDApath '\multivariate;'  FSDApath '\regression;',...
        FSDApath '\datasets\regression;' FSDApath '\datasets\multivariate;', ...
        FSDApath '\datasets\multivariate_regression;' FSDApath '\graphics;' FSDApath '\utilities;',...
        FSDApath '\examples;' FSDApath '\robust;' FSDApath '\combinatorial;', ...
        FSDApath '\clustering;' FSDApath '\datasets\clustering;' FSDApath '\FSDAdemos'];
else
    addp=[FSDApath ';' FSDApath '/multivariate;'  FSDApath '/regression;',...
        FSDApath '/datasets/regression;' FSDApath '/datasets/multivariate;', ...
        FSDApath '/datasets/multivariate_regression;' FSDApath '/graphics;' FSDApath '/utilities;',...
        FSDApath '/examples;' FSDApath '/robust;' FSDApath '/combinatorial;', ...
        FSDApath '/clustering;' FSDApath '/datasets/clustering;' FSDApath '/FSDAdemos'];
end
path(addp,path);

% Save current path for future sessions
savepath
end
