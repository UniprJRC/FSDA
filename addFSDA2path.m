function addFSDA2path(FSroot)
%Add FSDA toolbox to path
%
%  Required input arguments:
%
%  Optional input arguments:
%
%    FSroot:    path of the folder which contains FSDA toolbox. String.
%               A string containing the path which contains the root folder
%               of FSDA toolbox.
%               Example - 'D:\MATLAB\FSDA'   
%               Data Types - char
%               
%  Output:  
%           
%   Add a series of folders to the MATLAB search path
%
% More About:
%
% Function addFSDA2path adds to path the following subfolders of the main
% root of FSDA:
% \regression;
% \multivariate;
% \regression;
% \clustering;
% \graphics; 
% \datasets\regression;
% \datasets\multivariate;
% \datasets\multivariate_regression;
% \datasets\clustering
% \combinatorial;
% \utilities;
% \utilities_stat;
% \examples;
% \FSDAdemos
%
% In order to check that the previous folders have been added to path click
% on Home|Set Path
%
% See also path
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

if nargin<1
    FSroot= fileparts(which('docsearchFS.m'));
end

f=filesep;


    addp=[FSroot ';' FSroot f 'multivariate;' FSroot f 'regression;',...
        FSroot f 'clustering;' FSroot f 'graphics;',... 
        FSroot f 'datasets' f 'regression;' FSroot f 'datasets' f 'multivariate;', ...
        FSroot f 'datasets' f 'multivariate_regression;' FSroot f 'datasets' f 'clustering;'...
        FSroot f 'combinatorial;', ...
        FSroot f 'utilities;' FSroot f 'utilities_stat;' ...
        FSroot f 'examples;' FSroot f 'FSDAdemos'];

    path(addp,path);

% disp('REMARK: Remember to save the added folders to path in the MATLAB window')
% disp(['In the menu ' '''Home|Set path''' ' click on the button '  '''Save'''])
% Save current path for future sessions
savepath
%    

end
