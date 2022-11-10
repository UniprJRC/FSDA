function openMatlabFileFromHTML(MatlabFileName,LineToOpen)
% openMatlabFileFromHTML enables to put in HTML an hypertextual link to a specific MATLAB file
%
% <a href="matlab: docsearchFS('openMatlabFileFromHTML')">Link to the help function</a>
%
%
% Sometimes it is necessary to write HTML pages which contain link to
% MATLAB scripts or functions which have to be opened on user clik. This
% function enables to open the MATLAB script or function. For example, if
% when writing an HTML page, the user wants to add a link to MATLAB file
% named examples_categorical.m, the code to incorporate inside the HTML
% file is:
% <a href="matlab:linkToMatlabFileFromHTML('examples_categorical.m',1)">examples_categorical.m</a>
%
%
% Required input arguments:
%
%   MatlabFileName:   Name of Matlab file. Character. 
%               Filename  of the MATLAB function or script which has to be
%               opened when the user clicks on the hypertextual link. If it
%               is known whether the file to open is either a .m or a
%               .mlx file, please include the extension in the filename.
%
% Optional input arguments:
%
%   LineToOpen: Number of the line of the MATLAB file. Scalar.
%               Scalar which specifies the line where to position the
%               cursor when opening the file. 
%                 Example - 3
%
% Output:
%
%
% More About:
% 
% This function uses the undocumented MATLAB function opentoline. Finally
% note that the MATLAB file to open must be in MATLAB path because we use
% function 'which'.
%
% See also: findDir
%
% References:
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('openMatlabFileFromHTML')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    
    % openMatlabFileFromHTML with all default options.
    % In this case file examples_regression.m is opened.
    % editor is not available in terminal mode (i.e. Azure servers)
    if matlab.desktop.editor.isEditorAvailable ==true
        openMatlabFileFromHTML('examples_regression.m')
    end
%}
%

%{
    % openMatlabFileFromHTML in the command window.
    % If the user clicks on the link file examples_regression is opened on
    % line five.
    disp(['<a href="matlab:openMatlabFileFromHTML(''examples_regression.m'',5)">examples_regression.m</a>'])
%}

%% Beginning of code

% Position the cursor on line 1 if not specified otherwise
if nargin ==1
    LineToOpen=1;
end

% Find file on the MATLAB path
examp=which(MatlabFileName);

examp1=strrep(examp,filesep,[filesep, filesep]);

% Open the MATLAB file at a particular line
opentoline(examp1,LineToOpen)

end
%FScategory:UTIGEN