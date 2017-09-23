function linkToMatlabFileFromHTML(MatlabFileName,LineToOpen)
% LinkToMatlabFileFromHTML enables to put in HTML an hypertextual link to a specific MATLAB file
%
% <a href="matlab: docsearchFS('LinkToMatlabFileFromHTML')">Link to the help function</a>
%
%
% Required input arguments:
%
%   MatlabFileName:   Name of Matlab file. Character. 
%               String which contains the name of the MATLAB file which has
%               to be opened when user clicks on the hypertextual link.
%               Note that to avoid confusion it is better to write the name
%               with the extension. If the name is specified without
%               extension and there are both the .m and the .mlx file, the
%               .mlx file is opened. 
%               For example if when writing an HTML page the user wants to
%               add a link to MATLAB file named  examples_categorical.m
%               here is the code he has to write:
%               <a href="matlab:linkToMatlabFileFromHTML('examples_categorical.m',1)">examples_categorical.m</a>
%
% Optional input arguments:
%
%   LineToOpen: Number of line of the MATLAB file. Scalar.
%               Scalar which specifies the line number of the -m file which
%               has to be opened.
%
% Output:
%
%
% See also: findDir
%
% References:
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('LinkToMatlabFileFromHTML')">Link to the help page for this function</a>
% Last modified 31-05-2016
%
% Examples:
%
%{
    
    % linkToMatlabFileFromHTML with all default options.
    % In this case file examples_regression.m is opened.
    linkToMatlabFileFromHTML('examples_regression.m')
%}
%

%{
    % linkToMatlabFileFromHTML in the command window.
    % If the user clicks on the link file examples_regression is opened on
    % line five.
    disp(['<a href="matlab:linkToMatlabFileFromHTML(''examples_regression.m'',5)">examples_regression.m</a>'])
%}

%% Beginning of code

% Write in structure 'options' the options chosen by the user
if nargin ==1
    LineToOpen=1;
end

% Find file on the MATLAB path
examp=which(MatlabFileName);
examp1=strrep(examp,'\','\\');

% Open the MATLAB file at a particular line
opentoline(examp1,LineToOpen)

end
%FScategory:UTIGEN