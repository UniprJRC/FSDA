function out=pdfprotect(inputfile, varargin)
%pdfprotect protects pdf files against printing and copying content of a pdf file
%
%<a href="matlab: docsearchFS('pdfprotect')">Link to the help page for this function</a>
%
%   This function protect pdf files against printing and copying content.
%   Also, when needed, can add a watermark diagonally on all pages of the manuscript.
%   Please note that the output pdf file will be encrypted, so it is better to
%   save the original document for backup purposes.
%   pdfprotect assumes that the user installed miniconda (a Python distribution)
%   with all the default options. If this is not the case, the user should
%   edit the path to Python executable inside MATLAB code of this function
%   according to the custom setup.
%
%
% Required input arguments:
%
% inputfile:    Input pdf file. Character. Input pdf file specified as a
%               char vector is the original file to be encrypted protected.
%
%
%  Optional input arguments:
%
% watermark:  text of the watermark. Character. The text that will be printed
%           diagonally in grey and in a big size font (85 points) on
%           each page of the manuscript.
%           Example - 'watermark', '(C)FSDA toolbox'
%           Data Types - char
%
% outputfile: name of the outputfile. Character. The name of the pdf file to be
%           encrypted, password protected and permission edited.
%           default vslue is the inpufileENC.pdf that is we append the suffix ENC
%           to the inputfilename
%           Example - 'outputfile', 'myoutputfile.pdf'
%           Data Types - char
%
% print:    print permission flag. Boolean. This flag controls the user permission to print the
%           document, can be either true or false.
%           Example - 'print', true
%           Data Types - boolean
%
% edit:     edit permission flag. Boolean. This flag controls the user permission to edit and copy
%           any content of the document, can be either true or false .
%           Example - 'edit', false
%           Data Types - boolean
%
%
% passedit: owner encryption password. Character. This argument should
%           contain the pdf owner password for encrypting the pdf and access the
%           permission flags.
%           The default password is 'FSDA'.
%           Example - 'passedit', 'MyPassword1'
%           Data Types - char
%
% passread: user password. Character. This argument should
%           contain the pdf owner password for encrypting the pdf and access the
%           permission flags.
%           The default password is not set.
%           Example - 'passread', 'MyPassword2'
%           Data Types - char


%
% Output:
%
%         out:  return status. Scalar. Returns the status of calling the python pdf protectiong
%         function, -1 indicates failure.
%
% See also:
%
% References:
%
%
% Official documentation concerning the PDF 1.7 specifications,
% https://opensource.adobe.com/dc-acrobat-sdk-docs/pdfstandards/pdfreference1.7old.pdf
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('pdfprotect')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2024-10-28 00:27:12 #$: Date of the last commit

% Examples:

%{
   % Example with all default values.
   % encrypt and protect a file revoking all authorizations to copy and
   % print.
    pdfprotect('sourcefile.pdf')
%}

%{

   % Example with watermark option.
   % encrypt and protect a file revoking all authorizations to copy and
   % print and add a custom watermark to each page of the manuscript.
    pdfprotect('sourcefile','watermark','(C)FSDAToolbox')
%}


%{
    % Example with watermark and print options.
    % encrypt and protect a file revoking copy and paste authorizations but
    % allow manuscript printing and add a custom watermark to each page of 
    % the manuscript.
    pdfprotect('sourcefile','watermark','(C) FSDA Toolbox','print_flag',true)
%}

%{
    % Example with all the available options. encrypt and protect a file
    % revoking copy, paste and print authorizations add a custom watermark
    % to each page of the manuscript and specify the edit password and read
    % password as well.
    pdfprotect('sourcefile','watermark','FSDAToolbox','print',false, 'edit',false,'outputfile', ...
       'mypdf', 'passedit','bigsecret', 'passread','easyguess')

%}



%% Beginning of code

if nargin < 1
    error('FSDA:pdfprotect:missingInputs','A required input argument (input file name) is missing.')
end





% allows use of strings as arguments
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end



% add pdf extension to inputfile if needed
FileName=(inputfile);
[pathstrcf,name,ext]=fileparts(FileName);

if isempty(ext) & isempty(pathstrcf)
    % the input pdf file is in the same folder
    inputfile = [inputfile '.pdf'];
elseif isempty(ext) & ~isempty(pathstrcf)
    % the input pdf file is in a different folder
    inputfile = [pathstrcf filesep inputfile '.pdf'];
else
    inputfile = [name '.' ext];
end


% check if inputfile exists 
if exist(inputfile, 'file') == 0
    error('FSDA:pdfprotect:missingInputFile','A pdf input file should be specified.')
end


% default parameters values
watermark = 'F';
outputfile = 'outputfile.pdf';
print = 'F';
edit = 'F';
passedit = 'FSDA';
passread = 'F';





UserOptions=varargin(1:2:length(varargin));

if ~isempty(UserOptions)

    options = struct('watermark',watermark, 'outputfile', outputfile, ...
        'print', print, 'edit',edit, 'passedit', passedit, 'passread', passread);


    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:pdfprotect:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    aux.chkoptions(options,UserOptions)

    % We now overwrite inside structure options the default values with
    % those chosen by the user
    % Notice that in order to do this we use dynamic field names
    for j=1:2:length(varargin)
        options.(varargin{j})=varargin{j+1};
    end

    % sanitize the watermark text
    watermark=matlab.lang.makeValidName(options.watermark);
    outputfile=options.outputfile;


    % add pdf extension to outputfile if needed
    FileName=(outputfile);
    [~,name,ext]=fileparts(FileName);

    if isempty(ext) 
        % the input pdf file is in the same folder
        outputfile = [name '.pdf'];
    end


    if options.print == true
        print='T';
    else
        print='F';
    end

    if options.edit==true
        edit='T';
    else
        edit='F';
    end
    passedit=options.passedit;
    passread=options.passread;
end


% space
sp=' ';

if ismac
      % linux: TODO!
    disp('Sorry, MacOS version is not available at the moment....')
    out=-1;
    return

    % get the name of the MacOS current USER
    [~,curruser]=system('id -un');
    pythonpath=['/Users/' curruser '/miniconda3/bin/'];
    % compose the string
    str=[ pythonpath '/python pdf_encryption_wm_creation.py ' inputfile sp ...
        watermark sp outputfile sp print_flag sp edit_flag sp password_text];

elseif ispc
    % get the path to python
    pythonpath = fullfile(getenv('USERPROFILE'), 'miniconda3');
    pythoncode = which('pdfprotect.m');
    [pythoncode1]=fileparts(pythoncode);
    pythoncode2=[pythoncode1 filesep 'private' filesep 'pdf_encryption_wm_creation.py'];
    % compose the string
    str=[ pythonpath '\python ' pythoncode2 sp inputfile sp ...
        watermark sp outputfile sp print sp edit sp passedit sp passread];
else
    % linux: TODO!
    disp('Sorry, Linux version is not available at the moment....')
    out=-1;
    return
end
disp(str)
[status,~]=system(str);

if status ~=0
    error('FSDA:pdfprotect:callerror','The call to python function was unsuccessful, please check miniconda setup ans python path.');
end

% delete watermark temporary file
delete watermark_layer.pdf

% returns the status of calling the python pdf protectiong function
out=status;
end