function out=pdfprotect(inputfile, varargin)
%pdfprotect protects pdf files against printing and copying content of a pdf file
%
%<a href="matlab: docsearchFS('pdfprotect')">Link to the help page for this function</a>
%
%   This function protects pdf files against printing and copying content.
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
% inputfile:    Input pdf file. String scalar | Character vector. 
%               Input pdf file specified as a
%               char vector is the original file to be encrypted protected.
%               Can be specified with or without extension.
%           Data Types - char or string
%
%
%  Optional input arguments:
%
%
% edit:     edit permission flag. Boolean. This flag controls the user
%           permission to edit and copy any content of the document, can be
%           either true or false. The default value of edit is false, that
%           is you cannot edit the file. 
%           Example - 'edit', false
%           Data Types - boolean
%
% outputfile: name of the outputfile. String scalar | Character vector. 
%           Name of the pdf file to be encrypted, password protected and
%           permission edited. Can be specified with or without extension.
%           Default value is inpufileENC.pdf that is we append the
%           suffix ENC to the inputfilename 
%           Example - 'outputfile','myoutputfile.pdf' 
%           Data Types - char or string
%
% passedit: owner encryption password. String scalar | Character vector. If the password is set you need to supply the
%           password to modify and print the document.
%           The default password is 'FSDA'.
%           Example - 'passedit', 'MyPassword1'
%           Data Types - char or string
%
% passopen: user password. String scalar | Character vector. If the password is set you need to supply the
%           password to view the document.
%           If this option is not specified 
%           the password to open the file is not set.
%           Example - 'passread', 'MyPassword2'
%           Data Types - char or string
%
% print:    print permission flag. Boolean. 
%           This flag controls the user permission to print the
%           document, can be either true (allow) or false (deny).
%           The default value of print is false.            
%           Example - 'print', true
%           Data Types - boolean
%
%
% watermark:  text of the watermark. String scalar | Character vector. 
%           The text that will be printed diagonally in grey and in a big
%           size font (85 points) on each page of the manuscript. The
%           default value fo the watermark is '' that is no watermark is
%           added.
%           Example - 'watermark', '(C)FSDA toolbox'
%           Data Types - char or string
%
% Output:
%
%         out:  return status. Scalar. Returns the status of calling the python pdf protectiong
%               function, -1 indicates failure.
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
   % Create pdf file tmptmpENC.pdf 
   % create .m file==> transform to mlx ==> export to pdf
   FileName='addFSDA2path.m';
   FileNameMLX='tmptmp.mlx';
   matlab.internal.liveeditor.openAndSave(which(FileName),FileNameMLX);
   % Create file tmptmp.pdf
   export(FileNameMLX);
   % Extract file name without extension
   [~,NameWithoutExtension]=fileparts(FileNameMLX);
   % delete temporary mlx file
   delete(FileNameMLX);

   % Call to pdfprotect with all default options
   % encrypt and protect a file revoking all authorizations to copy and
   % print.
   pdfprotect(NameWithoutExtension);
   
   disp('Files named "tmptmpENC.pdf" has been created')
   disp('in the current folder. In this file editing and printing is disabled')
   disp('Original input file tmptmp.pdf on the other hand is not protected')
%}

%{
   %% Example with watermark option.
   % Encrypt and protect a file revoking all authorizations to copy and
   % print and add a custom watermark to each page of the manuscript.

   % Create pdf file tmptmpENC.pdf 
   % create .m file==> transform to mlx ==> export to pdf
   FileName='addFSDA2path.m';
   FileNameMLX='tmptmp.mlx';
   matlab.internal.liveeditor.openAndSave(which(FileName),FileNameMLX);
   % Create file tmptmp.pdf
   export(FileNameMLX);
   % Extract file name without extension
   [~,NameWithoutExtension]=fileparts(FileNameMLX);
   % delete temporary mlx file
   delete(FileNameMLX);

   % Call to pdfprotect with all default options
   % encrypt and protect a file revoking all authorizations to copy and
   % print.
   pdfprotect(NameWithoutExtension,'watermark','secret');
   
   disp('Files named "tmptmpENC.pdf" has been created')
   disp('in the current folder. In this file editing and printing is disabled')
   disp('and a watermark with text "Secret" has been added')
   disp('Original input file tmptmp.pdf on the other hand is not protected')

%}


%{
   %% Example with watermark and print option.
   % Encrypt and protect a file revoking copy and paste authorizations but
   % allow manuscript printing and add a custom watermark to each page of 
   % the manuscript.
   % Create pdf file tmptmpENC.pdf 
   % create .m file==> transform to mlx ==> export to pdf
   FileName='addFSDA2path.m';
   FileNameMLX='tmptmp.mlx';
   matlab.internal.liveeditor.openAndSave(which(FileName),FileNameMLX);
   % Create file tmptmp.pdf
   export(FileNameMLX);
   % Extract file name without extension
   [~,NameWithoutExtension]=fileparts(FileNameMLX);
   % delete temporary mlx file
   delete(FileNameMLX);

   % Call to pdfprotect with all default options
   % encrypt and protect a file revoking all authorizations to copy and
   % print.
   pdfprotect(NameWithoutExtension,'watermark','FSDA_Toolbox','print',true);
  
   disp('Files named "tmptmpENC.pdf" has been created')
   disp('in the current folder. In this file editing is disabled')
   disp('printing is allowed and a watermark with text "FSDA_Toolbox" has been added')
   disp('Original input file tmptmp.pdf on the other hand is not protected')
%}

%{
   %% Example with personalized passwords and name of output file.
   % Encrypt and protect a file
   % revoking copy, paste and print authorizations add a custom watermark
   % to each page of the manuscript, specify the edit password and read
   % password as well. Also specify name of output file
   % Create pdf file tmptmpENC.pdf 
   % create .m file==> transform to mlx ==> export to pdf
   FileName='addFSDA2path.m';
   FileNameMLX='tmptmp.mlx';
   matlab.internal.liveeditor.openAndSave(which(FileName),FileNameMLX);
   % Create file tmptmp.pdf
   export(FileNameMLX);
   % Extract file name without extension
   [~,NameWithoutExtension]=fileparts(FileNameMLX);
   % delete temporary mlx file
   delete(FileNameMLX);

   % Call to pdfprotect with options passedit, 
   % passopen and outputfile
   pdfprotect(NameWithoutExtension,'outputfile', ...
       'mypdfENC', 'passedit','bigsecret', 'passopen','easyguess')
  
   disp('File named "mypdfENC.pdf" has been created')
   disp('in the current folder. In this file editing is disabled')
   disp('printing is disables. There is also a pwd to open file.')
   disp('Original input file tmptmp.pdf on the other hand is not protected')
%}



%% Beginning of code

if nargin < 1
    error('FSDA:pdfprotect:missingInputs','A required input argument (input file name) is missing.')
else
inputfile=convertStringsToChars(inputfile);
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
    inputfile = [name ext];
end


% check if inputfile exists 
if exist(inputfile, 'file') == 0
    error('FSDA:pdfprotect:missingInputFile','A pdf input file should be specified.')
end


% default parameters values
watermark = 'F';
outputfile = [name 'ENC.pdf'];
print = 'F';
edit = 'F';
passedit = 'FSDA';
passopen = 'F';



UserOptions=varargin(1:2:length(varargin));

if ~isempty(UserOptions)

    options = struct('watermark',watermark, 'outputfile', outputfile, ...
        'print', print, 'edit',edit, 'passedit', passedit, 'passopen', passopen);


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

    if isempty(options.watermark)
        watermark='F';
    else
    % sanitize the watermark text
    watermark=matlab.lang.makeValidName(options.watermark);
    end

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
    passopen=options.passopen;
end


% space
sp=' ';

if ismac
      % linux: TODO!
    disp('Sorry, MacOS version is not available at the moment....')
    out=-1;
    return

    % % get the name of the MacOS current USER
    % [~,curruser]=system('id -un');
    % pythonpath=['/Users/' curruser '/miniconda3/bin/'];
    % % compose the string
    % str=[ pythonpath '/python pdf_encryption_wm_creation.py ' inputfile sp ...
    %     watermark sp outputfile sp print_flag sp edit_flag sp password_text];

elseif ispc
    % get the path to python
    pythonpath = fullfile(getenv('USERPROFILE'), 'miniconda3');
    pythoncode = which('pdfprotect.m');
    [pythoncode1]=fileparts(pythoncode);
    pythoncode2=[pythoncode1 filesep 'private' filesep 'pdf_encryption_wm_creation.py'];
% check if output file is already open
    [fstatus, errmsg]=fopen(outputfile,"w");
    if fstatus<0
        disp('The output file is already open, please close it.')
        disp(errmsg)
        return
    else
        fclose(fstatus);
        % compose the string
        str=[ pythonpath '\python ' pythoncode2 sp inputfile sp ...
            watermark sp outputfile sp print sp edit sp passedit sp passopen];
    end
else
    % linux: TODO!
    disp('Sorry, Linux version is not available at the moment....')
    out=-1;
    return
end
% disp(str)
[status,~]=system(str);

if status ~=0
    error('FSDA:pdfprotect:callerror','The call to python function was unsuccessful, please check miniconda setup ans python path.');
end

% delete watermark temporary file
delete watermark_layer.pdf

% returns the status of calling the python pdf protectiong function
out=status;
end