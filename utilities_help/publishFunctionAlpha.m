function fstring=publishFunctionAlpha(InputCell, varargin)
%publishFunctionAlpha enables to create web page which contains the alphabetical list of functions
%
% This routins uses as input the cell which is created with routine
% makecontentsfileFS.m and uses template file function-alphaEmpty.html
% which is contained inside:
% (main root of FSDA) filesep 'helpfiles' filesep 'FSDA;
% to create in a fully automatic way the alphabetical list of
% functions with automatic links for each alphabetic letter.
% The output file will be created inside
% (main root of FSDA) filesep 'helpfiles' filesep 'FSDA;
% and will have name function-alphaEmpty.html
% The automatic help file starting from structured .m file can be created
% using function publishFS.m
%
%
%<a href="matlab: docsearchFS('publishFunctionAlpha')">Link to the help function</a>
%
% Required input arguments:
%
%   InputCell: Cell created by function makecontentsfileFS.m. Cell. Cell
%              containing information about all files which have to be
%              included inside the alphabetical HTML file.
%
% Optional input arguments:
%
%     CreateTxtFile : create txt file. Boolean. If CreateTxtFile is false
%                     (default) no txt file is created, else, a txt file named
%                     function-alpha.txt is created which contains, the
%                     names of the files (separated by commas), inside
%                     folder (main root of FSDA)\helpfiles\FSDA. Every
%                     function which will be automatically created by
%                     publishFS looks into this file and (through a
%                     Javascript) adds the link on top and at the bottom of
%                     the page to the two file which come before and after
%                     current file in alphabetical order.
%                 Example - 'CreateTxtFile',false
%                 Data Types - boolean
%
% webhelp :   Option which specifies the default path to create html file
%             containing the alphabetical list of functions.
%             Logical.
%             If webhelp is true, the output is produced in the path
%             (FSDA root)\helpfiles\FSDAweb.
%             If webhelp is false (default), the output is produced in the path
%             (FSDA root)\helpfiles\FSDA.
%             Note that this option is valid just if outpuDir option below
%             is omitted.
%             Example - 'webhelp',true
%             Data Types - logical
%
% outputDir : Output folder. String.
%             Output folder to which the HTML document is saved and where
%             template file function-alphaEmpty.html is located, specified
%             as the comma-separated pair consisting of 'outputDir' and the
%             full path. You must specify the full path as a string, for
%             example 'C:\PublishedOutput'. Note that inside outputDir
%             there must be a file named "function-alphaEmpty.html" which
%             contains the template to create the alphabetical list of
%             functions. The defaults of 'outputDir' are as follows:
%             if input option webhelp is false  outputDir is
%             (FSDA root)\helpfiles\FSDA path,
%             else if input option webhelp is true  outputDir is
%             (FSDA root)\helpfiles\FSDAweb path.
%             Remark - outputDir must be a valid path.
%             Example - 'outputDir','C:'
%             Data Types - string
%
% Createnavbars : create files topscript.js and bottoscript.js. Boolean.
%             If Createnavbars is true (default), this routine assumes
%             there is a templsate input files inside (FSDA
%             root)\helpfiles\FSDA\includesFS named topscriptEMPTY.js and
%             bottomscriptEMPTY.js which contain top navigation bar to
%             include in each .html file. When this option is set to true,
%             file (FSDA root)\helpfiles\FSDA\includesFS\top(bottom)script.js
%             are created. This files contain an array variable that is the
%             list of files in input option InputCell and enables
%             the navigation bar functionality.
%             Example - 'Createnavbars',false
%             Data Types - boolean
%
% Output:
%
%       fstring:  string containing list of files in alphabetical order.
%                String. This string contains the full HTML files which all
%                hypertextual links to all HTML files for each alphabetical
%                letter. The HTML file function-alpha.html also produced inside
%                folder (main root of FSDA)\helpfiles\FSDA
%
% See also:    publishFunctionCate.m, publishFS.m
%
% References:
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('publishFunctionAlpha')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%
%{
    % Creation of file containing alphabetical list of functions.
    % Create contents file for each .m file
    % findDir with optional arguments 'InclDir' and 'ExclDir'.
    FileName='addFSDA2path';
    FullPath=which(FileName);
    FSDAroot=FullPath(1:end-length(FileName)-3);
    InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat'};
    ExclDir={'privateFS'  'datasets'};
    list = findDir(FSDAroot,'InclDir',InclDir,'ExclDir',ExclDir)
    out=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory','force',false);
    cd(fileparts(which('docsearchFS.m')))
    % Create HTML file containing alphabetical list of functions
    fstring=publishFunctionAlpha(out);
    % open outfile file in web browser
    FileWithFullPath=which('docsearchFS.m');
    [pathFSDAstr]=fileparts(FileWithFullPath);
    fsep=filesep;
    outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA\function-alpha.html'];
    web(outputOFHtmlHelpFile,'-browser');
%}

%{
    % Creation of txt file.
    % File function-alpha.txt is created which contains, the
    % names of the files (separated by commas), inside
    % folder (main root of FSDA)\helpfiles\FSDA.
    % Create contents file for each .m file
    % findDir with optional arguments 'InclDir' and 'ExclDir'.
    FileName='addFSDA2path';
    FullPath=which(FileName);
    FSDAroot=FullPath(1:end-length(FileName)-3);
    InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat'};
    ExclDir={'privateFS'  'datasets'};
    list = findDir(FSDAroot,'InclDir',InclDir,'ExclDir',ExclDir)
    out=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory','force',false);
    cd(fileparts(which('docsearchFS.m')))
    % Create HTML file containing alphabetical list of functions
    fstring=publishFunctionAlpha(out,'CreateTxtFile',true);
    % open outfile txt in web browser
    FileWithFullPath=which('docsearchFS.m');
    [pathFSDAstr]=fileparts(FileWithFullPath);
    fsep=filesep;
    outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA\function-alpha.txt'];
    web(outputOFHtmlHelpFile,'-browser');
%}


%% Beginning of code

% % Use file separator of current operating system
% % \ = Windows
% % / = Unix
fsep=filesep;
CreateTxtFile=false;
Createnavbars=true;

if nargin>1
    UserOptions=varargin(1:2:length(varargin));
    checklms2 = strcmp(UserOptions,'webhelp');
    if sum(checklms2)
        webhelp = varargin{2*find(checklms2)};
    else
        webhelp=false;
    end
else
    webhelp=false;
end

FileWithFullPath=which('docsearchFS.m');
[pathFSDAstr]=fileparts(FileWithFullPath);

if webhelp == false
    outputDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA'];
else
    outputDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDAweb'];
end


if nargin>1
    options=struct('CreateTxtFile',CreateTxtFile,'outputDir',...
        outputDir,'webhelp',webhelp,'Createnavbars', Createnavbars);
    
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:publishFunctionAlpha:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    CreateTxtFile=options.CreateTxtFile;
    outputDir=options.outputDir;
    Createnavbars=options.Createnavbars;
end


% Open input function-alphaEmpty.html file, put it in a string and do a series of preliminary operations
FileWithFullPath=[outputDir fsep 'function-alphaEmpty.html'];

fileID = fopen(char(FileWithFullPath), 'r');

if fileID==-1
    disp(['Output path: '''  char(outputDir) ''' must  contain a file named  ''function-alphaEmpty.html'''])
    error('FSDA:publishFunctionAlpha:WrongPath','Output path does not have the input file')
end


fileID = fopen(char(FileWithFullPath), 'r');

% Insert the file into fstring
fstring=fscanf(fileID,'%c');

seqAZ=('A':'Z')';
for i=1:length(seqAZ)
    letterINI=['--BEGOF-' seqAZ(i) '--'];
    if i<length(seqAZ)
        letterFIN=['--BEGOF-' seqAZ(i+1) '--'];
    else
        letterFIN='END OF FSDA TEXT';
    end
    % find position of letter inside empty file
    ini=regexp(fstring,letterINI);
    fin=regexp(fstring,letterFIN);
    
    iniHTMLTEXT=regexp(fstring,'HTMLTEXT');
    iniHTMLTEXT=iniHTMLTEXT(iniHTMLTEXT>ini);
    
    iniHTMLTEXT=iniHTMLTEXT(iniHTMLTEXT<fin);
    
    % Select all lines of InputCell starting with letter seqAZ(i)
    if strcmp(seqAZ(i),'M')
        boo=cellfun('length',regexpi(InputCell(:,1),['\<' seqAZ(i) '\w*']))==1;
    else
        boo=cellfun('isempty',regexpi(InputCell(:,1),['\<' seqAZ(i) '\w*']));
    end
    Inputi=InputCell(~boo,:);
    
    strInsert='';
    
    if ~isempty(Inputi)
        for j=1:size(Inputi,1)
            mfilename=Inputi{j,6};
            description=Inputi{j,7};
            strInsert= sprintf([strInsert '<tr>\r' ...
                '<td class="term"><a href="' mfilename '.html">' mfilename '</a></td>\r' ...
                '<td class="description">' description '</td>\r' ...
                '</tr>\r']);
        end
    end
    fstring=[fstring(1:iniHTMLTEXT-1) strInsert fstring(iniHTMLTEXT+8:end)];
end

% Write output file

% Print fstring in new HTML file named function-alpha.html
name='function-alpha';
file1ID=fopen([outputDir fsep name '.html'],'w');

if file1ID==-1
    
    if ismac || isunix
        errmsg= [' Path ' outputDir '/' name '.html does not exist or output file '  name '.html is not writable'];
    elseif ispc
        outputDir=strrep(outputDir,'\','\\');
        errmsg= [' Path ' outputDir '\\' name '.html does not exist or output file '  name '.html is not writable'];
    else
        errmsg= [' Path ' outputDir '/' name '.html does not exist or output file '  name '.html is not writable'];
    end
    
    error('FSDA:publishFunctionAlpha:WrngOutFolder',errmsg);
end
fprintf(file1ID,'%s',fstring);

if CreateTxtFile
    %now write inside file function-alpha.txt if requested
    file2ID=fopen([outputDir fsep name '.txt'],'w');
    
    if file2ID==-1
        
        if ismac || isunix
            errmsg= [' Path ' outputDir '/' name '.html does not exist or output file '  name '.html is not writable'];
        elseif ispc
            outputDir=strrep(outputDir,'\','\\');
            errmsg= [' Path ' outputDir '\\' name '.html does not exist or output file '  name '.html is not writable'];
        else
            errmsg= [' Path ' outputDir '/' name '.html does not exist or output file '  name '.html is not writable'];
        end
        
        error('FSDA:publishFunctionAlpha:WrngOutFolder',errmsg);
    end
    
    % now construct string which has to be written inside the file
    % Initialize the string.
    ReQstring='index.html';
    
    for i=1:size(InputCell,1)
        ReQstring =[ ReQstring ',' InputCell{i,6} '.html']; %#ok<AGROW>
    end
    ReQstring=[ ReQstring ',' 'function-cate.html'];
    
    
    fprintf(file2ID,'%s',ReQstring);
    
end

%% Navbar creation

if Createnavbars == true
    
    % first create topscript.js
    
    % Open input function-alphaEmpty.html file, put it in a string and do a series of preliminary operations
    FileTopBottomWithPath=[outputDir fsep 'includesFS' fsep 'topscriptEmpty.js'];
    
    %now write inside file function-alpha.txt if requested
    file3ID=fopen(FileTopBottomWithPath,'r');
    
    if file3ID==-1
        
        errmsg= ['Input file' FileTopBottomWithPath  'does not exist'];
        error('FSDA:publishFunctionAlpha:WrngOutInputFiler',errmsg);
    end
    
    % Insert the content of file topscriptEmpty into fstring
    fstring=fscanf(file3ID,'%c');
    
    strInsertTop='"index.html';
    
    for i=1:size(InputCell,1)
        strInsertTop =[ strInsertTop '","' InputCell{i,6} '.html']; %#ok<AGROW>
    end
    strInsertTop=[ strInsertTop '","' 'function-cate.html"'];
    
    iniHTMLTEXT=regexp(fstring,'LIST_OF_FILES');
    
    fstring=[fstring(1:iniHTMLTEXT-1) strInsertTop fstring(iniHTMLTEXT+13:end)];
    
    % Open file topscript.js for writing
    FileTopBottomWithPath=[outputDir fsep 'includesFS' fsep 'topscript.js'];
    
    file4ID=fopen(FileTopBottomWithPath,'w');
    % Write the string fstring into topscript.js
    fprintf(file4ID,'%s',fstring);
    
    % then create bottomscript.js
    
    % Open input function-alphaEmpty.html file, put it in a string and do a series of preliminary operations
    FileTopBottomWithPath=[outputDir fsep 'includesFS' fsep 'bottomscriptEmpty.js'];
    
    %now write inside file function-alpha.txt if requested
    file3ID=fopen(FileTopBottomWithPath,'r');
    
    if file3ID==-1
        
        errmsg= ['Input file' FileTopBottomWithPath  'does not exist'];
        error('FSDA:publishFunctionAlpha:WrngOutInputFiler',errmsg);
    end
    
    % Insert the content of file topscriptEmpty into fstring
    fstring=fscanf(file3ID,'%c');
    
    strInsertTop='"index.html';
    
    for i=1:size(InputCell,1)
        strInsertTop =[ strInsertTop '","' InputCell{i,6} '.html']; %#ok<AGROW>
    end
    strInsertTop=[ strInsertTop '","' 'function-cate.html"'];
    
    iniHTMLTEXT=regexp(fstring,'LIST_OF_FILES');
    
    fstring=[fstring(1:iniHTMLTEXT-1) strInsertTop fstring(iniHTMLTEXT+13:end)];
    
    % Open file topscript.js for writing
    FileTopBottomWithPath=[outputDir fsep 'includesFS' fsep 'bottomscript.js'];
    
    file4ID=fopen(FileTopBottomWithPath,'w');
    % Write the string fstring into topscript.js
    fprintf(file4ID,'%s',fstring);
    
end

fclose('all');
end
%FScategory:UTIHELP
