function fstring=publishFunctionCate(InputCell, varargin)
%publishFunctionCate enables to create web page which contains the categorical list of functions
%
%<a href="matlab: docsearchFS('publishFunctionCate')">Link to the help function</a>
%
% This routins uses as input the cell which is created with routine
% makecontentsfileFS.m and uses template file function-cateEmpty.html
% which is contained inside:
% (main root of FSDA) filesep 'helpfiles' filesep 'FSDA;
% to create in a fully automatic way the categorical list of
% functions with automatic links for each category.
% The output file will be created inside
% (main root of FSDA) filesep 'helpfiles' filesep 'FSDA;
% and will have name function-cate.html
% The automatic help file starting from structured .m file can be created
% using function publishFS.m
%
%
%
% Required input arguments:
%
%   InputCell: Cell created by function makecontentsfileFS.m. Cell. Cell
%              containing information about all files which have to be
%              included inside the categorical list of functions HTML file.
%
% Optional input arguments:
%
% webhelp :   Option which specifies the default path to create html file
%             containing the categorical list of functions.
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
%             template file function-cateEmpty.html is located, specified
%             as the comma-separated pair consisting of 'outputDir' and the
%             full path. You must specify the full path as a string, for
%             example 'C:\PublishedOutput'. Note that inside outputDir
%             there must be a file named "function-cateEmpty.html" which
%             contains the template to create the categorical list of
%             functions. The defaults of 'outputDir' are as follows:
%             if input option webhelp is false  outputDir is 
%             (FSDA root)\helpfiles\FSDA path,
%             else if input option webhelp is true  outputDir is 
%             (FSDA root)\helpfiles\FSDAweb path.
%             Remark - outputDir must be a valid path.
%             Example - 'outputDir','C:'
%             Data Types - string
%
%
% Output:
%
%       fstring:  string containing list of files in categorical order.
%                String. This string contains the full HTML files which all
%                hypertextual links to all HTML files for each category.
%                The HTML file function-cate.html also produced inside
%                input option "outputDir" folder.
%
% See also:    publishfunctionAlpha.m, publishFS.m
%
%
%
% References:
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('publishFunctionCate')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%
%{
    % Interactive_example
    % Creation of HTML file containing categorical list of functions.
    % Make sure you are inside the main folder of FSDA.
    % Create contents file for each .m file
    % findDir with optional arguments 'InclDir' and 'ExclDir'.
    FileName='addFSDA2path';
    FullPath=which(FileName);
    FSDAroot=FullPath(1:end-length(FileName)-3);
    InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat'};
    ExclDir={'privateFS'  'datasets'};
    list = findDir(FSDAroot,'InclDir',InclDir,'ExclDir',ExclDir)
    out=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',false);
    cd(fileparts(which('docsearchFS.m')))
    % Create HTML file containing categorical list of functions
    fstring=publishFunctionCate(out);
    % open output file in web browser
    FileWithFullPath=which('docsearchFS.m');
    [pathFSDAstr]=fileparts(FileWithFullPath);
    fsep=filesep;
    outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA\function-cate.html'];
    web(outputOFHtmlHelpFile,'-browser');
%}


%{
    % Interactive_example
    % Create output file in personalized folder
    % Create HTML file containing categorical list of functions in 
    % local path "D:\tmp".
    % Note that we assume that inside path D:\tmp there is the template
    % file named "function-cateEmpty.html"
    fstring=publishFunctionCate(out,'outputDir','D:\tmp');
    
%}

%% Beginning of code
% % Use file separator of current operating system
% % \ = Windows
% % / = Unix
fsep=filesep;


if nargin>1
    [varargin{:}] = convertStringsToChars(varargin{:});
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
    options=struct('outputDir',outputDir,'webhelp',webhelp);
    
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:publishFunctionCate:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    outputDir=options.outputDir;
 %    webhelp=options.webhelp;
end


%% Open input function-cateEmpty.html file, put it in a string and do a series of preliminary operations


% Open input function-cateEmpty.html file, put it in a string and do a series of preliminary operations
FileWithFullPath=[outputDir fsep 'function-cateEmpty.html'];
fileID = fopen(char(FileWithFullPath), 'r');

if fileID==-1
    disp(['Output path: '''  char(outputDir) ''' must  contain a file named  ''function-cateEmpty.html'''])
    error('FSDA:publishFunctionCate:WrongPath','Output path does not have the input file')
end
        
% Insert the file into fstring
fstring=fscanf(fileID,'%c');


seqCAT=unique(InputCell(:,8));
for i=1:length(seqCAT)
    letterINI=seqCAT{i};
    
    letterFIN='</table>';
    
    disp(letterINI)
    % find position of letter inside empty file
    ini=regexp(fstring,letterINI);
    
    assert(~isempty(ini),['Category: ' letterINI ' not found'])
    % pick the 2nd instance of the Category anchor because the 1st is in
    % the navigation left panel


    if strcmp(letterINI,'MULT-Categorical')
        ini=ini(1);
    else
        ini=ini(2);
        assert(length(ini)==1,['Category: ' letterINI ' is duplicated'])
    end
    
    
    
    fin=regexp(fstring,letterFIN);
    fin=fin(fin>ini);
    fin=fin(1);
    
    iniHTMLTEXT=regexp(fstring,'HTMLTEXT');
    iniHTMLTEXT=iniHTMLTEXT(iniHTMLTEXT>ini);
    
    iniHTMLTEXT=iniHTMLTEXT(iniHTMLTEXT<fin);
    
    boo=strcmp(seqCAT(i),InputCell(:,8));
    Inputi=InputCell(boo,:);
    
    
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

name='function-cate';

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
    
    error('FSDA:publishFS:WrngOutFolder',errmsg);
end

fprintf(file1ID,'%s',fstring);
fclose('all');

end
%FScategory:UTIHELP
