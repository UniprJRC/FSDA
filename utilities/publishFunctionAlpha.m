function fstring=publishFunctionAlpha(InputCell)
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
%
% Output:
%
%       fstring:  string containing list of files in alphabetical order.
%                String. This string contains the full HTML files which all
%                hypertextual links to all HTML files for each alphabetical
%                letter. The HTML file function-alpha.html also produced inside 
%                folder (main root of FSDA)\helpfiles\FSDA
%
% See also:    CreateHTLMfunctioncate.m, publishFS.m
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('CreateHTLMfunctionAlpha')">Link to the help function</a>
% Last modified 06-Feb-2015
%
% Examples:
%
% 
%{
    % Interactive_example
    % Creation of file containing alphabetical list of functions
    % make sure you are inside the main folder of FSDA.
    % create contents file for each .m file 
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

%% Beginning of code
% % Use file separator of current operating system
% % \ = Windows
% % / = Unix
fsep=filesep;


% Open input function-alphaEmpty.html file, put it in a string and do a series of preliminary operations
FileWithFullPath=which('docsearchFS.m');
pathFSDAstr=fileparts(FileWithFullPath);

FileWithFullPath=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA\function-alphaEmpty.html'];
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

% Write output file in subfolder \(FSDAroot)\helpfiles\FSDA

outputDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA\'];

% Pring fstring in new HTML file named function-alpha.html
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
    
    error('FSDA:publishFS:WrngOutFolder',errmsg);
end

fprintf(file1ID,'%s',fstring);
fclose('all');

end
%FScategory:UTIHELP
