function fstring=publishFunctionCate(InputCell)
%publishFunctionCate enables to create web page which contains the categorical list of functions
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
%<a href="matlab: docsearchFS('publishFunctionCate')">Link to the help function</a>
%
% Required input arguments:
%
%   InputCell: Cell created by function makecontentsfileFS.m. Cell. Cell
%              containing information about all files which have to be
%              included inside the categorical HTML file.
%
% Optional input arguments:
%
%
% Output:
%
%       fstring:  string containing list of files in categorical order.
%                String. This string contains the full HTML files which all
%                hypertextual links to all HTML files for each category.
%                The HTML file function-cate.html also produced inside 
%                folder (main root of FSDA)\helpfiles\FSDA
%
% See also:    publishfunctionAlpha.m, publishFS.m
%
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('publishFunctionCate')">Link to the help function</a>
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
    out=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory','force',false);
    cd(fileparts(which('docsearchFS.m')))
    % Create HTML file containing categorical list of functions
    fstring=publishFunctionCate(out);
    % open outfile file in web browser
    FileWithFullPath=which('docsearchFS.m');
    [pathFSDAstr]=fileparts(FileWithFullPath);
    fsep=filesep;
    outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA\function-cate.html'];
    web(outputOFHtmlHelpFile,'-browser');
%}

%% Beginning of code

% % Use file separator of current operating system
% % \ = Windows
% % / = Unix
fsep=filesep;

% Write output file in subfolder \(FSDAroot)\helpfiles\FSDA
FileWithFullPath=which('docsearchFS.m');
[pathFSDAstr]=fileparts(FileWithFullPath);

outputDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA\'];

%% Open input function-cateEmpty.html file, put it in a string and do a series of preliminary operations

FileWithFullPath=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA\function-cateEmpty.html'];
fileID = fopen(char(FileWithFullPath), 'r');

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
    if strcmp(letterINI,'GUI')
        ini=ini(1);
    else
    assert(length(ini)==1,['Category: ' letterINI 'is duplicated'])
    end
    
    fin=regexp(fstring,letterFIN);
    fin=fin(fin>ini);
    fin=fin(1);
    
    iniHTMLTEXT=regexp(fstring,'HTMLTEXT');
    iniHTMLTEXT=iniHTMLTEXT(iniHTMLTEXT>ini);
    
    iniHTMLTEXT=iniHTMLTEXT(iniHTMLTEXT<fin);
    
    % Select all lines of InputCell starting with letter seqAZ(i)
    %     if strcmp(seqCAT(i),'M')
    %         boo=cellfun('length',regexpi(InputCell(:,1),['\<' seqCAT(i) '\w*']))==1;
    %     else
    %         boo=cellfun('isempty',regexpi(InputCell(:,1),['\<' seqCAT(i) '\w*']));
    %     end
    
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
