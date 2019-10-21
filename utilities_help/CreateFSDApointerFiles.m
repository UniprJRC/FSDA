function outHELP=CreateFSDApointerFiles(InputCell,OUT)
%CreateFSDApointerFiles create HTML files which are just pointers to true HTML help files
%
% Given that:
% 1) buildocsearchdb enables to create the database just for
% the files which are inside the helpfile location specified in the
% info.xml file
% 2) all the entries which are found by buildocsearchdb are opened in the
% frame of the right (iframe)
% 3) FSDA team spent a lot of time to have proper HTML documentation which
% is not relagated to the right panel
% 4) true HTML documentation file of FSDA are inside (docroot)/FSDA otherwise
% the lucene MATLAB search engine does not work.
%
% The purpose of this function is to take as input all the files which are
% inside (FSDA root folder)/helpfiles/FSDA
% and for each of them to create a corresponding pointer HTML file inside path
% (FSDA root folder)/helpfiles/pointersHTML
% which contains
% a) the minimum necessary information to be indexed by lucene search engine
% b) a response.redirect which links to the true HTML file which is
% contained inside (docroot)/FSDA.
%
%
%{
    % Example of creation of pointer files.
    % It is necessary to call first makecontentsfileFS.m and
    % publishFSallFiles.m
    % 1) Create personalized contents file which will be the input the
    % procedure CreateFSDApointerFiles.
    % Find full path of the main root where FSDA is installed
    FileName='addFSDA2path';
    FullPath=which(FileName);
    root=FullPath(1:end-length(FileName)-3);
    InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat'};
    ExclDir={'privateFS'  'datasets'};
    % Create list of folders for which contents files has to be created
    list = findDir(root,'InclDir',InclDir,'ExclDir',ExclDir)
    % Create personalized contents file for main folder of FSDA
    % and required subfolders.
    [InputCell,Excluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',false)
    if verLessThan('matlab','8.1.0')==0
        % Publish all files
        [FilesWithProblems,OUT]=publishFSallFiles(InputCell,'write2file',false,'evalCode',false);
        % Create HTML pointer files
        CreateFSDApointerFiles(InputCell,OUT)
    else
        warning('At least MATLAB version 2013a is needed')
    end
%}

% Copyright 2008-2019.
% Written by FSDA team
%
%
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code
[FSDAroot]=fileparts(which('docsearchFS.m'));
fsep=filesep;
dirpathj=[FSDAroot fsep 'helpfiles' fsep 'pointersHTML'];

for i=1:size(InputCell,1)
    NameInputFile=InputCell{i,6};
    PurposeInputFile=InputCell{i,7};
    NameOutputFile=[NameInputFile '.html'];
    try
        DescrRef=OUT{i};
    catch
        dd=1;
    end
    [fileID,errMsg] = fopen([dirpathj filesep NameOutputFile],'w');
    if fileID < 0
        error(message('MATLAB:filebrowser:MakeContentsFileOpenError', errMsg))
    end
    % Create string to include the HTML pointer file
    try
        fstring=stringHTMLpointer(NameInputFile,PurposeInputFile,DescrRef);
    catch
        error('FSDA:CreateFSDApointerFiles:WrongInputFile',['Could not create pointer file: ' NameInputFile])
    end
    
    fprintf(fileID,'%s',fstring);
    fclose(fileID);
end
outHELP=true;


end

function fstring=stringHTMLpointer(NameInputFile,DescriptionInputFile,DescrRef)
NameInputFileHTML=['/FSDA/' NameInputFile '.html'];

Description=DescrRef.description;
MoreAbout=DescrRef.MoreAbout;
if isempty(DescrRef.References)
    References='';
else
    References=DescrRef.References{:};
end

try
   
    if ~isempty(Description)
        DescriptionHTML=['<h2>Description</h2>'  ...
            '<P>' Description '</P>'];
    else
        DescriptionHTML='';
    end
    
    if ~isempty(MoreAbout)
        MoreAboutHTML=[ sprintf('<h2>More About</h2>\r')  ...
            '<P>' MoreAbout '</P>'];
    else
        MoreAboutHTML='';
    end
    
    if ~isempty(References)
        ReferencesHTML=[sprintf('<h2>References</h2>\r')  ...
            '<P>' References '</P>'];
    else
        ReferencesHTML='';
    end
    
    fstring=[sprintf(['<!DOCTYPE HTML> \r'  ...
        '<html itemscope="" xmlns="http://www.w3.org/1999/xhtml"> \r'  ...
        '<head> \r'  ...
        '<title>' NameInputFile '</title> \r'  ...
        '<meta content="refpage" name="chunktype">\r'  ...
        '<meta content="function:' NameInputFile  ' " itemprop="refentity" name="refentity">\r'  ...
        '<meta content="fcn" itemprop="pagetype" name="toctype">\r'  ...
        '<meta content="ref/function" itemprop="infotype" name="infotype" />\r'  ...
        '<meta content="' NameInputFile ' ' DescriptionInputFile '" itemprop="description" name="description" />\r'  ...
        '<h1 itemprop="title">' NameInputFile '</h1>\r'  ...
        '<script type="text/javascript">\r'  ...
        '<!--\r'  ...
        '   function Redirect() {\r'  ...
        'var l = document.getElementById(''link'');\r'  ...
        'l.click();\r'  ...
        '   }\r'  ...
        '   setTimeout(''Redirect()'', 400);\r'  ...
        '//-->\r'  ...
        '</script>\r'  ...
        '</head>\r'  ...
        ' <a href="matlab:web([docrootFS ''' NameInputFileHTML '''])"; target="_top" id="link">Link to formatted HTML documentation in Mathworks style of ''' NameInputFileHTML '''</a> \r'  ...
        '<P>If redirecting does not work you can see the proper' ...
        ' HTML documentation of this page in Mathworks style at the web address of the Robust Statistics Academy of the University of Parma (RoSA)' ...
        '<P> <a href="http://rosa.unipr.it/FSDA/' NameInputFile '.html">http://rosa.unipr.it/FSDA/' NameInputFile '.html</a></P>'...
        '<hr />' ...
        '<p style="background-color:#A9CCE3 "><em>Syllabus page indexed by builddocsearchdb for function: ' NameInputFile '</em></p>\r'  ...
        '<P>' NameInputFile '</P>\r'  ...
        '<P>' DescriptionInputFile '</P>\r'])  ...
        DescriptionHTML ...
        MoreAboutHTML ...
        ReferencesHTML ...
        '</html>'];
catch
    error('FSDA:CreateFSDApointerFiles:noHTML',['It was impossible to create HTML file' NameInputFile])
end

end