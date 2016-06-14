function outHELP=CreateFSDApointerFiles(InputCell,OUT)
%CreateFSDApointerFiles create points to HTML help files
%
% Given that:
% 1) buildocsearchdb enables to create the database just for
% the files which are inside hte helpfile location specified in the
% info.xml file
% 2) all the entries which are found by buildocsearchdb are opened in the
% frame of the right (iframe)
% 3) FSDA team spent a lot of time to have proper HTML documentation which
% is not relagated to the right panel,
% 4) true HTML documentation file of FSDA are inside (docroot)/FSDA otherwise
% the lucene MATLAB search engine does not work.
%
% The purpose of this function is to take as input all the files which are
% inside (FSDA root folder)/helpfiles/FSDA
% and for each of them to create a corresponding pointer HTML file inside path
% (FSDA root folder)/helpfiles/pointersHTML
% which contains
% a) the minimum necessary information to be indexed by lucene search enegine
% b) a response.redirect which links to the true HTML file which is
% contained inside (docroot)/FSDA.
%
%
%{
    % Create personalized contents file which will be the input the
    % procedure CreateFSDApointerFiles.
    % Find full path of the main root where FSDA is installed
    FileName='addFSDA2path';
    FullPath=which(FileName);
    root=FullPath(1:end-length(FileName)-3);
    InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat'};
    ExclDir={'privateFS'  'datasets'};
    % Create list of folders which must be presonlized contents file
    list = findDir(root,'InclDir',InclDir,'ExclDir',ExclDir)
    % Crete personalized contents file for main folder of FSDA
    % and required subfolders.
    [InputCell,Excluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',false)
    % Create HTML pointer files
    CreateFSDApointerFiles(InputCell)
%}

% Copyright 2008-2016.
% Written by FSDA team
%

%% Beginning of code
[FSDAroot]=fileparts(which('docsearchFS.m'));
fsep=filesep;
dirpathj=[FSDAroot fsep 'helpfiles' fsep 'pointersHTML'];

for i=1:size(InputCell,1)
    NameInputFile=InputCell{i,6};
    PurposeInputFile=InputCell{i,7};
    NameOutputFile=[NameInputFile '.html'];
    DescrRef=OUT{i};
    [fileID,errMsg] = fopen([dirpathj filesep NameOutputFile],'w');
    if fileID < 0
        outHELP=false;
        error(message('MATLAB:filebrowser:MakeContentsFileOpenError', errMsg))
    end
    % Create string to include the HTML pointer file
    fstring=stringHTMLpointer(NameInputFile,PurposeInputFile,DescrRef);
    fprintf(fileID,'%s',fstring);
    fclose(fileID);
end
outHELP=true;


% open each file extract relevant information



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
    fstring=[sprintf(['<!DOCTYPE HTML> \r'  ...
        '<html itemscope="" xmlns="http://www.w3.org/1999/xhtml"> \r'  ...
        '<head> \r'  ...
        '<title>' NameInputFile '</title> \r'  ...
        '<meta content="refpage" name="chunktype">\r'  ...
        '<meta content="function:' NameInputFile  ' " itemprop="refentity" name="refentity">\r'  ...
        '<meta content="fcn" itemprop="pagetype" name="toctype">\r'  ...
        '<meta content="ref/function" itemprop="infotype" name="infotype" />\r'  ...
        '<meta content="' NameInputFile ' ' DescriptionInputFile '" itemprop="description" name="description" />\r'  ...
        '<script type="text/javascript">\r'  ...
        '<!--\r'  ...
        '   function Redirect() {\r'  ...
        '      window.location= "matlab:web([docrootFS ''' NameInputFileHTML '''])";\r'  ...
        '   }\r'  ...
        '   setTimeout(''Redirect()'', 10);\r'  ...
        '//-->\r'  ...
        '</script>\r'  ...
        '</head>\r'  ...
        '<p style="background-color:#A9CCE3 "><em>Syllabus page indexed by builddocsearchdb for function: ' NameInputFile '</em></p>\r'  ...
        '<h1 itemprop="title">' NameInputFile '</h1>\r'  ...
        '<h1 class="reftitle">' NameInputFile ' ' DescriptionInputFile '</h1>\r'  ...
        '<h2>Description</h2>\r'])  ...
        '<P>' Description '</P>'  ...
        sprintf('<h2>More About</h2>\r')  ...
        '<P>' MoreAbout '</P>'  ...
        sprintf('<h2>References</h2>\r')  ...
        '<P>' References '</P>'  ...
        '</html>'];
catch
    ddd=1;
end

end