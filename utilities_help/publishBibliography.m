function [fstring,citsCell]=publishBibliography(InputCell,OUT, varargin)
%publishBibliography enables to create web page which contains the references inside the input .m files
%
%<a href="matlab: docsearchFS('publishBibliography')">Link to the help function</a>
%
% This routine uses as input the cell which is created with routine
% makecontentsfileFS.m and the cell created by publishFSallFiles.m and uses
% template file bibliographyEmpty.html
% which is contained inside:
% (main root of FSDA) filesep 'helpfiles' filesep 'FSDA;
% to create in a fully automatic way the bibliography file containing all
% the references.
% The output file (if the path is not specified in input option outputDir)
% will be created inside:
% (main root of FSDA) filesep 'helpfiles' filesep 'FSDA;
% and will have name bibliography.html
%
%
% Required input arguments:
%
%   InputCell: Cell created by function makecontentsfileFS.m. Cell. Cell
%              containing the names of all files which have to be
%              processed to produce the HTML file contanining the bibliography.
%  OUT       : Cell created by function publishFSallFiles.m. Cell.
%              Cell of length size(InputCell,1). OUT{i} contains
%              the output of the application of routine publishFS.m
%              to i-th file, i=1, 2, ..., lengthInputCell.
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
%       fstring: HTML string with all the references. Character.
%                The HTML file bibliography.html is also produced inside
%                folder specified by input option "outputDir" .
%     citsCell : References, associated file and path. Cell.
%                Cell with three columns.
%                First column = citation.
%                Second column = FileName of associated file where citation
%                was found.
%                Third column = Path of associated file.
%
% More About:
%
% Guidelines to create the references.
% [1] Example of book citation (note that the publisher town is not given and
% that the title of the book must be written in ""):
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag.
% [2] Example of article citation (note that just the volume is given and that
% the name of the Journal is inside ""):
% Wong, C.K. and Easton, M.C. (1980), An Efficient Method for Weighted
% Sampling Without Replacement, "SIAM Journal of Computing", Vol. 9, pp.
% 111-113.
% [3] Example of InCollection (note that just the volume is given and that
% the name of the Journal is inside""):
% Riani, M., Cerioli, A., Atkinson, A.C., Perrotta, D. and Torti, F. (2008),
% Fitting Mixtures of Regression Lines with the Forward Search, in:
% "Mining Massive Data Sets for Security", Fogelman-Soulie F. et al. Eds.,
% pp. 271-286, IOS Press.
% [4] IMPORTANT.
% Extra information which is attached to a citation inside the MATLAB .m
% file but does not have to appear in the output file bibliography.html
% must be inserted inside square parenthesis.
% For example in the citation below inside the .m file:
%
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [quoted below as RPRH].
%
% The string [quoted below as RPRH] does not appear in the bibliography
%
%
% See also:    publishfunctionAlpha.m, publishfunctionCate.m, publishFS.m
%
%
% References:
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('publishBibliography')">Link to the help function</a>
%
%$LastChangedDate:: 2018-06-12 19:18:50 #$: Date of the last commit
%
% Examples:
%
%

%{
    % publishBibliography with all the default options.
    % Create the requested input arguments.
    FileName='addFSDA2path';
    FullPath=which(FileName);
    %Navigate to the main folder of FSDA
    FSDAroot=fileparts(FullPath);
    cd(FSDAroot)
    % Specify subfolders of main folders for which contents file has to be
    % created
    InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
        'examples' 'utilities' 'utilities_stat' 'utilities_help'};
    ExclDir={'privateFS'  'datasets'};
    % Create list of folders which must have the personalized contents file
    list = findDir(FSDAroot,'InclDir',InclDir,'ExclDir',ExclDir);
    % Crete personalized contents file for main folder of FSDA
    % and required subfolders.
    force=false;
    [FilesIncluded,FilesExcluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',force,'printOutputCell','Contents.m');
    [~,OUT]=publishFSallFiles(FilesIncluded, 'evalCode','false');
    % Create HTML file containing all the items which make up the bibliography
    fileBiblio=publishBibliography(FilesIncluded,OUT);
%}

%{
    % Example of use of option outputDir.
    % We assume that path outpuPath exist and that inside this path there
    % is template file
    outputPath='C:\temp';
    fileBiblio=publishBibliography(FilesIncluded,OUT,'outputDir',outputPath);
%}

%{
    % Example where second output is returned.
    % In this case output Citations contains the list of the citations
    % which have been found.
    [fileBiblio,Citations]=publishBibliography(FilesIncluded,OUT);
    disp(Citations)
%}


%% Beginning of code
% % Use file separator of current operating system
% % \ = Windows
% % / = Unix
fsep=filesep;


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

if nargin>3
    options=struct('outputDir',outputDir,'webhelp',webhelp);
    
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:publishBilbiography:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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


citsCell=extractRefs(OUT, InputCell);
% Sort all the rows by alphabetical order based on citation
citsCell=sortrows(citsCell,1);
% citsCell=citsCell;
% remove in any citation what is inside the square brackets
% The expression ' .\[.*\]' means find a string with any character (but just one) followed
% by [ containing any sequence of words up to ]
for i=1:length(citsCell)
    citsCell{i,1}=regexprep(citsCell{i,1},'.\[.*\]','');
end


% cits=unique(citsCell(:,1),'sorted');
[~,ia]=unique(citsCell(:,1));

cits=citsCell(ia,:);

% Remove duplicate references
uniquesTab=cell2table(cits);

% Remove reference which are virtually equal, that is the lines which have
% an edit distance https://en.wikipedia.org/wiki/Edit_distance smaller than
% 8
%
boo=false(size(cits,1),1);

for i=2:size(uniquesTab,1)
    if wfEdits(char(uniquesTab{i,1}),char(uniquesTab{i-1,1}))<8
        boo(i)=true;
    end
end
uniquesTab(boo,:)=[];

% Create 'bibliography.html' file for all uniques references inside all
% MATLAB source files.

% strInsert is the string which will contain all the citations in HTML format
strInsert='';

for ii = 2:size(uniquesTab,1)
    
    refEntry=char(uniquesTab{ii,1});
    refText=refEntry;
    
    % write in italic what is inside double quotes
    dbquotesVect=regexp(refText,'\"');
    if (~isempty(dbquotesVect))
        tmpText=[ refText(1:dbquotesVect(1)-1) '<i>' ...
            refText(dbquotesVect(1)+1:dbquotesVect(2)-1) '</i>' ...
            refText(dbquotesVect(2)+1:end)];
        refText=tmpText;
    end
    
    % delete every citation that is inside square brackets
    % also including brackets (inside html help files only!)
    quadStart=regexp(refText,'\[','once');
    if (~isempty(quadStart))
        quadEnd=regexp(refText,'\]','once');
        refText=[ refText(1:quadStart-1) refText(quadEnd+1:end)];
    end
    
    strInsert = [strInsert '<p>' refText '</p>' newline]; %#ok<AGROW>
end

 strInsert = [strInsert '<p><hr /> <center><i>This page has been automatically generated' ...
     ' by our routine <a href="publishBibliography.html">publishBibliography</a></i> </center></p>' newline];


%% Open input bibliographyEmpty.html file, put it in a string and do a series of preliminary operations
FileWithFullPath=[outputDir fsep 'bibliographyEmpty.html'];
fileID = fopen(char(FileWithFullPath), 'r');

if fileID==-1
    disp(['Output path: '''  char(outputDir) ''' must  contain a file named  ''bibliographyEmpty.html'''])
    error('FSDA:publishBibliography:WrongPath','Output path does not have the input file')
end

% Insert the file into fstring
fstring=fscanf(fileID,'%c');

[~,finHTMLTEXT] =regexp(fstring,'BEGINNING OF FSDA TEXT');

% Insert inside fstring all the citations
fstring=[fstring(1:finHTMLTEXT+3) strInsert fstring(finHTMLTEXT+4:end)];

name='bibliography';

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
    
    error('FSDA:publishBibliography:WrngOutFolder',errmsg);
end

fprintf(file1ID,'%s',fstring);
fclose('all');

end

function References=extractRefs(OUT, FilesIncluded)
% Creates a table of references from input struct OUT and list of files
% FilesIncluded

% Make sure that empty elements of OUT are removed
boo=cellfun(@isempty,OUT);
OUT=OUT(~boo);
FilesIncluded=FilesIncluded(~boo,:);

ri=size(OUT,1);
% References is a cell with three columns
% First column = citation
% Second column = FileName of associated file
% Third column = Path of associated file
References=cell(5000,3);
j=1;
for ii=1:ri
    
    fileRefs=OUT{ii,1}.References;
    refsRows=size(fileRefs,1);
    for jj=1:refsRows
        References{j,1}=removeExtraSpacesLF(fileRefs{jj});
        References(j,2)=FilesIncluded(ii,1);
        References(j,3)=FilesIncluded(ii,9);
        j=j+1;
    end
end
References=References(1:j-1,:);
end

function d = wfEdits(S1,S2)
% Wagner–Fischer algorithm to calculate the edit distance / Levenshtein distance.
%
N1 = 1+numel(S1);
N2 = 1+numel(S2);
%
D = zeros(N1,N2);
D(:,1) = 0:N1-1;
D(1,:) = 0:N2-1;
%
for r = 2:N1
    for c = 2:N2
        D(r,c) = min([D(r-1,c)+1, D(r,c-1)+1, D(r-1,c-1)+~strcmpi(S1(r-1),S2(c-1))]);
    end
end
d = D(end);
%
end

%FScategory:UTIHELP
