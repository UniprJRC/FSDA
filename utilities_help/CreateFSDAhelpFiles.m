function outHELP=CreateFSDAhelpFiles()
%CreateFSDAhelpFiles prepares all documentation files 

% The purpose of this function is to automatize all steps which lead to the
% HTML documentation
% The procedure is follows:
% 1) Create (using function makecontentsfileFS) a personalized contents
% file for each requested subfolder and the main folder of FSDA. Contents
% file in the main folder will contain the list of all requested files in
% the requested subfolders
%
% 2) Create for each .m file, filtered by previous procedure
% makecontentsfileFS.m  the corresponding HTML file (to accomplish
% this task use procedure publishFSallFiles).
% Check that all HTML files have been generated correctly.
%
% 3) Create files function-alpha.html and function-cate.html (list of
% functions in alphabetical and categorical order)
%
% 4) Create file function.alpha.txt which contains the list of all files
% separated by commas. This file is necessary in order to create the right
% and left buttons which enable us to navigate in alphabetical order inside the
% HTML navigation system
% 
% 5) Creare all pointer files using routine CreateFSDApointerFiles
%
% Remark: remember that your setup program must execaute command
% builddocsearchdb in folder [FSDAroot filesep 'helpfiles'  filesep 'pointersHTML'] 
% to create a searchable help
%
%
% Required input arguments:
%
%  Output:
%
%         outHELP:   structure which contains the following fields
%
%            outHELP.FilesIncluded  = list of files for which HTML help has
%                                       to be created 
%           outHELP.FilesIncluded   = list of files for which HTML help
%                                       does not have to be created
%            outHELP.FilesWithProblems = list of files with problems in the
%                                        generation of HTML help files
%              outHELP.fileAlpha    = string containing fucntion-alpha.html
%              outHELP.fileCate    = string containing fucntion-cate.html
%
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%$LastChangedDate::                      $: Date of the last commit


%% STEP 1 create personalized contents files
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
force=true;
[FilesIncluded,FilesExcluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',force,'printOutputCell','Contents.m');
disp('List of files which have been excluded (with path)')
disp(FilesExcluded(:,[1 9]))

%% STEP 2: create HTML for all files filtered using makecontentsFilesFS
[FilesWithProblems,OUT]=publishFSallFiles(FilesIncluded);

% Check correctness of HTML link inside each .m file
chkHTMLlink=cell2mat(FilesWithProblems(:,6));
if max(chkHTMLlink)>0
    disp('Files with wrong reference to HTML page')
    disp(FilesWithProblems{chkHTMLlink==1,1})
    error('FSDA:CreateFSDAhelpFiles','Files with wrong references in docsearchFS')
end

% Check correctness of HTML file creation
boo=~cellfun('isempty',FilesWithProblems(:,5));
seq=1:size(FilesWithProblems,1);
IndexesofFiles=seq(boo);
if ~isempty(IndexesofFiles)
    disp('Files whose HTML reference page could not be created')
    for i=1:length(IndexesofFiles)
        disp(FilesWithProblems{IndexesofFiles(i),1})
    end
    error('FSDA:CreateFSDAhelpFiles','Files without HTML web page')
end

%% STEP 3: create alphabetical list of functions and txt file
fsep=filesep;

% Make sure one more time you are inside main root of FSDA
cd(fileparts(which('docsearchFS.m')))
% Create HTML file containing alphabetical list of functions
fileAlpha=publishFunctionAlpha(FilesIncluded,'CreateTxtFile',true);
% open html file in web browser
outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA\function-alpha.html'];
web(outputOFHtmlHelpFile,'-browser');

fsep=filesep;
outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA\function-alpha.txt'];
% open outfile txt in web browser
disp('Check .txt file')
web(outputOFHtmlHelpFile,'-browser');


%% STEP 4: create categorical list of functions
fsep=filesep;

% Make sure one more time you are inside main root of FSDA
cd(fileparts(which('docsearchFS.m')))
% Create HTML file containing categorical list of functions
fileCate=publishFunctionCate(FilesIncluded);
% open outfile file in web browser
outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA\function-cate.html'];
web(outputOFHtmlHelpFile,'-browser');

%% STEP 5: create HTML pointer files

 h=CreateFSDApointerFiles(FilesIncluded,OUT);
if h
    disp('Successful creation of pointer files')
end

%% STEP 6: create searchable database with different versions of MATLAB
FileName='addFSDA2path';
FullPath=which(FileName);
%Navigate to the main folder of FSDA
FSDAroot=fileparts(FullPath);
% Navigate to subfolder which contains pointerHTML
pointersHTMLroot=[FSDAroot filesep 'helpfiles'  filesep 'pointersHTML'];
% Create searchable database
builddocsearchdb(pointersHTMLroot)

%% Now if all was well let us do the setup.exe
disp('Congratulations the FSDA package is ready to be deployed')

% Store all quantities inside structure outHELP
outHELP=struct;
outHELP.FilesIncluded=FilesIncluded;
outHELP.FilesExcluded=FilesExcluded;
outHELP.FilesWithProblems=FilesWithProblems;
outHELP.fileAlpha=fileAlpha;
outHELP.fileCate=fileCate;

end
