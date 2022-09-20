%% File which automatically creates toolbox file FSDA.mltbx for fileexchange

%% Beginning of code

% specify the version number, please use the format 'major.minor.revision'
newVersion = '8.5.34';

% Add the sentence which describes the new feature of the release
commentRelease='Added new datasets';

% Specify folder where to create the project
FSDAProjFolder='D:\tmp';

% Specify name of the file which will contain the project
ProjectFileName='FSDAproject.prj';

% If for examples FSDAProj is D:\tmp  This file will create:
% D:\tmp\FSDAproject.prj
% D:\tmp\ToolboxPackagingConfiguration.prj
% D:\tmp\FSDA.mltbx
% D:\tmp\FSDA which will contains a selection of files from git repository
% https://github.com/UniprJRC/FSDA

%% Preliminary operations

% Navigate into FSDAProjFolder
try
    cd(FSDAProjFolder)
catch
    disp(['Supplied path "' FSDAProjFolder '" does not exist'])
    error('FSDA:CreatetoolboxFile:WrgPath','Wrong input path')
end

% Get filesep
fsep=filesep;

% Check if subfolder resources exists
% If the answer is yes delete it in order to start from scratch
if isfolder('resources') == true
    rmdir('resources','s')
end

% Check if subfolder FSDA exists
% If the answer is yes delete it in order to start from scratch
if isfolder('FSDA') == true
    rmdir('FSDA','s')
end


% Check if ProjectFileName exists
% If the answer is yes delete it in order to start from scratch
if isfile(ProjectFileName) == true
    delete(ProjectFileName)
end

% Check if FSDA.mltbx exists
% If the answer is yes delete it in order to start from scratch
if isfile('FSDA.mltbx') == true
    delete('FSDA.mltbx')
end

% Check if ToolboxPackagingConfiguration.prj exists
% If the answer is yes delete it in order to start from scratch
if isfile('ToolboxPackagingConfiguration.prj') == true
    delete('ToolboxPackagingConfiguration.prj')
end

% Create the project inside FSDAProjFolder
% File  Blank_project.prj   will be created
FSDAproj = matlab.project.createProject(FSDAProjFolder);
% Rename file  Blank_project.prj  into ProjectFileName
movefile('Blank_project.prj',ProjectFileName)
% Label the project (before was "blank_project")
FSDAproj.Name = "FSDA (Flexible Statistics Data Analysis)";

%% CLONE FROM GIT
% Clone from github repo
try
    !git clone https://github.com/UniprJRC/FSDA.git --progress
catch
    disp('Due to network connection it was not possibe to download from')
    disp('https://github.com/UniprJRC/FSDA.git')
    error('FSDA:CreatetoolboxFile:WrgPath','Connection problem')
end


%% REMOVE UNNECESSARY FOLDERS inside project folder
FSroot='FSDA';

% remove subfolder _automation_tools
folder_to_remove=[FSroot fsep '_automation_tools'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder _development
folder_to_remove=[FSroot fsep '_development'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder _TODO
folder_to_remove=[FSroot fsep '_TODO'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder helpfiles/XML
folder_to_remove=[FSroot fsep 'helpfiles' fsep 'XML'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder .git
folder_to_remove=[FSroot fsep '.git'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder .github
folder_to_remove=[FSroot fsep '.github'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder .circleci
folder_to_remove=[FSroot fsep '.circleci'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder Univ
folder_to_remove=[FSroot fsep 'Univ'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder docker
folder_to_remove=[FSroot fsep 'docker'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

%% REMOVE UNNECESSARY FILES

% All png files  inside helpfiles\FSDA\images
% delete([FSroot fsep 'helpfiles' fsep 'FSDA' fsep 'images' fsep '*.png'])

% files mlx inside inside subfolder examples
delete([FSroot fsep 'examples' fsep 'examples_categorical.mlx'])
delete([FSroot fsep 'examples' fsep 'examples_multivariate.mlx'])
delete([FSroot fsep 'examples' fsep 'examples_regression.mlx'])
delete([FSroot fsep 'examples' fsep 'examples_MixSim.mlx'])

% file pptx which explains the flow chart
delete([FSroot fsep 'utilities_help' fsep 'FlowChart.pptx'])

% remove license files
delete([FSroot fsep 'eupllicense.pdf'])
delete([FSroot fsep 'Copyright notice.pdf'])

% remove readme and installation notes files
delete([FSroot fsep 'installationNotes.docx'])
delete([FSroot fsep 'installationNotes.pdf'])

% remove md files
delete([FSroot fsep 'readme.md'])
delete([FSroot fsep '404.md'])
delete([FSroot fsep 'CODE_OF_CONDUCT.md'])
delete([FSroot fsep 'CONTRIBUTING.md'])

delete([FSroot fsep 'requirements.txt'])
% delete([FSroot fsep 'package.json'])


%% Add files to project
% Add all files to the project which are inside folder FSDA
% and subfolders
addFolderIncludingChildFiles(FSDAproj,FSroot);

% Check that for example file addFSDA2path.m in the main folder of FSDA has
% been added
% findFile(FSDAproj,'FSDA/addFSDA2path.m')

%% Create searchable database

% save current path
oldpath = path;
path2add=[FSDAProjFolder fsep 'FSDA'];
addpath(path2add);

FileName=[FSroot filesep 'addFSDA2path'];
FullPath=which(FileName);
%Navigate to the main folder of FSDA
FSDAroot=fileparts(FullPath);
% FSDApointers = full path for folder which creates pointer files
FSDApointers=[FSDAroot fsep 'helpfiles' fsep 'pointersHTML'];

% Starting in R2022a, the builddocsearchdb function creates the subfolder
% helpsearch-v4 to contain the search database files. Previously,
% builddocsearchdb created a subfolder named helpsearch-v3.
% To ensure the documentation for FSDA toolbox is searchable in any version
% we have to run buildocsearchdb also on R2021b

% The following line assumes that the path of MATLAB 2021b is 
% C:\Program Files\MATLAB\R2021b

% bdocsearch2021b=['eval(''addpath(''''' FSDAroot ''''');builddocsearchdb(''''' FSDApointers ''''')'')'];
% run builddocsearchdb and quit
% bdocsearch2021b=['eval(''addpath(''''' FSDAroot ''''');builddocsearchdb(''''' FSDApointers ''''');quit'')'];
bdocsearch2021b=['addpath(''' FSDAroot ''');builddocsearchdb(''' FSDApointers ''');quit'];
[status]=system(['"C:\Program Files\MATLAB\R2021b\bin\matlab.exe" -nodesktop -nosplash -noFigureWindows -r ' bdocsearch2021b ]);

if status~=0
       error('FSDA:CreatetoolboxFile:WrgBdocSdb','Could not execute builddocsearchdb in MATLAB 2021b')
end

builddocsearchdb(FSDApointers);

% restore previous path
path(oldpath);

%% Add FSDA paths to the project
pt=cell(15,1);
pt{1}=FSroot;
pt{2}=[FSroot fsep 'multivariate'];
pt{3}= [FSroot fsep 'regression'];
pt{4}= [FSroot fsep 'clustering'];
pt{5}= [FSroot fsep 'graphics'];
pt{6}= [FSroot fsep 'datasets' fsep 'regression'];
pt{7}= [FSroot fsep 'datasets' fsep 'multivariate'];
pt{8}= [FSroot fsep 'datasets' fsep 'multivariate_regression'];
pt{9}= [FSroot fsep 'datasets' fsep 'clustering'];
pt{10}= [FSroot fsep 'combinatorial'];
pt{11}= [FSroot fsep 'utilities'];
pt{12}= [FSroot fsep 'utilities_stat'];
pt{13}= [FSroot fsep 'utilities_help'];
pt{14}= [FSroot fsep 'examples'];
pt{15}= [FSroot fsep 'FSDAdemos'];

for i=1:length(pt)
    folderonpath = addPath(FSDAproj,pt{i});
end

%% Run dependency analyzer
updateDependencies(FSDAproj);


%% Copy file ToolboxPackagingConfiguration.prj into FSDAProjFolder (current folder)
copyfile([FSDAroot fsep 'utilities_help' fsep 'ToolboxPackagingConfiguration.prj'],FSDAProjFolder)

%% Set release compatibility in ToolboxPackagingConfiguration.prj file
setToolboxStartEnd('ToolboxPackagingConfiguration.prj')

%% Publish contents file in the root inside subfolder html
% This instruction is necessary in order to display subfolder examples in
% Mathworks web site
publish([FSDAroot filesep 'Contents.m']);

%% Package toolbox and create file FSDA.mltbx

toolboxFile = 'ToolboxPackagingConfiguration.prj';

% set the version number
previousVersion = matlab.addons.toolbox.toolboxVersion(toolboxFile,newVersion);

outputFile ='FSDA.mltbx';
matlab.addons.toolbox.packageToolbox(toolboxFile, outputFile)

%% Close the project
close(FSDAproj)

% Open project
% FSDAproj = openProject(FSDAProjFolder);

%% Copy FSDA.mltbx to Github, create a new release and tag it

% copy d:\tmp\FSDA.mltbx ./bin/FSDA.mltbx

FSrootGitHub = fileparts(which('docsearchFS.m'));

% FSrootGitHub = 'C:\FSDA';

% create 'bin' subfolder
mkdir( [FSrootGitHub fsep 'bin']);

% mkdir( [FSrootGitHub fsep '.github/workflows/']);

[status,msg] = copyfile([FSDAProjFolder  fsep 'FSDA.mltbx'],[FSrootGitHub fsep 'bin']);

% Navigate to FSDA root (not project root)
cd(FSrootGitHub)

% Remove the cached key for 192.168.1.123 on the local machine to remove
% warning
!ssh-keygen -R 140.82.121.4
!git add ./bin/FSDA.mltbx
!git commit -m "added last build of FSDA.mltbx"
!git push

% tag the release and start the upload of FSDA.mltbx to release assets
% via the Github workflow and insert comment
eval(['!git tag -a ', newVersion ' -m "' commentRelease '"'])

!git push --tags

%% delete temporary <FSDAroot>/bin folder on GitHub
% Note that we delete the file after waiting 2 minutes just to make sure
% that github workflow has come to an end.

disp('Pausing for 2 minutes before deleting the file FSDA.mltbx from subfolder /bin on GitHub')
pause(60*2);

% Remove mltbx file
!git rm -r ./bin

% commit only the removal of the 'FSDA.mltbx' file just removed by 
% the command 'git rm -r ./bin' and comment this modificantion with the
% 3 numbers of the new version.
eval(['!git commit -m " FSDA version ' newVersion ' now released."'])

% push these modifications
!git push



