%% File which automatically creates toolbox file FSDA.mltbx for fileexchange

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
if exist('resources','dir') ==7
    rmdir('resources','s')
end

% Chceck if ProjectFileName exists
% If the answer is yes delete it in order to start from scratch
if exist(ProjectFileName,'file') >0
    delete(ProjectFileName)
end

% Create the project inside FSDAProjFolder
% File  Blank_project.prj   will be created
FSDAproj = matlab.project.createProject(FSDAProjFolder);
% Rename file  Blank_project.prj  into ProjectFileName
movefile('Blank_project.prj',ProjectFileName)
% Label the project (before was "blank_project")
FSDAproj.Name = "FSDA (Flexible Statistics Data Analysis)";

%% CLONE FROM GIT
% Check if a subfolder FSDA of current folder exists and if not download it from github
if exist(fullfile([pwd filesep 'FSDA' filesep 'addFSDA2path.m']),'file')~=2
    try
        !git clone https://github.com/UniprJRC/FSDA.git
    catch
        disp('Due to network connection it was not possibe to download from')
        disp('https://github.com/UniprJRC/FSDA.git')
        error('FSDA:CreatetoolboxFile:WrgPath','Connection problem')
    end
end


%% REMOVE UNNECESSARY FOLDERS
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

% remove subfolder .github
folder_to_remove=[FSroot fsep '.circleci'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

% remove subfolder .Univ
folder_to_remove=[FSroot fsep 'Univ'];
if exist(folder_to_remove,'dir') ==7
    rmdir(folder_to_remove,'s')
end

%% REMOVE UNNECESSARY FILES

% All png files  inside helpfiles
delete([FSroot fsep 'helpfiles' fsep 'FSDA' fsep 'images' fsep '*.png'])

% files mlx inside inside subfolder examples
delete([FSroot fsep 'examples' fsep 'examples_categorical.mlx'])
delete([FSroot fsep 'examples' fsep 'examples_multivariate.mlx'])
delete([FSroot fsep 'examples' fsep 'examples_regression.mlx'])
delete([FSroot fsep 'examples' fsep 'examples_MixSim.mlx'])

% file pptx which explains
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
delete([FSroot fsep 'package.json'])


%% Add files to project
% Add all files to the project which are inside folder FSDA
% and subfolders
addFolderIncludingChildFiles(FSDAproj,FSroot);

% Check that for example file addFSDA2path.m in the main folder of FSDA has
% been added
% findFile(FSDAproj,'FSDA/addFSDA2path.m')


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

%% Create searchable database
FileName=[FSroot filesep 'addFSDA2path'];
FullPath=which(FileName);
%Navigate to the main folder of FSDA
FSDAroot=fileparts(FullPath);
builddocsearchdb([FSDAroot fsep 'helpfiles' fsep 'pointersHTML'])

%% Copy file ToolboxPackagingConfiguration.prj into FSDAProjFolder (current folder)
copyfile([FSDAroot fsep 'utilities_help' fsep 'ToolboxPackagingConfiguration.prj'],FSDAProjFolder)

%% Publish contents file in the root inside subfolder html
% This instruction is necessary in order to display subfolder examples in
% Mathworks web site
publish([FSDAroot filesep 'Contents.m'])

%% Package toolbox and create file FSDA.mltbx
outputFile ='FSDA.mltbx';
matlab.addons.toolbox.packageToolbox('ToolboxPackagingConfiguration.prj', outputFile)

%% Close the project
close(FSDAproj)

% Open project
% FSDAproj = openProject(FSDAProjFolder);


