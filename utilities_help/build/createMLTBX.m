function createMLTBX 
%% File which automatically creates toolbox file FSDA.mltbx for fileexchange


%% Beginning of code

% specify the version number, please use the format 'major.minor.revision'
newVersion = getenv('GITHUB_RELEASE_TAG');

if isempty(newVersion)
    disp('Missing a GITHUB_RELEASE_TAG to build the toolbox : using version 10.0.0')
    newVersion = "10.0.0";
end

% Specify folder where to create the toolbox project
[FSDAroot, cleanup] = changeDirToRootWithCleanup; %#ok<ASGLU>

% Get filesep
fsep=filesep;

%% Publish contents file in the root inside subfolder html
% This instruction is necessary in order to display subfolder examples in
% Mathworks web site
publish('Contents.m');



%% Create toolbox project file

% uuid identified of FSDA
uuid = '20669fbc-61ca-4050-bc87-575422f4c0b8';

% Note that when matlab.addons.toolbox.ToolboxOptions the files attached to
% the project are automatically added inside options.ToolboxFiles 
options = matlab.addons.toolbox.ToolboxOptions(FSDAroot, uuid);

options = removeFoldersFromToolboxPackage(options, [ ...
    "_automation_tools"
    "_development"
    "_TODO"
    ".buildtool"
    ".circleci"
    ".git"
    ".github"
    "bin"
    "docker"
    "helpfiles" + fsep + "XML"
    "utilities_help" + fsep + "build"
    "Univ"
    ]);

options = removeFilesFromToolboxPackage(options, [...
    "examples" + fsep + "examples_categorical.mlx"
    "examples" + fsep + "examples_multivariate.mlx"
    "examples" + fsep + "examples_regression.mlx"
    "examples" + fsep + "examples_MixSim.mlx"
    "utilities_help" + fsep + "FlowChart.pptx"
    "eupllicense.pdf"
    "Copyright notice.pdf"
    "installationNotes.docx"
    "installationNotes.pdf"
    "buildfile.m"
    "FSDA.prj"
    "readme.md"
    "404.md"
    "CODE_OF_CONDUCT.md"
    "CONTRIBUTING.md"
    "helpfiles" + fsep + "FSDA" + fsep + "images" + fsep + "githubimgexamples.jpg"
    "helpfiles" + fsep + "FSDA" + fsep + "images" + fsep + "githubimgindex.jpg"
    "helpfiles" + fsep + "FSDA" + fsep + "images" + fsep + "githubimgtutorials.jpg"    
    ]);


% add FSDA paths
pathsToAdd = [ ...
    ""
    "multivariate"
    "regression"
    "clustering"
    "graphics"
    "datasets" + fsep + "regression"
    "datasets" + fsep + "multivariate"
    "datasets" + fsep + "multivariate_regression"
    "datasets" + fsep + "clustering"
    "combinatorial"
    "utilities"
    "utilities_stat"
    "utilities_help"
    "examples"
    "FSDAdemos"    
    ];

% NOTE scalar expansion here with vector of pathsToAdd
options.ToolboxMatlabPath = FSDAroot + fsep + pathsToAdd;

% add toolbox name
options.ToolboxName = "FSDA";

% add toolbox version
options.ToolboxVersion = newVersion;

% Detailed description of the toolbox.
options.Description="Flexible Statistics and Data Analysis (FSDA) extends MATLAB for " + ...
    "a robust analysis of data sets affected by different sources of heterogeneity. " + ...
    "It is open source software licensed under the European Union Public Licence (EUPL). " + ...
    "FSDA is a joint project by the University of Parma and the Joint Research Centre " + ...
    "of the European Commission."; 

% add summary description of the toolbox
options.Summary="Flexible Statistics Data Analysis Toolbox";

% add toolbox author.
options.AuthorName= "Marco Riani";

% add email address address
options.AuthorEmail = "FSDA@unipr.it";

% added company
options.AuthorCompany = "University of Parma (UNIPR) and Joint Research Centre of the " + ...
    "European Commission(JRC).";
 
% add architecture support (option name changed in the release)
options.SupportedPlatforms.Win64 = true;
options.SupportedPlatforms.Maci64 = true;
options.SupportedPlatforms.Glnxa64 = true;
options.SupportedPlatforms.MatlabOnline = true;

% version compatibility
options.MinimumMatlabRelease = 'R2018a';
options.MaximumMatlabRelease = '';

% add big logo
options.ToolboxImageFile = fullfile(FSDAroot, "logoblue.jpg");

% add getting startup file
options.ToolboxGettingStartedGuide = fullfile(FSDAroot, 'doc', 'GettingStarted.mlx');

% add gallery files
options.AppGalleryFiles=[];
% options.OutputFile='testFSDA';
mkdir(FSDAroot, 'bin')

options.OutputFile = fullfile(FSDAroot, "bin", "FSDA");

disp(options)

disp('Package toolbox and create file FSDA.mltbx')
matlab.addons.toolbox.packageToolbox(options);


end

function options = removeFoldersFromToolboxPackage(options, foldersToRemove)
root = options.ToolboxFolder;
foldersToRemove = root + filesep + foldersToRemove + filesep;
itemsToRemove = startsWith(options.ToolboxFiles, foldersToRemove);
options.ToolboxFiles(itemsToRemove) = [];
end

function options = removeFilesFromToolboxPackage(options, filesToRemove)
root = options.ToolboxFolder;
filesToRemove = root + filesep + filesToRemove;
itemsToRemove = matches(options.ToolboxFiles, filesToRemove);
options.ToolboxFiles(itemsToRemove) = [];
end