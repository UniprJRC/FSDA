function createMLTBX 
%% Create toolbox file for FSDA, both locally and in a CI action

% Specify the version number, please use the format 'major.minor.revision'
% During a github build this will be taken from github.ref_name and stored
% in GITHUB_ENV - and for local sandbox builds that don't specify this we
% simply default to a large 10.0.0 release
newVersion = getenv('GITHUB_RELEASE_TAG');

if isempty(newVersion)
    disp('Missing a GITHUB_RELEASE_TAG to build the toolbox : using version 10.0.0')
    newVersion = "10.0.0";
end

% Builds need to start from the root to cleanly define all the relative
% paths, but revert back to where ever we started from when the "cleanup"
% object goes out of scope
[FSDAroot, cleanup] = changeDirToRootWithCleanup; %#ok<ASGLU>

% Get filesep
fsep=filesep;

%% Publish contents file in the root inside subfolder html
% This instruction is necessary in order to display subfolder examples in
% Mathworks web site
publish('Contents.m');



%% Create toolbox project file

% Create our desired toolbox options from the default toolbox package
% configuration project fil in the root. This specifies as little as
% possible ... but it does define the following:
%
% 1. Our File Exchange UUID (currently 20669fbc-61ca-4050-bc87-575422f4c0b8)
% 2. By default ALL files should be packaged in the toolbox - we will
%   exclude those that are not needed below using some helpers.
options = matlab.addons.toolbox.ToolboxOptions(fullfile(FSDAroot, "defaultToolboxPackageConf.prj"));

% Firstly there are a set of folders in this repository that we do not want
% in the packaged toolbox - remove those using a helper function
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

% Secondly there are a set of files in this repository that we do not want
% in the packaged toolbox - remove those using a helper function
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


% Define the paths that we want to add to an installed MATLAB path
% NOTE - we need the root as well as some sub-folders so include the empty
% string at the beginning
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

% AND NOTE scalar expansion here with vector of pathsToAdd
options.ToolboxMatlabPath = FSDAroot + fsep + pathsToAdd;

% Define our TOOLBOX name, version and other metadata
options.ToolboxName = "FSDA";

options.ToolboxVersion = newVersion;

options.Description="Flexible Statistics and Data Analysis (FSDA) extends MATLAB for " + ...
    "a robust analysis of data sets affected by different sources of heterogeneity. " + ...
    "It is open source software licensed under the European Union Public Licence (EUPL). " + ...
    "FSDA is a joint project by the University of Parma and the Joint Research Centre " + ...
    "of the European Commission."; 

options.Summary="Flexible Statistics Data Analysis Toolbox";

options.AuthorName= "Marco Riani";

options.AuthorEmail = "FSDA@unipr.it";

options.AuthorCompany = "University of Parma (UNIPR) and Joint Research Centre of the " + ...
    "European Commission(JRC).";
 
% Define the final architectures that we will correctly work on
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

mkdir(FSDAroot, 'bin')
options.OutputFile = fullfile(FSDAroot, "bin", "FSDA");

% Display the options during build in case the github action fails and we
% need some debugging output.
disp(options)

% Finally - package the toolbox
disp('Package toolbox and create file FSDA.mltbx')
matlab.addons.toolbox.packageToolbox(options);

end

function options = removeFoldersFromToolboxPackage(options, foldersToRemove)
% This helper function removes all files in a a set of folders (defined
% relative to the root of the package) from the ToolboxFiles list. This is
% done using string comparison startsWith to find all the currently listed
% files in ToolboxFiles that start with the desired folders.
root = options.ToolboxFolder;
% NOTE - likely scalar expansion of root to accomodate multiple
% foldersToRemove - DO NOT replace with fullfile which does not support
% expansion.
% Also note addition of filesep on the end to ensure that ONLY files in
% folders are removed. Without this a file that began with the name
% folderToRemove and continued with more chars would ALSO be removed.
foldersToRemove = root + filesep + foldersToRemove + filesep;
% NOTE - startsWith takes multiple patterns for second option - we are in
% that case
itemsToRemove = startsWith(options.ToolboxFiles, foldersToRemove);
options.ToolboxFiles(itemsToRemove) = [];
end

function options = removeFilesFromToolboxPackage(options, filesToRemove)
% This helper function removes all files listed from the ToolboxFiles list.
% This is done using string comparison matches to find all the currently
% listed files in ToolboxFiles that match the desired list.
root = options.ToolboxFolder;
% NOTE likely scalar expansion of root to accomodate multiple
% filesToRemove - DO NOT replace with fullfile which does not support
% expansion
filesToRemove = root + filesep + filesToRemove;
% NOTE - matches takes multiple patterns for second option - we are in
% that case
itemsToRemove = matches(options.ToolboxFiles, filesToRemove);
options.ToolboxFiles(itemsToRemove) = [];
end