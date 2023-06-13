function plan = buildfile

plan = buildplan(localfunctions);

% Build doc before packaging toolbox as it is needed in the toolbox
plan("toolbox").Dependencies = "doc";

% Make the "toolbox" task the default task in the plan
plan.DefaultTasks = "toolbox";
end

function docTask(~)
% This task builds the doc search DB for the current version of MATLAB - the
% expected output will be in the folder ./helpfiles/pointersHTML
cleanup = iCdWithRevert(fullfile(iGetRootFolder, "utilities_help", "build")); %#ok<NASGU>
buildDocSearchForToolbox
end

function toolboxTask(~)
% This task packages the toolbox MLTBX file - the expected output will be
% in the ./bin/ folder (defined in the createMLTBX file)
cleanup = iCdWithRevert(fullfile(iGetRootFolder, "utilities_help", "build")); %#ok<NASGU>
createMLTBX
end

function root = iGetRootFolder()
% buildfile.m is always in the root of the toolbox - so use this definition
% to get the root of the toolbox from the location of buildfile.m
persistent THE_ROOT
if isempty(THE_ROOT)
    thisFile = mfilename("fullpath");
    THE_ROOT = fileparts(thisFile);
end
root = THE_ROOT;
end

function cleanup = iCdWithRevert(folder)
% Change to folder (probably build folder) and revert afterwards to pwd
cdir = pwd;
cleanup = onCleanup(@() cd(cdir));
cd(folder)
end