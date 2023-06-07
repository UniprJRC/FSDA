function plan = buildfile

plan = buildplan(localfunctions);

plan("toolbox").Dependencies = "doc";

% Make the "createMLTBX" task the default task in the plan
plan.DefaultTasks = "toolbox";
end

function docTask(~)
cleanup = iCdWithRevert(fullfile(iGetRootFolder, "utilities_help", "build")); %#ok<NASGU>
buildDocSearchForToolbox
end

function toolboxTask(~)
cleanup = iCdWithRevert(fullfile(iGetRootFolder, "utilities_help", "build")); %#ok<NASGU>
createMLTBX
end

function root = iGetRootFolder()
persistent THE_ROOT
if isempty(THE_ROOT)
    thisFile = mfilename("fullpath");
    THE_ROOT = fileparts(thisFile);
end
root = THE_ROOT;
end

function cleanup = iCdWithRevert(folder)
cdir = pwd;
cleanup = onCleanup(@() cd(cdir));
cd(folder)
end