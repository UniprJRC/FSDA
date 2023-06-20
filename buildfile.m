function plan = buildfile

plan = buildplan(localfunctions);

% Build doc before packaging toolbox as it is needed in the toolbox
plan("toolbox").Dependencies = "doc";

% Make the "toolbox" task the default task in the plan
plan.DefaultTasks = "toolbox";
end

function lintTask(~)
issues = codeIssues;
errors = issues.Issues(issues.Issues.Severity == "error", ...
    ["Location" "Severity" "Description"]);
assert(isempty(errors),formattedDisplayText(errors))
end

function testTask(~, cat2test, options)
arguments
    ~
    cat2test char = getenv('CATEGORY_TO_TEST')
    options.Performance (1,1) logical = false
end
runAllMyTestsFS(cat2test, Performance=options.Performance)
end

function docTask(context)
% This task builds the doc search DB for the current version of MATLAB - the
% expected output will be in the folder ./helpfiles/pointersHTML
cleanup = iCdWithRevert(fullfile(context.Plan.RootFolder, "utilities_help", "build")); %#ok<NASGU>
buildDocSearchForToolbox
end

function toolboxTask(context)
% This task packages the toolbox MLTBX file - the expected output will be
% in the ./bin/ folder (defined in the createMLTBX file)
cleanup = iCdWithRevert(fullfile(context.Plan.RootFolder, "utilities_help", "build")); %#ok<NASGU>
createMLTBX
end

function cleanup = iCdWithRevert(folder)
% Change to folder (probably build folder) and revert afterwards to pwd
cdir = pwd;
cleanup = onCleanup(@() cd(cdir));
cd(folder)
end