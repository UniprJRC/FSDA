function plan = buildfile

plan = buildplan(localfunctions);

% Build doc before packaging toolbox as it is needed in the toolbox
% removed for testing !!!!!
plan("toolbox").Dependencies = "doc";

% Make the "toolbox" task the default task in the plan
plan.DefaultTasks = "toolbox";
end

function checkTask(~)
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
cd toolbox
runAllMyTestsFS(cat2test, Performance=options.Performance)
end

function docTask(context)
% This task builds the doc search DB for the current version of MATLAB - the
% expected output will be in the folder ./helpfiles/pointersHTML
cleanup = iCdWithRevert(fullfile(context.Plan.RootFolder, "toolbox", "utilities_help", "build")); %#ok<NASGU>
buildDocSearchForToolbox
end

function toolboxTask(context)
% This task packages the toolbox MLTBX file - the expected output will be
% in the ./bin/ folder (defined in the createMLTBX file)
cleanup = iCdWithRevert(fullfile(context.Plan.RootFolder, "toolbox", "utilities_help", "build")); %#ok<NASGU>
createMLTBX
end

function releaseToGithubTask(~, opts)
% This task tags a new release and builds the toolbox MLTBX file -on GitHub-
% This task uses createMLTBX.m and GiHub Actions 
% e.g. buildtool releaseToGithub(Version="1.1.22",Comment="do not use just a test")
arguments
    ~
    opts.Version(1,1) string = ""
    opts.Comment(1,1) string = ""
end

if opts.Version == "" || opts.Comment == ""
    error('FSDA:ReleaseToGithub:IncorrectInputs',['You must specify ' ...
        '2 strings, Version and Comment, as input to the releaseToGithub built task'])
end

[OK, msg] = system(['git tag -a ' char(opts.Version)  ' -m "'  char(opts.Comment) '"']);
disp(msg);
if OK ~= 0
    error("FSDA:ReleaseToGithub:GitTagFailed", "git tag failed");
end
[OK, msg] = system("git push origin " + opts.Version);
disp(msg)
if OK ~= 0
    error("FSDA:ReleaseToGithub:GitPushFailed", "git push failed");
end
end

function cleanup = iCdWithRevert(folder)
% Change to folder (probably build folder) and revert afterwards to pwd
cdir = pwd;
cleanup = onCleanup(@() cd(cdir));
cd(folder)
end
