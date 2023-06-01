function plan = buildfile

plan = buildplan(localfunctions);
% Make the "createMLTBX" task the default task in the plan
plan.DefaultTasks = "createMLTBX";
end

function createMLTBXtask(~)

cdir = pwd;
cleanup = onCleanup(@() cd(cdir));
cd('utilities_help')
CreateFSDAtoolboxInGithubAction

end
