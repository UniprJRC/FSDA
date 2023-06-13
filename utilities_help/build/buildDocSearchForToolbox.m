function buildDocSearchForToolbox
% Select the root folder
[FSDAroot, cleanup] = changeDirToRootWithCleanup; %#ok<ASGLU>
% build the doc search database
builddocsearchdb(fullfile(FSDAroot, "helpfiles", "pointersHTML"))
end