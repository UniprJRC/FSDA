function buildDocSearchForToolbox
% Select the root folder
[FSDAroot, cleanup] = changeDirToRootWithCleanup; %#ok<ASGLU>

rehash
drawnow
pause(3)
which info.xml

% build the doc search database
builddocsearchdb(fullfile(FSDAroot, "helpfiles", "pointersHTML"))
end