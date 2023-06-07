function buildDocSearchForToolboxPre22a
% Select the root folder
[FSDAroot, cleanup] = changeDirToRootWithCleanup; %#ok<ASGLU>
% build the doc search database
% builddocsearchdb(fullfile(FSDAroot, "helpfiles", "pointersHTML"))
com.mathworks.mlwidgets.help.HelpPrefs.addProduct(FSDAroot); %#ok<JAPIMATHWORKS>
com.mathworks.mlwidgets.help.customdoc.CustomToolboxIndexer.index(fullfile(pwd, "helpfiles", "pointersHTML")); %#ok<JAPIMATHWORKS>
ls helpfiles\pointersHTML\helpsearch-v3\
end