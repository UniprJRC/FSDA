function buildDocSearchForToolboxPre22a
% Change to the root folder so that we can then pick up all the correct
% files.
[FSDAroot, cleanup] = changeDirToRootWithCleanup; %#ok<ASGLU>
% build the doc search database - Normally we would just call the following
% line of code - however prior to R2022b this doesn't work in nodisplay
% MATLAB so we need to look at its internals and do the bits needed.

% builddocsearchdb(fullfile(FSDAroot, "helpfiles", "pointersHTML"))

% Firstly we need to ensure that the help preference system knows about the
% info.xml file that is in the FSDA root - this is what registers the
% existence of the toolbox with help search
com.mathworks.mlwidgets.help.HelpPrefs.addProduct(FSDAroot); %#ok<JAPIMATHWORKS>
% Then we need to index the help using the customer toolbox indexer - this
% is what you can see in builddocsearchdb
com.mathworks.mlwidgets.help.customdoc.CustomToolboxIndexer.index(fullfile(pwd, "helpfiles", "pointersHTML")); %#ok<JAPIMATHWORKS>
% Finally - lets show that this worked
if exist(fullfile(FSDAroot, 'helpfiles', 'pointersHTML', 'helpsearch-v3'), 'dir')
    disp('docsearch index built in helpsearch-v3 folder')
end
end