function mHelpGUI(FileName)
%mHelpGUI calls the GUI which enables to modfify the Help for FileName
% create structure out from input XML file

% Copyright 2008-2021.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

InputStructForHelpGUI=xmlreadFS(FileName); %#ok<NASGU>
[FSDAroot]=fileparts(which('docsearchFS.m'));

fsep=filesep;
fullpathFileName=[FSDAroot fsep 'helpfiles' fsep 'XML' fsep 'InputStructForHelpGUI'];
save(fullpathFileName,'InputStructForHelpGUI')
HelpGUI
end
