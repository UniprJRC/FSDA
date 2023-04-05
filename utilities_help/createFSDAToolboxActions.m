% This file is meant to run on a Linux node in GitHub Actions
prjFileName=".."+filesep+"/fsdaToolboxPackaging.prj";
toolboxOption=matlab.addons.toolbox.ToolboxOptions(prjFileName);
% toolboxOption.Identifier="fsdaToolbox.UParma";

% Need to edit the line below before every release
toolboxOption.ToolboxVersion="0.0.22";

% Name for MLTBX file
toolboxOption.OutputFile=".."+filesep+"release"+filesep+"fsda.mltbx";

%Create toolbox MLTBX
matlab.addons.toolbox.packageToolbox(toolboxOption);