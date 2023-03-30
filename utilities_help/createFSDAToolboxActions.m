% This file is meant to run on GitHub Actions, the runner is always a Linux
prjFileName=".."+filesep+"/fsdaToolboxPackaging.prj";
toolboxOption=matlab.addons.toolbox.ToolboxOptions(prjFileName);
% toolboxOption.Identifier="fsdaToolbox.UParma";

% Need to edit the line below for every release
toolboxOption.ToolboxVersion="0.0.17";
toolboxOption.OutputFile=".."+filesep+"release"+filesep+"fsda.mltbx";

%Create toolbox MLTBX
matlab.addons.toolbox.packageToolbox(toolboxOption);