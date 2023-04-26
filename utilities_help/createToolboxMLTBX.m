function createToolboxMLTBX(verNumber)
% This file is meant to run on a Linux node in GitHub Actions
% Execute the code from project root
% The code also assumes that there is a fsdaToolboxPackaging.prj at the
% project root.
% Input:
%   verNumber- string in major.minor.bugfix format
% Copyright 2023 The MathWorks, Inc.

prjFileName = "fsdaToolboxPackaging.prj";
toolboxOption = matlab.addons.toolbox.ToolboxOptions(prjFileName);

% Update the version number
toolboxOption.ToolboxVersion = verNumber;

% Name for MLTBX file
toolboxOption.OutputFile = fullfile("release", "fsda.mltbx");
% Create toolbox MLTBX
matlab.addons.toolbox.packageToolbox(toolboxOption);