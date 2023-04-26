function createToolboxMLTBX(toolboxVersion)
%Package toolbox as MLTBX file
% This file is meant to run on a Linux node in GitHub Actions
% Execute the code from project root
% The code also assumes that there is a fsdaToolboxPackaging.prj at the
% project root.
% Input:
%   toolboxVersion- string in major.minor.bugfix format
% Copyright 2023 The MathWorks, Inc.

prjFileName = "fsdaToolboxPackaging.prj";
packagingData = matlab.addons.toolbox.ToolboxOptions(prjFileName);

% Update the version number
packagingData.ToolboxVersion = toolboxVersion;

% Name for MLTBX file
packagingData.OutputFile = fullfile("release", "fsda.mltbx");
% Create toolbox MLTBX
matlab.addons.toolbox.packageToolbox(packagingData);