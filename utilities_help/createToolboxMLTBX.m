function createToolboxMLTBX(prjFileName, toolboxVersion)
%Package toolbox as MLTBX file.
%   createToolboxMLTBX(toolboxVersion) builds the MLTBX file and saves it in the 
%   release folder. Input toolboxVersion is a string of the form Major.Minor.Bug.Build.

%   Copyright 2023 The MathWorks, Inc.

if isfile(prjFileName)
    packagingData = matlab.addons.toolbox.ToolboxOptions(prjFileName);
else
    msg="Unable to find file " + "'" + prjFileName+ "'";
    error(msg);
end

% Check if toolboxVersion is a string
% if ~isstring(toolboxVersion)
%     error("Toolbox version must be a string.");
% end

% Update the version number
packagingData.ToolboxVersion = toolboxVersion;

% Name for MLTBX file
packagingData.OutputFile = fullfile("release", "fsda.mltbx");

% Create toolbox MLTBX
matlab.addons.toolbox.packageToolbox(packagingData);