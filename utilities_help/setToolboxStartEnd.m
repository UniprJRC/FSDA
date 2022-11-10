function setToolboxStartEnd(toolboxFile, varargin)
%setToolboxStartEnd sets release compatibility in ToolboxPackagingConfiguration.prj file
%
%<a href="matlab: docsearchFS('setToolboxStartEnd')">Link to the help function</a>
%
% Required input arguments:
%
%    toolboxFile:  Toolbox packaging configuration file. Character or String.
%                  Absolute or relative path of the MATLAB file which
%                  contains the Toolbox Packaging Project. Note that this
%                  file must have a .prj extension.
%                  Example - ['tmp' filesep 'ToolboxPackagingConfiguration.prj']
%                  Data Types - Character or string scalar.
%
% Optional input arguments:
%
%   startVersion: Oldest release. Character or String or empty value (default).
%                 Oldest release for which compatibility is ensured in the
%                 format R[4digitsyear][a|b].
%                 For example if startVersion is 'R2017b' it means that
%                 compatibility is ensured from MATLAB 2017b. If this
%                 option is empty older release is automatically set to
%                 that of five years older than the currect release. For
%                 example if current release is R2021b startVersion is
%                 automatically set to 'R2016b'.
%               Example - 'startVersion','R2018a'
%               Data Types - Character or string scalar.
%
%   endVersion:   Newest releases. Character or String or empty value (default).
%                 Newest release for which compatibility is ensured in the
%                 format R[4digitsyear][a|b].
%                 For example if Newest is 'R2021b' it means that
%                 compatibility is ensured up to MATLAB 2021b. If this
%                 option is empty newest release is automatically set to
%                 currect release. For example if current release is R2021b endVersion is
%                 automatically set to 'R2021b'.
%               Example - 'endVersion','R2021b'
%               Data Types - Character or string scalar.
%
%
%
% Output:
%
%
% See also:  CreateFSDAtoolboxFile.m, CreateFSDAhelpFiles.m
%
% References:
%
% Acknowledgements:
%
% We would like to thank Jos Martin and the packaging team of Mathworks.
% Note that this routine is java based.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('setToolboxStartEnd')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %%  Call of setToolboxStartEnd with all the default options.
    FullPath=which('addFSDA2path');
    FSDAroot=fileparts(FullPath);
    fsep=filesep;
    fname='ToolboxPackagingConfiguration.prj';
    copyfile([FSDAroot fsep 'utilities_help' fsep fname],pwd)
    disp('Show release start and release end before calling setToolboxStartEnd')
    fileID = fopen(fname, 'r+');
    % Insert the file into fstring
    fstring=fscanf(fileID,'%c');
    disp(regexp(fstring,'<param.release.start>.*</param.release.start>','match'));
    disp(regexp(fstring,'<param.release.end>.*</param.release.end>','match'));
    fclose(fileID);
    % Set oldest and newest version inside ToolboxPackagingConfiguration.prj
    % In the case the oldest version is 5 years before the currect version and
    % end version is the current version on which MATLAB is running.
    setToolboxStartEnd(fname)
    disp('Show release start and release end after calling setToolboxStartEnd')
    fileID = fopen(fname, 'r+');
    % Insert the file into fstring
    fstring=fscanf(fileID,'%c');
    disp(regexp(fstring,'<param.release.start>.*</param.release.start>','match'));
    disp(regexp(fstring,'<param.release.end>.*</param.release.end>','match'));
    fclose(fileID);
%}

%{
    %% Call of setToolboxStartEnd with name/pairs.
    FullPath=which('addFSDA2path');
    FSDAroot=fileparts(FullPath);
    fsep=filesep;
    fname='ToolboxPackagingConfiguration.prj';
    copyfile([FSDAroot fsep 'utilities_help' fsep fname],pwd)
    % Set oldest and newest version inside ToolboxPackagingConfiguration.prj
    startVersion='R2012b';
    endVersion='R2020a';
    
    disp('Show release start and release end before calling setToolboxStartEnd')
    fileID = fopen(fname, 'r+');
    % Insert the file into fstring
    fstring=fscanf(fileID,'%c');
    disp(regexp(fstring,'<param.release.start>.*</param.release.start>','match'));
    disp(regexp(fstring,'<param.release.end>.*</param.release.end>','match'));
    fclose(fileID);
    % Set oldest and newest version inside ToolboxPackagingConfiguration.prj
    setToolboxStartEnd(fname, 'startVersion',startVersion,'endVersion',endVersion)
    disp('Show release start and release end after calling setToolboxStartEnd')
    fileID = fopen(fname, 'r+');
    % Insert the file into fstring
    fstring=fscanf(fileID,'%c');
    disp(regexp(fstring,'<param.release.start>.*</param.release.start>','match'));
    disp(regexp(fstring,'<param.release.end>.*</param.release.end>','match'));
    fclose(fileID);
%}

%% Beginning of code

% Verify the file exists
validateattributes(toolboxFile,{'char','string'},{'scalartext'}, ...
    'matlab.addons.toolbox.toolboxStartEnd','ToolboxFile',1)
toolboxFile = char(toolboxFile);
if exist(toolboxFile, 'file') ~= 2
    error(message('FSDA:setToolboxStartEnd:ToolboxFileNotFound',toolboxFile));
end

startVersion='';
endVersion='';
if nargin>1
    options=struct('startVersion',startVersion,'endVersion',endVersion);

    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:setToolboxStartEnd:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)

        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end

    end
    % set the options chosen by the user
    startVersion=options.startVersion;
    endVersion=options.endVersion;
end

if isempty(endVersion)
    endVersion=['R' version('-release')];
end

if isempty(startVersion)
    oldestVersion =num2str(str2num(endVersion(2:5))-5);
    startVersion=endVersion;
    startVersion(2:5)=oldestVersion;
end

% Validate oldestVersion & startVersion
validateattributes(startVersion,{'char','string'},{'scalartext'}, ...
    'startVersion','Version')
validateattributes(endVersion,{'char','string'},{'scalartext'}, ...
    'startVersion','Version')
cstartVersion=char(startVersion);
cendVersion=char(endVersion);

if length(cstartVersion)~=6 ||  cstartVersion(1)~=82
    error('FSDA:setToolboxStartEnd:WrongInputOpt',['startVersion must be in the format R[4digitsyear][a|b]...' ...
        ' For example ''R2020a''']);
elseif length(cendVersion)~=6 ||  cendVersion(1)~=82
    error('FSDA:setToolboxStartEnd:WrongInputOpt',['EndVersion must be in the format R[4digitsyear][a|b]...' ...
        ' For example ''R2020a''']);
else
end

% Get the absolute path to the file in case it was input as a relative path
if ~java.io.File(toolboxFile).isAbsolute
    toolboxFile = fullfile(pwd,toolboxFile);
end

[~,~,ext] = fileparts(toolboxFile);
switch lower(ext)
    case '.prj' % the input must be a .prj file
        service = com.mathworks.toolbox_packaging.services.ToolboxPackagingService;

        % Load the project
        try
            configKey = service.openProject(toolboxFile);
            c = onCleanup(@()service.closeProject(configKey));
        catch e
            error(message('MATLAB:toolbox_packaging:packaging:InvalidToolboxProjectFile',toolboxFile));
        end

        % Fetch the current version
        config = service.getConfiguration(configKey);

        % Set Oldest and Newest version
        try
            config.setParamAsString('param.release.start', startVersion);
            config.setParamAsString('param.release.end', endVersion);
        catch e
            error('FSDA:setToolboxStartEnd:NotRun','Unable to set start and end values');
        end
        if (~service.save(configKey))
            error(message('MATLAB:toolbox_packaging:packaging:SaveFailed', toolboxFile));
        end
    otherwise
        error(message('MATLAB:toolbox_packaging:packaging:InvalidFile',toolboxFile));
end

end

%FScategory:UTIHELP