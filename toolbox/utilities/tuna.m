function tuna(varargin)
%tuna Toolbox Update Notification Alert
%
%<a href="matlab: docsearchFS('tuna')">Link to the help function</a>
%
%
% This function:
%   1. Finds the installed version of the toolbox via Add-On Explorer;
%   2. Queries the latest release from GitHub;
%   3. Compares versions and, if the user has installed the latest version
%    simply shows a message in the "Command Window". On the other hand,
%   if the user has not installed the latest version notifies with a dialog
%   and a message, and with a message in the Command Window.
%   If no input argument is given tuna is applied to FSDA
%   and checks whether the user has installed the latest version of FSDA,
%
%  Required input arguments:
%
%
%  Optional input arguments:
%
%
%        toolboxName: name of toolbox. character vector or string.
%                       The name of the toolbox which appears in the first
%                       column of the table which is
%                       queried using matlab.addons.installedAddons.
%                       if this positinal argument is not specified or if
%                       it is empty we use 'FSDA'
%                       Example - 'PIVlab'
%                       Data Types - char or string scalar or empty.
%
%       gitHubOwner:   Owner of the GitHub repo. character vector or string.
%                      GitHub username or organization name that owns the
%                      repository. For example in the case of the web address
%                      https://github.com/uniprJRC/FSDA the fitHubOwner is
%                      'uniprJRC'
%                       Example - 'Shrediquette'
%                       Data Types - char or string scalar or empty.
%
%       gitHubRepo:   repository name. character vector or string.
%                      This is the second part of the GitHub URL.
%                      For example in the case of the web address
%                      https://github.com/uniprJRC/FSDA the gitHubRepo is
%                      'FSDA'
%                       Example - 'PIVlab'
%                       Data Types - char or string scalar or empty.
%
%           gui:   show a modal dialog. Boolean. If gui is true (default) and
%                      the installed version is not the latest, a modal dialog
%                      notifies the user and the installed/latest versions are
%                      echoed to the Command Window. If gui is false no dialog is
%                      shown: nothing is displayed when the toolbox is up to date,
%                      an update notice is written to the Command Window when it is
%                      not, and a failure to reach GitHub is non-fatal. Use
%                      gui=false for unattended checks (e.g. at session startup).
%                       Example - 'gui', false
%                       Data Types - logical
%
%  Output:
%
% See also: matlab.addons.installedAddons, uialert
%
% References:
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tuna')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    %% tuna with no input argument.
    % Check whether the latest version of FSDA is installed.
    tuna()
%}

%{
    % tuna with just two positional optional input arguments.
    tuna('PIVlab','Shrediquette')
%}

%{
    % tuna with three positional optional input arguments.
    tuna('PIVlab','Shrediquette','PIVlab')
%}

%{
    % tuna with three positional optional input arguments passed as strings.
    tuna("FSDA","uniprJRC","FSDA")
%}

%{
    % tuna in quiet mode (no modal dialog), suitable for session startup.
    % Nothing is shown if FSDA is up to date; a Command Window notice appears
    % otherwise, and a failure to reach GitHub is non-fatal.
    tuna('FSDA','uniprJRC','FSDA','gui',false)
%}

%% Beginning of code

% Extract the optional 'gui' name/value pair (default true) before parsing the
% positional arguments, so every existing positional call syntax is unchanged.
% gui=false is the quiet/startup mode: no modal dialog, nothing displayed when
% the toolbox is up to date, a Command Window notice when it is not, and a
% GitHub connection failure is non-fatal.
gui = true;
idx = find(cellfun(@(x) (ischar(x) || isstring(x)) && strcmpi(x, 'gui'), varargin), 1);
if ~isempty(idx)
    if idx == numel(varargin)
        error('FSDA:tuna','The ''gui'' option requires a logical value.');
    end
    gui = logical(varargin{idx+1});
    varargin([idx idx+1]) = [];
end

npos = numel(varargin);
if npos==0
    toolboxName='FSDA';
    gitHubOwner='uniprJRC';
    gitHubRepo='FSDA';
elseif npos==1
    error('FSDA:tuna','At least two input arguments must be given.');
elseif npos==2
    toolboxName=varargin{1};
    gitHubOwner=varargin{2};
    gitHubRepo=toolboxName;
elseif npos==3
    toolboxName=varargin{1};
    gitHubOwner=varargin{2};
    gitHubRepo=varargin{3};
else
    error('FSDA:tuna','Too many positional input arguments, they do not have to exceed 3.');
end



% Get installed add-ons
addons = matlab.addons.installedAddons;
row = strcmp(addons.Name, toolboxName);

if ~any(row)
    warning('Toolbox "%s" is not installed.', toolboxName);
    return;
    % error('FSDA:tuna:wrgToolboxName','Toolbox "%s" is not installed', toolboxName);
end

installedVersion = addons.Version(row);
id=addons.Identifier(row);
if gui
    fprintf('Installed version of "%s": %s\n', toolboxName, installedVersion{1});
end

% Query GitHub API for latest release
apiURL = sprintf('https://api.github.com/repos/%s/%s/releases/latest', ...
    gitHubOwner, gitHubRepo);
try
    if gui
        options = weboptions('ContentType', 'json', 'Timeout', 15);
    else
        options = weboptions('ContentType', 'json', 'Timeout', 5);
    end
    data = webread(apiURL, options);
catch ME
    if gui
        disp(['web site <a href="' apiURL '">"' apiURL '"</a> not reachable']);
        error('FSDA:tuna:GitHubUnreachable','Could not query GitHub API: %s', ME.message);
    else
        % Quiet mode: a connection failure must not break the session.
        warning('FSDA:tuna:GitHubUnreachable','Could not query GitHub API: %s', ME.message);
        return;
    end
end

if isfield(data, 'tag_name')
    latestVersion = data.tag_name;
else
    error('FSDA:tuna:GitHubNoTag','No tag_name found in GitHub API response.');
end

if gui
    fprintf('Latest available version: %s\n', latestVersion);
end

% Compare versions
if isUpdateAvailable(installedVersion{1}, latestVersion)

    Extra=['or click on the uninstall link below\n'....
        '\n'....
        '\n'....
        '\n'....
        '\n'];

    msg = sprintf(['A new version of "%s" is available!\n\n' ...
        'Installed: %s\nLatest: %s\n\n' ...
        'Remove current version from Add-Ons|Manage Add Ons\n' ...
        'before installing the new version\n' ...
        Extra ...
        'To install the new version from the HOME tab\n' ...
        '"Add-Ons|Explore Add-Ons".\n' ...
        'In the search textbox type: "%s"'], ...
        toolboxName, installedVersion{1}, latestVersion, ...
        toolboxName);

    if gui
        % Display popup
        try
            customAlert(msg, 'Toolbox Update Available',toolboxName,installedVersion{1},id)
        catch
            fprintf(2, '%s\n', msg); % Fallback to Command Window
        end
    else
        % Quiet mode: notify in the Command Window, no modal dialog.
        fprintf(2, '%s\n', msg);
    end
else
    if gui
        fprintf('Toolbox %s is up-to-date.\n',toolboxName);
    end
end
end

function flag = isUpdateAvailable(installed, latest)
% Normalize versions: remove leading 'v' and split into numbers
normalize = @(v) str2double(strsplit(regexprep(v, '^v', ''), '.'));
inst = normalize(installed);
lat  = normalize(latest);

% Pad with zeros to same length
len = max(numel(inst), numel(lat));
inst(end+1:len) = 0;
lat(end+1:len)  = 0;

% Compare sequentially
flag = false;
for i = 1:len
    if lat(i) > inst(i)
        flag = true;
        return;
    elseif lat(i) < inst(i)
        return;
    end
end
end

function customAlert(message, titleText,toolboxName,installedVersion,currentID)
% Create a small uifigure with no decorations
fig = uifigure('Name', titleText, 'Position', [100 100 400 350], ...
    'Resize', 'on', 'WindowStyle', 'modal');

% Remove the toolbar and menu
fig.MenuBar = 'none';
fig.ToolBar = 'none';

% Message label | [left bottom width height]
uilabel(fig,'Text', message,'Position', [20 50 360 330], ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 14, 'WordWrap', 'off');

pnl = uipanel(fig, ...
    'Position', [60 160 335 35], ...
    'Scrollable', 'off');
yPos=20;
uilabel(pnl, ...
    'Text', sprintf('%s (v%s)', toolboxName, installedVersion), ...
    'Position', [10 yPos-10 280 22]);
% The Hyperlink component
hl = uihyperlink(pnl, ...
    'Text', 'Uninstall', ...
    'FontColor', [0.85 0.33 0.1], ...
    'Position', [250 yPos-10 80 22]);

% This is the functional equivalent of the matlab: protocol.
% We use an anonymous function to pass the ID to the callback.
hl.HyperlinkClickedFcn = @(src, event) uninstallAction(fig, toolboxName, currentID,installedVersion);


% OK button
uibutton(fig, 'Text', 'OK','Position', [150 20 100 30], ...
    'ButtonPushedFcn', @(src, event) close(fig));
% Set warning icon (default MATLAB icon)
warnIcon = fullfile(matlabroot,'toolbox','matlab','icons','warning.gif');
% btn.Icon = warnIcon;
uiimage(fig, 'ImageSource', warnIcon, 'Position', [10 160 50 50]);
end


% Callback function for the link
function uninstallAction(fig, name, id,installedVersion)
msg = sprintf('Are you sure you want to uninstall\n "%s" version "%s"?', name, installedVersion);
selection = uiconfirm(fig, msg, 'Confirm Uninstall', 'Icon', 'warning');

if strcmp(selection, 'OK')
    try
        fprintf('Uninstalling: %s...\n', name);
        matlab.addons.uninstall(id);
        uialert(fig, 'Add-on uninstalled. Now manually install the new version.', 'Success', 'CloseFcn', @(src, event) close(fig));
    catch ME
        uialert(fig, ME.message, 'Uninstall Failed');
    end
end
end
