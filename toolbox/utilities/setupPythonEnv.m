function setupPythonEnv(varargin)
%setupPythonEnv configures a Python environment compatible with MATLAB.
%
%<a href="matlab: docsearchFS('setupPythonEnv')">Link to the help function</a>
%
% This function detects a Python installation compatible with the current
% MATLAB release, excludes unsupported system interpreters (such as the
% Python distributed with Xcode on macOS), configures pyenv accordingly,
% and optionally executes a pip command in the selected Python environment.
%
% Required input arguments:
%
%   None.
%
% Optional input arguments:
%
%   PipCommand : Pip command to be executed in the selected Python
%                environment. If PipCommand does not start with a standard
%                pip action such as 'install', 'uninstall', 'list' or
%                'show', and contains a single token only, it is
%                interpreted as a package name and automatically expanded
%                to 'install <package>'. In other words you can use both
%                'PipCommand','install yfinance' or 'PipCommand','yfinance'
%                Example - 'PipCommand','uninstall yfinance -y'
%                Data Types - char | string
 %
%   Verbose    : Logical scalar which controls the display of diagnostic
%                messages.
%                Default is true.
%                Example - 'Verbose',false
%                Data Types - logical
%
% Output:
%
%   This function does not return output arguments. It configures the
%   Python environment for the current MATLAB session and optionally
%   executes a pip command.
%
% See also: pyenv, pyrun, system
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('setupPythonEnv')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples:
%
%{
    % Just configure a compatible Python environment.
    setupPythonEnv
%}

%{
    % Install Python package yfinance in the selected environment.
    setupPythonEnv('PipCommand','install yfinance');
%}

%{
    % Equivalent compact syntax: if a single package name is supplied, it
    % is interpreted as an install request.
    setupPythonEnv('PipCommand','yfinance');
%}

%{
    % Show the list of installed packages.
    setupPythonEnv('PipCommand','list');
%}

%{
    % Show reduced output.
    setupPythonEnv('PipCommand','show yfinance','Verbose',false);
%}

%% Beginning of code

PipCommand = '';
Verbose    = true;

if nargin > 0

    options = struct('PipCommand',PipCommand,'Verbose',Verbose);
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));

    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:setupPythonEnv:WrongInputOpt', ...
                ['Number of supplied options is invalid. Probably values ' ...
                'for some parameters are missing.']);
        end
        aux.chkoptions(options,UserOptions);
    end

    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end

    PipCommand = options.PipCommand;
    Verbose    = options.Verbose;
end

% 1. Check operating system, detect a compatible Python interpreter and
% configure pyenv if needed.
try
    [setup_success, python_exe] = run_system_checks(Verbose);
catch ME
    error('FSDA:setupPythonEnv:SystemCheckError', ...
        'Critical error during system analysis: %s', ME.message);
end

if ~setup_success
    if Verbose
        disp('Operation aborted. Python requirements not met.');
    end
    return;
end

% 2. Optionally execute a pip command in the selected environment.
if ~isempty(PipCommand)

    PipCommand = strtrim(char(PipCommand));
    pipCmdPadded = [PipCommand ' '];

    isExplicitCmd = startsWith(pipCmdPadded,'install ')   || ...
                    startsWith(pipCmdPadded,'uninstall ') || ...
                    startsWith(pipCmdPadded,'list ')      || ...
                    startsWith(pipCmdPadded,'show ')      || ...
                    strcmp(PipCommand,'list')             || ...
                    startsWith(PipCommand,'--version');

    if ~isExplicitCmd && ~contains(PipCommand,' ')
        full_command = ['install ' PipCommand];
    else
        full_command = PipCommand;
    end

    sys_cmd = ['"' char(python_exe) '" -m pip ' full_command];

    if Verbose
        disp(['Internal pip execution: ' sys_cmd]);
    end

    [status, cmdout] = system(sys_cmd);

    if status == 0
        if Verbose && ~isempty(strtrim(cmdout))
            disp(cmdout);
        end
    else
        error('FSDA:setupPythonEnv:PipExecutionError', ...
            'Error during pip execution:\n%s', cmdout);
    end
end

end

% -------------------------------------------------------------------------
% Subfunctions
% -------------------------------------------------------------------------

function [success, valid_path] = run_system_checks(Verbose)

success = false;
valid_path = "";

% Extract current MATLAB release, for example '2024b'.
m_rel = version('-release');
pe = pyenv;

% Anonymous function to identify Python interpreters associated with Xcode.
is_xcode = @(p) contains(p,'Xcode.app') || strcmp(p,'/usr/bin/python3');

% If pyenv is already loaded, reuse it unless it points to the Xcode
% Python interpreter on macOS.
if pe.Executable ~= "" && pe.Status == "Loaded"
    if ismac && is_xcode(char(pe.Executable))
        error(['MATLAB has already loaded the Xcode Python interpreter. ' ...
            'Restart MATLAB to apply the bypass.']);
    else
        success = true;
        valid_path = char(pe.Executable);
        return;
    end
end

if Verbose
    disp('--- Starting Python Environment Analysis ---');
end

if ispc

    if Verbose
        disp('Operating System: Windows');
    end

    paths = get_system_paths('windows');
    [valid_path, python_exists] = find_compatible_python(paths,m_rel,false);

    if ~python_exists
        if Verbose
            disp('-> No Python installation detected on the system.');
            suggest_installation(m_rel);
        end
    elseif valid_path == ""
        if Verbose
            disp(['-> Python is installed, but the version is not ' ...
                'compatible with the MATLAB release.']);
            suggest_installation(m_rel);
        end
    else
        success = true;
    end

elseif ismac

    if Verbose
        disp('Operating System: macOS');
    end

    [xcode_status, ~] = system('xcode-select -p 2>/dev/null');
    xcode_installed = (xcode_status == 0);
    paths = get_system_paths('macos');

    if ~xcode_installed

        if Verbose
            disp('Xcode Status: Not installed.');
        end

        [valid_path, python_exists] = find_compatible_python(paths,m_rel,false);

        if ~python_exists
            if Verbose
                disp('-> Python Detection: No installation found.');
                suggest_installation(m_rel);
            end
        elseif valid_path == ""
            if Verbose
                disp('-> Python Detection: Versions found but not compatible with MATLAB.');
                suggest_installation(m_rel);
            end
        else
            success = true;
        end

    else

        if Verbose
            disp('Xcode Status: Installed (activating system dependencies bypass).');
        end

        
        [valid_path, python_exists] = find_compatible_python(paths,m_rel,true);

        if ~python_exists
            if Verbose
                disp('-> Python Detection: No valid non-Xcode installation found.');
                suggest_installation(m_rel);
            end
        elseif valid_path == ""
            if Verbose
                disp('-> Python Detection: Independent installations are not compatible.');
                suggest_installation(m_rel);
            end
        else
            success = true;
        end
    end

else
    error('FSDA:setupPythonEnv:UnsupportedOS', ...
        'Operating system not supported.');
end

if success && pe.Executable == ""
    pyenv('Version',valid_path);
    if Verbose
        disp(['-> Configuration completed. Executable in use: ' valid_path]);
    end
end

end

function paths = get_system_paths(os_type)

if strcmp(os_type,'windows')
    [~, cmdout] = system('where python python3 python3.13 2>nul');
    fallback = "";
else
    [~, cmdout] = system('/bin/zsh -lc "which -a python3 python3.13 2>/dev/null"');
    pythonEnv = pyenv;
    pythonUserPath=['\n' char(pythonEnv.Executable)];
    fallback = sprintf(['\n/opt/homebrew/bin/python3' ...
        pythonUserPath ...
        '\n/opt/homebrew/bin/python3.13' ...
        '\n/usr/local/bin/python3' ...
        '\n/Library/Frameworks/Python.framework/Versions/3.13/bin/python3.13']);
end

raw_string = string(cmdout) + string(fallback);
paths = unique(split(raw_string,newline),'stable');

end

function [valid_path, any_python_exists] = find_compatible_python(paths,m_rel,bypass_xcode)

valid_path = "";
any_python_exists = false;

for i = 1:length(paths)

    curr_str = strtrim(paths(i));

    if curr_str == "" || ~isfile(curr_str)
        continue
    end

    curr_char = char(curr_str);

    if bypass_xcode && ...
            (contains(curr_char,'Xcode.app') || strcmp(curr_char,'/usr/bin/python3'))
        continue
    end

    [status, ver_str] = system(['"' curr_char '" --version']);

    if status == 0
        any_python_exists = true;
        if check_version_compatibility(ver_str,m_rel)
            valid_path = curr_char;
            break
        end
    end
end

end

function is_compat = check_version_compatibility(ver_str,m_rel)

tokens = regexp(string(ver_str),'3\.\d+','match');

if isempty(tokens)
    is_compat = false;
    return
end

py_ver_short = tokens(1);

switch m_rel
    case {'2026a','2025b'}
        compat_list = ["3.10","3.11","3.12","3.13"];
    case {'2025a','2024b'}
        compat_list = ["3.9","3.10","3.11","3.12"];
    case {'2024a','2023b'}
        compat_list = ["3.9","3.10","3.11"];
    case {'2023a','2022b'}
        compat_list = ["3.8","3.9","3.10"];
    case {'2022a','2021b'}
        compat_list = ["3.7","3.8","3.9"];
    otherwise
        compat_list = [];
end

is_compat = ismember(py_ver_short,compat_list);

end

function suggest_installation(m_rel)

disp('---------------------------------------------------------');
disp(['ACTION REQUIRED: Install a Python version compatible with MATLAB R' m_rel]);
disp('Download: https://www.python.org/downloads/');

switch m_rel
    case {'2026a','2025b'}
        disp('Supported: 3.10, 3.11, 3.12, 3.13');
    case {'2025a','2024b'}
        disp('Supported: 3.9, 3.10, 3.11, 3.12');
    case {'2024a','2023b'}
        disp('Supported: 3.9, 3.10, 3.11');
    case {'2023a','2022b'}
        disp('Supported: 3.8, 3.9, 3.10');
    case {'2022a','2021b'}
        disp('Supported: 3.7, 3.8, 3.9');
    otherwise
        disp('MATLAB version not classified within the R2021b-R2026a range.');
end

disp('---------------------------------------------------------');

end