function list = findDir(root,varargin)
% findDir finds recursively all directories in root. 
% 
% <a href="matlab: docsearchFS('findDir')">Link to the help function</a>
%
%
% Required input arguments:
%
%   root:       Root directory. String. Case sensitive string that can indicate
%               relative or absolute path. Use '.' for current directory.
%               Data Types - string
%
%
% Optional input arguments:
%
%   InclDir:    Include directory pattern(s). String or cell arrays of
%               strings. User can use wildcards. Do not use regular expression,
%               for examples 'abc' and 'ab*de*'. File separator (i.e. '\' in 
%               Windows or '/' in Unix) is not allowed. Case-sensitive.
%               This filter is skipped if InclDir='', InclDir={},
%               InclDir='*' or if one of the element of the cell array
%               InclDir is '*'. Default: InclDir={''}.
%               Example - 'InclDir','dirname'
%               Data Types - (cell array of) string
%   ExclDir:    Exclude directory pattern(s). String or cell arrays of
%               strings. User can use wildcards. Do not
%               use regular expression. Examples: 'abc' and 'ab*de*'. Use ''
%               or {} to skip this stage. File separator (i.e. '\' in 
%               Windows or '/' in Unix) is not allowed. Case-sensitive. Default: ExclDir={''}.
%               Example - 'ExclDir','dirname'
%               Data Types - (cell array of) string
%
% Output:
%
%   list:       All subdirectories. Cell array of strings. List of all 
%               subdirectories under root directory, sorted in alphabetical
%               and ascending order. Always return absolute path. 
%
%
% See also: findFile
%
% References:
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('findDir')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:
%
%{
    %% findDir with all default options.
    FullPath=which('findDir');
    root=FullPath(1:end-length('findDir.m')-1);
    list = findDir(root)
%}   
%
%{
    % findDir with otpional arguments.
    FullPath=which('findDir');
    root=FullPath(1:end-length('findDir.m')-1);
    list = findDir(root,'InclDir','datasets')
%}   

% Assign input arguments.
options=struct('InclDir',{''},'ExclDir',{''});

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSR:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

% Force root to be an absolute path.
if ~ischar(root)
    error('root is not a string.');
end
tmp = pwd;
cd(root);
root = pwd;
cd(tmp);

% Force InclDir and ExclDir to be cell array.
InclDir=options.InclDir;
ExclDir=options.ExclDir;
if isempty(InclDir)
    InclDir = {};
else
    if ~iscell(InclDir)
        InclDir = {InclDir};
    end
end
if isempty(ExclDir)
    ExclDir = {};
else
    if ~iscell(ExclDir)
        ExclDir = {ExclDir};
    end
end

% Make sure that 'InclDir' and 'ExclDir' do not contain file separator.
c = strfind(InclDir, filesep);
c = [c{:}];
if ~isempty(c)
    error('One of the InclDir patterns contains file separator.');
end
c = strfind(ExclDir, filesep);
c = [c{:}];
if ~isempty(c)
    error('One of the ExclDir patterns contains file separator.');
end

% Find subdirectories under root (non-recursive).
files = dir(root);
if isempty(files)
    list = {};
    return;
end
isdir = logical(cat(1,files.isdir));
dirs = files(isdir); % select only directory entries from the current listing

% Remove "." and ".." from dirs. 
dirs(strcmp('.', {dirs.name})) = [];
dirs(strcmp('..', {dirs.name})) = [];

% Remove directories from 'dirs' according to 'ExclDir'.
if ~isempty(ExclDir)
    midx = [];
    for n = 1:length(dirs)
        c = regexp(dirs(n).name, regexptranslate('wildcard', ExclDir), 'match');
        c = [c{:}];
        if any(strcmp(dirs(n).name, c))
            midx = [midx, n];
        end
    end
    dirs(midx) = [];
end

% Recursively descend through all directories.
list = {root};
for i=1:length(dirs)
   dirname = dirs(i).name;
   list = [list, findDir(fullfile(root, dirname),'InclDir','','ExclDir',ExclDir)];   % Take out InclDir in all recursive calls.
end

% Check InclDir.
if ~isempty(InclDir) && ~any(strcmp(InclDir, '*'))
    
    % Modify 'InclDir' by adding '\' at the beginning and at the end.
    for n = 1:length(InclDir)
        InclDir{n} = [filesep, InclDir{n}, filesep];
    end
    
    % Create 'p_new' by adding '\' at the end.
    p_new = cell(1, length(list));
    for n = 1:length(list)
        p_new{n} = [list{n}, filesep];
    end

    % Only return the path whose subdirectory matches InclDir.
    midx = [];
    for n = 1:length(list)
        c = regexp(p_new{n}, regexptranslate('wildcard', InclDir), 'match');
        c = [c{:}];
        if ~isempty(c)
            midx = [midx, n];
        end
    end
    list = list(midx);
end

% Sort the list
list = sort(list);