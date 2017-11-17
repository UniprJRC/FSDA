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
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('findDir')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
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
    % Example of the use of InclDir. find all subfolders inside matlab help
    % documentaion which contain the string optim
    pathdocroot=docroot;
    FoldersWithOptim=findFile(pathdocroot,'InclDir','*optim*');
%}
%
%{
    % findDir with optional arguments.
    FullPath=which('findDir');
    root=FullPath(1:end-length('findDir.m')-1);
    list = findDir(root,'InclDir','datasets')
%}

%{
    % findDir with optional arguments 'InclDir' and 'ExclDir'.
    FileName='addFSDA2path';
    FullPath=which(FileName);
    root=FullPath(1:end-length(FileName)-3);
    InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat'};
    ExclDir={'privateFS'  'datasets'};
    list = findDir(root,'InclDir',InclDir,'ExclDir',ExclDir)
%}

%% Beginning of code

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
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

% Force root to be an absolute path.
if ~ischar(root)
    error('FSDA:findDir:WrongInputOpt','root is not a string.');
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
        error('FSDA:findDir:WrongInputOpt','One of the InclDir patterns contains file separator.');
end
c = strfind(ExclDir, filesep);
c = [c{:}];
if ~isempty(c)
        error('FSDA:findDir:WrongInputOpt','One of the ExclDir patterns contains file separator.');
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
    % booToDelete = Boolean vector which contains a true in
    % correspondence of the elements of list which have to be removed
    booToDelete = false(length(dirs),1);
    
    for i = 1:length(dirs)
        c = regexp(dirs(i).name, regexptranslate('wildcard', ExclDir), 'match');
        c = [c{:}];
        if any(strcmp(dirs(i).name, c))
            booToDelete(i)=true;
        end
    end
    dirs(booToDelete) = [];
end

% Recursively descend through all directories.
list=cell(length(dirs)+1,1);
list{1}=root;
ij=2;
for i=1:length(dirs)
    dirname = dirs(i).name;
    fdir=findDir(fullfile(root, dirname),'InclDir','','ExclDir',ExclDir);
    if length(fdir)==1
        list{ij} =char(fdir) ;   % Take out InclDir in all recursive calls.
    else
        list(ij:ij+length(fdir)-1)=fdir;
    end
    ij=ij+length(fdir);
end

% list = {root};
% for i=1:length(dirs)
%     dirname = dirs(i).name;
%     list = [list, findDir(fullfile(root, dirname),'InclDir','','ExclDir',ExclDir)];   % Take out InclDir in all recursive calls.
% end


% Check InclDir.
if ~isempty(InclDir) && ~any(strcmp(InclDir, '*'))
    
    % Modify 'InclDir' by adding '\' at the beginning and at the end.
    for i = 1:length(InclDir)
        InclDir{i} = [filesep, InclDir{i}, filesep];
    end
    
    % Create 'p_new' by adding '\' at the end.
    p_new = cell(1, length(list));
    for i = 1:length(list)
        
        p_new{i} = [list{i}, filesep];
        
    end
    
    % Only return the path whose subdirectory matches InclDir.
    booToKeep = false(length(list),1);
    
    for i = 1:length(list)
        c = regexp(p_new{i}, regexptranslate('wildcard', InclDir), 'match');
        c = [c{:}];
        if ~isempty(c)
            booToKeep(i)=true;
        end
    end
    list=list(booToKeep);
end

% Sort the list
list = sort(list);
end
%FScategory:UTIGEN