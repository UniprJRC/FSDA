function list = findFile(root,varargin)
% findFile finds recursively all files in root.
%
% <a href="matlab: docsearchFS('findFile')">Link to the help function</a>
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
%
%   ExclDir:    Exclude directory pattern(s). String or cell arrays of
%               strings. User can use wildcards. Do not
%               use regular expression. Examples: 'abc' and 'ab*de*'. Use ''
%               or {} to skip this stage. File separator (i.e. '\' in
%               Windows or '/' in Unix) is not allowed. Case-sensitive. Default: ExclDir={''}.
%               Example - 'ExclDir','dirname'
%               Data Types - (cell array of) string
%
%   InclFiles:  Include file pattern(s). String or cell arrays of
%               strings. User can use wildcards. Do not use regular expression.
%               Use '*' to include all files. Note that '*' and '*.*' give
%               the same behaviour and return all files. File separator
%               (i.e. '\' in Windows or '/' in Unix) is not allowed.
%               Case-sensitive. Default: InclFile={'*'}.
%               Example - 'InclFiles','filename'
%               Data Types - (cell array of) string
%
%   ExclFiles:  Exclude file pattern(s). String or cell arrays of
%               strings. User can use wildcards. Do not use regular expression.
%               Use '' or {} to skip this check. File separator (i.e. '\'
%               in Windows or '/' in Unix) is not allowed.
%               Case-sensitive. Default: ExclFiles={''}.
%               Example - 'ExclFiles','filename'
%               Data Types - (cell array of) string
%
% Output:
%
%   list:       All files. Cell array of strings. List of all
%               files matched under root directory, sorted in alphabetical
%               and ascending order. Always return absolute path.
%
%
% See also: findDir
%
% References:
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('findFile')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% findFile with all default options.
    %find the location of findFile.m
    FullPath=which('findFile');
    %extract the directory containing findFile.m
    root=fileparts(FullPath);
    %list the content of the directory containing findFile.m
    list = findFile(root)
%}
%
%{
    % findFile with optional arguments.
    % Make sure that the FSDA path is loaded (call routine addFSDA2path).
    addFSDA2path
    % Find the location of findFile.m
    FullPath=which('addFSDA2path');
    % extract the root directory of FSDA
    root=fileparts(FullPath);
    %list the content of the directory under FSDA named 'graphics'
    list = findFile(root,'InclDir','graphics')
%}

%{
    % List the content of the directory under FSDA named 'graphics'
    % and exclude all those which start with res.
    % find the location of findFile.m
    FullPath=which('findFile');
    % extract the root directory of FSDA
    root=FullPath(1:strfind(FullPath,'FSDA')+3);
    list = findFile(root,'InclDir','graphics','ExclFiles','res*')
%}

%{
    % find the location of help file gplotmatrix.html.
    pathdocroot=docroot;
    pathExtHelpFile=findFile(pathdocroot,'InclFiles','gplotmatrixee.html');
%}

%% Beginning of code

% Assign input arguments.
options=struct('InclDir',{''},'ExclDir',{''},'InclFiles',{'*'},'ExclFiles',{''});

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:findFile:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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

InclDir=options.InclDir;
ExclDir=options.ExclDir;
InclFiles=options.InclFiles;
ExclFiles=options.ExclFiles;

% Check root.
if ~ischar(root)
    error('FSDA:findFile:WrongInputOpt','root is not a string.');
end

if isempty(root)
    warning('FSDA:findFile:WrongInputOpt','root is empty.');
end

% Force InclFiles and ExclFiles to be a cell array of string.
if isempty(InclFiles)
    InclFiles = {};
else
    if ~iscell(InclFiles)
        InclFiles = {InclFiles};
    end
end
if isempty(ExclFiles)
    ExclFiles = {};
else
    if ~iscell(ExclFiles)
        ExclFiles = {ExclFiles};
    end
end

% Make sure that InclFiles and ExclFiles do not contain file separator.
c = strfind(InclFiles, filesep);
c = [c{:}];
if ~isempty(c)
    error('FSDA:findFile:WrongInputOpt','One of the InclFiles file patterns contains file separator.');
end
c = strfind(ExclFiles, filesep);
c = [c{:}];
if ~isempty(c)
    error('FSDA:findFile:WrongInputOpt','One of the ExclFiles patterns contains file separator.');
end


% Get all suddirectories under rootdir. Contain absolute paths.
alldirs = findDir(root,'InclDir',InclDir,'ExclDir',ExclDir);


% Find all files that match.
list=cell(10000,1);
ij=1;
for m = 1:length(alldirs)
    for i = 1:length(InclFiles)
        filelist = dir(fullfile(alldirs{m}, InclFiles{i}));
        filelist = filelist(~[filelist.isdir]);
        for k = 1:length(filelist)
            list{ij}=fullfile(alldirs{m}, filelist(k).name);
            ij=ij+1;
        end
    end
end
if ij>1
    list=list(1:ij-1);
end


% Remove files from match according to excludefile.
if ~isempty(ExclFiles)
    % booToDelete = Boolean vector which contains a true in
    % correspondence of the elements of list which have to be removed
    booToDelete = false(length(list),1);
    for i = 1:length(list)
        [~, name, ext] = fileparts(list{i});
        filename = [name, ext];
        c = regexp(filename, regexptranslate('wildcard', ExclFiles), 'match');
        c = [c{:}];
        if any(strcmp(filename, c))
            booToDelete(i)=true;
        end
    end
    list(booToDelete) = [];
end

% REMARK.
% Note that
% cellfun('isempty', list);
% is much faster than
% cellfun(@isempty,list)
% as documented in
% http://undocumentedmatlab.com/blog/cellfun-undocumented-performance-boost/

% Sort files in list, if list is not empty.
if sum(cellfun(@isempty,list))<length(list)
    list = sort(list);
else
    disp('No file which matches the criteria has been found')
    list=[];
end

end
%FScategory:UTIGEN