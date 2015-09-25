function p = searchdir(varargin)

%
%     FUNCTION: searchdir - Search Directories Recursively.
%
%       SYNTAX: p = searchdir(root);
%               p = searchdir(root, include, exclude);
%
%  DESCRIPTION: Search directories recursively under the root directory.
%               Searches are performed in two stages:
%
%                   Stage 1: Recursively get all directories. 
%
%                   Stage 2: Exclude all directories (and their subdirectories)
%                            that match one of the patterns specified in
%                            'exclude'. If this stage is skipped, then no
%                            directory will be excluded by this stage. If
%                            exclude = '*', then this stage will return only 
%                            the root directory.
%
%                   Stage 3: Include only the directories (returned by Stage 2)
%                            that match of the patterns specified in 'include'.
%                            If this stage is skipped, then all directories
%                            (returned by Stage 2) will be returned in output
%                            cell array.
%
%               Example 1: Search all sub-directories under directory 'dlog'
%                               >> p = searchdir('dlog');
%
%               Example 2: Search all sub-directories under directory 'dlog' but
%                          exclude directory '.svn' and all subdirectories under
%                          '.svn'.
%                               >> p = searchdir('dlog', '', '.svn');
%                                                   OR
%                               >> p = searchdir('dlog', '*', '.svn');
%
%               Example 3: Search all sub-directories under directory 'dlog'.
%                          Exclude directory '.svn' and all subdirectories under
%                          '.svn'. Also only include directory 'asm' and its
%                          subdirectories. However 'dlog\asm\.svn' and its
%                          subdirectories will be excluded.
%                               >> p = searchdir('dlog', 'asm', '.svn');
%
%               Example 4: This is a special case where user does not want this
%                          function to transverse down into any subdirectory
%                          (therefore this function will run fast). In this
%                          case, this function will always return the root
%                          directory.
%                               >> p = searchdir('dlog', '', '*');
%                                                   OR
%                               >> p = searchdir('dlog', '*', '*');
%
%        INPUT: - root (string)
%                   Root directory. Can be relative path or absolute path. Use
%                   '.' for current directory. Case-sensitive.
%
%               - include (string or 1-D row/col cell array of strings)
%                   Include directory pattern(s). User can use wildcards. Do not
%                   use regular expression. Examples: 'abc' and 'ab*de*'. File
%                   separator (i.e. '\' in Windows or '/' in Unix) is not
%                   allowed. Case-sensitive. Note that if include = {''}, then
%                   there will be no directory name that matches '' and so
%                   output cell array p will be {}. This function will skip this
%                   stage if one of the following conditions is true:
%                       (1) include = ''
%                       (2) include = {}
%                       (3) include = '*'
%                       (4) include = a cell array and one of its elements is
%                           '*'.
%
%               - exclude (string or 1-D row/col cell array of strings)
%                   Exclude directory pattern(s). User can use wildcards. Do not
%                   use regular expression. Examples: 'abc' and 'ab*de*'. Use ''
%                   or {} to skip this stage. File separator (i.e. '\' in 
%                   Windows or '/' in Unix) is not allowed. Case-sensitive.
%
%       OUTPUT: - p (1-D row cell array of string)
%                   Cell array of all subdirectories under root directory.
%                   Always return absolute path. Sorted in alphabetical
%                   and ascending order.
%
%    $Revision: 1610 $
%
%        $Date: 2007-03-05 11:03:15 -0800 (Mon, 05 Mar 2007) $
%
%      $Author: khung $


%
% Assign input arguments.
%
switch nargin
case 1
    root = varargin{1};
    include = {};
    exclude = {};
case 3
    [root, include, exclude] = deal(varargin{:});
otherwise
    error('Invalid number of input arguments.');
end


%
% Force root to be an absolute path.
%
if ~ischar(root)
    error('root is not a string.');
end
tmp = pwd;
cd(root);
root = pwd;
cd(tmp);


%
% Force include and exclude to be 1-D cell array.
%
if isempty(include)
    include = {};
else
    if ~iscell(include)
        include = {include};
    end
end
if isempty(exclude)
    exclude = {};
else
    if ~iscell(exclude)
        exclude = {exclude};
    end
end


%
% Make sure that 'include' and 'exclude' do not contain file separator.
%
c = strfind(include, filesep);
c = [c{:}];
if ~isempty(c)
    error('One of the include patterns contains file separator.');
end
c = strfind(exclude, filesep);
c = [c{:}];
if ~isempty(c)
    error('One of the exclude patterns contains file separator.');
end


%
% Find subdirectories under root (non-recursive).
%
files = dir(root);
if isempty(files)
    p = {};
    return;
end
isdir = logical(cat(1,files.isdir));
dirs = files(isdir); % select only directory entries from the current listing


%
% Remove "." and ".." from dirs. 
%
dirs(strcmp('.', {dirs.name})) = [];
dirs(strcmp('..', {dirs.name})) = [];


%
% Remove directories from 'dirs' according to 'exclude'.
%
if ~isempty(exclude)
    midx = [];
    for n = 1:length(dirs)
        c = regexp(dirs(n).name, regexptranslate('wildcard', exclude), 'match');
        c = [c{:}];
        if any(strcmp(dirs(n).name, c))
            midx = [midx, n];
        end
    end
    dirs(midx) = [];
end


%
% Recursively descend through all directories.
%
p = {root};
for i=1:length(dirs)
   dirname = dirs(i).name;
   p = [p, searchdir(fullfile(root, dirname), '', exclude)];   % Take out include in all recursive calls.
end


%
% Check include.
%
if ~isempty(include) && ~any(strcmp(include, '*'))
    
    % Modify 'include' by adding '\' at the beginning and at the end.
    for n = 1:length(include)
        include{n} = [filesep, include{n}, filesep];
    end
    
    % Create 'p_new' by adding '\' at the end.
    p_new = cell(1, length(p));
    for n = 1:length(p)
        p_new{n} = [p{n}, filesep];
    end

    % Only return the path whose subdirectory matches include.
    midx = [];
    for n = 1:length(p)
        c = regexp(p_new{n}, regexptranslate('wildcard', include), 'match');
        c = [c{:}];
        if ~isempty(c)
            midx = [midx, n];
        end
    end
    p = p(midx);
    
end


%
% Sort p.
%
p = sort(p);


%
% Exit function.
%
return;


