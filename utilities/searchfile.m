function match = searchfile(varargin)

%
%     FUNCTION: searchfile - Search Files Recursively.
%
%       SYNTAX: match = searchfile(rootdir, includefile);
%               match = searchfile(rootdir, includefile, excludefile);
%               match = searchfile(rootdir, includedir, excludedir, includefile, excludefile);
%
%  DESCRIPTION: Search files recursively under root directory.
%
%               Example 1: Search all files under directory 'dlog'. This will
%                          return files with or without extension.
%                               >> match = searchfile('dlog', '*');
%
%               Example 2: Search all files under directory 'dlog' but exclude
%                          files with extension = *.svn-work or *.svn-base.
%                               >> match = searchfile('dlog', '*', {'*.svn-work', '*.svn-base'});
%
%               Example 3: Search all *.c files under directory 'dlog' and
%                          'sbmp'. Exclude the following directories: '.svn',
%                          'asm' and 'build'. Exclude generated files.
%                               >> match = searchfile('.', {'dlog', 'sbmp'}, {'.svn', 'asm', 'build'}, '*.c', '*gen*');
%
%               Example 4: Search all files under root directory 'dlog' without
%                          transversing into any of its subdirectories
%                          (therefore this function will run fast).
%                               >> match = searchfile('dlog', '', '*', '*', '');
%
%        INPUT: - rootdir (string)
%                   Root directory. Can be relative path or absolute path. Use
%                   '.' for current directory.
%
%               - includefile (string or 1-D row/col cell array of string)
%                   Include file pattern(s). Use wildcard (recognized by DIR).
%                   Do not use regular expression. Use '*' to include all files.
%                   Note that '*' and '*.*' give the same behaviour and return
%                   all files. File separator (i.e. '\' in Windows or '/' in
%                   Unix) is not allowed. Case-sensitive.
%
%               - excludefile (string or 1-D row/col cell array of string)
%                   Exclude file pattern(s). User can use wildcards. Do not
%                   use regular expression. Use '' or {} to skip this check.
%                   File separator (i.e. '\' in Windows or '/' in Unix) is not
%                   allowed. Case-sensitive.
%
%               - includedir (string or 1-D row/col cell array of strings)
%                   Include directory pattern(s). Refer to input argument
%                   'include' in searchdir.m for details.
%
%               - excludedir (string or 1-D row/col cell array of strings)
%                   Exclude directory pattern(s). Refer to input argument
%                   'exclude' in searchdir.m for details.
%
%       OUTPUT: - match (1-D row cell array of string)
%                   Cell array of all matched files. Always return absolute 
%                   path. Sorted in alphabetical and ascending order.
%
%    $Revision: 1600 $
%
%        $Date: 2007-03-02 10:31:20 -0800 (Fri, 02 Mar 2007) $
%
%      $Author: khung $


%
% Assign input arguments.
%
switch nargin
case 2
    [rootdir, includefile] = deal(varargin{:});
    includedir = '';
    excludedir = '';
    excludefile = '';
case 3
   [rootdir, includefile, excludefile] = deal(varargin{:});
    includedir = '';
    excludedir = '';   
case 5
    [rootdir, includedir, excludedir, includefile, excludefile] = deal(varargin{:});
otherwise
    error('Invalid number of input arguments.');
end


%
% Check rootdir.
%
if ~ischar(rootdir)
    error('rootdir is not a string.');
end


%
% Force includefile and excludefile to be a 1-D cell array of string.
%
if isempty(includefile)
    includefile = {};
else
    if ~iscell(includefile)
        includefile = {includefile};
    end
end
if isempty(excludefile)
    excludefile = {};
else
    if ~iscell(excludefile)
        excludefile = {excludefile};
    end
end


%
% Make sure that 'includefile' and 'excludefile' do not contain file separator.
%
c = strfind(includefile, filesep);
c = [c{:}];
if ~isempty(c)
    error('One of the include file patterns contains file separator.');
end
c = strfind(excludefile, filesep);
c = [c{:}];
if ~isempty(c)
    error('One of the exclude file patterns contains file separator.');
end


%
% Get all suddirectories under rootdir. Contain absolute paths.
%
alldirs = searchdir(rootdir, includedir, excludedir);


%
% Find all files that match.
%
match = {};
for m = 1:length(alldirs)
    for n = 1:length(includefile)
        filelist = dir(fullfile(alldirs{m}, includefile{n}));
        filelist = filelist(~[filelist.isdir]);
        for k = 1:length(filelist)
            match{end+1} = fullfile(alldirs{m}, filelist(k).name);
        end
    end
end
            

%
% Remove files from match according to excludefile.
%
if ~isempty(excludefile)
    midx = [];
    for n = 1:length(match)
        [pathstr, name, ext] = fileparts(match{n});
        filename = [name, ext];
        c = regexp(filename, regexptranslate('wildcard', excludefile), 'match');
        c = [c{:}];
        if any(strcmp(filename, c))
            midx = [midx, n];
        end
    end
    match(midx) = [];
end


%
% Sort matched files.
%
match = sort(match);
            
            
%
% Exit function.
%
return;





