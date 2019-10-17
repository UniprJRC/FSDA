function [out, Excluded]=makecontentsfileFS(varargin)
%makecontentsfileFS  extends Matlab function makecontentsfile
%
%<a href="matlab: docsearchFS('makecontentsfileFS')">Link to the help function</a>
%
%   makecontentsfileFS creates a file named contents which contains
%   information about selected files present in folder and/or subfolders
%
% Required input arguments:
%
% Optional input arguments:
%
%    dirpath:       path to use.
%                   String or cell array of strings.
%                   Absolute path of the folder(s) for which summary file must
%                   be created.
%                   If dirpath is not specified dirpath is the one of the
%                   current folder and it is found using MATLAB instruction
%                   pwd. If dirpath is a cell array of string then NameOutputFile
%                   is created for all entry of dirpath.
%                   Example - 'dirpath',pwd
%                   Data Types - string
%                   Remark: dirpath can be conveniently created with
%                   function findDir.m
%  NameOutputFile : Name of output file. String.
%                   String containing the name of the file which has to be
%                   created. If this option is not specified Contents.m is used.
%                   Example - 'NameOutputFile','Mycontents.m'
%                   Data Types - string
%    force  :       Overwrite existing file.
%                   Boolean.
%                   Boolean which specifies whether existing output file must
%                   be overwritten. If force is true (default) FileName is
%                   overwritten else just out cell is produced
%                   Example - 'force',false
%                   Data Types - string
% FilterOutFileName  :  filter files depending on their name. String.
%                   String which specifies which .m files do not have to be
%                   included inside NameOutputFile. All files whose name contains
%                   FilterOutFileName will not be included inside NameOutputFile.
%                   If this optional argument is not specified all files containing
%                   string [OlD] (lower case or uppercase) will not be
%                   listed inside NameOutputFile
%                   Example - 'FilterOutFileName','veryold'
%                   Data Types - string
%FilterFileContent :  filter .m files depending on their content. String.
%                   String which specifies the string files must contain
%                   to be included inside NameOutputFile.
%                   If this optional argument is not specified or it is empty
%                   all files .m file (whose name does not contain
%                   FilterOutFileName) will be listed inside
%                   NameOutputFile. In the final column of output cell dout, the letters after
%                   FilterFileContent will be shown. For example suppose
%                   that FilterFileContent is '%TAGFDSA' and that two
%                   files respectively contain the words '%TAGFDSAregression' and
%                   '%TAGFDSAmultivariate' then the words 'regression' and
%                   'multivariate' will be given inside last column of dout
%                   and inside NameOutputFile under column 'category'
%                   Example - 'FilterFileContent','%FScategory'
%                   Data Types - string
% printOutputCell : print out cell is a file. String.
%                   String contaning the name of the file of the current
%                   folder which contains the summary of all files present
%                   in findDir folders. The default value of
%                   printOutputCell is '' that is the overall content file
%                   is not created
%                   Example - 'printOutputCell','ContentsAll.m'
%                   Data Types - string
%
% Output:
%
%          out:   structured information of filtered files containing
%                 selected Tags. Cell.
%                 Cell of size r-times-9 containing detailed information about
%                 the files present in required folder or folders if
%                 dirpath is a cell array of strings.
%               The columns of dout contain the following information:
%               out(:,1)= name of the file (with extension);
%               out(:,2)= date (in local format);
%               out(:,3)= size of the files in bytes;
%               out(:,4)= boolean 1 if name is a folder and 0 if name is a file;
%               out(:,5)= Modification date as serial date number;
%               out(:,6)= matlab file name (without .m extension);
%               out(:,7)= file description (the so called H1 line of the file);
%               out(:,8)= string which identifies file category. File
%                           category is the string after
%                           'FilterFileContent' in each file.
%               out(:,9)= string which contains file path.
%               Remark: if the file contains:
%               FilterFileContentword1 in row 34;
%               FilterFileContentAnotherName in row 456;
%               FilterFileContentObama, in row 1243;
%               than category is associated to
%               the last instance which is found (in this example 'Obama').
%    Excluded:  structured information of m files not included.
%                 Cell.
%                 Cell of size r-times-9 containing detailed information about
%                 the files present in required folder or folders if
%                 dirpath is a cell array of strings but which have been
%                 excluded by the filters.
%               The columns of Excluded contain the following information:
%               Excluded(:,1)= name of the file (with extension);
%               Excluded(:,2)= date (in local format);
%               Excluded(:,3)= size of the files in bytes;
%               Excluded(:,4)= boolean 1 if name is a folder and 0 if name is a file;
%               Excluded(:,5)= Modification date as serial date number;
%               Excluded(:,6)= matlab file name (without .m extension);
%               Excluded(:,7)= file description (the so called H1 line of the file);
%               Excluded(:,8)= string which identifies file category. File
%                           category is the string after
%                           'FilterFileContent' in each file.
%               Excluded(:,9)= string which contains file path.
%
% See also: publishFunctionAlpha, publishFunctionCate, publishFS
%
% References:
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('makecontentsfileFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % makecontentsfileFS with all default options.
    out=makecontentsfileFS('force',false);
%}

%{
    %Just create output cell out but not the output file if it already exists.
    % Create personalized content file of current folder.
    out=makecontentsfileFS('force',false);
%}

%{
    % Create personalized content file of selected subfolders.
    % Create contents file for subfolder '(main root of FSDA)\utilities' and list
    % just the files which contain string '%FScategory:'
    FileName='addFSDA2path';
    FullPath=which(FileName);
    FSDAroot=FullPath(1:end-length(FileName)-3);
    dirpathsel=[FSDAroot filesep 'utilities'];
    %  Note that force is false and therefore just output cell is created.
    out=makecontentsfileFS('dirpath',dirpathsel,'FilterFileContent','%FScategory:','force',false)
%}

%{
    % Example of the use of options dirpath and  FilterFileContent.
    % Preliminary step: create a list of the subfolders which have to be
    % included, using routine findDir with options 'InclDir' and 'ExclDir'.
    % Find full path of the main root where FSDA is installed
    FileName='addFSDA2path';
    FullPath=which(FileName);
    root=FullPath(1:end-length(FileName)-3);
    InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat'};
    ExclDir={'privateFS'  'datasets'};
    % Create list of folders which must have a personalized contents file
    list = findDir(root,'InclDir',InclDir,'ExclDir',ExclDir)
    % Crete personalized contents file for main folder of FSDA
    % and required subfolders.
    [out,Excluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',false)
%}


%% Beginning of code
NameOutputFile='Contents.m'; % default MATLAB name of output file
dirpath=pwd;   % default path (current folder)
force=true;   % overwrite existing NameOutputFile
FilterOutFileName='old'; % Remove all files whose name contains string old
FilterFileContent=''; % Do not look into the content of the file
printOutputCell=''; % specifies if the overall output has to be written in a file

options=struct('NameOutputFile',NameOutputFile,'dirpath',dirpath,...
    'force',force,'FilterOutFileName',FilterOutFileName,...
    'FilterFileContent',FilterFileContent,'printOutputCell',printOutputCell);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    dirpath=options.dirpath;
    force=options.force;
    FilterOutFileName=options.FilterOutFileName;
    FilterFileContent=options.FilterFileContent;
    printOutputCell=options.printOutputCell;
    
    % Check if dirpath is a string or a cell array of strings
    if iscell(dirpath)
        for i=1:length(dirpath)
            assert(exist(dirpath{i},'dir')==7,['Supplied path ' dirpath{i} ' does not exist'])
        end
        ldirpath=length(dirpath);
    else
        assert(ischar(dirpath),'''supplied path'' should be a charater array (string)')
        % Check if dirpath exists
        assert(exist(dirpath,'dir')==7,['Supplied path ' dirpath ' does not exist'])
        ldirpath=1;
    end
else
    ldirpath=1;
end

% lineSep = character which defines line separator
lineSep = char(java.lang.System.getProperty('line.separator'));

% Initialize out with a large number of rows

out=cell(1000,9);
iout=1;
Excluded=out;
iExcluded=1;


% Find for all file in all required subfolder:
% maximal name length of a file
% maximal length of a file category
% maximal length of a file description
maxNameLenAll=0;
maxCategoryLenAll = 0;
maxDescriptionLenAll=0;


for j=1:ldirpath
    
    if iscell(dirpath)
        dirpathj=dirpath{j};
    else
        % Find all files which have extension .m
        dirpathj=dirpath;
    end
    
    d = dir([dirpathj filesep '*.m']);
    % Act just if in the folder there are .m files
    if ~isempty(d)
        if max(strcmp(fieldnames(d),'folder'))>0
            d=rmfield(d,'folder');
        end
        
        % if file contains the string FilterOutFileName then it is not
        % listed
        if ~isempty(FilterOutFileName)
            boo=cellfun('isempty',regexpi({d.name},FilterOutFileName));
            lExcl=sum(~boo);
            if lExcl>0
                Excluded(iExcluded:iExcluded+lExcl-1,1:5)=struct2cell(d(~boo))';
                iExcluded=iExcluded+lExcl;
            end
            d=d(boo);
        end
        
        % Sort files in alphabetical order
        [~,sortIndex] = sort(lower({d.name}));
        d = d(sortIndex);
        
        
        % For each particulat folder find:
        % maximal name length of a file
        % maximal length of a file category
        % maximal length of a file description
        % maxNameLen is linked to the maximal name length of a file
        maxNameLen = 0;
        % maxDescriptionLen is linked to the maximal length of a file category
        maxCategoryLen = 0;
        % maxDescriptionLen is linked to the maximal length of a file description
        maxDescriptionLen=0;
        
        
        killIndex = [];
        noContentsFlag = 1;
        for i = 1:length(d)
            % create mfilename (that is name of the file without .m extension) from
            % file name
            d(i).mfilename = regexprep(d(i).name,'\.m$','');
            if strcmp(d(i).name,NameOutputFile)
                % Special case: remove the NameOutputFile.m file from the list
                % NameOutputFile.m should not list itself.
                killIndex = i;
                noContentsFlag = 0;
            else
                disp(d(i).name)
                try
                    [description,category]=get_H1lineandCategory(d(i).name,FilterFileContent);
                catch
                    description='';
                    category='';
                    % error('FSDA:normBoxCox:WrongInputFile',['error while parsing file :' d(i).name])
                end
                d(i).description=description;
                d(i).category=category;
                % Store also the path
                d(i).path=dirpathj;
                
                % Check what is the name of the file with the maximum length
                maxNameLen = max(length(d(i).mfilename), maxNameLen);
                maxCategoryLen= max(length(d(i).category), maxCategoryLen);
                maxDescriptionLen= max(length(d(i).description), maxDescriptionLen);
            end
        end
        
        maxNameLenAll=max(maxNameLenAll,maxNameLen);
        maxCategoryLenAll = max(maxCategoryLenAll,maxCategoryLen);
        maxDescriptionLenAll=max(maxDescriptionLenAll,maxDescriptionLen);
        
        
        d(killIndex) = [];
        
        maxNameLenStr = num2str(maxNameLen);
        maxCategoryLenStr = num2str(maxCategoryLen);
        maxDescriptionLenStr = num2str(maxDescriptionLen);
        
        maxNameLenAllStr = num2str(maxNameLenAll);
        maxCategoryLenAllStr = num2str(maxCategoryLenAll);
        maxDescriptionLenAllStr = num2str(maxDescriptionLenAll);
       
        dout=struct2cell(d)';
        % dout(:,1)= name of the file (with extension)
        % dout(:,2)= date (in local format)
        % dout(:,3)= size ()
        % dout(:,4)= boolean (1 if a folder)
        % dout(:,5)= date (in numeric format)
        % dout(:,6)= matlab file name (wthout .m extension)
        % dout(:,7)= file description
        % dout(:,8)= string which identifies file category
        % dout(:,9)= file path
        
        if ~isempty(FilterFileContent)
            boo=~cellfun(@isempty,dout(:,8));
            % Extract just the rows of d for which category does not the string
            % specified inside input option FilterFileContent
            
            lExcl=sum(~boo);
            if lExcl>0
                Excluded(iExcluded:iExcluded+lExcl-1,:)=dout(~boo,:);
                iExcluded=iExcluded+lExcl;
            end
            dout=dout(boo,:);
            d=d(boo);
        end
    else
        dout='';
    end
    
    %Include inside dout output of j-th folder which has been analyzed
    if ~isempty(dout)
        out(iout:(iout+size(dout,1)-1),:)=dout;
        iout=iout+size(dout,1);
    else
    end
    
    % Output file is created if force is true or Outputfile does not exist
    % provided there is something to write. In other words, if dout is
    % empty this means that there are not .m file for which summary has to
    % be created.
    if (noContentsFlag || force == true) && ~isempty(dout)
        [fid,errMsg] = fopen([dirpathj filesep NameOutputFile],'w');
        if fid < 0
            error(message('MATLAB:filebrowser:MakeContentsFileOpenError', errMsg))
        end
        % nm = name of the folder which will be printed in the contents filte
        [~,nm] = fileparts(dirpathj);
        % Print in uppercase the name of the folder and then leave two lines
        fprintf(fid,'%% %s%s%%%s',upper(nm), lineSep, lineSep);
        % fprintf(fid,'%% Files%s', lineSep);
        fprintf(fid,'%% File names, description, category and date last modified%s', lineSep);
        
        fprintf(fid,['%%' lineSep]);
        
        fprintf(fid,['%%   %-' maxNameLenStr 's - %-' maxDescriptionLenStr 's - %-' maxCategoryLenStr 's- %s%s'], ...
            'Name', 'Description','Category', 'Date last modified', lineSep);
        fprintf(fid,'%%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%s', lineSep);
        
        % date of last modified file in required format
        date=datestr(cell2mat(dout(:,5)),'yyyy mmm dd');
        
        for i = 1:size(dout,1)
            fprintf(fid,['%%   %-' maxNameLenStr 's - %-' maxDescriptionLenStr 's - %-' maxCategoryLenStr 's- %s%s'], ...
                d(i).mfilename, d(i).description,d(i).category, date(i,:), lineSep);
        end
        fclose(fid);
    else
        warning(message('MATLAB:filebrowser:MakeContentsFileExists'))
    end
end

out=out(1:iout-1,:);
Excluded=Excluded(1:iExcluded-1,:);

% sort output in alphabetical order
[~,sortIndex] = sort(lower(out(:,6)));
out = out(sortIndex,:);

if ~isempty(printOutputCell)
    dirpathmain=pwd;
    if (noContentsFlag || force == true) && ~isempty(out)
        [fid,errMsg] = fopen([dirpathmain filesep printOutputCell],'w');
        if fid < 0
            error(message('MATLAB:filebrowser:MakeContentsFileOpenError', errMsg))
        end
        % nm = name of the folder which will be printed in the contents filte
        [~,nm] = fileparts(dirpathmain);
        % Print in uppercase the name of the folder and then leave two lines
        fprintf(fid,'%% %s%s%%%s',upper(nm), lineSep, lineSep);
        % fprintf(fid,'%% Files%s', lineSep);
        fprintf(fid,'%% File names, description, category and date last modified%s', lineSep);
        
        fprintf(fid,['%%' lineSep]);
        
        fprintf(fid,['%%   %-' maxNameLenAllStr 's - %-' maxDescriptionLenAllStr 's - %-' maxCategoryLenAllStr 's- %s%s'], ...
            'Name', 'Description','Category', 'Date last modified', lineSep);
        fprintf(fid,'%%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%s', lineSep);
        
        % date of last modified file in required format
        date=datestr(cell2mat(out(:,5)),'yyyy mmm dd');
        
        for i = 1:size(out,1)
            %            fprintf(fid,['%%   %-' maxNameLenStr 's - %-' maxDescriptionLenStr 's - %-' maxCategoryLenStr 's- %s%s'], ...
            %                d(i).mfilename, d(i).description,d(i).category, date(i,:), lineSep);
            fprintf(fid,['%%   %-' maxNameLenAllStr 's - %-' maxDescriptionLenAllStr 's - %-' maxCategoryLenAllStr 's- %s%s'], ...
                out{i,6}, out{i,7},out{i,8}, date(i,:), lineSep);
        end
        fclose(fid);
    else
        warning(message('MATLAB:filebrowser:MakeContentsFileExists'))
    end
    
end
end


function [H1line,category] = get_H1lineandCategory(filename,searchTag)
%GET H1 LINE and file category through input option searchTag

if nargin<2
    searchTag='';
end

[~,name,ext] = fileparts(filename);
H1line = ''; % default output
if strcmp(ext,'.m')
    fid = fopen(filename); % open file
    tline = fgetl(fid); % read first line
    while ischar(tline)
        k = strfind(tline,'%'); % find comment
        if ~isempty(k)
            k = k(1);
            ispercents = false(size(tline(k:end)));
            ispercents(strfind(tline(k:end),'%'))=true;
            start = k+find(~(isspace(tline(k:end)) | ispercents),1,'first')-1;
            if ~isempty(start)
                tline = tline(start:end); % remove leading space/percent
                IX = strfind(lower(tline),lower(name));
                if ~isempty(IX)
                    if IX(1)==1
                        tline = tline(length(name)+1:end); % remove function name
                    end
                    tline = strtrim(tline); % remove any leading/trailing space
                end
                H1line = tline;
                H1line = strtrim(H1line);
                if ~isempty(H1line)
                    if strcmp(H1line(end),'.') % remove trailing period
                        H1line = H1line(1:end-1);
                    end
                    H1line(1) = upper(H1line(1)); % capitalize first letter
                end
            end
            tline = -1; % set tline to numeric
        else
            tline = fgetl(fid); % read next line
        end
    end
    
    % now get category
    fstring=fscanf(fid,'%c');
    
    if ~isempty(searchTag)
        % Note that in order to avoid extracting all the files which at the
        % end of the searchTag symbol ' we added as condition '[^'']' that
        % is not as subsequent character symbol '
        FScatPos=regexp(fstring,[searchTag '[^'']']);
        fincatPos=regexp(fstring,'\s');
        
        if isempty(FScatPos)
            category='';
        else
            % check the position where search category terminates
            % If more than one instance is found than just take the final one
            fincatPos=fincatPos(fincatPos>FScatPos(end));
            
            if isempty(fincatPos)
                category=strtrim(fstring(FScatPos(end)+length(searchTag):end));
            else
                category=strtrim(fstring(FScatPos(end)+length(searchTag):fincatPos(1)));
            end
        end
    else
        category='';
    end
    fclose(fid);
end

end
%FScategory:UTIHELP