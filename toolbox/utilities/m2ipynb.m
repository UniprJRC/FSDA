function [Incl, Excluded]=m2ipynb(varargin)
%m2ipynb  convert selected m files into Jupyter notebook files
%
%<a href="matlab: docsearchFS('m2ipynb')">Link to the help function</a>
%
%   m2ipynb transforms m files which have a predefined label first
%   into mlx files and then into Jupiter notebook files. This file also
%   automatically appends inside README.md files a table written in markup
%   language. To understand how this table looks like see
%   https://github.com/UniprJRC/FigMonitoringBook/tree/main/cap1 or
%   https://github.com/UniprJRC/FigMonitoringBook/tree/main/cap2 or
%   https://github.com/UniprJRC/FigMonitoringBook/tree/main/cap3 or every
%   other chapter in the book
%
% Required input arguments:
%
% Optional input arguments:
%
%
%   append2README: append or not the list of filtered files to README.md
%                   file. Boolean. By default the list of filtered files in
%                   the specififed folder will be appended to the README.md
%                   file. If this table already exists inside the
%                   README.md, it will be replaced with the current one.
%                   Example - 'append2README',false
%                   Data Types - logical
%
% CatchError:       Whether to catch errors. Boolean.
%                   Whether to catch errors when running the live script or
%                   function during conversion, specified as a numeric or
%                   logical 1 (true) or 0 (false). If CatchError is true
%                   and an error occurs, export includes the error in the
%                   converted file. If CatchError is false and an error
%                   occurs, export displays the error in the Command Window
%                   and does not create a converted file.
%                   Example - 'CatchError',false
%                   Data Types - logical
%
%  category :       label inside the files which have to be translated into
%                   ipynb format. Charater or string. As default all .m
%                   files which contain the label '%InsideREADME' will be
%                   converted into mlx format and later into ipynb format.
%                   They will also be included  in the table inside the
%                   README.md file
%                   Example - 'category','##myPersonalLabel##'
%                   Data Types - character or string
%
%
%
% deleteMLXfiles : delete or not the .mlx files after their conversion to
%                   jupiter notebook format.
%                   Boolean. The default is false, that is .mlx are note
%                   deleted after their conversion to ipynb format.
%                   Example - 'deleteMLXfiles',true
%                   Data Types - logical
%
%    dirpath:       path to use or file to convert.
%                   Cell array of characters or character. Absolute path of
%                   the folder(s) for which m2ipynb files must be created.
%                   If dirpath is not specified or it is empty all .m files
%                   in the current folder with the category label will be
%                   converted. If dirpath is a cell array of characters
%                   then .m file are converted for all specified subfolders.
%                   If dirpath is a charater containing a single file
%                   in the current folder, just this file will be
%                   converted.
%                   Example - 'dirpath',pwd
%                   Data Types - cell array of characters or char
%                   Remark: dirpath can be conveniently created
%                   with function findDir.m
%
%
% FilterOutFileName  :  filter files depending on their name. Character or String.
%                   Character or string which specifies which .m files do not have to be
%                   included inside NameOutputFile. All files whose name contains
%                   FilterOutFileName will not be included inside NameOutputFile.
%                   If this optional argument is not specified all files containing
%                   string [OlD] (lower case or uppercase) will not be
%                   considered
%                   Example - 'FilterOutFileName','veryold'
%                   Data Types - string
%
%        msg  :     It controls whether to display or not messages on the
%                   screen. Boolean.
%                   If msg==false (default) messages are not displayed
%                   on the screen about all .m files which are considered for
%                   translation to ipynb
%                   Example - 'msg',false
%                   Data Types - logical
%
% printOutputCell : print output cell in the screen. Boolean.
%                   If printOutputCell
%                   folder which contains the summary of all files present
%                   in findDir folders. The default value of
%                   printOutputCell is '' that is the overall content file
%                   is not created
%                   Example - 'printOutputCell','ContentsAll.m'
%                   Data Types - string
%
%
%     repoName :    GitHub repository name. Character or string.
%                   String which the GitHub repository address
%                   The default is to use 'UniprJRC/FigMonitoringBook'
%                   Example - repoName 'github/awesome-matlab'
%                   Data Types - Character or string
%
%            run :  Whether to run the code in the .mlx files include
%                   outputs. Boolean.
%                   If this option is set to true  (default) the code
%                   inside the MLX files is run and the output is included
%                   in the converted file. Through option runExcluded it is
%                   possible to control the files which have to be
%                   translated to ipynb format with run set to false.
%                   Example - 'run',false
%                   Data Types - logical
%
%   runExcluded  :  String which identifies files with run set to false.
%                   'Character' or 'string'. The default of runExcluded is
%                   'Interactive', in the sense that all files whose name
%                   contains the string Interactive have to be translated to
%                   ipynb format with option run set to false.
%                   Example - 'runExcluded','UserInteraction'
%                   Data Types - logical
%
% Output:
%
%          Incl:   structured information of translated files Cell.
%                 Cell of size r-times-4 containing detailed information about
%                 the files present in required folder or subfolders for which
%                 jpynb format has been created.
%               The columns of Incl contain the following information:
%               Incl(:,1)= name of the m file (with extension);
%               Incl(:,2:3)= Information which has been added in the second column
%                          of the table inside README.md
%               out(:,4)= file paths;
%
%    Excluded:  structured information of m files not included.
%                 Cell.
%                 Cell of size r-times-5 containing detailed information about
%                 the files present in required folder or folders if
%                 dirpath is a cell array of strings but which have been
%                 excluded by the filters.
%               The columns of Excluded contain the following information:
%               Excluded(:,1)= name of the file (with extension);
%               Excluded(:,2)= date (in local format);
%               Excluded(:,3)= size of the files in bytes;
%               Excluded(:,4)= boolean 1 if name is a folder and 0 if name is a file;
%               Excluded(:,5)= Modification date as serial date number;
%
%
% See also: export, publishFunctionAlpha, publishFunctionCate, publishFS
%
% References:
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('m2ipynb')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % m2ipynb with all default options.
    % Convert first to .mlx and then to .ipynn all .m files in the current
    % folder which contain '%InsideREADME' and append the table to
    % README.md
    out=m2ipynb();
%}

%{
    % As before but do not run the mlx files.
    out=m2ipynb('run',false);
%}

%{
    % Example 1 of the use of option dirpath.
    % In this case dipath is a character associated with a single .m file
    % Just convert file whose name is 'rescaleFS.m'
    InclDir='rescaleFS.m';
    try
        out=m2ipynb('dirpath',InclDir);
    catch
        %disp(lasterror)
        emsg='Function m2ipynb has been called on a file which does not exist';
        MException('FSDA:m2ipynb:WrongInputOpt',emsg)
    end
%}

%{
    % Example 2 of the use of option dirpath.
    % In this case dipath is a cell which identifies a particular subfolder whose
    % .m files have to be converted to .ipynb.
    InclDir={'cap1'};
    try
        out=m2ipynb('dirpath',InclDir);
    catch
        %disp(lasterror)
        emsg='Function m2ipynb has been called on a subfolder which does not exist';
        MException('FSDA:m2ipynb:WrongInputOpt',emsg)
    end
%}

%{
    % Example 3 of the use of option dirpath.% Preliminary call to findDir
    % Include all subfolderw whose name starts with cap or whose name is
    % solutionsEX
    InclDir={'cap*' 'solutionsEX'};
    % Remove all subfolders whose names are 'AnalysisByDataset' or
    %  'IncomeDatasets' whose name starts with fig or with Fig5 or with
    %  exercises
    ExclDir = {'AnalysisByDataset' 'IncomeDatasets' 'fig*' 'Fig5*' 'exercises*'};
    % dirpath contains 
    dirpath=findDir(pwd,'InclDir',InclDir,'ExclDir',ExclDir);
    [Incl,Excl]=m2ipynb('run',false,'dirpath',dirpath);
%}

%% Beginning of code

dirpath=pwd;   % default path (if it is empty it is current folder)
FilterOutFileName='old'; % Do not translate files whose name contains string old
printOutputCell=false; % specifies whether to print on the screen the output cell
msg=true;
runMLXfile=true;
deleteMLXfiles=false;
append2README=true;
category='%InsideREADME';
repoName='UniprJRC/FigMonitoringBook';
CatchError=true;
runExcluded='Interactive';

options=struct('dirpath',dirpath,...
    'run',runMLXfile,'FilterOutFileName',FilterOutFileName,'msg',msg,...
    'category',category,'printOutputCell',printOutputCell, ...
    'repoName',repoName,'append2README',append2README', ...
    'deleteMLXfiles',deleteMLXfiles,'CatchError',CatchError, ...
    'runExcluded',runExcluded);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:m2ipynb:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    dirpath=options.dirpath;
    FilterOutFileName=options.FilterOutFileName;
    category=options.category;
    printOutputCell=options.printOutputCell;
    runMLXfile=options.run;
    msg=options.msg;
    deleteMLXfiles=options.deleteMLXfiles;
    append2README=options.append2README;
    CatchError=options.CatchError;
    runExcluded=options.runExcluded;

    Folder=true;
    % Check if dirpath is a cell array of strings or a char associated with
    % a single file or a single folder or an empty value
    if isempty(dirpath)
        ldirpath=1;
        dirpath=pwd;
    elseif iscell(dirpath)
        for i=1:length(dirpath)
            assert(exist(dirpath{i},'dir')==7,['Supplied path ' dirpath{i} ' does not exist'])
        end
        ldirpath=length(dirpath);
    else
        assert(ischar(dirpath),'''supplied path'' should be an empty value a  char or a cell array of chars')
        % Check if the char which is supplied is a file or a folder which
        % exists
        if exist(dirpath,'dir')==7
            % In this case dirpath is char identifying a subfolder which
            % exists
        elseif exist(dirpath,'file')==2
            % In this case dirpath is a file in the current folder
            Folder=false;
        else
            errmsg= ['Supplied path ' dirpath ' does not exist and/or finename not found in current folder'];
            error('FSDA:m2ipynb:WrongInputOpt',errmsg)
        end
        ldirpath=1;
    end
else
    ldirpath=1;
    Folder=true;
end

prerepo='[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=';


% Initialize Incl with a large number of rows
% List of the files which have been translated to ipynbb format
Incl=cell(1000,4);
iout=1;
Excluded=Incl;
iExcluded=1;
currentPath=pwd;

for j=1:ldirpath
    cd(currentPath);

    if iscell(dirpath)
        dirpathj=dirpath{j};
        d = dir([dirpathj filesep '*.m']);
    elseif Folder==true % isempty(dirpath)
        % Find all files which have extension .m in the current folder
        dirpathj=dirpath;
        d = dir([dirpathj filesep '*.m']);
    else
        % input is a single m file
        assert(isfile(dirpath),'If dirpath is a char it must an .m file in the current folder')
        d = dir(dirpath);
        dirpathj=pwd;
    end

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
        filel=lower({d.name});
        % If FileNames contain number an attempt will be made to order the
        % file according to numbers. For examples if FileNames are
        % ex5_2.m ex10_1.m ex10_1 cames later than ex5_2.m
        personOrder=true;
        if personOrder==true
            %names=extract(string(filel(4)),digitsPattern);
            %str2num(strjoin(names,''))
            names=zeros(length(d),2);
            for ii=1:length(d)
                namii=extract(string(filel(ii)),digitsPattern);
                if length(namii)>=2
                    names(ii,:)=double(namii(1:2));
                end
            end

            [~,sortIndex]=sortrows(names,[1 2]);
        else
            [~,sortIndex] = sort(lower({d.name}));
        end

        d = d(sortIndex);
        dout=cell(100,3);
        ij=1;
        cd(dirpathj)
        for i = 1:length(d)
            % create mfilename (that is name of the file without .m extension) from
            % file name
            d(i).mfilename = regexprep(d(i).name,'\.m$','');
            FileName=d(i).name;



            if strcmp(d(i).name,'m2ipynb.m')
                catBoolean=false;
            else
                try
                    [H1line,H2line,catBoolean]=get_H1lineandCategory(FileName,category);
                    H2line=regexprep(H2line,'\n',' ');
                catch
                    catBoolean=false;
                    % error('FSDA:m2ipynb:WrongInputFile',['error while parsing file :' d(i).name])
                end
            end
            % Convert file to mlx and subsequently to ipynb format
            if catBoolean ==true
                if msg==true
                    disp(['Converting file: ' FileName])
                end

                FileNameMLX=[FileName(1:end-2) '.mlx'];
                matlab.internal.liveeditor.openAndSave(FileName,FileNameMLX);
                %  try
                if contains(FileName,runExcluded)
                    export(FileNameMLX,'Format','ipynb','Run',false,'CatchError',CatchError);
                else
                    export(FileNameMLX,'Format','ipynb','Run',runMLXfile,'CatchError',CatchError);
                end
                FileNameipynb=[FileName(1:end-2) '.ipynb'];
                disp(['File '  num2str(FileNameipynb) ' created'])
                % catch
                %     warning('FSDA:m2ipynb:WrongInputFile',['CatchError in file ' FileName])
                %     % error('FSDA:m2ipynb:WrongInputFile','Source code error in original .m file')
                % disp(['CatchError in file ' FileName])
                % end
                dout{ij,1}=FileName;
                dout{ij,2}=H1line;
                dout{ij,3}=H2line;
                % Delete mlx files if option  deleteMLXfiles is true
                if deleteMLXfiles ==true
                    delete(FileNameMLX);
                end
                ij=ij+1;
            end
        end
    end

    dout=dout(1:ij-1,:);
    if ~isempty(dout)
        boo=~cellfun(@isempty,dout(:,1));
        % Extract just the rows of d for which category does not the string
        % specified inside input option FilterFileContent

        lExcl=sum(~boo);
        if lExcl>0
            Excluded(iExcluded:iExcluded+lExcl-1,:)=dout(~boo,:);
            iExcluded=iExcluded+lExcl;
        end
        dout=dout(boo,:);
        % d=d(boo);
    end

    % Below is the typical file which has to be automatically created
    % | FileName | Description | Open in MATLAB on line | Jupiter notebook |
    % |---|---|---|---|
    % |MentalIllness.m |Contaminated illness data.<br/> This file creates Figure 4.33.|[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/FigMonitoringBook&file=/cap4/MentalIllness.m)| [[ipynb](MentalIllness.ipynb)]
    % |Stars.m|Stars data.<br/> This file creates Figures 4.1-4.4, 4.9-4.11|[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/FigMonitoringBook&file=/cap4/Stars.m)| [[ipynb](Stars.ipynb)]
    % |SurgicalUnit.m|Surgical Unit data.<br/> This file creates Figures 4.30-4.33.|[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/FigMonitoringBook&file=/cap4/SurgicalUnit.m)| [[ipynb](SurgicalUnit.ipynb)]
    % |Wool.m|Wool data.<br/> This file creates Figures 4.5-4.8.|[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/FigMonitoringBook&file=/cap4/Wool.m)| [[ipynb](Wool.ipynb)]

    if append2README ==true
        % Create a table in markup language and append it at the end of the
        % README.md file
        fid = fopen('README.md'); % open file
        if fid==-1
            fstring1='';
        else
            fstring=fscanf(fid,'%c');
            fclose(fid);

            % findTBL=regexp(fstring,'(|*)\s*FileName\s*|\s*Description\s*|\s*Script\s*\s*Jupiter');
            findTBL=regexp(fstring,'(|*)\s*FileName\s*|\s*Description');

            if ~isempty(findTBL)
                fstring1=fstring(1:findTBL(1)-2);
            else
                fstring1=fstring;
            end
        end

        %Include inside dout output of j-th folder which has been analyzed
        if ~isempty(dout)

            % Create table
            TBL='| FileName | Description | Open in MATLAB on line | Jupiter notebook | \r |---|---|---|---| \r ';

            for i = 1:size(dout,1)

                % | Income1 |[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/FigMonitoringBook&file=/cap1/Income1main.m) |
                FileName=dout{i,1};
                posSep=regexp(dirpathj,filesep);
                folderName=[dirpathj(posSep(end)+1:end) '/'];
                % folderName='';
                postrepo=['&file=' folderName '/' FileName ')'];
                FileNameNoExt=FileName(1:end-2);
                folderName='';
                ipy=['| [[ipynb](' folderName FileNameNoExt '.ipynb)]'];
                Row2add=[prerepo repoName postrepo ipy];
                TBL=[TBL '|' FileName '|' dout{i,2} '<br/> ' dout{i,3} '|' Row2add '\r']; %#ok<AGROW>
                %  [[ipynb](/cap4/ARregression.ipynb)]

            end
            TBLmarkup=sprintf(TBL);
            file1ID=fopen([dirpathj filesep 'README.md'],'w');

            outstring=[fstring1 TBLmarkup];
            fprintf(file1ID,'%s',outstring);
            fclose('all');


            Incl(iout:(iout+size(dout,1)-1),1:3)=dout;
            Incl(iout:(iout+size(dout,1)-1),4)=repelem({dirpathj},size(dout,1),1);
            iout=iout+size(dout,1);

        end
    end
end
cd(currentPath);
Incl=Incl(1:iout-1,:);
Excluded=Excluded(1:iExcluded-1,:);

if printOutputCell==true
    disp(Incl)
end

end


function [H1line,H2line, catBoolean] = get_H1lineandCategory(filename,searchTag)
%GET H1 LINE and file category through input option searchTag

if nargin<2
    searchTag='';
end

[~,~,ext] = fileparts(filename);
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
                H1line = tline;
                H1line = strtrim(H1line);
                if ~isempty(H1line)
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

    findPercentage=strfind(fstring,'%%');
    fstringSel=fstring(1:findPercentage(1)-1);
    fstringSel=strrep(fstringSel,'%','');
    H2line=removeExtraSpacesLF(fstringSel);


    if  ~contains(fstring,searchTag)
        catBoolean=false;
    else
        catBoolean=true;
    end
    fclose(fid);
end


end
%FScategory:UTIGEN