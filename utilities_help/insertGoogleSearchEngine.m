function [FilesWithProblems,OUT]=insertGoogleSearchEngine(InputCell,varargin)
%insertGoogleSearchEngine insert instructions to include google search engine
%
% Required input arguments:
%
%   InputCell : n-by-1 cell containing list of files
%
% Optional input arguments:
%
% outputDir : Output folder. String.
%             Output folder to which the HTML document is saved, specified
%             as the comma-separated pair consisting of 'outputDir' and the
%             full path. You must specify the full path as a string, for
%             example 'C:\PublishedOutput'.
%             The default value, '', specifies the (FSDA root)\helpfiles\FSDA
%             path.
%             Remark - outputDir must be a valid path.
%             Example - 'outputDir','C:'
%             Data Types - string
%
%
%
% Output:
%
%    FilesWithProblems :
%
%

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
    [outTest,Excluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',false);
    [FilesWithProblems,OUT]=publishFSallFiles(outTest,'evalCode',false,'write2file',false);
%}
%

%% Beginning of code
% Write output file in subfolder \(FSDAroot)\helpfiles\FSDA
FileWithFullPath=which('docsearchFS.m');
[pathFSDAstr]=fileparts(FileWithFullPath);

fsep=filesep;

outputDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDAweb'];


if nargin>1
    options=struct('outputDir',outputDir,'write2file',true);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:regressB:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
        
        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
        
    end
    
    write2file=options.write2file;
    outputDir=options.outputDir;
    
    
end



searchTagIni='includesFS/headJS.js';
searchTagFin='/headJS.js';

InputCell=InputCell';

FilesWithProblems=[];


for i=1:size(InputCell,1)
    dirpathi=[pathFSDAstr '\helpfiles\FSDA\' InputCell{i,end}];
    disp(['Processing file: ' dirpathi ])
    out= get_SearchEngine(InputCell{i,end},dirpathi,searchTagIni,searchTagFin,outputDir);
    if (out==-1)
       FilesWithProblems=[FilesWithProblems  InputCell{i,end} ', ']; 
    end
end
OUT=out;

end

function [file1ID] = get_SearchEngine(filename,filenamefullpath,letterINI,letterFIN,outputDir)
%insert search enegines lines
%    letterINI='DOVE INIZIA A RIMPIAZZARE';
%
%     letterINI='DOVE FIISCE A RIMPIAZZARE';



% [~,name,ext] = fileparts(filename);

fid = fopen(filenamefullpath); % open file

if fid > 1
    
    % now get the part which contains seacrh engine
    fstring=fscanf(fid,'%c');
    
    
    
    strInsert='web';
    %                strInsert= sprintf([strInsert '<tr>\r' ...
    %                 '<td class="term"><a href="' mfilename '.html">' mfilename '</a></td>\r' ...
    %                 '<td class="description">' description '</td>\r' ...
    %                 '</tr>\r']);
    
    idxStart=strfind(fstring, letterINI);
    offsetINI=length(letterINI);
    offsetFIN=length(letterFIN);
    idxEnd=strfind(fstring, letterFIN);
    
    fstring=[fstring(1:idxStart+offsetINI-4) strInsert fstring(idxStart+offsetINI-3:end)];
    
    fclose(fid);
    
    file1ID=fopen([outputDir filesep filename],'w');
    
    if file1ID==-1
        
        if ismac || isunix
            errmsg= [' Path ' outputDir '/' name '.html does not exist or output file '  name '.html is not writable'];
        elseif ispc
            outputDir=strrep(outputDir,'\','\\');
            errmsg= [' Path ' outputDir '\\' name '.html does not exist or output file '  name '.html is not writable'];
        else
            errmsg= [' Path ' outputDir '/' name '.html does not exist or output file '  name '.html is not writable'];
        end
        
        error('FSDA:publishFS:WrngOutFolder',errmsg);
    end
    
    fprintf(file1ID,'%s',fstring);
    fclose('all');
else
    fclose('all');
    file1ID=-1;
end
end