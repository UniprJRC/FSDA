function [FilesWithProblems,OUT]=xmlcreateFSallFiles(InputCell,varargin)
%xmlcreateFSallFiles passes routine xmlcreateFS to all files found with makecontentsfileFS
%
% Required input arguments:
% 
%   InputCell : Cell created using routine makecontentsfileFS
%
% Optional input arguments:
%
% write2file: Option to write HTML file. Logical. Option which specifies
%             whether XML file must be created or if just structure out
%             must be created. The default value of write2file is true,
%             that is xml file is created
%             Example - 'write2file','false'
%
% Output:
%
%    FilesWithProblems : information about files whose automatic XML
%                        creation created probelms. cell of size k-by-6
%                        where k is the number of files with problems.
%                        1st column = file name
%                        2nd column = file path 
%        OUTchr         :  output of xmlcreateFS for each file. Cell of length
%                        size(InputCell,1). OUTchr{i} contains the output of
%                        the application of routine xmlcreateFS to i-th
%                        file (that is the the serialized DOM node as
%                       it appears in an XML file for the i-th file)
%
%
% Copyright 2008-2023.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit


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
    [FilesWithProblems,OUTchr]=xmlcreateFSallFiles(outTest);
%}
%

%% Beginning of code

write2file=true;

if nargin>1
    options=struct('write2file',write2file);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:xmlcreateFSallFiles:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
        
        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
        
    end
    
    write2file=options.write2file;
end 
    
FilesWithProblems=cell(1000,6);
OUT=cell(size(InputCell,1),1);
ij=1;
for i=1:size(InputCell,1)
    dirpathi=InputCell{i,end};
    disp(['Processing file: ' dirpathi filesep InputCell{i,1}])
    try
        % call xmlcreateFS
        [~,docNodechr]=xmlcreateFS(InputCell{i,1},'write2file',write2file);
        
        % Store output cell out inside OUT
        OUT{i}=docNodechr;
        
    catch
        FilesWithProblems{ij,1}=InputCell{i,1};
        % Store file path
        FilesWithProblems{ij,2}=InputCell{i,9};
        % Store catch message
        FilesWithProblems{ij,3}='notrun';
        
        ij=ij+1;
        errmsg=['Could not parse file: ' InputCell{i,1}];
        warning('FSDA:xmlcreateFSallFiles:WrongInput',errmsg)
        
    end
end
FilesWithProblems=FilesWithProblems(1:ij-1,:);
end
