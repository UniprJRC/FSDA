function [FilesWithProblems,OUT,InputCell]=publishFSallFiles(InputCell,varargin)
%publishFSallFiles passes routine publishFS to all files found with makecontentsfileFS
%
% Required input arguments:
% 
%   InputCell : Cell created using routine makecontentsfileFS
%
% Optional input arguments:
%   
% evalCode :  Option to run code. Logical. Option to evaluate code of the
%             examples in the input .m files enclosed in tags "%{" "%}" whose
%             first line starts with symbols "%%".
%             If evalcode is true the code associated with the examples
%             which start with symbols '%%' will be run and the output will
%             be automatically included into the HTML output file. The
%             images will be save in subfolder images_help of the
%             outputDir. The default value of evalCode is true.
%             Example - 'evalCode','false'
%             Data Types - Boolean
% write2file: Option to write HTML file. Logical. Option which specifies
%             whether HTML file must be created or if just structure out
%             must be created. The default value of write2file is true,
%             that is html file is created
%             Example - 'write2file','false'
%             Data Types - Boolean
% outputDir : Output folder. String.
%             Output folder to which the HTML document is saved, specified
%             as the comma-separated pair consisting of 'outputDir' and the
%             full path. You must specify the full path as a string, for
%             example 'C:\PublishedOutput'.
%             if input option webhelp is false the default value, '',
%             specifies the (FSDA root)\helpfiles\FSDA path,
%             else if the input option webhelp is true the default value,
%             '', specifies the (FSDA root)\helpfiles\FSDAweb path,
%             Remark - outputDir must be a valid path.
%             Example - 'outputDir','C:'
%             Data Types - string
% imagesDir : Output folder of png images. String.
%             Output folder to which the images attached to the HTML
%             document are saved, specified as the comma-separated pair
%             consisting of 'outputDir' and the full path. You must specify
%             the full path as a string, for example
%             'C:\PublishedOutput'.
%             if input option webhelp is false the default value, '',
%             specifies the "(FSDA root)\helpfiles\FSDA\images" path,
%             else if the input option webhelp is true the default value,
%             '', specifies the "(FSDA root)\helpfiles\FSDAweb\images"
%             path.
%             Remark - if imageDir is not specified but outputDir is
%             specified images will be saved into the same folder of the
%             HTML output file
%             Remark - imagesDir must be a valid path.
%             Example - 'imagesDir','C:'
%             Data Types - string
% webhelp :   Option to create web page version of local html help file. Logical.
%             This option substitutes the MATLAB search form component with 
%             a Google local search form, enabling Google to create an index
%             of online documentation help.
%             In this case 'outputDir' and 'imagesDir' should be different
%             from  those previously used wiht option webhelp=false to 
%             avoid overwriting local help files with WEB help files.
%             Default values for 'outputDir' and 'imagesDir' are these namepairs: 
%             'outputDir','(FSDA root) filesep helpfiles filesep FSDAweb'
%             'imagesDir','(FSDA root) filesep helpfiles filesep FSDAweb filesep images'
%             It is worth noting that setting up folders different form default
%             values, implies that the user should manually create the path.
%             The default value of webhelp is false.
%             Example - 'webhelp',true
%             Data Types - logical
%ErrWrngSeeAlso: Option to check links in the see also part. Logical.
%            If ErrWrngSeeAlso is true publishFS checks whether the strings
%            inside see also are valid files and puts an hyperlink to the
%            corresponding HTML file. If publishFS cannot find the files
%            exits the procedure with an error. If ErrWrngSeeAlso is false
%            it just produces a warning. Default value of ErrWrngSeeAlso is
%            true.
%             Example - 'ErrWrngSeeAlso',false
%             Data Types - logical
%
% Output:
%
%    FilesWithProblems : information about files whose automatic HTML
%                        creation created probelms. cell of size k-by-6
%                        where k is the number of files with problems.
%                        1st column = file name
%                        2nd column = file path 
%                        3rd column = information about Input Arguments mismatch
%                        4th column = information about Output Arguments mismatch
%                        5th column = information about lack of execution
%                        of examples inside m files
%                        6th column = information about mismatch in
%                        docsearchFS. This column is true if docsearchFS
%                        does not contain correct information
%        OUT          :  output of publishFS for each file. Cell of length
%                        size(InputCell,1). OUT{i} contains the output of
%                        the application of routine publishFS to i-th file
%
%       InputCell     :  What is given is InputCell is also shown in the
%                        output. This is done because workspace is cleared
%                        and InputCell is lost if it is not an output
%                        argument
%
% Copyright 2008-2019.
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
    [FilesWithProblems,OUT]=publishFSallFiles(outTest,'evalCode',false,'write2file',false);
%}
%

%% Beginning of code

% % Use file separator of current operating system
% % \ = Windows
% % / = Unix
fsep=filesep;

% Write output file in subfolder \(FSDAroot)\helpfiles\FSDA
FileWithFullPath=which('docsearchFS.m');
[pathFSDAstr]=fileparts(FileWithFullPath);

evalCode=true;
write2file=true;
Display='none';
ErrWrngSeeAlso=true;

if nargin>1
    UserOptions=varargin(1:2:length(varargin));
    checklms2 = strcmp(UserOptions,'webhelp');
    if sum(checklms2)
        webhelp = varargin{2*find(checklms2)};
    else
        webhelp=false;
    end
else
    webhelp=false;
end

if webhelp == false
    outputDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA'];
    imagesDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA' fsep 'images'];
else
    outputDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDAweb'];
    imagesDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDAweb' fsep 'images'];
end

if nargin>1
    options=struct('evalCode',evalCode,'write2file',write2file,'Display',Display,...
       'webhelp',webhelp,'outputDir',outputDir,'imagesDir',imagesDir,...
       'ErrWrngSeeAlso',ErrWrngSeeAlso);
    
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
    % set the options chosen by the user
    evalCode=options.evalCode;
    write2file=options.write2file;
    Display=options.Display;
    webhelp=options.webhelp;
    outputDir=options.outputDir;
    imagesDir=options.imagesDir;
end 
    
FilesWithProblems=cell(1000,6);
OUT=cell(size(InputCell,1),1);
ij=1;
for i=1:size(InputCell,1)
    dirpathi=InputCell{i,end};
    disp(['Processing file: ' dirpathi filesep InputCell{i,1}])
    try
        % call publishFS
        %out=publishFS(InputCell{i,1},'evalCode',evalCode,'write2file',write2file);
        
        out=publishFS(InputCell{i,1},'evalCode',evalCode,'Display',Display,...
        'webhelp',webhelp,'outputDir',outputDir,...
        'imagesDir',imagesDir,'write2file',write2file,'ErrWrngSeeAlso',ErrWrngSeeAlso);
        % Store output cell out inside OUT
        OUT{i}=out;
        
        if  (size(out.InpArgsMisMatch,1)+size(out.OutArgsStructMisMatch,1))>2 || ~isempty(out.laste) || out.linkHTMLMisMatch==1
            % Store information about files with problems
            % Store file name
            FilesWithProblems{ij,1}=InputCell{i,1};
            % Store file path
            FilesWithProblems{ij,2}=InputCell{i,9};
            % Store InputMismatch
            FilesWithProblems{ij,3}=out.InpArgsMisMatch;
            % Store OutputMismatch
            FilesWithProblems{ij,4}=out.OutArgsStructMisMatch;
            % Store laste
            FilesWithProblems{ij,5}=out.laste;
            FilesWithProblems{ij,6}= out.linkHTMLMisMatch;
            
            ij=ij+1;
        end
        
    catch
        FilesWithProblems{ij,1}=InputCell{i,1};
        % Store file path
        FilesWithProblems{ij,2}=InputCell{i,9};
        % Store catch message
        FilesWithProblems{ij,3}='Not Parsed';
        
        ij=ij+1;
        errmsg=['Could not parse file: ' InputCell{i,1}];
        warning('FSDA:publishFSallFiles:WrongInput',errmsg)
        
    end
end
FilesWithProblems=FilesWithProblems(1:ij-1,:);
end
