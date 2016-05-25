function FilesWithProblems=publishFSallFiles(InputCell)
%publishFSallFiles passes routine publishFS to all files found with makecontentsfileFS
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
    [outTest,Excluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',false)
    FilesWithProblems=publishFSallFiles(outTest);
%}
%

%% Beginning of code
FilesWithProblems=cell(1000,6);
ij=1;
for i=1:size(InputCell,1)
    dirpathi=InputCell{i,end};
    disp(['Processing file: ' dirpathi filesep InputCell{i,1}])
    try
        % call publishFS
        out=publishFS(InputCell{i,1});
        
        
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
        FilesWithProblems{ij,3}='notrun';
        
        ij=ij+1;
        errmsg=['%s','Could not parse file: ' InputCell{i,1}];
        warning('FSDA:publishFSallFiles:WrongInput',errmsg)
        
    end
end
FilesWithProblems=FilesWithProblems(1:ij-1,:);
end