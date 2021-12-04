% This routine exports all dataset file of FSDA into Excel format
% (FSDAroot)/datasetExcel is created with subfolders
% clustering 
% multivariate
% multivariate_regression
% regression
% Every subfoler contains all FSDA datasets in Excel (.xlsx) format

%% Beginning of code

%% Find FSDA root
FileName='addFSDA2path';
FullPath=which(FileName);
%Navigate to the main folder of FSDA
FSDAroot=fileparts(FullPath);
cd(FSDAroot)

% create subfolder datasetsExcel
mkdir('datasetsExcel')
cd('datasetsExcel\')
datasetExcelFolder=pwd;
cd('..\')

%% Export datasets inside subfolder clustering

cd(datasetExcelFolder)
folders={'clustering' 'multivariate' 'regression' 'multivariate_regression'};
for ii=1:length(folders)
    folderi=folders{ii};

    mkdir(folderi)
    pathFSDAfolder=[FSDAroot '\datasets\' folderi '\'];
    pathFSDAfolderExcel=[datasetExcelFolder '\' folderi '\'];
    d=dir([pathFSDAfolder '*.mat']);
    n=length(d);
    if n==0
        error('FSDA:exportFSDAdatasetstoExcel','Files not found')
    end
    for i=1:n
        FileName=d(i).name;
        [~,FileNameNoext]=fileparts(FileName);
        load(FileName)
        writetable(eval(['' FileNameNoext '']),[pathFSDAfolderExcel FileNameNoext '.xlsx'],'WriteRowNames',true)
    end

end