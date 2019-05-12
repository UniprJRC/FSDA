%% Load necessary elements for performance test
% load OUT

FileName='addFSDA2path';
FullPath=which(FileName);
%Navigate to the main folder of FSDA
FSDAroot=fileparts(FullPath);
cd(FSDAroot)
% Specify subfolders of main folders for which contents file has to be
% created
InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat' 'utilities_help'};
ExclDir={'privateFS'  'datasets'};
% Create list of folders which must have the personalized contents file
list = findDir(FSDAroot,'InclDir',InclDir,'ExclDir',ExclDir);
% Crete personalized contents file for main folder of FSDA
% and required subfolders.
force=false;
[FilesIncluded,FilesExcluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',force,'printOutputCell','Contents.m');
disp('List of files which have been excluded (with path)')
disp(FilesExcluded(:,[1 9]))
[~,OUT]=publishFSallFiles(FilesIncluded, 'evalCode','false','write2file',false);

ij=1;
nfiles=length(OUT);
sz=[5000, 7];
TotSummary = table('Size',sz,'VariableTypes',{'cellstr' 'cellstr' 'cellstr' 'double' 'double' 'cellstr' 'cellstr'},...
    'VariableNames',{'FileName' 'Category', 'Identifier' 'MeanTime' 'MedianTime'  'Code' 'TestActivity'});

%% Performance part
% nfiles = number of files

for i=1:nfiles
    try
    Ex=OUT{i,1}.Ex;
    Extra=OUT{i,1}.ExtraEx;
    Excomb=[Ex;Extra];
    for iEx=1:size(Excomb,1)
        close all
        Exi=Excomb{iEx,3};
        
        if ~isempty(Exi) && isempty(strfind(Excomb{iEx,1},'Interactive example'))
            
            Exi=regexprep(Exi,'&lt;','<');
            Exi=regexprep(Exi,'&gt;','>');
            
            % Write Exi to a file
            file1ID=fopen('tempfile.m','w');
            fprintf(file1ID,'%s',Exi);
            fclose('all');
            outp=runperf('tempfile.m');
            
            MeanS=outp.sampleSummary.Mean;
            FindNaN=isnan(MeanS);
            MeanS=MeanS(~FindNaN);
            TotSummary{ij,'MeanTime'}= MeanS;

            MedianS=outp.sampleSummary.Median;
            MedianS=MedianS(~FindNaN);
            TotSummary{ij,'MedianTime'}= MedianS;
            
            % TotSummary{ij,'MedianTime'}= outp.sampleSummary.Median;
            TotSummary(ij,'Code')=Ex(1,3);
            TotSummary(ij,'TestActivity')={outp.TestActivity};
            TotSummary(ij,'FileName')=FilesIncluded(i,1);
            TotSummary(ij,'Category')=FilesIncluded(i,8);
            TotSummary(ij,'Identifier')={['Ex' num2str(iEx)]};
            ij=ij+1;
        end
    end
    catch
        disp(['Error on the following function: (fx .num:)' int2str(i)])
        OUT{i,1}.titl;   
    end
end
