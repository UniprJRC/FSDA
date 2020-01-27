
%% Load necessary elements for performance test
% load OUT
run addFSDA2path
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
[~,OUT]=publishFSallFiles(FilesIncluded, 'evalCode','false',...
    'write2file',false,'ErrWrngSeeAlso',false);

%% Category to test
cat2test=getenv('CATEGORY_TO_TEST');
disp('---------------')
disp('Test for category:')
disp(cat2test)
disp('---------------')

% VIS GUI     MULT CLUS REG UTI
if strcmp(cat2test,'graphics')
    str=regexp(FilesIncluded(:,8),'VIS*');
    boo1=~cellfun(@isempty,str);
    str=regexp(FilesIncluded(:,8),'GUI');
    boo2=~cellfun(@isempty,str);
    boo=boo1 | boo2;
    
elseif strcmp(cat2test,'multivariate')
    str=regexp(FilesIncluded(:,8),'MULT*');
    strnot1=regexp(FilesIncluded(:,8),'CLUS-RobClaMULT');
    boo=~cellfun(@isempty,str) & cellfun(@isempty,strnot1);
    
elseif strcmp(cat2test,'multivariate-clustering')
    str=regexp(FilesIncluded(:,8),'CLUS-RobClaMULT');
    boo=~cellfun(@isempty,str);
    
elseif strcmp(cat2test,'regression-clustering')
    str=regexp(FilesIncluded(:,8),'CLUS-RobClaREG');
    boo=~cellfun(@isempty,str);
    
elseif strcmp(cat2test,'mixsim')
    str=regexp(FilesIncluded(:,8),'CLUS-MixSim');
    boo=~cellfun(@isempty,str);
    
elseif strcmp(cat2test,'regression')
    str=regexp(FilesIncluded(:,8),'REG-Regression');
    boo=~cellfun(@isempty,str);
    OUT=OUT(boo,:);
    FilesIncluded=FilesIncluded(boo,:);
    
    str=regexp(FilesIncluded(:,1),'(TS)|(ts)*');
    boo=cellfun(@isempty,str);
    
elseif strcmp(cat2test,'regressionTS')
    str=regexp(FilesIncluded(:,8),'REG-Regression');
    boo=~cellfun(@isempty,str);
    OUT=OUT(boo,:);
    FilesIncluded=FilesIncluded(boo,:);
    str=regexp(FilesIncluded(:,1),'(TS)|(ts)*');
    booTS=~cellfun(@isempty,str);
    % Remove LTSts
    str=regexp(FilesIncluded(:,1),'LTSts.m');
    booNotLTSts=cellfun(@isempty,str);
    boo=booTS & booNotLTSts;
    
elseif strcmp(cat2test,'regressionLTSts')
    str=regexp(FilesIncluded(:,1),'LTSts.m');
    boo=~cellfun(@isempty,str);
    
elseif strcmp(cat2test,'regressionEXT')
    str=regexp(FilesIncluded(:,8),'REG-*');
    strnot1=regexp(FilesIncluded(:,8),'REG-Regression');
    strnot2=regexp(FilesIncluded(:,8),'CLUS-RobClaREG');
    boo=~cellfun(@isempty,str) & cellfun(@isempty,strnot1) & cellfun(@isempty,strnot2);
    
elseif strcmp(cat2test,'utilities')
    str=regexp(FilesIncluded(:,8),'UTI*');
    boo=~cellfun(@isempty,str);
    
else
    error('FSDA:runTests:WrgCLS','Wrong class')
end

OUT=OUT(boo,:);
FilesIncluded=FilesIncluded(boo,:);


ij=1;
nfiles=length(OUT);
sz=[5000, 7];
TotSummary = table('Size',sz,'VariableTypes',{'cellstr' 'cellstr' 'cellstr' 'double' 'double' 'cellstr' 'cellstr'},...
    'VariableNames',{'FileName' 'Category', 'Identifier' 'MeanTime' 'MedianTime'  'Code' 'TestActivity'});

% [tmp]=publishFS('existFS', 'evalCode','false',...
%     'write2file',false,'ErrWrngSeeAlso',false);
% disp(tmp)


%% Performance and Testing part
% make sure to be in the FSDAroot
cd(FSDAroot)

% Use perf = true if for each example you want run runperf.m
% Use perf = true if for each example you want run runtests.m
perf = false;
testpath= ['tests-' cat2test];
mkdir(testpath);


for i=1:nfiles
    clc
    disp(['Filename ' FilesIncluded{i,1}])
    disp(['Executing file ' FilesIncluded{i,1} '  Number  ' num2str(i) ' of ' num2str(nfiles)])
    
    Ex=OUT{i,1}.Ex;
    
    Extra=OUT{i,1}.ExtraEx;
    Excomb=[Ex;Extra];
    for iEx=1:size(Excomb,1)
        disp(['Running example number ' num2str(iEx)  ' of '   num2str(size(Excomb,1)) '  File: '  FilesIncluded{i,1}]);
        
        close all
        Exi=Excomb{iEx,3};
        if ~isempty(Exi) && isempty(strfind(Excomb{iEx,1},'Interactive example'))
            
            Exi=regexprep(Exi,'&lt;','<');
            Exi=regexprep(Exi,'&gt;','>');
            close all
            if iEx==1
                Exif=[Exi,newline,'close all',newline 'save tempfileWS'];
            else
                Exif=['load tempfileWS',newline,Exi,newline,'close all',newline, 'save tempfileWS'];
            end
            
            % Write Exif to a file which name begins with 'text'
            filename2open=[testpath '/test' FilesIncluded{i,1}(1:end-2) '_' num2str(iEx) ...
                'of'  num2str(size(Excomb,1)) FilesIncluded{i,1}(end-1:end)];
            file1ID=fopen(filename2open,'w');
            %file1ID=fopen('tempfile.m','w');
            fprintf(file1ID,'%s',Exif);
            fclose('all');
            try
                if perf==1
                    outp=runperf('tempfile.m');
                    MeanS=outp.sampleSummary.Mean;
                    FindNaN=isnan(MeanS);
                    MeanS=MeanS(~FindNaN);
                    TotSummary{ij,'MeanTime'}= MeanS;
                    
                    MedianS=outp.sampleSummary.Median;
                    MedianS=MedianS(~FindNaN);
                    TotSummary{ij,'MedianTime'}= MedianS;
                    
                    % TotSummary{ij,'MedianTime'}= outp.sampleSummary.Median;
                    TotSummary(ij,'TestActivity')={outp.TestActivity};
                else
                    %tic
                    %outp=runtests(filename2open);
                    %time=toc;
                    %TotSummary{ij,'MeanTime'}=time;
                end
%                 TotSummary(ij,'Code')=Ex(1,3);
%                 TotSummary(ij,'FileName')=FilesIncluded(i,1);
%                 TotSummary(ij,'Category')=FilesIncluded(i,8);
%                 TotSummary(ij,'Identifier')={['Ex' num2str(iEx)]};
                ij=ij+1;
            catch
                disp(['Error on example ' num2str(iEx)])
                disp(['Name of the file: '  FilesIncluded{i,1}])
                warning('Stop here')
            end
        else
            disp('Interactive example')
        end
    end
    % disp(['Error on the following function: (fx .num:)' int2str(i)])
    % OUT{i,1}.titl;
end

cd(testpath);
testResults = runtests;
cd ..
save ([cat2test '_testResults.mat'], 'testResults')
%rmdir('tests', 's')

% TotSummary1=TotSummary(1:ij-1,:);
% disp(TotSummary1)
% cfol=pwd
% FSDAroot 
% filename = [FSDAroot '/test-results/' cat2test '_test.xlsx'];
% writetable(TotSummary1,filename,'Sheet',1,'Range','A1');


