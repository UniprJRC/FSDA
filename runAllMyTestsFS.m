
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
warning('off')
[FilesIncluded,FilesExcluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',force,'printOutputCell','Contents.m');
[filesWithProblems,OUT]=publishFSallFiles(FilesIncluded, 'evalCode','false',...
    'write2file',false,'ErrWrngSeeAlso',false,'msg',false);
FilesIncludedAll= FilesIncluded;
warning('on')
disp(filesWithProblems)
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
    % multivariate clustering excluding tclust*
    str=regexp(FilesIncluded(:,8),'CLUS-RobClaMULT');
    booMULT=~cellfun(@isempty,str);
    
    str=regexp(FilesIncluded(:,1),'tclust*');
    booTCLUST=~cellfun(@isempty,str);
    boo=booMULT & ~booTCLUST;
    
elseif strcmp(cat2test,'tclustMULT')
    % multivariate clustering just tclust*
    str=regexp(FilesIncluded(:,8),'CLUS-RobClaMULT');
    booMULT=~cellfun(@isempty,str);
    str=regexp(FilesIncluded(:,1),'tclust*');
    booTCLUST=~cellfun(@isempty,str);

    str=regexp(FilesIncluded(:,1),'(gpcm)|(GPCM)*');
    notGPCM=cellfun(@isempty,str);

    boo=booMULT & booTCLUST & notGPCM;

 elseif strcmp(cat2test,'tclustMULTgpcm')
    % multivariate clustering just tclust*
    str=regexp(FilesIncluded(:,8),'CLUS-RobClaMULT');
    booMULT=~cellfun(@isempty,str);
    str=regexp(FilesIncluded(:,1),'tclust*');
    booTCLUST=~cellfun(@isempty,str);

    str=regexp(FilesIncluded(:,1),'(gpcm)|(GPCM)*');
    GPCM=~cellfun(@isempty,str);
    boo=booMULT & booTCLUST & GPCM;   
    
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
    
    % disp(['Filename ' FilesIncluded{i,1}])
    disp(['File ' FilesIncluded{i,1} '  Number  ' num2str(i) ' of ' num2str(nfiles)])
    
    Ex=OUT{i,1}.Ex;
    
    Extra=OUT{i,1}.ExtraEx;
    Excomb=[Ex;Extra];
    for iEx=1:size(Excomb,1)
        disp(['Writing to file  example number ' num2str(iEx)  ' of '   num2str(size(Excomb,1)) '  contained in file: '  FilesIncluded{i,1}]);
        
        close all
        Exi=Excomb{iEx,3};
        if ~isempty(Exi) && isempty(strfind(Excomb{iEx,1},'Interactive example'))
            
            Exi=regexprep(Exi,'&lt;','<');
            Exi=regexprep(Exi,'&gt;','>');
            close all
            if iEx==1
                %Exif=[Exi,newline,'close all',newline 'save tempfileWS'];
                Exif=[Exi,newline,'close all'];
            else
                %Exif=['load tempfileWS',newline,Exi,newline,'close all',newline, 'save tempfileWS'];
                Exif=[Exi,newline,'close all'];
            end
            
            % Write Exif to a file which name begins with 'text'
            filename2open=[testpath '/test' FilesIncluded{i,1}(1:end-2) '_' num2str(iEx) ...
                'of'  num2str(size(Excomb,1)) FilesIncluded{i,1}(end-1:end)];
            file1ID=fopen(filename2open,'w');
            %file1ID=fopen('tempfile.m','w');
            fprintf(file1ID,'%s',Exif);
            fclose('all');
        else
            disp('Interactive example')
        end
    end
    % disp(['Error on the following function: (fx .num:)' int2str(i)])
    % OUT{i,1}.titl;
end


cd(testpath);
import matlab.unittest.plugins.CodeCoveragePlugin;


% Clear last warning
lastwarn('');
warn1 = lastwarn;

suite = testsuite(pwd);

warn2 = lastwarn;
if ~strcmp(warn1, warn2)
    disp(warn2)
    error('FSDA:runAllMyTestsFS:WrongExample','A warning occurred during test suite creation!')
end

if perf==false
    
    % Create a TestRunner
    runner = matlab.unittest.TestRunner.withTextOutput();
    
    % Add a plugin to produce a JUnit-style test report
    runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat(['test-' cat2test '-report.xml']));
    
    % Get file paths of source code being tested
    filePaths = fullfile(FilesIncluded(:,9), FilesIncluded(:,1));
    % Indicate where the Cobertura coverage report should be created
    covFile = matlab.unittest.plugins.codecoverage.CoberturaFormat(['coverage-' cat2test '-report.xml']);
    % Add the CodeCoveragePlugin
    runner.addPlugin(matlab.unittest.plugins.CodeCoveragePlugin.forFile(filePaths, 'Producing', covFile));
    
    % Run the test suite
    disp('Run the test suite')
    runner.run(suite);
    
    % create a new subfolder for CircleCI
    [status, ~, ~] = mkdir(['test-' cat2test '-report']);
    copyfile(['test-' cat2test '-report.xml'], ['test-' cat2test '-report']);
    
    
else
    import matlab.perftest.TimeExperiment
    experiment = TimeExperiment.limitingSamplingError('NumWarmups',1,...
        'MaxSamples',4,'RelativeMarginOfError',0.08,'ConfidenceLevel',0.97);
    resultsTE = run(experiment,suite);
    
end

%% CircleCI 

%cd(testpath);
%testResults = runtests([FSDAroot '/' testpath]);
%cd ..
%save ([cat2test '_testResults.mat'], 'testResults')
%rmdir('tests', 's')

% TotSummary1=TotSummary(1:ij-1,:);
% disp(TotSummary1)
% cfol=pwd
% FSDAroot
% filename = [FSDAroot '/test-results/' cat2test '_test.xlsx'];
% writetable(TotSummary1,filename,'Sheet',1,'Range','A1');


