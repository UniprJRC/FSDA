
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
FilesIncluded=FilesIncluded(1:10,:);
OUT=OUT(1:10);

FilesIncludedAll= FilesIncluded;
warning('on')
disp(filesWithProblems)


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
testpath= ['tests-Jenkins'];
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

cat2test='Jenkins';

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
   
