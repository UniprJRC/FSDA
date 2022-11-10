function outHELP=CreateFSDAhelpFiles()
%CreateFSDAhelpFiles prepares all documentation files

% The purpose of this function is to automatize all steps which lead to the
% HTML documentation
% The procedure is follows:
% 1) Create (using function makecontentsfileFS) a personalized contents
% file for each requested subfolder and the main folder of FSDA. Contents
% file in the main folder will contain the list of all requested files in
% the requested subfolders
%
% 2) Create for each .m file, filtered by previous procedure
% makecontentsfileFS.m  the corresponding HTML file (to accomplish
% this task use procedure publishFSallFiles).
% Check that all HTML files have been generated correctly.
%
% 3) Create files function-alpha.html and function-cate.html (list of
% functions in alphabetical and categorical order)
%
% 4) Create file function.alpha.txt which contains the list of all files
% separated by commas. This file is necessary in order to create the right
% and left buttons which enable us to navigate in alphabetical order inside the
% HTML navigation system
%
% 5) create bibliography file
%
% 6) Creare all pointer files using routine CreateFSDApointerFiles
%
% Remark: remember that your setup program must execaute command
% builddocsearchdb in folder [FSDAroot filesep 'helpfiles'  filesep 'pointersHTML']
% to create a searchable help
%
%
% Required input arguments:
%
%  Output:
%
%         outHELP:   structure which contains the following fields
%
%            outHELP.FilesIncluded  = list of files for which HTML help has
%                                       to be created
%           outHELP.FilesIncluded   = list of files for which HTML help
%                                       does not have to be created
%            outHELP.FilesWithProblems = list of files with problems in the
%                                        generation of HTML help files
%              outHELP.fileAlpha    = string containing function-alpha.html
%              outHELP.fileCate    = string containing function-cate.html
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code 

%% STEP 1 create personalized contents files
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

%% STEP 2: create HTML for all files filtered using makecontentsFilesFS

% Make sure that the format is short
format short

[FilesWithProblems,OUT,FilesIncluded]=publishFSallFiles(FilesIncluded);

% If it is necessary to call publishFSallFiles without execution of the
% examples and without creating HTML files
% [~,OUT]=publishFSallFiles(FilesIncluded, 'evalCode',false,'write2file',false);

if ~isempty(FilesWithProblems)
    % Check correctness of HTML link inside each .m file
    chkHTMLlink=cell2mat(FilesWithProblems(:,6));
    if max(chkHTMLlink)>0
        disp('Files with wrong reference to HTML page')
        disp([FilesWithProblems{chkHTMLlink==1,1}])
        error('FSDA:CreateFSDAhelpFiles','Files with wrong references in docsearchFS')
    end
    
    % Check correctness of HTML file creation
    boo=~cellfun('isempty',FilesWithProblems(:,5));
    seq=1:size(FilesWithProblems,1);
    IndexesofFiles=seq(boo);
    if ~isempty(IndexesofFiles)
        disp('Files whose HTML reference page could not be created')
        for i=1:length(IndexesofFiles)
            disp(FilesWithProblems{IndexesofFiles(i),1})
        end
        error('FSDA:CreateFSDAhelpFiles','Files without HTML web page')
    end
end

%% STEP 2B some documentation about FSDA html help pages
% Help pages are modular and contain portion of html code that are stored
% in files={'headjs.js', 'engine.js','lbar.js', 'lbarsimple.js', '...'} and
% then injected in the target html file by the document.write() JS
% function. Any other known method would fail because these methods are not
% synchronous! The issue here is that the html pages are fully dynamic and
% dense of javascript methods that must be run in the exact same
% hierarchical order and timing in which they appear in the html file.
%
% headJS.js, for example, contains instructions like 
% '<script src='../includes/shared/scripts/jquery.highlight.js'></script>";'
% which points to an '/include' folder that contains all the new and version
% bound JS functions related to the MATLAB most recent build.
% This is a convenient solution to keep the pages updated without changing
% the html page code.
%
% headJS.js also contains all the references to subsequent JS packed html
% code sections like {'engine.js','lbar.js', 'lbarsimple.js', '...'} that
% are later called in the code, so all the declarations are stored in this
% JS file and NOT inside the <body> tag of the html functions.
% 
% headJS.js is organised like a giant string variable as follows: headJS=
% "<script type='text/x-mathjax-config'>"; headJS=headJS +
% "MathJax.Hub.Config({"; 
% that unfortunatley lowers the code readability.
% Please note that due to the constraints of JS string concatenation, ALL
% the instances of double quotes should be carefully removed! 
% headJS.js is declared like: 
% "<script src='includesFS/headJS.js' type='text/javascript'></script>"; 
% and then injected in the right spot like this: 
% <script type="text/javascript">document.write(headJS);</script>
%
% static pages are different from regular dynamic pages and are created
% using a template called 'template_static_page.html' that contains the
% needed 'blueprints'.
%
% function_alpha.html and function_cate.html contain a JS file named
% 'OnThisPagefcate.js' and 'OnThisPagefalpha.js' that point to a common
% variable with the same name called by document.write(fcate); also they
% need the 'lbarsimple.js' file unlike all the other help pages.
%
% The very same pages are used in the web server (ROSA
% http://rosa.unipr.it), provided that a file named 'headJSweb.js' is used
% in place of 'engine.js'. This file replaces the MATLAB internal search
% engine with a reference to Google Search which works mostly in the same
% way and will not need any further modification of the html code.


%% STEP 3: create alphabetical list of functions and txt file
% rerun step 1 to regenerate FilesIncluded
fsep=filesep;

% Make sure one more time you are inside main root of FSDA
cd(fileparts(which('docsearchFS.m')))
% Create HTML file containing alphabetical list of functions
fileAlpha=publishFunctionAlpha(FilesIncluded,'CreateTxtFile',true);
% open html file in web browser
outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA' fsep 'function-alpha.html'];
web(outputOFHtmlHelpFile,'-browser');

fsep=filesep;
outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA' fsep 'function-alpha.txt'];
% open outfile txt in web browser
disp('Check .txt file')
web(outputOFHtmlHelpFile,'-browser');


%% STEP 4: create categorical list of functions
fsep=filesep;

% Make sure one more time you are inside main root of FSDA
cd(fileparts(which('docsearchFS.m')))
% Create HTML file containing categorical list of functions
fileCate=publishFunctionCate(FilesIncluded);
% open outfile file in web browser
outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA' fsep 'function-cate.html'];
web(outputOFHtmlHelpFile,'-browser');


%% STEP 5: create bibliography file
fsep=filesep;
% Make sure one more time you are inside main root of FSDA
cd(fileparts(which('docsearchFS.m')))

% Note that input OUT comes from function publishFSallFiles
% cell OUT can also be more easily created using option 'evalCode',false
% and 'write2file',false
% [~,OUT]=publishFSallFiles(FilesIncluded, 'evalCode',false,'write2file',false);

% Create HTML file containing all the items which make up the bibliography
[fileBiblio,Cits]=publishBibliography(FilesIncluded,OUT);

% open outfile file in web browser
outputOFHtmlHelpFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA' fsep 'bibliography.html'];
web(outputOFHtmlHelpFile,'-browser');

%% STEP 6: create HTML pointer files
h=CreateFSDApointerFiles(FilesIncluded,OUT);
if h
    disp('Successful creation of pointer files')
end

%% STEP not compulsory: create searchable database
% Use as folder the one which contains all pointers files
FileName='addFSDA2path';
FullPath=which(FileName);
%Navigate to the main folder of FSDA
FSDAroot=fileparts(FullPath);
% Navigate to subfolder which contains pointerHTML
pointersHTMLroot=[FSDAroot filesep 'helpfiles'  filesep 'pointersHTML'];
% Create searchable database
builddocsearchdb(pointersHTMLroot)


%% Step 7 Create zip file containing images of HTML files
% zip all images files from (FSDAroot)/helpfiles/FSDA/images/
% into (FSDAroot)/helpfiles/FSDA/FSDAweb/images.zip
% This file will have to be sent to rosa web site
zipfilename=[FSDAroot fsep 'helpfiles' filesep 'FSDAweb' filesep 'images.zip'];
filenames=[FSDAroot fsep 'helpfiles' filesep 'FSDA' filesep 'images' filesep];
zip(zipfilename,filenames);


%% DOCUMENTATION FOR THE WEB SITE http://rosa.unipr.it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE STEPS WHICH FOLLOW ARE NECESSARY TO CREATE THE DOCUMENTATION ON THE WEB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 1web: create HTML to publish in rosa.unipr.it for all files filtered using makecontentsFilesFS

[FilesWithProblemsweb,OUTweb]=publishFSallFiles(FilesIncluded,'webhelp',true);

% Check correctness of HTML file creation
boo=~cellfun('isempty',FilesWithProblemsweb(:,5));
seq=1:size(FilesWithProblemsweb,1);
IndexesofFiles=seq(boo);
if ~isempty(IndexesofFiles)
    disp('Files whose HTML reference page could not be created')
    for i=1:length(IndexesofFiles)
        disp(FilesWithProblemsweb{IndexesofFiles(i),1})
    end
    error('FSDA:CreateFSDAhelpFiles','Files without HTML web page for WEB')
end

%% STEP 2web: Insert in all files not parsed by publishFSallFiles Google search engine
%  Patch with Google Search module all static files.
%  Then, copy external static files: acknowledgments.html, developers.html, group.html
%  license.txt, links_relevant.html, poster_fsda.pdf to FSDAweb folder.
%  In the end create a sitemap for Google Search module.
%
FileName='addFSDA2path';
FullPath=which(FileName);
%Navigate to the main folder of FSDA
FSDAroot=fileparts(FullPath);

fsep=filesep;
% List of static files
ListofFiles={'bibliography.html' 'cluster_intro.html' 'datasets.html' ...
    'datasets_clu.html' 'datasets_mv.html' 'datasets_reg.html'...
    'examples.html',  'function-alpha.html', 'function-cate.html', ...
    'getting-started.html' 'index.html' 'introFS.html' 'introrob.html' ...
    'introrobmulttech.html' 'introrobregtech.html' 'mult_fsm.html' ...
    'mult_fsmeda.html' 'mult_fsmfan.html' 'mult_fsmtra.html' 'mult_mcd.html'...
    'mult_sandmm.html' 'mult_unibiv.html' 'multivariate_intro.html'...
    'multivariatetransf_intro.html' ...
    'regression_fsr.html' 'regression_fsreda.html' 'regression_intro.html'...
    'regression_lxs.html' 'regression_mms.html' 'regression_mscp.html'...
    'regression_mst.html' 'regressionms_intro.html' 'release_notes.html'...
    'statistical_visualizationFS.html' 'statistical_visualization_cds.html'...
    'statistical_visualization_fan.html' 'statistical_visualization_intro.html'...
    'statistical_visualization_mdr.html' 'statistical_visualization_monres.html'...
    'statistical_visualization_resindex.html' 'statistical_visualization_yx.html'...
    'transf_fsrfan.html' 'transf_intro.html' 'transf_score.html' 'tutorials.html'};

FilesNotFound=insertGoogleSearchEngine(ListofFiles);

if ~isempty(FilesNotFound)
    disp('Files for which google search engine could not be inserted')
    for i=1:length(FilesNotFound)
        disp(FilesNotFound{i,1})
    end
    error('FSDA:CreateFSDAhelpFiles','Files not found')
end

extFiles= {'acknowledgments.html' 'developers.html' 'group.html' 'license.txt' ...
    'links_relevant.html' 'poster_fsda.pdf'};

for i=1:length(extFiles)
    inputFile=[FSDAroot fsep 'helpfiles' fsep 'FSDA' fsep extFiles{i}];
    status=copyfile(inputFile, ...
        [FSDAroot fsep 'helpfiles' fsep 'FSDAweb' fsep extFiles{i}]);
    if status==0
        disp(['File: ' inputFile ' not found'])
        error('FSDA:CreateFSDAhelpFiles','File: not found')
    end
    
end
% copy various files ....
inputFile=[FSDAroot fsep 'InstallationNotes.pdf'];
status=copyfile(inputFile, ...
    [FSDAroot fsep 'helpfiles' fsep 'FSDAweb' fsep 'InstallationNotes.pdf']);
if status==0
    disp(['File: ' inputFile ' not found'])
    error('FSDA:CreateFSDAhelpFiles','File: not found')
end


% Create a simple URL Sitemap in txt format
outputOFHtmlHelpFileWeb=[FSDAroot fsep 'helpfiles' fsep 'FSDAweb'];
allWebFiles=getWebMatlabFiles(outputOFHtmlHelpFileWeb);
allHttpUrl=cell(length(allWebFiles),1);

for i=1:length(allWebFiles)
    allHttpUrl{i,1}=['http://rosa.unipr.it/FSDA/' allWebFiles{i}];
end

fileSitemap=cell2table(allHttpUrl);
writetable(fileSitemap,[FSDAroot fsep 'helpfiles' fsep 'FSDAweb' fsep 'sitemap.txt'],...
    'WriteVariableNames',false);

%% Step 3 Create zip file containing images of HTML files
% zip all images files from (FSDAroot)/helpfiles/FSDA/images/
% Note that last / otherwise the files inside the zip archive will have
% subfolder images
% into (FSDAroot)/helpfiles/FSDA/FSDAweb/images.zip
zipfilename=[FSDAroot fsep 'helpfiles' filesep 'FSDAweb' filesep 'images.zip'];
filenames=[FSDAroot fsep 'helpfiles' filesep 'FSDA' filesep 'images' filesep];
zip(zipfilename,filenames);

%% Now if all was well let us do the setup.exe
disp('Congratulations the FSDA package is ready to be deployed')

% Store all quantities inside structure outHELP
outHELP=struct;
outHELP.FilesIncluded=FilesIncluded;
outHELP.FilesExcluded=FilesExcluded;
outHELP.FilesWithProblems=FilesWithProblems;
outHELP.fileAlpha=fileAlpha;
outHELP.fileCate=fileCate;
outHELP.fileBiblio=fileBiblio;
outHELP.OUT=OUT;
outHELP.OUTweb=OUTweb;
outHELP.Cits=Cits;

end
