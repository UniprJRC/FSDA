
ListofFiles={'bibliography.html' 'cluster_intro.html' 'datasets.html' ...
    'datasets_clu.html' 'datasets_mv.html' 'datasets_reg.html' 'empty.html'...
    'examples.html' ...
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

%ListofFiles=ListofFiles';


[issues,OUT]=insertGoogleSearchEngine(ListofFiles);

% rthin subsets
%'function-cateEmpty.html'
ListofFiles2=getWebMatlabFiles();
%ListofFiles2=ListofFiles2';
ListofFiles2=sort(ListofFiles2);
%ListofFiles2{1,274}='index.html'
[htmlDinamicFiles,ia]=setdiff(ListofFiles2, ListofFiles);
[htmlStaticFiles,ia]=setdiff(ListofFiles2, htmlDinamicFiles);

htmlStaticFiles=htmlStaticFiles';

[~,nfiles]=size(htmlStaticFiles);
%%
FileName='addFSDA2path';
FullPath=which(FileName);
root=FullPath(1:end-length(FileName)-3);
InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat' 'utilities_help'};
ExclDir={'privateFS'  'datasets'};
% Create list of folders which must have a personalized contents file
list = findDir(root,'InclDir',InclDir,'ExclDir',ExclDir)
% Crete personalized contents file for main folder of FSDA
% and required subfolders.
[outTest,Excluded]=makecontentsfileFS('dirpath',list,'FilterFileContent','%FScategory:','force',false);
%%

%files2convert={'htmlwriteFS','makecontentsfileFS', 'mreadFS', 'publishFS', 'publishFunctionAlpha', 'publishFunctionCate', 'xmlcreateFS'};
%filePaths={'C:\FSDA\utilities_help','C:\FSDA\utilities_help', 'C:\FSDA\utilities_help', 'C:\FSDA\utilities_help', 'C:\FSDA\utilities_help', 'C:\FSDA\utilities_help', 'C:\FSDA\utilities_help'};


nfiles=size(outTest,1);

for i=176 % 4:nfiles
    disp(['Iterazione' num2str(i)])
    try
        tmp=publishFS([outTest{i,9} '\' outTest{i,1}],'evalCode',true,'write2file',true);
        % tmp=publishFS([outTest{i,9} '\' outTest{i,1}],'evalCode',true,'Display','iter-detailed','webhelp',false,'outputDir','C:\FSDA\webhelpfiles\FSDA','imagesDir','C:\FSDA\webhelpfiles\FSDA\images','write2file',true);
        %tmp=publishFS([filePaths{i} '\' files2convert{i} '.m'],'evalCode',true,'Display','iter-detailed','webhelp',true,'outputDir','C:\FSDA\webhelpfiles\FSDA','imagesDir','C:\FSDA\webhelpfiles\FSDA\images');
    catch
        
        jjj=1;
        error(['jjjjjjjjjjjjj' num2str(i)])
    end
end
