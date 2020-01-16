%% This file performs the necessary tasks to properly run FSDA toolbox
% More specifically, this file
% 1) copies all the .html files inside (docroot/FSDA). In order to
% succesfully perform this operation it is necessary to run MATLAB as
% administrator;
% 2) runs file addFSDA2path.m. This files adds to the MATLAB path all the
% required FSDA folders.
% 3) launches buildocsearchdb in subfolder (main FSDA folder)/helpfiles/pointersHTML. 
% 4) Install the apps. 3 apps brushRES, brushFAN and brushROB which help
% the user to familiarize with the dynamic interactive routines of FSDA.


%% Preliminary step: make sure that the user is inside main FSDA folder
FileName='addFSDA2path';
FullPath=which(FileName);
if isempty(FullPath)
    error('FSDA:setup:WrongLocation','In order to properly run this file please navigate to the main folder of FSDA')
else
    %Navigate to the main folder of FSDA
    FSDAroot=fileparts(FullPath);
    cd(FSDAroot)
end

%% Copy all .html files inside docrootFSDA
fsep=filesep;
try
    source=['helpfiles' filesep 'FSDA'];
    destination=[docroot filesep 'FSDA'];
    [status,msg]=copyfile(source,destination,'f');
    if status ==1
        disp('HTML files correctly copied')
    else
        disp('Due to write permission problems HTML files in:')
        disp([pwd filesep source])
        disp('could not be copied inside folder')
        disp(destination)
        disp('To solve the problem, please run MATLAB as administrator')
        disp('or manually copy the files')
    end
catch
    disp('Unknown error when trying to copy the HTML files from')
    disp([pwd filesep source])
    disp('to:')
    disp(destination)
end

%% 2) ADD relevant FSDA paths to MATLAB path
try
    addFSDA2path
    disp('FSDA added to the MATLAB path')
    disp('The following folders are added to MATLAB path:')
    disp([pwd filesep 'clustering'])
    disp([pwd filesep 'regression'])
    disp([pwd filesep 'multivariate'])
    disp([pwd filesep 'utilities_stat'])
    disp([pwd filesep 'utilities_help'])
    disp([pwd filesep 'FSDAdemos'])
    disp([pwd filesep 'graphics'])
    disp([pwd filesep 'utilities'])
    disp([pwd filesep 'examples'])
    disp([pwd filesep 'combinatorial'])
    disp([pwd filesep 'datasets' filesep 'clustering'])
    disp([pwd filesep 'datasets' filesep 'multivariate'])
    disp([pwd filesep 'datasets' filesep 'regression'])
    disp([pwd filesep 'datasets' filesep 'multivariate_regression'])
    disp('')
    disp('Please do not remove them otherwise the FSDA toolbox will not work.')
catch
    disp('Unknown error when trying to add FSDA folders to MATLAB path')
    disp('File: addFSDA2path could not run')
end

%% 3) Launch buildocsearchdb
folderwithSearchableDatabase=[pwd filesep 'helpfiles' filesep 'pointersHTML'];
try
    builddocsearchdb(folderwithSearchableDatabase)
    disp('FSDA searchable database correctly added')
catch
    disp('Unknown error when trying to run MATLAB routine builddocsearchdb')
    disp(['in folder:  '  folderwithSearchableDatabase])
end

%% 4) Install the apps
try
    matlab.apputil.install('brushRES');
catch
    disp('Unknown error when trying to install brushRES app')
end
try
    matlab.apputil.install('brushFAN')
catch
    disp('Unknown error when trying to install brushFAN app')
end

try
    matlab.apputil.install('brushROB');
catch
    disp('Unknown error when trying to install brushROB app')
end

