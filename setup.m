%% This file performs the tasks to run
% More specifically
% 1) It copies all the .html files inside (docroot/FSDA)
% 2) Run file addFSDA2path.m
% 3) Launch buildocsearchdb
% 4) Install the apps

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

