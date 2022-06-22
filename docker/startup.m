function startup

% Make sure we don't end up in another folder after running this!    
cdir = pwd;
c = onCleanup(@() cd(cdir));

FSDAroot = "/opt/fsda/FSDA-" + version('-release');
cd(FSDAroot)

% ALWAYS DO THE FOLLOWING
%% 2) ADD relevant FSDA paths to MATLAB path
try
    addFSDA2path
 catch
    disp('Unknown error when trying to add FSDA folders to MATLAB path')
    disp('File: addFSDA2path could not run')
end



%% 3) Launch buildocsearchdb - this just doesn't seem to work!
% folderwithSearchableDatabase=[pwd filesep 'helpfiles' filesep 'pointersHTML'];
% try
%     builddocsearchdb(folderwithSearchableDatabase)
%     disp('FSDA searchable database correctly added')
% catch
%     disp('Unknown error when trying to run MATLAB routine builddocsearchdb')
%     disp(['in folder:  '  folderwithSearchableDatabase])
% end

oneTimeFile = string(prefdir) + filesep + "FSDA-one-time-file";

if exist(oneTimeFile, 'file')
    return
end

%% 4) Install the apps
try
    matlab.apputil.install('brushRES');
catch
    disp('Unknown error when trying to install brushRES app')
end
try
    matlab.apputil.install('brushFAN');
catch
    disp('Unknown error when trying to install brushFAN app')
end

try
    matlab.apputil.install('brushROB');
catch
    disp('Unknown error when trying to install brushROB app')
end

system("touch " + oneTimeFile);