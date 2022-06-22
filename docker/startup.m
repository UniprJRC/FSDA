function startup

% Setup default proxy settings based on the environment variables that
% we will have set in the run.sh script
host = getenv('MW_PROXY_HOST');
port = getenv('MW_PROXY_PORT');
if ~isempty(host) && ~isempty(port)
    % Replace the deprecated JAVA API with a wrapper
    matlab.net.internal.copyProxySettingsFromEnvironment();
end

% Make sure we do not end up in another folder after running this!    
cdir = pwd;
c = onCleanup(@() cd(cdir));

[~, FSDAroot] = system("cat /opt/fsda/fsda-location.txt");
cd(FSDAroot)

% ALWAYS DO THE FOLLOWING
%% 2) ADD relevant FSDA paths to MATLAB path
try
    addFSDA2path
 catch
    disp('Unknown error when trying to add FSDA folders to MATLAB path')
    disp('File: addFSDA2path could not run')
end



%% 3) Launch buildocsearchdb - this just does not seem to work!
if  matlab.ui.internal.hasDisplay && ...
    matlab.internal.lang.capability.Capability.isSupported(matlab.internal.lang.capability.Capability.ModalDialogs)
    
    folderwithSearchableDatabase=[pwd filesep 'helpfiles' filesep 'pointersHTML'];
    try
        builddocsearchdb(folderwithSearchableDatabase)
        disp('FSDA searchable database correctly added')
    catch
        disp('Unknown error when trying to run MATLAB routine builddocsearchdb')
        disp(['in folder:  '  folderwithSearchableDatabase])
    end
end

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

end