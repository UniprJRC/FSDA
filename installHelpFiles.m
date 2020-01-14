% This file copies the HTML files which are in subfolder
%       (FSDA path)/helpfiles/FSDA
% inside
%       (MATLAB docroot)/help/FSDA
% 
% Note that to properly copy these file under windows, it may be necessary to have
% administrator privileges (or to run MATLAB with administrator privileges)
%


%% Beginning of code
% Store current folder so that after the execution we go back to this
% folder
CurrentFolder=pwd;

% Navigate to FSDA main folder
FileName='addFSDA2path';
FullPath=which(FileName);
if isempty(FullPath)
    error('FSDA:setup:WrongLocation','In order to properly run this file please navigate to the main folder of FSDA')
else
    %Navigate to the main folder of FSDA
    FSDAroot=fileparts(FullPath);
    cd(FSDAroot)
end

%% Copy all FSDA .html files inside (MATLAB docroot)/FSDA
fsep=filesep;

try
    source=['helpfiles' filesep 'FSDA'];
    destination=[docroot filesep 'FSDA'];
    disp('Copying all FSDA documentation files which are in folder')
    disp([FSDAroot filesep source])
    disp('into')
    disp(destination)
    disp('------------------------')
    
    
    if ispc
        % create a shell application
        uac = actxserver('Shell.Application');
        
        
        % Create temporary file named copy_FSDA_help_files.bat
        stringToIncludeInBatFile=['robocopy /E "' FSDAroot filesep source '"  "' destination '"'];
        
        file1ID=fopen('copy_FSDA_help_files.bat','w');
        fprintf(file1ID,'%s',stringToIncludeInBatFile);
        fclose('all');
        
        % run copy_FSDA_help_files,bat file with admin privileges
        uac.ShellExecute('copy_FSDA_help_files.bat', 'ELEV', '', 'runas', 1);
       
        pause(5);
        % delete temporary file copy_FSDA_help_files.bat
        delete copy_FSDA_help_files.bat
        
    else
        % Mac users
        [status,msg]=copyfile(source,destination,'f');
        if status ==1
            disp('HTML FSDA documentation files correctly copied')
        else
            disp('Due to write permission problems HTML files in:')
            disp([pwd filesep source])
            disp('could not be copied inside folder')
            disp(destination)
            disp('To solve the problem, please run MATLAB as administrator')
            disp('or manually copy the files, otherwise the HTML FSDA help files will not be visible')
            warning('FSDA:startup:NotCopied','Impossible to copy FSDA HTML documentation files')
        end
    end
catch
    disp('Unknown error when trying to copy the HTML FSDA documentation files from')
    disp([pwd filesep source])
    disp('to:')
    disp(destination)
end

% Create searchable databaset (not necessary)
% subPointersFolder=[filesep 'helpfiles' filesep 'pointersHTML'];
% folderwithSearchableDatabase=[FSDAroot subPointersFolder];
%
%
% % Launch buildocsearchdb
% % run builddocsearchdb in subfolder pointersHTML
% try
%     builddocsearchdb(folderwithSearchableDatabase)
%     disp('FSDA searchable database correctly added')
% catch
%     disp('Unknown error when trying to run MATLAB routine builddocsearchdb')
%     disp(['in folder:  '  folderwithSearchableDatabase])
% end

cd(CurrentFolder)
clearvars
