% This file makes sure that the HTML files which are in subfolder
% (FSDA path)/helpfiles/FSDA
% are copied inside (docroot)/help/FSDA
% Note that to properly copy these file it is necessary to have
% administrator privileges (or to run MATLAB with adminitrator privileges)


%% Copy all FSDA .html files inside docroot/FSDA
fsep=filesep;
if exist([docroot filesep 'FSDA'],'dir')~=7
    
    try
        source=['helpfiles' filesep 'FSDA'];
        destination=[docroot filesep 'FSDA'];
        [status,msg]=copyfile(source,destination,'f');
        if status ==1
            disp('HTML FSDA documentation files correctly copied')
            % Launch buildocsearchdb
            % FIrst navigate to FSDA main folder
            FileName='addFSDA2path';
            FullPath=which(FileName);
            if isempty(FullPath)
                error('FSDA:setup:WrongLocation','In order to properly run this file please navigate to the main folder of FSDA')
            else
                %Navigate to the main folder of FSDA
                FSDAroot=fileparts(FullPath);
                cd(FSDAroot)
            end
            % run builddocsearchdb in subfolder pointersHTML
            folderwithSearchableDatabase=[pwd filesep 'helpfiles' filesep 'pointersHTML'];
            try
                builddocsearchdb(folderwithSearchableDatabase)
                disp('FSDA searchable database correctly added')
            catch
                disp('Unknown error when trying to run MATLAB routine builddocsearchdb')
                disp(['in folder:  '  folderwithSearchableDatabase])
            end
            
        else
            disp('Due to write permission problems HTML files in:')
            disp([pwd filesep source])
            disp('could not be copied inside folder')
            disp(destination)
            disp('To solve the problem, please run MATLAB as administrator')
            disp('or manually copy the files, otherwise the HTML FSDA help files will not be visible')
            warning('FSDA:startup:NotCopied','Impossible to copy FSDA HTML documentation files')
        end
    catch
        disp('Unknown error when trying to copy the HTML FSDA documentation files from')
        disp([pwd filesep source])
        disp('to:')
        disp(destination)
    end
    
end
