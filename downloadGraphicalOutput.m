
%% Preamble
% The purpose of this file is to download and extract the graphical output
% of the html examples of FSDA contained in file
% http://rosa.unipr.it/FSDA/images.zip into
% (MATLAB docroot)/FSDA/images
% if this path exists
% or to
% (FSDA root)/helpfiles/FSDA/images.
% if the previous path does not exist.
%
% REMARK: this file is necessary just from those who have downloaded FSDA
% from Mathworks file exchange due to the limitation of 20MB size of the
% files inside fileexchange. Those who have downloaded FSDA from github do
% not need to run this file.

%% Beginning of code

% url from which to download file (including file name)
url = 'http://rosa.unipr.it/FSDA/images.zip';

% Set the timeout value to Inf so that the connection does not time out.
options = weboptions('Timeout',Inf);

% Define output folder
fsep=filesep;
% generic outputfolder: works on all OSes
outputfolder=[docroot filesep 'FSDA' filesep 'images'];



if exist(outputfolder,'dir') == 0
    % check  that output folder. If output folder already exists this means
    % that user has already copied there all the html documentation files
    % and simply needs to extract the graphical output
    % in this case the HTML documentation files have not been copied yet
    % and therefore the graphical output will not be extracted into the
    % folder
    
    %disp('the HTML documentation files have not been copied yet ')
    %disp('and therefore the graphical output will not be extracted into the folder.')
    %disp('please run ')
    %error('FSDA:setup:WrongLocation','In order to properly run this file please navigate to the main folder of FSDA')
    installHelpFiles;
    
else
    
end

try
    if ispc
        % save currebn folder location to restore later
        oldFolder = pwd;
        
        % locate main FSDA folder
        FileName='addFSDA2path';
        FullPath=which(FileName);
        
        % Navigate to the main folder of FSDA
        FSDAroot=fileparts(FullPath);
        cd(FSDAroot)
        
        % create a temp folder inside main folder of FSDA
        if exist('imgtemp','dir') == 0
            mkdir 'imgtemp';
        end
        % path to extract all contents of zip file
        % before copying its content in outputfolder dir
        % with admin privileges
        IMGtmpfolder= [pwd filesep 'imgtemp'];
        
        % create a shell application
        uac = actxserver('Shell.Application');
        
        % unzip file into the (temp) outputfolder
        unzip(url, IMGtmpfolder)
        
        % create a batch file to copy all image files
        stringToIncludeInBatFile=['robocopy /E "' IMGtmpfolder '"  "' outputfolder '"'];
        
        file1ID=fopen('copy_FSDA_img_files.bat','w');
        fprintf(file1ID,'%s',stringToIncludeInBatFile);
        fclose('all');
        
        % run batch file with admin privileges
        uac.ShellExecute('copy_FSDA_img_files.bat', 'ELEV', '', 'runas', 1);
        
        pause(5);
        % delete temporary file copy_FSDA_img_files.bat
        delete copy_FSDA_img_files.bat
        
        % remove the temp folder with all image files
        [status, message] = rmdir('imgtemp', 's');
        
        % go back
        cd(oldFolder);
        
    elseif ismac
        % MACOS, no admin privileges issues, just
        % unzip file into outputfolder
        unzip(url,outputfolder)
        disp('Files correctly extracted')
        disp('Now graphical output of the FSDA files is also visible locally')
    else
        % Linux OS
        % from the terminal, type: sudo chown -R $LOGNAME: ~/MATLAB/help/
        % under '/home/user/MATLAB/help/FSDA/images'
        % unzip file into outputfolder
        try
            system('cd ~')
            system('sudo chown -R $LOGNAME: ~/MATLAB/help/')
            unzip(url,outputfolder)
        catch
            disp('Due to write permission problems zip file in:')
            disp(url)
            disp('could not be extracted inside folder')
            disp(outputfolder)
            disp('To solve the problem, please run MATLAB as administrator')
            disp('or manually from the terminal, type: "sudo chown -R $LOGNAME: ~/MATLAB/help/", otherwise the HTML FSDA graphical output will not be visible locally')
            warning('FSDA:dowloadGraphicalOutput:NotExtracted','Impossible to extract FSDA graphical output')
        end
        disp('Files correctly extracted')
        disp('Now graphical output of the FSDA files is also visible locally')
    end
    disp('Files correctly extracted')
    disp('Now graphical output of the FSDA files is also visible locally')
catch
    disp('Due to write permission problems zip file in:')
    disp(url)
    disp('could not be extracted inside folder')
    disp(outputfolder)
    disp('To solve the problem, please run MATLAB as administrator')
    disp('or manually upzip the files, otherwise the HTML FSDA graphical output will not be visible locally')
    warning('FSDA:dowloadGraphicalOutput:NotExtracted','Impossible to extract FSDA graphical output')
    
end




