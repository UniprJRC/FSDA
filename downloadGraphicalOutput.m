
%% Preamble
% The purpose of this file is to download and extract the graphical output
% of the html examples of FSDA contained in file
% http://rosa.unipr.it/FSDA/outputHTMLfiles.zip into
% (MATLAB docroot)/FSDA//images
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
outputfolder=[docroot filesep 'FSDA' filesep 'images'];

if exist(outputfolder,'dir')==7
    % check  that output folder. If output folder already exists this means
    % that user has already copied there all the html documentation files
    % and simply needs to extract the graphical output
    
    try
        % unzip file into outputfolder
        unzip(url,outputfolder)
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
else
    % in this case the HTML documentation files have not been copied yet
    % and therefore the graphical output will be extracted into
    
    FileName='addFSDA2path';
    FullPath=which(FileName);
    if isempty(FullPath)
        error('FSDA:setup:WrongLocation','In order to properly run this file please navigate to the main folder of FSDA')
    else
        %Navigate to the main folder of FSDA
        FSDAroot=fileparts(FullPath);
        %  cd(FSDAroot)
    end
    outputfolder=[FSDAroot filesep 'helpfiles' filesep 'FSDA' filesep 'images'];
    
    try
        % unzip file into outputfolder
        unzip(url,outputfolder)
        disp(['Files correctly extracted into '  outputfolder])
        disp('Now it is necessary to copy manually folder')
        disp([FSDAroot filesep 'helpfiles' filesep 'FSDA'])
        disp('in')
        disp([docroot filesep help 'FSDA'])
        disp('-----------------------------')
        disp('Alternatively, run routine checkFSDAsetup.m')
        disp('which is located in the main folder of FSDA')
    catch
        disp('Unknown error when trying to unzip file into folder')
        disp(outputfolder)
    end
    
end


