function webFS(varargin)
% We check whether the user has chosen as
% Documentation Location "web" or "locally", depending on that we show the
% FSDA documentation on rosa.unipr.it/FSDA or locally

%% Beginning of code

s=settings;

if verLessThanFS([9 10]) % < R2021a
    web(varargin{:});
elseif verLessThanFS([9 14]) % From R2021a to R2022b
    try
        PersonalvalueDOC=s.matlab.help.DocCenterLocation21a.PersonalValue;
    catch
        PersonalvalueDOC='WEB';
    end

else % From 2023a
    try
        PersonalvalueDOC=s.matlab.help.DocLocation.PersonalValue;
    catch
        PersonalvalueDOC='WEB';
    end
end

    [~,nameFile,ext]=fileparts(varargin{:});

if strcmp(PersonalvalueDOC,'WEB')
    disp('In order to be able to see the FSDA documentation locally,')
    disp('it is necessary in ''Home|Preferences'' to set ''Documentation Location''')
    disp('to ''Installed Locally''')
    varargin{:}=['http://rosa.unipr.it/FSDA/' nameFile ext];
else
    % Make sure that the file exists in the destination path
    % If the destination file does not exist, call routine installHelpFiles
    if exist(varargin{:},'file') ==0
        status=installHelpFiles;
        if status==0
            % If routine installHelpFiles failed than redirect to rosa
            % website for the documentation
            varargin{:}=['http://rosa.unipr.it/FSDA/' nameFile ext];
        end
    end
end

web(varargin{:});
end


