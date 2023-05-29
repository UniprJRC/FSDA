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
if strcmp(PersonalvalueDOC,'WEB')
    disp('In order to be able to see the FSDA documentation locally,')
    disp('it is necessary in ''Home|Preferences'' to set ''Documentation Location''')
    disp('to ''Installed Locally''')
    [~,nameFile,ext]=fileparts(varargin);
    varargin{:}=['http://rosa.unipr.it/FSDA/' nameFile ext];
end
web(varargin{:});
end
