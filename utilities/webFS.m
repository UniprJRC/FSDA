function webFS(varargin)
% In case of MATLAB 2021A or greater if the user has chosen as
% Documentation Location "web on Mathworks" and the user is searching for a
% FSDA function temporarily change the Domain in order to prevent to search
% in Mathworks web site.
if verLessThanFS([9 10])
    web(varargin{:});
else
    s=settings;
    try
        PersonalvalueDOC=s.matlab.help.DocCenterLocation21a.PersonalValue;
    catch
        PersonalvalueDOC='WEB';
    end

    if strcmp(PersonalvalueDOC,'WEB')
        disp('In order to be able to see locally the FSDA documentation')
        disp('It is necessary in Home|Preferences to set ''Documentation Location''')
        disp('to ''Installed Locally or restart MATLAB''')
        warning('FSDA:webFS:Wrongloc','Please change Documentation Location in order to view HTML FSDA documentation')
        %     % Change in a temporary way the user option to INSTALLED
        s.matlab.help.DocCenterLocation21a.PersonalValue='INSTALLED';
        %     web(varargin{:});
        %     % Restore previous option
        %     s.matlab.help.DocCenterLocation21a.PersonalValue=PersonalvalueDOC;
        %  com.mathworks.mlservices.MLHelpServices.setDocCenterDomain('');
    end
    web(varargin{:});
end
