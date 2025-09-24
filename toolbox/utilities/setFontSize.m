function setFontSize(fontSize)
% setFontSize Sets the font size for MATLAB code in the editor and Command Window
%
% Syntax:
%   setFontSize()           % Uses default size (10)
%   setFontSize(fontSize)   % Uses specified size
%
% Input:
%   fontSize - (optional) Desired font size (positive number)
%              Default value: 10
%
% Example:
%   setFontSize(14)         % Sets font to 14pt
%   setFontSize()           % Sets font to 10pt (default)
%

    % Default value if not specified
    if nargin == 0
        fontSize = 10;
    end
    
    % Input validation
    if ~isnumeric(fontSize) || ~isscalar(fontSize) || fontSize <= 0
        error('fontSize must be a positive number');
    end
    
    % Access MATLAB settings and set code font size
    s = settings;
    s.matlab.fonts.codefont.Size.TemporaryValue = fontSize;
    
    % Display confirmation message
    fprintf('Code font size set to: %.0f\n', fontSize);
    
end