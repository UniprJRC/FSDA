function status = closereqFS(~,~)
%closereqFS is a user-defined close request function which displays a question dialog box
%
%
% Copyright 2008-2014.
% Written by FSDA team
%
%

status = 0;
selection = questdlg('Close this figure and all the other linked figures and exit brushing mode',...
    'Close figure',...
    'Yes','No','Yes');
switch selection,
    case 'Yes',
        delete(get(0,'CurrentFigure')); %delete(src);
        status = 1;
    case 'No'
        return
end
end