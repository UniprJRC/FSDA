function status = closereqFS(~,~)
%closereqFS is a user-defined close request function which displays a question dialog box
%
%
% Copyright 2008-2011.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti 
%            and Vytis Kopustinskas (2009-2010)
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