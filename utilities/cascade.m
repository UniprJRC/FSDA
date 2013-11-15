function cascade
%cascade is a third party function used by some FSDA demo/example functions
%(brushRES.m and example_regression.m in folder examples) that
%produce plots which should not direcly overlap. 
%
%
% Code, comments and authorship follows.
%
%
%CASCADE Cascade existing figures so that they don't directly overlap
%   CASCADE takes and returns no arguments.  This function will cascade as
%   many figures as will fit the height/width of the screen.  If there are
%   more figures than can cascade in a screen, those additional figures are
%   left in their original position.
%
%   Author: Isaac Noh
%   Copyright The MathWorks, Inc.
%   November 2007

% Existing Figures
figs = findobj(0,'Type','figure');
figs = sort(figs);

% Demos are not subject to repositioning and are removed from figure list
demo=((findobj('type','figure','Tag','demo')));
if ~isempty(demo)
    figs(intersect(figs,demo))=[]; 
end

% Size of Entire Screen
ss = get(0,'ScreenSize');

units = get(figs,'Units');
set(figs,'Units','pixels')

for n = 1:length(figs)
    pos = get(figs(n),'Position');
    if n == 1
        bot = ss(4) - pos(4) - 140;
        set(figs(n),'Position',[0 bot pos(3:4)]);
    else
        pPos = get(figs(n-1),'Position');
        left = pPos(1) + 150;
        bot = pPos(2) - 70;
        if ((left + pos(3)) > ss(3)) || (bot < 0)
            break
        end
        set(figs(n),'Position',[left bot  pos(3:4)]);
    end
    figure(figs(n));
end
if length(figs)>1
    set(figs,{'Units'},units);
end