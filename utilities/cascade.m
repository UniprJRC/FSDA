function cascade()
%cascade is a third party function used in FSDA demos and examples
% producing plots that should not direcly overlap.
%
%<a href="matlab: docsearchFS('cascade')">Link to the help page for this function</a>
%
% function cascade cascades existing figures so that they don't directly overlap
%   cascade takes and returns no arguments.  This function will cascade as
%   many figures as will fit the height/width of the screen.  If there are
%   more figures than can cascade in a screen, those additional figures are
%   left in their original position.
%
% Required input arguments:
%
%
% Optional input arguments:
%
% Output:
%
%
%
%  See also brushRES.m 
%
% References:
%
%   Tufte E.R. (1983), "The visual display of quantitative information",
%   Graphics Press, Cheshire.

%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('cascade')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples
%
%{
    % Cascade multiple figures.
    close all;
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);
    yXplot(y,X);
    [out]=LXS(y,X,'nsamp',10000);
    [out]=FSReda(y,X,out.bs);
    out1=out;
    out1.RES=out.RES.^2;
    resfwdplot(out1);
    levfwdplot(out1);
    resindexplot(out1.RES);
    plot(out.S2(:,1),out.S2(:,2)); title('Plot of R^2');
    cascade;
%}

%% Beginning of code

% Existing Figures
figs = findobj(0,'Type','figure');

% Demos are not subject to repositioning and are removed from figure list
finddemo=strcmp(get(figs,'Tag'),'demo');

if sum(finddemo)>0
    figs(finddemo)=[];
end

% Size of Entire Screen
ss = get(0,'ScreenSize');

units = get(figs,'Units');
set(figs,'Units','pixels')

for j = 1:length(figs)
    pos = get(figs(j),'Position');
    if j == 1
        bot = ss(4) - pos(4) - 140;
        set(figs(j),'Position',[0 bot pos(3:4)]);
    else
        pPos = get(figs(j-1),'Position');
        left = pPos(1) + 150;
        bot = pPos(2) - 70;
        if ((left + pos(3)) > ss(3)) || (bot < 0)
            break
        end
        set(figs(j),'Position',[left bot  pos(3:4)]);
    end
    figure(figs(j));
end
if length(figs)>1
    set(figs,{'Units'},units);
end
end
%FScategory:UTIGEN