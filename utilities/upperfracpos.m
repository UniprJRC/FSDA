function upperfracpos(hfigl , hfigr , fraction)
%upperfracpos positions two figures on the upper part of the screen.
%
%
%upperfracpos positions two figures with handles hfigl and hfigr
%respectively on the mid half-right and mid half-left of the upper part of
%the screen. The fraction of the screen which will be occupied is
%determined by option parameter 'fraction', which is a number between 0 and
%1. For example, with fraction=0.5 the two figures will occupy the upper
%half of the screen and with  fraction=0.3 they will occupy the upper third
%of the screen.
%
% Copyright 2008-2015.
% Written by FSDA team
%
% Last modified 06-Feb-2015
%
% Examples
%{
%% Example of use of upperfracpos
close all;

% create two figures, rescale and position them
hfigl = figure; plot(sin(rand(10,1)),'r'); title('goes on left');
hfigr = figure; plot(cos(rand(10,1)),'b'); title('goes on right');
upperfracpos(hfigl , hfigr , 0.5);

% now rescale the figures to a smaller proportion
upperfracpos(hfigl , hfigr , 0.2);

% this is just to bring the rescaled figures in the screen foreground
figure(hfigl); figure(hfigr);
%}

%% Beginning of code
if nargin < 3
    fraction = 0.5;
end

%Ensure root units are pixels and get the size of the screen
set(0,'Units','pixels') ;
scnsize = get(0,'ScreenSize');

%The figure Position property only includes the drawable window part,
%without borders. Obtain the full window size from OuterPosition
position = get(hfigl,'Position');
outerpos = get(hfigl,'OuterPosition');
borders = outerpos - position;

% Define the desired size and location of the figures. Leave a space equal
% to their border width between them

edge = -borders(1)/2;
pos1 = [edge,...
    scnsize(4) * (1-fraction),...
    scnsize(3)/2 - edge,...
    scnsize(4)*fraction];
pos2 = [scnsize(3)/2 + edge,...
    pos1(2),...
    pos1(3),...
    pos1(4)];

% Reposition the two figures by changing their OuterPosition properties
set(hfigl,'OuterPosition',pos1);
set(hfigr,'OuterPosition',pos2);

end