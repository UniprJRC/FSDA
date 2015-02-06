function position(plmain)
%position controls the position of the open figures
%
% Required input arguments:
%
% Optional input arguments:
%
%   plmain, the handle of a 'main' figure to be positioned at the top-left
%   side of the screen, which is supposed to be the position attracting
%   first the attention of a user. In absence of the option, plmain is set
%   to be the smaller handle, which usually corresponds to the first
%   created figure. 
%
% A number of relevant FSDA plots are positioned according to a predefined
% layout. 
%
% Reminder: the position property format is [left, bottom, width, height]
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
% Last modified 06-Feb-2015
%


%% Beginning of code

% find all open figures and if there are no figures do nothing. This check
% is necessary because the figure with handle plmain could have been closed
% before the call of position.
openfigs = findobj(allchild(0),'flat','Visible','on');
if isempty(openfigs),
    return
else
    if nargin == 0 
        plmain = min(openfigs);
    end
    %Ensure root units are pixels and get the size of the screen
    set(0,'Units','pixels') ;
    scnsize = get(0,'ScreenSize');
    scrwidth = scnsize(3); scrheight = scnsize(4); halfscrheight = scrheight/2;
    % The size and coordinates of the plot under brushing are not used,
    % because the width/height ratio is fixed.
    % plmainpos = get(plmain,'Position');plmainwidth = plmainpos(3);
    % plmainheight = plmainpos(4);
    whratio=560/420;
    % Set an extra space equal to the border width between the main figure
    % and other figures
    edge = 2;
    % Fraction of the upper screen which will be occupied by the main plot
    fraction = 0.5; 
end

%% Reposition the main figure. 
% The main figure is usually a plot under brushing
% posmain = [1 halfscrheight halfscrheight*whratio halfscrheight];
posmain = [edge,...
            scrheight * (1-fraction),...
            halfscrheight*whratio - edge,...
            scrheight * fraction];
set(plmain,'WindowStyle','normal','OuterPosition',posmain);

%% Reposition all other figures

% ensure first that they are not docked, otherwise re-positioning will not work
dockedfigs = findobj(openfigs,'flat','WindowStyle','docked');
if ~isempty(dockedfigs),
    set(dockedfigs,'WindowStyle','normal');
end

% Set the new position of a GUI or demo
demo=((findobj(openfigs,'flat','Tag','demo')));
if isempty(demo)
    demo=[];
else
%   newdemopos = [halfscrheight*whratio halfscrheight halfscrheight*whratio halfscrheight-50];
    newdemopos = [posmain(3) + edge,...
        posmain(2),...
        posmain(3),...
        posmain(4)-50];
    set(demo,'Position',newdemopos);
end

% Check the *linkable* FS plots that are open and get their handles
% residuals. In findobj, "'type','figure'" is replaced with
% "openfigs,'flat'" to accelerate the search of the right handle.

% REGRESSION PLOTS
% yXplot
plyX=((findobj(openfigs,'flat','Tag','pl_yX')));
if isempty(plyX), plyX=[]; end;
% resfwdplot
plres=((findobj(openfigs,'flat','Tag','pl_resfwd'))); 
if isempty(plres), plres=[]; end;
% resfwdplot
plresindex=((findobj(openfigs,'flat','Tag','pl_resindex'))); 
if isempty(plresindex), plresindex=[]; end;
% levfwdplot
pllev=((findobj(openfigs,'flat','Tag','pl_levfwd')));
if isempty(pllev), pllev=[]; end;
% mdrplot
plmdr=((findobj(openfigs,'flat','Tag','pl_mdr')));
if isempty(plmdr), plmdr=[]; end;
% fanplot
plfan=((findobj(openfigs,'flat','Tag','pl_fan')));
if isempty(plfan), plfan=[]; end;

% MULTIVARIATE PLOTS  
% malfwdplot
plmalfwd=((findobj(openfigs,'flat','Tag','pl_malfwd')));
if isempty(plmalfwd), plmalfwd=[]; end;
% malindexplot
plmalindex=((findobj(openfigs,'flat','Tag','pl_malindex')));
if isempty(plmalindex), plmalindex=[]; end;
% mmdplot
plmmd=((findobj(openfigs,'flat','Tag','pl_mmd')));
if isempty(plmmd), plmmd=[]; end;
% spmplot
plspm=((findobj(openfigs,'flat','Tag','pl_spm')));
if isempty(plspm), plspm=[]; end;

%Set the new position of the FS plots linked to that under brushing
lkplots = horzcat(plmdr,plres,plresindex,plyX,pllev,plfan,...
            plmalindex,plmalfwd,plmmd,plspm);
lkplots(lkplots==plmain)=[];
nlkplots = length(lkplots);
width = scrwidth / (nlkplots + 1);
left = zeros(nlkplots,1);
for i=1:nlkplots
    left(i) = floor(scrwidth-(i*width));
end
fb = 30 ; %position from bottom
switch nlkplots;
    case 0
    case 1
        poslkplots = cat(2,left,[fb width scrheight/2-fb]);
        set(lkplots,'OuterPosition',poslkplots);
    otherwise
        poslkplots = cat(2,left,repmat([fb width scrheight/2-fb],nlkplots,1));
        cposlkplots=num2cell(poslkplots,2);
        set(lkplots,{'OuterPosition'},cposlkplots);
end

%Set the new position and size of any other open plot
otherfigs = openfigs';
otherfigs(ismember(otherfigs,lkplots)) = [];
otherfigs(otherfigs==plmain)=[];
if ~isempty(demo)
    otherfigs(otherfigs==demo)=[];
end

notherfigs = length(otherfigs);
switch notherfigs;
    case 0
    case 1
        posotherfigs = [scrwidth-(scrwidth/3) scrheight-scrheight/3 scrwidth/3 scrheight/3];
        set(otherfigs,'OuterPosition',posotherfigs);
    otherwise
        displacements = 20*ones(notherfigs,notherfigs).*repmat((0:notherfigs-1)',1,notherfigs);
        posotherfigs = [scrwidth-(scrwidth/3) scrheight-scrheight/3 scrwidth/3 scrheight/3];
        posotherfigs = repmat(posotherfigs,notherfigs,1);
        posotherfigs(:,1:notherfigs) = posotherfigs(:,1:notherfigs) - displacements;
        cposotherfigs=num2cell(posotherfigs,2);
        set(otherfigs,{'OuterPosition'},cposotherfigs);
end

end

