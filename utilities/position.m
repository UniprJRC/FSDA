function position(plmain)
%position controls the position of the open figures
%
%<a href="matlab: docsearchFS('position')">Link to the help function</a>
%
% A number of relevant FSDA plots are positioned according to a predefined
% layout.
% Reminder: the position property format is [left, bottom, width, height]
%
% Required input arguments:
%
%   plmain:     Figure handle. Scalar. The handle of a 'main' figure to be
%               positioned at the top-left side of the screen, which is
%               supposed to be the position attracting first the attention
%               of a user.
%                -  If plmain not given, or it is empty, it is set to be
%                   the smaller handle, which normally is the handle of the
%                   first created figure.
%                -  If plmain is set to zero (0), then function cascade is
%                   applied.
%
%
% Optional input arguments:
%
%  Output:
%
% See also: cascade
%
% References:
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('position')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Three plots, with a recognized tag. They are organized properly in
    % the screen area.

    Y1=load('geyser2.txt');
    Y2=load('fishery.txt');
    Y3=load('M5data.txt');
    figure('Tag','pl_spm'); spmplot(Y1); hmain=gcf;
    figure('Tag','pl_spm'); spmplot(Y2);
    figure('Tag','pl_spm'); spmplot(Y3);
    position(hmain);
%}
%{
    % The three plots have now an unknown tag. In this case, we assume that 
    % the plots are not relevant and are therefore put in a non-interesting
    % screen area (top-right).
    
    close all
    figure('Tag','aaaaa'); spmplot(Y1); hmain=gcf;
    figure('Tag','bbbbb'); spmplot(Y2);
    figure('Tag',''); spmplot(Y3);
    position(hmain);
%}


%% Beginning of code

% find all open figures
openfigs = findobj(allchild(0),'flat','Visible','on');

% in case of 0 or 1 figures, do nothing
if isempty(openfigs) || length(double(openfigs))==1
    return
end

%  if the input argument 'plmain' is zero, there is no 'main'
%  figure and simple cascade is used.
if nargin == 1 && plmain == 0
    cascade;
    return
end

if nargin == 0
    %  If 'plmain' is not indicated (i.e. no input arguments) the 'main'
    %  figure positioned at the top-left side of the screen is set to be
    %  the one with smaller handle, which usually corresponds to the first
    %  created figure.
    plmain = min(double(openfigs));
end

%% set the geometry of the space where figures will be positioned

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
if ~isempty(dockedfigs)
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

% Check the *linkable* FS plots that are open and get their handles. In
% findobj, "'type','figure'" is replaced with "openfigs,'flat'" to
% accelerate the search of the right handle.

% REGRESSION PLOTS
% yXplot
plyX=((findobj(openfigs,'flat','Tag','pl_yX')));
if isempty(plyX), plyX=[]; end
% resfwdplot
plres=((findobj(openfigs,'flat','Tag','pl_resfwd')));
if isempty(plres), plres=[]; end
% resfwdplot
plresindex=((findobj(openfigs,'flat','Tag','pl_resindex')));
if isempty(plresindex), plresindex=[]; end
% levfwdplot
pllev=((findobj(openfigs,'flat','Tag','pl_levfwd')));
if isempty(pllev), pllev=[]; end
% mdrplot
plmdr=((findobj(openfigs,'flat','Tag','pl_mdr')));
if isempty(plmdr), plmdr=[]; end
% fanplot
plfan=((findobj(openfigs,'flat','Tag','pl_fan')));
if isempty(plfan), plfan=[]; end

% MULTIVARIATE PLOTS
% malfwdplot
plmalfwd=((findobj(openfigs,'flat','Tag','pl_malfwd')));
if isempty(plmalfwd), plmalfwd=[]; end
% malindexplot
plmalindex=((findobj(openfigs,'flat','Tag','pl_malindex')));
if isempty(plmalindex), plmalindex=[]; end
% mmdplot
plmmd=((findobj(openfigs,'flat','Tag','pl_mmd')));
if isempty(plmmd), plmmd=[]; end
% spmplot
plspm=((findobj(openfigs,'flat','Tag','pl_spm')));
if isempty(plspm), plspm=[]; end

%Set the new position of the FS plots linked to that under brushing
lkplots = horzcat(plmdr,plres,plresindex,plyX,pllev,plfan,...
    plmalindex,plmalfwd,plmmd,plspm);
lkplots(lkplots==plmain)=[];

nlkplots = length(lkplots);
if nlkplots > 0
    % the width of the plots
    switch nlkplots
        case 1
            lkplots_pos = get(lkplots,'OuterPosition');
            width = lkplots_pos(3);
        case {2,3,4}
            lkplots_pos = cell2mat(get(lkplots,'OuterPosition'));
            lkplots_width = min(min(lkplots_pos(:,3)),ceil(scrwidth / nlkplots));
            width = lkplots_width;
        otherwise
            width = ceil(scrwidth / 4);
    end
    heigth = ceil(width/whratio);
    
    % the left position of the plots depends on their width
    left = ((1:nlkplots)'-1)*width;
    left = mod(left,scrwidth);
    
    % set the displacement from the bottom of the screen
    fb = 30 ;
    if nlkplots > 1
        fb = fb*ones(nlkplots,1);
        for i=1:nlkplots
            nrows = ceil(i/4);
            fb(i) = min(nrows*fb(i)+heigth*(nrows-1) , heigth*3);
        end
    end
    
    % set the new positions
    switch nlkplots
        case 1
            poslkplots = cat(2,left,[fb width scrheight/2-fb]);
            set(lkplots,'OuterPosition',poslkplots);
        otherwise
            poslkplots = cat(2,left,fb,repmat(width,nlkplots,1),repmat(heigth,nlkplots,1));
            cposlkplots=num2cell(poslkplots,2);
            set(lkplots,{'OuterPosition'},cposlkplots);
    end
    
end

%Set the new position and size of any other open plot
otherfigs = openfigs';
otherfigs(ismember(otherfigs,lkplots)) = [];
otherfigs(otherfigs==plmain)=[];
if ~isempty(demo)
    otherfigs(otherfigs==demo)=[];
end

notherfigs = length(otherfigs);
switch notherfigs
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
%FScategory:UTIGEN
