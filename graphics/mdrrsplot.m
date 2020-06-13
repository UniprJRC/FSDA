function brushedUnits=mdrrsplot(out,varargin)
%mdrrsplot plots the trajectory of minimum deletion residual from random starts
%
%<a href="matlab: docsearchFS('mdrrsplot')">Link to the help function</a>
%
% Required input arguments:
%
%  out :  Structure containing the following fields.
%
%   out.mdrrs = a matrix of size (n-ninit)-by-(nsimul+1)containing the
%               monitoring of minimum deletion
%               residual in each step of the forward search for each of the nsimul random starts.
%               The first
%               column of mdr must contain the fwd search index
%               This matrix can be created using function FSRmdrrs
% out.BBrs   =  3D array of size n-by-n-(init)-by-nsimul containing units
%               forming subset for each random start.
%               This field is necessary if datatooltip is true or databrush
%               is not empty.
%     out.y   =   a vector containing the response (necessary only if
%               option databrush is not empty)
%     out.X   =   a matrix containing the explanatory variables
%               (necessary only if option databrush is not empty)
%
%
% Optional input arguments:
%
%       quant   :   vector containing quantiles for which envelopes have
%                   to be computed. Vector or scalar. The default is to
%                   produce 1%, 50% and 99% envelopes. In other words the
%                   default is quant=[0.01;0.5;0.99];
%               Example - 'quant',[0.01 0.99]
%               Data Types - double
%
%       envm    :   Size of the sample which is
%                   used to superimpose the envelope. Scalar. The default is to add
%                   an envelope based on all the observations (size n
%                   envelope).
%               Example - 'envm',n
%               Data Types - double
%
%       xlimx   :   vector with two elements controlling minimum and
%                   maximum on the x axis. Default value is mdr(1,1)-3 and
%                   mdr(end,1)*1.3
%                   Example - 'xlimx',[20 100]
%                   Data Types - double
%
%       ylimy   :   min and max on the y axis. Vector. Vector with two
%                   elements controlling minimum and
%                   maximum on the y axis. Default value is min(mdr(:,2))
%                   and max(mdr(:,2));
%                   Example - 'ylimy',[2 6]
%                   Data Types - double
%
%       lwdenv  :   Line width of the envelopes. Scalar. Scalar which
%                   controls the width of the lines associated
%                   with the envelopes. Default is lwdenv=1.
%                   Example - 'lwdenv',2
%                   Data Types - double
%
%       tag     :   tag of the plot. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag 'pl_mdrrs'. Notice that if the
%                   program finds a plot which has a tag equal to the one
%                   specified by the user, then the output of the new plot
%                   overwrites the existing one in the same window else a
%                   new window is created
%                   Example - 'tag','mymdrrs'
%                   Data Types - char
%
%   datatooltip :   empty value or structure. The default is datatooltip=''
%                   If datatooltip is not empty the user can use the mouse
%                   in order to have information about the unit seected,
%                   the step in which the unit enters the search and the
%                   associated label. If datatooltip is a structure, it is
%                   possible to control the aspect of the data cursor (see
%                   function datacursormode for more details or the
%                   examples below). The default options of the structure
%                   are DisplayStyle='Window' and SnapToDataVertex='on'
%                   Example - 'datatooltip',1
%                   Data Types - empty value, numeric or structure
%
%    databrush :    interactive mouse brushing. Empty value (default),
%                   scalar or structure.
%                   DATABRUSH IS AN EMPTY VALUE .
%                   If databrush is an empty
%                   value (default), no brushing is done. The activation of
%                   this option (databrush is a scalar or a structure) enables
%                   the user  to select a set of trajectories in the
%                   current plot and to see them highlighted in the spm
%                   (notice that if the spm does not exist it is automatically created).
%                   In addition, units forming subset in the selected steps
%                   selected trajectories can be highlighted in the
%                   monitoring MD plot Note that the window style of the
%                   other figures is set equal to that which contains the
%                   monitoring residual plot. In other words, if the
%                   monitoring residual plot is docked all the other
%                   figures will be docked too.
%                  DATABRUSH IS A SCALAR.
%                   If databrush is a scalar the default selection tool is a
%                   rectangular brush and it is possible to brush only once
%                   (that is persist='').
%                  DATABRUSH IS A STRUCTURE.
%                   If databrush is a structure, it is
%                   possible to use all optional arguments
%                   of function selectdataFS.m and the following optional
%                   argument:
%                  persist. Persist is an empty value or a scalar
%                   containing the strings 'on' or 'off' If persist = 'on'
%                   or 'off' brusing can be done as many time as the user
%                   requires. If persist='on' then the unit(s) currently
%                   brushed are added to those previously brushed. If
%                   persist='off' every time a new brush is performed units
%                   previously brushed are removed. The default value of
%                   persist is '' that is brushing is allowed only once. If
%                   persist is 'on' it is possible, every time a new
%                   brushing is done, to use a different color for the
%                   brushed units
%                  labeladd. If this option is '1', we label the units
%                     of the last selected group with the unit row index in
%                     matrices X and y. The default value is labeladd='',
%                     i.e. no label is added.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%
%       FontSize:   Label font size. Scalar. Scalar which controls the
%                   fontsize of the labels of the
%                   axes. Default value is 12
%                   Example - 'FontSize',14
%                   Data Types - single | double
%
%       ColorTrj:   Color of trajectories. Scalar. Integer which controls 
%                   the color of the trajectories. Default value is 1,
%                   which rotates fixed colors. ColorTrj = 0 produces a
%                   colormap proportional to sum or mdr along a relevant
%                   part of the trajectory. ColorTrj > 1 can be used for
%                   rotating fixed colors for the ColorTrj trajectories
%                   with larger mdr; no more than 7 trajectories will be
%                   considered. ColorTrj > 1 also adds a marker every 10
%                   steps. Note that if the largest mdr are in the final
%                   part of the search (due to a group of outliers), the
%                   peak is not informative and it is therefore not
%                   considered because.
%                   Example - 'ColorTrj',0
%                   Data Types - single | double
%
%    SizeAxesNum:   Size of axes numbers. Scalar. Scalar which controls the
%                   fontsize of the numbers of the axes.
%                   Default value is 10.
%                   Example - 'SizeAxesNum',14
%                   Data Types - single | double
%
%       nameX   :   cell array of strings of length p containing the labels
%                   of the varibles of the regression dataset. If it is empty
%                 	(default) the sequence X1, ..., Xp will be created
%                   automatically
%                   Example - 'nameX',{'Age','Income','Married','Profession'}
%                   Data Types - cell
%
%       namey   :   response label. Character. Character containing the
%                   label of the response.
%                   Example - 'namey','mylabel'}
%                   Data Types - character
%
%       lwd     :   Trajectories line width. Scalar. Scalar which controls
%                   linewidth of the curve which contains the monitoring of
%                   minimum deletion residual.
%                   Default line width=2
%                   Example - 'lwd',3
%                   Data Types - single | double
%
%       titl    :   titel label. Charater. A label for the title (default: '')
%                   Example - 'namey','Plot title'
%                   Data Types - char
%
%       labx    :   x axis title. Character.
%                   A label for the x-axis (default: 'Subset size m').
%                   Example - 'labx','Subset size m'
%                   Data Types - char
%
%       laby    :   y axis title. Character. A label for the y-axis
%                  (default: 'Minimum deletion residual').
%                   Example - 'laby','mmd'
%                   Data Types - char
%
%
%  Output:
%
% brushedUnits  : brushed units. Vector. Vector containing the list of the
%                 units which are inside subset in the trajectories which
%                 have been brushed using option databrush. If option
%                 databrush has not been used brushedUnits will be an empty
%                 value.
%
% See also: mmdrsplot, FSRmdrrs, tclustreg
%
%
% References:
%
%   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
%   Springer Verlag, New York.
%
% Copyright 2008-2020.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mdrrsplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of the use of function mdrrsplot with all the default
    % options.
    % The X data have been introduced by Gordaliza, Garcia-Escudero & Mayo-Iscar (2013).
    % The dataset presents two parallel components without contamination.
    X  = load('X.txt');
    y = X(:,end);
    X =X(:,1:end-1);
    [out]=FSRmdrrs(y,X,'bsbsteps',0,'nsimul',300);
    mdrrsplot(out);
    title('ColorTrj = 1 (default)','FontSize',16,'Interpreter','Latex');
    % now with trajectories with colormap
    mdrrsplot(out,'ColorTrj',0,'tag','ColorTrj_0');
    title('ColorTrj = 0 -- color gradient proportional to sum of mdr','FontSize',16,'Interpreter','Latex');
    % now with 2 trajectories with marker symbols
    mdrrsplot(out,'ColorTrj',2,'tag','ColorTrj_2');
    title('ColorTrj = 2 -- highlight the 2 trj with max mdr','FontSize',16,'Interpreter','Latex');
    % now with trajectories with marker symbols
    mdrrsplot(out,'ColorTrj',3,'tag','ColorTrj_3');
    title('ColorTrj = 3 -- highlight the 3 trj with max mdr','FontSize',16,'Interpreter','Latex');
    cascade;
%}

%{
    % Example of the use of function mdrrsplot with personalized envelopes.
    % tclustreg of contaminated X data using all default options.
    % The X data have been introduced by Gordaliza, Garcia-Escudero & Mayo-Iscar (2013).
    % The dataset presents two parallel components without contamination.
    X  = load('X.txt');
    y = X(:,end);
    X =X(:,1:end-1);
    [out]=FSRmdrrs(y,X,'bsbsteps',0);
    mdrrsplot(out,'quant',[0.99;0.9999]);
%}

%{
    % Example of option datatooltip.
    % Example of the use of function mdrrsplot with datatooltip passed as
    % scalar (that is using default options for datacursor (i.e.
    % DisplayStyle =window)
    load('loyalty.txt','loyalty');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    [out]=FSRmdrrs(y,X);
    mdrrsplot(out,'datatooltip',1);
%}

%{
    % Example of the use option datatooltip passed as structure.
    load('loyalty.txt','loyalty');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    [out]=FSRmdrrs(y,X);
    tooltip=struct;
    tooltip.SnapToDataVertex='on'
    tooltip.DisplayStyle='datatip'
    mdrrsplot(out,'datatooltip',tooltip);
%}


%{
   % Example of the use of option envm.
   % Example of the use of function mdrplot with personalized envelopes.
    load('loyalty.txt','loyalty');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    [out]=FSRmdrrs(y,X);
   %In this case the resuperimposed envelope is based on n-2 observations
   mdrrsplot(out,'envm',length(out.y)-2);

%}

%{
    % Interactive_example
    %Example of the use of function mdrplot with option databrush.
    load('loyalty.txt','loyalty');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    % Note that in this case it is necessary to call function FSRmdrrs in
    % order to save the units belonging to subset for all the steps of the
    % forward search.
    [out]=FSRmdrrs(y,X,'bsbsteps',0);
    brushedUnits=mdrrsplot(out,'databrush',1);
%}

%{
    % Interactive_example
    % Example where databrush is a structure and option labeladd is used
    X  = load('X.txt');
    y = X(:,end);
    X =X(:,1:end-1);
    [out]=FSRmdrrs(y,X,'bsbsteps',0);
    databrush=struct
    databrush.selectionmode='Lasso'
    databrush.labeladd='1';
    mdrrsplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    %Example of the use of option databrush using brush mode
    X  = load('X.txt');
    y = X(:,end);
    X =X(:,1:end-1);
    [out]=FSRmdrrs(y,X);
    databrush=struct
    databrush.selectionmode='Brush'
    mdrrsplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    % Example of the use of persistent non cumulative brush. Every time a
    % brushing action is performed previous highlights are removed but the
    % labels are not removed from the scatterplot matrix.
    % If persist is 'off' after each selection, all trajectories except
    % those selected in the current iteration are plotted in
    % greysh color. If persist is 'on' after each selection, trajectories
    % never selected in any iteration are plotted in greysh color.
    X  = load('X.txt');
    y = X(:,end);
    X =X(:,1:end-1);
    [out]=FSRmdrrs(y,X);
    databrush=struct
    databrush.persist='off'
    databrush.labeladd='1';
    databrush.RemoveLabels='off';
    mdrrsplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    % Example of the use of persistent cumulative brush. Every time a
    % brushing action is performed current highlights are added to
    % previous highlights.
    % If persist is 'off' after each selection, all trajectories except
    % those selected in the current iteration are plotted in
    % greysh color. If persist is 'on' after each selection, trajectories
    % never selected in any iteration are plotted in greysh color.
    X  = load('X.txt');
    y = X(:,end);
    X =X(:,1:end-1);
    [out]=FSRmdrrs(y,X);
    databrush=struct
    databrush.persist='on';
    databrush.selectionmode='Rect'
    mdrrsplot(out,'databrush',databrush)
%}


%% Beginning of code

brushedUnits=[];
mdr=out.mdrrs;
ntrajectories=size(mdr,2)-1;

BBrs=out.BBrs;
nsimul=size(mdr,2)-1;
[n,p]=size(out.X);

% Default limits for x axis
xl1=mdr(1,1)-3;
xl2=mdr(end,1)*1.1;
xlimx=[xl1 xl2];

% Default limits for y axis
yl1=min(min(mdr(:,2:end)));
yl2=max(max(mdr(:,2:end)))*1.1;

thresh=6;
if yl2 >thresh
    mdrrstmp=mdr(:,2:end);
    mdrrstmp(mdrrstmp>thresh)=NaN;
    yl2=max(max(mdrrstmp))*1.1;
else
end
ylimy=[yl1 yl2];

% Default quantiles to compute the envelopes
quant=[0.01;0.5;0.99];

labx='Subset size m';
laby='Minimum deletion residual';


%% User options

options=struct('quant', quant,...
    'envm',n,'xlimx',xlimx,'ylimy',ylimy,'lwd',1.5,'lwdenv',1,...
    'FontSize',12,'SizeAxesNum',10,'tag','pl_mdrrs',...
    'datatooltip','','databrush','',...
    'titl','','labx',labx,'laby',laby,'nameX','','namey','',...
    'ColorTrj',1);

if nargin<1
    error('FSDA:mdrrsplot:missingInputs','A required input argument is missing.')
end

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mdrrsplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

if nargin>1
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

databrush=options.databrush;

if isstruct(databrush)
    fdatabrush=fieldnames(databrush);
    % labeladd option
    d=find(strcmp('labeladd',fdatabrush));
    if d>0
        labeladd=databrush.labeladd;
        % labeladd=labeladd{1};
        % options.databrush(d:d+1)=[];
        databrush=rmfield(databrush,'labeladd');
        fdatabrush=fieldnames(databrush);
    else
        labeladd='';
    end
    
    % bivarfit option
    d=find(strcmp('bivarfit',fdatabrush));
    if d>0
        bivarfit=databrush.bivarfit;
        databrush=rmfield(databrush,'bivarfit');
        fdatabrush=fieldnames(databrush);
    else
        bivarfit='';
    end
    
    % multivarfit option
    d=find(strcmp('multivarfit',fdatabrush));
    if d>0
        multivarfit=databrush.multivarfit;
        databrush=rmfield(databrush,'multivarfit');
        fdatabrush=fieldnames(databrush);
    else
        multivarfit='';
    end
    
    % persist option
    dpers=find(strcmp('persist',fdatabrush));
    if dpers>0
        persist=databrush.persist;
        databrush=rmfield(databrush,'persist');
        
        ColorOrd=[1 0 0;0 1 1; 1 0 1; 1 1 0; 0 0 0; 0 1 0; 0 0 1];
        ColorOrd=repmat(ColorOrd,4,1);
    else
        persist='';
        ColorOrd=[1 0 0];
    end
    
    % FlagColor option Initialize colors for the brushing option: default
    % colors are blue (unbrushed unit) and red (brushed units)
    d=find(strcmp('FlagColor',fdatabrush));
    if d>0
        % flagcol=options.databrush{d+1};
        flagcol=databrush.FlagColor;
        databrush=rmfield(databrush,'FlagColor');
        fdatabrush=fieldnames(databrush);
        
        clr=['b' flagcol 'cmykgbrcmykg'];
    else
        clr='brcmykgbrcmykgbrcmykg';
    end
else
    bivarfit='';
    multivarfit='';
    labeladd='';
    persist='';
    ColorOrd=[1 0 0];
    clr='brcmykgbrcmykgbrcmykg';
end

% Symbol types for y|X plot (in case of brushing)
styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};

% Quantiles associated with the envelope based on all the observations.
quant=options.quant;
nenvel=length(quant);

% envm = number of units to use in order to compute the envelope
envm=options.envm;

% lwd = line width of the curve which contains mdr
lwd=options.lwd;

% lwdenv = line width of the curves associated with the envelopes
lwdenv=options.lwdenv;

% FontSize = font size of the axes labels
FontSize =options.FontSize;

% FontSizeAxes = font size for the axes numbers
SizeAxesNum=options.SizeAxesNum;

% ColorTrj determines the color of the trajectories
% - ColorTrj = 0 for colormap proportional to sum or mdr along the trajectory; 
% - ColorTrj = 1 for rotation of fixed colors;
% - ColorTrj = 2:7 for rotation of fixed colors for the trajectories with 
%              larger mdr (no more than 7 allowed). Markers are also added.
ColorTrjUser = options.ColorTrj;
ColorTrj     = ColorTrjUser;
maxcolors = 7; % this is the number of colors set later in variable ColorSet
if ColorTrj > maxcolors
    ColorTrj = maxcolors;
end

init=mdr(1,1);
if init-p==0
    init=init+1;
end

%% Display the minimum deletion residual plot

% Create a figure to host the plot or clear the existing one
h=findobj('-depth',1,'tag',options.tag);
if (~isempty(h))
    clf(h);
    figure(h)
    axes;
else
    % create a new figure
    h=figure;
end
set(h,'Name', 'Monitoring of Minimum deletion residual', 'NumberTitle', 'off');
hold('all');

% displays the boundary of the current axes.
box('on');

% Compute the theoretical envelopes for minimum deletion residual
[gmin] = FSRenvmdr(envm,p,'prob',quant,'init',init);

% sel = extract the elements of matrix mdr which have to be plotted
sel = 1:size(mdr,1)-n+envm;

% x coordinates where to put the messages about envelopes
xcoord = max([xlimx(1) init]);
for i=1:length(quant)
    % Superimpose chosen envelopes. Note that if internationaltrade==1 the
    % ordering of residuals that determine the FS progression changes, and
    % the envelopes do not apply anymore; therefore, they are not plotted.
    if out.internationaltrade==0
        c05   = 'g';
        c099  =  [0.2 0.8 0.4];
        c0999 =  [0  0 0];
    else
        c05   = 'w';
        c099  = 'w';
        c0999 = 'w';
    end
    if quant(i)==0.5
        % Superimpose 50% envelope
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color',c05,'tag','env');
    elseif quant(i)<=0.99
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color',c099,'tag','env');
    else
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color',c0999,'tag','env');
    end
    % [figx figy] = dsxy2figxy(gca, max([xlimx(1) init]), gmin(1,i+1));
    [figx, figy] = dsxy2figxy(gca, xcoord,gmin(gmin(:,1)==xcoord,i+1));
    kx=0; ky=0;
    
    if isempty(figy) || figy<0
        figy=0;
    else
        if figy>1
            figy=1;
        end
    end
    if isempty(figx) || figx<0
        figx=0;
    else
        if figx>1
            figx=1;
        end
    end
    if out.internationaltrade==0
        annotation(gcf,'textbox',[figx figy kx ky],...
            'String',{[num2str(100*quant(i)) '%']},...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'EdgeColor','none',...
            'BackgroundColor','none',...
            'FitBoxToText','off',...
            'FontSize',FontSize);
    end
end

if out.internationaltrade==0
    % Write an extra message on the plot
    annotation(gcf,'textbox',[0.2 0.8 0.1 0.1],...
        'EdgeColor','none',...
        'String',['Envelope based on ' num2str(envm) ' obs.'],...
        'FontSize',FontSize);
end

% plot minimum deletion residual, but make trajectories invisible
tagstat = 'rs_data_mdr';
plot1   = plot(mdr(sel,1),mdr(sel,2:end),'tag',tagstat,...
               'LineWidth',0.1,'LineStyle','-','Color','w');

% set the x and y axis
xlimx=options.xlimx;
ylimy=options.ylimy;
xlim(xlimx);
ylim(ylimy);
labx=options.labx;
laby=options.laby;
titl=options.titl;
title(titl);

% Add to the plot the labels for values of la
xlabel(labx,'Fontsize',FontSize);
ylabel(laby,'Fontsize',FontSize);

set(gca,'FontSize',SizeAxesNum)

hold('off');

% include specified tag in the current plot
set(gcf,'tag',options.tag);

% %%% Managment of the rainbow %%%

% used to skip end part of FS to avoid typical peaks due to outliers
skipafter    = round(find(mdr(:,1)==floor(n*0.9)));
% used to skip a fixed percentage of the initial part of FS to avoid spurious peaks
skipbefore   = max(init,floor(n*0.1));

% used to skip a the initial extreme spurious peaks outside the bands
%skipbefore   = find(max(mdr(:,2:end),[],2) <= gmin(:,end) , 1);

% this is an adaptive estimateion of the part to skip
% ia = -1; skipbefore = 1; 
% while or(length(ia) > floor(n*0.25) , ia == -1)
%     skipbefore = skipbefore + floor(n*0.05);
%     resmax     = max(mdr(skipbefore:skipafter,2:end),[],1);
%     [~,ia,~]   = unique(resmax);
% end
     
% markers to use each lm steps; by default markers are not used ('none')
lm          = min(10 , floor((n-skipbefore)/10));
MarkerSet   = {'none' ; 'o' ; '+' ; '*'  ; 'x' ; ...
               'square' ; 'diamond' ; '.' };
slinmkr     = repmat(MarkerSet(1),nsimul,1); 

% linestyle; by default, one color per trajectory
LinestyleSet  = {'-';'--';':';'-.';'-';'--';':'};
slinsty       = repmat(LinestyleSet,ceil(nsimul/length(LinestyleSet)),1);
slinsty(nsimul+1:end,:)=[];

% define the selected colors in RGB form
%ColorSet={'b';'g';'r';'c';'m';'y';'k'};
ColorSet = {FSColors.black.RGB;...
            FSColors.blue.RGB;...
            FSColors.red.RGB;...
            FSColors.magenta.RGB;...
            FSColors.green.RGB;...
            FSColors.cyan.RGB;...
            FSColors.yellow.RGB;...
            };

% Line width default
slinwdt  = lwd*ones(nsimul,1);

switch ColorTrj
    case 0  % use a blue colormap
        % colors for the trajectories
        ressum    = sum(mdr(skipbefore:end,(2:end)),1);
        A         = rescaleFS(ressum,1,0);
        [B , iA]  = sort(A,'descend');
        bgcolors  = [zeros(nsimul,1) , B' , ones(nsimul,1)];
        fcol      = num2cell(bgcolors,2);
        
        % colors for the reference bar
        bgcolmap = [zeros(nsimul,1) , linspace(1,0,nsimul)',  ones(nsimul,1) ];
        colormap(bgcolmap);
        c = colorbar;
        c.Label.String = ['Standardised sum of mdr from step ' num2str(skipbefore)];
        c.FontSize = SizeAxesNum;
        caxis([0 1]);
        
        % use only the main line style '-'
        slinsty(2:end) = slinsty(1); 
        
        % line width a bit smaller than the standard size
        slinwdt = slinwdt*0.5;
        
    case 1 % color rotation if ColorTrj > 0, or more in general if ColorTrj ~= 0
        fcol=repmat(ColorSet,ceil(nsimul/size(ColorSet,1)),1);
        fcol(nsimul+1:end,:)=[];
        iA = (1:nsimul)';
        
    otherwise % ColorTrj > 1: use colors and symbols for the main ColorTrj trajectories
        % start with a single line style ':'
        slinsty(1:end) = LinestyleSet(3);
        
        % start with a greish color for all trajectories
        bgcolors  = repmat(FSColors.darkgrey.RGB,nsimul,1);
        fcol      = num2cell(bgcolors,2);
        
        % find the trajectories to highlight
        % resmax    = max(mdr(skip:end,(2:end)),[],1);
        resmax    = max(mdr(skipbefore:skipafter,2:end),[],1);
        [~,ia,ic] = unique(resmax);
        % the next while statement is in case the max is due to a peak in
        % the last part of the search, because of a set of outliers
        while length(ia)==1
            [row,col] = find(mdr(:,2:end)==resmax(1));
            if sum(col) == nsimul*(nsimul+1)/2 %this is a double check
                mdrtmp = mdr;
                mdrtmp(unique(row),:)=[];
                resmax    = max(mdrtmp(skipbefore:skipafter,(2:end)),[],1);
                [~,ia,ic] = unique(resmax);
            end
        end
        
        % keep the best 'maxcolors' indices of the modes
        nmodes    = length(ia);
        ia = flipud(ia);
        if length(ia) > maxcolors
            ia = ia(1:maxcolors);
        end
        ntrj = min(length(ia),ColorTrj);
        
        % set the trajectories
        seq = 1:nsimul;        
        for ii=1:ntrj 
            selii            = seq(ic==nmodes-ii+1);
            slinsty(selii)   = LinestyleSet(ii);            
            fcol(selii,:)    = ColorSet(ii,:);
            slinmkr(selii,:) = MarkerSet(ii+1);
        end
        iA = (1:nsimul)';
        idarkgrey    = find(cellfun(@(x) isequal(x,FSColors.darkgrey.RGB), fcol(:), 'UniformOutput', 1));
        inotdarkgrey = setdiff(iA,idarkgrey);
        iA = [idarkgrey ; inotdarkgrey];
        slinsty = slinsty(iA);
        fcol    = fcol(iA);
        slinmkr = slinmkr(iA);
        slinwdt = [0.5*slinwdt(idarkgrey) ; slinwdt(inotdarkgrey)];
        if ntrj < ColorTrjUser
            disp(['Note that you asked to highlight ColorTrj = ' num2str(ColorTrjUser) ...
                ' groups, but only ' num2str(ntrj) ' are displayed']);
        end
end

% set Color, Linestyle and Marker of the trajectories
set(plot1(iA),{'LineStyle'},slinsty,{'Color'},fcol,{'LineWidth'},num2cell(slinwdt)); 
set(plot1(iA),{'Marker'}   ,slinmkr,'MarkerIndices',1:lm:length(sel));

% Store the handle of the mdrplot inside handle hmin
hmin=gcf;

%% Set the datatooltip for the mdrrsplot
if ~isempty(options.datatooltip)
    try
        chkgpu=gpuDevice; %#ok<NASGU>
        hdt = datacursormode;
        set(hdt,'Enable','on');
        % If options.datatooltip is not a struct then use our default options
        if ~isstruct(options.datatooltip)
            set(hdt,'DisplayStyle','window','SnapToDataVertex','on');
        else
            % options.datatooltip contains a structure where the user can set the
            % properties of the data cursor
            set(hdt,options.datatooltip);
        end
        % Declare a custom datatooltip update function to display additional
        % information about the selected unit
        set(hdt,'UpdateFcn',{@mdrplotLbl,out})
    catch
        disp('No graphical device, interactive datatooltip not enabled')
    end
end

%% Brush mode (call to function selectdataFS)

if ~isempty(options.databrush) || isstruct(options.databrush)
    % Check that BBrs has the proper size
    if size(BBrs,2)~=n-init+1
        warning('FSDA:mdrrsplot:WrongInputOpt','Input 3Darray BBrs inside input structure out \n does not contain information about all subsets for each subset size')
        error('FSDA:mdrrsplot:WrongInputOpt','Please call input function FSRmdrrs with option ''bsbsteps'',0 \n to save the units beloging to subset for all steps of the forward search')
    end
    
    % Preliminary check that BBrs has stored all the subsets for for all
    % the trajectories.
    % BBrs
    
    if isstruct(options.databrush)
        
        
        % sele={options.databrush{:} 'Ignore' findobj(gcf,'tag','env')}; % old code
        % we need to transform the input structure databrush into a cell array
        
        cv=[fieldnames(databrush) struct2cell(databrush)]';
        
        % sele=[options.databrush {'Ignore'} {findobj(gcf,'tag','env')}];
        sele=[cv(:)' 'Ignore' {findobj(gcf,'tag','env')}];
        
        % Add the FlagSize of the brushed points if it has not been previously
        % specified by the user
        d=find(strcmp('FlagSize',fdatabrush));
        if d==0
            sele=[sele 'FlagSize' '3'];
        end
        
    else
        sele={'selectionmode' 'Rect' 'Ignore' findobj(gcf,'tag','env') };
    end
    
    % sele={sele{:} 'Tag' options.tag}; OLD inefficient code
    sele=[sele 'Tag' {options.tag}];
    
    % group = VECTOR WHICH WILL CONTAIN THE IDENTIFIER OF EACH GROUP
    % e.g. group(14)=3 means that unit 14 was selected at the third brush
    group=ones(n,1);
    
    % some local variables
    but=0; brushcum=[]; ij=1; ijk=1;
    
    % Check if X includes the constant term for the intercept.
    X=out.X;
    y=out.y;
    % p=size(X,2);
    
    intcolumn = find(max(X,[],1)-min(X,[],1) == 0);
    
    if intcolumn==1
        p1=1:(p-numel(intcolumn));
        Xsel=X;
        Xsel(:,intcolumn)=[];
    else
        p1=1:p;
        Xsel=X;
    end
    
    % Set the labels of the axes.
    if isempty(options.nameX)
        nameX=cellstr(num2str(p1','X%d'));
    else
        nameX=options.nameX;
    end
    
    if isempty(options.namey)
        namey=char('y');
    else
        namey=options.namey;
    end
    
    
    % add to cell sele option FlagColor (color of selection) and
    % FlagMarker (symbol to be used for selection)
    sele=[sele 'FlagColor' ColorOrd(ij,:) 'FlagMarker' char(styp(ij+1))];
    
    % seleTracumboo = boolean vector which will contain true is a
    % trajectorry has been selected at least once.
    seleTracumboo=false(ntrajectories,1);
    
    % loop brushing
    while but<=1
        
        figure(hmin);
        
        % Remark: function selectdataFS cannot be used on the current figure if
        % the "selection mode" or the "zoom tool" are on. Setting the
        % plotedit mode initially to on and then to off, has the effect to
        % deselect plotedit mode.
        plotedit on
        plotedit off
        
        if strcmp(persist,'off')
            % Remove from the current plot the yellow selection left by
            % selectdataFS, if present.
            a=findobj(gcf,'Tag','selected');
            delete(a);
        elseif strcmp(persist,'on')
            if ij>1
                chkexist=find(strcmp('FlagColor',sele)==1);
                if ~isempty(chkexist)
                    sele{chkexist+1}=ColorOrd(ij,:);
                    sele{chkexist+3}=char(styp(ij+1));
                end
            end
            
        end
        % FROM HERE
        % seqsim = sequence from 1 to the number of random starts
        seqsim=1:nsimul;
        
        % CALL TO FUNCTION selectdataFS
        disp('Select trajectories to brush in the current plot');
        % xs and ys are two cells with size number of random starts-by-1
        % Each element of xs contains the selected steps associated to a
        % particular random starts
        % Each element of ys contains the values of mmd selected associated to a
        % particular random starts
        [pointslist,xs,ys] = selectdataFS(sele{:});
        
        % exit from function if the figure  was closed before selection
        if ~isempty(pointslist) && isnumeric(pointslist) && (min(pointslist) == -999)
            return
        end
        
        % The selected steps
        % The selected x-points corresponding to the selected steps
        % We should have: selxsteps = selsteps-init+1
        selxsteps = unique(cell2mat(pointslist))+init-1;
        
        % nselTraj = total number of trajectories which have been selected
        nselTraj= nsimul-sum(cellfun('isempty', xs));
        disp([num2str(nselTraj) ' trajectories selected' ]);
        
        %% For each brushing operation, do the following:
        if nselTraj>0
            
            % brushj = vector which will contain the units belonging to
            % subset for at least one of the  selected trajectories
            brushj=zeros(n,1);
            
            % Given that the ordering of the trajectories given by selectdataFS
            % may be different from the ordering in which the trajectories are
            % stored inside gca we get the Ydata from gca(hmin) and we find the
            % positions associated with the values of mmd which have been
            % selected
            htraj = findall(gca(hmin),'tag',tagstat);
            a=get(htraj,'Ydata');
            a=cell2mat(a)';
            
            % sumys is a vector with length equal to nsimsul
            % Given that each element of the cell ys may have a different
            % length we use as a summary of each element of ys the sum
            % The elements of sumys which are greater then 0 are the
            % random starts which are selected
            sumys=cellfun(@sum,ys);
            sumys=round(sumys,8);
            
            % UniqueValuesSummarymmd = summary value of unique searches for each
            % selection. The length of UniqueValuesSummarymmd is equal to the number of uniqeu searches selected
            UniqueValuesSummarymmd=unique(sumys);
            UniqueValuesSummarymmd=UniqueValuesSummarymmd(UniqueValuesSummarymmd>0);
            
            % seleUnitjboo is a Boolean matrix of size
            % n-by-length(UniqueValuesSummarymmd) which will contain a
            % true in ith-position if unit i is inside the subset in
            % the min selected stesp of the selected trajectory
            seleUnitjboo=false(n,length(UniqueValuesSummarymmd));
            
            % seleTrajboo is a Boolean vector of size nsimul-by-1
            % which will contain a true in ith-position if trajectory i
            % has been selected
            seleTrajboo=false(nsimul,1);
            
            for j=1:length(UniqueValuesSummarymmd)
                selysteps=[ys{sumys==UniqueValuesSummarymmd(j)}]';
                selxsteps=xs{sumys==UniqueValuesSummarymmd(j)};
                
                % Now find the searches associated with the values of mmd which have
                % been selected
                % The last input element of find is 1 because it is enough to take one
                % search whose value of  mmd is equal to selysteps(1) at step selxsteps(1)
                seleUnitjboo(abs(mdr(selxsteps(1)-init+1,2:end)-selysteps(1))<1e-10,j)=true;
                
                % selTraj contains the reverse indexes of the selected trajectories
                % using the ordering in which they are stored inside the figure
                % For examples if selTraj =[ 6 14 ] and there are 15
                % trajectrories 6 means 10 and 14 means 2. Note that a is
                % not out.mmdrs but changes in every iteration becase
                % selected trajectories are put on top.
                selTraj=find(abs(a(selxsteps(1)-init+1,:)-selysteps(1))<1e-10);
                
                % disp('Indexes of the selected trajectories as they are stored inside the plot')
                % disp(selTraj)
                seleTrajboo(selTraj)=true;
                
                % seleTracumboo = boolean vector which contains true in
                % correspondence of the selected trajectories for all
                % brushing
                seleTracumboo(selTraj)=true;
                su=sum(seleTracumboo);
                
                % Now find the units belonging to subset at step selxsteps(1) for
                % random start number searchtoExtr
                searchtoExtr=find(seleUnitjboo(:,j),1);
                BB_selTraj = BBrs(:,:,searchtoExtr);
                bsbstepjboo = BB_selTraj(:,selxsteps(1)-init+1);
                
                % bsbjboo = units forming subset in selected step
                bsbjboo=~isnan(bsbstepjboo);
                brushj(bsbstepjboo(bsbjboo))=bsbstepjboo(bsbjboo);
            end
            
            % nbrush= complete set of brushed units in the last selection
            nbrush=brushj(brushj>0);
            
            if strcmp(persist,'on')
                % Set new colors and new line style for the
                % trajectories never selected in any iteration
                % Set new colors and new line style for the unselected trajectories
                set(htraj(~seleTracumboo),'Color',FSColors.greysh.RGB);
                set(htraj(~seleTracumboo),'LineStyle',':');
            else
                % Set new colors and new line style for the
                % trajectories not selcted in the current iteration
                set(htraj(~seleTrajboo),'Color',FSColors.greysh.RGB);
                set(htraj(~seleTrajboo),'LineStyle',':');
            end
            
            % seleTracumboo contains the trajectories which have been
            % selected at least once. Note that selected trajectories are
            % always put on the first columns of matrix a where
            % a is defined as
            % htraj = findall(gca(hmin),'tag',tagstat);
            % a=get(htraj,'Ydata');
            % a=cell2mat(a)';
            seleTracumboo=false(ntrajectories,1);
            seleTracumboo(1:su)=true;
            
            % Set the color of selected trajectories
            set(htraj(seleTrajboo),'Color',clr(ij+1));
            
            % selTraj = indices of the selected trajectories
            selTraj=seqsim(seleTrajboo);
            
            %  unselTraj   = indices of the unselected trajectories
            unselTraj   = seqsim(~seleTrajboo);
            
            % Now reorder the lines inside the plot in such a way that those which are selected
            % appear on top of the others
            chH=get(gca,'Children');
            % select the lines which have the tag rs_data_mmd (that is the
            % trajectories associated with the random starts)
            lineindexes=strcmp(get(chH,'Tag'),'rs_data_mdr');
            
            % inlch = number of children which do not have the tag rs_data_mmd, (i.e.
            % length of the legends and of the lines associated with the envelopes)
            inlch=sum(~lineindexes);
            
            % The if which follows is necessary because after the first
            % selection the lines associated with the envelopes have
            % already been put at the beginning. Note that it is necessary
            % to use index ijk instead of ij because it may happen that the
            % user does selections which do not contain anything
            if ijk==1
                set(gca,'Children',[chH(~lineindexes);chH(selTraj+inlch-nenvel);chH(unselTraj+inlch-nenvel)])
                ijk=ijk+1;
            else
                set(gca,'Children',[chH(~lineindexes);chH(selTraj+inlch);chH(unselTraj+inlch)])
            end
            
            % brushcum =
            % - the list of selected observations in all iterations if
            %   persist=on
            % - the list of selected observations in the current iteration
            %   if persist=off
            if strcmp(persist,'on')
                brushcum=unique([brushcum; nbrush(:)]);
            else
                brushcum=nbrush;
                group=ones(n,1);
            end
            
            % group=vector of length(Xsel) observations taking values
            % from 1 to the number of groups selected.
            % unigroup= list of selected groups.
            group(nbrush)=ij+1;
            unigroup=unique(group);
            
            %% - display the yXplot with the corresponding groups of units highlighted
            
            h=findobj('-depth',1,'Tag','pl_yX');
            
            if (~isempty(h))
                % delete from the current figure all graphics objects whose
                % handles are not hidden
                clf(h);
                % Make the figure identified by handle h become the current
                % figure make it visible, and raise it above all other
                % figures on the screen.
                figure(h);
            else
                % create a new figure and set its style equal to that of
                % the mdrplot.
                figure;
                set(gcf,'WindowStyle',get(hmin,'WindowStyle'));
            end
            
            [H,AX,BigAx] = gplotmatrix(Xsel,y,group,clr(unigroup),char(styp{unigroup}),[],'on',[],nameX,namey);
            
            % Assign to this figure a name and a tag=pl_yX
            set(gcf,'Name','Scatter plot matrix y|X with selected groups highlighted');
            set(gcf,'tag','pl_yX');
            
            % Set markers
            for mfc=1:length(unigroup)
                set(findobj(gcf,'marker',char(styp(unigroup(mfc)))),'MarkerFaceColor',clr(unigroup(mfc)));
            end
            
            % Set the legend properties of the gplotmatrix
            set(H(:,:,1),'DisplayName','Unbrushed units');
            for brugrp = 2:length(unigroup)
                set(H(:,:,brugrp),'DisplayName',['Brushed units ' num2str(brugrp-1)]);
            end
            
            % save the indices of the last selected units (nbrush) to the
            % 'UserData' field of the last selected group of H(:,:,end)
            set(H(:,:,end), 'UserData' , nbrush);
            
            % The following line adds objects to the panels of the yX
            % add2yX(H,AX,BigAx,out,group,nbrush,bivarfit,multivarfit,labeladd)
            intercept=1;
            add2yX(H,AX,BigAx,'intercept',intercept,'bivarfit',bivarfit,'multivarfit',multivarfit,'labeladd',labeladd);
            
            disp('Brushed steps');
            disp(selxsteps);
            
            disp('Associated brushed units');
            disp(nbrush);
            
            %             disp('Steps of entry of brushed units');
            %             disp(Un(selindex,:));
            
        end
        
        
        %% check condition to exit from the brush mode
        % If the option persistent is not equal off or on than get out of
        % the loop
        if strcmp('on',persist) || strcmp('off',persist)
            
            % variable ij is linked to the highlighting color
            if strcmp('on',persist)
                ij=ij+1;
                % Set but=1 so that previous highlightments in other
                % figures are not deleted
                but=1;
            end
            
            
            % Before waitforbuttonpress:
            % - the mdrplot is highlighted again
            figure(hmin);
            
            % Lay down the plots before continuing
            position(hmin);
            
            % - and a function to be executed on figure close is set
            set(gcf,'CloseRequestFcn',@closereqFS);
            
            disp('Press a mouse key to continue brushing, a keyboard key to stop')
            ss=waitforbuttonpressFS;
            disp('------------------------')
            
            %TODO:mdrrsplot:closereq
            set(gcf,'CloseRequestFcn','closereq');
            Open_yX = findobj(0, 'type', 'figure','tag','pl_yX');
            Open_res = findobj(0, 'type', 'figure','tag','pl_resfwd');
            Open_mdr = findobj(0, 'type', 'figure','tag','pl_mdrrs');
            if isempty(Open_mdr)  % User closed the main brushing window
                if ~isempty(Open_yX); delete(Open_yX); end    % yX plot is deleted
                if ~isempty(Open_res); delete(Open_res); end  % monitoring residual plot is deleted
                delete(get(0,'CurrentFigure')); % deletes Figure if still one left open
            end
            
            % - and the 'but' variable is set if keyboard key was pressed
            if ss==1
                but=2;
            end
        else
            but=2;
        end
        
    end % close loop associated with but
    brushedUnits=brushcum;
end % close options.databrush

    function output_txt = mdrplotLbl(~,event_obj,out)
        %Provides information about the selected point in the trajectory of the minimum deletion residual
        %(once selected with the mouse)
        %
        % Required input arguments:
        %
        %       obj =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the sense that
        %               these arguments are automatically passed to the function when it executes.
        %     out   =   a structure containing the following fields
        %       y   =   the response of the regressione model
        %      Un   =   a matrix containing the list of the units which entered the subset
        %               in each step of the search
        %   label   =  (optional argument) if it is present it must be
        %               a cell array of strings containing the labels of
        %               the rows of the regression dataset
        %
        % Output:
        %
        %   output_txt=  Datatip text (string or string cell array) which informs
        %                about the value of the mdr at the selected step of the
        %                search, the selected step of the search and the unit which
        %                will enter the search at the next step.
        %
        % REMARK: this function is called by function mdrplot
        %
        % References:
        %
        %   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
        %   Springer Verlag, New York.
        %
        % Written by FSDA team
        
        %%
        pos = get(event_obj,'Position');
        
        % x and y, plot coordinates of the mouse
        x = pos(1); y = pos(2);
        
        % output_txt is what it is shown on the screen
        output_txt = {['mdr=',num2str(y,4)]};
        
        % Add information abou the step of the search which is under investigation
        output_txt{end+1} = ['Step m=' num2str(x)];
        
        % If structure out does not contain labels for the rows then
        % labels row1....rown are added automatically
        if isempty(intersect('label',fieldnames(out)))
            out.label=cellstr(num2str((1:length(out.y))','row%d'));
        end
        
        
        % Number of searches associated to a particular value selected
        seltraj=abs(mdr(x-init+1,2:end)-y(1))<1e-10;
        
        % Add information about the number of trajectories selected
        % search
        output_txt{end+1} = ['Number of trajectories='  num2str(sum(seltraj))];
        
    end

end
%FScategory:CLUS-RobClaREG

