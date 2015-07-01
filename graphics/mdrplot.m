function mdrplot(out,varargin)
%mdrplot plots the trajectory of minimum deletion residual (mdr)
%
%<a href="matlab: docsearchFS('mdrplot')">Link to the help function</a>
%
% Required input arguments:
%
%  out :  Structure containing monitoring of mdr. Structure. 
%               Structure containing the following fields.
%    out.mdr =  minimum deletion residual. A matrix containing the monitoring of minimum deletion
%               residual in each step of the forward search. The first
%               column of mdr must contain the fwd search index
%               This matrix can be created using function FSReda
%               (compulsory argument)
%       out.Un  =   order of entry of each unit. Matrix containing the order of entry of each unit
%               (necessary if datatooltip is true or databrush is not empty)
%       out.y   =   response. Vector containing the response (necessary only if
%               option databrush is not empty)
%       out.X   =   Regressors. A matrix containing the explanatory variables
%               (necessary only if option databrush is not empty)
%     out.Bols  =   Monitoring of beta coefficients. (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward search
%               (necessary only if option databrush is not empty and
%               suboption multivarfit is not empty)
%
%
% Optional input arguments:
%
%       quant   :   Quantiles for which envelopes have
%                   to be computed. Vector. The default is to produce 1%, 50% and
%                   99% envelopes. In other words the default is
%                   quant=[0.01;0.5;0.99];
%                   Example - 'quant',[0.05;0.5;0.95]
%                   Data Types - double
%       exact:      Exact of approximate cdf for envelope calculation. Scalar. If it is equal to 1 the calculation of the
%                   quantiles of the T and F distribution is based on
%                   functions finv and tinv from the Matlab statistics
%                   toolbox, otherwise the calculations of the former
%                   quantiles is based on functions invcdff and invcdft.
%                   The solution has a tolerance of 1e-8 (change variable
%                   tol in files invcdff.m and invcdft.m if required.
%                   Example - 'exact',1
%                   Data Types - double
%                   Remark: the use of functions tinv and finv is more
%                   precise but requires more time. The default value of
%                   exact is 0 (approximate solution).
%       sign    :   mdr with sign. Scalar. If it is equal 1 (default) we distinguish steps
%                   for which minimum deletion residual was associated with
%                   positive or negative value of the residual. Steps
%                   associated with positive values of mdr are plotted in
%                   black, while other steps are plotted in red
%                   Example - 'sign',1
%                   Data Types - double
%       mplus1  :   plot of (m+1)th order statistic. Scalar. Scalar, which specifies if it is necessary to plot the
%                   curve associated with (m+1)th order statistic
%                   Example - 'mplus1',1
%                   Data Types - double
%       envm    :   sample size for drawing enevlopes. Scalar. Scalar which specifies the size of the sample which is
%                   used to superimpose the envelope. The default is to add
%                   an envelope based on all the observations (size n
%                   envelope).
%                   Example - 'envm',100
%                   Data Types - double
%       xlimx   :   min and max for x axis. Vector. vector with two elements controlling minimum and
%                   maximum on the x axis. Default value is mdr(1,1)-3 and
%                   mdr(end,1)*1.3
%                   Example - 'xlimx',[20 100]
%                   Data Types - double
%       ylimy   :   min and max for x axis. Vector. Vector with two elements controlling minimum and
%                   maximum on the y axis. Default value is min(mdr(:,2))
%                   and max(mdr(:,2));
%                   Example - 'ylimy',[2 6]
%                   Data Types - double
%       lwdenv  :   Line width. Scalar. Scalar which controls the width of the lines associated
%                   with the envelopes. Default is lwdenv=1
%                   Example - 'lwdenv',2
%                   Data Types - double
%       tag     :   plot handle. String. String which identifies the handle of the plot which
%                   is about to be created. The default is to use tag
%                   'pl_mdr'. Notice that if the program finds a plot which
%                   has a tag equal to the one specified by the user, then
%                   the output of the new plot overwrites the existing one
%                   in the same window else a new window is created
%                   Example - 'tag','mymdr'
%                   Data Types - char 
%   datatooltip :   interactive clicking. Empty value (default) or
%                   structure. The default is datatooltip=''.
%                   If datatooltip is not empty the user can use the mouse
%                   in order to have information about the unit selected,
%                   the step in which the unit enters the search and the
%                   associated label. If datatooltip is a structure, it is
%                   possible to control the aspect of the data cursor (see
%                   function datacursormode for more details or the
%                   examples below). The default options of the structure
%                   are DisplayStyle='Window' and SnapToDataVertex='on'
%                   Example - 'datatooltip',''
%                   Data Types - char 
%       label   :   row labels. Cell of strings. Cell containing the labels
%                   of the units (optional argument used when
%                   datatooltip=1. If this field is not present labels
%                   row1, ..., rown will be automatically created and
%                   included in the pop up datatooltip window)
%                   Example - 'label',{'Smith','Johnson','Robert','Stallone'}
%                   Data Types - cell 
%    databrush :    interactive mouse brushing. Empty value (default),
%                   scalar or structure.
%                   DATABRUSH IS AN EMPTY VALUE 
%                   If databrush is an empty
%                   value (default), no brushing is done. The activation of
%                   this option (databrush is a scalar or a structure) enables
%                   the user  to select a set of trajectories in the
%                   current plot and to see them highlighted in the y|X
%                   plot (notice that if the plot y|X does not exist it is
%                   automatically created). In addition, brushed units can
%                   be highlighted in the monitoring residual plot
%                   Remark. the window style of the
%                   other figures is set equal to that which contains the
%                   monitoring residual plot. In other words, if the
%                   monitoring residual plot is docked all the other
%                   figures will be docked too.
%                  DATABRUSH IS A SCALAR
%                   If databrush is a scalar the default selection tool is a
%                   rectangular brush and it is possible to brush only once
%                   (that is persist='').
%                  DATABRUSH IS A STRUCTURE If databrush is a structure, it is
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
%                  bivarfit. This option is to add one or more least
%                     square lines to the plots of y|X, depending on the
%                     selected groups. bivarfit = ''
%                       is the default: no line is fitted.
%                     bivarfit = '1'
%                       fits a single ols line to all points of each
%                       bivariate plot in the scatter matrix y|X.
%                     bivarfit = '2'
%                       fits two ols lines: one to all points and another
%                       to the last selected group. This is useful when
%                       there are only two groups, of which one refers to a
%                       set of potential outliers.
%                     bivarfit = '0'
%                       fits one ols line for each selected group. This is
%                       useful for the purpose of fitting mixtures of
%                       regression lines.
%                     bivarfit = 'i1' or 'i2' or 'i3' etc.
%                       fits a ols line to a specific group, the one with
%                       index 'i' equal to 1, 2, 3 etc.
%                   multivarfit. If this option is '1', we add to each scatter
%                     plot of y|X a line based on the fitted hyperplane
%                     coefficients. The line added to the scatter plot y|Xi
%                     is mean(y)+Ci*Xi, being Ci the coefficient of Xi. The
%                     default value of multivarfit is '', i.e. no line is
%                     added.
%                   labeladd. If this option is '1', we label the units
%                     of the last selected group with the unit row index in
%                     matrices X and y. The default value is labeladd='',
%                     i.e. no label is added.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct 
%                   Remark: if databrush is a cell, it is possible to
%                   specify all optional arguments of function selectdataFS
%                   and LXS inside the curly brackets of option databrush.
%       Fontsize:   Size of axes labels. Scalar. Scalar which controls the
%                   fontsize of the labels of the axes. Default value is 12
%                   Example - 'Fontsize',14
%                   Data Types - single | double
%    SizeAxesNum:   Size of axes numbers. Scalar which controls the fontsize of the numbers of
%                   the axes. Default value is 10
%                   Example - 'SizeAxesNum',14
%                   Data Types - single | double
%       nameX   :   Regressors names. Cell array of strings. Cell array of strings of length p containing the labels
%                   of the varibles of the regression dataset. If it is empty
%                 	(default) the sequence X1, ..., Xp will be created
%                   automatically
%                   Example - 'nameX',{'Age','Income','Married','Profession'}
%                   Data Types - cell 
%       namey   :   response label. Character. Character containing the label of the response
%                   Example - 'namey','response label'
%                   Data Types - char 
%       lwd     :   Curves line width. Scalar. Scalar which controls linewidth of the curve which
%                   contains the monitoring of minimum deletion residual.
%                   Default line width=2
%                   Example - 'lwd',3
%                   Data Types - single | double 
%       titl    :   main title. Character. A label for the title (default: '')
%                   Example - 'namey','Plot title'
%                   Data Types - char 
%       labx    :   x axis title. Character. A label for the x-axis (default: 'Subset size m')
%                   Example - 'labx','Subset size m'
%                   Data Types - char 
%       laby    :   y axis title. Character. A label for the y-axis (default: 'Minimum deletion residual')
%                   Example - 'laby','mdr'
%                   Data Types - char 
%
%
%
%
% Output: 
%
% See also: resfwdplot
%
% References:
%
%   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
%   Springer Verlag, New York.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mdrplot')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    %Steps common to all the examples
    load('loyalty.txt','loyalty');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
%}

%{
    [outLXS]=LXS(y,X,'nsamp',1000);
    [out]=FSReda(y,X,outLXS.bs);
    % Example of the use of function mdrplot with all the default options
    mdrplot(out);
%}

%{
    %Example of the use of function mdrplot with personalized envelopes
    mdrplot(out,'quant',[0.99;0.9999]);
%}

%{
    %Example of the use of function mdrplot with datatooltip passed as
    %scalar (that is using default options for datacursor (i.e.
    %DisplayStyle =window)
     mdrplot(out,'datatooltip',1);
%}

%{
    %Example of the use of function mdrplot with datatooltip passed as
    %structure

    clear tooltip
    tooltip.SnapToDataVertex='on'
    tooltip.DisplayStyle='datatip'
    mdrplot(out,'datatooltip',tooltip);
%}

%{
    %Example of the information which can be extracted from option sign=1
    %(default). If the data come from a distribution which has positive
    %asymmetry generally the last part of the search is associated with
    %positive values of the minimum residual among the units belonging to
    %subset The example below shows that the last 60 steps of the curve are
    %associated with a black curve (positive value of mdr)

    state = 137; state1=4567;
    rand('state', state);
    randn('state', state1);
    X=randn(200,3);
    y=chi2rnd(8,200,1);
    [outLXS]=LXS(y,X,'nsamp',1000);
    [out]=FSReda(y,X,outLXS.bs);
    mdrplot(out,'sign',1);
%}

%{
   %Example of the use of option envm
   %In this case the resuperimposed envelope is based on n-2 observations
   mdrplot(out,'envm',length(out.y)-2);

%}

%{
    % Interactive_example
    %Example of the use of function mdrplot with databrush
     mdrplot(out,'databrush',1);
%}

%{
    % Interactive_example
    %Example where databrush is a structure
    databrush=struct
    databrush.selectionmode='Lasso'
     mdrplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    %Example of the use of brush using brush mode
    databrush=struct
    databrush.selectionmode='Brush'
    databrush.Label='on';
    mdrplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    %Example of the use of persistent non cumualtive brush. Every time a
    %brushing action is performed previous highlightments are removed
    databrush=struct
    databrush.persist='off'
    databrush.RemoveLabels='off';
    mdrplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    %Example of the use of persistent cumulative brush. Every time a
    %brushing action is performed current highlightments are added to
    %previous highlightments
    databrush=struct
    databrush.persist='on';
    databrush.selectionmode='Rect'
    mdrplot(out,'databrush',databrush)
%}


%% Initialization

% Extract the absolute value of minimum deletion residual
mdr=abs(out.mdr);

[n,p]=size(out.X);

% seq= column vector containing the sequence 1 to n
% seq= (1:n)';

% numtext= a cell of strings used to labels the units with their position
% in the dataset.
% numtext=cellstr(num2str(seq,'%d'));

% Default limits for x axis
xl1=mdr(1,1)-3;
xl2=mdr(end,1)*1.1;
xlimx=[xl1 xl2];

% Default limits for y axis
yl1=min(mdr(:,2));
yl2=max(mdr(:,2))*1.1;
ylimy=[yl1 yl2];

% Default quantiles to compute the envelopes
quant=[0.01;0.5;0.99];

labx='Subset size m';
laby='Minimum deletion residual';


%% User options

options=struct('quant', quant,'exact',0,'sign',0,'mplus1',0,...
    'envm',n,'xlimx',xlimx,'ylimy',ylimy,'lwd',2,'lwdenv',1,...
    'FontSize',12,'SizeAxesNum',10,'tag','pl_mdr',...
    'datatooltip','','databrush','',...
    'titl','','labx',labx,'laby',laby,'nameX','','namey','','label','');

if nargin<1
    error('FSDA:mdrplot:missingInputs','A required input argument is missing.')
end

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mdrplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

if nargin>1
    for i=1:2:length(varargin);
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
        %         bivarfit=options.databrush(d+1);
        %         bivarfit=bivarfit{1};
        %         options.databrush(d:d+1)=[];
        bivarfit=databrush.bivarfit;
        databrush=rmfield(databrush,'bivarfit');
        fdatabrush=fieldnames(databrush);
    else
        bivarfit='';
    end
    
    % multivarfit option
    d=find(strcmp('multivarfit',fdatabrush));
    if d>0
        %         multivarfit=options.databrush(d+1);
        %         multivarfit=multivarfit{1};
        %         options.databrush(d:d+1)=[];
        multivarfit=databrush.multivarfit;
        databrush=rmfield(databrush,'multivarfit');
        fdatabrush=fieldnames(databrush);
    else
        multivarfit='';
    end
    
    % persist option
    dpers=find(strcmp('persist',fdatabrush));
    if dpers>0
        %         persist=options.databrush(d+1);
        %         options.databrush(d:d+1)=[];
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
        %if dpers>0
        %    ColorOrd=
    else
        clr='brcmykgbrcmykgbrcmykg';
        flagcol='r';
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
exact=options.exact;
mplus1=options.mplus1;
sign=options.sign;

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

init=mdr(1,1);
if init-p==0;
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

% Theoretical envelopes for minimum deletion residual
[gmin] = FSRenvmdr(envm,p,'prob',quant,'init',init,'exact',exact);

box('on');

% sel= extract the elements of matrix mdr which have to be plotted
sel=1:size(mdr,1)-n+envm;

% set the x and y axis
xlimx=options.xlimx;
ylimy=options.ylimy;
xlim(xlimx);
ylim(ylimy);

% x coordinates where to put the messages about envelopes
xcoord=max([xlimx(1) init]);
for i=1:length(quant)
    % Superimpose chosen envelopes
    if quant(i)==0.5;
        % Superimpose 50% envelope
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
    elseif quant(i)<=0.99
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
    else
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color',[0  0 0],'tag','env');
    end
    
    % [figx figy] = dsxy2figxy(gca, max([xlimx(1) init]), gmin(1,i+1));
    [figx, figy] = dsxy2figxy(gca, xcoord,gmin(gmin(:,1)==xcoord,i+1));
    kx=0; ky=0;
    
    if isempty(figy) || figy<0;
        figy=0;
    else
        if figy>1;
            figy=1;
        end
    end
    if isempty(figx) || figx<0;
        figx=0;
    else
        if figx>1;
            figx=1;
        end
    end
    
    
    
    
    annotation(gcf,'textbox',[figx figy kx ky],'String',{[num2str(100*quant(i)) '%']},...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','off',...
        'FontSize',FontSize);
end

% plot minimum deletion residual
stdColor = 'b'; stdLineStyle = '-';
plot(mdr(sel,1),mdr(sel,2),'tag','data_mdr',...
    'Color',stdColor,'LineStyle',stdLineStyle,'LineWidth',lwd);



% Write an extra message on the plot
annotation(gcf,'textbox',[0.2 0.8 0.1 0.1],'EdgeColor','none','String',['Envelope based on ' num2str(envm) ' obs.'],'FontSize',FontSize);

% If mplus1=1 add the line associated with (m+1) ordered deletion residual
if mplus1;
    line(mdr(sel,1),mdr(sel,3),'LineWidth',lwd,'LineStyle','.','Color','cyan','tag','env');
end

labx=options.labx;
laby=options.laby;
titl=options.titl;
title(titl);
% Add to the plot the labels for values of la
xlabel(labx,'Fontsize',FontSize);
ylabel(laby,'Fontsize',FontSize);

set(gca,'FontSize',SizeAxesNum)

set(gca,'FontSize',SizeAxesNum)

if sign==1;
    mdrori=out.mdr(sel,1:3);
    mdrori(:,3)=NaN;
    % select steps with positive values of mdr or steps where a negative
    % value is followed by a positive value of mdr
    selstep=(mdrori(:,2)>0) | [(mdrori(1:end-1,2)<0 & mdrori(2:end,2)>0);0];
    mdrori(selstep>0,3)=mdr(selstep>0,2);
    
    % Find steps where it is necessary to include NaN These steps are those
    % in which a positive value of mdr is followed by a negative value
    stepna=mdrori(mdrori(1:end-1,2).*mdrori(2:end,2)<0 & mdrori(2:end,2)<0,1)+1;
    % add a column of NaN
    stepna=cat(2,stepna,NaN*zeros(length(stepna),1));
    
    % Vertically concatenate matrix mdrori and matrix stepna
    mdrt=cat(1,mdrori(:,[1 3]),stepna);
    
    % Sort all rows of matrix mdrt and make sure that in correspondence of
    % the duplicated steps NaN is in the first value of the two duplicated
    % steps
    mdrsor=sortrows(mdrt,-1);
    mdrsor=mdrsor(length(mdrsor):-1:1,:);
    
    % Add with blue color the line associated with positive values of mdr
    line(mdrsor(:,1),mdrsor(:,2),'LineWidth',2,'Color','black','tag','env');
end

hold('off');

% displays the boundary of the current axes.
box on

% include specified tag in the current plot
set(gcf,'tag',options.tag);

% Store the handle of the mdrplot inside handle hmin
hmin=gcf;

%% Set the datatooltip for the mdrplot
if ~isempty(options.datatooltip)
    % datacursormode on;
    hdt = datacursormode;
    set(hdt,'Enable','on'); % DDD
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
end

%% Brush mode (call to function selectdataFS)
% if (~isempty(options.databrush) || iscell(options.databrush))
%
%     if isscalar(options.databrush)
%         sele={'selectionmode' 'Rect' 'Ignore' findobj(gcf,'tag','env') };
%     else
%         sele=[options.databrush {'Ignore' findobj(gcf,'tag','env')}];
%     end
%

if ~isempty(options.databrush) || isstruct(options.databrush)
    
    
    if isstruct(options.databrush)
        
        % If option Label is 'on' then matrix Un is added to UserData
        d=max(strcmp('Label',fieldnames(databrush)));
        if d==1 && strcmp(databrush.Label,'on')
            set(gcf,'UserData',out.Un)
        end
        
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
    but=0; brushcum=[]; ij=1;
    
    % Check if X includes the constant term for the intercept.
    X=out.X;
    y=out.y;
    % p=size(X,2);
    
    intcolumn = find(max(X,[],1)-min(X,[],1) == 0);
    
    if intcolumn==1;
        intercept=1;
        p1=1:(p-numel(intcolumn));
        Xsel=X;
        Xsel(:,intcolumn)=[];
    else
        intercept=0;
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
    
    %Before brushing, check if there is a resplot open and save the colors,
    % the type and the width of its lines.
    hplres = findobj('-depth',1,'Tag','pl_resfwd');
    if ~isempty(hplres)
        hLplres = findobj(hplres,'Type','line');
        lineresC = get(hLplres,'Color');
        lineresLwd = get(hLplres(1),'LineWidth');
    end
    
    % loop brushing
    while but<=1
        figure(hmin);
        
        % Remark: function selectdataFS cannot be used on the current figure if
        % the "selection mode" or the "zoom tool" are on. Setting the
        % plotedit mode initially to on and then to off, has the effect to
        % deselect plotedit mode.
        plotedit on
        plotedit off
        
        if strcmp(persist,'off');
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
        
        % CALL TO FUNCTION selectdataFS
        disp('Select steps to brush in the current plot');
        % [unused,xs] = selectdataFS(sele{:}, 'Label', 'off');
        [unused,xs] = selectdataFS(sele{:});
        
        % exit from function if the mdrplot was closed before selection
        if ~isempty(unused) && isnumeric(unused) && (min(unused) == -999)
            return
        end
        
        % selsteps= the list of the steps which have been brushed
        selsteps=xs+1;
        
        % selindex= the indexes of selsteps in matrix un
        lselsteps=length(selsteps);
        selindex=zeros(lselsteps,1);
        Un=out.Un;
        for i=1:lselsteps
            selindex(i)=find(Un(:,1)==selsteps(i));
        end
        
        % Unsel= units which entered the search in the brushed steps
        Unsel=Un(selindex,2:end);
        
        % nbrush = vector which contains the list of the selected units
        nbrush=Unsel(~isnan(Unsel));
        
        %% For each brushing operation, do the following:
        if ~isempty(nbrush)
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
            
            [H,~,~] = gplotmatrix(Xsel,y,group,clr(unigroup),char(styp{unigroup}),[],'on',[],nameX,namey);
            
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
            add2yX('intercept',intercept,'bivarfit',bivarfit,'multivarfit',multivarfit,'labeladd',labeladd);
            
            
            %% - highlight brushed trajectories also in the resplot, if it is open
            
            % Now check if the figure which monitors the residuals is open.
            % If it is, then also in that figure highlight the trajectories
            % of the brushed units
            h=findobj('-depth',1,'Tag','pl_resfwd');
            if (~isempty(h))
                
                % make figure which contains monitoring scaled residuals
                % become the current figure
                figure(h);
                
                % Condition || but==0 if but=0 then it is necessary to
                % remove previous highlightments (even if persist='on')
                if strcmp(persist,'off') || but==0;
                    lines = findobj(h,'Type','line');
                    set(lines,{'Color'},lineresC,'LineWidth',lineresLwd);
                    
                    %                     % If set of values has already been highlighted in the
                    %                     % mdr plot, remove it
                    %                     a=findobj(h,'Tag','brush_res');
                    %                     delete(a);
                    %
                    %                     % Remove the yellow selection in this plot if present
                    %                     a=findobj(h,'Tag','selected');
                    %                     delete(a);
                end
                
                if strcmp('on',persist)
                    set(lines((n+1-nbrush)),'Color',ColorOrd(ij,:),'LineWidth',3);
                else
                    set(lines((n+1-nbrush)),'Color',flagcol,'LineWidth',3);
                end
                
                %                 % get the x and y coordinates of the monitoring of the
                %                 % scaled residuals
                %                 a=findobj(h,'tag','data_res'); %a=lines;
                %                 abrush=get(a(n+1-nbrush),'ydata');
                %                 if iscell(abrush)
                %                     ycoord=cell2mat(abrush); % y coordinates of scaled residuals (values)
                %                 else
                %                     ycoord=abrush; % y coordinate of scaled residual (value)
                %                 end
                %                 xcoord=get(a(1),'Xdata'); % x coordinates of mdr (steps)
                %
                %                 hold('on');
                %                 if strcmp('on',persist)
                %                     plot(gca,xcoord,ycoord,'LineWidth',4,'color',ColorOrd(ij,:),'tag','brush_res');
                %                 else
                %                     plot(gca,xcoord,ycoord,'LineWidth',4,'color',flagcol,'tag','brush_res');
                %                 end
                %                 hold('off');
                
            end
            
            disp('Brushed steps');
            disp(selsteps);
            
            disp('Associated brushed units');
            disp(nbrush);
            
            disp('Steps of entry of brushed units');
            disp(Un(selindex,:));
            
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
            %             % After waitforbuttonpress:
            %             % - the standard MATLAB function to be executed on figure
            %             %   close is recovered
            %             set(gcf,'CloseRequestFcn','closereq');
            
            % After waitforbuttonpress:
            % - the standard MATLAB function to be executed on figure
            %   close is recovered
            set(gcf,'CloseRequestFcn','closereq');
            Open_yX = findobj(0, 'type', 'figure','tag','pl_yX');
            Open_res = findobj(0, 'type', 'figure','tag','pl_resfwd');
            Open_mdr = findobj(0, 'type', 'figure','tag','pl_mdr');
            if isempty(Open_mdr)  % User closed the main brushing window
                if ~isempty(Open_yX); delete(Open_yX); end    % yX plot is deleted
                if ~isempty(Open_res); delete(Open_res); end  % monitoring residual plot is deleted
                delete(get(0,'CurrentFigure')); % deletes Figure if still one left open
            end
            
            
            
            % - and the 'but' variable is set if keyboard key was pressed
            if ss==1;
                but=2;
            end
        else
            but=2;
        end
        
    end % close loop associated with but
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
        
        % Add information about the corresponding row label of what has been
        % selected
        %output_txt{end+1} = ['Unit: ' cell2mat(out.label(row))];
        
        % Add information about the next step in which the selected unit entered the
        % search
        Un=out.Un;
        idx = find(Un(:,1)==x,1)+1;
        sel=Un(idx,2:end);
        
        output_txt{end+1} = ['Unit(s) entered in step ' num2str(x+1) '='  num2str(sel(~isnan(sel)))];
    end

end

