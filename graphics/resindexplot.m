function resindexplot(residuals,varargin)
%resindexplot plots the residuals from a regression analysis versus index number or any other variable
%
%
%<a href="matlab: docsearchFS('resindexplot')">Link to the help function</a>
%
% Required input arguments:
%
%  residuals :  residuals to plot. Numeric vector or structure. If
%               residuals is a vector it contains the n residuals.
%               If residuals is a structure it contains the following fields
%               residuals.residuals = vector of residuals (compulsory
%               field)
%               residuals.y = response (compulsory field if interactive
%               brushing is used)
%               residuals.X = n-by-p matrix containing explanatory
%               variables(compulsory field if interactive
%               brushing is used)
%                Data Types - single|double
%
% Optional input arguments:
%               h : the axis handle of a figure where to send the resindexplot.
%                   This can be used to host the resindexplot in a subplot of a
%                   complex figure formed by different panels (for example a panel
%                   with residuals from a classical ols estimator and another
%                   with residuals from a robust regression: see example
%                   below).
%                   Example -'h',h1 where h1=subplot(2,1,1)
%                   Data Types - Axes object (supplied as a scalar)
%              x :  the vector to be plotted on the x-axis. Numeric vector.
%                   As default the sequence 1:length(residuals) will be
%                   used
%                   Example -'x',1:100
%                   Data Types - double
%           labx :  a label for the x-axis. Character.  (default: '')
%                   Example -'labx','row index'
%                   Data Types - char
%           laby :  a label for the y-axis.  Character.  (default: '')
%                   Example -'laby','scaled residuals'
%                   Data Types - char
%          title :  a label containing the title of the plot.  Character. Default value is
%                   'Index plot of residuals'
%                   Example -'title','scaled residuals'
%                   Data Types - char
%          numlab:  number of points to be identified in plots.
%                   [] | cell ({5}) default) | numeric vector | structure.
%                   NUMLAB IS A CELL.
%                   If numlab is a cell containing scalar k, the units
%                   with the k largest residuals are labelled in the plots.
%                   The default value of numlab is {5}, that is the units
%                   with the 5 largest residuals are labelled.
%                   For no labelling leave it empty.
%                   NUMLAB IS A VECTOR.
%                   If numlab is a vector, the units inside vector numlab are
%                   labelled in the plots. 
%                   NUMLAB IS A STRUCTURE.
%                   If numlab is a struct it is possible to control the
%                   size of the points identified. It contains the
%                   following fields:
%                   numlab.numlab = number of points to be identified (cell
%                   or vector, see above);
%                   numlab.FontSize = fontsize of the labels of the
%                   points. The default value is 12.
%                   Example -'numlab',[3,10,35]
%                   Data Types - double
%        conflev :  confidence interval for the horizontal bands. Numeric
%                   vector.
%                   It can be a vector of different confidence level
%                   values.
%                   Example -'conflev',[0.95,0.99,0.999]
%                   Data Types - double
%                   Remark: confidence interval is based on the chi^2 distribution
%        FontSize:  Scalar which controls the fontsize of the labels of the
%                   axes. Default value is 12.
%                   Example -'Fontsize',10
%                   Data Types - double
%     SizeAxesNum:  Scalar which controls the fontsize of the numbers of
%                   the axes. Default value is 10.
%                   Example -'SizeAxesNum',10
%                   Data Types - double
%           ylimy:  Vector with two elements which controla minimum and maximum
%                   value of the y axis. Default is '', automatic scale.
%                   Example -'ylimy',[-5 5]
%                   Data Types - double
%           xlimx:  Vector with two elements controlling minimum and maximum
%                   on the x axis. Default value is '' (automatic scale).
%                   Example -'xlimx',[-5 5]
%                   Data Types - double
%          lwdenv:  width of the lines associated
%                   with the envelopes. Scalar. Default is lwdenv=1.
%                   Example -'lwdenv',2
%                   Data Types - double
%      MarkerSize:  size of the marker in points. Scalar.
%                   The default value for MarkerSize is 6 points (1 point =
%                   1/72 inch).
%                   Example -'MarkerSize',10
%                   Data Types - double
% MarkerFaceColor:  Marker fill color.
%                   'none' | 'auto' | RGB triplet | color string.
%                   Fill color for markers that are closed shapes
%                   (circle, square, diamond, pentagram, hexagram, and the
%                   four triangles).
%                   Example -'MarkerFaceColor','b'
%                   Data Types - char
%    databrush  :   interactive mouse brushing. Empty value, scalar or structure.
%                   If databrush is an empty value (default), no brushing
%                   is done.
%                   The activation of this option (databrush is a scalar or
%                   a cell) enables the user  to select a set of
%                   trajectories in the current plot and to see them
%                   highlighted in the y|X plot, i.e. a matrix of scatter
%                   plots of y against each column of X, grouped according
%                   to the selection(s) done by brushing. If the plot y|X
%                   does not exist, it is automatically created.
%                   Please, note that the window style of the other figures is set
%                   equal to that which contains the monitoring residual
%                   plot. In other words, if the monitoring residual plot
%                   is docked all the other figures will be docked too
%                   DATABRUSH IS A SCALAR.
%                   If databrush is a scalar the default selection tool is
%                   a rectangular brush and it is possible to brush only
%                   once (that is persist='').
%                   DATABRUSH IS A STRUCTURE.
%                   If databrush is a structure, it is possible to use all
%                   optional arguments of function selectdataFS.m and the
%                   following fields
%                   - databrush.persist = repeated brushng enabled. Persist is an empty value or a scalar
%                     containing the strings 'on' or 'off'.
%                     The default value of persist is '', that is brushing
%                     is allowed only once.
%                     If persist is 'on' or 'off' brushing can be done as
%                     many time as the user requires.
%                     If persist='on' then the unit(s) currently brushed
%                     are added to those previously brushed. it is
%                     possible, every time a new brushing is done, to use a
%                     different color for the brushed units.
%                     If persist='off' every time a new brush is performed
%                     units previously brushed are removed.
%                   - databrush.labeladd = add labels of brushed units.
%                     Character. [] (default) | '1'.
%                     If databrush.labeladd='1', we label the units
%                     of the last selected group with the unit row index in
%                     matrices X and y. The default value is labeladd='',
%                     i.e. no label is added.
%                   - databrush.bivarfit = this option adds one or more least
%                     square lines based on SIMPLE REGRESSION to the plots
%                     of y|X, depending on the selected groups.
%                     bivarfit = ''
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
%                   - databrush. multivarfit = this option adds one or more least square
%                       lines, based on MULTIVARIATE REGRESSION of y on X,
%                       to the plots of y|Xi.
%                     multivarfit = ''
%                       is the default: no line is fitted.
%                     multivarfit = '1'
%                       fits a single ols line to all points of each
%                       bivariate plot in the scatter matrix y|X.
%                       The line added to the scatter plot y|Xi
%                       is avconst +Ci*Xi, where Ci is the
%                       coefficient of Xi in the multivariate regression
%                       and avconst is the effect of all the other
%                       explanatory variables different from Xi evaluated
%                       at their centroid (that is overline{y}'C))
%                     multivarfit = '2'
%                       exactly equal to multivarfit ='1' but this time we
%                       add the line based on the group of unselected
%                       observations.
%                   - databrush.labeladd = if this option is '1', we label the units
%                     of the last selected group with the unit row index in
%                     matrices X and y. The default value is labeladd='',
%                     i.e. no label is added.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%       nameX   :   regressor labels. Cell array of strings of length p containing the labels
%                   of the variables of the regression dataset. If it is
%                   empty (default) the sequence X1, ..., Xp will be created
%                   automatically
%                   Example - 'nameX',{'Age','Income','Married','Profession'}
%                   Data Types - cell
%       namey   :   response label. Character. Character containing the
%                   label of the response. If it is
%                   empty (default) label 'y' will be used.
%                   Example - 'namey','response'
%                   Data Types - char
%           tag  :  Figure tag. Character.
%                   Tag of the figure which will host the malindexplot. The
%                   default tag is pl_resindex
%                   Example - 'tag','indexPlot'
%                   Data Types - character
%
% Output:
%
%
% See also resfwdplot.m
%
% References:
%
%   Rousseeuw P.J., Leroy A.M. (1987), "Robust regression and outlier
%   detection", Wiley.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('resindexplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Residual plot of 100 random numbers.
    resindexplot(randn(100,1))
%}

%{
    % Compare OLS residuals with robust residuals for the stack loss data.
    load('stack_loss.txt');
    y=stack_loss(:,4);
    X=stack_loss(:,1:3);
    % Define confidence level
    conflev=[0.95,0.99];
    figure;
    h1=subplot(2,1,1);
    % Compute studentized residuals (deletion residuals)
    stats=regstats(y,X,'linear',{'standres','studres'});
    resindexplot(stats.studres,'h',h1,'conflev',conflev,'labx','Index number','laby','Deletion residuals');

    % Compute robust residuals
    [out]=LXS(y,X,'nsamp',0,'rew',1,'lms',0);
    h2=subplot(2,1,2);
    resindexplot(out.residuals,'h',h2,'conflev',conflev,'labx','Index number','laby','Robust LTS reweighted residuals');
%}

%{
    % Just plot robust residuals.
    [out]=LXS(y,X,'nsamp',0,'rew',1,'lms',0);
    bonfconf = 1-0.01/size(y,1);    % 99% Bonferronised
    resindexplot(out.residuals,'conflev',[0.95,0.99,bonfconf],'labx','Index number','laby','Robust LTS reweighted residuals');
%}

%{
    % Interactive_example.
    databrush=struct;
    databrush.selectionmode='Brush'; % Brush selection
    databrush.persist='on'; % Enable repeated mouse selections
    databrush.Label='on'; % Write labels of the units while selecting
    databrush.RemoveLabels='on'; % Remove labels after selection
    databrush.RemoveTool    = 'on'; % Remove yellow tool after selection
    databrush.RemoveFlagged = 'on'; % Remove filled red color for selected points after selection

    [out]=LXS(y,X,'rew',1,'lms',0,'yxsave',1);
    resindexplot(out,'databrush',databrush)

    [outFS]=FSReda(y,X,out.bs);
    resfwdplot(outFS,'databrush',databrush)
%}


%{
    % Example of usage of option numlab.
    % Write the row number for the units which have the 3 largest
    % residuals (in absolute value)
    [out]=LXS(y,X,'nsamp',1000);
    resindexplot(out.residuals,'numlab',{3});
%}

%{
    % First example in which numlab is passed as structure.
    % In this case we control the FontSize of the associated labels.
    numlab=struct;
    % Set a font size for the labels equal to 20
    numlab.FontSize=20;
    resindexplot(randn(100,1),'numlab',numlab)
%}

%{
    % Second example in which numlab is passed as structure.
    % In this case we control both the number of units to label and
    % also the FontSize of the associated labels.
    numlab=struct;
    % Show just the two most important residuals.
    numlab.numlab={2}; 
    % Set a font size for the labels equal to 20
    numlab.FontSize=20;
    resindexplot(randn(100,1),'numlab',numlab)
%}

%% Initialization

if nargin<1
    error('FSDA:resindexplot:missingInputs','In order to run this function a vector of residuals has to be supplied')
end

if isstruct(residuals)
    out=residuals;
    residuals=out.residuals;
end

% The following line is to make sure residuals is always a column vector
residuals = residuals(:);

% Close existing pl_residual figure.
% Remark: existing figures containing subplots generated with resindexplot
% (using option 'h') will not be closed.
if ~isempty(findobj('type','figure','Tag','pl_resindex'))
    close(findobj('type','figure','Tag','pl_resindex'));
end

% n is the number of observations, as usual;
n=length(residuals);
numlabdef={{5}};
% Set standard options
options=struct('h','','x',1:n,'labx','','laby','','numlab',numlabdef,'conflev',0.975,...
    'title','Index plot of residuals','FontSize',12,'SizeAxesNum',10,...
    'xlimx','','ylimy','','lwdenv',1,'MarkerSize',6,'MarkerFaceColor','w',...
    'databrush','','tag','pl_resindex','nameX','','namey','');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:resindexplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
end


[h, x, labx, laby, titl, numlab, conflev, FontSize, SizeAxesNum, lwdenv,MarkerSize,MarkerFaceColor] = ...
    deal(options.h, options.x, options.labx, options.laby,...
    options.title, options.numlab, options.conflev,...
    options.FontSize, options.SizeAxesNum, options.lwdenv,options.MarkerSize,options.MarkerFaceColor);

% conflev
numconflev = length(conflev);
conflev = sort(conflev,'descend');

numlabFontSizedef=10;
if isstruct(numlab)
    if isfield(numlab,'FontSize')
        numlabFontSize=numlab.FontSize;
    else
        numlabFontSize=numlabFontSizedef;
    end
    
    if isfield(numlab,'numlab')
        numlab=numlab.numlab;
    else
        numlab=numlabdef{:};
    end
    
else
    numlabFontSize=numlabFontSizedef;
end


% numlab: if it is a cell, extract the number
if iscell(numlab)
    numlab=numlab{:};
    ord=abs(residuals);
    [~,ind]=sort(ord);
    ind=ind(n-numlab+1:n);
else
    ind=numlab;
    if size(ind,2)>1
        ind=ind';
    end
end

% Create the figure that will host the resindexplot
hfig = figure('Name', 'Residual plot', 'NumberTitle', 'off',...
    'Tag','pl_resindex');

% Get figure's axis
afig = axes('Parent',hfig);

% Plot the resindexplot and add relevant labels
plot(afig,x,residuals,'bo','MarkerFaceColor',MarkerFaceColor,...
    'MarkerSize',MarkerSize,'tag','data_res');

% Set the font size for the axes numbers
set(afig,'FontSize',SizeAxesNum);

text(x(ind),residuals(ind),int2str(ind),'VerticalAlignment','Baseline', 'FontSize',numlabFontSize);
% dx=(max(x)-min(x))/80; dy=(max(residuals)-min(residuals))/80;
% Displacement: baseline does already well the job of displacing the text.

% set the colors for the confidence bands
ColorOrd=[1 0 0 ; 0 1 1 ; 1 0 1; 1 1 0; 0 0 0; 0 1 0; 0 0 1];
ColorOrd1=ColorOrd(1:numconflev,:);
ColorOrd2=repmat(ColorOrd1,numconflev,1);
set(afig,'ColorOrder',ColorOrd2);

% set the string legend for the confidence bands
legendstring = cell(numconflev,1);
legendstring(:) = cellstr('% band');
legendstring2 = strcat(num2str(((conflev)*100)'),legendstring);

% set the confidence bands values
rangeaxis=axis;
quant = sqrt(chi2inv(conflev,1));
V=repmat([rangeaxis(1);rangeaxis(2)],1,2*numconflev);
QUANT=[[quant;quant],[ -quant;-quant]];

% plot the confidence bands
hline = line(V, QUANT,'LineWidth',lwdenv,'Tag','conflevline');
set(hline(1:numconflev),{'Displayname'},legendstring2);
set(hline(numconflev+1 : numconflev*2),{'Displayname'},legendstring2);

% make the legend for the confidence bands clickable
clickableMultiLegend(hline(1:numconflev),legendstring2);

axis(axis);
% Remark: the previous line is used to fix the y-axis, otherwise it may happen
% that the figure changes if one removes the bands by clicking on the legend

if ~isempty(h)
    % Eventually send the resindexplot into a different figure/subplot
    hfigh = get(h,'Parent');
    
    set(hfigh,'Name','Residual plots','NumberTitle','off');
    set(h,'Tag','res_subplot');
    copyobj(allchild(afig),h);
    pause(0.0000001);
    delete(hfig);
    hline2 = findobj(h, 'Tag','conflevline');
    hlineh = flipud(hline2);
    if length(findobj(get(h,'Parent'),'Tag','res_subplot'))==1
        clickableMultiLegend(hlineh(1:numconflev),legendstring2);
    else
        legend_h = legend(hlineh(1:numconflev),legendstring2);
        %         verMatlab=verLessThan('matlab','8.4.0');
        %         if verMatlab
        %             legend(legend_h,'hide');
        %         else
        %             legend_h.Visible='off';
        %         end
        set(legend_h,'Visible','off')
        
    end
    % Fix the y-axis
    set(h,'YLimMode', 'manual');
    % Add title and axis labels for the figure with subplots, and set their FontSize
    title(gca,titl);
    xlabel(gca,labx,'Fontsize',FontSize);
    ylabel(gca,laby,'Fontsize',FontSize);
    % Set the font size for the axes numbers
    set(gca,'FontSize',SizeAxesNum);
    
else
    % If the resindexplot has not to be sent in a different figure/subplot
    % add the figure title and axis labels, and set their FontSize
    title(afig,titl);
    xlabel(afig,labx,'Fontsize',FontSize);
    ylabel(afig,laby,'Fontsize',FontSize);
end

% set the axis limits
if ~isempty(options.ylimy)
    ylim(options.ylimy);
end
if ~isempty(options.xlimx)
    xlim(options.xlimx);
else
    xlim([0 n+1])
end

% The options labeladd, bivarfit, persist and multivarfit, must be removed
% from cell options.databrush, because they are not valid options for
% selectdataFS.
% Remark: "d=find(strcmp(..." is equivalent to "d=strmatch(...",
% but strmatch is slower and will be obsolete from 2011.

databrush=options.databrush;

if isstruct(databrush)
    fdatabrush=fieldnames(databrush);
    % labeladd option
    d=find(strcmp('labeladd',fdatabrush));
    if d>0
        labeladd=databrush.labeladd;
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
        fdatabrush=fieldnames(databrush);
        
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

% Initialize stype
styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
styp=repmat(styp,ceil(n/13),1);


%% Brush mode (call to function selectdataFS)
if ~isempty(options.databrush) || isstruct(options.databrush)
    
    
    outnames=fieldnames(out);
    d=strcmp('y',outnames)+strcmp('X',outnames);
    if max(d)<1
        mess=['Missing data in input structure.\n',...
            'With databrush option, this function requires an input structure \n',...
            'including two fields y and X containing the original data \n',...
            'For example, if your input structure is INP and your data matrix X and y  are MYDATAX MYDATAy\n',...
            'run first the following instruction: INP.y = MYDATAy; INP.X = MYDATAX; \n',...
            'Note that the estimator which generated the Mahalanobis distances has an option \n',...
            'to automatically save the input data in the output structure.   \n',...
            ''];
        error('FSDA:resindexplot:InvalidArg1',mess);
    else
        y=out.y;
        X=out.X;
    end
    
    
    if isstruct(options.databrush)
        % transform the input structure databrush into a cell array
        cv=[fdatabrush struct2cell(databrush)]';
        sele=[cv(:)' 'Ignore' {findobj(gcf,'tag','env')}];
        % add the FlagSize of the brushed points if it has not been
        % previously specified by the user
        d=find(strcmp('FlagSize',fdatabrush));
        if d==0
            sele=[sele 'FlagSize' '3'];
        end
    else
        sele={'selectionmode' 'Rect' 'Ignore' findobj(gcf,'tag','env') };
    end
    
    sele=[sele 'Tag' {options.tag}];
    
    
    % Check if X includes the constant term for the intercept.
    p=size(X,2);
    intcolumn = find(max(X,[],1)-min(X,[],1) == 0);
    
    if intcolumn==1
        intercept = 1;
        p1=1:(p-numel(intcolumn));
        Xsel=X;
        Xsel(:,intcolumn)=[];
    else
        intercept = 0;
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
    
    % group = vector which will contain the identifier of each group e.g.
    % group(14)=3 means that unit 14 was selected at the third brush
    group=ones(n,1);
    
    % some local variables
    but=0; brushcum=[]; ij=1;
    
    
    sele=[sele 'FlagColor' ColorOrd(ij,:) 'FlagMarker' char(styp(ij+1))];
    
    plot1=gcf;
    % set(plot1,'tag','data_res')
    
    % loop brushing
    while but<=1
        
        figure(plot1);
        
        % Remark: function selectdataFS cannot be used on the current
        % figure if the "selection mode" or the "zoom tool" are on. Setting
        % the plotedit mode initially to on and then to off, has the effect
        % to deselect plotedit mode.
        plotedit on
        plotedit off
        
        if strcmp(persist,'on')
            % add to cell sele option FlagColor (color of selection) and
            % FlagMarker (symbol to be used for selection)
            
            if ij>1
                chkexist=find(strcmp('FlagColor',sele)==1);
                sele{chkexist+1}=ColorOrd(ij,:);
                sele{chkexist+3}=char(styp(ij+1));
            end
        end
        
        % call to selectdataFS
        disp('Select a region to brush in the index plot of residuals');
        pl = selectdataFS(sele{:});
        
        % exit if the resfwdplot figure was closed before selection
        if isnumeric(pl) && ~isempty(pl) && (pl == -999)
            return
        end
        
        if ~isempty(cell2mat(pl))
            
            nbrush=pl{1};
            
            disp('Brushed units, yvalue and X values');
            disp([nbrush out.y(nbrush) out.X(nbrush,:)]);
            
        else
            disp('Wrong selection: Try again');
            disp('Select a region to brush in the residual plot');
            figure(plot1);
            nbrush='';
        end
        
        %% For each brushing operation, do the following:
        if ~isempty(nbrush)
            % brushcum = - the list of selected observations in all
            % iterations if
            %   persist=on
            % - the list of selected observations in the current iteration
            %   if persist=off
            if strcmp(persist,'on')
                brushcum=unique([brushcum; nbrush]);
            else
                brushcum=nbrush;
                group=ones(n,1);
            end
            
            % group=vector of length(Xsel) observations taking values from
            % 1 to the number of groups selected. unigroup= list of
            % selected groups.
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
                % the resfwdplot.
                figure;
                set(gcf,'WindowStyle',get(plot1,'WindowStyle'));
            end
            [H,AX,BigAx] = gplotmatrix(Xsel,y,group,clr(unigroup),char(styp{unigroup}),[],'on',[],nameX,namey);
            
            % Assign to this figure a name and a tag=pl_yX
            set(gcf,'Name','Scatter plot matrix y|X with selected groups highlighted');
            set(gcf,'tag','pl_yX');
            
            % Set markers
            for mfc=1:length(unigroup)
                set(findobj(gcf,'marker',char(styp(unigroup(mfc)))),'MarkerFaceColor',clr(unigroup(mfc)));
            end
            
            % Set the legenda properties of the gplotmatrix
            set(H(:,:,1),'DisplayName','Unbrushed units');
            for brugrp = 2:length(unigroup)
                set(H(:,:,brugrp),'DisplayName',['Brushed units ' num2str(brugrp-1)]);
            end
            
            % save the indices of the last selected units (nbrush) to the
            % 'UserData' field of the last selected group of H(:,:,end)
            set(H(:,:,end), 'UserData' , nbrush);
            
            % add objects to the panels of the yX
            % add2yX(H,AX,BigAx,out,group,nbrush,bivarfit,multivarfit,labeladd);
            add2yX(H,AX,BigAx,'intercept',intercept,'bivarfit',bivarfit,'multivarfit',multivarfit,'labeladd',labeladd);
            
            %% - check condition to exit from the brush mode
            % If the option persistent is not equal off or on than get out
            % of the loop
            if strcmp('on',persist) || strcmp('off',persist)
                if strcmp('on',persist)
                    ij=ij+1;
                    % but=1 makes that previous highlightments in other
                    % figures are not deleted
                    but=1;
                end
                
                % Before waitforbuttonpress:
                % - the resfwdplot is highlighted again
                figure(plot1);
                % - a function to be executed on figure close is set
                set(gcf,'CloseRequestFcn',@closereqFS);
                
                % - and lay down the plots before continuing
                position(plot1);
                disp('Highlight the index plot of residuals then: click on it to continue brushing or press a keyboard key to stop');
                ss=waitforbuttonpressFS;
                disp('------------------------');
                
                % After waitforbuttonpress: - the standard MATLAB function
                % to be executed on figure
                % close is recovered
                set(gcf,'CloseRequestFcn','closereq');
                Open_yX = findobj(0, 'type', 'figure','tag','pl_yX');
                Open_resindex = findobj(0, 'type', 'figure','tag','pl_resindex');
                if isempty(Open_resindex)  % User closed the main brushing window
                    if ~isempty(Open_yX); delete(Open_yX); end    % yX plot is deleted
                    delete(get(0,'CurrentFigure')); % deletes Figure if still one left open
                end
                
                % - and the 'but' variable is set if keyboard key was
                % pressed
                if ss==1
                    but=2;
                end
            else
                but=2;
            end % close loop associated with persist 'on' or persist 'off'
            position(plot1);
        end % for each brushing operation do ...
    end % close loop associated with but (loop brushing)
end % close options.databrush


end
%FScategory:VIS-Reg

