function MCDenv=malindexplot(md,v,varargin)
%malindexplot plots the Mahalanobis distances versus a selected variable.
%The selected variable is typically a generic index number.
%
%<a href="matlab: docsearchFS('malindexplot')">Link to the help function</a>
%
% Required input arguments:
%
%   md : Mahalanobis distances. Vector or structure. Vector of Mahalanobis
%        distances (in squared units) or a structure containing fields md
%        and Y. In this second case md is a structure with the following
%        fields:
%        md.md = contains the Mahalanobis distances (this field is compulsory);
%        md.Y = contains the original data matrix whose Mahalanobis
%        distances have been computed (this field is compulsory is option
%           databrush is used).
%        md.class = this field is not compulsory. In the case of
%        md.class='mcdCorAna' simulated envelopes are used to define the
%           empirical quantiles. Note that if the simulated bands have been
%           precalculated they can be passed through the seconf input
%           argument v.
%                Data Types - single|double
%
%  v : Number of variables or matrix of size n-by-k containing empirical envelope.
%       Scalar or matrix with the same rows of length(md).
%       If v is a scalar, it contains the number of variables of the
%       original data matrix which have been used to compute md. The
%       threshold in this case is based on the Chi^2 distribution with v
%       degrees of freedom. If v is a matrix with size(v,1)=length(md)
%       the empirical precalculated envelope in v are used to obtain the
%       confidence bands.
%                Data Types - single|double
%
% Optional input arguments:
%
%
%               h : Where to plot. Axis hadle.
%                   The axis handle of the Figure where to send the
%                   malindexplot. This can be used to host the malindexplot
%                   in a subplot of a complex figure formed by different
%                   panels (e.g. a panel with malindexplot from a classical
%                   mle estimator and another with Mahalanobis distances
%                   from a robust analysis, see example below).
%                   Example - 'h',gca
%                   Data Types - graphics handle
%
%              x :  x-axis index. Vector. The vector to be plotted on the
%                   x-axis.
%                   Default is the sequence 1:length(md).
%                   Example - 'x','1:100'
%                   Data Types - numeric
%
%           labx :  x label. Character. A label for the x-axis (default: '').
%                   Example - 'labx','unit number'
%                   Data Types - character
%
%           laby :  y label. Character. A label for the y-axis (default: '').
%                   Example - 'laby','MD'
%                   Data Types - character
%
%          title :  plot title. Character. A label containing the title of the plot.
%                   Default is 'Index plot of Mahalanobid distances'.
%                   Example - 'title','Index plot of MD'
%                   Data Types - character
%
%          numlab:  number of points to be labelled in the plot.
%                   vector or cell. If numlab is a cell containing scalar k, the units
%                   with the k largest md are labelled in the plots.
%                   If numlab is a vector, the units indexed by the vector
%                   are labelled in the plot.
%                   Default is numlab={5}, that is units with the 5
%                   largest md are labelled.
%                   Use numlab='' for no labelling.
%                   Example - 'numlab',{3}
%                   Data Types - numeric vector or cell.
%
%        conflev :  confidence interval for the horizontal bands. Vector.
%                   It can be a vector of different confidence level values,
%                   e.g. [0.95,0.99,0.999]. Confidence interval is based on
%                   the chi^2 distribution.
%                   Example - 'conflev',0.99
%                   Data Types - numeric
%
%        FontSize:  Labels font size. Scalar. Scalar which controls the
%                   font size of the labels of the axes.
%                   Default value is 12.
%                   Example - 'FontSize',12
%                   Data Types - numeric
%
%     SizeAxesNum:  Numbers font size. Scalar. Scalar which controls the
%                   fontsize of the numbers of  the axes.
%                   Default value is 10.
%                   Example - 'SizeAxesNum',12
%                   Data Types - numeric
%
%           ylimy:  ylimits. Vector. Vector with two elements controlling minimum and
%                   maximum value of the y axis.
%                   Default is '' (automatic scale).
%                   Example - 'ylimiy',[-3 3]
%                   Data Types - numeric
%
%           xlimx:  xlimits. Vector. Vector with two elements controlling minimum and
%                   maximum value of the x axis.
%                   Default is '' (automatic scale).
%                   Example - 'xlimix',[1 30]
%                   Data Types - numeric
%
%          lwdenv:  Envelope line width. Scalar. Scalar which controls the
%                   width of the lines associated  with the envelopes.
%                   Default is lwdenv=1.
%                   Example - 'lwdenv',4
%                   Data Types - numeric
%
%      MarkerSize:  Marker size of points. Scalar. Scalar specifying the
%                   size of the marker in points (1 point = 1/72 inch).
%                   Default is MarkerSize = 6.
%                   Example - 'MarkerSize',4
%                   Data Types - numeric
%
% MarkerFaceColor:  Marker fill color of points. Character or length 3 RGB
%                   numeric vector. The fill color for markers that are closed shapes
%                   (circle, square, diamond, pentagram, hexagram, and the
%                   four triangles).
%                   Example - 'MarkerFaceColor','b'
%                   Data Types - numeric | character
%
%             tag:  Figure tag. Character.
%                   Tag of the figure which will host the malindexplot.
%                   The default tag is pl_malindex.
%                   Example - 'tag','indexPlot'
%                   Data Types - character
%
%    databrush  :   interactive mouse brushing. Empty value, scalar or structure.
%                   If databrush is an empty value (default),
%                   no brushing is done. The activation of this option
%                   (databrush is a scalar or a structure) enables the user
%                   to select a set the points in the current plot and to see them
%                   highlighted in the scatter plot matrix (spm). If spm
%                   does not exist it is automatically created.
%                   DATABRUSH IS A SCALAR.
%                   If databrush is a scalar the default
%                   selection tool is a rectangular brush and it is
%                   possible to brush only once (that is persist='').
%                   DATABRUSH IS A STRUCTURE.
%                   If databrush is a structure,
%                   it is possible to use all optional arguments of
%                   function selectdataFS.m and the following optional
%                   arguments:
%                   databrush.persist = persisent brushing.
%                     Persist is an empty value or a scalar
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
%                   databrush.labeladd = add labels. If this option is '1',
%                     we label in the scatter plot matrix the units of the
%                     last selected group with the unit row index in matrix
%                     Y. The default value is labeladd='', i.e. no label is
%                     added.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%
%                   REMARK: the options which follow work in connection
%                   with previous option databrush and produce their effect
%                   on the scatter plot matrix of the original data.
%
%           nameY: variables labels of the original data matrix. Cell. Cell
%                array of strings containing the labels of the
%                variables. As default value, the labels which are added
%                are Y1, ..., Yv. This option is used just if previous
%                option databrush is not empty.
%                   Example - 'nameY',{'Y_1' Y_2'}
%                   Data Types - character
%
%      label : row labels. Cell.
%               Cell of length n containing the labels of the rows.
%                   Example - 'label',{'UK' ...  'IT'}
%                   Data Types - cell
%
%
%
%
%  Output:
%
%           MCDenv : Empirical envelopes. Array.
%                   Matrix with size n-by-length(conflev) which contains
%                   the empirical confidence envelopes or vector of length
%                   length(conflev) containing teh quantiles of the
%                   reference distribution.
%
%
% See also resfwdplot.m, resindexplot.m
%
% References:
%
%   Rousseeuw P.J., Leroy A.M. (1987), "Robust regression and outlier
%   detection", Wiley.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('malindexplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %   Mahalanobis distance plot of 100 random numbers.
    %   Numbers are from from the chi2 with 5 degrees of freedom
     MCDenv=malindexplot(chi2rnd(5,100,1),5);
%}

%{
    %   Compare traditional md with robust md for the stack loss data.
    load('stack_loss.txt');
    X=stack_loss(:,1:3);
    [n,v]=size(X);
    % Define confidence level
    conflev=[0.95,0.99];
    figure;
    h1=subplot(2,1,1);
    % Compute traditional Mahalanobis distances
    mdtrad=mahal(X,X);
    malindexplot(mdtrad,v,'h',h1,'conflev',conflev,'labx','Index number','laby','Traditional md');

    % Compute robust md
    [out]=FSM(X,'init',5,'plots',0);
    seq=1:size(X,1);
    good=setdiff(seq,out.outliers);
    mdrob=mahal(X,X(good,:));
    
    h2=subplot(2,1,2);
    malindexplot(mdrob,v,'h',h2,'conflev',conflev,'labx','Index number','laby','Robust md','title','');

%}


%{
    % Interactive_example
    %   Index plot Mahalanobis distance with databrush option.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [RAW,REW]=mcd(Ycont);
    RAW.Y=Ycont;
    malindexplot(RAW,v,'databrush',1)
%}

%{
    % Interactive_example
    %   Index plot Mahalanobis distance with personalized databrush option.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [RAW,REW]=mcd(Ycont);
    RAW.Y=Ycont;
    databrush=struct;
    databrush.selectionmode='Brush'; % Brush selection
    databrush.persist='on'; % Enable repeated mouse selections
    databrush.Label='on'; % Write labels of the units while selecting
    databrush.RemoveLabels='on'; % Remove labels after selection
    databrush.RemoveTool    = 'off'; % Do not remove yellow tool after selection
    databrush.RemoveFlagged = 'off'; % Do not remove filled red color for selected points after selection
    databrush.labeladd = '1'; % Write number of seleceted units in the scatter plot matrix
    malindexplot(RAW,v,'databrush',databrush)
%}



%% Beginning of code

if nargin<1
    error('FSDA:malindexplot:missingInputs','To run this function a vector of Mahalanobis distances has to be supplied')
end

if nargin<2
    error('FSDA:malindexplot:missingInputs','To run this function the number of variables which have been used to construct md has to be supplied')
end

if isempty(v) || any(isnan(v))
    error('FSDA:malindexplot:Wrongv','v must be scalar (or vector with the same length of md) non empty and non missing!!!');
end


% Close existing pl_mahal figure. Remark: existing figures containing
% subplots generated with malindexplot (using option 'h') will not be closed.
if ~isempty(findobj('type','figure','Tag','pl_malindex'))
    close(findobj('type','figure','Tag','pl_malindex'));
end

if isstruct(md)
    
    if isfield(md,'N') && isfield(md,'h')
        % mcdCorAnatype
        out=md;
        hh=md.h;
        N=md.N;
        md=md.md;
        class='mcdCorAna';
    elseif isfield(md,'N')
        % non robust distances
        out=md;
        N=md.N;
        class='mcdCorAna';
        hh=sum(N,'all');
    else
        out=md;
        md=out.md;
        class='';
    end
else
    class='';
end

% The following line is to ensure md is always a column vector
md = md(:);

% n = number of observations;
n=length(md);

% Set standard options
options=struct('h','','x',1:n,'labx','','laby','','numlab',{{5}},'conflev',0.975,...
    'title','Index plot of Mahalanobis distances','FontSize',12,'SizeAxesNum',10,...
    'xlimx','','ylimy','','lwdenv',1,'MarkerSize',6,'MarkerFaceColor','w',...
    'databrush','','tag','pl_malindex','nameY','','label','');


UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:malindexplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
end

[h, x, labx, laby, titl, numlab, conflev, FontSize, SizeAxesNum, ...
    lwdenv,MarkerSize,MarkerFaceColor] = deal(...
    options.h, options.x, options.labx, options.laby,...
    options.title, options.numlab, options.conflev,...
    options.FontSize, options.SizeAxesNum, options.lwdenv,...
    options.MarkerSize,options.MarkerFaceColor);

% conflev
numconflev = length(conflev);
conflev = sort(conflev,'descend');

% numlab: if it is a cell, extract the number
if iscell(numlab)
    numlab=numlab{:};
    ord=abs(md);
    [~,ind]=sort(ord);
    ind=ind(n-numlab+1:n);
else
    ind=numlab;
    if size(ind,2)>1
        ind=ind';
    end
end

Tag=options.tag;

% Create the figure that will host the malindexplot
hfig = figure('Name', 'Residual plot', 'NumberTitle', 'off',...
    'Tag',Tag);


% Get figure's axis
afig = axes('Parent',hfig);
% Set the font size for the axes numbers
set(afig,'FontSize',SizeAxesNum);

% Plot the malindexplot and add relevant labels
plot(afig,x,md,'bo','MarkerFaceColor',MarkerFaceColor,...
    'MarkerSize',MarkerSize);

label=options.label;
if isempty(label)
    text(x(ind),md(ind),int2str(ind),'VerticalAlignment','Baseline');
else
    text(x(ind),md(ind),label(ind),'VerticalAlignment','Baseline');
end
% dx=(max(x)-min(x))/80; dy=(max(md)-min(md))/80;
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

V=[rangeaxis(1);rangeaxis(2)];

if isempty(class)
    if isscalar(v)
        quant = chi2inv(conflev,v);
        QUANT=[quant;quant];
        
        % plot the confidence bands
        hline = line(V, QUANT,'LineWidth',lwdenv,'Tag','conflevline');
        MCDenv=quant;
    else
        % plot the empirical confidence band(s)
        hline = line(x, v,'LineWidth',lwdenv,'Tag','conflevline');
        MCDenv=v;
    end
    
elseif strcmp(class,'mcdCorAna')
    % if the class is mcdCorAna simulated envelopes are used for each row
    [I,J]=size(N);
    % nrowt = column vector containing row marginal totals
    nrowt=sum(N,2);
    % ncolt = row vector containing column marginal totals
    ncolt=sum(N,1);
    
    nsimul=200;
    mmdStore=zeros(I,nsimul);
    
    parfor j=1:nsimul
        % Generate the contingency table
        out1=rcontFS(I,J,nrowt,ncolt);
        Nsim=out1.m144;
        
        RAW=mcdCorAna(Nsim,'plots',0,'msg',0,'bdp',hh);
        mmdStore(:,j)=RAW.md;
    end
    
    % Sort rows of matrix mmdStore
    mmdStore=sort(mmdStore,2);
    
    
    MCDenv=mmdStore(:,round(nsimul*conflev));
    
    % plot the confidence bands
    hline = line(1:I, MCDenv,'LineWidth',lwdenv,'Tag','conflevline');
else
    % Should never be here
    error('wrong class')
end

if isscalar(v)
    set(hline(1:numconflev),{'Displayname'},legendstring2);
    
    % make the legend for the confidence bands clickable
    clickableMultiLegend(hline(1:numconflev),legendstring2);
end

% fix the y-axis, otherwise the figure may change if one hides the bands
% by clicking on the legend
axis(axis);

if ~isempty(h)
    % Eventually send the malindexplot into a different figure/subplot
    hfigh = get(h,'Parent');
    
    set(hfigh,'Name','Mahalanobis distance plot','NumberTitle','off');
    set(h,'Tag','md_subplot');
    copyobj(allchild(afig),h);
    delete(hfig);
    hline2 = findobj(h, 'Tag','conflevline');
    hlineh = flipud(hline2);
    if length(findobj(get(h,'Parent'),'Tag','md_subplot'))==1
        clickableMultiLegend(hlineh(1:numconflev),legendstring2);
    else
        legend_h = legend(hlineh(1:numconflev),legendstring2);
        %         verMatlab=verLessThan('matlab','8.4.0');
        %         if verMatlab
        %             legend(legend_h,'hide');
        %         else
        %             legend_h.Visible='off';
        %         end
        set(legend_h,'Visible','off');
    end
    drawnow;
    % Fix the y-axis
    set(h,'YLimMode', 'manual');
    % Add title and axis labels for the figure with subplots
    title(gca,titl);
    xlabel(gca,labx,'Fontsize',FontSize);
    ylabel(gca,laby,'Fontsize',FontSize);
    % Set the font size for the axes numbers
    set(gca,'FontSize',SizeAxesNum);
    
else
    % If the plot has not to be sent in a different figure/subplot
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
end

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
    d=strcmp('Y',outnames);
    if max(d)<1
        mess=['Missing data in input structure.\n',...
            'With databrush option, this function requires an input structure \n',...
            'including a field Y containing the original data \n',...
            'For example, if your input structure is INP and your data matrix is MYDATA \n',...
            'run first the following instruction: INP.Y = MYDATA; \n',...
            'Note that the estimator which generated the Mahalanobis distances has an option \n',...
            'to automatically save the input data in the output structure.   \n',...
            ''];
        error('FSDA:malindexplot:InvalidArg1',mess);
    else
        Y=out.Y;
        v=size(Y,2);
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
    
    % Adding the instruction is necessary to make sure that the selection
    % is done in the  current plot
    sele=[sele 'Tag' {Tag}];
    
    
    % Create the labels
    p1=1:v;
    
    % Set the labels of the axes.
    if isempty(options.nameY)
        nameY=cellstr(num2str(p1','Y%d'));
    else
        nameY=options.nameY;
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
        disp('Select a region to brush in the index plot of Mahalanobis distances');
        pl = selectdataFS(sele{:});
        
        % exit if the malfwdplot figure was closed before selection
        if isnumeric(pl) && ~isempty(pl) && (pl == -999)
            return
        end
        
        if ~isempty(cell2mat(pl))
            
            nbrush=pl{1};
            
            disp('Brushed units, Yvalues ');
            disp([nbrush out.Y(nbrush,:)]);
            
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
            
            
            
            %% - display the spm with the corresponding groups of units highlighted
            
            h=findobj('-depth',1,'Tag','pl_spm');
            
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
                figure('Tag','pl_spm');
                set(gcf,'WindowStyle',get(plot1,'WindowStyle'));
            end
            
            plo=struct; plo.nameY=nameY; plo.labeladd=labeladd;
            H = spmplot(Y,group,plo);
            
            % Assign to this figure a name and a tag=pl_spm
            set(gcf,'Name','Scatter plot matrix with selected groups highlighted');
            
            % Set markers
            for mfc=1:length(unigroup)
                set(findobj(gcf,'marker',char(styp(unigroup(mfc)))),'MarkerFaceColor',clr(unigroup(mfc)));
            end
            
            % save the indices of the last selected units (nbrush) to the
            % 'UserData' field of the last selected group of H(:,:,end)
            set(H(:,:,end), 'UserData' , nbrush);
            
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
                disp('Highlight the index plot of MD then: click on it to continue brushing or press a keyboard key to stop');
                ss=waitforbuttonpressFS;
                disp('------------------------');
                
                % - the standard MATLAB function to be executed on figure
                %   close is recovered
                set(gcf,'CloseRequestFcn','closereq');
                Open_spm = findobj(0, 'type', 'figure','tag','pl_spm');
                Open_malindex = findobj(0, 'type', 'figure','tag','pl_malindex');
                if isempty(Open_malindex)  % User closed the main brushing window
                    if ~isempty(Open_spm); delete(Open_spm); end    % yX plot is deleted
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
%FScategory:VIS-Mult
