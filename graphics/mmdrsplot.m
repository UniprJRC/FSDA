function mmdrsplot(out,varargin)
%mmdrsplot plots the trajectories of minimum Mahalanobis distances from different starting points
%
%<a href="matlab: docsearchFS('mmdrsplot')">Link to the help function</a>
%
% Required input arguments:
%
%  out :  Structure containing the following fields.
%
%   out.mmdrs = a matrix of size (n-ninit)-by-(nsimul+1)containing the
%               monitoring of minimum Mahalanobis distance
%               in each step of the forward search for each of the nsimul random starts.
%               The first column of mmdrs must contain the fwd search
%               index.
%               This matrix can be created using function FSMmmdrs.
%
% out.BBrs   =  3D array of size n-by-n-(init)-by-nsimul containing units
%               forming subset for rach random start.
%               This field is necessary if datatooltip is true or databrush
%               is not empty.
%
%     out.Y   =   n-by-v matrix containing the original datamatrix
%               This field is necessary if datatooltip is true or databrush
%               is not empty.
%
% Optional input arguments:
%
%  quant:    quantiles for which envelopes have
%               to be computed. Vector.
%               1 x k vector containing quantiles for which envelopes have
%               to be computed. The default is to produce 1%, 50% and 99%
%               envelopes.
%               Example - 'quant',[0.01 0.99]
%               Data Types - double
%
%       envm    :  sample size to use. Scalar. Scalar which specifies the
%                   size of the sample which is used to superimpose the
%                   envelope. The default is to add an envelope based on
%                   all the observations (size n envelope).
%               Example - 'envm',n
%               Data Types - double
%
%       xlimx   :   min and max for x axis. Vector. vector with two
%                   elements controlling minimum and maximum on the x axis.
%                   Default value is mmd(1,1)-3 and mmd(end,1)*1.3
%                   Example - 'xlimx',[20 100]
%                   Data Types - double
%
%       ylimy   :   min and max for x axis. Vector. Vector with two
%                   elements controlling minimum and maximum on the y axis.
%                   Default value is min(mmd(:,2)) and max(mmd(:,2));
%                   Example - 'ylimy',[2 6]
%                   Data Types - double
%
%       lwdenv  :   Line width. Scalar. Scalar which controls the width of
%                   the lines associated with the envelopes. 
%                   Default is lwdenv=1
%                   Example - 'lwdenv',2
%                   Data Types - double
%
%       tag     :   plot handle. String. String which identifies the handle
%                   of the plot which is about to be created. The default
%                   is to use tag 'pl_mmdrs'. Notice that if the program
%                   finds a plot which has a tag equal to the one specified
%                   by the user, then the output of the new plot overwrites
%                   the existing one in the same window else a new window
%                   is created.
%                   Example - 'tag','mymmdrs'
%                   Data Types - char
%
%   datatooltip :   interactive clicking. Empty value (default) or
%                   structure.
%                   If datatooltip is not empty the user can use the mouse
%                   in order to have information about the unit selected,
%                   the step in which the unit enters the search and the
%                   associated label. If datatooltip is a structure, it is
%                   possible to control the aspect of the data cursor (see
%                   function datacursormode for more details or the
%                   examples below). The default options of the structure
%                   are DisplayStyle='Window' and SnapToDataVertex='on'.
%                   Example - 'datatooltip',1
%                   Data Types - empty value, numeric or structure
%
%       label   :   rwo labels. Cell. Cell containing the labels of the
%                   units (optional argument used when datatooltip=1. If
%                   this field is not present labels row1, ..., rown will
%                   be automatically created and included in the pop up
%                   datatooltip window).
%                   Example - 'label',{'Smith','Johnson','Robert','Stallone'}
%                   Data Types - cell
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
%                   Remark: if databrush is a struct, it is possible to
%                   specify all optional arguments of function selectdataFS
%                   inside the curly brackets of option databrush.
%
%       FontSize:   Size of axes labels. Scalar. Scalar which controls the
%                   fontsize of the labels of the axes. 
%                   Default value is 12.
%                   Example - 'FontSize',14
%                   Data Types - single | double
%
%    SizeAxesNum:   Size of axes numbers. Scalar which controls the
%                   fontsize of the numbers of the axes.
%                   Default value is 10.
%                   Example - 'SizeAxesNum',14
%                   Data Types - single | double
%
%       nameY   :   Regressors names. Cell array of strings. Cell array of
%                   strings of length v containing the labels
%                   of the varibales of the original data matrix. If it is empty
%                 	(default) the sequence Y1, ..., Yp will be created
%                   automatically.
%                   Example - 'nameY',{'Age','Income','Married','Profession'}
%                   Data Types - cell
%
%       lwd     :   Curves line width. Scalar. Scalar which controls
%                   linewidth of the curve which contains the monitoring of
%                   minimum Mahalanobis distance.
%                   Default line width=2
%                   Example - 'lwd',3
%                   Data Types - single | double
%
%       titl    :   main title. Character.
%                   A label for the title (default: '').
%                   Example - 'namey','Plot title'
%                   Data Types - char
%
%       labx    :   x axis title. Character.
%                   A label for the x-axis (default: 'Subset size m').
%                   Example - 'labx','Subset size m'
%                   Data Types - char

%       laby    :   y axis title. Character. A label for the y-axis
%                  (default: 'Minimum Mahalnobis distance').
%                   Example - 'laby','mmd'
%                   Data Types - char
%
%       labenv  :   label the envelopes. Boolean. 
%                   If labenv is true (default) labels of the confidence
%                   envelopes which are used are added on the y axis.
%                   Example - 'labenv',false
%                   Data Types - boolean
%
%        scaled :   scaled or unscaled envelopes. Boolean. 
%                   Use reference envelopes scaled or unscaled).
%                   If scaled=1 the envelopes are produced for
%                   scaled Mahalanobis distances (no consistency factor is
%                   applied) else the traditional consistency factor is applied
%                   Default is to use unscaled envelopes.
%                   Example - 'scaled',0
%                   Data Types - char
%
% Output:
%
% See also: FSMmmdrs.m
%
% References:
%
% Atkinson, A.C., Riani, M. and Cerioli, A. (2004), "Exploring multivariate
% data with the forward search", Springer Verlag, New York.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mmdrsplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of mmdrsplot with all the default options.
    %Steps common to all the examples
    load('swiss_banknotes');
    Y=swiss_banknotes.data;
    [fre]=unibiv(Y);
    %create an initial subset with the 4 observations, which fell the smallest
    %number of times outside the robust bivariate ellipses, and with the
    %lowest Mahalanobis Distance.
    fre=sortrows(fre,[3 4]);
    m0=20;
    bs=fre(1:m0,1);
    [outeda]=FSMeda(Y,bs);
    [out]=FSMmmdrs(Y,'bsbsteps',0,'cleanpool',0,'nsimul',80);
    mmdrsplot(out)

%}

%{
    %Example of the use of function mmdrsplot with personalized envelopes.
    mmdrsplot(out,'quant',[0.99;0.9999]);
%}

%{
    % Interactive_example
    % mmdrsplot with option dataooltip.
    %Example of the use of function mmdrsplot with datatooltip passed as
    %scalar (that is using default options for datacursor (i.e.
    %DisplayStyle =window)
     mmdrsplot(out,'datatooltip',1);
%}

%{
    % Interactive_example
    % mmdrsplot with option dataooltip passed as structure.
    %Example of the use of function mmdrsplot with datatooltip passed as
    %structure
    clear tooltip
    tooltip.SnapToDataVertex='on'
    tooltip.DisplayStyle='datatip'
    mmdrsplot(out,'datatooltip',tooltip);
%}

%{
    % Interactive_example
    %Example of the use of option databrush.
    mmdrsplot(out,'databrush',1);
%}

%{
    % Interactive_example
    % Example of the use of option databrush. 
    % Selected units are also highlighted in the malfwdplot.
    load('swiss_banknotes');
    Y=swiss_banknotes.data;

    out=FSMmmdrs(Y,'bsbsteps',0,'cleanpool',0,'nsimul',80);
    outEDA=FSMeda(Y,1:10,'init',20,'scaled',1);
    malfwdplot(outEDA)
    databrush=struct;
    databrush.persist='on';
    mmdrsplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    %% Two groups with approximately the same number of units.
    close all
    rng('default')
    rng(10);
    n1=100;
    n2=100;
    v=3;
    Y1=rand(n1,v);
    Y2=rand(n2,v)+1;
    Y=[Y1;Y2];
    group=[ones(n1,1);2*ones(n2,1)];
    spmplot(Y,group);
    title('Two simulated groups')
    Y=[Y1;Y2];
    close all
    % parfor of Parallel Computing Toolbox is used (if present in current computer)
    % and pool is not cleaned after the execution of the random starts
    % The number of workers which is used is the one specified
    % in the local/current profile
    [out]=FSMmmdrs(Y,'nsimul',100,'init',10,'plots',1,'cleanpool',0);
    ylim([2 5])
    disp('The two peaks in the trajectories of minimum Mahalanobis distance (mmd).')
    disp('clearly show the presence of two groups.')
    disp('The decrease after the peak in the trajectories of mmd is due to the masking effect.')
    mmdrsplot(out,'databrush',1)
   
%}

%% Beginning of code 

% Initialization
if nargin<1
    error('FSDA:mmdrsplot:missingInputs','A required input argument is missing.')
end

mmdrs=out.mmdrs;
ntrajectories=size(mmdrs,2)-1;

BBrs=out.BBrs;
nsimul=size(mmdrs,2)-1;

[n,v]=size(out.Y);

% seq= column vector containing the sequence 1 to n
% seq= (1:n)';

% numtext= a cell of strings used to labels the units with their position
% in the dataset.
% numtext=cellstr(num2str(seq,'%d'));

% Default limits for x axis
xl1=mmdrs(1,1)-3;
xl2=mmdrs(end,1)*1.1;
xlimx=[xl1 xl2];

% Default limits for y axis
yl1=min(min(mmdrs(:,2:end)));
yl2=max(max(mmdrs(:,2:end)))*1.1;

thresh=6;
if yl2 >thresh
    mdrrstmp=mmdrs(:,2:end);
    mdrrstmp(mdrrstmp>thresh)=NaN;
    yl2=max(max(mdrrstmp))*1.1;
else
end
ylimy=[yl1 yl2];

% Default quantiles to compute the envelopes
quant=[0.01;0.5;0.99];

labx='Subset size m';
laby='Minimum Mahalanobis distance';


%% User options

options=struct('quant', quant,...
    'xlimx',xlimx,'ylimy',ylimy,'lwd',2,'lwdenv',1,...
    'FontSize',12,'SizeAxesNum',10,'tag','pl_mmdrs',...
    'datatooltip','','databrush','',...
    'titl','','labx',labx,'laby',laby,'nameY','','label','','envm',n,'scaled',0,'labenv',true);

if nargin<1
    error('FSDA:mmdrsplot:missingInputs','A required input argument is missing.')
end

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mmdrsplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    labeladd='';
    persist='';
    ColorOrd=[1 0 0];
    clr='brcmykgbrcmykgbrcmykg';
    flagcol='r';
end

% Symbol types for spm plot (in case of brushing)
styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};

% Quantiles associated with the envelope based on all the observations.
quant=options.quant;
nenvel=length(quant);

% scaled = scaled or unscaled envelopes (default is unscaled envelopes)
scaled=options.scaled;

envm=options.envm;

% lwd = line width of the curve which contains mdr
lwd=options.lwd;

% lwdenv = line width of the curves associated with the envelopes
lwdenv=options.lwdenv;

% FontSize = font size of the axes labels
FontSize=options.FontSize;

% FontSizeAxes = font size for the axes numbers
SizeAxesNum=options.SizeAxesNum;

% labenv = add labels to the nevelopes
labenv=options.labenv;


init=mmdrs(1,1);
if init-v==0
    init=init+1;
end

%% Display the min MD random start plot

%Create a figure to host the plot or clear the existing one
h=findobj('-depth',1,'tag',options.tag);
if (~isempty(h))
    clf(h(end));
    figure(h(end))
    axes;
else
    % create a new figure
    h=figure;
end
set(h,'Name', 'Monitoring of Minimum Mahalnobis distance', 'NumberTitle', 'off');
hold('all');

% Theoretical envelopes for minimum Mahalnobis distance
[gmin] = FSMenvmmd(envm,v,'prob',quant,'init',init,'scaled',scaled);

box('on');

% set the x and y axis
xlimx=options.xlimx;
ylimy=options.ylimy;
xlim(xlimx);
ylim(ylimy);

% x coordinates where to put the messages about envelopes
xcoord=max([xlimx(1) init]);
for i=1:length(quant)
    % Superimpose chosen envelopes
    if quant(i)==0.5
        % Superimpose 50% envelope
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
    elseif quant(i)<=0.99
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
    else
        line(gmin(:,1),gmin(:,i+1),'LineWidth',lwdenv,'LineStyle','--','Color',[0  0 0],'tag','env');
    end
    
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
    
    if labenv == true
    annotation(gcf,'textbox',[figx figy kx ky],'String',{[num2str(100*quant(i)) '%']},...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'EdgeColor','none',...
        'BackgroundColor','none',...
        'FitBoxToText','off',...
        'FontSize',FontSize);
    end
end


box('on');

% sel= extract the elements of matrix mdr which have to be plotted
sel=1:size(mmdrs,1)-n+envm;

tagstat='rs_data_mmd';
% plot minimum deletion residual
plot1=plot(mmdrs(sel,1),mmdrs(sel,2:end),'tag',tagstat,...
    'LineWidth',lwd);

LineStyle={'-';'--';':';'-.'};
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));

% Write an extra message on the plot
annotation(gcf,'textbox',[0.2 0.8 0.1 0.1],'EdgeColor','none','String',['Envelope based on ' num2str(envm) ' obs.'],'FontSize',FontSize);


labx=options.labx;
laby=options.laby;
titl=options.titl;
title(titl);
% Add to the plot the labels for values of la
xlabel(labx,'Fontsize',FontSize);
ylabel(laby,'Fontsize',FontSize);

set(gca,'FontSize',SizeAxesNum)

set(gca,'FontSize',SizeAxesNum)

hold('off');

% displays the boundary of the current axes.
box on

% include specified tag in the current plot
set(gcf,'tag',options.tag);

% Store the handle of the mdrplot inside handle hmin
hmin=gcf;



%% Set the datatooltip for the mmdplot
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
    set(hdt,'UpdateFcn',{@mmdplotLbl,out})
end

%% Brush mode (call to function selectdataFS)

if ~isempty(options.databrush) || isstruct(options.databrush)
    
    if isstruct(options.databrush)
        
        %         % If option Label is 'on' then matrix Un is added to UserData
        %         d=max(strcmp('Label',fieldnames(databrush)));
        %         if d==1 && strcmp(databrush.Label,'on')
        %             set(gcf,'UserData',Un)
        %         end
        
        
        cv=[fieldnames(databrush) struct2cell(databrush)]';
        
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
    
    sele=[sele 'Tag' {options.tag}];
    
    
    
    % group = VECTOR WHICH WILL CONTAIN THE IDENTIFIER OF EACH GROUP
    % e.g. group(14)=3 means that unit 14 was selected at the third brush
    group=ones(n,1);
    
    % some local variables
    but=0; brushcum=[]; ij=1; ijk=1;
    
    % Extract Y
    Y=out.Y;
    % Set the labels of the axes.
    d=find(strcmp('nameY',fieldnames(out)),1);
    if  isempty(d)
        v=size(Y,2);
        p1=1:v;
        nameY=cellstr(num2str(p1','y%d'));
    else
        nameY=options.nameY;
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
                seleUnitjboo(abs(mmdrs(selxsteps(1)-init+1,2:end)-selysteps(1))<1e-10,j)=true;
                
                % selTraj contains the true indexes of the selected trajectories
                % using the ordering in which they are stored inside the figure
                selTraj=find(abs(a(selxsteps(1)-init+1,:)-selysteps(1))<1e-10);
                
                % disp('Indexes of the selected trajectories as they are stored inside the plot')
                % disp(selTraj)
                seleTrajboo(selTraj)=true;
                
                % seleTracumboo = boolean vector which contains true in
                % correspondence of the selected trajectories for all
                % brushing
                seleTracumboo(selTraj)=true;
                
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
            
            % Set new colors and new line style for the unselected trajectories
            set(htraj(~seleTracumboo),'Color',FSColors.greysh.RGB);
            set(htraj(~seleTracumboo),'LineStyle',':');
            
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
            lineindexes=strcmp(get(chH,'Tag'),'rs_data_mmd');
            
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
                brushcum=unique([brushcum; nbrush]);
            else
                brushcum=nbrush;
                group=ones(n,1);
            end
            
            % group=vector of length(Xsel) observations taking values
            % from 1 to the number of groups selected.
            % unigroup= list of selected groups.
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
                % the mdrplot.
                figure('Tag','pl_spm');
                set(gcf,'WindowStyle',get(hmin,'WindowStyle'));
            end
            
            plo=struct; plo.nameY=nameY; plo.labeladd=labeladd;
            plo.clr=clr(unigroup);
            plo.sym= char(styp{unigroup});
            H = spmplot(Y,group,plo);
            
            % Assign to this figure a name and a tag=pl_spm
            set(gcf,'Name','Scatter plot matrix Y with selected groups highlighted');
            
            % Set markers
            for mfc=1:length(unigroup)
                set(findobj(gcf,'marker',char(styp(unigroup(mfc)))),'MarkerFaceColor',clr(unigroup(mfc)));
            end
            
            % save the indices of the last selected units (nbrush) to the
            % 'UserData' field of the last selected group of H(:,:,end)
            set(H(:,:,end), 'UserData' , nbrush);
            
            %% - highlight brushed trajectories also in the malfwdplot, if it is open
            
            % Now check if the figure which monitors the Mahalanobis distances is open.
            % If it is, then also in that figure highlight the trajectories
            % of the brushed units
            
            h=findobj('-depth',1,'Tag','pl_malfwd');
            
            if (~isempty(h))
                
                % make figure which contains monitoring scaled residuals
                % become the current figure
                figure(h);
                
                % Condition || but==0 if but=0 then it is necessary to
                % remove previous highlightments (even if persist='on')
                if strcmp(persist,'off') || but==0
                    % If set of values has already been highlighted in the
                    % mdr plot, remove it
                    a=findobj(h,'Tag','brush_res');
                    delete(a);
                    
                    % Remove the yellow selection in this plot if present
                    a=findobj(h,'Tag','selected');
                    delete(a);
                end
                
                % get the x and y coordinates of the monitoring of the
                % scaled residuals
                a=findobj(h,'tag','data_res');
                abrush=get(a(n+1-nbrush),'ydata');
                if iscell(abrush)
                    ycoord=cell2mat(abrush); % y coordinates of scaled residuals (values)
                else
                    ycoord=abrush; % y coordinate of scaled residual (value)
                end
                xcoord=get(a(1),'Xdata'); % x coordinates of mdr (steps)
                
                hold('on');
                if strcmp('on',persist)
                    plot(gca,xcoord,ycoord,'LineWidth',4,'color',ColorOrd(ij,:),'tag','brush_res');
                else
                    plot(gca,xcoord,ycoord,'LineWidth',4,'color',flagcol,'tag','brush_res');
                end
                hold('off');
                
            end
            
            disp('Brushed steps');
            disp(selxsteps');
            
            disp('Units forming subset at least once in the selected steps of the selected trajectories ');
            disp(nbrush');
            
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
            
            
            
            
            % After waitforbuttonpress:
            % - the standard MATLAB function to be executed on figure
            %   close is recovered
            %      set(gcf,'CloseRequestFcn','closereq');
            
            
            set(gcf,'CloseRequestFcn','closereq');
            Open_yX = findobj(0, 'type', 'figure','tag','pl_spm');
            Open_res = findobj(0, 'type', 'figure','tag','pl_malfwd');
            Open_mdr = findobj(0, 'type', 'figure','tag','pl_mmdrs');
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
end % close options.databrush

    function output_txt = mmdplotLbl(~,event_obj,out)
        %Provides information about the selected point in the trajectory of the mmd
        %(once selected with the mouse)
        %
        % Required input arguments:
        %
        %       obj =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the sense that
        %               these arguments are automatically passed to the function when it executes.
        %     out   =   a structure containing the following fields
        %       Y   =   the response of the regressione model
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
        % REMARK: this function is called by function mmdrsplot
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
        output_txt = {['mmd=',num2str(y,4)]};
        
        % Add information abou the step of the search which is under investigation
        output_txt{end+1} = ['Step m=' num2str(x)];
        
        % If structure out does not contain labels for the rows then
        % labels row1....rown are added automatically
        if isempty(intersect('label',fieldnames(out)))
            out.label=cellstr(num2str((1:n)','row%d'));
        end
        
        % Number of searches associated to a particular value selected
        seltraj=abs(mmdrs(x-init+1,2:end)-y(1))<1e-10;
        
        % Add information about the number of trajectories selected
        % search
        output_txt{end+1} = ['Number of trajectoris='  num2str(sum(seltraj))];
    end

end
%FScategory:VIS-Mult