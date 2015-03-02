function mmdplot(out,varargin)
%mmdplot plots the trajectory of minimum mhalanobis distance (mmd)
%
%<a href="matlab: docsearchFS('mmdplot')">Link to the help function</a>
%
% Required input arguments:
%
%  out :  structure containing the following fields
%        mmd =  a matrix containing the monitoring of minimum Mahalanobis 
%               distance in each step of the forward search. The first
%               column of mdr must contain the fwd search index This matrix
%               can be created using function FSMeda (compulsory argument).
%               If this matrix has three columns in the third colum there
%               is the monitoring of the (m+1)-th Mahalanobis distance
%       Un  =   matrix containing the order of entry of each unit
%               (necessary if datatooltip is true or databrush is not empty)
%       Y   =   the original n x v data matrix (necessary only if
%               option databrush is not empty)
%
%
% Optional input arguments:
%
%       quant   :   vector containing quantiles for which envelopes have
%                   to be computed. The default is to produce 1%, 50% and
%                   99% envelopes. In other words the default is
%                   quant=[0.01;0.5;0.99];
%       exact:      scalar, if it is equal to 1 the calculation of the
%                   quantiles of the F distribution is based on
%                   functions finv and tinv from the Matlab statistics
%                   toolbox, otherwise the calculations of the former
%                   quantiles is based on function invcdff invcdft.
%                   The solution has a tolerance of 1e-8 (change variable
%                   tol in file invcdff.m required.
%                   Remark: the use of function finv is more
%                   precise but requires more time. The default value of
%                   exact is 1 (exact solution).
%       mplus1  :   Scalar, if mplus1=1 it is also possible to plot the
%                   curve associated with (m+1)th order statistic.
%                   The default is mplus1=0
%       envm    :   Scalar which specifies the size of the sample which is
%                   used to superimpose the envelope. The default is to add
%                   an envelope based on all the observations (size n
%                   envelope)
%       xlimx   :   vector with two elements controlling minimum and
%                   maximum on the x axis. Default value is mdr(1,1)-3 and
%                   mdr(end,1)*1.3
%       ylimy   :   vector with two elements controlling minimum and
%                   maximum on the y axis. Default value is min(mdr(:,2))
%                   and max(mdr(:,2));
%       lwdenv  :   Scalar which controls the width of the lines associated
%                   with the envelopes. Default is lwdenv=1
%       tag     :   string which identifies the handle of the plot which
%                   is about to be created. The default is to use tag
%                   'pl_mmd'. Notice that if the program finds a plot which
%                   has a tag equal to the one specified by the user, then
%                   the output of the new plot overwrites the existing one
%                   in the same window else a new window is created
%   datatooltip :   empty value or structure. The default is datatooltip=''
%                   If datatooltip is not empty the user can use the mouse
%                   in order to have information about the unit seected,
%                   the step in which the unit enters the search and the
%                   associated label. If datatooltip is a structure, it is
%                   possible to control the aspect of the data cursor (see
%                   function datacursormode for more details or the
%                   examples below). The default options of the structure
%                   are DisplayStyle='Window' and SnapToDataVertex='on'
%       label   :   cell containing the labels of the units (optional
%                   argument used when datatooltip=1. If this field is not
%                   present labels row1, ..., rown will be automatically
%                   created and included in the pop up datatooltip window)
%    databrush :    empty value, scalar or structure.
%                   DATABRUSH IS AN EMPTY VALUE If databrush is an empty
%                   value (default), no brushing is done. The activation of
%                   this option (databrush is a scalar or a structure) enables
%                   the user  to select a set of trajectories in the
%                   current plot and to see them highlighted in the spm
%                   (notice that if the spm does not exist it is
%                   automatically created). In addition, brushed units can
%                   be highlighted in the monitoring MD plot
%                   Remark: the window style of the
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
%                  labeladd. If this option is '1', we label the units
%                     of the last selected group with the unit row index in
%                     matrices X and y. The default value is labeladd='',
%                     i.e. no label is added.
%       Fontsize:   Scalar which controls the fontsize of the labels of the
%                   axes. Default value is 12
%    SizeAxesNum:   Scalar which controls the fontsize of the numbers of
%                   the axes. Default value is 10
%       nameY   :   cell array of strings of length p containing the labels
%                   of the varibles of the regression dataset. If it is empty
%                 	(default) the sequence Y1, ..., Yp will be created
%                   automatically
%       lwd     :   Scalar which controls linewidth of the curve which
%                   contains the monitoring of minimum Mahalanobis distance.
%                   Default line width=2
%       titl    :   a label for the title (default: '')
%       labx    :   a label for the x-axis (default: 'Subset size m')
%       laby    :   a label for the y-axis (default: 'Minimum Mahalnobis distance')
%
%
% See also:
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
%<a href="matlab: docsearchFS('mmdplot')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    %Steps common to all the examples
    Y=load('head.txt');
    [fre]=unibiv(Y);
    %create an initial subset with the 4 observations, which fell the smallest
    %number of times outside the robust bivariate ellipses, and with the
    %lowest Mahalanobis Distance.
    fre=sortrows(fre,[3 4]);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs);
    mmdplot(out)
%}

%{
    %Example of the use of function mdrplot with personalized envelopes
    mmdplot(out,'quant',[0.99;0.9999]);
%}

%{
    % Interactive_example
    %Example of the use of function mmdplot with datatooltip passed as
    %scalar (that is using default options for datacursor (i.e.
    %DisplayStyle =window)
     mmdplot(out,'datatooltip',1);
%}

%{
    % Interactive_example
    %Example of the use of function mmdplot with datatooltip passed as
    %structure

    clear tooltip
    tooltip.SnapToDataVertex='on'
    tooltip.DisplayStyle='datatip'
    mmdplot(out,'datatooltip',tooltip);
%}


%{
   %Example of the use of option envm
   %In this case the resuperimposed envelope is based on n-2 observations
   mmdplot(out,'envm',size(Y,1)-2);

%}

%{
    % Interactive_example
    %Example of the use of function mdrplot with databrush
    mmdplot(out,'databrush',1);
%}

%{
    % Interactive_example
    %Example where databrush is a structure
    databrush=struct
    databrush.selectionmode='Lasso'
    mmdplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    %Example of the use of brush using brush mode
    databrush=struct
    databrush.selectionmode='Brush'
    databrush.Label='on';
    mmdplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    %Example of the use of persistent non cumulative brush. Every time a
    %brushing action is performed previous highlightments are removed
    databrush=struct
    databrush.persist='off'
    databrush.Label='on';
    databrush.RemoveLabels='off';
    mmdplot(out,'databrush',databrush)
%}

%{
    % Interactive_example
    %Example of the use of persistent cumulative brush. Every time a
    %brushing action is performed current highlightments are added to
    %previous highlightments
    databrush=struct
    databrush.persist='on';
    databrush.selectionmode='Rect'
    mmdplot(out,'databrush',databrush)
%}


%% Initialization

% Extract the absolute value of minimum deletion residual
% or minimum Mahalanobis distance

mdr=abs(out.mmd);

[n,v]=size(out.Y);

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
laby='Minimum Mahalanobis distance';


%% User options

options=struct('quant', quant,'exact',0,'sign',0,'mplus1',0,...
    'envm',n,'xlimx',xlimx,'ylimy',ylimy,'lwd',2,'lwdenv',1,...
    'FontSize',12,'SizeAxesNum',10,'tag','pl_mmd',...
    'datatooltip','','databrush','',...
    'titl','','labx',labx,'laby',laby,'nameY','','label','','scaled',0);

if nargin<1
    error('FSDA:mmdplot:missingInputs','A required input argument is missing.')
end

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mmdplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
exact=options.exact;
mplus1=options.mplus1;

% scaled = scaled or unscaled envelopes (default is unscaled envelopes)
scaled=options.scaled;

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
if init-v==0;
    init=init+1;
end

%% Display the min MD plot

%Create a figure to host the plot or clear the existing one
h=findobj('-depth',1,'tag',options.tag);
if (~isempty(h))
    clf(h);
    figure(h(end))
    axes;
else
    % create a new figure
    h=figure;
end
set(h,'Name', 'Monitoring of Minimum Mahalnobis distance', 'NumberTitle', 'off');
hold('all');

% Theoretical envelopes for minimum Mahalnobis distance
% [gmin] = FSRenvmdr(envm,p,'prob',quant,'init',init,'exact',exact);
[gmin] = FSMenvmmd(envm,v,'prob',quant,'init',init,'exact',exact,'scaled',scaled);

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

% plot minimum MD
stdColor = 'b'; stdLineStyle = '-';
plot(mdr(sel,1),mdr(sel,2),'tag','data_mmd',...
    'Color',stdColor,'LineStyle',stdLineStyle,'LineWidth',lwd);



% Write an extra message on the plot
annotation(gcf,'textbox',[0.2 0.8 0.1 0.1],'EdgeColor','none','String',['Envelope based on ' num2str(envm) ' obs.'],'FontSize',FontSize);

% If mplus1=1 add the line associated with (m+1) ordered MD
if mplus1;
    line(mdr(sel,1),mdr(sel,3),'LineWidth',lwd,'LineStyle',':','Color','cyan','tag','env');
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

hold('off');

% displays the boundary of the current axes.
box on

% include specified tag in the current plot
set(gcf,'tag',options.tag);

% Store the handle of the mdrplot inside handle hmin
hmin=gcf;

Un=out.Un;

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
        
        % If option Label is 'on' then matrix Un is added to UserData
        d=max(strcmp('Label',fieldnames(databrush)));
        if d==1 && strcmp(databrush.Label,'on')
            set(gcf,'UserData',Un)
        end
        
        
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
    but=0; brushcum=[]; ij=1;
    
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
            
            % Now check if the figure which monitors the residuals is open.
            % If it is, then also in that figure highlight the trajectories
            % of the brushed units
            
            h=findobj('-depth',1,'Tag','pl_malfwd');
            
            if (~isempty(h))
                
                % make figure which contains monitoring scaled residuals
                % become the current figure
                figure(h);
                
                % Condition || but==0 if but=0 then it is necessary to
                % remove previous highlightments (even if persist='on')
                if strcmp(persist,'off') || but==0;
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
            % After waitforbuttonpress:
            % - the standard MATLAB function to be executed on figure
            %   close is recovered
            set(gcf,'CloseRequestFcn','closereq');
            % - and the 'but' variable is set if keyboard key was pressed
            if ss==1;
                but=2;
            end
        else
            but=2;
        end
        
    end % close loop associated with but
end % close options.databrush

    function output_txt = mmdplotLbl(~,event_obj,out)
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
        output_txt = {['mmd=',num2str(y,4)]};
        
        % Add information abou the step of the search which is under investigation
        output_txt{end+1} = ['Step m=' num2str(x)];
        
        % If structure out does not contain labels for the rows then
        % labels row1....rown are added automatically
        if isempty(intersect('label',fieldnames(out)))
            out.label=cellstr(num2str((1:n)','row%d'));
        end
        
        % Add information about the next step in which the selected unit entered the
        % search
        idx = find(Un(:,1)==x,1)+1;
        sel=Un(idx,2:end);
        
        output_txt{end+1} = ['Unit(s) entered in step ' num2str(x+1) '='  num2str(sel(~isnan(sel)))];
        
    end

end

