function covplot(out,varargin)
%covplot plots the trajectories of the elements of the covariance (correlation) matrix monitored
%along the forward search
%
%<a href="matlab: docsearchFS('covplot')">Link to the help function</a>
%
%  Required input arguments:
%  out :  Data to plot. Structure. Structure containing the following fields
%  out.S2cov=   (n-init+1) x (v*(v+1)/2+1) matrix containing the monitoring
%               of the elements of the covariance matrix in each step
%               of the forward search:
%               1st col = fwd search index (from init to n);
%               2nd col = monitoring of S(1,1);
%               3rd col = monitoring of S(1,2);
%               ...;
%               (compulsory argument).
%     out.Un=   matrix containing the order of entry of each unit
%               (necessary if datatooltip is true).
%     out.Y =   n x v data matrix; n observations
%               and v variables.
%      Data Types - struct
%
%
%  Optional input arguments:
%
%       standard:   Appearance of the plot. Structure. Structure which
%                   defines the appearance of the plot in terms of xlim,
%                   ylim, axes labels and their font size style, color of
%                   the lines, etc.
%                   The structure contains the following fields:
%                   standard.SizeAxesNum = scalar specifying the fontsize of the
%                       axes numbers. Default value is 10.
%                   standard.xlim = two elements vector with minimum and maximum of
%                       the x axis. Default value is '' (automatic scale).
%                   standard.ylim = two elements vector with minimum and maximum of
%                       the y axis. Default value is '' (automatic scale).
%                   standard.titl = a label for the title (default: '').
%                   standard.labx = a label for the x-axis (default: 'Subset size m').
%                   standard.laby = a label for the y-axis (default: 'Elements of
%                       the correlation (covariance) matrix)
%                   standard.SizeAxesLab = Scalar specifying the fontsize of the
%                       labels of the axes. Default value is 12.
%                   standard.LineWidth = scalar specifying line width for the
%                       trajectories.
%                   standard.Color = cell array of strings containing the colors to
%                       be used for the standard units.
%                       If length(Color)=1 the same color will be used for
%                       all units.
%                       If length(Color)=2 half of the trajectories will
%                       appear with Color{1} and the other half with
%                       Color{2}. And so on with 3 cell elements, etc.
%                   standard.LineStyle = cell containing the line types.
%                   Example - 'standard',standard
%                   Data Types - structure
%
%         fground : Trajectories in foregroud. Structure.
%                   Structure which defines the trajectories in foregroud,
%                   that is which trajectories are highlighted and how
%                   they are plotted to be distinguishable from the others.
%                   It is possible to control the label, the width, the
%                   color, the line type and the marker of the highlighted
%                   covariances. The structure fground contains the following
%                   fields:
%                   fground.fthresh = (alternative to funit) numeric vector of
%                       length 1 or 2 which specifies the highlighted
%                       trajectories.
%                       If length(fthresh)=1 the highlighted trajectories
%                       are those of units that throughtout the search had
%                       at leat once a covariance greater (in absolute value)
%                       than thresh. The default value of fthresh is Inf.
%                       If length(fthresh)=2 the highlighted trajectories
%                       are those of units that throughtout the search had
%                       a covariance at leat once bigger than fthresh(2) or
%                       smaller than fthresh(1).
%                   fground.funit = (alternative to fthresh) scalar containing the
%                       number of trajectories of the covariances to be
%                       highlighted. For example if funit=5 the
%                       trajectories with the 5 highest values of the
%                       covariances are highlighted.
%                       Notice that if funit is supplied,
%                       fthresh is ignored. The default value of funit is
%                       5.
%                   fground.flabstep = numeric vector which specifies the steps of
%                       the search where to put labels for the highlighted
%                       trajectories (units). The default is to put the
%                       labels at the initial and final steps of the
%                       search, that is fground.flabstep=[m0 n].
%                       flabstep='' means no label.
%                   fground.LineWidth = scalar specifying line width for the
%                       highlighted trajectories (units). Default is 1.
%                   fground.Color = cell array of strings containing the colors to
%                       be used for the highlighted trajectories (units).
%                       If length(scolor)==1 the same color will be used
%                       for all highlighted units. If for example
%                       length(scolor)=2 and there are 6 highlighted units,
%                       3 highlighted trajectories appear with
%                       selunitcolor{1} and 3 highlighted trajectories with
%                       selunitcolor{2}.
%                   fground.LineStyle = cell containing the line type of
%                                       the highlighted trajectories.                      
%                   fground.fmark  = scalar controlling whether to plot highlighted
%                       trajectories as markers. If 1 each trajectory is
%                       plotted using a different marker else (default) no
%                       marker is used.
%                   Example - 'fground',fground
%                   Data Types - structure
%
%       tag     :   Handle of the plot. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag 'pl_cov'. Note that if the
%                   program finds a plot which has a tag equal to the one
%                   specified by the user, then the output of the new plot
%                   overwrites the existing one in the same window else a
%                   new window is created.
%                   Example - 'tag','pl_mycov'
%                   Data Types - string
%
%   datatooltip :   Information about the unit selected. Empty value or structure.
%                   The default is datatooltip=1.
%                   If datatooltip is not empty the user can use the mouse
%                   in order to have information about the unit selected,
%                   the step in which the unit enters the search and the
%                   associated label.
%                   If datatooltip is a structure, it is possible to
%                   control the aspect of the data cursor (see function
%                   datacursormode for more details or the examples below).
%                   The default options of the structure are
%                   DisplayStyle='Window' and SnapToDataVertex='on'.
%                   Example - 'datatooltip',''
%                   Data Types - Empty value or structure.
%
% Output:
%
% See also:
%
% References:
%
% Atkinson, A.C., Riani, M. and Cerioli, A. (2004), "Exploring multivariate
% data with the forward search", Springer Verlag, New York.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('covplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% covplot with all default options.
    % generate input structure for malfwdplot
    n=100;
    p=4;
    state1=141243498;
    randn('state', state1);
    Y=randn(n,p);
    kk=[1:10];
    Y(kk,:)=Y(kk,:)+4;
    [fre]=unibiv(Y);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'init',30);
    % Produce monitoring covariances plot with all the default options
    covplot(out)
%}
%
%{
    % covplot with optional arguments.
    % Example of the use of some options inside structure standard.
    n=100;
    p=4;
    state1=141243498;
    randn('state', state1);
    Y=randn(n,p);
    kk=[1:10];
    Y(kk,:)=Y(kk,:)+4;
    [fre]=unibiv(Y);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'init',30);
    % Initialize structure standard
    standard=struct;
    standard.LineStyle={'-';'-.';'--'};
    % Specify the line width
    standard.LineWidth=0.5;
    covplot(out,'standard',standard)
%}
%
%{
    % Example of the use of some options inside structure fground.
    n=100;
    p=4;
    state1=141243498;
    randn('state', state1);
    Y=randn(n,p);
    kk=[1:10];
    Y(kk,:)=Y(kk,:)+4;
    [fre]=unibiv(Y);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'init',30);
    % Initialize structure fground
    fground = struct;
    % Specify the number of trajectories which have to be highlighted
    fground.funit=2;
    % Specify the steps in which labels have to be put
    n=size(Y,1);
    fground.flabstep=[n/2 n*0.75 n+0.5];
    % Specify the line width of the highlighted trajectories
    fground.LineWidth=3;
    % Produce a monitoring residuals plot in which labels are put for units
    % [2 5 20 23 35 45] in steps [n/2 n*0.75 n+0.5] of the search
    covplot(out,'fground',fground)
%}
%
%{
    %% Example of the use of option tag.
    n=100;
    p=4;
    state1=141243498;
    randn('state', state1);
    Y=randn(n,p);
    kk=[1:10];
    Y(kk,:)=Y(kk,:)+4;
    [fre]=unibiv(Y);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'init',30);
     % Initialize structure fground
    fground = struct;
    % Specify the number of trajectories which have to be highlighted
    fground.funit=2;
    % Specify the steps in which labels have to be put
    n=size(Y,1);
    fground.flabstep=[n/2 n*0.75 n+0.5];
    % Specify the line width of the highlighted trajectories
    fground.LineWidth=3;
   covplot(out,'fground',fground,'tag','pl_mycov')
%}
%

%% Beginning of code 

% Initialization

[n,v]=size(out.Y);

% seq= column vector containing the sequence 1 to n
seq= (1:n)';


m0=out.S2cov(1,1);


S2cov=out.S2cov;
residuals=out.S2cov(:,2:end)';
vv=size(residuals,1);

% the following lines construct the cell array of strings
% 1,2; 1,3; ....  v-1,v;
aco=triu(ones(v,v));
ind = find(abs(aco)>0);
[I,J]=ind2sub(v,ind);
numtext=cellstr([num2str(I) repmat(',',vv,1) num2str(J)]);

% laby= label used for the y-axis of the malfwdplot.
labx='Subset size m';

if max(max(abs(out.S2cov(:,2:end))))<1
    laby='Elements of correlation matrix';
else
    laby='Elements of covariance matrix';
end

LineStyle={'-';  '--'; ':'; '-.'};

% Default options for all trajectories
standarddef = struct(...
    'xlim','','ylim','','titl','','labx',labx,'laby',laby,...
    'Color',{{'b'}},'LineStyle',{LineStyle},...
    'LineWidth',1,'SizeAxesLab',12,'SizeAxesNum',10);

% Default options for the trajectories in foreground
fgrounddef = struct('fthresh','','flabstep',[m0 n],'funit',5,...
    'fmark',0,'LineWidth','','Color','','LineStyle','');

options=struct('standard',standarddef,'fground',fgrounddef,'tag','pl_cov','datatooltip',1);


%% Preliminary checks

if nargin<1
    error('FSDA:covplot:missingInputs','A required input argument is missing.')
end

%get optional user options
if nargin>1
    UserOptions=varargin(1:2:length(varargin));
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:covplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end


%% Prepate the figure to display the covplot

% Create a figure to host the plot or clear the existing one
h=findobj('-depth',1,'tag',options.tag);
if (~isempty(h))
    clf(h);
    figure(h);
    axes;
else
    h=figure;
end
set(h,'Name', 'Monitoring of cov elements', 'NumberTitle', 'off');
hold('all');


%% standard options
% get the option names: structure option was initialised with standarddef
% and updated with optional user's options. For the options not set by the
% user, use their default value


if ~isequal(options.standard,standarddef)
    fld=fieldnames(options.standard);
    
    % Check if user options inside options.fground are valid options
    chkoptions(standarddef,fld)
    for i=1:length(fld)
        standarddef.(fld{i})=options.standard.(fld{i});
    end
end

standard=standarddef;


plot1=plot(S2cov(:,1),S2cov(:,2:end),'LineWidth',standard.LineWidth);

% Apply color
scol=standard.Color;

if ~isempty(scol)
    if size(scol,2)>1
        scol=scol';
    end
    scol=repmat(scol,ceil(vv/length(scol)),1);
    set(plot1,{'Color'},scol(1:vv));
end

% Apply Line Style
slintyp=standard.LineStyle;
if ~isempty(slintyp)
    if size(slintyp,2)>1
        slintyp=slintyp';
    end
    slintyp=repmat(slintyp,ceil(vv/length(slintyp)),1);
    set(plot1,{'LineStyle'},slintyp(1:vv));
end

% control minimum and maximum for x and y axis
if ~isempty(standard.xlim)
    xlim(standard.xlim);
end
if ~isempty(standard.ylim)
    ylim(standard.ylim);
end

% Main title of the plot and labels for the axes
labx=standard.labx;
laby=standard.laby;
titl=standard.titl;
title(titl);

% Add the x and y labels to the plot.
SizeAxesLab=standard.SizeAxesLab;
xlabel(labx,'Fontsize',SizeAxesLab);
ylabel(laby,'Fontsize',SizeAxesLab);

% FontSizeAxes = font size for the axes numbers
SizeAxesNum=standard.SizeAxesNum;
% Specify the FontSize of the number on the axes
set(gca,'FontSize',SizeAxesNum)

% displays the boundary of the current axes.
box on

%% fground options
if ~isempty(options.fground)
    
    % Control the appearance of the trajectories to be highlighted
    if ~isequal(options.fground,fgrounddef)
        
        fld=fieldnames(options.fground);
        
        % Check if user options inside options.fground are valid options
        chkoptions(fgrounddef,fld)
        for i=1:length(fld)
            fgrounddef.(fld{i})=options.fground.(fld{i});
        end
    end
    
    % For the options not set by the user use their default value
    fground=fgrounddef;
    
    % fground.flabstep option and check if the choice of flabsteps is valid
    if ~isempty(fground.flabstep)
        steps=floor(fground.flabstep);
        if max(steps)>n || min(steps)<m0
            mess=sprintf(['Warning: steps that you have chosen outside the range [m0 ... n]\n',...
                'are not considered']);
            fprintf('%s\n',mess);
            steps=steps(steps>=m0 & steps<=n);
        end
    end
    
    
    % fthresh= threshold to define units which have to be displayed in foreground
    % (highlighted)
    fthresh=fground.fthresh;
    % funit= List of the units to be displayed in foreground (highlighted)
    funit=fground.funit;
    
    if ~isempty(funit)
        selmaxabs=max(abs(residuals),[],2);
        [~,IX]=sort(selmaxabs,'descend');
        
        % Some checks on minimum and maximum of vector funit
        if max(funit)>vv || min(funit)<1
            sprintf(['Warning: the number of trajectories to label',...
                'must be >1 and smaller than ' num2str(vv) ]);
        end
        
        funit=IX(1:min(round(funit),length(IX)));
    else
        selmax=max(residuals,[],2);
        selmin=min(residuals,[],2);
        if length(fthresh)>1
            funit=seq(selmax>fthresh(2) | selmin<fthresh(1));
        else
            funit=seq(selmax>fthresh | selmin<-fthresh);
        end
    end
    
    
    % lunits = number of units which must be highlighted
    lunits=length(funit);
    
    % Specify the line type for the highlighted units
    % (those forming vector funit)
    slintyp=fground.LineStyle;
    
    if ~isempty(slintyp)
        slintyp=repmat(slintyp,ceil(n/length(slintyp)),1);
        set(plot1(funit),{'Line'},slintyp(funit));
    end
    
    if ~isempty(fground.flabstep)
        % lsteps = number of steps for which it is necessary to add the labels
        lsteps=length(steps);
        lall=lunits*lsteps;
        
        text(reshape(repmat(steps,lunits,1),lall,1),reshape(residuals(funit,steps-m0+1),lall,1),reshape(repmat(numtext(funit),1,lsteps),lall,1))
    end
    
    % if requested, set the color of the selected trajectories note that if
    % scolor contains more than one color, e.g. options.scolor = {'b';'g';'r'},
    % then the colors of the trajectories alternate.
    fcol=fground.Color;
    if ~isempty(fcol)
        fcol=repmat(fcol,ceil(lunits/length(fcol)),1);
        set(plot1(funit),{'Color'},fcol(1:lunits));
    end
    
    % if requested, set the selected trajectories in LineWidth
    if isnumeric(fground.LineWidth)
        set(plot1(funit),'LineWidth',fground.LineWidth);
    end
    
    % If requested, add markers to all the trajectories
    if fground.fmark==1
        styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
        styp=repmat(styp,ceil(n/13),1);
        set(plot1(funit),{'Marker'},styp(funit))
    end
    
end

hold('off')
% include specified tag in the current plot
set(gcf,'tag',options.tag)
set(gca,'Position',[0.1 0.1 0.85 0.85])

Un=out.Un;

%% Set the datatooltip for the covplot
if ~isempty(options.datatooltip)
    hTarget=[];
    hTargetlwd=[];
    hTargetcol=[];
    try
        % chkgpu=gpuDeviceCount; 
        % datacursormode on;
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
        set(hdt,'UpdateFcn',{@covplotLbl,out})
    catch
        disp('No graphical device, interactive datatooltip not enabled')
    end
end


    function output_txt = covplotLbl(~,event_obj,~)
        %% covplotLbl provides information about the selected covariances (correlations)
        %
        % Required input arguments:
        %
        %   obj     =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object
        %               (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the
        %               sense that these arguments are automatically passed
        %               to the function when it executes.
        %      label=  (optional argument) if it is present it must be
        %               a cell array of strings containing the labels of
        %               the rows of the regression dataset
        %
        % Output:
        %
        %   output_txt=  Datatip text (string or string cell array) which
        %                informs
        %                about the step of the search which has been
        %                selected, the unit(s) selected and its (their)
        %                entry during the fwd search
        %
        %
        % References:
        %
        %   Atkinson and Riani (2000), Robust Diagnostic Regression
        %   Analysis, Springer Verlag, New York.
        %
        % Written by FSDA team
        
        if ~isempty(hTarget)
            % set old line width and old color for old selection
            set(hTarget,'LineWidth',hTargetlwd,'Color',hTargetcol);
        else
        end
        
        % Store line width and color of selected trajectory
        % Notice that changing event_obj.Target in subsequent lines seems
        % to affect also hTarget
        hTarget=event_obj.Target;
        hTargetlwd=get(hTarget,'LineWidth');
        hTargetcol=get(hTarget,'Color');
        
        % Increase Line width and set color to red (or to blue if previous
        % color was red) of selected trajectory
        if sum(get(hTarget,'Color')==[1 0 0])==3
            set(hTarget,'LineWidth',hTargetlwd+1.5,'Color','b');
        else
            set(hTarget,'LineWidth',hTargetlwd+1.5,'Color','r');
        end
        
        
        pos = get(event_obj,'Position');
        
        % x and y, plot coordinates of the mouse
        x1 = pos(1); y1 = pos(2);
        
        % Find index to retrieve obs. name Consider that find return the
        % linear indexing of matrix xydata
        % residuals=out.RES;
        idx = find(residuals == y1,1);
        % Linear indexing is transformed into normal indexing using
        % function ind2sub row and column contain the column and row
        % indexed of the observation which has been selected with the mouse
        [row,~] = ind2sub(size(residuals),idx);
        
        covij=numtext{row};
        
        
        if isempty(covij)
            output_txt{1}=['no covariance has coordinates x,y' num2str(x1) '' num2str(y1)] ;
        else
            
            % If structure out does not contain labels for the rows then
            % labels row1....rown are added automatically
            
            output_txt=cell(3,1);
            
            % output_txt is what it is shown on the screen
            output_txt{1,1} = ['Cov(' covij ')=' num2str(y1,4)];
            
            % Add information about the step of the search which is under
            % investigation
            output_txt{2,1} = ['Step m=' num2str(x1)];
            
            
            % Add information about the units entering at the selcted step
            % of the search
            idx = find(Un(:,1)==x1,1);
            sel=Un(idx,2:end);
            
            output_txt{3,1} = ['Unit(s) entered in step ' num2str(x1) '='  num2str(sel(~isnan(sel)))];
            
            
        end
    end


end

%FScategory:VIS-Mult
