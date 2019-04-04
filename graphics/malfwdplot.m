function plotopt=malfwdplot(out,varargin)
%malfwdplot plots the trajectories of scaled Mahalanobis distances along the search
%
%<a href="matlab: docsearchFS('malfwdplot')">Link to the help function</a>
%
%  Required input arguments:
%
%  out :  Structure containing monitoring of Mahalanobis distance. Structure.
%               Structure containing the following fields.
%       out.MAL =   matrix containing the squared Mahalanobis distances monitored in each
%               step of the forward search. Every row is associated with a
%               unit (row of data matrix Y).
%               This matrix can be created using function FSMeda
%               (compulsory argument)
%       out.Un  =   matrix containing the order of entry of each unit
%               (necessary if datatooltip is true or databrush is not empty)
%       out.Y   =   n x v data matrix; n observations
%               and v variables
%                Data Types - single|double
%
%
%  Optional input arguments:
%           standard : plot layout. Structure. Appearance of the plot
%                   in terms of xlim, ylim, axes labels and their font size
%                   style, color of the lines, etc. 
%                   Structure standard contains the following fields:
%                   standard.SizeAxesNum = scalar specifying the fontsize of the
%                       axes numbers. Default value is 10.
%                   standard.xlim = two elements vector with minimum and maximum of
%                       the x axis. Default value is '' (automatic scale).
%                   standard.ylim = two elements vector with minimum and maximum of
%                       the y axis. Default value is '' (automatic scale).
%                   standard.titl = a label for the title (default: '').
%                   standard.labx = a label for the x-axis (default: 'Subset size m').
%                   standard.laby = a label for the y-axis (default: 'Mahalanobis distances'
%                       or 'Scaled Mahalanobis distances').
%                   standard.SizeAxesLab = Scalar specifying the fontsize of the
%                       labels of the axes. Default value is 12.
%                   standard.subsize = numeric vector containing the subset size
%                       with length equal to the number of columns of
%                       matrix of MD. The default value of subsize is
%                       size(MAL,1)-size(MAL,2)+1:size(MAL,1)
%                   standard.LineWidth = scalar specifying line width for the
%                       trajectories.
%                   standard.Color= cell array of strings containing the colors to
%                       be used for the highlighted units.
%                       If length(Color)=1 the same color will be used for
%                       all units.
%                       If length(Color)=2 half of the trajectories will
%                       appear with Color{1} and the other half with
%                       Color{2}. And so on with 3 cell elements, etc.
%                   standard.LineStyle = cell containing the line types.
%
%                   Remark. The default values of structure standard are:
%                   standard.SizeAxesNum=10;
%                   standard.SizeAxesLab=12;
%                   standard.xlim='' (Automatic scale);
%                   standard.ylim='' (Automatic scale);
%                   standard.titl='' (empty title);
%                   standard.labx='Subset size m';
%                   standard.laby='Mahalanobis distances';
%                   standard.LineWidth=1;
%                   standard.Color={'b'};
%                   standard.LineStyle={'-'}
%
%                   Example - 'standard.LineWidth','1'
%                   Data Types - struct
%         fground : trajectories in foregroud.
%                   Structure. Structure which controls which trajectories
%                   are highlighted and how they are plotted to be
%                   distinguishable from the others.
%                   It is possible to control the label, the width, the
%                   color, the line type and the marker of the highlighted
%                   units.
%                   Structure fground contains the following fields:
%                   fground.fthresh = (alternative to funit) numeric vector
%                       of length 1 or 2 which specifies the criterion to
%                       select the trajectories to highlight.
%                       If length(fthresh)=1 the highlighted trajectories
%                       are those units that throughtout the search had
%                       at leat once a MD greater (in absolute value)
%                       than fthresh. The default value of fthresh is 2.5.
%                       If length(fthresh)=2 the highlighted trajectories
%                       are those of units that throughtout the search had
%                       a MD at least once bigger than fthresh(2) or
%                       smaller than fthresh(1).
%                   fground.funit = (alternative to fthresh) vector
%                       containing the list of the units to be highlighted.
%                       If funit is supplied, fthresh is ignored.
%                   fground.flabstep = numeric vector which specifies the steps of
%                       the search whre to put labels for the highlighted
%                       trajectories (units). The default is to put the
%                       labels at the initial and final steps of the search.
%                       flabstep='' means no label.
%                   fground.LineWidth = scalar specifying line width for the
%                       highlighted trajectories (units). Default is 1.
%                   fground.Color = cell array of strings containing the colors to
%                       be used for the highlighted trajectories (units).
%                       If length(Color)==1 the same color will be used for
%                       all highlighted units Remark: if for example
%                       length(Color)=2 and there are 6 highlighted units,
%                       3 highlighted trajectories appear with
%                       Color{1} and 3 highlighted trajectories with
%                       Color{2}.
%                   fground.LineStyle = cell containing the line type of
%                       the highlighted trajectories.
%                   fground.fmark  = scalar controlling whether to plot
%                       highlighted trajectories as markers. if 1 each line
%                       is plotted using a different marker else no marker
%                       is used (default).
%                   fground.FontSize  = scalar controlling font size of the
%                       labels of the trajectories in foreground
%
%                   Remark. The default values of structure fground are:
%                    fground.fthresh=2.5;
%                    fground.flabstep=[m0 n];
%                    fground.LineWidth=1;
%                    fground.LineStyle={'-'};
%                    fground.FontSize=12;
%                   Remark. if fground='' no unit is highlighted and no
%                   label is inserted into the plot.
%
%                   Example - 'fground.LineWidth','1'
%                   Data Types - struct
%         bground : characterictics of the trajectories in background.
%                   Structure.
%                    Structure which specifies the trajectories in background,
%                   i.e. the trajectories corresponding to "unimmportant"
%                   units in the central part of the data. The structure
%                   also specifies the style used in the plot to give them
%                   less emphasis, so that to not distract the eye of the
%                   analyst from the trajectories of the relevant units.
%                   Structure bground contains the following fields:
%                   bground.bthresh = numeric vector of length 1 or 2 which
%                       specifies how to define the unimmportant trajectories.
%                       Unimmportant trajectories will be plotted using a
%                       colormap, in greysh or will be hidden.
%                       - If length(bthresh)=1 the irrelevant units are
%                         those which always had a MD smaller
%                         (in absolute value) than thresh.
%                       - If length(bthresh)=2 the irrelevant units
%                         are those which always had a MD greater
%                         than bthresh(1) and smaller than bthresh(2).
%                       The default is:
%                           bthresh=2.5 if n>100 and bthresh=-inf if n<=100
%                       i.e. to treat all trajectories as important if
%                       n<=100 and, if n>100, to reduce emphasis only to
%                       trajectories having in all steps of the search a
%                       value of scaled MD smaller than 2.5.
%                   bground.bstyle = specifies how to plot the unimportant
%                       trajectories as defined in option bthresh.
%                       'faint': unimportant trajectories are plotted using
%                           a colormap.
%                       'hide': unimportant trajectories are hidden.
%                       'greysh': unimportant trajectories are displayed in
%                           a faint grey.
%                       When n>100 the default option is 'faint'.
%                       When n<=100 and bthresh = -Inf option bstyle is
%                       ignored.
%                   Example - 'bground.bstyle','faint'
%                   Data Types - struct
%                   Remark: bground='' is equivalent to bground.bthresh=-Inf
%                   that is all trajectories are considered relevant.
%
%       tag     :   Personalized plot tag. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag 'pl_mal'.
%                   Note that if the program finds a plot which has a tag
%                   equal to the one specified by the user, then the output
%                   of the new plot overwrites the existing one in the same
%                   window else a new window is created.
%                   Example - 'tag','myplot'
%                   Data Types - char
%   datatooltip :   interactive clicking. Empty value (default) or
%                   structure. The default is datatooltip=''.
%                   If datatooltip = 1, the user can select with the
%                   mouse an individual MD trajectory in order to
%                   have information about the corresponding unit, the
%                   associated label and the step of the search in which
%                   the unit enters the subset.
%                   If datatooltip is a structure it may contain the
%                   following the fields
%                   datatooltip.DisplayStyle = Determines how the data
%                       cursor displays. datatip | window.
%                       - datatip displays data cursor
%                       information in a small yellow text box attached to a
%                       black square marker at a data point you interactively
%                       select.
%                       - window displays data cursor information for the
%                       data point you interactively select in a floating
%                       window within the figure.
%                   datatooltip.SnapToDataVertex=  Specifies whether the
%                       data cursor snaps to the nearest data value or is
%                       located at the actual pointer position.  on | off.
%                       - on data cursor snaps to the nearest data value
%                       - off data cursor is located at the actual pointer
%                       position.
%                   (see the MATLAB function datacursormode or the examples
%                   below). Default values are datatooltip.DisplayStyle =
%                   'Window' and datatooltip.SnapToDataVertex = 'on'.
%                   datatooltip.LineColor = controls the line color of the
%                       selected trajectory. Vector with three elements
%                       specifying RGB color.
%                       By default, the MD trajectory selected with the
%                       mouse is highlighted in red. (a RGB vector can be
%                       conveniently chosen with our MATLAB class FSColor,
%                       see documentation).
%                       Example -
%                       datatooltip.LineColor=[155/255,190/255,61/255];
%                       or using class FScolors:
%                       datatooltip.LineColor='greenish'
%                   datatooltip.SubsetLinesColor = highlights the
%                       trajectories of the units
%                       that are in the subset at a given step of the
%                       search. Vector with three elements
%                       specifying RGB color.
%                       Example - datatooltip.SubsetLinesColor=[0 0 0];
%                       or using class FScolors:
%                       datatooltip.SubsetLinesColor='black'; highlights in
%                       black the trajectories of the units which are
%                       inside subset in correpondence of the selected
%                       steps.
%                   Remark. This can be done (repeatedly) with a left mouse click
%                       on the x axis ('subset size m') in proximity of the
%                       step of interest. A right mouse click will
%                       terminate the selection by marking with a up-arrow
%                       the step corresponding to the highlighted lines.
%                       The highlighted lines by default are in blue, but
%                       different colors can be specified as RGB vectors in
%                       the field SubsetLinesColor. By default
%                       SubsetLinesColor = '', i.e. the modality is not
%                       active. Any initialization for SubsetLinesColor
%                       which cannot be interpreted as RGB vector will be
%                       converted to blue, i.e. SubsetLinesColor will be
%                       forced to be [0 0 1].
%                   Example - 'datatooltip',''
%                   Data Types - empty value, scalar or struct
%       label   :   row labels. Cell of strings. Cell containing the labels
%                   of the units (optional argument used when
%                   datatooltip=1. If this field is not present labels
%                   row1, ..., rown will be automatically created and
%                   included in the pop up datatooltip window).
%                   Example - 'label',{'Smith','Johnson','Robert','Stallone'}
%                   Data Types - cell
%    databrush  :   interactive mouse brushing. Empty value, scalar or structure.
%                   If databrush is an empty value (default), no brushing
%                   is done.
%                   The activation of this option (databrush is a scalar or
%                   a structure) enables the user  to select a set of
%                   trajectories in the current plot and to see them
%                   highlighted in the scatter plot matrix (spm).
%                   If spm does not exist it is automatically created.
%                   In addition, brushed units are automatically highlighted in the
%                   minimum MD plot if it is already open.
%                   Please, note that the window style of the other figures is set
%                   equal to that which contains the monitoring MD
%                   plot. In other words, if the MD plot
%                   is docked all the other figures will be docked too.
%                   DATABRUSH IS A SCALAR.
%                   If databrush is a scalar the default selection tool is
%                   a rectangular brush and it is possible to brush only
%                   once (that is persist='').
%                   DATABRUSH IS A STRUCTURE.
%                   If databrush is a structure, it is possible to use all
%                   optional arguments of function selectdataFS.m and the
%                   following fields:
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
%                   - databrush.Label = add labels of brushed units in the
%                     malfwdplot.
%                   - databrush.labeladd = add labels of brushed units in spm.
%                     Character. [] (default) | '1'.
%                     If databrush.labeladd='1', we label in the scatter
%                     plot matrix the units of the last selected group with
%                     the unit row index in matrices X and y. The default
%                     value is labeladd='', i.e. no label is added.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%       nameY   :   variable labels. Cell array of strings.
%                   Cell array of strings of length v containing the labels
%                   of the variables of the dataset. Cell of strings. If it
%                   is empty (default) the sequence Y1, ..., Yv will be used
%                   automatically
%                   Example - 'nameY',{'var1', var2, 'var3'}
%                   Data Types - cell of strings
%       msg     :   display or save used options. Scalar. Scalar which
%                   controls whether to display or to save
%                   as output the options inside structures standard,
%                   fground and bground which have been used to draw the
%                   plot.
%                   plotopt=malfwdplot(out,'msg',1) enables to save inside
%                   cell  plotopt the options which have been used to draw
%                   the three types of trajectories (standard, foreground
%                   and background)
%                   plotopt=malfwdplot(out,'msg',2) saves inside cell plotopt
%                   the options which have been used and prints them on the
%                   screen
%                   Example - 'msg',1
%                   Data Types - single or double
%        conflev :  confidence interval for the horizontal bands. Vector.
%                   It can be a vector of different confidence level values,
%                   e.g. [0.95,0.99,0.999]. Confidence interval is based on
%                   the chi^2 distribution.
%                   Example - 'conflev',0.99
%                   Data Types - single or double
%
% Output:
%
%   plotopt : options which have been used to create the plot. Cell array
%               of strings. Store all options which have been used to
%               generate the plot inside cell plotopt.
%
%
% See also: resfwdplot
%
% References:
%
% Atkinson, A.C., Riani, M. and Cerioli, A. (2004), "Exploring multivariate
% data with the forward search", Springer Verlag, New York.
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('malfwdplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Produce monitoring MD plot with all the default options.
    % Generate input structure for malfwdplot
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
    [out]=FSMeda(Y,bs,'plots',1,'init',30);
    % Produce monitoring MD plot with all the default options
    malfwdplot(out)
%}
%
%{
    %% Example of the use of some options inside structure standard.
    % Initialize structure standard
    standard=struct;
    standard.LineStyle={'-';'-.';':'};
    % Specify the line width
    standard.LineWidth=0.5;
    malfwdplot(out,'standard',standard)
%}
%
%{
    % Example of the use of some options inside structure fground.
    % Initialize structure fground
    fground = struct;
    % Specify which trajectories have to be highlighted
    fground.funit=[2 5 20 23 35 45];
    % Specify the steps in which labels have to be put
    n=size(Y,1);
    fground.flabstep=[n/2 n*0.75 n+0.5];;
    % Specify the line width of the highlighted trajectories
    fground.LineWidth=3;
    % Produce a monitoring MD plot in which labels are put for units
    % [2 5 20 23 35 45] in steps [n/2 n*0.75 n+0.5] of the search
    malfwdplot(out,'fground',fground)
%}
%
%{
    % Example of the use of some options inside structure bground.
    bground = struct;
    % Specify a threshold to define the "background" trajectories
    bground.bthresh=2;
    % Trajectories whose MD is always between -btresh and +bthresh
    % are shown as specified in bground.bstyle
    bground.bstyle='faint';
    malfwdplot(out,'bground',bground)
%}
%
%
%{
    % Interactive_example
    %   Example of the use of option databrush.
    %   (brushing is done only once using a rectangular selection tool)
    malfwdplot(out,'databrush',1)
    %   An equivalent statement is
    databrush=struct;
    databrush.selectionmode='Rect';
    malfwdplot(out,'databrush',databrush)
%}

%{
    % Example of the use of some options inside structure fground.
    % load Swiss banknotes
    Y=load('swiss_banknotes.txt');
    [fre]=unibiv(Y);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'init',30);

    % Initialize structure fground
    fground = struct;
    % Specify which trajectories have to be highlighted
    fground.funit=out.Un(end-15:end,2);
    % Specify the steps in which labels have to be put
    n=size(Y,1);
    fground.flabstep=[n/2 n*0.75 n+0.5];;
    % Specify the line width of the highlighted trajectories
    fground.LineWidth=3;
    % Produce a monitoring MD plot in which labels are put for units
    % out.Un(end-15:end,2)in steps [n/2 n*0.75 n+0.5] of the search
    malfwdplot(out,'fground',fground)
%}
%
%{
    %   Example of the use of option datatooltip.
    %   Gives the user the possibility of clicking on the different points
    %   and have information about the unit selected, the step of entry
    %   into the subset and the associated label
    malfwdplot(out,'datatooltip',1);
%}

%{
    % Example of the use of option datatooltip personalized.
    % Gives the user the possibility of clicking on the different points
    % and have information about the unit selected, the step of entry
    % into the subset and the associated label.
    datatooltip = struct;
    % In this example the style of the datatooltip is 'datatip'. Click on a
    % trajectory when the resfwdplot is displayed.
    %
    datatooltip.DisplayStyle = 'datatip';
    malfwdplot(out,'datatooltip',datatooltip);
    %
    % Now we use the default style, which is 'window'.
    datatooltip.DisplayStyle = 'window';
    malfwdplot(out,'datatooltip',datatooltip);

    % Here we specify the RGB color used to highlight the selected trajectory.
    % Note that we can obtain the RGB vector with our MATLAB class FSColors.
    %
    datatooltip = struct;

    datatooltip.LineColor = FSColors.yellowish.RGB;
    malfwdplot(out,'datatooltip',datatooltip);
    % now LineColor is not a valid RGB vector, but red (default) will be used
    datatooltip.LineColor = [123 41 12 32 1];
    malfwdplot(out,'datatooltip',datatooltip);
%}
%
%{
    % Interactive_example
    % Another example of the use of option datatooltip.
    % The user can highlight the trajectories of the units that are in
    % the subset at a given step with a mouse click in proximity
    % of that step. A right click will terminate the exercise.
    % To activate this modality, we set the field SubsetLinesColor,
    % which specifies the color used to highlight the trajectories.
    datatooltip = struct;
    datatooltip.SubsetLinesColor = FSColors.purplish.RGB;
    malfwdplot(out,'datatooltip',datatooltip);

    % Here we show that the modality is also activated when
    % SubsetLinesColor is not a valid RGB vector.
    % In this case the default highlight color (blue) is used.
    datatooltip = struct;
    datatooltip.SubsetLinesColor = 999;
    malfwdplot(out,'datatooltip',datatooltip);
%}
%
%{
    % Interactive_example
    %   Example of the use of option databrush.
    %   (brushing is done only once using a rectangular selection tool)
    malfwdplot(out,'databrush',1)
    %   An equivalent statement is
    databrush=struct;
    databrush.selectionmode='Rect';
    malfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    %   Example of the use of brush using a rectangular selection tool and
    %   a cyan colour.
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.FlagColor='c';
    malfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    %   Example of the use of brush using multile selection circular tools.
    databrush=struct;
    databrush.selectionmode='Brush';
    malfwdplot(out,'databrush',databrush);
%}
%
%{
    % Interactive_example
    %   Example of the use of brush using lasso selection tool and fleur
    %   pointer.
    databrush=struct;
    databrush.selectionmode='lasso';
    databrush.Pointer='fleur';
    malfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    %   Example of the use of rectangular brush.
    %    Labels are added for the brushed units. Persistent labels appear in the plot which has
    %   been brushed
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    malfwdplot(out,'databrush',databrush)
%}
%
%   All previuos examples used a non persistent brushing (that is brushing
%   could be done only once). The examples below use persistent brushing
%   (that is brushing can be done multiple times)
%
%{
    % Interactive_example
    %   Example of the use of persistent non cumulative brush. Every time a
    %   brushing action is performed previous highlightments are removed
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.persist='off';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    malfwdplot(out,'databrush',databrush);
%}
%
%{
    % Interactive_example
    %   Example of the use of persistent cumulative brush. Every time a
    %   brushing action is performed current highlightments are added to
    %   previous highlightments
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.persist='on';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    malfwdplot(out,'databrush',databrush)
%}
%
%{
    % Example of the use of some options inside structure fground.
    % load Swiss banknotes
    Y=load('swiss_banknotes.txt');
    [fre]=unibiv(Y);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'init',30);

    % Initialize structure fground
    fground = struct;
    % Specify which trajectories have to be highlighted
    fground.funit=out.Un(end-15:end,2);
    % Specify the steps in which labels have to be put
    n=size(Y,1);
    fground.flabstep=[n/2 n*0.75 n+0.5];;
    % Specify the line width of the highlighted trajectories
    fground.LineWidth=3;
    % Produce a monitoring MD plot in which labels are put for units
    % out.Un(end-15:end,2)in steps [n/2 n*0.75 n+0.5] of the search
    % and store the options to produce the plot inside plotopt
    plotopt=malfwdplot(out,'fground',fground,'msg',2)
    % In order to reuse the options which have been stored inside plotopt
    % use the following sintax
    % malfwdplot(out,plotopt{:})
%}

%{
    % malfwdplot starting from MM estimators.
    n=100;
    v=3;
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out]=MMmulteda(Ycont);
    malfwdplot(out,'conflev',0.99)
%}

%{
    %   Example of use of databrush with RowNames labels shown inside spmplot. 
    close all
    load carsmall
    x1 = Weight;
    x2 = Horsepower;    % Contains NaN data
    y = MPG;    % response
    Y=[x1 x2 y];
    % Remove Nans
    boo=~isnan(y);
    Y=Y(boo,:);
    %   RowLabelsMatrixY is the cell which contains the names of the rows.
    RowLabelsMatrixY=cellstr(Model(boo,:));
    [fre]=unibiv(Y);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'init',30);
    % Row names of Y are added to structure out
    out.label=RowLabelsMatrixY;
    databrush=struct;
    databrush.labeladd='1';
    malfwdplot(out,'databrush',databrush)
%}

%
%   REMARK: note that at present options.databrush and options.datatooltip
%   cannot be used at the same time. These options will be replaced by a
%   single option (options.interact) in future versions of the toolbox.
%

%% Initialization
%

% The number of rows of matrix MDmat is associated with the number of
% units. The number of columns is associated with the number of steps of
% the fwd search.
MDvalues=out.MAL;
[n,nsteps]=size(MDvalues);

% Extract original matrix Y
Y=out.Y;

% seq= column vector containing the sequence 1 to n
seq= (1:n)';

% numtext= a cell of strings used to labels the units with their position
% in the dataset.
numtext=cellstr(num2str(seq,'%d'));

% x= vector which contains the subset size (numbers on the x axis)
x=(n-nsteps+1):n;

% selthdef= threshold to select the MDs labelled in the malfwdplot.
% selline=  threshold to select the MDs highlighted in the malfwdplot.
% unselline= to select how to represent the unselected units ('faint' 'hide' 'greish');
% laby=     label used for the y-axis of the malfwdplot.
if isfield(out,'class')
    if strcmp(out.class,'Smulteda') || strcmp(out.class,'mveeda') || strcmp(out.class,'mcdeda')
        labx= 'Break down point';
    elseif strcmp(out.class,'MMmulteda')
        labx= 'Efficiency';
    else
        labx= 'Subset size m';
    end
else
    labx= 'Subset size m';
end

laby='Mahalanobis distances';
fthresh=2.5^2;
if n>100
    bthresh=2.5^2;
    bstyle='faint';
else
    bthresh=-inf;
    bstyle='';
end

% maximum and minimum MD along the search for each unit
selmax=max(MDvalues,[],2);
selmin=min(MDvalues,[],2);

styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
styp=repmat(styp,ceil(n/13),1);

% Default options for all trajectories
standarddef = struct(...
    'subsize',x,'xlim','','ylim','','titl','','labx',labx,'laby',laby,...
    'Color',{{'b'}},'LineStyle',{{'-'}},...
    'LineWidth',1,'SizeAxesLab',12,'SizeAxesNum',10);

% Default options for the trajectories in foreground
fgrounddef = struct(...
    'fthresh',fthresh,'funit','','flabstep',[],...
    'fmark',0,'LineWidth','','Color','','LineStyle','','FontSize',12);

% Default options for the trajectories in background
bgrounddef=struct('bthresh',bthresh, 'bstyle',bstyle);

conflev='';
options=struct(...
    'standard',standarddef,'fground',fgrounddef,'bground',bgrounddef,...
    'tag','pl_malfwd','datatooltip',1,'label','','databrush','',...
    'nameY','','msg','','conflev',conflev);


%% Preliminary checks

if nargin<1
    error('FSDA:malfwdplot:missingInputs','A required input argument is missing.')
end

%get optional user options
if nargin>1
    UserOptions=varargin(1:2:length(varargin));
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:malfwdplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end


% if LineColor and SubsetLinesColor are not specified, set to default.
LineColor=[1 0 0];      %[1 0 0] is red.
SubsetLinesColor = '';
% Option SubsetLinesColor and LineColor must be removed from cell
% options.datatooltip, because it is not a valid property of the
% datacursormanager class.
datatooltip=options.datatooltip;
if isstruct(datatooltip)
    fdatatooltip=fieldnames(datatooltip);
    % SubsetLinesColor option
    d=find(strcmp('SubsetLinesColor',fdatatooltip));
    if d>0
        SubsetLinesColor=datatooltip.SubsetLinesColor;
        datatooltip=rmfield(datatooltip,'SubsetLinesColor');
        fdatatooltip=fieldnames(datatooltip);
        [a, b]=size(SubsetLinesColor);
        if ~(a==1 && b==3 && sum(SubsetLinesColor<=1)==3)
            % if it is not a valid RGB vector, set to default (blue)
            SubsetLinesColor = [0 0 1];
        end
    end
    % LineColor option
    d=find(strcmp('LineColor',fdatatooltip));
    if d>0
        LineColor=datatooltip.LineColor;
        datatooltip=rmfield(datatooltip,'LineColor');
        %fdatatooltip=fieldnames(datatooltip);
        [a, b]=size(LineColor);
        if ~(a==1 && b==3 && sum(LineColor<=1)==3)
            % if LineColor is not a valid RGB vector, set to default (red)
            LineColor=[1 0 0];
        end
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
        % flagcol=options.databrush{d+1};
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


%% Prepate the figure to display the malfwdplot

% Create a figure to host the plot or clear the existing one
h=findobj('-depth',1,'tag',options.tag);
if (~isempty(h))
    clf(h);
    figure(h);
    axes;
else
    h=figure;
end
set(h,'Name', 'Monitoring of Mahalanobis distances', 'NumberTitle', 'off');
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
% extract the vector associated with the subset size (x)
x=standard.subsize;




% If structure out contains fielname class we check whether input structure
% out comes from MMeda or Seda
if any(strcmp(fieldnames(out),'class'))
    if strcmp(out.class,'MMmulteda')
        x=out.eff;
        out.Un='';
    elseif strcmp(out.class,'Smulteda') || strcmp(out.class,'mveeda') || strcmp(out.class,'mcdeda')
        x=out.bdp;
        out.Un='';
    end
    
end

plot1=plot(x,MDvalues,'tag','data_res','LineWidth',standard.LineWidth);

% Make sure that in the case the length of x is not large the labels on the
% x axis correspond to the values of bdp or of eff
if length(x)<10
    if x(end)<x(1)
        set(gca,'XTick',fliplr(x))
    else
        set(gca,'XTick',x)
    end
end

if any(strcmp(fieldnames(out),'class'))
    if strcmp(out.class,'Smulteda') || strcmp(out.class,'mveeda') || strcmp(out.class,'mcdeda')
        set(gca,'XDir','reverse')
    end
end

% Apply color
scol=standard.Color;

if ~isempty(scol)
    if size(scol,2)>1
        scol=scol';
    end
    scol=repmat(scol,ceil(n/length(scol)),1);
    set(plot1,{'Color'},scol(1:n));
end

% Apply Line Style
slintyp=standard.LineStyle;
if ~isempty(slintyp)
    slintyp=slintyp(:);
    slintyp=repmat(slintyp,ceil(n/length(slintyp)),1);
    
    set(plot1,{'LineStyle'},slintyp(1:n));
end

% save the resfwdplot lines handles, for subsequent use with option persist
plot1lines=plot1;

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
        if max(steps)>n || min(steps)<x(1)
            mess=sprintf(['Warning: steps that you have chosen outside the range [m0 ... n]\n',...
                'are re-assigned to m0 or to n']);
            fprintf('%s\n',mess);
            steps(steps<x(1)) = x(1);
            steps(steps>n) = n;
            steps = sort(unique(steps));
            % before, the steps outside range were not considered
            %steps=steps(steps>=x(1) & steps<=n);
        end
    elseif isempty(fground.flabstep)&& ischar(fground.flabstep)
        steps=[];
        fground.flabstep=steps;
    elseif isempty(fground.flabstep)&& ~ischar(fground.flabstep)
        steps=[x(1) x(end)];
        fground.flabstep=steps;
    else
        steps=[x(1) x(end)];
        fground.flabstep=steps;        
    end
    
    % fthresh= threshold to define units which have to be displayed in foreground
    % (highlighted)
    fthresh=fground.fthresh;
    % funit= List of the units to be displayed in foreground (highlighted)
    funit=fground.funit;
    if ~isempty(funit)
        % Some checks on minimum and maximum of vector funit
        if max(funit)>n || min(funit)<1
            mess=sprintf(['Warning: units that you have chosen outside the range [1 ... n]\n',...
                'are not considered']);
            fprintf('%s\n',mess);
            funit=funit(funit>0 & funit<=n);
        end
    else
        selmax=max(MDvalues,[],2);
        selmin=min(MDvalues,[],2);
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
        set(plot1(funit),{'LineStyle'},slintyp(funit));
    else
        slintyp={'-'};
        slintyp=repmat(slintyp,ceil(n/length(slintyp)),1);
        set(plot1(funit),{'LineStyle'},slintyp(funit));
    end
    
    if ~isempty(fground.flabstep)
        % lsteps = number of steps for which it is necessary to add the labels
        lsteps=length(steps);
        lall=lunits*lsteps;
        
        % indsteps = indexes of the columns of the matrix of the
        % Mahalanobis distances
        % which have to be taken
        % For FS indsteps is simply
        % indsteps=steps-steps(1)+1; however we want to write it in very
        % general terms to cope with S and MM
        indsteps=zeros(lsteps,1);
        for i=1:lsteps
            indsteps(i)=find(x==steps(i));
        end
        
        % HA = the HorizontalAlignment of the labels
        nflabstep = lunits*numel(steps);
        HA = repmat({'center'},nflabstep,1);
        if sum(steps==x(1))
            HA(1:lunits) = {'right'};
        end
        if sum(steps==n)
            HA(nflabstep-lunits+1:nflabstep) = {'left'};
        end
        
        % strings = the labels supplied by the user if they
        % exist, otherwise we simply use the sequence 1 to n
        if isempty(options.label)
            strings = numtext(funit);
        else
            out.label=options.label;
            strings = out.label(funit);
        end
        
        % Label the units
        %         h=text(reshape(repmat(steps,lunits,1),lall,1),...
        %             reshape(MDvalues(funit,steps-x(1)+1),lall,1),...
        %             reshape(repmat(strings,1,lsteps),lall,1),...
        %             'FontSize',fground.FontSize);
        
        % Label the units
        h=text(reshape(repmat(steps,lunits,1),lall,1),...
            reshape(MDvalues(funit,indsteps),lall,1),...
            reshape(repmat(strings,1,lsteps),lall,1),...
            'FontSize',fground.FontSize);
        
        
        
        set(h,{'HorizontalAlignment'},HA)
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
        set(plot1(funit),{'Marker'},styp(funit))
    end
    
end

%% bground options
if ~isempty(options.bground)
    
    if ~isequal(options.bground,bgrounddef)
        
        fld=fieldnames(options.bground);
        
        % Check if user options inside options.fground are valid options
        chkoptions(bgrounddef,fld)
        
        for i=1:length(fld)
            bgrounddef.(fld{i})=options.bground.(fld{i});
        end
    end
    
    % For the options not set by the user use their default value
    bground=bgrounddef;
    
    
    % units = the units which do not have to be modified
    % backunits = the other units which must be plotted using a colormap or
    % which must be hidden or which have to be plotted in greysh
    bthresh=bground.bthresh;
    
    if ~isempty(bthresh) && ischar(bthresh)
        error('FSDA:malfwdplot:WrongBthresh','Specify bthresh as a numeric vector');
    else
        if length(bthresh)>1
            units=seq(selmax>bthresh(2) | selmin<bthresh(1));
        elseif length(bthresh)==1
            units=seq(selmax>bthresh | selmin<-bthresh);
        end
    end
    
    % backunits are defined as the trajectories not belonging to units
    % backunits are associated with unimportant trajectories
    backunits = setdiff(seq,units);
    
    % set line style for trajectories associated with "backunits"
    bstyle=bground.bstyle;
    switch bstyle
        case 'faint'
            % Set to degrading faint blue the color of the 'unimportant
            % trajectories'. Note that stdColor above is 'b', i.e. [0 0 1],
            % and that cyan is [0 1 1].
            Z = rescale(nanmean(abs(MDvalues),2),1,0);
            colormapres = num2cell(colormap([zeros(n,1) Z ones(n,1)]),2);
            set(plot1(backunits),{'Color'},colormapres(backunits,:));
        case 'hide'
            % hide the curves not selected in vector units
            set(plot1(backunits),'Visible','off');
        case 'greysh'
            set(plot1(backunits),'Color',[0.9 0.9 0.9]);
        otherwise
            % do nothing, i.e. leave the default color (blue).
    end
    
end

hold('off')

% include specified tag in the current plot
set(gcf,'tag',options.tag)
set(gca,'Position',[0.1 0.1 0.85 0.85])

% Store the handle of the resfwdplot figure
% Remark: so far plot1 was the vector of the lines in the resfwdplot figure;
%         from now on, plot1 is the handle of the resfwdplot figure.
plot1=gcf;

Un=out.Un;

if ~isempty(options.conflev)
    conflev=options.conflev;
    v=size(out.Loc,2)-1;
    quant = chi2inv(conflev,v);
    rangeaxis=axis;
    V=[rangeaxis(1);rangeaxis(2)];
    QUANT=[quant;quant];
    lwdenv=2;
    numconflev = length(conflev);
    
    % set the string legend for the confidence bands
    legendstring = cell(numconflev,1);
    legendstring(:) = cellstr('% band');
    legendstring2 = strcat(num2str(((conflev)*100)'),legendstring);
    % plot the confidence bands
    hline = line(V, QUANT,'LineWidth',lwdenv,'Tag','conflevline','color','r');
    
    % make the legend for the confidence bands clickable
    clickableMultiLegend(hline(1:numconflev),legendstring2);
    
    % fix the y-axis, otherwise the figure may change if one hides the bands
    % by clicking on the legend
    axis(axis);
end

%% Datatooltip mode (call to function ginputFS)
% This is to highlight trajectories of unit in the subset at given step
if ~isempty(datatooltip) && isstruct(datatooltip) && ~isempty(SubsetLinesColor)
    butt=1;
    % Notice that changing event_obj. Target in subsequent lines
    % seems to affect also hTarget
    hTarget=[];
    
    % save the original XTick values
    xtick_ori=get(gca,'XTick');
    % butt=1 if the left  button is pressed
    % butt=3 if the right button is pressed
    while butt<3
        % Remember to check compatibility of ginput with INUX and MACOSX.
        [x1,~,butt] = ginputFS(1,'right');
        if butt == 99
            return;
        else
            x1=round(x1);
            if x1>n, x1=n; end
            if ~isempty(hTarget)
                % set old line width and old color for old selection
                set(hTarget,{'LineWidth'},hTargetlwd2,{'Color'},hTargetcol2);
            end
            if x1>=out.S2cov(1,1)
                set(gca,'XTick',x1);
                % find the units belonging to subset at selected step
                seqsel=seq(~isnan(out.BB(:,end-n+x1)));
                % find corresponding trajectories (lines)
                aa=get(gca,'Children');
                hTarget=aa(end-seqsel+1);
                
                % Store line width and color of selected trajectory.
                hTargetlwd2=get(hTarget,'LineWidth');
                hTargetcol2=get(hTarget,'Color');
                hTargetlwd22=cell2mat(hTargetlwd2);
                %hTargetcol22=cell2mat(hTargetcol2);
                
                % Increase Line width of trajectories belonging to subset at
                % that particular step. The color is set according to
                % the value of SubsetLinesColor.
                set(hTarget,{'LineWidth'},num2cell(hTargetlwd22+1.5),{'Color'},{SubsetLinesColor});
            end
            if x1 <= min(xtick_ori) || x1 == n
                set(gca,'XTick',xtick_ori);
                if x1 == n
                    set(hTarget,{'LineWidth'},hTargetlwd2,{'Color'},hTargetcol2);
                end
            end
        end
    end
    % set the original XTick values
    set(gca,'XTick',xtick_ori);
    % put an arrow in correspondence of the selected step
    text(x1,min(get(gca,'YLim')),'\uparrow','FontSize',18,'HorizontalAlignment','center','Color',SubsetLinesColor); % or hTargetcol2{end}
end

%% Datatooltip mode (call to function malfwdplotLbl)
if ~isempty(datatooltip)
    hTarget=[];
    hTargetlwd=[];
    hTargetcol=[];
    % datacursormode on;
    hdt = datacursormode;
    set(hdt,'Enable','on');
    % If options.datatooltip is not a struct then use our default options
    if ~isstruct(datatooltip)
        set(hdt,'DisplayStyle','window','SnapToDataVertex','on');
    else
        % options.datatooltip contains a structure where the user can set
        % the properties of the data cursor
        set(hdt,datatooltip);
    end
    
    % Declare a custom datatooltip update function to display additional
    % information about the selected unit
    set(hdt,'UpdateFcn',{@malfwdplotLbl,out,LineColor});
end

%% Datatooltip mode (call to function ginputFS)
% This is to highlight trajectories of unit in the subset at given step

% if ~isempty(options.datatooltip)
%
%     hTarget=[];
%     hTargetlwd=[];
%     hTargetcol=[];
%     % datacursormode on;
%
%
%     % datacursormode on;
%     hdt = datacursormode;
%     set(hdt,'Enable','on');
%     % If options.datatooltip is not a struct then use our default options
%     if ~isstruct(options.datatooltip)
%         set(hdt,'DisplayStyle','window','SnapToDataVertex','on');
%     else
%         % options.datatooltip contains a structure where the user can set the
%         % properties of the data cursor
%         set(hdt,options.datatooltip);
%     end
%     % Declare a custom datatooltip update function to display additional
%     % information about the selected unit
%     set(hdt,'UpdateFcn',{@malfwdplotLbl,out})
% end

%% Brush mode (call to function selectdataFS)
if ~isempty(options.databrush) || isstruct(options.databrush)
    
    if isstruct(options.databrush)
        % transform the input structure databrush into a cell array
        cv=[fdatabrush struct2cell(databrush)]';
        sele=[cv(:)' 'Ignore' {findobj(gcf,'tag','env')}];
        % Add the FlagSize of the brushed points if it has not been
        % previously specified by the user
        d=find(strcmp('FlagSize',fdatabrush));
        if d==0
            sele=[sele 'FlagSize' '3'];
        end
    else
        sele={'selectionmode' 'Rect' 'Ignore' findobj(gcf,'tag','env') };
    end
    
    % sele={sele{:} 'Tag' options.tag}; OLD inefficient code
    sele=[sele 'Tag' {options.tag}];
    
    % Check if X includes the constant term for the intercept.
    v=size(Y,2);
    p1=1:v;
    
    % Set the labels of the axes.
    d=find(strcmp('nameY',fieldnames(out)),1);
    if  isempty(d)
        nameY=cellstr(num2str(p1','y%d'));
    else
        nameY=out.nameY;
    end
    
    % group = vector which will contain the identifier of each group e.g.
    % group(14)=3 means that unit 14 was selected at the third brush
    group=ones(n,1);
    
    % some local variables
    but=0; brushcum=[]; ij=1;
    
    % specify line type for the set of units which are brushed each time
    stypebrushed={'--';':';'-.';'-'};
    stypebrushed=repmat(stypebrushed,10,1);
    
    sele=[sele 'FlagColor' ColorOrd(ij,:) 'FlagMarker' char(styp(ij+1))];
    
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
        
        
        % CALL to selectdataFS
        disp('Select a region to brush in the monitoring MD plot');
        pl = selectdataFS(sele{:});
        
        % exit if the resfwdplot figure was closed before selection
        if isnumeric(pl) && ~isempty(pl) && (pl == -999)
            return
        end
        
        if ~isempty(cell2mat(pl))
            
            % sel = vector which contains the list of brushed units
            sel=seq(cellfun('isempty',pl(1:n)));
            % nbrush = vector which contains the list of the selected units
            nbrush=setdiff(seq,sel);
            disp('Brushed units, yvalue and X values');
            disp([nbrush out.Y(nbrush,:) ]);
            
            % if persist='off', before highlighting the new selection
            if strcmp(persist,'off')
                % set back to default style the previous selection
                set(findobj(plot1,'tag','data_mal'),...
                    'Color',standard.Color{:},'LineStyle',standard.LineStyle{:},'LineWidth',standard.LineWidth);
                % set the standard color style of the 'unimportant units'.
                switch bstyle
                    case 'faint'
                        set(plot1lines(backunits),{'Color'},colormapres(backunits,:));
                    case 'hide'
                        % hide the curves not selected in vector units
                        set(plot1lines(backunits),'Visible','off');
                    case 'greysh'
                        set(plot1lines(backunits),'Color',[0.9 0.9 0.9]);
                    otherwise
                        % do nothing, i.e. leave the default color (blue).
                end
            end
            
            % highlight all trajectories which have been brushed using the
            % same color adopted for brushing figure(plot1);
            linext=findobj(plot1,'Type','line');
            for i=1:length(nbrush)
                set(linext(length(linext)-nbrush(i)+1),'Color',ColorOrd(ij,:));
            end
            % Change in the current plot the style of the lines which have
            % been brushed, using the same line style
            set(linext(length(linext)-nbrush+1),'LineStyle',stypebrushed{ij})
            % Change the line width of the brushed trajectories the width
            % is 1 points bigger than the one of the normal trajectories.
            % linewidthStd=get(linext(length(linext)-nbrush+1),'LineWidth');
            % if iscell(linewidthStd),
            % linewidthStd=max(cell2mat(linewidthStd)); end;
            set(linext(length(linext)-nbrush+1),'LineWidth',standard.LineWidth+1);
            % This line would set a fixed width (3 points, in the example).
            %set(linext(length(linext)-nbrush+1),'LineWidth',3)
            % position(plot1);
        else
            disp('Wrong selection: Try again');
            disp('Select a region to brush in the monitoring MD plot');
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
            
            % nbrush= vector which contains the brushed steps selstesp=
            % vector which will contain the steps in which the brushed
            % units enter the search Given nbrush, find selstesp
            selsteps=zeros(n,11);
            ii=0;
            for i=1:length(nbrush)
                idx = find(Un(:,2:end) == nbrush(i));
                
                % Find the required row(s) of matrix Un
                row = ind2sub(size(Un(:,2:end)),idx);
                % when isempty(row) is true the selected unit entered the
                % subset at step before Un(1,:), that is before the first
                % step which has been monitored
                if ~isempty(row)
                    % Note that mod(row,size(Un,1)) is used rather than
                    % Un(row,:) because the unit could have entered the
                    % subset together with other units (i.e. could be in a
                    % column  of matrix Un different from the second).
                    % Finally note that the expression row-1e-12 is
                    % necessary otherwise when row is= size(Un,1) then
                    % mod(row-1e-12,size(Un,1))is equal to zero
                    selsteps(ii+1:ii+length(row),:)=Un(ceil(mod(row-1e-12,size(Un,1))),:);
                    ii=ii+length(row);
                end
            end
            
            selsteps=sortrows(selsteps(1:ii,:),1);
            % m1 contains the indexes of the unique steps;
            [~, m1]=unique(selsteps(:,1));
            selsteps=selsteps(m1,:);
            % Remark: selsteps=unique(selsteps,'rows') does not seem to
            % work
            
            disp('Steps of entry of brushed units');
            disp(selsteps);
            
            %% - highlight brushed units also in the minimum MD, if it is open
            
            h=findobj('-depth',1,'Tag','pl_mmd');
            
            if (~isempty(h))
                % Remove unnecessary rows from vector selsteps -1 is
                % necessary because we are considering minimum outside
                selsteps=selsteps(:,1)-1;
                
                % make figure which contains mdr become the current figure
                figure(h);
                
                % Condition || but==0 if but=0 then it is necessary to
                % remove previous highlightments (even if persist='on')
                if strcmp(persist,'off') || but==0
                    % If set of values has already been highlighted in the
                    % mmd plot, remove it
                    a=findobj(h,'Tag','brush_mmd');
                    delete(a);
                    
                    % Remove the yellow selection in this plot if present
                    a=findobj(gcf,'Tag','selected');
                    delete(a);
                end
                
                % get the x and y coordinates of mdr
                a=findobj(h,'tag','data_mmd');
                xdata=get(a,'Xdata'); % x coordinates of mmd (steps)
                ydata=get(a,'ydata'); % y coordinates of mmd (values)
                
                [~, ~, ib]=intersect(selsteps, xdata);
                % Stack together x and y coordinates
                xx=[xdata; ydata];
                
                % Just in case the first step of mmd is selected remove it
                % because we also consider ib-1
                ib=ib(ib>1);
                % For each of the brushed units extract coordinates of mmd
                % referred to the step before their entry and the step
                % before
                xxsel=xx(:,[ib-1 ib])';
                % Sort all steps
                xxselr=sortrows(xxsel,1);
                % xxlim=length(nbrush);
                xxlim=length(ib);
                % Reshape previous matrix in such a way that the first
                % length(nbrush) columns refer to the steps which have to
                % be plotted and the remining columns refer to their
                % corresponding values of mdr
                xy=reshape(xxselr,2,2*xxlim);
                % Add to the previous matrix a row of missing values This
                % operation is necessary if the steps are not contiguous
                xy=cat(1,xy,NaN*zeros(1,2*xxlim));
                
                % Reshape the set of x and y coordinates in two column
                % vectors Note the NaN between the steps which are not
                % consecutive
                xcoord=reshape(xy(:,1:xxlim),3*xxlim,1);
                ycoord=reshape(xy(:,xxlim+1:end),3*xxlim,1);
                hold('on');
                if strcmp('on',persist)
                    % If necessary it isalso possible to specify a line
                    % style for the brushed steps
                    % 'LineStyle',stypebrushed{ij},
                    plot(gca,xcoord,ycoord,'LineWidth',1.5,'color',ColorOrd(ij,:),'tag','brush_mdr');
                    set(gca,'Position',[0.1 0.1 0.85 0.85])
                else
                    plot(gca,xcoord,ycoord,'LineWidth',1.5,'color',ColorOrd(1,:),'tag','brush_mdr');
                end
                hold('off');
            end
            
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
            
            % generate the scatterplot matrix
            plo=struct; plo.nameY=nameY; plo.labeladd=labeladd; 
            if max(strcmp('label',fieldnames(out)))>0 && ~isempty(out.label)
                plo.label=out.label(:);
            end
            H = spmplot(Y,group,plo);
            
            % Assign to this figure a name
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
                % - and a function to be executed on figure close is set
                set(gcf,'CloseRequestFcn',@closereqFS);
                
                % Lay down the plots before continuing
                position(plot1);
                disp('Highlight the monitoring MD plot then: click on it to continue brushing or press a keyboard key to stop');
                ss=waitforbuttonpressFS;
                disp('------------------------');
                
                % After waitforbuttonpress:
                % - the standard MATLAB function to be executed on figure
                %   close is recovered
                set(gcf,'CloseRequestFcn','closereq');
                Open_spm = findobj(0, 'type', 'figure','tag','pl_spm');
                Open_mal = findobj(0, 'type', 'figure','tag','pl_malfwd');
                Open_mmd= findobj(0, 'type', 'figure','tag','pl_mmd');
                if isempty(Open_mal)  % User closed the main brushing window
                    if ~isempty(Open_spm), delete(Open_spm); end % spmplot is deleted
                    if ~isempty(Open_mmd), delete(Open_mmd); end % mmdplot is deleted
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
            position(plot1)
        end % for each brushing operation do ...
    end % close loop associated with but (loop brushing)
end % close options.databrush


    function output_txt = malfwdplotLbl(~,event_obj,out,color)
        %% malfwdplotLbl provides information about the selected fwd scaled MD
        %
        % Required input arguments:
        %
        %   obj     =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object
        %               (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the
        %               sense that these arguments are automatically passed
        %               to the function when it executes.
        %       out =   a structure containing the following fields
        %               MAL = a matrix containing the MD monitored
        %               along the search with n rows and n-init+1 columns
        %               Each row of matrix out.MAL is associated with a unit
        %       Y   =   the data matrix
        %       Un  =   a matrix containing the list of the units
        %               which entered the subset in each step of the search
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
        % REMARK: this function is called by function resfwdplot
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
        
        % Increase Line width and set color to 'color' (default is red) or
        % to blue if previous color was 'color') of selected trajectory
        if sum(get(hTarget,'Color')==color)==3
            set(hTarget,'LineWidth',hTargetlwd+1.5,'Color','b');
        else
            set(hTarget,'LineWidth',hTargetlwd+1.5,'Color',color);
        end
        
        pos = get(event_obj,'Position');
        
        % x and y, plot coordinates of the mouse
        x1 = pos(1); y1 = pos(2);
        
        % Find index to retrieve obs. name Consider that find return the
        % linear indexing of matrix xydata
        idx = find(MDvalues == y1,1);
        
        
        % Linear indexing is transformed into normal indexing using
        % function ind2sub row and column contain the column and row
        % indexed of the observation which has been selected with the mouse
        [row,col] = ind2sub(size(MDvalues),idx);
        
        if isempty(row)
            output_txt{1}=['no MD has coordinates x,y' num2str(x1) '' num2str(y1)] ;
        else
            
            % If structure out does not contain labels for the rows then
            % labels row1....rown are added automatically
            if isempty(intersect('label',fieldnames(out)))
                out.label=cellstr(num2str((1:length(out.Y))','row %d'));
            end
            
            
            
            output_txt=cell(5,1);
            
            % output_txt is what it is shown on the screen
            output_txt(1,1) = {['MD equal to: ',num2str(y1,4)]};
            
            % Add information about the corresponding row label of what has
            % been selected
            output_txt{3,1} = ['Unit: ' num2str(cell2mat(out.label(row)))];
            
            if any(strcmp(fieldnames(out),'class'))
                if strcmp(out.class,'MMmulteda')
                    output_txt{2,1} = ['eff=' num2str(x1)];
                    output_txt{4,1} = ['weight=' num2str(out.Weights(row,col))];
                    
                elseif strcmp(out.class,'Smulteda') || strcmp(out.class,'mveeda') || strcmp(out.class,'mcdeda')
                    output_txt{2,1} = ['bdp=' num2str(x1)];
                    output_txt{4,1} = ['weight=' num2str(out.Weights(row,col))];
                else
                    % Add information about the step of the search which is under
                    % investigation
                    output_txt{2,1} = ['Step $m$=' num2str(x1)];
                    
                    % Add information about the step (to be precise the last three
                    % steps) in which the selected unit entered the search
                    idx = find(Un(:,2:end) == row,3,'last');
                    [row,~] = ind2sub(size(Un(:,2:end)),idx);
                    
                    if isempty(row)
                        
                        output_txt{4,1} = ['Unit entered before step $m$=' num2str(Un(1,1))];
                    elseif length(row)<2
                        output_txt{4,1} = ['Unit entered in step $m$=' num2str(Un(row,1))];
                    elseif length(row)==2
                        output_txt{4,1} = ['Unit entered in step $m$=' num2str(Un(row(1),1)) ' and then in step $m$=' num2str(Un(row(2),1))];
                    else
                        output_txt{4,1} = ['Unit entered in steps $m$=' num2str(Un(row(1),1)) ', $m$=' num2str(Un(row(2),1)) ' and $m$=' num2str(Un(row(3),1))];
                    end
                end
            end
            set(0,'ShowHiddenHandles','on');    % Show hidden handles
            hText = findobj('Type','text','Tag','DataTipMarker');
            set(hText,'Interpreter','latex');
            
        end
        
    end

if ~isempty(options.msg)
    if options.msg==1 || options.msg==2
        plotopt=cell(6,1);
        plotopt{1}='standard';
        plotopt{2}=standard;
        plotopt{3}='fground';
        plotopt{4}=fground;
        plotopt{5}='bground';
        plotopt{6}=bground;
    else
    end
    
    if options.msg==2
        disp('standard')
        disp(standard)
        disp('fground')
        disp(fground)
        disp('bground')
        disp(bground)
    end
end

end
%FScategory:VIS-Mult