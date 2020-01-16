function plotopt=resfwdplot(out,varargin)
%resfwdplot plots the trajectories of the monitored scaled (squared) residuals
%
%
%<a href="matlab: docsearchFS('resfwdplot')">Link to the help function</a>
%
%  Required input arguments:
%
%  out :  Structure containing monitoring of scaled residuals. Structure.
%               Structure containing the following fields.
%   out.RES =   matrix containing the residuals monitored in each step of
%               the forward search or any other robust procedure. Every row
%               is associated with a residual (unit). This matrix can be
%               created using function FSReda, Sregeda, MMregeda.
%   out.Un  =   matrix containing the order of entry in the subset of each
%               unit (required only when datatooltip is true or databrush
%               is not empty).
%     out.y  =   a vector containing the response (required only when option
%               databrush is not empty).
%     out.X  =   a matrix containing the explanatory variables (required only
%               when option databrush is not empty).
%   out.Bols =   (n-init+1) x (p+1) matrix containing the estimated beta
%               coefficients monitored in each step of the robust procedure
%               (required only when option databrush is not empty and
%               suboption multivarfit is not empty).
%                Data Types - single|double
%
%
%  Optional input arguments:
%
%           standard : appearance of the plot
%                   in terms of xlim, ylim, axes labels and their font size
%                   style, color of the lines, etc. Structure.
%                   Structure standard contains the following fields
%                   standard.SizeAxesNum  = scalar specifying the fontsize of the
%                       axes numbers. Default value is 10.
%                   standard.xlim = two elements vector with minimum and maximum of
%                       the x axis. Default value is '' (automatic scale).
%                   standard.ylim = two elements vector with minimum and maximum of
%                       the y axis. Default value is '' (automatic scale).
%                   standard.titl = a label for the title (default: '').
%                   standard.labx = a label for the x-axis (default: 'Subset size m').
%                   standard.laby = a label for the y-axis (default:'Scaled residuals'
%                       or 'Scaled squared residuals').
%                  standard. SizeAxesLab = Scalar specifying the fontsize of the
%                       labels of the axes. Default value is 12.
%                   standard.xvalues = vector. x axis values.
%                       Numeric vector of length(size(out.RES,2)) controlling the x
%                       axis coordinates. The default value of xvalues is
%                       size(out.RES,1)-size(out.RES,2)+1:size(out.RES,1)
%                   standard.LineWidth = scalar specifying line width for the
%                       trajectories.
%                   standard.Color = cell array of strings containing the colors to
%                       be used for the highlighted units.
%                       If length(Color)=1 the same color will be used for
%                       all units.
%                       If length(Color)=2 half of the trajectories will
%                       appear with Color{1} and the other half with
%                       Color{2}. And so on with 3 cell elements, etc.
%                   standard.LineStyle = cell containing the line types.
%
%                   The default values of structure standard are:
%                   standard.SizeAxesNum=10
%                   standard.SizeAxesLab=12
%                   standard.xlim='' (Automatic scale)
%                   standard.ylim='' (Automatic scale)
%                   standard.titl='' (empty title)
%                   standard.labx='Subset size m'
%                   standard.laby='Scaled residuals'
%                   standard.LineWidth=1
%                   standard.Color={'b'}
%                   standard.LineStyle={'-'}
%                   standard.xvalues=size(out.RES,1)-size(out.RES,2)+1:size(out.RES,1)
%                   Example - 'standard.LineWidth','1'
%                   Data Types - struct
%
%         fground : trajectories in foregroud.
%                   Structure. Structure which controls which trajectories
%                   are highlighted and how
%                   they are plotted to be distinguishable from the others.
%                   It is possible to control the label, the width, the
%                   color, the line type and the marker of the highlighted
%                   units. The structure fground contains the following
%                   fields:
%                   fground.fthresh = (alternative to funit) numeric vector of
%                       length 1 or 2 which specifies the highlighted
%                       trajectories.
%                       If length(fthresh)=1 the highlighted trajectories
%                       are those of units that throughtout the search had
%                       at leat once a residual greater (in absolute value)
%                       than thresh. The default value of fthresh is 2.5.
%                       If length(fthresh)=2 the highlighted trajectories
%                       are those of units that throughtout the search had
%                       a residual at leat once bigger than fthresh(2) or
%                       smaller than fthresh(1).
%                   fground.funit = (alternative to fthresh) vector containing the
%                       list of the units to be highlighted. If funit is
%                       supplied, fthresh is ignored.
%                   fground.flabstep = numeric vector which specifies the steps of
%                       the search where to put labels for the highlighted
%                       trajectories (units). The default is to put the
%                       labels at the initial and final steps of the search.
%                       flabstep='' means no label.
%                   fground.LineWidth = scalar specifying line width for the
%                       highlighted trajectories (units). Default is 1.
%                   fground.Color = cell array of strings containing the colors to
%                       be used for the highlighted trajectories (units).
%                       If length(scolor)==1 the same color will be used for
%                       all highlighted units Remark: if for example
%                       length(scolor)=2 and there are 6 highlighted units,
%                       3 highlighted trajectories appear with
%                       selunitcolor{1} and 3 highlighted trajectories with
%                       selunitcolor{2}
%                   fground.LineStyle = cell containing the line type of the highlighted
%                       trajectories.
%                   fground.fmark  = scalar controlling whether to plot highlighted
%                       trajectories as markers.
%                       if 1 each line is plotted using a different marker
%                       else no marker is used (default).
%                fground.FontSize = scalar controlling font size of the labels of
%                       the trajectories in foreground.
%
%                   The default values of structure fground are:
%                    fground.fthresh=2.5
%                    fground.flabstep=[m0 n]
%                    fground.LineWidth=1
%                    fground.LineStyle={'-'}
%                    fground.FontSize=12
%
%                   Example - 'fground.LineWidth','1'
%                   Data Types - struct
%                   Remark: if fground='' no unit is highlighted and no
%                   label is inserted into the plot.
%
%         bground : characterictics of the trajectories in background.
%                   Structure.
%                    Structure which specifies the trajectories in background,
%                   i.e. the trajectories corresponding to "unimmportant"
%                   units in the central part of the data. The structure
%                   also specifies the style used in the plot to give them
%                   less emphasis, so that to not distract the eye of the
%                   analyst from the trajectories of the relevant units.
%                   The structure contains the following fields:
%                   bground.bthresh = numeric vector of length 1 or 2 which
%                       specifies how to define the unimmportant trajectories.
%                       Unimmportant trajectories will be plotted using a
%                       colormap, in greysh or will be hidden.
%                       - If length(thresh)=1 the irrelevant units are
%                         those which always had a residual smaller
%                         (in absolute value) than thresh.
%                       - If length(bthresh)=2 the irrelevant units
%                         are those which always had a residual greater
%                         than bthresh(1) and smaller than bthresh(2).
%                       The default is:
%                           bthresh=2.5 if n>100 and bthresh=-inf if n<=100
%                       i.e. to treat all trajectories as important if
%                       n<=100 and, if n>100, to reduce emphasis only to
%                       trajectories having in all steps of the search a
%                       value of scaled residual smaller than 2.5.
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
%                   Remark: bground='' is equivalent to bground.thresh=-Inf
%                   that is all trajectories are considered relevant.
%
%       tag     :    Personalized plot tag. String. String which identifies
%                   the handle of the plot which
%                   is about to be created. The default is to use tag
%                   'pl_resfwd'. Note that if the program finds a plot which
%                   has a tag equal to the one specified by the user, then
%                   the output of the new plot overwrites the existing one
%                   in the same window else a new window is created.
%                   Example - 'tag','myplot'
%                   Data Types - char
%
%   datatooltip :   interactive clicking. Empty value, scalar or
%                   structure.
%                   It is inactive if it is an empty value. The default is
%                   datatooltip = 1, i.e. the user can select with the
%                   mouse an individual residual trajectory in order to
%                   have information about the corresponding unit. The
%                   information displayed depends on the estimator in use.
%                   For example for class FSReda the information concerns
%                   the label and the step of the search in which the unit
%                   enters the subset. If datatooltip is a
%                   structure it may contain the following fields
%                   datatooltip.DisplayStyle = determines how the data cursor displays.
%                       Possible values are 'datatip' and 'window'
%                       (default).
%                       'datatip' displays data cursor information in a
%                       small yellow text box attached to a black square
%                       marker at a data point you interactively select.
%                       'window' displays data cursor information for the
%                       data point you interactively select in a floating
%                       window within the figure.
%                   datatooltip.SnapToDataVertex = specifies whether the
%                       data cursor snaps to the nearest data value or is
%                       located at the actual pointer position. Possible
%                       values are 'on' (default) and 'off'.
%                   datatooltip.LineColor = controls the color of the
%                       trajectory selected with the mouse. It can be an
%                       RGB triplet of values between 0 and 1, or character
%                       vector indicating a color name. Note that a RGB
%                       vector can be conveniently chosen with our MATLAB
%                       class FSColor, see documentation.
%                   datatooltip.SubsetLinesColor = enables to control the
%                       color of the trajectories of the units that are in
%                       the subset at a given step of the search (if
%                       resfwdplot is applied to an object of class
%                       FSReda) or have a weight greater than 0.9 (if
%                       resfwdplot is applied to an object of class LXSeda,
%                       Sregeda and MMregeda).
%                       This can be done (repeatedly) with a left mouse
%                       click in proximity of the step of interest. A right
%                       mouse click will terminate the selection by marking
%                       with a up-arrow the step corresponding to the
%                       highlighted lines. The highlighted lines by default
%                       are in red, but a different color can be specified
%                       as RGB triplet or character of color name.
%                       Note that a RGB vector can be conveniently chosen with
%                       our MATLAB class FSColor, see documentation.
%                       By default SubsetLinesColor ='', i.e. the modality
%                       is not active.
%                       Any initialization for SubsetLinesColor which
%                       cannot be interpreted as RGB vector will be
%                       converted to blue, i.e. SubsetLinesColor will be
%                       forced to be [0 0 1].
%                       If SubsetLinesColor is not empty previous option
%                       LineColor is overlooked.
%                   Example - 'datatooltip',''
%                   Data Types - empty value, scalar or struct
%
%       label   :   row labels. Cell of strings.
%                   Cell containing the labels of the units (optional
%                   argument used when datatooltip=1. If this field is not
%                   present labels row1, ..., rown will be automatically
%                   created and included in the pop up datatooltip window).
%                   Example - 'label',{'Smith','Johnson','Robert','Stallone'}
%                   Data Types - cell
%
%    databrush  :   interactive mouse brushing. Empty value, scalar or structure.
%                   If databrush is an empty value (default), no brushing
%                   is done.
%                   The activation of this option (databrush is a scalar or
%                   a cell) enables the user  to select a set of
%                   trajectories in the current plot and to see them
%                   highlighted in the y|X plot, i.e. a matrix of scatter
%                   plots of y against each column of X, grouped according
%                   to the selection(s) done by brushing. If the plot y|X
%                   does not exist it is automatically created.
%                   In addition, brushed units are automatically highlighted in the
%                   minimum deletion residual plot if it is already open.
%                   The extension to the following plots will be available
%                   in future versions of the toolbox:
%                   - monitoring leverage plot;
%                   - maximum studentized residual;
%                   - s^2 and R^2;
%                   - Cook distance and modified Cook distance;
%                   - deletion t statistics.
%                   Note that the window style of the other figures is set
%                   equal to that which contains the monitoring residual
%                   plot. In other words, if the monitoring residual plot
%                   is docked all the other figures will be docked too
%                   DATABRUSH IS A SCALAR
%                   If databrush is a scalar the default selection tool is
%                   a rectangular brush and it is possible to brush only
%                   once (that is persist='').
%                   DATABRUSH IS A STRUCTURE.
%                   If databrush is a structure, it is possible to use all
%                   optional arguments of function selectdataFS.m and the
%                   following optional argument:
%                   - persist. Persist is an empty value or a scalar
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
%                   - bivarfit. This option adds one or more least
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
%                   - multivarfit: This option adds one or more least square
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
%                       at their centroid (that is $\overline{y}'C$)).
%                     multivarfit = '2'
%                       exactly equal to multivarfit ='1' but this time we
%                       add the line based on the group of unselected
%                       observations.
%                   - labeladd. If this option is '1', we label the units
%                     of the last selected group with the unit row index in
%                     matrices X and y. The default value is labeladd='',
%                     i.e. no label is added.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%
%       nameX   :  labels of the explanatory variables. Cell. Cell array of
%                   strings of length p containing the labels
%                   of the variables of the regression dataset. If it is
%                   empty the sequence X1, ..., Xp will be created
%                   automatically.
%                   Example - 'nameX',{'var1', var2, 'var3'}
%                   Data Types - cell of strings
%
%       namey   :   response label. Character. Character containing the
%                   label of the response.
%                   Example - 'namey','response'
%                   Data Types - character
%
%       msg     :   display  used options. Scalar.
%                   Scalar which controls whether to display
%                   as output the options inside structures standard,
%                   fground and bground which have been used to draw the
%                   plot.
%                   plotopt=resfwdplot(out,'msg',1) prints on the screen
%                   the options which have been used to draw
%                   the three types of trajectories (standard, foreground
%                   and background).
%                   Example - 'msg',1
%                   Data Types - single or double
%
% Output:
%
%   plotopt : options which have been used to create the plot. Cell array
%               of strings. Store all options which have been used to
%               generate the plot inside cell plotopt.
%
% See also: malfwdplot
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('resfwdplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    % Monitoring residuals plot with all the default options.
    % generate input structure for the resfwdplot
    % In this example FS estimator is used
    n=100;
    y=randn(n,1);
    X=randn(n,4);
    [out]=LXS(y,X,'nsamp',1000);
    [out]=FSReda(y,X,out.bs);
    % Produce  resfwdplot
    resfwdplot(out);
%}
%

%
%{
    % Example of the use of some options inside structure standard.
    % Initialize structure standard
    % Specify the steps in which labels have to be put
    standard=struct;
    standard.LineStyle={'-';'-.';':'};
    % Specify the line width
    standard.LineWidth=0.5;
    resfwdplot(out,'standard',standard)
%}
%
%{
    % Example of the use of some options inside structure fground.
    % Initialize structure fground
    fground = struct;
    % Specify which trajectories have to be highlighted
    fground.funit=[2 5 20 23 35 45];
    % Specify the steps in which labels have to be put
    n=length(y);
    fground.flabstep=[n/2 n*0.75 n+0.5];;
    % Specify the line width of the highlighted trajectories
    fground.LineWidth=3;
    % Produce a monitoring residuals plot in which labels are put for units
    % [2 5 20 23 35 45] in steps [n/2 n*0.75 n+0.5] of the search
    resfwdplot(out,'fground',fground)
%}
%
%
%{
    % Example of the use of some options inside structure bground.
    bground = struct;
    % Specify a threshold to define the "background" trajectories
    bground.bthresh=2;
    % Trajectories whose residual is always between -btresh and +bthresh
    % are shown as specified in bground.bstyle
    bground.bstyle='hide';
    resfwdplot(out,'bground',bground)
%}
%
%
%{
    % Example of the use of option datatooltip.
    % Gives the user the possibility of clicking on the different points
    % and have information about the unit selected, the step of entry
    % into the subset and the associated label.
    datatooltip = struct;
    % In this example the style of the datatooltip is 'datatip'. Click on a
    % trajectory when the resfwdplot is displayed.
    %
    datatooltip.DisplayStyle = 'datatip';
    resfwdplot(out,'datatooltip',datatooltip);
    %
    % Now we use the default style, which is 'window'.
    datatooltip.DisplayStyle = 'window';
    resfwdplot(out,'datatooltip',datatooltip);

    % Here we specify the RGB color used to highlight the selected trajectory.
    % Note that we can obtain the RGB vector with our MATLAB class FSColors.
    %
    datatooltip = struct;

    datatooltip.LineColor = FSColors.yellowish.RGB;
    resfwdplot(out,'datatooltip',datatooltip);
    % now LineColor is not a valid RGB vector, but red (default) will be used
    datatooltip.LineColor = [123 41 12 32 1];
    resfwdplot(out,'datatooltip',datatooltip);
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
    resfwdplot(out,'datatooltip',datatooltip);

    % Here we show that the modality is also activated when
    % SubsetLinesColor is not a valid RGB vector.
    % In this case the default highlight color (blue) is used.
    datatooltip = struct;
    datatooltip.SubsetLinesColor = 999;
    resfwdplot(out,'datatooltip',datatooltip);
%}
%
%{
    % Interactive_example
    %   Example of the use of option databrush.
    %   (brushing is done only once using a rectangular selection tool)
    resfwdplot(out,'databrush',1)
    %   An equivalent statement is
    databrush=struct;
    databrush.selectionmode='Rect';
    resfwdplot(out,'databrush',databrush);
%}
%
%{
    % Interactive_example
    %   Example of the use of brush using a rectangular selection tool.
    % Use a cyan colour.
    resfwdplot(out,'databrush',{'selectionmode' 'Rect' 'FlagColor' 'c'})
%}
%
%{
    % Interactive_example
    % Example of the use of brush using multiple selection circular tools.
    databrush=struct;
    databrush.selectionmode='Brush';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    %   Example of the use of brush using lasso selection tool and fleur pointer.
    databrush=struct;
    databrush.selectionmode='lasso';
    databrush.Pointer='fleur';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Example of the use of rectangular brush.
    %  We havesuperimposed labels
    %   for the brushed units and persistent labels in the plot which has
    %   been brushed.
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    resfwdplot(out,'databrush',databrush)
%}
%
%   All previous examples used a non persistent brushing (that is brushing
%   could be done only once). The examples below use persistent brushing
%   (that is brushing can be done multiple times)
%
%{
    % Interactive_example
    % Example of persistent brushing.
    %   Example of the use of persistent non cumulative brush. Every time a
    %   brushing action is performed previous highlightments are removed
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    databrush.persist='off';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Example of persistent cumulative brush.
    %   Example of the use of persistent cumulative brush. Every time a
    %   brushing action is performed current highlightments are added to
    %   previous highlightments
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='on';
    databrush.persist='on';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Example of persistent cumulative brush with ols line.
    %   The same as before, but also fit one ols line to each selected group
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='on';
    databrush.persist='on';
    databrush.bivarfit='0';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Example of persistent cumulative brush with a single ols line.
    %   The same but now fit a single ols line to all data.
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    databrush.persist='on';
    databrush.bivarfit='1';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Example of persistent cumulative brush with two ols lines.
    %   The same but now fit a first ols line to all data and a second line
    %   on the group of observations which remain unselected.
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    databrush.persist='on';
    databrush.bivarfit='2';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Example of persistent cumulative brush with ols line for a particular group.
    %   The same but now fit a single ols line to the group with index 4.
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    databrush.persist='on';
    databrush.bivarfit='i4';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Example of persistent cumulative brush with personalized line.
    %   The same but now add the line mean(y)+Ci*Xi.
    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    databrush.persist='on';
    databrush.multivarfit='1';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Example of persistent cumulative brush with option multivarfit.
    % Example of use of option databrush suboptions multivarfit 2, persist on
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);
    [out]=LXS(y,X,'nsamp',10000);
    [out]=FSReda(y,X,out.bs);
    out1=out;
    out1.RES=out.RES.^2;
    mdrplot(out,'ylimy',[1 5]);

    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='off';
    databrush.persist='on';
    databrush.multivarfit='2';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Interactive_example
    % Brushing in the loyalty cards data example.
    % Fidelity cards data
    XX=load('loyalty.txt');
    namey='Sales'
    y=XX(:,end);
    y=y.^0.4;
    nameX={'Number of visits', 'Age', 'Number of persons in the family'};

    X=XX(:,1:end-1);
    [out]=LXS(y,X,'nsamp',10000);
    [out]=FSReda(y,X,out.bs);
    mdrplot(out,'ylimy',[1.6 4.5],'xlimx',[420 510],'lwd',2,'FontSize',16,'SizeAxesNum',16,'lwdenv',2);
    mdrplot(out,'ylimy',[1.6 4.5],'envm',489,'xlimx',[420 510],'lwd',2,'FontSize',16,'SizeAxesNum',16,'lwdenv',2);

    databrush=struct;
    databrush.selectionmode='Rect';
    databrush.Label='on';
    databrush.RemoveLabels='on';
    databrush.persist='on';
    databrush.multivarfit='2';
    resfwdplot(out,'databrush',databrush)
%}
%
%{
    % Example of the  use of nargout.
    % Fidelity cards data
    XX=load('loyalty.txt');
    namey='Sales'
    y=XX(:,end);
    y=y.^0.4;
    X=XX(:,1:end-1);
    [out]=LXS(y,X,'nsamp',10000);
    [out]=FSReda(y,X,out.bs);
    % msg is 1 therefore used options are also shown on the screen.
    plotopt=resfwdplot(out,'msg',1)
    % In order to reuse the options which have been stored inside plotopt
    % use the following sintax
    % resfwdplot(out,plotopt{:})
%}

%{
    % Example of use of option xvalues inside structure standard.
    XX=load('loyalty.txt');
    namey='Sales'
    y=XX(:,end);
    y=y.^0.4;
    X=XX(:,1:end-1);
    [out]=LXS(y,X,'nsamp',10000);
    [out]=FSReda(y,X,out.bs);
    standard=struct;
    % xlabels start from 400
    standard.xvalues=400:1:(size(out.RES,2)+400-1);
    resfwdplot(out,'standard',standard)
%}

%{
    % Monitoring residuals plot using S estimator.
    n=100;
    y=randn(n,1);
    X=randn(n,4);
    [out]=Sregeda(y,X);
    % Produce  resfwdplot
    resfwdplot(out);
%}

%{
    % Interactive_example
    % First example of brushing with row names.
    % The labels of the units appear in the resfwdplot
    % and automatically in the yXplot with the same color.
    % In this example brushing is done just once. To
    % check what happens when persist is 'on', see the next example.
    close all
    load carsmall
    x1 = Weight;
    x2 = Horsepower;    % Contains NaN data
    y = MPG;    % response
    X=[x1 x2];
    % Remove Nans
    boo=~isnan(y);
    y=y(boo,:);
    X=X(boo,:);
    RowLabelsy=cellstr(Model(boo,:));
    [out]=LXS(y,X,'nsamp',10000);
    [out]=FSReda(y,X,out.bs);

    % Write labels of trajectories inside the resfwdplot while brushing
    databrush.Label='on';
    % Do not remove labels after selection in the resfwdplot
    databrush.RemoveLabels='off';
    % Add the labels of the units in the associated yXplot matrix
    databrush.labeladd='1'; %
    resfwdplot(out,'databrush',databrush,'label',RowLabelsy)
%}

%{
    % Interactive_example
    % Second example of brushing with row names.
    % This example is exactly equal as before but now persist is 'on'.
    % In this case each set of brushed units appears with a particular
    % color (both the points and the associated labels)
    close all
    load carsmall
    x1 = Weight;
    x2 = Horsepower;    % Contains NaN data
    y = MPG;    % response
    X=[x1 x2];
    % Remove Nans
    boo=~isnan(y);
    y=y(boo,:);
    X=X(boo,:);
    RowLabelsy=cellstr(Model(boo,:));
    [out]=LXS(y,X,'nsamp',10000);
    [out]=FSReda(y,X,out.bs);

    % Write labels of trajectories inside the resfwdplot while brushing
    databrush.Label='on';
    % Do not remove labels after selection in the resfwdplot
    databrush.RemoveLabels='off';
    % Add the labels of the units in the associated yXplot matrix
    databrush.labeladd='1';
    % Enable repeated brushing actions
    databrush.persist='on';
    resfwdplot(out,'databrush',databrush,'label',RowLabelsy)
%}

%
%   REMARK: note that at present options.databrush and options.datatooltip
%   cannot be used at the same time. These options will be replaced by a
%   single option (options.interact) in future versions of the toolbox.
%

%% Beginning of code 

% Initialization

% Close existing yXplot and resfwdplot, unless the fanplot is open.
% if isempty(findobj('type','figure','Tag','pl_fan'))
%     close(findobj('type','figure','Tag','pl_yX'));
%     close(findobj('type','figure','Tag','pl_resfwd'));
% end


% The rows of matrix residuals are associated with the units in the dataset.
% The columns are associated with the steps of the fwd search.
residuals = out.RES;
[n,nsteps]= size(residuals);

% seq: column vector containing the sequence 1 to n
seq = (1:n)';

% numtext= a cell of strings used to labels the units (that is the rownames of the units)
numtext=cellstr(num2str(seq,'%d'));

% x: vector which contains the subset size (numbers on the x axis)
x = (n-nsteps+1):n;

% maximum and minimum residual along the search for each unit
selmax=max(residuals,[],2);
selmin=min(residuals,[],2);

% default values for fthresh, bthresh, bstyle, labx and laby.
if strcmp(out.class,'Sregeda')
    labx= 'Break down point';
elseif strcmp(out.class,'MMregeda')
    labx= 'Efficiency';
else
    labx= 'Subset size m';
end


if min(min(residuals))<0
    fthresh=2.5;
    laby='Scaled residuals';
    if n>100
        bthresh=2.5;
        bstyle='faint';
    end
else
    laby='Squared scaled residuals';
    fthresh=2.5^2;
    if n>100
        bthresh=2.5^2;
        bstyle='faint';
    end
end
if n<=100
    bthresh=-inf;
    bstyle='';
end

styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
styp=repmat(styp,ceil(n/13),1);

% Default options for all trajectories
standarddef = struct(...
    'xvalues',x,'xlim','','ylim','','titl','','labx',labx,'laby',laby,...
    'Color',{{'b'}},'LineStyle',{{'-'}},...
    'LineWidth',1,'SizeAxesLab',12,'SizeAxesNum',10);

% Default options for the trajectories in foreground
fgrounddef = struct(...
    'fthresh',fthresh,'funit','','flabstep','',...
    'fmark',0,'LineWidth','','Color','','LineStyle','','FontSize',12);

% Default options for the trajectories in background
bgrounddef=struct('bthresh',bthresh, 'bstyle',bstyle);

options=struct(...
    'standard',standarddef,'fground',fgrounddef,'bground',bgrounddef,...
    'tag','pl_resfwd','datatooltip',1,'label','','databrush','',...
    'nameX','','namey','','msg','');

%% Preliminary checks

if nargin<1
    error('FSDA:resfwdplot:missingInputs','A required input argument is missing.')
end

%get optional user options
if nargin>1
    UserOptions=varargin(1:2:length(varargin));
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:resfwdplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
%% Prepate the figure to display the resfwdplot

% Create a figure to host the plot or clear the existing one
h=findobj('-depth',1,'tag',options.tag);
if (~isempty(h))
    clf(h);
    figure(h);
    axes;
else
    h=figure;
end
set(h,'Name', 'Monitoring of Scaled Residuals', 'NumberTitle', 'off');
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
x=standard.xvalues;

% If structure out contains fielname class we check whether input structure
% out comes from MMeda or Seda
if any(strcmp(fieldnames(out),'class'))
    if strcmp(out.class,'MMregeda')
        x=out.eff;
        out.Un='';
    elseif strcmp(out.class,'Sregeda')
        x=out.bdp;
        out.Un='';
    elseif strcmp(out.class,'MPDPeda')
        x=out.alpha;
        out.Un='';
    end
    
end

plot1=plot(x,residuals,'tag','data_res','LineWidth',standard.LineWidth);

if strcmp(out.class,'Sregeda') || strcmp(out.class,'MPDPeda')
    set(gca,'XDir','reverse')
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
    if size(slintyp,2)>1
        slintyp=slintyp';
    end
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
        if max(steps)>x(end) || min(steps)<x(1)
            mess=sprintf(['Warning: steps that you have chosen outside the range [m0 ... n]\n',...
                'are re-assigned to m0 or to n']);
            fprintf('%s\n',mess);
            steps(steps<x(1)) = x(1);
            steps(steps>x(end)) = x(end);
            steps = sort(unique(steps));
        end
    else
        steps=[x(1) x(end)];
        fground.flabstep=steps;
    end
    
    % fthresh= threshold to define units which have to be displayed in
    % foreground (highlighted)
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
    
    % Specify the line type for the highlighted units (those forming vector
    % funit)
    slintyp=fground.LineStyle;
    slintyp=slintyp(:);
    
    if ~isempty(slintyp)
        slintyp=repmat(slintyp,ceil(n/length(slintyp)),1);
        set(plot1(funit),{'LineStyle'},slintyp(1:length(funit)));
        
    else
        slintyp={'-';'--';':';'-.'};
        % slintyp={'-'};
        slintyp=repmat(slintyp,ceil(n/length(slintyp)),1);
        
        set(plot1(funit),{'LineStyle'},slintyp(1:length(funit)));
        
    end
    
    if ~isempty(fground.flabstep)
        % lsteps = number of steps for which we add the labels
        lsteps=length(steps);
        lall=lunits*lsteps;
        % indsteps = indexes of the columns of the matrix named residuals
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
        %         if isempty(options.label)
        %
        %             % numtext: cell of strings used to label the units and their position in
        %             % the dataset
        %             numtext = cellstr(num2str(seq,'%d'));
        %             strings = numtext(funit);
        %         else
        %             %             out.label=options.label;
        %             %             strings = out.label(funit);
        %             %             numtext=out.label;
        %             strings = options.label(funit);
        %             numtext=options.label;
        %         end
        
        %%%%%%%%%%%%%%%%%%%%
        if isempty(options.label)
            % In old releases of FSDA it was possible to supply row names
            % directly from input structure out, so for compatibility we
            % leave the instruction below
            % If structure out does not contain labels for the rows then
            % labels row1....rown are added automatically
            if isempty(intersect('label',fieldnames(out)))
                out.label=numtext;
            end
            strings = numtext(funit);
        else
            numtext=options.label;
            out.label=options.label;
            strings = numtext(funit);
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%
        
        
        % Label the units
        h=text(reshape(repmat(steps,lunits,1),lall,1),...
            reshape(residuals(funit,indsteps),lall,1),...
            reshape(repmat(strings,1,lsteps),lall,1),...
            'FontSize',fground.FontSize);
        set(h,{'HorizontalAlignment'},HA)
    end
    
    % if requested, set the color of the selected trajectories note that if
    % scolor contains more than one color, e.g. options.scolor =
    % {'b';'g';'r'}, then the colors of the trajectories alternate.
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
    
    % units = the units which do not have to be modified backunits = the
    % other units which must be plotted using a colormap or which must be
    % hidden or which have to be plotted in greysh
    bthresh=bground.bthresh;
    
    if ~isempty(bthresh) && ischar(bthresh)
        error('FSDA:resfwdplot:WrongBthresh','Specify bthresh as a numeric vector');
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
            Z = rescaleFS(nanmean(abs(residuals),2),1,0);
            colormapres = num2cell(colormap([zeros(n,1) Z ones(n,1)]),2);
            set(plot1(backunits),{'Color'},colormapres(backunits,:));
        case 'hide'
            % hide the curves not selected in vector units
            set(plot1(backunits),'Visible','off');
        case 'greysh'
            % set(plot1(backunits),'Color',[0.9 0.9 0.9]);
            set(plot1(backunits),'Color',[0.7 0.9 1]);
        otherwise
            % do nothing, i.e. leave the default color (blue).
    end
    
end

hold('off')

% include specified tag in the current plot
set(gcf,'tag',options.tag)
set(gca,'Position',[0.1 0.1 0.85 0.85])

% Store the handle of the resfwdplot figure.
% Remark: so far plot1 was the vector of the lines in the resfwdplot figure;
%         from now on, plot1 is the handle of the resfwdplot figure.
plot1=gcf;

Un=out.Un;
X=out.X;
y=out.y;

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
            if x1>=out.S2(1,1)
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

%% Datatooltip mode (call to function resfwdplotLbl)
if ~isempty(datatooltip)
    hTarget=[];
    hTargetlwd=[];
    hTargetcol=[];
    try
        chkgpu=gpuDevice; %#ok<NASGU>
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
        set(hdt,'UpdateFcn',{@resfwdplotLbl,out,LineColor});
    catch
        disp('No graphical device, interactive datatooltip not enabled')
    end
end

%% Brush mode (call to function selectdataFS)
if ~isempty(options.databrush) || isstruct(options.databrush)
    
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
    
    sele=[sele 'Tag' {options.tag}  'RowNamesLabels' {numtext}];
    
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
        
        % call to selectdataFS
        disp('Select a region to brush in the monitoring residual plot');
        pl = selectdataFS(sele{:});
        
        % exit if the resfwdplot figure was closed before selection
        if isnumeric(pl) && ~isempty(pl) && (pl == -999)
            return
        end
        
        if ~isempty(cell2mat(pl))
            
            % sel = vector which contains the list of unbrushed units
            sel=seq(cellfun('isempty',pl(1:n)));
            % nbrush = vector which contains the list of brushed units
            nbrush=setdiff(seq,sel);
            disp('Brushed units, yvalue and X values');
            disp([nbrush out.y(nbrush) out.X(nbrush,:)]);
            
            % if persist='off', before highlighting the new selection
            if strcmp(persist,'off')
                % set back to default style the previous selection
                set(findobj(plot1,'tag','data_res'),...
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
            disp('Select a region to brush in the monitoring residual plot');
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
            
            %% - highlight brushed units also in the minimum deletion residual, if it is open
            
            h=findobj('-depth',1,'Tag','pl_mdr');
            
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
                    % mdr plot, remove it
                    a=findobj(h,'Tag','brush_mdr');
                    delete(a);
                    
                    % Remove the yellow selection in this plot if present
                    a=findobj(gcf,'Tag','selected');
                    delete(a);
                end
                
                % get the x and y coordinates of mdr
                a=findobj(h,'tag','data_mdr');
                xdata=get(a,'Xdata'); % x coordinates of mdr (steps)
                ydata=get(a,'ydata'); % y coordinates of mdr (values)
                
                [c, ia, ib]=intersect(selsteps, xdata); %#ok<ASGLU>
                % Stack together x and y coordinates
                xx=[xdata; ydata];
                
                % Just in case the first step of mdr is selected remove it
                % because we also consider ib-1
                ib=ib(ib>1);
                % For each of the brushed units extract coordinates of mdr
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
            
            
            %             % generate the scatterplot matrix
            %             plo=struct;
            %             % plo.namey=namey; plo.nameX=nameX;
            %             plo.clr=clr(unigroup);
            %             plo.sym=char(styp{unigroup});
            %             plo.labeladd=labeladd;
            %             if max(strcmp('label',fieldnames(out)))>0 && ~isempty(out.label)
            %                 plo.label=out.label(:);
            %             end
            %             [H,AX,BigAx]=yXplot(y,Xsel,group,plo,'namey',namey,'nameX',nameX);
            
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
            
            % add objects to the panels of the yX
            % add2yX(H,AX,BigAx,out,group,nbrush,bivarfit,multivarfit,labeladd);
            add2yX(H,AX,BigAx,'intercept',intercept,'bivarfit',bivarfit,...
                'multivarfit',multivarfit,'labeladd',labeladd,'RowNamesLabels',numtext);
            
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
                disp('Highlight the monitoring residuals plot then: click on it to continue brushing or press a keyboard key to stop');
                ss=waitforbuttonpressFS;
                disp('------------------------');
                
                % After waitforbuttonpress:
                % - the standard MATLAB function to be executed on figure
                %   close is recovered
                set(gcf,'CloseRequestFcn','closereq');
                Open_yX = findobj(0, 'type', 'figure','tag','pl_yX');
                Open_res = findobj(0, 'type', 'figure','tag','pl_resfwd');
                Open_mdr = findobj(0, 'type', 'figure','tag','pl_mdr');
                if isempty(Open_res)  % User closed the main brushing window
                    if ~isempty(Open_yX); delete(Open_yX); end    % yX plot is deleted
                    if ~isempty(Open_mdr); delete(Open_mdr); end  % mdr plot is deleted
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


    function output_txt = resfwdplotLbl(~,event_obj,out,color)
        %% resfwdplotLbl provides information about the selected fwd scaled residual
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
        %               RES = a matrix containing the residuals monitored
        %               along the search with n rows and n-init+1 columns
        %               Each row of matrix
        %               out.RES is associated with a unit
        %               y   =   the response of the regressione model
        %               Un  =   a matrix containing the list of the units
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
        
        % Find index to retrieve obs. name. Consider that find returns the
        % linear indexing of matrix xydata
        % residuals=out.RES;
        idx = find(residuals == y1,1);
        
        % Linear indexing is transformed into normal indexing using
        % function ind2sub row and column contain the column and row
        % indexed of the observation which has been selected with the mouse
        [row,col] = ind2sub(size(residuals),idx);
        
        if isempty(row)
            output_txt{1}=['no residual has coordinates x,y' num2str(x1) '' num2str(y1)] ;
        else
            
            output_txt=cell(5,1);
            
            % output_txt is what is shown on the screen
            output_txt(1,1) = {['Residual equal to: ',num2str(y1,4)]};
            
            % Add information about the corresponding row label of what has
            % been selected
            output_txt{3,1} = ['Unit name: ' num2str(cell2mat(out.label(row)))];
            
            if any(strcmp(fieldnames(out),'class'))
                if strcmp(out.class,'MMregeda')
                    output_txt{2,1} = ['eff=' num2str(x1)];
                    output_txt{4,1} = ['weight=' num2str(out.Weights(row,col))];
                    
                elseif strcmp(out.class,'Sregeda')
                    output_txt{2,1} = ['bdp=' num2str(x1)];
                    output_txt{4,1} = ['weight=' num2str(out.Weights(row,col))];
                    
                elseif strcmp(out.class,'MPDPeda')
                    output_txt{2,1} = ['alpha=' num2str(x1)];
                    
                else
                    % Add information about the step of the search which is under
                    % investigation or the value of bdp or the value of eff
                    output_txt{2,1} = ['Step m=' num2str(x1)];
                    
                    
                    % Add information about the step (to be precise the last three
                    % steps) in which the selected unit entered the search
                    idx = find(Un(:,2:end) == row,3,'last');
                    [row,~] = ind2sub(size(Un(:,2:end)),idx);
                    if isempty(row)
                        output_txt{4,1} = ['Unit entered before step m=' num2str(Un(1,1))];
                    elseif length(row)<2
                        output_txt{4,1} = ['Unit entered in step m=' num2str(Un(row,1))];
                    elseif length(row)==2
                        output_txt{4,1} = ['Unit entered in step m=' num2str(Un(row(1),1)) ' and then in step m=' num2str(Un(row(2),1))];
                    else
                        output_txt{4,1} = ['Unit entered in steps m=' num2str(Un(row(1),1)) ', m=' num2str(Un(row(2),1)) ' and m=' num2str(Un(row(3),1))];
                    end
                end
            end
        end
    end

if options.msg==1
    disp('standard')
    disp(standard)
    disp('fground')
    disp(fground)
    disp('bground')
    disp(bground)
end

% Return the handles from boxplot only if the user requested them.
if nargout>0
    plotopt=cell(6,1);
    plotopt{1}='standard';
    plotopt{2}=standard;
    plotopt{3}='fground';
    plotopt{4}=fground;
    plotopt{5}='bground';
    plotopt{6}=bground;
end

end
%FScategory:VIS-Reg