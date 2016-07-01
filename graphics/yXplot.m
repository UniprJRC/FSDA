function [plot1]=yXplot(y,X,varargin)
%yXplot produces an interactive scatterplot of y against each variable of X in the input dataset.
%
%<a href="matlab: docsearchFS('yXplot')">Link to the help function</a>
%
%  Required input arguments:
%
%    y: Response variable. Vector. A vector with n elements that contains
%       the response variable. y can be either a row or a column vector.
%    X: Data matrix of explanatory variables (also called 'regressors') of
%       dimension. Rows of X represent observations, and columns
%       represent variables.
%
%      or a structure out containing the following fields
%
%               yXplot(out,varargin)
%
%       y   =   a vector containing the response
%       X   =   a matrix containing the explanatory variables
%
%               If out contains just the two above fields the yXplot will
%               be immediately created. On the other hand, if out also
%               contains information about the search, it is possible to
%               exploit the brushing (i.e. the automatic
%               interaction with the other plots) and datatooltip
%               possibilities
%
%       RES =   matrix containing the residuals monitored in each
%               step of the forward search. Every row is associated with a
%               residual (unit).
%               This matrix can be created using function FSReda
%               (compulsory argument)
%       Un  =   matrix containing the order of entry of each unit
%               (necessary if datatooltip is true or databrush is not empty)
%
%
%  Optional input arguments:
%
%          group:   identifier vector. Vector. 
%                   Vector with n elements, grouping variable that
%                   determines the marker and color assigned to each point.
%                   It can be a categorical variable, vector, string
%                   matrix, or cell array of strings.
%               Example - group, 
%               Data Types - single or double
%       nameX   :   explanatory variables names. Cell. Cell array of strings of length p containing the labels
%                   of the varibles of the regression dataset. If it is empty
%                 	(default) the sequence X1, ..., Xp will be created
%                   automatically
%               Example - group, 
%               Data Types - single or double
%       namey   :   response variable name. Character. Character containing the label of the response
%                   Example - 'namey','response'
%                   Data Types - character
%       ylim    :   y limits. Vector. vector with two elements controlling
%                   minimum and maximum on the y axis. Default value is ''
%                   (automatic scale)
%                   Example - 'ylim',[0 4]
%                   Data Types - single or double
%       xlim    :   x limits. Vector. vector with two elements controlling minimum and maximum
%                   on the x axis. Default value is '' (automatic scale)
%                   Example - 'xlim',[0 4]
%                   Data Types - single or double
%       tag     :   plot tag. String. string which identifies the handle of the plot which
%                   is about to be created. The default is to use tag
%                   'pl_yX'. Notice that if the program finds a plot which
%                   has a tag equal to the one specified by the user, then
%                   the output of the new plot overwrites the existing one
%                   in the same window else a new window is created
%                   Example - 'tag','myplot'
%                   Data Types - character
%       selunit :   unit labelling. Cell array of strings, string, or numeric vector for
%                   labelling units. If out is a structure the threshold is
%                   associated with the trajectories of the residuals
%                   monitored along the search else it refers to the values
%                   of the response variable.
%                   If it is a cell array of strings, only
%                   the lines associated with the units that in at least
%                   one step of the search had a residual smaller than
%                   selunit{1} or greater than selline{2} will have a
%                   textbox.
%                   If it is a string it specifies the threshold
%                   above which labels have to be put. For example
%                   selunit='2.6' means that the text labels are written
%                   only for the units which have in at least one step of
%                   the search a value of the scaled residual greater than
%                   2.6 in absolute value.
%                   If it is a numeric vector it
%                   contains the list of the units for which it is
%                   necessary to put the text labels.
%                   The default value of
%                   selunit is string '2.5' if the input is a structure
%                   else it is an empty value if if the input is matrices y
%                   and X.
%                   Example - 'selunit','3'
%                   Data Types - numeric or character
%
%
%       The options which follow can only be used if the input is a
%       structure which contains information about the fwd search (i.e. the
%       two fields RES and Un)
%
%   datatooltip :   personalized tooltip. Empty value or structure. 
%                    The default is datatooltip=''
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
%                   Data Types - char
%       label   :   uni labels. Cell. Cell containing the labels of the units (optional
%                   argument used when datatooltip=1. If this field is not
%                   present labels row1, ..., rown will be automatically
%                   created and included in the pop up datatooltip window)
%                   Example - 'label',{'Row 1' 'Row 2' 'Row 3'}
%                   Data Types - cell
%     databrush :   empty value, scalar or cell.
%                   DATABRUSH IS AN EMPTY VALUE
%                   If databrush is an empty value (default), no brushing
%                   is done.
%                   The activation of this option (databrush is a scalar or
%                   a cell) enables the user  to select a set of
%                   observations in the current plot and to see them
%                   highlighted in the resfwdplot, i.e. the plot of the
%                   trajectories of all observations, grouped according
%                   to the selection(s) done by brushing. If the resfwdplot
%                   does not exist it is automatically created.
%                   In addition, brushed units can be highlighted in the
%                   other following plots (only if they are already open):
%                   - minimum deletion residual plot;
%                   - monitoring leverage plot;
%                   - maximum studentized residual;
%                   - s^2 and R^2;
%                   - Cook distance and modified Cook distance;
%                   - deletion t statistics.
%                   Remark. The window style of the other figures is set
%                   equal to that which contains the monitoring residual
%                   plot. In other words, if the scatterplot matrix plot
%                   is docked all the other figures will be docked too
%                   DATABRUSH IS A SCALAR
%                   If databrush is a scalar the default selection tool is
%                   a rectangular brush and it is possible to brush only
%                   once (that is persist='')
%                   DATABRUSH IS A CELL
%                   If databrush is a cell, it is possible to use all
%                   optional arguments of function selectdataFS.m and the
%                   following optional argument:
%                   - persist. Persist is an empty value or a scalar
%                     containing the strings 'on' or 'off'.
%                     The default value of persist is '', that is brushing
%                     is allowed only once.
%                     If persist is 'on' or 'off' brushing can be done as
%                     many time as the user requires.
%                     If persist='on' then the unit(s) currently brushed
%                     are added to those previously brushed. It is
%                     possible, every time a new brushing is done, to use a
%                     different color for the brushed units.
%                     If persist='off' every time a new brush is performed
%                     units previously brushed are removed.
%                   - bivarfit. This option is to add one or more least
%                     square lines to the plots of y|X, depending on the
%                     selected groups.
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
%                   - multivarfit. If this option is '1', we add to each scatter
%                     plot of y|X a line based on the fitted hyperplane
%                     coefficients. The line added to the scatter plot y|Xi
%                     is mean(y)+Ci*Xi, being Ci the coefficient of Xi.
%                     The default value of multivarfit is '', i.e. no line is
%                     added.
%                   - labeladd. If this option is '1', we label the units
%                     of the last selected group with the unit row index in
%                     matrices X and y. The default value is labeladd='',
%                     i.e. no label is added.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%
%       The options which follow work in connection with previous option
%       databrush and produce their effect on the monitoring residuals plot
%
%       subsize :   x axis control. Numeric vector. Numeric vector
%                   containing the subset size with length
%                   equal to the number of columns of matrix residuals.
%                   If it is not specified it will be set equal to
%                   size(residuals,1)-size(residuals,2)+1:size(residuals,1)
%                   Example - 'subsize',10:100
%                   Data Types - single | double
%       selstep :   text in selected steps. Numeric vector. Numeric vector
%                   which specifies for which steps of the
%                   forward search textlabels are added in the monitoring
%                   residual plot after a brushing action in the yXplot.
%                   The default is to write the labels at the initial and
%                   final step. The default is selstep=[m0 n] where m0 and
%                   n are respectively the first and final step of the
%                   search.
%                   Example - 'selstep',100
%                   Data Types - single | double
%
%  Output:
%
%        plot1  :   handle to the figure. Graphic handle. Handle to the
%        yXplot.
%
% See also: spmplot, mdrplot, fanplot, resfwdplot
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('yXplot')">Link to the help function</a>
% Last modified 14-06-2016

% Examples:

%{
    %First example of yXplot.
    % In the first example as input there are two matrices
    % y and X respectively
    % A simple yX plot is created
    n=100;
    p=3;
    X=randn(n,p);
    y=100+randn(n,1);
    % Example of the use of function yXplot with all the default options
    yXplot(y,X);
%}
%
%
%
%{
    % Example of the use of function yXplot with option group.
    n=100;
    p=3;
    X=randn(n,p);
    y=100+randn(n,1);
    sel=51:100;
    y(sel)=y(sel)+2;
    group=ones(n,1);
    group(sel)=2;
    yXplot(y,X,'group',group);
%}
%
%{
    % Example of option selunit.
    % Example of the use of function yXplot putting the text for the units
    % which have a value of y smaller than 98 and greater than 102.
    yXplot(y,X,'selunit',{'98' '102'});
%}
%
%{
    % yXplot with personalized labelling.
    % Example of the use of function yXplot putting the text for the units
    % which have a value of y smaller than 1% percentile and greater than
    % 99% percentile of y
    selth={num2str(prctile(y,1)) num2str(prctile(y,99))};
    yXplot(y,X,'selunit',selth);
%}

%{
    % Example of use of yXplot when input is a structure.
    % In the following example the input is a strucure which also contains
    % information about the forward search.
    [out]=LXS(y,X,'nsamp',1000);
    [out]=FSReda(y,X,out.bs);
    % Example of the use of function yXplot with all the default options
    yXplot(out);
%}
%
%   In all the examples which follow, we assume that the input structure
%   contains information about the fwd search (i.e. contains the two fields
%   RES and Un)
%
%{
    % Interactive_example
    %   Example of the use of options selunit and selstep.
    yXplot(out,'selunit',[2 5 20 23 35 45],'selstep',[20 22 27 36],...
            'databrush',{'persist','off','selectionmode' 'Rect'})
    %   produces a resfwdplot in which labels are put for units
    %   [2 5 20 23 35 45] in steps [20 22 27 36] of the search
%}
%
%{
    % Interactive_example  
    %   Example of the use of options selstep, selunit, selunitbold and
    %   selunitcolor.
   yXplot(out,'selstep',[40 21 80],'selunit','1.5',...
           'databrush',{'persist','off','selectionmode' 'Rect'})

    %   produces a yXplot plot in which labels are put for units
    %   which have a residual greater and 1.5. When a set of units is brushed in the yXplot
    %   in the monitoring residuals plot the labels are added in steps selsteps
%}
%
%{
    % Interactive_example
    %   Example of the use of option selunit (notice that in this
    %   case selunit is a cell array of strings.
    yXplot(out,'selunit',{'-3';'2'},...
            'databrush',{'selectionmode' 'Rect'});
    %   highlight only the trajectories which in at least one step of the
    %   search had a value smaller than -3 or greater than 2 and label
    %   them at the beginning and at the end of the search.
%}
%
%{
    %   Example of the use of option datatooltip.
   yXplot(out,'datatooltip',1);
    %   gives the possibility of clicking on the different points and have
    %   information about the unit selected, the step of entry into the
    %   subset and the associated label
%}
%{
    % Interactive_example
    %   Example of the use of option databrush
    %   (brushing is done only once using a rectangular selection tool).
       yXplot(out,'databrush',1)
    %   An equivalent statement is
       yXplot(out,'databrush',{'selectionmode' 'Rect'})
%}
%
%{
    % Interactive_example
    %   Example of the use of brush using a rectangular selection tool and
    %   a cyan colour.
   yXplot(out,'databrush',{'selectionmode' 'Rect' 'FlagColor' 'c'})
%}
%
%{
    % Interactive_example
    %  Example of the use of brush using multiple selection circular tools.
    yXplot(out,'databrush',{'selectionmode' 'Brush'})
%}
%
%{
    % Interactive_example
    %   Example of the use of brush using lasso selection tool and fleur pointer.
    yXplot(out,'databrush',{'selectionmode' 'lasso','Pointer','fleur'})
%}
%
%{
    % Interactive_example
    %   Example of the use of rectangular brush. Superimposed labels for
    %   the brushed units and persistent labels in the plot which has been
    %   brushed
    yXplot(out,'databrush',{'selectionmode' 'Rect' 'Label' 'on' 'RemoveLabels' 'off'})
%}
%
%   All previuos examples used a non persistent brushing (that is brushing
%   could be done only once). The examples below use persistent brushing
%   (that is brushing can be done multiple times)
%
%   Example of the use of persistent non cumulative brush. Every time a
%   brushing action is performed previous highlightments are removed
%{
    % Interactive_example
    % Example of persistent cumulative brushing (with persist off).
    %   Every time a brushing action is performed
    %   current highlightments replace previous highlightments
    yXplot(out,'databrush',{'selectionmode','Rect','persist' 'off' ...
                            'Label' 'on' 'RemoveLabels' 'off'})
%}

%{
    % Interactive_example
    % Example of persistent cumulative brushing (with persist on).
    %   Every time a brushing action is performed
    %   current highlightments are added to previous highlightments
    yXplot(out,'databrush',{'selectionmode','Rect','persist' 'on' ...
                            'Label' 'off' 'RemoveLabels' 'on'})
%}


%{
    % Interactive_example
    % Example of persistent cumulative brushing (with persist on and labeladd '1').
    %   Now option labeladd '1'. In this case the row numbers of the
    %   selected units is displayed in the monitoring residuals plot
    yXplot(out,'databrush',{'selectionmode','Rect','persist' 'on' ...
                            'labeladd' '1'})
%}

%% Beginning of code
% (common with resfwdplot)
%
% Close existing yXplot and resfwdplot
% if isempty(findobj('type','figure','Tag','pl_fan'))
%     close(findobj('type','figure','Tag','pl_yX'));
%     close(findobj('type','figure','Tag','pl_resfwd'));
% end

if nargin<1
    error('FSDA:yXplot:missingInputs','A required input argument is missing.')
end

% Check if the first argument is a structure or not
if ~isstruct(y)
    [n,p]=size(X);
    onlyyX=1;
else
    % If the first argument is a structure and the number of supplied
    % arguments is greater than 1 it is necessary to add to varargin the
    % second input (which is X)
    if nargin>1
        varargin=[X varargin];
    end
    out=y;
    % Retrieve y and X from the input structure
    y=out.y;
    X=out.X;
    [n,p]=size(X);
    
    if find(strcmp('RES',fieldnames(out)))
        % The number of rows of matrix residual is associated with the number of
        % units. The number of columns is associated with the number of steps of
        % the fwd search.
        residuals=out.RES;
        [n,nsteps]=size(residuals);
        onlyyX=0;
    else
        onlyyX=1;
    end
end
% seq= column vector containing the sequence 1 to n
seq= (1:n)';

% numtext= a cell of strings used to labels the units with their position
% in the dataset.
numtext=cellstr(num2str(seq,'%d'));

% Initialize line width
linewidthStd = 0.5;

if onlyyX==0
    % x= vector which contains the subset size (numbers on the x axis)
    x=(n-nsteps+1):n;
    
    % selthdef= threshold used to decide which residuals are labelled in the resfwdplot.
    % laby= label used for the y-axis of the resfwdplot.
    labx='Subset size m';
    if min(min(residuals))<0
        laby='Scaled residuals';
        selthdef=num2str(2.5);
    else
        laby='Squared scaled residuals';
        selthdef=num2str(2.5^2);
    end
    
    % maximum and minimum residual along the search for each unit
    selmax=max(residuals,[],2);
    selmin=min(residuals,[],2);
else
    x=1;
    % The default is not to add textlabels to any unit
    selthdef='';
end

%% User options
one=ones(n,1);
options= struct('subsize',x,'selstep',x([1 end]),'selunit',selthdef,...
    'xlim','','ylim','','tag','pl_yX',...
    'datatooltip',0,'label','','databrush','','nameX','','namey','','group',one);

if nargin>2 || (nargin>1 && onlyyX==0)
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:yXplot:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

group=options.group;

% labeladd option
d=find(strcmp('labeladd',options.databrush));
if d>0
    labeladd=options.databrush(d+1);
    labeladd=labeladd{1};
    % This option must be removed from cell options.databrush because it is
    % not a valid option for the function selectdataFS.
    options.databrush(d:d+1)=[];
else
    labeladd='';
end

% bivarfit option
d=find(strcmp('bivarfit',options.databrush));
if d>0
    bivarfit=options.databrush(d+1);
    bivarfit=bivarfit{1};
    % This option must be removed from cell options.databrush because it is
    % not a valid option for the function selectdataFS.
    options.databrush(d:d+1)=[];
else
    bivarfit='';
end

% multivarfit option
d=find(strcmp('multivarfit',options.databrush));
if d>0
    multivarfit=options.databrush(d+1);
    multivarfit=multivarfit{1};
    % This option must be removed from cell options.databrush because it is
    % not a valid option for the function selectdataFS.
    options.databrush(d:d+1)=[];
else
    multivarfit='';
end

% persist option
d=find(strcmp('persist',options.databrush));
if d>0
    persist=options.databrush(d+1);
    % This option must be removed from cell options.databrush because it is
    % not a valid option for the function selectdataFS.
    options.databrush(d:d+1)=[];
    
    ColorOrd=[1 0 0;0 1 1; 1 0 1; 1 1 0; 0 0 0; 0 1 0; 0 0 1];
    ColorOrd=repmat(ColorOrd,4,1);
else
    persist='';
    ColorOrd=[1 0 0];
end

% FlagColor option
% Initialize colors: default colors are blue (unbrushed unit)
% and red (brushed units)
d=find(strcmp('FlagColor',options.databrush));
if d>0
    flagcol=options.databrush{d+1};
    clr=['b' flagcol 'cmykgbrcmykg'];
else
    clr='brcmykgbrcmykgbrcmykg';
    
end

% selstep option

% extract the vector associated with the subset size (x)
x=options.subsize;
% and check if the choice of selsteps is valid
steps=options.selstep;
if max(steps)>n
    mess=sprintf(['One of the steps which has beeen chosen is greater than n. \n',...
        'It is deleted.']);
    fprintf('%s\n',mess);
    steps=steps(steps<=n);
end
if min(steps)<x(1)
    mess=sprintf(['One of the steps which has beeen chosen is smaller than m0. \n',...
        'It is deleted.']);
    fprintf('%s\n',mess);
    steps=steps(steps>=x(1));
end

% selunit option

% units= the list of the units which must be labelled.
% It can be a cell array of strings (defining lower and upper threhold),
% a string (defining just one threshold) or a numeric vector.

units=options.selunit;

if onlyyX
    
    if ~isempty(options.databrush)
        disp('It is not possible to use option databrush without supplying structure out produced by FSReda')
        return
    end
    if iscellstr(units)
        selunit=str2double(units);
        units=seq(y>selunit(2) | y<selunit(1));
    elseif ischar(units)==1
        % if units is a character than the user has specified a threshold
        % convert character to numeric
        thresh=str2double(units);
        % Label the units whose maximum residual along the search is greater
        % than the required threshold or smalleR than -threshold
        units=seq(y>thresh | y<-thresh);
    else
        % Some checks on minimum and maximum of vector units
        if max(units)>n
            mess=sprintf(['One of the units which has beeen chosen is greater than n. \n',...
                'It is deleted.']);
            fprintf('%s\n',mess);
            units=units(units<=n);
        end
        if min(units)<1
            mess=sprintf(['One of the units which has beeen chosen is smaller than 1. \n',...
                'It is deleted.']);
            fprintf('%s\n',mess);
            units=units(units>0);
        end
    end
    
else
    if iscellstr(units)
        selunit=str2double(units);
        selmax=max(residuals,[],2);
        selmin=min(residuals,[],2);
        units=seq(selmax>selunit(2) | selmin<selunit(1));
    elseif ischar(units)==1
        % if units is a character than the user has specified a threshold
        % convert character to numeric
        thresh=str2double(units);
        % Label the units whose maximum residual along the search is greater
        % than the required threshold or smalleR than -threshold
        units=seq(selmax>thresh | selmin<-thresh);
    else
        % Some checks on minimum and maximum of vector units
        if max(units)>n
            mess=sprintf(['One of the units which has beeen chosen is greater than n. \n',...
                'It is deleted.']);
            fprintf('%s\n',mess);
            units=units(units<=n);
        end
        if min(units)<1
            mess=sprintf(['One of the units which has beeen chosen is smaller than 1. \n',...
                'It is deleted.']);
            fprintf('%s\n',mess);
            units=units(units>0);
        end
    end
end

% lunits = number of units which must be labelled
lunits=length(units);
% lsteps = number of steps for which it is necessary to add the labels
lsteps=length(steps);
lall=lunits*lsteps;

%% Display the yXplot

% Create a figure to host the gplotmatrix or clear the existing one
h=findobj('-depth',1,'tag',options.tag);
if (~isempty(h))
    clf(h);
    figure(h);
    axes;
else
    h=figure;
end
set(h,'Name', 'yXplot: plot of y against each column of matrix X', 'NumberTitle', 'off');
hold('all');

% Set default value for potential groups of selected units
unigroup=unique(group);
styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};

% Check if X includes the constant term for the intercept.

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

% Display the initial gplotmatrix
[H,AX,BigAx] = gplotmatrix(Xsel,y,group,clr(unigroup),char(styp{unigroup}),[],'on',[],nameX,namey);

%[H,AX,BigAx] = gplotmatrix(Xsel,Y,group,clr(unigroup),charsym,siz,doleg,'hist',nameY,nameY);


% default legenda
set(H,'DisplayName','Units');

% if ~isempty(units)
%     for i = 1:length(AX)
%         set(gcf,'CurrentAxes',AX(i));
%         xlimits = get(AX(i),'Xlim'); ylimits = get(AX(i),'Ylim');
%         dx = (xlimits(2)-xlimits(1))*0.01*length(AX); dy = (ylimits(2)-ylimits(1))*0.01*length(AX)/2; % displacement
%         text(Xsel(units,i)+dx,y(units)+dy,numtext(units),'HorizontalAlignment', 'Left');
%     end
% end


for i = 1:length(AX)
    set(gcf,'CurrentAxes',AX(i));
    % If the user has specified the min and max for y limit
    
    if ~isempty(options.ylim)
        ylim(AX(i),options.ylim);
    end
    % if selunit option has been defined label the units.
    
    if ~isempty(units)
        xlimits = get(AX(i),'Xlim'); ylimits = get(AX(i),'Ylim');
        dx = (xlimits(2)-xlimits(1))*0.01*length(AX); dy = (ylimits(2)-ylimits(1))*0.01*length(AX)/2; % displacement
        text(Xsel(units,i)+dx,y(units)+dy,numtext(units),'HorizontalAlignment', 'Left');
    end
end



%%  Prepare the yXplot for brushing

% the handle of the figure including the gplotmatrix
% (i.e. the closest ancestor of BigAx).
fig = ancestor(BigAx,'figure');
plot1 = fig;
%plot1=gcf;

%when the graphics objects are axes children, MATLAB clears the axes and
%resets some of their properties to default values. In order to add new
%graphics objects without clearing or resetting the current figure we set
%the figure and axes NextPlot properties to 'add'.
set(fig,'NextPlot','add');
set(AX,'NextPlot','add');

styp=repmat(styp,ceil(n/13),1);

% displays the boundary of the current axes.
box on

% set the specified tag in the current plot
set(gcf,'tag',options.tag)

% set the options.datatooltip (enable/disable interactive data cursor mode)
if options.datatooltip
    hdt = datacursormode;
    set(hdt,'Enable','on');
    
    % If options.datatooltip is not a struct then use our default options
    if ~isstruct(options.datatooltip)
        set(hdt,'DisplayStyle','window','SnapToDataVertex','on');
    else
        % options.databrush contains a structure where the user can set the
        % properties of the data cursor
        set(hdt,options.datatooltip);
    end
    % Declare a custom datatooltip update function to display additional
    % information about the selected unit
    set(hdt,'UpdateFcn',{@yXplotLbl,out})
end

%% Brush mode (call to function selectdataFS)
if ~isempty(options.databrush) || iscell(options.databrush)
    if isscalar(options.databrush)
        sele={'selectionmode' 'Rect' 'Ignore' findobj(gcf,'tag','env') };
    else
        sele={options.databrush{:} 'Ignore' findobj(gcf,'tag','env')}; %#ok<*CCAT>
    end
    
    sele={sele{:} 'Tag' options.tag};
    
    % group = VECTOR WHICH WILL CONTAIN THE IDENTIFIER OF EACH GROUP
    % e.g. group(14)=3 means that unit 14 was selected at the third brush
    group=ones(n,1);
    
    % some local variables
    i=0; but=0; brushcum=[]; ij=1;
    
    % loop brushing
    while but<=1
        i=i+1;
        figure(plot1);
        
        % Remark: function selectdataFS cannot be used on the current figure if
        % the "selection mode" or the "zoom tool" are on. Setting the
        % plotedit mode initially to on and then to off, has the effect to
        % deselect plotedit mode.
        plotedit on
        plotedit off
        
        if strcmp(persist,'off')
            % Remove from the current plot what has alredy been selected
            a=findobj(gcf,'Tag','selected');
            delete(a);
        elseif strcmp(persist,'on')
            % add to cell sele option FlagColor (color of selection) and
            % FlagMarker (symbol to be used for selection)
            chkexist=find(strcmp('FlagColor',sele)==1);
            if ~isempty(chkexist)
                sele=sele(1:chkexist(1)-1);
            end
            sele={sele{:} 'FlagColor' ColorOrd(ij,:)  'FlagMarker' char(styp(ij+1))};
            % add to sele the name of the file which calls selectdataFS in
            % order to avoid that, when persist = on, the user can start
            % brushing from a plot different to the yXplot
        end
        
        %% - select an area in the yXplot
        
        % A function to be executed on figure close is set
        set(gcf,'CloseRequestFcn',@closereqFS);
        
        disp('Select a region to brush in the yXplot');
        disp('Left mouse button picks points.');
        ss=waitforbuttonpressFS;
        
        set(get(0,'CurrentFigure'),'CloseRequestFcn','closereq'); %
        
        if ss==1
            but=2;
        else
            but=1;
        end
        
        %         button = get(fig, 'SelectionType');
        %         if strcmp(button,'open'), but = 1;
        %         elseif strcmp(button,'normal'), but = 1;
        %         elseif strcmp(button,'extend'), but = 2;
        %         % elseif strcmp(button,'alt'), % but = 3;
        %         else error('MATLAB:ginput:InvalidSelection', 'Invalid mouse selection.')
        %         end
        
        if but==1 && find(AX==gca) > 0 %if the User has made a selection in one of the scatterplots
            
            %% - identify the axes of the area
            
            % fig is the closest ancestor figure of BigAx. Set property
            % CurrentAxes equal to the axes of the scatterplot that we have
            % selected
            set(fig,'CurrentAxes',gca);
            %indice is a scalar identyfing the selected axes
            indice=find(AX==gca);
            %otherAxes is the list of the not selected scatterplot axes
            otherAxes = AX;
            otherAxes(indice)=[];
            
            %During the selection, not selected axes must have properties HandleVisibility and
            %HitTest set to off.
            set(otherAxes,'HandleVisibility','off');
            set(otherAxes,'HitTest','off');
            
            %% - call selectdataFS
            
            if ij>1
                legend(hLegend(length(AX)),'hide');
                %               set(hLegend(length(AX)),'Location', 'EastOutside');
            end;
            set(fig,'UserData',num2cell(seq));
            [pl,xselect,yselect] = selectdataFS(sele{:}, 'Label', 'off');
            
            % exit from yXplot if the yXplot figure was closed before selection
            if ~isempty(pl) && isnumeric(pl) && (min(pl) == -999)
                return
            end
            
            if ij>1; set(hLegend(length(AX)),'Location', 'Best'); end;
            
            %When the selection has been completed, axes properties
            %'HandleVisibility' and 'HitTest' must be set to on for an eventual
            %possible future selection.
            set(otherAxes,'HandleVisibility','on');
            set(otherAxes,'HitTest','on');
            
            %% - identify the selected units
            if iscell(xselect)
                %xselect is a cell if more than one selection have been
                %already done. If there is no intersection between the
                %current and the previous selection, the list of the X of
                %the current selection is the last cell of xselect, since
                %the other cells are empty. Otherwise, if one (more)
                %unit(s), that was (were) previously selected, is (are)
                %selected again, it is the conatenation of the last cell
                %with the previous cell where the unit (units) was (were)
                %previously selected.
                for z=1:length(xselect)
                    if z==1
                        xselect_all=xselect{z};
                        yselect_all=yselect{z};
                    else
                        xselect_all=cat(1,xselect_all,xselect{z} ) ;
                        yselect_all=cat(1,yselect_all,yselect{z} ) ;
                    end
                end
            else    %xselect is a vector if there is only one selection.
                xselect_all=xselect;
                yselect_all=yselect;
            end
            %coi= list of row identifiers of selected units.
            % Initialize coi assuming n rows (just in case there are
            % duplicate units)
            coi=NaN(n,1);
            jk=1;
            
            for k=1:length(xselect_all)
                %uno=vector of zeros and ones of the same length of X and
                %y. The ones are reported in the rows where the X values
                %are equal to the X selected
                uno=Xsel(:,indice)==(xselect_all(k));
                %due=same idea of uno: vector of zeros and ones of the same
                %length of X and y. The ones are reported in the rows where
                %the y values are equal to the y selected
                due=y(:)==(yselect_all(k));
                %find observations that have satisfied both (AND)
                %conditions uno and due, i.e. observations with both X and
                %y equal to the X and y selected.
                tre=uno + due;
                %coi=list of selected units with that particular
                %combination of X and y. Note that in presence of
                %duplicate units listbra can have have more than one row
                listbra=find(tre==2);
                % allocate inside vector coi brushed units
                coi(jk:jk+length(listbra)-1)=listbra;
                jk=jk+length(listbra);
            end
            % Resize coi with the number of selected units. 
            % jk-1 = number of selected units
            coi=coi(1:jk-1);
            
            %pl=vector with n rows. Ones represent observations selected at
            %the current iteration, zeros all the others.
            pl=zeros(length(Xsel),1);
            pl(coi)=1;
            if ~isempty(pl)
                % sel = vector which contains the list of units not selected at
                % the current iteration.
                sel=seq(pl<1);
                
                % nbrush = vector which contains the list of the units selected
                % at the current iteration.
                nbrush=setdiff(seq,sel);
                disp('Brushed units are');
                disp([nbrush out.y(nbrush) out.X(nbrush,:)]);
            else
                disp('Wrong selection: Try again');
                disp('Select a region to brush in the monitoring residual plot');
                figure(plot1);
                nbrush='';
            end
            
            if ~isempty(pl)
                disp('Brushed units are');
                disp([nbrush out.y(nbrush) Xsel(nbrush,indice)]);
            else
                disp('Wrong selection: Try again');
                disp('Select a region to brush in the monitoring residual plot');
                figure(plot1);
                nbrush='';
            end
            
            %% For each brushing operation, do the following:
            if ~isempty(nbrush)
                
                %% - Update the yXplot and keep it always in foreground
                
                %brushcum is the list of selected observations in all
                %iterations if persist=on or the list of selected
                %observations in the current iteration if persist=off
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
                %open the figure containing the scatterplot matrix
                figure(plot1);
                
                % set in this figure, containing the gplotmatrix, the
                % property currentAxes equal to the axes of the scatterplot
                % in which the selection has been done.
                set(plot1,'CurrentAxes',AX(indice));
                
                %replace the gplotmatrix which has come out from
                %selectdataFS
                %(the selection was done ONLY on one scatterplot) with a
                %new gplotmatrix which takes into consideration the
                %selection coming out from selectdataFS, in all scatterplots
                [H,AX,BigAx] = gplotmatrix(Xsel,y,group,clr(unigroup),char(styp{unigroup}),[],'on',[],nameX,namey);
                
                % Set the legenda properties of the gplotmatrix
                if nbrush>0
                    set(H(:,:,1),'DisplayName','Unbrushed units');
                    for brugrp = 2:length(unigroup)
                        set(H(:,:,brugrp),'DisplayName',['Brushed units ' num2str(brugrp-1)]);
                    end
                else
                    set(H,'DisplayName','Units');
                end
                
                %fig is the handle of the closest ancestor figure of BigAx.
                fig = ancestor(BigAx,'figure');
                
                %In order to add new graphics objects without clearing or
                %resetting the current figure we set the figure and axes
                %NextPlot properties to 'add'.
                set(fig,'NextPlot','add');
                set(AX,'NextPlot','add');
                
                % Now, for each plot of y|X do:
                for i = 1:length(AX)
                    % Make the axes of the scatterplot identified by the handle
                    % AX(i) the current axes.
                    set(fig,'CurrentAxes',AX(i));
                    % Remark: axes(AX(i)) would also do the job, but would restack
                    % the axes above all other axes in the figure.
                    
                    % Fit least square line(s) to the scatterplot AX(i).
                    switch bivarfit
                        case ''
                            %do nothing: no line is fit.
                        case '0'
                            h=olsline(0); % fit a line to each group
                            if length(h)==1
                                set(h,'DisplayName','fit on unbrushed units of y|Xi');
                            else
                                for brugrp = 1:size(h,1)-1
                                    set(h(brugrp),'DisplayName',['fit on brushed units ' num2str(size(h,1)-brugrp)]);
                                end
                                %set(h(2),'DisplayName','fit on unbrushed units of y|Xi');
                                set(h(size(h,1)),'DisplayName','fit on unbrushed units of y|Xi');
                            end
                        case '1'
                            h=olsline;    % fit 1 line to all data, regardless the groups
                            set(h,'DisplayName','fit on all units of y|Xi');
                        case '2'
                            h=olsline;    % fit a line to all data
                            set(h,'DisplayName','fit on all units of y|Xi');
                            h=olsline(size(H,3)); % fit a line to the unselected data, i.e. the last group among the handles H
                            set(h,'DisplayName','fit on unbrushed units of y|Xi');
                        otherwise
                            if strncmp('i',bivarfit,1)
                                token = strtok(bivarfit, 'i');
                                selgroup=str2double(token);
                                if ~isnan(selgroup) && selgroup > 0.5 && selgroup <= size(H,3)
                                    h=olsline(size(H,3)-round(selgroup)+1); % fit one group only: the one with index round(bivarfit)
                                    set(h,'DisplayName','fit on a group of y|Xi');
                                else
                                    %do nothing: no line is fit.
                                end
                            else
                                error('FSDA:yXplot:WrongBivarfit','''Valid values for option ''bivarfit'' are: '''', ''0'', ''1'', ''2'', ''i1'', ''i2'', ... , ''ig'', ... being ''g'' the index of a selected group.')
                                %do nothing: no line is fit
                            end
                    end
                    
                    % Add the labels of the last selected group.
                    if strcmp('1',labeladd)
                        xlimits = get(AX(i),'Xlim'); ylimits = get(AX(i),'Ylim');
                        dx = (xlimits(2)-xlimits(1))*0.01*length(AX); dy = (ylimits(2)-ylimits(1))*0.01*length(AX)/2; % displacement
                        text(Xsel(brushcum,i)+dx,y(brushcum)+dy,numtext(brushcum),'HorizontalAlignment', 'Left');
                    end
                    % Add the (fix) labels eventually set in selunit.
                    if ~isempty(units)
                        xlimits = get(AX(i),'Xlim'); ylimits = get(AX(i),'Ylim');
                        dx = (xlimits(2)-xlimits(1))*0.01*length(AX); dy = (ylimits(2)-ylimits(1))*0.01*length(AX)/2; % displacement
                        text(Xsel(units,i)+dx,y(units)+dy,numtext(units),'HorizontalAlignment', 'Left');
                    end
                    
                    % Add to each scatterplot AX(i) the line(s) based on the hyperplane
                    % fit to y|X.
                    switch multivarfit
                        case {'1' , '2'}
                            %int = out.beta(intcolumn,1); % the intercept value
                            coef =regress(y,X);
                            
                            % indcoef = vector which contains the indexes of columns of X
                            % except that which is about to be plotted
                            if intercept==1
                                indcoef=setdiff(1:p,i+1);
                            else
                                indcoef=setdiff(1:p,i);
                            end
                            % Find the mean value of the columns of matrix X using all
                            % the units
                            meaot=mean(X(:,indcoef))*coef(indcoef);
                            
                            coef(intcolumn)=[]; % the other coefficients
                            xlimits = get(AX(i),'Xlim');
                            hline1 = line(xlimits , meaot + coef(i).*xlimits);
                            datacolor = [.3 .3 .3]; %Dark grey;  [1 .62 .40] == %Copper
                            set(hline1,'Color',datacolor,'LineWidth',2,...
                                'DisplayName','fit on all units of y|X');
                            if strcmp(multivarfit, '2')
                                % Notice that group==1 is the group of the
                                % unselected units
                                Xgood = X(logical(group==1),:);
                                ygood = y(logical(group==1));
                                
                                coef = regress(ygood,Xgood);
                                
                                % Find the mean of all the other variables
                                % considering just good units
                                meaot=mean(Xgood(:,indcoef))*coef(indcoef);
                                
                                coef(intcolumn)=[];
                                hline2 = line(xlimits , meaot  + coef(i).*xlimits);
                                datacolor = get(H(1,i,1),'Color');
                                %datacolor = [1 .40 .80]; %Copper verso rosso
                                set(hline2,'Color',datacolor,'LineWidth',2,...
                                    'DisplayName','fit on unbrushed units of y|X');
                            end
                        otherwise
                            %do nothing
                    end
                    
                end
                
                % Now update the legends and make them clickable.
                
                hLines = findobj(AX(1), 'type', 'line');
                hLegend = zeros(size(AX));
                for iAxes = 1:length(AX)
                    % Find all lines in current axes and get their DisplayName
                    hLines = findobj(AX(iAxes), 'type', 'line');
                    eLegend = cell(length(hLines), 1);
                    for iLines = 1:length(hLines)
                        eLegend{iLines} = get(hLines(iLines), 'DisplayName');
                    end
                    %%hLegend(iAxes) = clickableLegend(hLines, eLegend{:});
                    legend_h = legend(hLines);
                    if iAxes < length(AX)
                        %                         if verMatlab
                        %                             legend(legend_h,'hide');
                        %                         else
                        %                             legend_h.Visible='off';
                        %                         end
                        set(legend_h,'Visible','off')
                        
                    end
                end
                hLegend(iAxes) = clickableMultiLegend(hLines, eLegend{:});
                
                %% - display the resfwdplot with the corresponding groups of trajectories highlighted.
                
                % creates the resfwdplot
                fres=findobj('-depth',1,'Tag','data_res');
                if isempty(fres)
                    figure;
                    % The window style of the new figure which has been
                    % created is set equal to that which contained the
                    % yXplot
                    set(gcf,'WindowStyle',get(plot1,'WindowStyle'));
                    set(gcf,'Name','resfwdplot','NumberTitle', 'off');
                    
                    % store the LineWidth property of the selected trajectories
                    
                    plot(x,residuals,'Tag','data_res','Color','b','LineWidth',0.5);
                    
                    set(gcf,'tag','data_res');
                    hold('off')
                    % control minimum and maximum for x and y axis
                    if ~isempty(options.xlim)
                        xlim(options.xlim);
                    end
                    
                    if ~isempty(options.ylim)
                        ylim(options.ylim);
                    end
                    
                    % displays the boundary of the current axes.
                    box on
                end
                
                
                %check if the resfwdplot is already open
                hh=findobj('-depth',1,'tag','data_res');
                figure(hh);
                set(gcf,'WindowStyle',get(plot1,'WindowStyle'));
                
                if strcmp('off',persist)
                    %if persist=off the previous selection must be replaced
                    set(gcf,'NextPlot','replace');
                    plot(x,residuals,'Tag','data_res','Color','b');
                    text(reshape(repmat(steps,lunits,1),lall,1),reshape(residuals(units,steps-x(1)+1),lall,1),reshape(repmat(numtext(units),1,lsteps),lall,1));
                    set(gcf,'tag','data_res');
                    set(gcf,'Name','resfwdplot');
                    set(gcf,'NextPlot','add');
                end
                hold on
                if strcmp('on',persist)
                    %find the handles of the residual trajectories
                    a=get(gcf,'Children');
                    aa=get(a,'Children');
                    %sort the handles of the residual trajectories
                    aa1=sort(aa);
                    %select the handles of the selected residual trajectories
                    aa1=aa1(nbrush);
                    %delete the selected trajectories
                    delete(aa1);
                end
                % plot on the residual plot previously created the same
                % tajectories that we have just deleted, but in another
                % color and with a width that is 2 points bigger than the
                % one of the standard trajectories.
                plot(x,residuals(nbrush,:),...
                    strcat(clr(unigroup(length(unigroup))),...
                    char(styp{unigroup(length(unigroup))})),...
                    'LineStyle','-','MarkerSize',0.3,'tag','data_res',...
                    'LineWidth',linewidthStd+2);
                
                %add labels, if necessary.
                if strcmp('1',labeladd)
                    if strcmp('off',persist)
                        text(reshape(repmat(steps,length(nbrush),1),length(nbrush)*length(steps),1),reshape(residuals(nbrush,steps-x(1)+1),length(nbrush)*length(steps),1),reshape(repmat(numtext(nbrush),1,length(steps)),length(nbrush)*length(steps),1));
                    end
                    if strcmp('on',persist)
                        text(reshape(repmat(steps,length(brushcum),1),length(brushcum)*length(steps),1),reshape(residuals(brushcum,steps-x(1)+1),length(brushcum)*length(steps),1),reshape(repmat(numtext(brushcum),1,length(steps)),length(brushcum)*length(steps),1));
                    end
                end
                %add (fix) labels eventually set in selunit.
                if ~isempty(units)
                    text(reshape(repmat(steps,length(units),1),length(units)*length(steps),1),reshape(residuals(units,steps-x(1)+1),length(units)*length(steps),1),reshape(repmat(numtext(units),1,length(steps)),length(units)*length(steps),1));
                end
                
                xlabel(labx);
                ylabel(laby);
                % Initialize vector selstesp. It will contains the steps in
                % which brushed units enter the search
                selsteps=zeros(n,11);
                ii=0;
                % Find steps where selected (brushed) units entered the
                % search vector nbrush contined the brushed steps
                % Now find the steps where the brushed units entered the
                % search
                Un=out.Un;
                for i=1:length(nbrush)
                    idx = find(Un(:,2:end) == nbrush(i));
                    
                    % Find the required row(s) of matrix Un
                    row = ind2sub(size(Un(:,2:end)),idx);
                    % if isempty(row) is true, it means that the selected
                    % unit entered the subset at step before Un(1,:), that
                    % is before the first step which has been monitored
                    if ~isempty(row)
                        % REMARK: mod(row,size(Un,1)) is used rather than
                        % Un(row,:) because the unit could have entered the
                        % subset together with other units (i.e. could be
                        % in a column  of matrix Un different from the
                        % second)
                        % Finally notice that the expression row-1e-12 is
                        % necessary otherwise when row is= size(Un,1) then
                        % mod(row-1e-12,size(Un,1))is equal to zero
                        selsteps(ii+1:ii+length(row),:)=Un(ceil(mod(row-1e-12,size(Un,1))),:);
                        ii=ii+length(row);
                    end
                end
                
                selsteps=sortrows(selsteps(1:ii,:),1);
                % m1 contains the indexes of the unique steps;
                [~, m1]=unique(selsteps(:,1));
                % WARNING: selsteps=unique(selsteps,'rows') does not seem
                % to work
                selsteps=selsteps(m1,:);
                
                disp('Steps of entry of brushed units');
                disp(selsteps);
                
                % Check if the figure containing minimum deletion residual
                % is open If it is, units are brushed in this plot too
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
                        % If set of values has already been highlighted in the mdr plot,
                        % remove it
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
                    
                    [~, ~, ib]=intersect(selsteps, xdata);
                    % Stack together x and y coordinates
                    xx=[xdata; ydata];
                    
                    % Just in case the first step of mdr is selected remove
                    % it because we also consider ib-1
                    ib=ib(ib>1);
                    % For each of the brushed units extract coordinates of
                    % mdr referred to the step before their entry and the
                    % step before
                    xxsel=xx(:,[ib-1 ib])';
                    % Sort all steps
                    xxselr=sortrows(xxsel,1);
                    % xxlim=length(nbrush);
                    xxlim=length(ib);
                    % Reshape previous matrix in such a way that the first
                    % length(nbrush) columns refer to the steps which have
                    % to be plotted and the remining columns refer to their
                    % corresponding values of mdr
                    xy=reshape(xxselr,2,2*xxlim);
                    % Add to the previous matrix a row of missing values
                    % This operation is necessary if the steps are not
                    % contiguous
                    xy=cat(1,xy,NaN*zeros(1,2*xxlim));
                    
                    % Reshape the set of x and y coordinates in two column
                    % vectors. Notice the NaN between the steps which are not
                    % consecutive
                    xcoord=reshape(xy(:,1:xxlim),3*xxlim,1);
                    ycoord=reshape(xy(:,xxlim+1:end),3*xxlim,1);
                    hold('on');
                    if strcmp('on',persist)
                        plot(gca,xcoord,ycoord,'LineWidth',4,'color',ColorOrd(ij,:),'tag','brush_mdr');
                    else
                        plot(gca,xcoord,ycoord,'LineWidth',4,'color',clr(ij),'tag','brush_mdr');
                    end
                    hold('off');
                end
                
                if strcmp('on',persist) || strcmp('off',persist)
                    
                    % variable ij is linked to the bridhing color
                    if strcmp('on',persist)
                        ij=ij+1;
                        % Set but=1 so that previous highlightments in
                        % other figures are not deleted
                        but=1;
                    end
                    
                    disp('Highlight the yXplot, then click on it to continue brushing or press a keyboard key to stop')
                    % The monitoring residual plot is highlighted again
                    % before continuing with waitforbuttonpress
                    figure(plot1);
                    set(gcf,'CloseRequestFcn',@closereqFS);
                    ss=waitforbuttonpressFS;
                    set(get(0,'CurrentFigure'),'CloseRequestFcn','closereq');
                    disp('------------------------')
                    if ss==1
                        but=2;
                    end
                else
                    but=2;
                end
            end
        end
    end
end
    function output_txt = yXplotLbl(~,event_obj,out)
        %% yXplotLbl provides information about the selected points
        %
        % Required input arguments:
        %
        %   obj     =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the sense that
        %               these arguments are automatically passed to the function when it executes.
        %       out =   a structure containing the following fields
        %       y   =   a vector containing the y values
        %       Un  =   a matrix containing the list of the units which entered the subset
        %               in each step of the search
        %      label=  (optional argument) if it is present it must be
        %               a cell array of strings containing the labels of
        %               the rows of the regression dataset
        %
        % Output:
        %
        %   output_txt=  Datatip text (string or string cell array) which informs
        %                about the unit(s) selected, its (their) X and y values
        %                and its (their) entry during the fwd search
        %
        % REMARK: this function is called by function yXplot
        %
        % References:
        %
        %   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
        %   Springer Verlag, New York.
        %
        % Written by FSDA team
        
        %% - identify the axes of the area
        
        % fig is the closest ancestor figure of BigAx. Set property CurrentAxes
        % equal to the axes of the scatterplot that we have selected
        set(fig,'CurrentAxes',gca);
        
        %indice is a scalar identyfing the selected axes
        indice=find(AX==gca);
        
        %otherAxes is the list of the not selected scatterplot axes
        otherAxes = AX;
        otherAxes(indice)=[];
        
        %During the selection, not selected axes must have properties
        %HandleVisibility and HitTest set to off.
        set(otherAxes,'HandleVisibility','off');
        set(otherAxes,'HitTest','off');
        
        pos = get(event_obj,'Position');
        
        % x and y, plot coordinates of the mouse
        x = pos(1); y = pos(2);
        
        % Find index to retrieve obs. name
        % Consider that find return the linear indexing of matrix xydata
        ord=out.y;
        idx = find(ord == y,1);
        
        % Linear indexing is transformed into normal indexing using function
        % ind2sub
        % row and column contain the column and row indexed of the observation
        % which has been selected with the mouse
        [row,~] = ind2sub(size(ord),idx);
        
        if isempty(row)
            output_txt{1}=['no residual has coordinates x,y' num2str(x) '' num2str(y)] ;
        else
            
            % If structure out does not contain labels for the rows then labels
            % row1....rown are added automatically
            if isempty(intersect('label',fieldnames(out)))
                out.label=cellstr(num2str((1:length(out.y))','row %d'));
            end
            
            % output_txt is what it is shown on the screen
            output_txt = {['y value equal to: ',num2str(y,4)]};
            
            % Add information about the step of the search which is under
            % investigation
            output_txt{end+1} = ['x value equal to:' num2str(x)];
            
            % Add information about the corresponding row label of what has been
            % selected
            output_txt{end+1} = ['Unit: ' cell2mat(out.label(row))];
            
            % Add information about the step in which the selected unit entered
            % the search
            Un=out.Un;
            idx = find(Un(:,2:end) == row,1);
            [row,~] = ind2sub(size(Un(:,2:end)),idx);
            output_txt{end+1} = ['Unit entered in step m=' num2str(Un(row,1))];
        end
        
        %When the selection has been completed, axes properties
        %'HandleVisibility' and 'HitTest' must be set to on for an eventual
        %possible future selection.
        set(otherAxes,'HandleVisibility','on');
        set(otherAxes,'HitTest','on');
        
    end
end
%FScategory:VIS-Reg