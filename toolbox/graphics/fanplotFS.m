function brushedUnits=fanplotFS(out,varargin)
%fanplotFS plots the fan plot for transformation in linear regression or deletion t stat
%
%<a href="matlab: docsearchFS('fanplotFS')">Link to the help function</a>
%
% Required input arguments:
%
%  out :  Data to plot. Structure. Structure containing the following fields
%     out.Score  =  (n-init) x length(la)+1 matrix:
%               1st col = fwd search index or bpd or eff
%               2nd col = value of the statistic in each step
%               of the fwd search or each bdp or each eff,
%               For example if for la(1);
%               ...;
%               last col  =  value of the score test in each step
%               of the fwd search for la(end).
%               Remark: note that out.Score can be replaced by any other
%               statistic monitored as function of breakdown point,
%               efficiency or subset szie
%               out.Tdel if
%               the input comes from routine FSRaddt.
%       out.la   =  vector containing the values of transformation
%               parameter lambda which have been used inside routine FSRfan
%               or the numbers associated to the columns of matrix out.X
%               for which deletion t stats have computed by routine
%               FSRaddt. Alternatively out.la can be a string array or cell
%               arry of characters containing the names of the variables
%               associated with the deletion t stats.
%       out.bs   =  matrix of size p x length(la) containing the units
%               forming the initial subset for each value of lambda.
%      out. Un   =  cell of size length(la). out.Un{i} is a (n-init) x 11
%               matrix which contains the unit(s) included in the subset
%               at each step of the fwd search (necessary only if option
%               datatooltip or databrush are not empty).
%         out.y  = a vector containing the response (necessary only if option
%              databrush is true).
%         out.X  = a matrix containing the explanatory variables (necessary
%              only if option databrush is not empty).
%      Data Types - struct
%
% Optional input arguments:
%
%
%      addxline : add to the plot a set of vertical lines.
%                 Numeric vector.
%                 It is possible to add to the overall plot or to each
%                 panel vertical line(s) with constant x-value. The default
%                 value is [], that is no vertical line is added.
%                   Example - 'addxline', [30 50 200]
%                   Data Types - double
%
%       conflev :   Confidence level. Scalar or vector. Confidence level
%                   for the bands (default is 0.99, that is we plot two
%                   horizontal lines in correspondence of value -2.58 and
%                   2.58).
%                   Example - 'conflev',[0.9 0.95 0.99]
%                   Data Types - double
%
%     flabstep  :   numeric vector which specifies the steps of
%                   where to put labels for the trajectories of the
%                   quantities which are monitored. The default is to put
%                   the labels in the final steps of the monitoring
%                   procedure. flabstep=[] means no label.
%                   Example - 'flabstep',[20 30 120]
%                   Data Types - double
%
%
%       FontSize:   Font size of the labels. Scalar. Scalar which controls
%                   the font size of the labels of the axes and of the
%                   labels inside the plot. Default value is 12.
%                   Example - 'FontSize',12
%                   Data Types - double
%
%    highlight : units to highlight in the fanplot plot. Vector or 2D array or empty
%               (default). If highlight is a vector it contains numbers
%               associated to the units whose entry step has to be shown in
%               trajectories of the fanplot. More specifically, the steps
%               in which the units inside highlight join the subset is
%               shown in red in the trajectories of the fanplot. The
%               default of  highlight is the empty vector []. If highlight
%               is a matrix it contains in column j the units whose entry
%               step has to be shown in the j-th trajectory of the fanplot.
%                 Example - 'highlight',[2 4 20]
%                 Data Types - double
%
%         label :   Labels. Cell array of strings. Cell containing the
%                   labels of the units (optional argument used when
%                   datatooltip=1). If this field is not present labels
%                   row1, ..., rown will be automatically created and
%                   included in the pop up datatooltip window.
%                   Example - 'corr',1
%                   Data Types - Cell array of strings
%
%       labx    :   x-axis label. String. A label for the x-axis (default:
%                   'Subset size m').
%                   Example - 'labx','Subset size m'
%                   Data Types - char or string
%
%       laby    :   y-axis label. String. a label for the y-axis
%                   (default:'Score test statistic').
%                   Example - 'laby','Score test statistic'
%                   Data Types - char or string
%
%       lwd     :   Linewidth. Scalar. Scalar which controls linewidth of
%                   the curves which contain the score test.
%                   Default line width=2.
%                   Example - 'lwd',2
%                   Data Types - double
%
%       lwdenv  :   Width of the envelope lines. Scalar. Scalar which
%                   controls the width of the lines associated with the
%                   envelopes. Default is lwdenv=1.
%                   Example - 'lwdenv',1
%                   Data Types - double
%
%       nameX   :   Labels of the X variables. Cell array of strings. Cell
%                   array of strings of length p containing the labels of
%                   the varibles of the regression dataset. If it is empty
%                 	(default) the sequence X1, ..., Xp will be created
%                   automatically.
%                   Example - 'nameX',''
%                   Data Types - Cell array of strings
%
%       namey   :   Labels of the y variable. String. String containing the
%                   label of the response variable.
%                   Example - 'namey',''
%                   Data Types - char
%
%multiPanel :  plots on a single or multiple panels. Boolean. If
%                   multiPanel is false each trajectory appears on a
%                   seperate subplot. If multiPanel is true (defaul) all the
%                   trajectories appear in single plot.
%                   Example - 'multiPanel',true
%                   Data Types - logical
%
%    SizeAxesNum:   Size of the numbers of the axis. Scalar. Scalar which
%                   controls the size of the numbers of the axes.
%                   Default value is 10.
%                   Example - 'SizeAxesNum',10
%                   Data Types - double
%
%       tag     :   Handle of the plot. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag pl_fan. Notice that if the
%                   program finds a plot which has a tag equal to the one
%                   specified by the user, then the output of the new plot
%                   overwrites the existing one in the same window else a
%                   new window is created.
%                   Example - 'tag','pl_mycov'
%                   Data Types - char
%
%       titl    :   Title. String. A label for the title (default: 'Fan plot')
%                   Example - 'titl','Fan plot'
%                   Data Types - char or string
%
%       xlimx   :   Min and Max of the x axis. Vector. Vector with two
%                   elements controlling minimum and maximum of the x axis.
%                   Default value is [init n].
%                   Example - 'xlimx',[init n]
%                   Data Types - double
%
%       ylimy   :   Min and Max of the y axis. Vector. Vector with two
%                   elements controlling minimum and maximum of the y axis.
%                   Default value for
%                   ylimy(1)=max(min(score_test),-20). Default value for
%                   ylimy(2)=min(max(score_test),20).
%                   Example - 'ylimy',[0 100]
%                   Data Types - double
%
%   datatooltip :   Information about the unit selected. Empty value,
%                   scalar or structure. The default is datatooltip=''. If
%                   datatooltip is not empty the user can use the mouse in
%                   order to have information about the unit selected, the
%                   step in which the unit enters the search and the
%                   associated label. If datatooltip is a structure, it is
%                   possible to control the aspect of the data cursor (see
%                   function datacursormode for more details or the
%                   examples below). The default options of the structure
%                   are DisplayStyle='Window' and SnapToDataVertex='on'.
%                   Example - 'datatooltip',1
%                   Data Types - Empty value, scalar or structure
%
%    databrush :    Databrush options. Empty value, scalar or cell.
%                   DATABRUSH IS AN EMPTY VALUE: If databrush is an empty
%                   value (default), no brushing is done. The activation of
%                   this option (databrush is a scalar or a cell) enables
%                   the user  to select a set of trajectories in the
%                   current plot and to see them highlighted in the y|X
%                   plot (notice that if the plot y|X does not exist it is
%                   automatically created). In addition, brushed units can
%                   be highlighted in the other following plots
%                   (only if they are already open): monitoring residual plot
%                   monitoring leverage plot maximum studentized residual
%                   $s^2$ and $R^2$ Cook distance and modified Cook
%                   distance deletion t statistics. The window style of the
%                   other figures is set equal to that which contains the
%                   monitoring residual plot. In other words, if the
%                   monitoring residual plot is docked all the other
%                   figures will be docked too. DATABRUSH IS A SCALAR: If
%                   databrush is a scalar the default selection tool is a
%                   rectangular brush and it is possible to brush only once
%                   (that is persist=''). DATABRUSH IS A CELL: If databrush
%                   is a cell, it is possible to use all optional arguments
%                   of function selectdataFS.m and LXS.m inside the curly
%                   brackets of option databrush and the following optional
%                   argument:
%                  persist = Persist is an empty value or a scalar containing
%                   the strings 'on' or 'off'. If persist = 'on' or 'off'
%                   brusing can be done as many time as the
%                   user requires. In case persist='off', every time a new
%                   brush is performed, units previously brushed are
%                   removed. In case persist='on' the unit(s)
%                   currently brushed are added to those previously
%                   brushed. However in both cases, if the user brushes a
%                   different trajectory from the one previously brushed,
%                   the previos brushed plots are stored in a figure in the
%                   background. The default value of persist is '' that is
%                   brushing is allowed only once.
%                 bivarfit = This option adds one or more least square
%                   lines, based on SIMPLE REGRESSION of y on Xi, to the
%                   plots of y|Xi.
%                   If bivarfit = ''
%                   is the default: no line is fitted.
%                   If bivarfit = '1'
%                   fits a single ols line to all points of each bivariate
%                   plot in the scatter matrix y|X.
%                   If bivarfit = '2'
%                   fits two ols lines: one to all points and another to
%                   the group of the genuine observations. The group of the
%                   potential outliers is not fitted.
%                   If bivarfit = '0'
%                   fits one ols line to each group. This is useful for the
%                   purpose of fitting mixtures of regression lines.
%                   If bivarfit = 'i1' or 'i2' or 'i3' etc
%                   fits an ols line to a specific group, the one with
%                   index 'i' equal to 1, 2, 3 etc. Again, useful in case
%                   of mixtures.
%                 multivarfit = This option adds one or more least square lines,
%                   based on MULTIVARIATE REGRESSION of y on X, to the
%                   plots of y|Xi.
%                   If multivarfit = ''
%                   is the default: no line is fitted.
%                   If multivarfit = '1'
%                   fits a single ols line to all points of each bivariate
%                   plot in the scatter matrix y|X. The line added to the
%                   scatter plot y|Xi is avconst +Ci*Xi, where Ci is the
%                   coefficient of Xi in the multivariate regression and
%                   avconst is the effect of all the other explanatory
%                   variables different from Xi evaluated at their centroid
%                   (that is overline{y}'C))
%                   If multivarfit = '2'
%                   exactly equal to multivarfit ='1' but this time we add the
%                   line based on the group of unselected observations.
%                   Example - 'databrush',1
%                   Data Types - Empty value, scalar or cell.
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
%
% See also: FSRfan, fanBIC, FSRaddt
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
%
% Atkinson, A.C. and Riani, M. (2002a), Tests in the fan plot for robust,
% diagnostic transformations in regression, "Chemometrics and Intelligent
% Laboratory Systems", Vol. 60, pp. 87-100.
%
% Atkinson, A.C. Riani, M., Corbellini A. (2019), The analysis of
% transformations for profit-and-loss data, Journal of the Royal
% Statistical Society, Series C, "Applied Statistics",
% https://doi.org/10.1111/rssc.12389
%
% Atkinson, A.C. Riani, M. and Corbellini A. (2021), The Box–Cox
% Transformation: Review and Extensions, "Statistical Science", Vol. 36,
% pp. 239-255, https://doi.org/10.1214/20-STS778
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('fanplotFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% fanplotFS with all default options.
    % load the wool data
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % FSRfan and fanplotFS with all default options
    [out]=FSRfan(y,X);
    fanplotFS(out);
%}
%
%{
    %%fanplotFS with optional arguments.
    %FSRfan and fanplotFS with specified lambda
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    % la = vector contanining the most common values of the transformation
    % parameter
    la=[-1 -0.5 0 0.5 1];
    [out]=FSRfan(y,X,'la',la);
    fanplotFS(out);
%}
%
%{
    %Interactive_example
    %FSRfan and fanplotFS with databrush option.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    la=[-1 -0.5 0 0.5 1];
    [out]=FSRfan(y,X,'la',la);
    fanplotFS(out,'databrush','1');
%}
%
%{
    %Interactive_example
    %FSRfan and fanplotFS with databrush, persist, label and RemoveLabels options.
    %Removelabels is a parameter of SelectdataFS function
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    la=[-1 -0.5 0 0.5 1];
    [out]=FSRfan(y,X,'la',la);
    fanplotFS(out,'databrush',{ 'persist' 'on' 'Label' 'on' 'RemoveLabels' 'off'});
%}
%
%{
    %Interactive_example
    %FSRfan and fanplotFS with databrush, bivarfit, label and  RemoveLabels options.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    la=[-1 -0.5 0 0.5 1];
    [out]=FSRfan(y,X,'la',la);
    fanplotFS(out,'databrush',{ 'bivarfit' '2' 'Label' 'on' 'RemoveLabels' 'off'});
%}

%{
    %Interactive_example
    %FSRfan and fanplotFS with databrush  and selectionmode options.
    %Example of the use of persistent cumulative brush.
    %Every time a brushing action is performed
    %current highlightments are added to previous highlightments

    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    la=[-1 -0.5 0 0.5 1];
    [out]=FSRfan(y,X,'la',la);
    fanplotFS(out,'databrush',{'selectionmode','Brush'});
    fanplotFS(out,'databrush',{'selectionmode' 'Lasso' 'persist' 'off'})
    fanplotFS(out,'databrush',{'selectionmode' 'Rect' 'persist' 'on'})
%}

%{
    %fanplotFS with datatooltip passed as scalar.
    %That is using default options for datacursor (i.e. DisplayStyle=window).
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    la=[-1 -0.5 0 0.5 1];
    [out]=FSRfan(y,X,'la',la);
    fanplotFS(out,'datatooltip',1);
%}

%{
   % Construct fan plot specifying the confidence level and the xlimits.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    la=[0 1/3 0.4 0.5];
    [out]=FSRfan(y,X,'la',la,'init',size(X,2)+2,'nsamp',20000);
    fanplotFS(out,'xlimx',[100 300],'conflev',0.95);
%}

%{
    %Interactive_example
    %Example of the use of multivarfit and xlimx.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    la=[0 1/3 0.4 0.5];
    [outs]=FSRfan(y,X,'la',la,'init',size(X,2)+2,'nsamp',20000);
    fanplotFS(outs,'xlimx',[10 520],'databrush',{'selectionmode' 'Brush' 'multivarfit' '2'})
%}

%{
    %Interactive_example
    %Example of the use of FlagSize, namey, namex, lwd,FontSize, SizeAxesNum.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    la=[-1 -0.5 0 0.5 1];
    [out]=FSRfan(y,X,'la',la);
     namey='Sales';
     nameX={'Number of visits', 'Age', 'Number of persons in the family'};
    %FlagSize controls how large must be the highlighted points. It is a
    %parametr of selectdataFS.
     fanplotFS(out,'xlimx',[10 520],'lwd',1.5,'FontSize',11,'SizeAxesNum',11,'nameX',nameX,'namey',namey,'databrush',{'selectionmode' 'Brush'...
     'multivarfit' '2' 'FlagSize' '5'})
%}

%{
    %Interactive_example
    % Only one brush specifying labels for y and X.
    load loyalty.mat
    y=loyalty{:,4};
    X=loyalty{:,1:3};
    la=[-1 -0.5 0 0.5 1];
    [out]=FSRfan(y,X,'la',la);
    namey=loyalty.Properties.VariableNames{4};
    nameX=loyalty.Properties.VariableNames(1:3);
    fanplotFS(out,'databrush',{'selectionmode' 'Brush' 'FlagSize' '5'},'nameX',nameX,'namey',namey)
%}

%{
    % fanplotFS with option highlight (first example).
    % load the wool data
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % FSRfan and fanplotFS with all default options
    [out]=FSRfan(y,X);
    % Units 24 27 and 23 enter the last three steps in the search with la=0.
    % Option highlight enables us to understand when these 3 units join the
    % subset for the other values of lambda
    fanplotFS(out,'highlight',[24 27 23]);
%}

%{
    %% fanplotFS with option highlight (second example).
    load gasoline.mat
    y=gasoline{:,2};
    X=gasoline{:,1};
    [out]=FSRfan(y,X,'la', [-1 0 1]);
    % In the code below, option highlight enables us to understand when units
    % [21 18 33 77] join the subset for the difference values of lambda
    fanplotFS(out,'highlight',[ 21 18 33 77]);
%}

%{
    % Example of option highlight when highlight is a 2D array.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    % la = vector contanining the most common values of the transformation
    % parameter lambda
    la=[0 0.5 1];
    [outFSRfan]=FSRfan(y,X,'la',la);
    
    % Store the units declares as outliers for the 3 values of lambda
    Highl=NaN(250,5);
    
    for j=1:length(la)
        if abs(la(j))>1e-06
            outlaj=FSR(y.^la(j),X,'msg',0,'plots',0);
        else
            outlaj=FSR(log(y),X,'msg',0,'plots',0);
        end
        oultlj=outlaj.outliers;
        Highl(1:length(oultlj),j)=oultlj;
    end
    % Highlight in red inside the fanplot the steps of entry of the units
    % declared as outliers for the different trajectories of lambda
    fanplotFS(outFSRfan,'highlight',Highl,'ylimy',[-50 20],'xlimx',[30 510]);
%}

%{
    %% fanplotFS based on the output of FSRaddt.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    nameX={'F1','F2','F3'};
    [out]=FSRaddt(y,X,'plots',0);
    % out.la contains the names of the variables which have to be shown
    out.la=nameX;
    fanplotFS(out);
%}

%{
    %% fanplotFS based on the output of Sregeda.
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);
    % Specify the type of estimator of t statistic
    covrob=2; 
    % Specify the rho function
    rhofunc='optimal'; 
    % Specify line width of the envelopes
    lwdenv=2;
    [outS]=Sregeda(y,X,'covrob',covrob,'rhofunc',rhofunc);
    fanplotFS(outS,'conflev',0.90,'lwdenv',lwdenv);
%}

%{
    %% fanplotFS based on the output of MMregeda.
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);
    % Specify the type of estimator of t statistic
    covrob=5; 
    % Specify the rho function
    rhofunc='optimal'; 
    % Specify line width of the envelopes
    lwdenv=2;
    [outMM]=MMregeda(y,X,'covrob',covrob,'rhofunc',rhofunc);
    fanplotFS(outMM,'conflev',0.90,'lwdenv',lwdenv);
%}

%{
    %% fanplotFS based on the output of FSReda.
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);
    % Specify line width of the envelopes
    lwdenv=2;
    [out]=LXS(y,X,'nsamp',1000);
    [outFS]=FSReda(y,X,out.bs);
    fanplotFS(outFS,'conflev',0.90,'lwdenv',lwdenv);
%}

%% Beginning of code
brushedUnits=[];


%% Input parameters checking

y=out.y;
X=out.X;
[n,p]=size(X);
if range(X(:,1))==0
    intercept=1;
else
    intercept=0;
end

if strcmp(out.class,'FSRfan')
    fanplotScore=true;
    Sco=out.Score;
    laby='Score test statistic';
    labx= 'Subset size m';

    la=out.la(:);
    las=string(la);
    dfvary=false;
    estimatorS=false;
    finalvalue=n;
    intercept=0;
    maxy=20;
    statName="X";
    xrange=[Sco(1,1) Sco(end,1)+1];


elseif strcmp(out.class,'FSRaddt')
    fanplotScore=false;
    Sco=out.Tdel;
    laby='Deletion t statistics';
    labx= 'Subset size m';
    dfvary=true;
    estimatorS=false;
    finalvalue=n;
    intercept=0;
    maxy=50;
    statName="X";
    xrange=[Sco(1,1) Sco(end,1)+1];


elseif strcmp(out.class,'Sregeda')
    fanplotScore=false;
    Sco=out.tStat;
    laby='t statistics from S estimator';
    labx= 'Break down point';

    dfvary=false;
    estimatorS=true;
    finalvalue=out.bdp(end)-0.005;
    maxy=50;
    statName="t_";
    xrange=[out.bdp(end) out.bdp(1)];


elseif strcmp(out.class,'MMregeda')
    fanplotScore=false;
    Sco=out.tStat;
    laby='t statistics from MM estimator';
    labx= 'Efficiency';
    dfvary=false;
    estimatorS=false;
    finalvalue=out.eff(end)+0.01;
    maxy=50;
    statName="t_";
    xrange=[out.eff(1) out.eff(end)+0.01];

elseif  strcmp(out.class,'FSReda')
    fanplotScore=false;
    Sco=out.Tols;
    laby='t statistics from FS';
    labx= 'Subset size m';
    dfvary=true;
    estimatorS=false;
    finalvalue=n;
    maxy=50;
    statName="t_";
    xrange=[Sco(1,1) Sco(end,1)+1];

else
    error('Unknown class')
end

if fanplotScore==false
    if isfield(out,'la')
        % if out.la is a cell array of characters convert it to string array.
        if iscell(out.la)
            out.la=string(out.la);
        end
        % Extract the variables for which deletion tstat have been computed
        % Note that there is -1 because intercept is present in out.X
        if isstring(out.la)
            las=out.la(:);
        else
            las=statName+"{"+string(out.la(:)-1)+"}";
        end
        la=las;
    else
        las=statName+"{"+(1:(p-intercept))+"}";
        la=las;
    end
end

%% User options
options=struct('conflev',0.99,'titl','Fan plot','labx',labx,...
    'laby',laby,'xlimx','','ylimy','','lwd',2,'lwdenv',1, ...
    'FontSize',12,'SizeAxesNum',10,'highlight',[],...
    'tag','pl_fan','datatooltip','','databrush','', ...
    'nameX','','namey','','label','','multiPanel', ...
    false,'addxline',[],'flabstep',[]);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:fanplotFS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    aux.chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end


% find and store the value of option labeladd inside cell options.databrush
d=find(strcmp('labeladd',options.databrush));
if d>0
    labeladd=options.databrush(d+1);
    labeladd=labeladd{1};
    % This option must be removed because is not an option of function
    % selectdataFS
    options.databrush(d:d+1)=[];
else
    labeladd='';
end

% find and store the value of option bivarfit inside cell options.databrush
d=find(strcmp('bivarfit',options.databrush));
if d>0
    bivarfit=options.databrush(d+1);
    bivarfit=bivarfit{1};
    % This option must be removed because is not an option of function
    % selectdataFS
    options.databrush(d:d+1)=[];
else
    bivarfit='';
end

% find and store the value of option multivarfit inside cell options.databrush
d=find(strcmp('multivarfit',options.databrush));
if d>0
    multivarfit=options.databrush(d+1);
    multivarfit=multivarfit{1};
    % This option must be removed because is not an option of function
    % selectdataFS
    options.databrush(d:d+1)=[];
else
    multivarfit='';
end

% find and store the value of option persist inside cell options.databrush
d=find(strcmp('persist',options.databrush));
if d>0
    persist=options.databrush(d+1);
    % This option must be removed because is not an option of function
    % selectdataFS
    options.databrush(d:d+1)=[];

    % Start highlighting with red color
    ColorOrd=[1 0 0;0 1 1; 1 0 1; 1 1 0; 0 0 0; 0 1 0; 0 0 1];
    ColorOrd=repmat(ColorOrd,4,1);
else
    persist='';
end

% Initialize colors for plot y|X
d=find(strcmp('FlagColor',options.databrush));
if d>0
    flagcol=options.databrush{d+1};
    clr=['b' flagcol 'cmykg'];
else
    % default color for brushing is red
    flagcol='r';
    % default colors for y|X are blue (unbrushed unit) and red (brushed
    % units)
    clr='brcmykgbrcmykg';
end

% Symbol types for y|X plot (in case of brushing)
styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};

%% Display the fanplot

% Specify where to send the output of the current procedure
h=findobj('-depth',1,'tag',options.tag);
if (~isempty(h))
    clf(h);
    figure(h)
    axes;
else
    figure;
end

if fanplotScore == true
    set(gcf,'Name',['fanplot for lambda=' mat2str(out.la) ]);
end

lla=length(la);

% lwd = line width of the trajectories which contain the score test
lwd=options.lwd;

% set the x and y axis
xlimx=options.xlimx;
ylimy=options.ylimy;

% lwdenv = line width of the curves associated with the envelopes
lwdenv=options.lwdenv;
conflev=options.conflev;
% FontSize = font size of the axes labels
FontSize =options.FontSize;
addxline = options.addxline;

% Use default limits for y axis
ylim1=max(-maxy,min(min(Sco(:,2+intercept:end))));
ylim1=min([ylim1 -3]);
ylim2=min(maxy,max(max(Sco(:,2+intercept:end))));
ylim2=max([ylim2 3]);

% Main title of the plot and labels for the axes
labx=options.labx;
laby=options.laby;
titl=options.titl;
highlight=options.highlight;
multiPanel=options.multiPanel;
flabstep=options.flabstep;

if ~isempty(highlight)

    if isvector(highlight)
        nBrush=repmat(highlight(:),1,lla);
    else
        nBrush=highlight;
    end
end

% SizeAxesNum = font size for the axes numbers
SizeAxesNum=options.SizeAxesNum;

% Line Width of vertical line
lwdVerticalLine=1.5;

if multiPanel==false

    plot1=plot(Sco(:,1),Sco(:,2+intercept:end),'LineWidth',lwd);
    set(gcf,'Tag',options.tag)

    % Specify the line type for the units inside vector units
    slin={'-';'--';':';'-.'};
    slin=repmat(slin,ceil(lla/4),1);
    set(plot1,{'LineStyle'},slin(1:lla));


    if ~isempty(xlimx)
        xlim(xlimx);
    else
        xlim(xrange)
    end

    if isempty(ylimy)
        ylim([ylim1 ylim2]);
    else
        % Use limits specified by the user
        ylim(ylimy);
    end

    % PLOT THE CONFIDENCE BANDS OF THE SCORE TESTS
    rangeaxis=axis;

    if dfvary==true
        conflev=(1+conflev)/2;
        % df = m-p
        Tdelenv=tinv(repmat(conflev,length(Sco(:,1)),1),repmat(Sco(:,1),1,length(conflev))-p);

        for i=1:length(conflev)
            % Add adaptive confidence bands
            line(Sco(:,[1 1]),[-Tdelenv(:,i) Tdelenv(:,i)],'LineWidth',lwdenv,'color','r','Tag','env');
        end
    else
        quant = sqrt(chi2inv(conflev,1));
        numconflev=length(conflev);
        V=repmat([rangeaxis(1);rangeaxis(2)],1,2*numconflev);
        QUANT=[[quant;quant],[ -quant;-quant]];
        % Assign to the confidence lines Tag env so that they cannot be selected
        % with options databrush
        line(V, QUANT,'LineWidth',lwdenv,'color','r','Tag','env');
    end

    % Add labels at the end of the search


    if isempty(flabstep)
        % Add labels at the end of the search
        text(finalvalue*ones(lla,1),Sco(end,2+intercept:end)',las,'FontSize',FontSize);
    else
        % Add labels in the steps specified by flabstep
        for ii=1:length(flabstep)
            posii=find(Sco(:,1)>=flabstep(ii),1);
            if ~isempty(posii)
                text(Sco(posii,1)*ones(lla,1),Sco(posii,2+intercept:end)'+0.1,las,'FontSize',FontSize);
            end
        end
    end

    if estimatorS==true
        set(gca,'XDir','reverse');
    end

    title(titl);

    % Add to the plot the labels for values of la
    % Add the horizontal lines representing asymptotic confidence bands
    xlabel(labx,'Fontsize',FontSize);
    ylabel(laby,'Fontsize',FontSize);

    % Add the vertical line(s)
    if ~isempty(addxline)
        xline(addxline,'LineWidth',lwdVerticalLine,'Color','k')
    end
    set(gca,'FontSize',SizeAxesNum)

    box on

    if ~isempty(highlight)

        hold('on');

        for j=1:lla
            nbrush=rmmissing(nBrush(:,j));

            if iscell(out.Un)
                Un=out.Un{j};
            else
                Un=out.Un;
            end
            xdata=Sco(:,1)'; % x coordinates of score (steps)
            ydata=Sco(:,j+intercept+1)'; % y coordinates of score (values)

            % Find step of entry of selected units and the associated
            % x,y coordinates
            [selsteps, xcoord,ycoord]= computeHG(n,nbrush,Un,xdata,ydata);

            if fanplotScore==true
                disp("Steps of entry of selected units when la= "+las(j));
            else
                disp(['Steps of entry of selected units in ' laby '  for ' char(las(j))]);
            end
            disp(selsteps);
            plot(gca,xcoord,ycoord,'LineWidth',3.5,'color',[1 0 0],'tag','brush_mdr');
        end

        hold('off');
    end
else % else plot on separate panels
    ij=0;

    if dfvary==true
        conflev=(1+conflev)/2;
        % df = m-p
        Tdelenv=tinv(repmat(conflev,length(Sco(:,1)),1),repmat(Sco(:,1),1,length(conflev))-p);
    else
        quant = sqrt(chi2inv(conflev,1));
        numconflev=length(conflev);

    end

    ht = tiledlayout('flow');
    for j=2+intercept:size(Sco,2)
        nexttile
        ij=ij+1;
        plot(Sco(:,1),Sco(:,j),'LineWidth',lwd);
        set(gcf,'Tag',options.tag)


        if ~isempty(xlimx)
            xlim(xlimx);
        else
            xlim(xrange)
        end

        if isempty(ylimy)
            ylim([ylim1 ylim2]);
        else
            % Use limits specified by the user
            ylim(ylimy);
        end

        % PLOT THE CONFIDENCE BANDS OF THE SCORE TESTS

        if dfvary==true
            for i=1:length(conflev)
                % Add adaptive confidence bands
                line(Sco(:,[1 1]),[-Tdelenv(:,i) Tdelenv(:,i)],'LineWidth',lwdenv,'color','r','Tag','env');
            end
        else
            rangeaxis=axis;
            V=repmat([rangeaxis(1);rangeaxis(2)],1,2*numconflev);
            QUANT=[[quant;quant],[ -quant;-quant]];      % Assign to the confidence lines Tag env so that they cannot be selected
            % with options databrush
            line(V, QUANT,'LineWidth',lwdenv,'color','r','Tag','env');
        end

        if isempty(flabstep)
            % Add labels at the end of the search
            text(finalvalue,Sco(end,j)',las(ij),'FontSize',FontSize);
        else
            for ii=1:length(flabstep)
                posii=find(Sco(:,1)>=flabstep(ii),1);
                if ~isempty(posii)
                    text(Sco(posii,1),Sco(posii,j)',las(ij),'FontSize',FontSize);
                end
            end
        end

        % Add the vertical line(s)
        if ~isempty(addxline)
            xline(addxline,'LineWidth',lwdVerticalLine,'Color','k')
        end

        if estimatorS==true
            set(gca,'XDir','reverse');
        end

        set(gca,'FontSize',SizeAxesNum)
        box on
        hold('on')
        if ~isempty(highlight)
            nbrush=rmmissing(nBrush(:,ij));

            if iscell(out.Un)
                Un=out.Un{ij};
            else
                Un=out.Un;
            end
            xdata=Sco(:,1)'; % x coordinates of score (steps)
            ydata=Sco(:,ij+intercept+1)'; % y coordinates of score (values)

            % Find step of entry of selected units and the associated
            % x,y coordinates
            [selsteps, xcoord,ycoord]= computeHG(n,nbrush,Un,xdata,ydata);

            if fanplotScore==true
                disp("Steps of entry of selected units when la= "+las(ij));
            else
                disp(['Steps of entry of selected units in ' laby ' for ' char(las(ij))]);
            end
            disp(selsteps);
            plot(gca,xcoord,ycoord,'LineWidth',3.5,'color',[1 0 0],'tag','brush_mdr');
        end

    end
    title(ht,titl);
    % Add to the plot the x and y labels
    xlabel(ht,labx,'Fontsize',FontSize);
    ylabel(ht,laby,'Fontsize',FontSize);

end
%% Set the datatooltip for the fanplot

if ~isempty(options.datatooltip)
    try
        hdt = datacursormode;
        if ~isstruct(options.datatooltip)
            set(hdt,'DisplayStyle','window');
        else
            % options.databrush contains a structure where the user can set the
            % properties of the data cursor
            set(hdt,options.datatooltip);
        end
        % Declare a custom datatip update function to display additional
        % information about the units
        set(hdt,'UpdateFcn',{@fanplotLbl,out})
    catch
        disp('No graphical device, interactive datatooltip not enabled')
    end
end

% Store all graphics elements of the current figure inside handle hmin
hfan=gcf;

%% Brush mode (call to function selectdataFS)

if (~isempty(options.databrush) || iscell(options.databrush))
    close(findobj('type','figure','Tag','pl_yX'));
    close(findobj('type','figure','Tag','pl_resfwd'));
    close(findobj('type','figure','Tag','pl_mdr'));

    if isscalar(options.databrush)
        sele={'selectionmode' 'Rect' 'Ignore' findobj(gcf,'tag','env') };
    else
        sele=[options.databrush {'Ignore' findobj(gcf,'tag','env')}];
    end

    % Check how the user set the option Label inside cell options.databrush
    d=find(strcmp('Label',options.databrush));
    if d>0
        Label=options.databrush{d+1};
    else
        Label='nolabel';
    end

    % This option is necessary to make sure that the user selects data from
    % the plot and not from other plots
    sele=[sele 'Tag' {options.tag}];

    % Check if X includes the constant term for the intercept. The variable
    % 'intercept' will be used later for plotting.
    intcolumn = find(max(X,[],1)-min(X,[],1) == 0);

    if intcolumn==1 && size(X,2)>1
        intercept=1;
        p1=1:(p-numel(intcolumn));
        Xsel=X;
        Xsel(:,intcolumn)=[];
    else
        intercept=0;
        p1=1:p;
        Xsel=X;
    end

    if isempty(options.nameX)
        nameX=cellstr(num2str(p1','X%d'));
    else
        nameX=options.nameX;
    end

    % group = VECTOR WHICH WILL CONTAIN THE IDENTIFIER OF EACH GROUP
    % e.g. group(14)=3 means that unit 14 was selected at the third brush
    group=ones(n,1);

    % some local variables
    but=0; brushcum=[]; ij=1;

    sele1=sele;

    %initialization of variable lasel. It will be a scalar associated with
    %the value of lambda which has been brushed.
    lasel='';

    %initialization of variable it. It will be a counter identifying the
    %iteration linked to the following while statement
    it=0;
    while but<=1
        it=it+1;
        figure(hfan);

        % Remark: function selectdataFS cannot be used on the current figure if
        % the "selection mode" or the "zoom tool" are on. Setting the
        % plotedit mode initially to on and then to off, has the effect to
        % deselect plotedit mode.

        plotedit on
        plotedit off

        if strcmp(persist,'off')
            % then remove from the current plot what has alredy been
            % selected
            % Remove the yellow selection in this plot if present
            a=findobj(gcf,'Tag','selected');
            delete(a);

        elseif strcmp(persist,'on')
            %then add to the option list sele, the list of colors for the
            %brushing and selections.
            sele1=[sele {'FlagColor' ColorOrd(ij,:)}];
        end
        if it>1
            disp('To continue brushing, with the mouse select additional steps to brush in the FAN plot');
            disp('To stop brushing press any key');
        else
            disp('To brush, with the mouse select steps in the FAN plot');
        end

        [pl,xs,ys] = aux.selectdataFS(sele1{:},'Label','off');

        % exit from levfwdplot if the levfwdplot figure was closed before selection
        if isnumeric(pl) && ~isempty(pl) && (max(pl) == -999)
            return
        end

        % selsteps= vector which contains the list of the steps which have
        % been brushed
        % sely = corresponding values of the score test associated
        % to the brushed steps
        if iscell(xs)
            selsteps=cell2mat(xs);
            sely=cell2mat(ys);
        else
            selsteps=xs;
            sely=ys;
        end

        lselsteps=length(selsteps);
        selindex=zeros(lselsteps,1);

        %  find to which value of \lambda the selected steps refer to
        if ~isempty(sely)

            % col1 is associated to the columns of \lambda which have been
            % brushed
            col1=zeros(length(sely),1);
            for i=1:length(sely)
                [~,col] = find(Sco(:,2:end)==sely(i));
                col1(i)=col;
            end

            %list_lambda is the list of lambda trajectories that the User
            %has selected
            list_lambda=unique(col1);
            %length_list_lambda is the number of lambda trajectories that
            %the User has selected
            length_list_lambda=length(list_lambda);

            %the for loop has been introduced because the ueser could
            %select more than one lambda trajectories
            for j =1: length_list_lambda

                %change_lambda=yes if the user selects a trajectory
                %different to the one selected at step j-1
                %change_lambda=no if the user selects the same trajectory
                %that he has selected in the step j-1

                if list_lambda(j)==lasel
                    change_lambda='no';
                else
                    change_lambda='yes';

                    %if the selected lambda has changed, all the existing
                    %plots will be stored in a new figure
                    h=findobj('-depth',1,'Tag','pl_yX');
                    if ~isempty(h)

                        %creation of a second figure which will contain a
                        %copy of the gplotmatrix, just created
                        hf2=figure;
                        figure(hf2);
                        set(hf2,'Name',['Results for lambda=' num2str(la(lasel)) ]);

                        %creation of subplot in figure hf2 where I want to
                        %copy the gplotmatrix of figure f1
                        s1=subplot(2,2,1);

                        %pos is the position of the subplot in hf2
                        pos=get(s1,'Position');

                        %I copy the gplotmatrix in figure hf2
                        delete(s1);
                        hax2=copyobj(AX,hf2);

                        sizeAX = size(AX,2);
                        if sizeAX==1
                            set(hax2, 'Position', pos);
                        else
                            space = .02;
                            for i=1:sizeAX
                                pos1=get(hax2(i),'Position');
                                pos1(1,1) = pos(1,1)+(i-1)*(pos(1,3)/sizeAX);
                                pos1(1,2) = pos(1,2);
                                pos1(1,3) = (1-space)*pos(1,3)/sizeAX;
                                pos1(1,4) = pos(1,4);
                                set(hax2(i), 'Position', pos1);
                            end
                        end
                    end


                    h=findobj('-depth',1,'Tag','pl_resfwd');
                    if ~isempty(h)
                        figure(h)
                        hax1=gca;
                        figure(hf2);
                        set(hf2,'Name',['Results for lambda=' num2str(la(lasel)) ]);
                        pos123 =strcat('2','3',num2str(length(AX)+1));
                        s1=subplot(pos123);
                        pos=get(s1,'Position');
                        delete(s1);
                        hax2=copyobj(hax1,hf2);
                        set(hax2, 'Position', pos);

                    end
                    h=findobj('-depth',1,'Tag','pl_mdr');
                    if ~isempty(h)
                        figure(h)
                        hax1=gca;
                        figure(hf2);
                        set(hf2,'Name',['Results for lambda=' num2str(la(lasel)) ]);
                        pos123 =strcat('2','3',num2str(length(AX)+2));
                        s1=subplot(pos123);
                        pos=get(s1,'Position');
                        delete(s1);
                        hax2=copyobj(hax1,hf2);
                        set(hax2, 'Position', pos);

                    end
                end

                % lasel = scalar associated with the j-th value of lambda
                % selected (brushed).
                lasel=list_lambda(j);

                Un=out.Un{lasel};

                for i=1:lselsteps
                    selindex(i)=find(Un(:,1)==selsteps(i));
                end


                % Find the units which entered the search in the brushed steps
                Unsel=Un(selindex,2:end);

                % nbrush= vector which contains the list of the selected units
                nbrush=Unsel(~isnan(Unsel));

                h=findobj('-depth',1,'Tag','pl_yX');

                if (~isempty(h))
                    % delete from the current figure all graphics objects
                    % whose handles are not hidden
                    clf(h);
                    % Make the figure identified by handle h become the
                    % current figure make it visible, and raise it above
                    % all other figures on the screen.
                    figure(h);

                else
                    % create a new figure
                    hh=figure;
                    % The window style of the new figure which has been
                    % created is set equal to that which contained the
                    % monitoring residual plot
                    set(hh,'WindowStyle',get(hfan,'WindowStyle'));
                end

                if strcmp(persist,'on') && strcmp(change_lambda,'no')
                    %if persist is on and the selected trajectory is the
                    %same selected at step j-1 then brushcum is the list of
                    %units selected from steps 1 to j.
                    brushcum=unique([brushcum; nbrush]);
                else
                    %if:
                    % 1. persist is on but the selected trajectory is
                    %    different from the one selected at step j-1.
                    % 2. persist is off and the selected trajectory is the
                    %    same selected at step j-1
                    % 3. persist is off and the selected trajectory is
                    %    different from the one selected at step j-1
                    %then brushcum is the list of units selected at step j.
                    brushcum=nbrush;
                    group=ones(n,1);
                end

                %group is a vector of length equal to the number of input
                %onservations. The unselected observations are represented
                %with ones. The selected observations with ij+1.
                group(nbrush)=ij+1;

                %unigroup is the list of groups.
                unigroup=unique(group);

                if fanplotScore==true
                    % la(lasel) is the value of lambda selected
                    if la(lasel) ==0
                        % Create label for transformed response
                        if isempty(options.namey)
                            namey=char('log(y)');
                        else
                            namey=char(['log(' options.namey ')']);
                        end
                        % Create transformed response
                        y1=log(y);

                    else
                        % Create label for transformed response
                        if isempty(options.namey)
                            namey=char(['y^{' num2str(la(lasel)) '}']);
                        else
                            namey=char([options.namey '^{' num2str(la(lasel)) '}']);
                        end

                        % Create transformed response
                        y1=y.^la(lasel);
                    end
                else
                    y1=y;
                    if isempty(options.namey)
                        namey='y';
                    else
                        namey=char(options.namey);
                    end
                end
                %% Run the forward search using the selected value of \lambda

                % i.e. use transformed response which is contained in
                % vector y1
                [oute] = FSReda(y1,X,out.bs(:,lasel),'nocheck',1);

                %% Display the yXplot with the corresponding groups of units highlighted
                [H,AX,BigAx] = gplotmatrix(Xsel,y1,group,clr(unigroup),char(styp{unigroup}),[],'on',[],nameX,namey);

                if fanplotScore==true
                    set(gcf,'Name',['Scatter plot matrix y^' num2str(la(lasel))  '|X with different symbols for brushed units']);
                else
                    set(gcf,'Name','Scatter plot matrix y|X with different symbols for brushed units');
                end
                set(gcf,'tag','pl_yX');
                %set(gcf,'Name',['yXplot: results for lambda=' num2str(la(lasel)) ]); %DDD

                if nbrush>0
                    set(H(:,:,1),'DisplayName','Unbrushed units');
                    set(H(:,:,2),'DisplayName','Brushed units');
                else
                    set(H,'DisplayName',' Units');
                end

                set(AX,'FontSize',FontSize)
                % Extract Xlabel and set the font

                xlab=get(AX,'Xlabel');
                if iscell(xlab)
                    for jj=1:length(xlab)
                        set(xlab{jj},'FontSize',SizeAxesNum)
                    end
                else
                    set(xlab,'FontSize',SizeAxesNum)
                end

                % Extract ylabel
                ylab=get(AX,'Ylabel');
                if iscell(ylab)
                    for jj=1:length(xlab)
                        set(ylab{jj},'FontSize',SizeAxesNum)
                    end
                else
                    set(ylab,'FontSize',SizeAxesNum)
                end

                % save the indices of the last selected units (nbrush) to the
                % 'UserData' field of the last selected group of H(:,:,end)
                set(H(:,:,end), 'UserData' , nbrush);

                % The following line adds objects to the panels of the yX
                % add2yX(H,AX,BigAx,oute,group,nbrush,bivarfit,multivarfit,labeladd)
                add2yX(H,AX,BigAx,'intercept',intercept,'bivarfit',bivarfit,'multivarfit',multivarfit,'labeladd',labeladd);


                for mfc=1:length(unigroup)
                    set(findobj(gcf,'marker',char(styp(unigroup(mfc)))),'MarkerFaceColor',clr(unigroup(mfc)));
                end

                set(gcf,'Name',['Scatter plot matrix y^' num2str(la(lasel))  '|X with different symbols for brushed units']);


                % Assign to this figure tag pl_yX
                %set(gcf,'tag','pl_yX');


                %% Highlight brushed trajectories also in the resfwdplot, if it is open

                if (ij==1 ||  strcmp(change_lambda,'yes'))
                    % the plot of scaled residuals has to be completely
                    % recreated when the User has selected a lambda trajectory
                    % different from the one selected at step j-1

                    fground=struct;
                    fground.funit=nbrush;
                    fground.Color={'b'};
                    fground.LineStyle={'-'};
                    standard=struct;
                    standard.SizeAxesLab=FontSize;
                    standard.SizeAxesNum=SizeAxesNum;
                    resfwdplot(oute,'standard',standard,'fground',fground)

                    %                         % Check if options label is on or off
                    %                     if strcmp(Label,'on')==1
                    %                         resfwdplot(oute,'selunit',nbrush,'selunitcolor',{'b'},'selunitbold',1.5,'FontSize',FontSize,'SizeAxesNum',SizeAxesNum, 'selunittype',{'-'});
                    %                     else
                    %                         resfwdplot(oute,'selstep',0,'selunit',nbrush,'selunitcolor',{'b'},'selunitbold',1.5,'FontSize',FontSize,'SizeAxesNum',SizeAxesNum, 'selunittype',{'-'});
                    %                     end
                    %
                end

                % Now check if the figure which monitors the residuals is open.
                h=findobj('-depth',1,'Tag','pl_resfwd');
                set(h,'Name',strcat('Monitoring of Scaled Residuals for lambda=',num2str(la(lasel))));
                if (~isempty(h))

                    % make figure which contains monitoring scaled residuals
                    % become the current figure
                    figure(h);

                    % if but=0 then it is necessary to remove previous
                    % highlightments (even if persist='on')
                    if (strcmp(persist,'off') || but==0 )

                        % If set of values has already been highlighted in the
                        %monitoring residuals plot, remove it
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
                        ycoord=abrush; % y coordinates of scaled residual (value)
                    end

                    xcoord=get(a(1),'Xdata'); % x coordinates of mdr (steps)

                    hold('on');
                    % Highlight the curve in the monitoring residuals plot
                    if (strcmp('on',persist))
                        %if persist is on use ColorOrd(ij,:)
                        plot(gca,xcoord,ycoord,'LineWidth',2,'color',ColorOrd(ij,:),'tag','brush_res');
                    else
                        %if persist is off use flagcol
                        plot(gca,xcoord,ycoord,'LineWidth',2,'color',flagcol,'tag','brush_res');
                    end
                    if strcmp(Label,'on')==1
                        % add the labels of brushed units to the fan plot
                        text(repmat(xcoord(1),length(nbrush),1),ycoord(:,1),cellstr(num2str(nbrush)));
                        text(repmat(xcoord(end),length(nbrush),1),ycoord(:,end),cellstr(num2str(nbrush)));
                    else
                    end
                    hold('off');
                end

                disp('Brushed steps');
                disp(selsteps);

                disp('Associated brushed units');
                disp(nbrush);

                disp('Steps of entry of brushed units');
                disp(Un(selindex,:));


                %% Highlight brushed trajectories also in the mdrplot, if it is open

                % Default limits for x axis
                h=findobj('-depth',1,'Tag','pl_mdr');

                %the mdrplot has to be completely recreated if
                % 1 persist is off
                % 2 persist is on but we are in the first iteration
                if strcmp('off',persist) || (strcmp('on',persist)&& it==1)
                    mdrplot(oute,'quant', [0.01;0.5;0.99],'sign',0,'mplus1',0,'envm',n,'tag','pl_mdr','databrush','');

                    % Now check if the figure containing minimum deletion
                    % residual is open. If it is, units are brushed in this
                    % plot too
                    h=findobj('-depth',1,'Tag','pl_mdr');
                    set(h,'Name',['Monitoring of Minimum deletion residual for lambda=' num2str(la(lasel)) ])
                else

                    % the mdrplot has also to be completly created if persist is
                    % on and we are not in the first iteration but the user has
                    % selected a new lambda trajectory
                    if strcmp(change_lambda,'yes') &&  it>1
                        mdrplot(oute,'quant', [0.01;0.5;0.99],'sign',0,'mplus1',0,'envm',n,'tag','pl_mdr','databrush','');

                        % Now check if the figure containing minimum deletion
                        % residual is open.
                        % If it is, units are brushed in this plot too
                        h=findobj('-depth',1,'Tag','pl_mdr');
                        set(h,'Name',['Monitoring of Minimum deletion residual where lambda=' num2str(la(lasel)) ])
                    end
                end
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

                    [~,~, ib]=intersect(selsteps, xdata);

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

                    % Add to the previous matrix a row of missing values
                    % This operation is necessary if the steps are not contiguous
                    xy=cat(1,xy,NaN*zeros(1,2*xxlim));

                    % Reshape the set of x and y coordinates in two column
                    % vectors
                    % Notice the NaN between the steps which are not consecutive
                    xcoord=reshape(xy(:,1:xxlim),3*xxlim,1);
                    ycoord=reshape(xy(:,xxlim+1:end),3*xxlim,1);
                    hold('on');
                    if strcmp('on',persist)

                        % If necessary it isalso possible to specify a line
                        %style for the brushed steps
                        % 'LineStyle',stypbrushed{ij},
                        plot(gca,xcoord,ycoord,'LineWidth',4,'color',ColorOrd(ij,:),'tag','brush_mdr');
                        set(gca,'Position',[0.1 0.1 0.85 0.85])
                    else
                        plot(gca,xcoord,ycoord,'LineWidth',4,'color',flagcol,'tag','brush_mdr');
                    end
                    hold('off');
                end
            end

            % If the option persistent is not equal off or on than get out
            % of the loop
            if strcmp('on',persist) || strcmp('off',persist)

                % variable ij is linked to the highlighting color
                if strcmp('on',persist)
                    ij=ij+1;

                    % Set but=1 so that previous highlightments in other
                    % figures are not deleted
                    but=1;
                end

                % Before waitforbuttonpress:
                % - the fanplotFS is highlighted again
                figure(hfan);

                % Lay down the plots before continuing
                position(hfan);
                % position;

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
                if ss==1
                    but=2;
                end
            else
                but=2;
            end % close condition associated with persist 'on' or persist 'off'
        end % close condition associated with sely
    end % close brushing loop
    brushedUnits=brushcum;
end % close options.databrush

    function output_txt = fanplotLbl(~,event_obj,out)
        %%Provides information about the selected trajectory of the score test
        %(once selected with the mouse)
        %
        % Required input arguments:
        %
        %   obj     =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the sense that
        %               these arguments are automatically passed to the function when it executes.
        %   out     =   a structure coming from FSReda, Sregeda, MMregeda
        %               or FSRaddt
        %               label =(optional argument) if it is present it must be
        %                       a cell array of strings containing the labels of
        %                       the rows of the regression dataset
        % Output:
        %
        %   output_txt=  Datatip text (string or string cell array) which informs
        %                1) about the value of the statistic which is monitored
        %                2) the step of the search (just if FS) or bdp (if
        %                Sregeda) or eff (if MMregeda)
        %
        % REMARK: this function is called by function fanplotFS
        %
        % Written by FSDA team
        %

        %%
        pos = get(event_obj,'Position');

        % x and y, plot coordinates of the mouse
        x = pos(1); y = pos(2);

        % If structure out does not contain labels for the rows then
        % labels row1....rown are added automatically
        if isempty(intersect('label',fieldnames(out)))

            % Find n=number of units
            % cUn=cell2mat(out.Un(1));
            % n=cUn(end,1);
            out.label=cellstr(num2str((1:n)','row%d'));
        end

        % Add information about the unit(s) and the selected steps they entered the
        % search


        % Add information about the corresponding row label of what has been
        % selected
        if strcmp(out.class,'FSRfan')
            % output_txt is what it is shown on the screen
            output_txt = {['Value of the score test=',num2str(y)]};

            % Add information abou the step of the search which is under investigation
            output_txt{end+1} = ['Step m=' num2str(x)];

            [~,col] = find(Sco(:,2:end)==y);
            Un=out.Un{col};
            idx = find(Un(:,1)==x,1);
            sel=Un(idx,2:end);


            output_txt{end+1} = ['H_0: \lambda =' num2str(out.la(col))];
            output_txt{end+1} = ['Unit(s) entered in step ' num2str(x) '='  num2str(sel(~isnan(sel)))];


        elseif strcmp(out.class,'FSRaddt')
            output_txt = {['Value of the deletion t test=',num2str(y)]};
            % Add information abou the step of the search which is under investigation
            output_txt{end+1} = ['Step m=' num2str(x)];

            [~,col] = find(Sco(:,2:end)==y);
            Un=out.Un{col};
            idx = find(Un(:,1)==x,1);
            sel=Un(idx,2:end);
            output_txt{end+1} = ['H_0: \beta_' num2str(out.la(col)-1) '=0'];
            output_txt{end+1} = ['Unit(s) entered in step ' num2str(x) '='  num2str(sel(~isnan(sel)))];

        elseif strcmp(out.class,'FSReda')
            output_txt = {['Value of the t test=',num2str(y)]};
            % Add information abou the step of the search which is under investigation
            output_txt{end+1} = ['Step m=' num2str(x)];

            [~,col] = find(Sco(:,2+intercept:end)==y);
            Un=out.Un;
            idx = find(Un(:,1)==x,1);
            sel=Un(idx,2:end);
            output_txt{end+1} = ['H_0: \beta_' num2str(col) '=0'];
            output_txt{end+1} = ['Unit(s) entered in step ' num2str(x) '='  num2str(sel(~isnan(sel)))];

        elseif strcmp(out.class,'Sregeda')

            output_txt = {['Value of the t test=',num2str(y)]};
            % Add information about bdp
            output_txt{end+1} = ['bdp=' num2str(x)];
            [~,col] = find(Sco(:,2+intercept:end)==y);
            output_txt{end+1} = ['H_0: \beta_' num2str(col) '=0'];
        elseif  strcmp(out.class,'MMregeda')

            output_txt = {['Value of the t test=',num2str(y)]};
            % Add information about eff
            output_txt{end+1} = ['eff=' num2str(x)];
            [~,col] = find(Sco(:,2+intercept:end)==y);
            output_txt{end+1} = ['H_0: \beta_' num2str(col) '=0'];
        else

        end


    end


end

function [selsteps, xcoord,ycoord]=computeHG(n,nbrush,Un,xdata,ydata)
% nbrush= vector which contains the brushed steps selstesp=
% vector which will contain the steps in which the brushed
% units enter the search Given nbrush,Un,xdata,ydata, find selstesp
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
selsteps=rmmissing(selsteps,2,'MinNumMissing',size(selsteps,1));

% remove the columns of selsteps which contain just NaN


%% - highlight units in the fanplot

% Remove unnecessary rows from vector selsteps -1 is
% necessary because we are considering minimum outside
% selsteps=selsteps(:,1)-1;


[c, ia, ib]=intersect(selsteps(:,1), xdata); %#ok<ASGLU>
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
end

%FScategory:VIS-Reg