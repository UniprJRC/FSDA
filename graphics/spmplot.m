function [H,AX,BigAx] = spmplot(Y,varargin)
%spmplot produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal
%
%<a href="matlab: docsearchFS('spmplot')">Link to the help function</a>
%
%  Required input arguments:
%
%     Y : data matrix (2D array) containing n observations on v variables
%         or a structure 'out' coming from function FSMeda. Matrix or
%         structure.
%
%  Optional input arguments if Y is a 2D array:
%
%     If Y is a 2D array, varargin can be either a sequence of name/value
%     pairs, detailed below, or one of the following explicit assignments:
%
%       spmplot(Y,group);
%
%       spmplot(Y,group,plo);
%
%       spmplot(Y,group,plo,dispopt);
%
%     where group, plo and dispopt have the meaning described below.
%
%     If varargin{1} is a n-elements vector, then it is interpreted
%     as a grouping variable vector 'group'. In this case, it can only be
%     followed by 'plo' and 'dispopt'. Otherwise, the program expects a
%     sequence of name/value pairs.
%
%  List of optional input arguments if input Y is a matrix:
%
%  group: vector with n elements. It is a grouping variable that determines
%         the marker and color assigned to each point. It can be a categorical
%         variable, vector, string matrix, or cell array of strings.
%         Remark: if 'group' is used to distinguish a set of outliers from
%         a set of good units, the id number for the outliers should be the
%         larger (see optional field 'labeladd' of option 'plo' for details).
%
%
%    plo: empty value, scalar of structure which controls the names which
%         are displayed in the margins of the scatter-plot matrix and the
%         labels of the legend.
%
%         If plo is the empty vector [], then nameY and labeladd are
%           both set to the empty string '' (default), and no label and
%           no name is added to the plot.
%
%         If plo = 1 the names Y1,..., Yv are added to the margins of the
%           the scatter plot matrix else nothing is added.
%
%         If plo is a structure it may contain the following fields:
%         - labeladd: if it is '1', the elements belonging to the max(group)
%                in the spm are labelled with their unit row index.
%                The default value is labeladd = '', i.e. no label is added.
%         - nameY: cell array of strings containing the labels of the
%                variables. As default value, the labels which are added
%                are Y1, ..., Yv.
%         - clr: a string of color specifications. By default, the colors
%                are 'brkmgcy'.
%         - sym: a string or a cell of marker specifications. For example,
%                if sym = 'o+x', the first group will be plotted with a
%                circle, the second with a plus, and the third with a 'x'.
%                This is obtained with the assignment plo.sym = 'o+x'
%                or equivalently with plo.sym = {'o' '+' 'x'}.
%                By default the sequence of marker types is:
%                '+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'
%         - siz: scalar, a marker size to use for all plots. By default the
%                marker size depends on the number of plots and the size of
%                the figure window. Default is siz = '' (empty value).
%       - doleg: a string to control whether legends are created or not.
%                Set doleg to 'on' (default) or 'off'.
%
% dispopt: string which controls how to fill the diagonals in a plot of
%       Y vs Y (main diagonal of the scatter plot matrix). Set dispopt to
%       'hist' (default) to plot histograms, or 'box' to plot boxplots.
%
%       REMARK 1: the style which is used for univariate boxplots is
%       'traditional' if the number of groups is <=5, else it is 'compact'.
%
%       REMARK 2: to set dispopt without changing the defaults for plo use,
%       e.g., spmplot(Y,group,[],'box');
%
%
%  Y is a structure.
%
%  If first input Y is a structure (generally created by function FSMeda),
%  then this structure must have the following fields:
%
%       REQUIRED FIELD IN INPUT STRUCTURE Y.
%
%       Y.Y   = a data matrix of size n-by-v.
%
%               If the input structure Y contains just the data matrix, a
%               standard static scatter plot matrix will be created.
%
%               On the other hand, if Y also contains information on
%               statistics monitored along a search, then the scatter plots
%               will be linked with other (forward) plots with interaction
%               possibilities, enabled via brushing and datatooltip. More
%               precisely, with option databrush it is possible to create
%               an automatic interaction with the other plots, while with
%               option datatooltip it is possible to retrieve information
%               about a particular unit once selected with the mouse).
%
%       OPTIONAL FIELDS IN INPUT STRUCTURE Y.
%
%       Y.MAL = matrix containing the Mahalanobis distances monitored in each
%               step of the forward search. Every row is associated with a
%               unit (this is a necessary field if the user wants to brush
%               the scatter plot matrix).
%       Y.Un  = matrix containing the order of entry of each unit
%               (necessary if datatooltip is true or databrush is not
%               empty).
%       Y.label = cell of length n containing the labels of the units
%               (optional argument used when datatooltip=1. If this
%               field is not present labels row1, ..., rown will be
%               automatically created and included in the pop up
%               datatooltip window)
%
%       The options which follow can only be used if the first argument
%       of spmplot is a structure containing information about the fwd
%       search (i.e. the fields MAL, Un and eventually label)
%
%   Optional input arguments: (if the first argument of spmplot is a
%                               structure)
%
%
%   datatooltip :   interactive clicking. Empty value (default) or
%                   structure.
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
%    databrush :    interactive mouse brushing. Empty value (default),
%                   scalar or structure.
%                   DATABRUSH IS AN EMPTY VALUE
%                   If databrush is an empty value (default), no brushing
%                   is done.
%                   The activation of this option (databrush is a scalar or
%                   a cell) enables the user  to select a set of
%                   observations in the current plot and to see them
%                   highlighted in the malfwdplot, i.e. the plot of the
%                   trajectories of all observations, grouped according
%                   to the selection(s) done by brushing. If the malfwdplot
%                   does not exist it is automatically created.
%                   In addition, brushed units can be highlighted in the
%                   other following plots (only if they are already open):
%                   - minimum Mahalanobis distance plot;
%                   Remark. the window style of the other figures is set
%                   equal to that which contains the spmplot. In other
%                   words, if the scatterplot matrix plot is docked all the
%                   other figures will be docked too
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
%                   - labeladd. If this option is '1', we label the units
%                     of the last selected group with the unit row index in
%                     matrices X and y. The default value is labeladd='',
%                     i.e. no label is added.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%
%       The options which follow work in connection with previous option
%       databrush and produce their effect on the monitoring MD plot
%       (malfwdplot)
%
%       subsize :   x axis control in malfwdplot. Vector. numeric vector containing the subset size with length
%                   equal to the number of columns of matrix Y.MAL.
%                   If it is not specified it will be set equal to
%                   size(Y.MAL,1)-size(Y.MAL,2)+1:size(Y.MAL,1)
%                   Example - 'subsize',10:100
%                   Data Types - single | double
%       selstep :   add text labels of brushed units in malfwdplot. Vector. Numeric vector
%                   which specifies for which steps of the
%                   forward search textlabels are added in the monitoring
%                   MD plot after a brushing action in the spmplot.
%                   The default is to write the labels at the initial and
%                   final step. The default is selstep=[m0 n] where m0 and
%                   n are respectively the first and final step of the
%                   search.
%                   Example - 'selstep',100
%                   Data Types - single | double
%
%  Output:
%
%        H      :   array of handles H to the plotted points. 3D array. See
%                   gplotmatrix for further details
%        AX     :   handles to the individual subaxes. Matrix. See
%                   gplotmatrix for further details
%      BigAx    :   handle to big (invisible) axes framing the subaxes.
%                   Scalar. See gplotmatrix for further details
%
% More About:
%
%   spmplot has the same output of gplotmatrix in the statistics toolbox:
%   [H,AX,BigAx] = spmplot(...) returns an array of handles H to the
%   plotted points; a matrix AX of handles to the individual subaxes; and a
%   handle BIGAX to big (invisible) axes framing the subaxes.  The third
%   dimension of H corresponds to groups in G. AX contains one extra row of
%   handles to invisible axes in which the histograms are plotted. BigAx is
%   left as the CurrentAxes so that a subsequent TITLE, XLABEL, or YLABEL
%   will be centered with respect to the matrix of axes.
%
%
% See also: gplotmatrix
%
% Copyright 2008-2016.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('spmplot')">Link to the help function</a>
% Last modified 31-05-2016

% Examples:

%{
    % Call of spmplot without name/value pairs.
    % Iris data: scatter plot matrix with univariate boxplots on the main
    % diagonal
    load fisheriris;
    plo=struct;
    plo.nameY={'SL','SW','PL','PW'};
    spmplot(meas,species,plo,'hist');
%}

%{
    %% Call of spmplot with name/value pairs.
    % Iris data: scatter plot matrix with univariate boxplots on the main
    % diagonal
    load fisheriris;
    plo=struct;
    plo.nameY={'SL','SW','PL','PW'};
    spmplot(meas,'group',species,'plo',plo,'dispopt','box');
%}

%{
    % Iris data: scatter plot matrix with univariate boxplots on the main
    % diagonal and personalized options for symbols, colors and symbol
    % size and no legend
    load fisheriris;
    plo=struct;
    plo.nameY={'SL','SW','PL','PW'}; % Name of the variables
    plo.clr='kbr'; % Colors of the groups
    plo.sym={'+' '+' 'v'}; % Symbols of the groups (inside a cell)
    % Symbols can also be specified as characters
    % plo.sym='++v'; % Symbols of the groups
    plo.siz=3.4; % Symbol size
    plo.doleg='off'; % Remove the legend
    spmplot(meas,species,plo,'box');
%}
%
%
%{
    % Example of spmpot called by routine FSM.
    % Generate contaminated data.
    state=100;
    randn('state', state);
    n=200;
    Y=randn(n,3);
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;

    % spmplot is called automatically by all outlier detection methods, e.g. FSM
    [out]=FSM(Ycont,'plots',1);
%}

%{
    % Now test the direct use of FSM. Set two groups, e.g. those obtained
    % from FSM
    % Generate contaminated data
    state=100;
    randn('state', state);
    n=200;
    Y=randn(n,3);
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
   
    close all;
    [out]=FSM(Ycont,'plots',1);

    group = zeros(n,1);
    group(out.outliers)=1;
    plo=struct;
    plo.labeladd='1'; % option plo.labeladd is used to label the outliers

    % By default, the legend identifies the groups with the identifiers
    % given in vector 'group'.
    spmplot(Ycont,group,plo,'box');

%}

%{
    % spm with personalized tags.
    % With two groups, and if the Tag of the figure contains the word
    % 'outlier', the legend will identify one group for outliers and the
    % other for normal units. The largest number in the 'group' variable
    % identifies the group of outliers.
    figure('tag','This is a scatterplot with ouTliErs');
    spmplot(Ycont,group);

    % If the Tag of the Figure contains the string 'group', then the
    % legend identifies the groups with 'Group 1', Group 2', etc.
    figure('tag','This scatterplot contains groups');
    spmplot(Ycont,group,plo,'box');

    % If the tag figure includes the word 'brush', the legend will identify
    % one group for 'Unbrushed units' and the others for 'Brushed units 1',
    % 'Brushed units 2', etc.
    figure('Tag','Scatterplot with brushed units');
    spmplot(Ycont,group,plo);
%}


%{
    % An example with 5 groups.
    rng('default')
    rng(2); n1=100;
    n2=80;
    n3=50;
    n4=80;
    n5=70;
    v=5;
    Y1=randn(n1,v)+5;
    Y2=randn(n2,v)+3;
    Y3=rand(n3,v)-2;
    Y4=rand(n4,v)+2;
    Y5=rand(n5,v);

    group=ones(n1+n2+n3+n4+n5,1);
    group(n1+1:n1+n2)=2;
    group(n1+n2+1:n1+n2+n3)=3;
    group(n1+n2+n3+1:n1+n2+n3+n4)=4;
    group(n1+n2+n3+n4+1:n1+n2+n3+n4+n5)=5;


    Y=[Y1;Y2;Y3;Y4;Y5];
    spmplot(Y,group,[],'box');
%}

%{
    % spmplot called with name/pairs.
    % In all previous examples spmplot was called without the
    % name/value pairs arguments
    % The example which follow make use of the name/value pairs arguments
    load fisheriris;
    plo=struct;
    plo.nameY={'SL','SW','PL','PW'}; % Name of the variables
    plo.clr='kbr'; % Colors of the groups
    plo.sym={'+' '+' 'v'}; % Symbols of the groups (inside a cell)
    % Symbols can also be specified as characters
    % plo.sym='++v'; % Symbols of the groups
    plo.siz=3.4; % Symbol size
    spmplot(meas,'group',species,'plo',plo,'dispopt','box');
%}

% In the previous examples the first argument of spmplot was a matrix. In
% the two examples below the first argument is a structure which contains
% the fields Y and Un

%{
    % Example of use of option datatooltip.
    % First input argument is a structure.
    n=100;
    v=3;
    m0=3;
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:10,:)=5;
    [fre]=unibiv(Ycont);
    %create an initial subset with the 3 observations with the lowest
    %Mahalanobis Distance
    fre=sortrows(fre,4);
    bs=fre(1:m0,1);
    [out]=FSMeda(Ycont,bs,'plots',1);
    % mmdplot(out);
    figure
    plo=struct;
    plo.labeladd='1';
    spmplot(out,'datatooltip',1)
%}

%{
    % Interactive_example.
    % Example when first input argument is a structure.
    % Example of use of option databrush
    close all
    rng(841,'shr3cong');
    n=100;
    v=3;
    m0=v+1;
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [fre]=unibiv(Y);
    %create an initial subset with the 3 observations with the lowest
    %Mahalanobis Distance
    fre=sortrows(fre,4);
    bs=fre(1:m0,1);
    [out]=FSMeda(Ycont,bs,'plots',1);
    % mmdplot(out);
    figure
    plo=struct;
    plo.labeladd='1';
    % Please note the difference between plo.labeladd='1' and option labeladd
    % '1' inside databrush.
    % plo.labeladd enables the user to label the units in the scatterplot
    % matrix once selected. Option labeladd '1' inside databrush enables to add
    % the labels of the selected units in the linked plots
    spmplot(out,'databrush',{'persist','on','selectionmode' 'Rect','labeladd','1'},'plo',plo,'dispopt','hist')
%}

%% Beginning of code
% if length(varargin{1})==n then we are in the old format of the function
% spmplot(Y,group)
% or
% spmplot(Y,group,plo)
% or
% spmplot(Y,group,plo,dispopt)

if nargin<1
    error('FSDA:spmplot:missingInputs','A required input argument is missing.')
end

% Check if the first argument is a structure or not
if ~isstruct(Y);
    [n,v]=size(Y);
    onlyyX=1;
else
    % The first argument (Y) is a structure
    out=Y;
    % Retrieve matrix Y from the input structure
    Y=out.Y;
    [n,v]=size(Y);
    
    if find(strcmp('MAL',fieldnames(out)));
        % The number of rows of matrix MAL (squared Mahalanobis distnaces
        % along the search) is associated with the number of units. The
        % number of columns is associated with the number of steps of the
        % fwd search.
        residuals=out.MAL;
        [n,nsteps]=size(residuals);
        % onlyyX is a scalar. If onlyyX =1 this means that the user can do
        % brushing on the spmplot, else just a static spmplot with
        % clickable multilegends is produced
        onlyyX=0;
    else
        onlyyX=1;
    end
end

if nargin>1
    if length(varargin{1})==n
        % In this case the user has called function spmplot with the
        % old format, that is
        % spmplot(Y,group,plo,dispopt), without name/value pairs
        
        group=varargin{1};
        
        if length(varargin)<2
            plo=1;
        else
            plo=varargin{2};
        end
        
        if length(varargin)<3
            dispopt='hist';
        else
            dispopt=varargin{3};
        end
        databrush='';
        datatooltip=0;
        tag='pl_spm';
    else
        % In the case the user has called function spmplot with the new
        % format name/value pairs
        namevaluepairs=1;
        
        if onlyyX==0;
            % x= vector which contains the subset size (numbers on the x axis)
            x=(n-nsteps+1):n;
            
            % selthdef= threshold used to decide which residuals are labelled in the resfwdplot.
            % laby= label used for the y-axis of the resfwdplot.
            labx='Subset size m';
            laby='Mahalanobis distances';
            selthdef=num2str(4*v);
            
            % maximum and minimum residual along the search for each unit
            selmax=max(residuals,[],2);
            selmin=min(residuals,[],2);
        else
            x=1;
            % The default is not to add textlabels to any unit
            selthdef='';
        end
        one=ones(n,1);
        options=struct('group',one,'plo',[],'subsize',x,'selstep',x([1 end]),...
            'selunit',selthdef,'datatooltip',0,'label','',...
            'dispopt','hist','databrush','','tag','pl_spm');
        
        UserOptions=varargin(1:2:length(varargin));
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:spmplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:spmplot:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
        
        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin);
            options.(varargin{i})=varargin{i+1};
        end
        group=options.group;
        plo=options.plo;
        dispopt=options.dispopt;
        databrush=options.databrush;
        tag=options.tag;
        datatooltip=options.datatooltip;
    end
else
    group=ones(n,1);
    plo='';
    dispopt='hist';
    tag='pl_spm';
    datatooltip=0;
    databrush='';
    namevaluepairs=1;
end

ngroups=length(unique(group));
% Specify default values for colors, symbols, size of symbols and presence
% of legend
clrdef='brkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcy';
symdef={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.';'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
sizdef='';
dolegdef='on';

% seq= column vector containing the sequence 1 to n
seq= (1:n)';

if isstruct(plo)
    fplo=fieldnames(plo);
    
    d=find(strcmp('nameY',fplo));
    if d>0
        nameY=plo.nameY;
    else
        nameY=cellstr(num2str((1:v)','Y%d'));
    end
    d=find(strcmp('labeladd',fplo));
    if d>0
        labeladd=plo.labeladd;
    else
        labeladd='';
    end
    d=find(strcmp('clr',fplo));
    if d>0
        clr=plo.clr;
        
        if length(clr) ~= ngroups
            warning('FSDA:spmplot:WrongNumColors','Number of colors which have been supplied is not equal to the number of groups')
            disp(['Number of groups =' num2str(ngroups)])
            disp(['Number of colors =' num2str(length(clr))])
            if length(clr)< ngroups
                disp('Supplied colors will be duplicated')
                if isrow(clr)
                    clr=repmat(clr,1,ngroups);
                else
                    clr=repmat(clr,ngroups,1);
                end
            else
                disp(['Just the first ' num2str(ngroups) ' colors will be used'])
            end
        end
        
    else
        clr=clrdef;
    end
    d=find(strcmp('sym',fplo));
    if isempty(d)
        d=0;
    end
    if d>0 % && (ngroups == numel(plo.sym))
        sym=plo.sym;
        if length(sym) ~= ngroups
            warning('FSDA:spmplot:WrongNumSymb','Number of symbols which have been supplied is not equal to the number of groups')
            disp(['Number of groups =' num2str(ngroups)])
            disp(['Number of symbols =' num2str(length(sym))])
            if length(sym)< ngroups
                disp('Supplied symbols will be duplicated')
                if isrow(sym)
                    sym=repmat(sym,1,ngroups);
                else
                    sym=repmat(sym,ngroups,1);
                end
            else
                disp(['Just the first ' num2str(ngroups) ' symbols will be used'])
            end
        end
    else
        sym=symdef;
    end
    d=find(strcmp('siz',fplo));
    if d>0
        siz=plo.siz;
    else
        siz=sizdef;
    end
    d=find(strcmp('doleg',fplo));
    if d>0
        doleg=plo.doleg;
    else
        doleg='on';
    end
    
else
    
    if ischar(plo) && namevaluepairs==0
        error('FSDA:spmplot:InvalidArg3',' Third argument must be a structure, or a scalar or an empty value []')
    end
    
    if plo==1
        nameY = cellstr(num2str((1:v)','Y%d'));
    else
        nameY = '';
    end
    labeladd='';
    
    clr=clrdef;
    sym=symdef;
    siz=sizdef;
    doleg=dolegdef;
end

if iscell(group)
    groupv = zeros(numel(group),1);
    guni = unique(group,'stable');
    for ii=1:numel(guni)
        groupv(strcmp(group,guni(ii))) = ii;
    end
    %guniv = cellstr(num2str(unique(groupv,'stable')));
else
    groupv = group;
end

unigroup = 1:ngroups;

%% The scatterplot matrix with histograms or boxplots (on the main diagonal) generalised to groups

% sym can be either a cell array or a character
if iscell(sym)
    charsym=char(sym{unigroup});
else
    charsym=sym(unigroup);
end

% The scatter matrix is generated with histograms on the main diagonal and
% the axes matrix AX will therefore contain handles for these histograms;
% later we superimpose new group histograms using the FSDA function histFS;
% we do not use gplotmatrix with option 'grpbars' because this would not
% work in MATLAB releases previous to R2015a.

[H,AX,BigAx] = gplotmatrix(Y,[],group,clr(unigroup),charsym,siz,doleg,'hist',nameY,nameY);

% Not necessary anymore: lines will be removed after additional testing.
% if strcmp(dispopt,'box')==1
%     % ahist=findobj(AX,'type','patch') ;
%     %set(ahist,'MarkerEdgeColor','none','EdgeColor','none','FaceColor','none')
% ahist=findobj(AX,'type','histogram');
% set(ahist,'Visible','off')
% end

% The third dimension of H distinguishes the groups. If there are no groups
% then ndims(H) = 2.

for i=1:size(AX,2)
    hold('on');
    
    if strcmp(dispopt,'hist')==1
        % Add the histograms generalised to groups
        
        ax = AX(i,i);
        Xlim=get(ax,'Xlim');
        
        ax = AX(end,i);
        
        % Not necessary anymore: lines will be removed after additional testing.
        % Ylim=get(ax,'Ylim');
        % set(gcf,'CurrentAxes',ax);
        
        % the strings used to label the tick marks
        XTickLabel = get(ax,'XTickLabel');
        YTickLabel = get(ax,'YTickLabel');
        XLabel = get(get(ax,'XLabel'),'String');
        YLabel = get(get(ax,'YLabel'),'String');
        
        [countfreq, ~] = histFS(Y(:,i),10,groupv,'',ax,clr(unigroup)); %'br'
        
        % Prevent from changing the limits when the figure is resized:
        % 1.Freeze the current limits
        set(ax,'XLimMode','manual','YLimMode','manual');
        set(ax,'xlim',Xlim)
        % set(ax,'ylim',Ylim)
        set(ax,'ylim',[0 max(sum(countfreq,2))])
        % 2.Freeze the current tick values
        set(ax,'XTickMode','manual','YTickMode','manual');
        %Now restore the labels of the gplotmatrix
        set(ax,'XTickLabel',XTickLabel,'YTickLabel',YTickLabel);
        set(get(ax,'XLabel'),'String',XLabel);
        set(get(ax,'YLabel'),'String',YLabel);
        
    else % if strcmp(dispopt,'box')==1
        % Add the boxplots generalised to groups. Note that we use AX(i,i)
        % just to find the position of the required panel on the main
        % diagonal to superimpose boxplots
        ax=AX(i,i);
        
        axPosition = get(ax,'position');
        % Now we create an axes object using axPosition.
        ax = axes('Position',axPosition);
        
        if length(unigroup) <= 5
            plotstyle = 'traditional';
        else
            plotstyle = 'compact';
        end
        
        % Boxplots are superimposed on object ax
        hbp = boxplot(ax,Y(:,i),groupv,'plotstyle',plotstyle,'colors',clr(unigroup),'labelverbosity','minor','symbol','+');
        
        % Remove the x tick labels from the graph containing boxplots
        set(ax,'XTickLabel',{' '});
        
        % Remove the y tick labels from the graph containing boxplots
        set(ax,'YTickLabel',{' '});
        
        % Put the graph containing boxplots in the correct position
        set(ax,'position',axPosition);
        
        % Not necessary anymore: lines will be removed after additional testing.
        % Adjust the vertical scale (using the y scale of a scatter in
        % the same row of the scatter plot matrix)
        %         if i < size(Y,2);
        %             ax1 = AX(i,end);
        %             Ylim=get(ax1,'Ylim');
        %         else
        %             ax1 = AX(i,1);
        %             Ylim=get(ax1,'Ylim');
        %         end
        %         set(ax,'Ylim',Ylim)
        
        % The tag is set for later use in add2spm by clickableMultiLegend
        for gg=1:numel(unigroup)
            set(hbp(:,gg),'Tag',['boxplot' num2str(gg)]);
        end
    end
end

if ndims(H) == 3
    H=double(H);
    H(:,:,end) = ~eye(size(H,1)).*H(:,:,end);
    
    % Put in the figure UserData field the list of units in the last group,
    % i.e. (depending on context) the outliers or the units last brushed.
    if isnumeric(group)
        set(H(:,:,end), 'UserData' , seq(group==max(group)));
    else
        set(H(:,:,end), 'UserData' , seq(groupv==max(groupv)));
    end
    
    if strcmp(doleg,'on')
        % Add to the spm the clickable multilegend and eventually the text labels
        % of the selections
        if isnumeric(group)
            % use context sensitive legends
            add2spm(H,AX,BigAx,'labeladd',labeladd,'userleg','1');
        else
            % use legends in guni
            add2spm(H,AX,BigAx,'labeladd',labeladd,'userleg',guni);
        end
    end
end

% the handle of the figure including the gplotmatrix (i.e. the closest ancestor of BigAx).
fig = ancestor(BigAx,'figure');

% Make BigAx the CurrentAxes
% The instruction below is necessary otherwise for example the instruction
% title will not be put on top of the screen but on top of the last panel
% on the bottom right
set(gcf,'CurrentAx',BigAx)

% Also set Title and X/YLabel visibility to on and strings to empty
set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
    'String','','Visible','on')


% set the options.datatooltip (enable/disable interactive data cursor mode)
if datatooltip;
    
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
    set(hdt,'UpdateFcn',{@spmplotLbl,out})
end

%% Brush mode (call to function selectdataFS)

% BEGINNING OF DATABRUSH OPTION
if ~isempty(databrush) || iscell(databrush)
    % Initialize line width
    linewidthStd = 0.5;
    units=options.selunit;
    
    % extract the vector associated with the subset size (x)
    x=options.subsize;
    % and check if the choice of selsteps is valid
    steps=options.selstep;
    if max(steps)>n
        mess=sprintf(['One of the steps which has beeen chosen is greater than n. \n',...
            'It is deleted.']);
        disp(mess);
        steps=steps(steps<=n);
    end
    if min(steps)<x(1)
        mess=sprintf(['One of the steps which has beeen chosen is smaller than m0. \n',...
            'It is deleted.']);
        disp(mess);
        steps=steps(steps>=x(1));
    end
    
    if onlyyX;
        
        if ~isempty(options.databrush)
            disp('It is not possible to use option databrush without supplying structure out produced by FSReda')
            return
        end
        if iscellstr(units)
            selunit=str2double(units);
            units=seq(y>selunit(2) | y<selunit(1));
        elseif ischar(units)==1;
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
                disp(mess);
                units=units(units<=n);
            end
            if min(units)<1
                mess=sprintf(['One of the units which has beeen chosen is smaller than 1. \n',...
                    'It is deleted.']);
                disp(mess);
                units=units(units>0);
            end
        end
        
    else
        
        if iscellstr(units)
            selunit=str2double(units);
            selmax=max(residuals,[],2);
            selmin=min(residuals,[],2);
            units=seq(selmax>selunit(2) | selmin<selunit(1));
        elseif ischar(units)==1;
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
                disp(mess);
                units=units(units<=n);
            end
            if min(units)<1
                mess=sprintf(['One of the units which has beeen chosen is smaller than 1. \n',...
                    'It is deleted.']);
                disp(mess);
                units=units(units>0);
            end
        end
    end
    
    % lunits = number of units which must be labelled
    lunits=length(units);
    % lsteps = number of steps for which it is necessary to add the labels
    lsteps=length(steps);
    lall=lunits*lsteps;
    
    % numtext= a cell of strings used to labels the units with their position
    % in the dataset.
    numtext=cellstr(num2str(seq,'%d'));
    
    %%  Prepare the spmplot for brushing
    
    plot1 = fig;
    %plot1=gcf;
    
    %when the graphics objects are axes children, MATLAB clears the axes and
    %resets some of their properties to default values. In order to add new
    %graphics objects without clearing or resetting the current figure we set
    %the figure and axes NextPlot properties to 'add'.
    set(fig,'NextPlot','add');
    set(AX,'NextPlot','add');
    
    % Set default value for potential groups of selected units
    styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
    styp=repmat(styp,ceil(n/13),1);
    
    % displays the boundary of the current axes.
    box on
    
    % set the specified tag in the current plot
    set(gcf,'tag',tag)
    
    % labeladd option
    d=find(strcmp('labeladd',databrush));
    if d>0
        labeladdDB=databrush(d+1);
        labeladdDB=labeladdDB{1};
        % This option must be removed from cell options.databrush because it is
        % not a valid option for the function selectdataFS.
        databrush(d:d+1)=[];
    else
        labeladdDB='';
    end
    
    % persist option
    d=find(strcmp('persist',databrush));
    if d>0
        persist=databrush(d+1);
        % This option must be removed from cell options.databrush because it is
        % not a valid option for the function selectdataFS.
        databrush(d:d+1)=[];
        ColorOrd=[1 0 0;0 1 1; 1 0 1; 1 1 0; 0 0 0; 0 1 0; 0 0 1];
        ColorOrd=repmat(ColorOrd,4,1);
    else
        persist='';
        ColorOrd=[1 0 0];
    end
    
    % FlagColor option
    % Initialize colors: default colors are blue (unbrushed unit) and red
    % (brushed units)
    d=find(strcmp('FlagColor',databrush));
    if d>0
        flagcol=databrush{d+1};
        clr=['b' flagcol 'cmykgbrcmykg'];
    else
        clr='brcmykgbrcmykgbrcmykg';
    end
    
    if isscalar(databrush)
        sele={'selectionmode' 'Rect' 'Ignore' findobj(gcf,'tag','env') };
    else
        sele={databrush{:} 'Ignore' findobj(gcf,'tag','env')}; %#ok<*CCAT>
    end
    
    sele={sele{:} 'Tag' tag};
    
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
        
        if strcmp(persist,'off');
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
        
        %% - select an area in the spmplot
        
        % A function to be executed on figure close is set
        set(gcf,'CloseRequestFcn',@closereqFS);
        
        disp('Select a region to brush in the spmplot');
        disp('Left mouse button picks points.');
        ss=waitforbuttonpressFS;
        
        % After waitforbuttonpress:
        % - the standard MATLAB function to be executed on figure
        %   close is recovered
        set(gcf,'CloseRequestFcn','closereq');
        Open_YY = findobj(0, 'type', 'figure','tag','pl_spm');
        Open_res = findobj(0, 'type', 'figure','tag','data_res');
        Open_mmd = findobj(0, 'type', 'figure','tag','pl_mmd');
        if isempty(Open_YY)  % User closed the main brushing window
            if ~isempty(Open_res); delete(Open_res); end    % spm plot is deleted
            if ~isempty(Open_mmd); delete(Open_mmd); end  % mmd plot is deleted
            delete(get(0,'CurrentFigure')); % deletes Figure if still one left open
        end
        
        if ss==1;
            but=2;
        else but=1;
        end
        
        if but==1 && ~isempty(find(AX==gca,1))  %if the User has made a selection in one of the scatterplots
            
            %% - identify the axes of the area
            
            % fig is the closest ancestor figure of BigAx. Set property
            % CurrentAxes equal to the axes of the scatterplot that we have
            % selected
            set(fig,'CurrentAxes',gca);
            %indice is a scalar identyfing the selected axes
            indice=find(AX==gca);
            
            % indicer and indicec respectively contain the row and colum
            % indexes of the scatter in which points have been selected
            [indicer,indicec]=ind2sub(size(AX),indice);
            
            otherAxes = AX;
            
            %% - call selectdataFS
            if ij>1;
                set( hLegend(1,end),'Visible','off')
            end
            
            set(fig,'UserData',num2cell(seq));
            [pl,xselect,yselect] = selectdataFS(sele{:}, 'Label', 'off');
            
            % exit from yXplot if the yXplot figure was closed before selection
            if ~isempty(pl) && isnumeric(pl) && (min(pl) == -999)
                return
            end
            
            if ij>1;
                set(hLegend(1,end),'Location', 'Best');
            end
            
            %When the selection has been completed, axes properties
            %'HandleVisibility' and 'HitTest' must be set to on for an eventual
            %possible future selection.
            set(otherAxes,'HandleVisibility','on');
            set(otherAxes,'HitTest','on');
            
            %% - identify the selected units
            if iscell(xselect)
                %xselect is a cell if more than one selection has been
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
            coi=NaN(1,length(xselect_all));
            
            for k=1:length(xselect_all)
                %uno=vector of zeros and ones of the same length of X and
                %y. The ones are reported in the rows where the X values
                %are equal to the X selected
                uno=Y(:,indicec)==(xselect_all(k));
                %due=same idea of uno: vector of zeros and ones of the same
                %length of X and y. The ones are reported in the rows where
                %the y values are equal to the y selected
                due=Y(:,indicer)==(yselect_all(k));
                %find observations that have satisfied both (AND)
                %conditions uno and due, i.e. observations with both X and
                %y equal to the X and y selected.
                tre=uno + due;
                %coi=list of selected units
                coi(k)=find(tre==2);
                if k==1
                    coi=coi(k);
                end
            end
            
            %pl=vector with n rows. Ones represent observations selected at
            %the current iteration, zeros all the others.
            pl=zeros(n,1);
            pl(coi)=1;
            if ~isempty(pl)
                % sel = vector which contains the list of units not
                % selected in the current iteration.
                sel=seq(pl<1);
                
                % nbrush = vector which contains the list of the units selected
                % in the current iteration.
                nbrush=setdiff(seq,sel);
                disp('Brushed units are');
                disp([nbrush out.Y(nbrush)]);
            else
                disp('Wrong selection: Try again');
                disp('Select a region to brush in one of the scatters out of the diagonal ');
                figure(plot1);
                nbrush='';
            end
            
            if ~isempty(pl)
                disp('Brushed units are');
                disp([nbrush Y(nbrush,:)]);
            else
                disp('Wrong selection: Try again');
                disp('Select a region to brush in the scatter plot matrix');
                figure(plot1);
                nbrush='';
            end
            
            %% For each brushing operation, do the following:
            if ~isempty(nbrush)
                
                %% - Update the spmplot and keep it always in foreground
                
                %brushcum is the list of selected observations in all
                %iterations if persist=on or the list of selected
                %observations in the current iteration if persist=off
                if strcmp(persist,'on')
                    brushcum=unique([brushcum; nbrush]);
                else
                    brushcum=nbrush;
                    group=ones(n,1);
                end
                
                % group=vector of length size(Y,1) taking values
                % from 1 to the number of groups selected.
                % unigroup= list of selected groups.
                group(nbrush)=ij+1;
                unigroup=unique(group);
                groupv=group;
                %open the figure containing the scatterplot matrix
                figure(plot1);
                
                % set in this figure, containing the gplotmatrix, the
                % property currentAxes equal to the axes of the scatterplot
                % in which the selection has been done.
                set(plot1,'CurrentAxes',AX(indicer,indicec));
                
                %replace the gplotmatrix which has come out from
                %selectdataFS
                %(the selection was done ONLY on one scatterplot) with a
                %new gplotmatrix which takes into consideration the
                %selection coming out from selectdataFS, in all scatterplots
                [H,AX,BigAx] = gplotmatrix(Y,[],group,clr(unigroup),char(styp{unigroup}),[],'on','hist',nameY,nameY);
                
                % Set the legend properties of the gplotmatrix
                if nbrush>0
                    set(H(:,:,1),'DisplayName','Unbrushed units');
                    for brugrp = 2:length(unigroup)
                        set(H(1,end,brugrp),'DisplayName',['Brushed units ' num2str(brugrp-1)]);
                    end
                else
                    set(H,'DisplayName','Units');
                end
                
                % Make the histograms on the main diagonal invisible if the
                % user has decided to put boxplots
                if strcmp(dispopt,'box')==1
                    ahist=findobj(AX,'type','patch') ;
                    set(ahist,'MarkerEdgeColor','none','EdgeColor','none','FaceColor','none')
                end
                
                % beginning of updating boxplots or histograms
                % on the main diagonal
                for i=1:size(AX,2)
                    hold('on');
                    
                    if strcmp(dispopt,'hist')==1
                        % Add the histograms generalised to groups
                        
                        ax = AX(i,i);
                        Xlim=get(ax,'Xlim');
                        
                        ax = AX(end,i);
                        Ylim=get(ax,'Ylim');
                        
                        % the strings used to label the tick marks
                        XTickLabel = get(ax,'XTickLabel');
                        YTickLabel = get(ax,'YTickLabel');
                        XLabel = get(get(ax,'XLabel'),'String');
                        YLabel = get(get(ax,'YLabel'),'String');
                        
                        [countfreq, ~] = histFS(Y(:,i),10,groupv,'',ax,clr(unigroup)); %'br'
                        % Prevent from changing the limits when the figure is resized:
                        % 1.Freeze the current limits
                        set(ax,'XLimMode','manual','YLimMode','manual');
                        set(ax,'xlim',Xlim)
                        set(ax,'ylim',Ylim)
                        
                        set(ax,'ylim',[0 max(sum(countfreq,2))])
                        
                        % 2.Freeze the current tick values
                        set(ax,'XTickMode','manual','YTickMode','manual');
                        %Now restore the labels of the gplotmatrix
                        set(ax,'XTickLabel',XTickLabel,'YTickLabel',YTickLabel);
                        set(get(ax,'XLabel'),'String',XLabel);
                        set(get(ax,'YLabel'),'String',YLabel);
                        
                    else
                        % Modify the boxplots generalised to groups Note
                        % that instead of using ax = AX(i,i); we use
                        % AX(end,i); because the last row of AX contains
                        % the handles to the invisible axes in which the
                        % histograms or boxplots are plotted
                        ax=AX(end,i);
                        
                        % get the position of AX(i,i)
                        axPosition = get(ax,'position');
                        
                        if length(unigroup) <= 5
                            plotstyle = 'traditional';
                        else
                            plotstyle = 'compact';
                        end
                        
                        hbp = boxplot(ax,Y(:,i),groupv,'plotstyle',plotstyle,'colors',clr(unigroup),'labelverbosity','minor','symbol','+');
                        
                        % Remove the x tick labels from the graph containing boxplots
                        set(ax,'XTickLabel',{' '});
                        
                        % Remove the y tick labels from the graph containing boxplots
                        set(ax,'YTickLabel',{' '});
                        
                        % Put the graph containing boxplots in the correct position
                        set(ax,'position',axPosition);
                        
                        % Adjust the vertical scale (using the y scale of a scatter in
                        % the same row of the scatter plot matrix)
                        if i < size(Y,2);
                            ax1 = AX(i,end);
                            Ylim=get(ax1,'Ylim');
                        else
                            ax1 = AX(i,1);
                            Ylim=get(ax1,'Ylim');
                        end
                        set(ax,'Ylim',Ylim)
                        
                        % The tag is set for later use in add2spm by
                        % clickableMultiLegend
                        for gg=1:numel(unigroup)
                            set(hbp(:,gg),'Tag',['boxplot' num2str(gg)]);
                        end
                        
                    end
                end
                % end of updating boxplots or histograms on the main
                % diagonal
                
                %fig is the handle of the closest ancestor figure of BigAx.
                fig = ancestor(BigAx,'figure');
                
                %In order to add new graphics objects without clearing or
                %resetting the current figure we set the figure and axes
                %NextPlot properties to 'add'.
                set(fig,'NextPlot','add');
                set(AX,'NextPlot','add');
                
                % Now update the legends and make them clickable.
                if strcmp(doleg,'on')
                    H=double(H);
                    H(:,:,end) = ~eye(size(H,1)).*H(:,:,end);
                    
                    if ~isempty(labeladd)
                        % Put in the figure UserData field the list of units in the last group,
                        % i.e. the units last brushed.
                        if isnumeric(group)
                            set(H(:,:,end), 'UserData' , seq(group==max(group)));
                        else
                            set(H(:,:,end), 'UserData' , seq(groupv==max(groupv)));
                        end
                    end
                    
                    % Add to the spm the clickable multilegend and eventually the text labels
                    % of the selections
                    % use context sensitive legends
                    add2spm(H,AX,BigAx,'labeladd',labeladd,'userleg','1');
                end
                
                hLines = findobj(AX(1,end), 'type', 'line');
                eLegend = cell(length(hLines), 1);
                for iLines = 1:length(hLines)
                    eLegend{iLines} = get(hLines(iLines), 'DisplayName');
                end
                hLegend=zeros(size(AX));
                hLegend(1,end) = clickableMultiLegend(hLines, eLegend{:});
                
                %% - display the malfwdplot with the corresponding groups of trajectories highlighted.
                
                % creates the resfwdplot
                fres=findobj('-depth',1,'Tag','data_res');
                if isempty(fres)
                    figure;
                    % The window style of the new figure which has been
                    % created is set equal to that which contained the
                    % yXplot
                    set(gcf,'WindowStyle',get(plot1,'WindowStyle'));
                    set(gcf,'Name','malfwdplot','NumberTitle', 'off');
                    
                    % store the LineWidth property of the selected trajectories
                    plot(x,residuals,'Tag','data_res','Color','b','LineWidth',0.5);
                    
                    set(gcf,'tag','data_res');
                    hold('off')
                    
                    % displays the boundary of the current axes.
                    box on
                end
                
                %check if the malfwdplot is already open
                hh=findobj('-depth',1,'tag','data_res');
                figure(hh);
                set(gcf,'WindowStyle',get(plot1,'WindowStyle'));
                
                if strcmp('off',persist)
                    %if persist=off the previous selection must be replaced
                    set(gcf,'NextPlot','replace');
                    plot(x,residuals,'Tag','data_res','Color','b');
                    text(reshape(repmat(steps,lunits,1),lall,1),reshape(residuals(units,steps-x(1)+1),lall,1),reshape(repmat(numtext(units),1,lsteps),lall,1));
                    set(gcf,'tag','data_res');
                    set(gcf,'Name','malfwdplot');
                    set(gcf,'NextPlot','add');
                end
                hold on
                if strcmp('on',persist);
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
                if strcmp('1',labeladdDB)
                    if strcmp('off',persist)
                        text(reshape(repmat(steps,length(nbrush),1),length(nbrush)*length(steps),1),reshape(residuals(nbrush,steps-x(1)+1),length(nbrush)*length(steps),1),reshape(repmat(numtext(nbrush),1,length(steps)),length(nbrush)*length(steps),1));
                    end
                    if strcmp('on',persist)
                        text(reshape(repmat(steps,length(brushcum),1),length(brushcum)*length(steps),1),reshape(residuals(brushcum,steps-x(1)+1),length(brushcum)*length(steps),1),reshape(repmat(numtext(brushcum),1,length(steps)),length(brushcum)*length(steps),1));
                    end
                end
                %add (fix) labels eventually set in selunit.
                if ~isempty(units)
                    text(reshape(repmat(steps,length(units),1),length(units)*length(steps),1),...
                        reshape(residuals(units,steps-x(1)+1),length(units)*length(steps),1),...
                        reshape(repmat(numtext(units),1,length(steps)),length(units)*length(steps),1));
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
                        % REAMRK: mod(row,size(Un,1)) is used rather than
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
                
                % Check if the figure containing minimum MD
                % is open If it is, units are brushed in this plot too
                h=findobj('-depth',1,'Tag','pl_mmd');
                
                if (~isempty(h))
                    
                    % Remove unnecessary rows from vector selsteps -1 is
                    % necessary because we are considering minimum outside
                    selsteps=selsteps(:,1)-1;
                    
                    % make figure which contains mdr become the current figure
                    figure(h);
                    
                    % Condition || but==0 if but=0 then it is necessary to
                    % remove previous highlightments (even if persist='on')
                    if strcmp(persist,'off') || but==0;
                        % If set of values has already been highlighted in the mdr plot,
                        % remove it
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
                    
                    disp('Highlight the spmplot, then click on it to continue brushing or press a keyboard key to stop')
                    % The monitoring residual plot is highlighted again
                    % before continuing with waitforbuttonpress
                    figure(plot1);
                    set(gcf,'CloseRequestFcn',@closereqFS);
                    ss=waitforbuttonpressFS;
                    set(get(0,'CurrentFigure'),'CloseRequestFcn','closereq');
                    disp('------------------------')
                    
                    set(gcf,'CloseRequestFcn','closereq');
                    Open_YY = findobj(0, 'type', 'figure','tag','pl_spm');
                    Open_res = findobj(0, 'type', 'figure','tag','data_res');
                    Open_mmd = findobj(0, 'type', 'figure','tag','pl_mmd');
                    if isempty(Open_YY)  % User closed the main brushing window
                        if ~isempty(Open_res); delete(Open_res); end    % spm plot is deleted
                        if ~isempty(Open_mmd); delete(Open_mmd); end  % mmd plot is deleted
                        delete(get(0,'CurrentFigure')); % deletes Figure if still one left open
                    end
                    
                    if ss==1;
                        but=2;
                    end
                else
                    but=2;
                end
            end
        end
    end
end

    function output_txt = spmplotLbl(~,event_obj,out)
        %% spmplotLbl provides information about the selected points
        %
        % Required input arguments:
        %
        %   obj     =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the sense that
        %               these arguments are automatically passed to the function when it executes.
        %       out =   a structure containing the following fields
        %         Y =   Data matrix (2D array) containing n observations on v variables.
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
        % REMARK: this function is called by function spmplot
        %
        % References:
        %
        %   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
        %   Springer Verlag, New York.
        %
        % Written by FSDA team
        
        %% - identify the axes of the area
        output_txt{1}='oooo';
        
        % fig is the closest ancestor figure of BigAx. Set property CurrentAxes
        % equal to the axes of the scatterplot that we have selected
        set(fig,'CurrentAxes',gca);
        
        %TODO:spmplot:ControlAxes
        
        pos = get(event_obj,'Position');
        
        % x and y, plot coordinates of the mouse
        x = pos(1); y = pos(2);
        
        if isempty(x)
            output_txt{1}=['no observation has coordinates x,y' num2str(x) '' num2str(y)] ;
        else
            
            % Find index to retrieve obs. name
            % Consider that find return the linear indexing of matrix xydata
            ord=out.Y;
            
            set(fig,'CurrentAxes',gca);
            %indice is a scalar identyfing the selected axes
            ind=find(AX==gca);
            
            % indicer and indicec respectively contain the row and colum
            % indexes of the scatter in which points have been selected
            [indr,indc]=ind2sub(size(AX),ind);
            
            %uno=vector of zeros and ones of the same length of X and
            %y. The ones are reported in the rows where the X values
            %are equal to the X selected
            uno=ord(:,indc)==x;
            %due=same idea of uno: vector of zeros and ones of the same
            %length of X and y. The ones are reported in the rows where
            %the y values are equal to the y selected
            due=ord(:,indr)==y;
            %find observations that have satisfied both (AND)
            %conditions uno and due, i.e. observations with both X and
            %y equal to the X and y selected.
            tre=uno + due;
            %coi=list of selected units
            row=find(tre==2);
            
            % If structure out does not contain labels for the rows then labels
            % row1....rown are added automatically
            if isempty(intersect('label',fieldnames(out)))
                out.label=cellstr(num2str((1:size(out.Y,1))','row %d'));
            end
            
            output_txt=cell(length(row)*2+2,1);
            % output_txt is what it is shown on the screen
            output_txt(1) = {['Y(,:' num2str(indr) ') value equal to: ',num2str(y,4)]};
            
            output_txt(2) = {['Y(,:' num2str(indc) ') value equal to: ',num2str(x,4)]};
            
            % Add information about the corresponding row label of
            % what has been selected
            
            % Add information about the step in which the selected unit entered
            % the search
            ij=3;
            Un=out.Un;
            for kk=1:length(row)
                output_txt{ij} = ['Unit: ' cell2mat(out.label(row(kk)))];
                
                idx = find(Un(:,2:end) == row(kk));
                [rw,~] = ind2sub(size(Un(:,2:end)),idx);
                
                if isempty(rw)
                    output_txt{ij+1} = ['Unit entered before step m=' num2str(Un(1,1))];
                elseif length(rw)<2;
                    output_txt{ij+1} = ['Unit entered in step m=' num2str(Un(rw,1))];
                elseif length(rw)==2;
                    output_txt{ij+1} = ['Unit entered in step m=' num2str(Un(rw(1),1)) ' and then in step m=' num2str(Un(rw(2),1))];
                else
                    output_txt{ij+1} = ['Unit entered in steps m=' num2str(Un(rw(1),1)) ', m=' num2str(Un(rw(2),1)) ' and m=' num2str(Un(rw(3),1))];
                end
                ij=ij+2;
            end
            
            %         %When the selection has been completed, axes properties
            %         %'HandleVisibility' and 'HitTest' must be set to on for an eventual
            %         %possible future selection.
            %         set(otherAxes,'HandleVisibility','on');
            %         set(otherAxes,'HitTest','on');
            
        end
    end
end
%FScategory:VIS-Mult
