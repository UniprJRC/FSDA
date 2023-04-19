function [h,Ntable] = balloonplot(N,varargin)
%balloonplot creates a balloon plot of a contingency table
%
%<a href="matlab: docsearchFS('balloonplot')">Link to the help function</a>
%
% A balloonplot is an alternative to bar plot for visualizing
% a large categorical data. It draws a graphical matrix of a contingency
% table, where each cell contains a dot whose size reflects the relative
% magnitude of the corresponding component. ballonplot is a bubble chart
% for contingency tables.
%
%  Required input arguments:
%
%       N    :    Contingency table (default) or n-by-2 input dataset.
%                 2D Array or Table.
%                 2D array or table which contains the input contingency
%                 table (say of size I-by-J) or the original data matrix X.
%                 In this last case N=crosstab(X(:,1),X(:,2)). As default
%                 procedure assumes that the input is a contingency table.
%
%  Optional input arguments:
%
%      ax     :  displays the bubble chart in the target axes ax.
%                Axes object.
%                If you do not specify the axes, MATLAB plots into the
%                current axes, or it creates an Axes object if one does
%                not exist.
%               Example - 'ax',myaxes
%               Data Types - Axes object
%
% datamatrix  : Data matrix or contingency table. Boolean. If
%               datamatrix is true the first input argument N is forced to
%               be interpreted as a data matrix, else if the input argument
%               is false N is treated as a contingency table. The default
%               value of datamatrix is false, that is the procedure
%               automatically considers N as a contingency table (in array
%               or table format). If datamatrix is true, N can be an array
%               or a table of size n-by-2. Note that if N has more than two
%               columns balloonplot is based on the first two
%               columns of N (and a warning is produced).
%               Example - 'datamatrix',true
%               Data Types - logical
%
%       Lr   :  Vector of row labels. Cell.
%               Cell containing the labels of the rows of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table, because in this case  Lr=N.Properties.RowNames;
%               Example - 'Lr',{'a' 'b' 'c'}
%               Data Types - cell array of strings
%
%       Lc   :  Vector of column labels. Cell.
%               Cell containing the labels of the columns of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table, because in this case Lc=N.Properties.VariableNames;
%               Example - 'Lc',{'c1' c2' 'c3' 'c4'}
%               Data Types - cell array of strings
%
%
%  Output:
%
%      h :    returns the BubbleChart object. Use h to modify properties of
%             the chart after creating it. For a list of properties, see
%             BubbleChart Properties.
%
%    Ntable:  $I$-by-$J$ contingency table.
%             Table. This is tha table which has been used to build the
%             balloonplot.
%
% See also bubblesize, bubblelegend
%
% References:
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('balloonplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% balloonplot with table input.
    % Load the Housetasks data (a contingency table containing the
    % frequency of execution of 13 house tasks in the couple).
    N=[156	14	2	4;
        124	20	5	4;
        77	11	7	13;
        82	36	15	7;
        53	11	1	57;
        32	24	4	53;
        33	23	9	55;
        12	46	23	15;
        10	51	75	3;
        13	13	21	66;
        8	1	53	77;
        0	3	160	2;
        0	1	6	153];
    rowslab={'Laundry' 'Main-meal' 'Dinner' 'Breakfeast' 'Tidying' 'Dishes' ...
        'Shopping' 'Official' 'Driving' 'Finances' 'Insurance'...
        'Repairs' 'Holidays'};
    colslab={'Wife'	'Alternating'	'Husband'	'Jointly'};
    tableN=array2table(N,'VariableNames',colslab,'RowNames',rowslab);
    % In this case a table is supplied
    balloonplot(tableN);
%}

%{
    % balloonplot with array input.
    % Load the Housetasks data (a contingency table containing the
    % frequency of execution of 13 house tasks in the couple).
    N=[156	14	2	4;
        124	20	5	4;
        77	11	7	13;
        82	36	15	7;
        53	11	1	57;
        32	24	4	53;
        33	23	9	55;
        12	46	23	15;
        10	51	75	3;
        13	13	21	66;
        8	1	53	77;
        0	3	160	2;
        0	1	6	153];
    rowslab={'Laundry' 'Main-meal' 'Dinner' 'Breakfast' 'Tidying' 'Dishes' ...
        'Shopping' 'Official' 'Driving' 'Finances' 'Insurance'...
        'Repairs' 'Holidays'};
    colslab={'Wife'	'Alternating'	'Husband'	'Jointly'};
    tableN=array2table(N,'VariableNames',colslab,'RowNames',rowslab);
    % In this case row and columns labels are supplied through options Lr
    % and Lc
    balloonplot(N,'Lr',rowslab,'Lc',colslab);
%}

%{
    %% balloonplot with original data matrix as input.
    load smoke
    [h,Ntable]=balloonplot(smoke,'datamatrix',true);
    disp(Ntable)
%}

%% Beginning of code

% Check whether N is a contingency table or a n-by-p input dataset (in this
% last case the contingency table is built using the first two columns of the
% input dataset).
if ~isempty(varargin)
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    checkdatamatrix = strcmp(UserOptions,'datamatrix');
    if sum(checkdatamatrix)
        datamatrix = varargin{2*find(checkdatamatrix)};
    else
        datamatrix=false;
    end
else
    datamatrix=false;
end

% If input is a datamatrix it is necessary to construct the contingency
% table
if datamatrix == true
    if size(N,2)>2
        warning('FSDA:CorAna:TooManyVars','Input array or table has more than 2 columns. CorAna uses the first two');
    end
    if istable(N)
        [N,~,~,labelsTab] =crosstab(N{:,1},N{:,2});
    else
        [N,~,~,labelsTab] =crosstab(N(:,1),N(:,2));
    end
    [I,J]=size(N);
    % default labels for rows of contingency table
    Lr=labelsTab(1:I,1);
    % default labels for columns of contingency table
    Lc=labelsTab(1:J,2);
    % Make valid names
    Lr=matlab.lang.makeValidName(Lr);
    Lc=matlab.lang.makeValidName(Lc);
else
    [I,J]=size(N);
    % Size of N
    % default labels for rows of contingency table
    Lr=cellstr(strcat('r',num2str((1:I)')));
    % default labels for columns of contingency table
    Lc=cellstr(strcat('c',num2str((1:J)')));
end
ax='';
options=struct('Lr',{Lr},'Lc',{Lc},'datamatrix',false,...
    'ax',ax);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:balloonplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    Lr  = options.Lr;
    Lc  = options.Lc;
    ax=options.ax;
    Lr=matlab.lang.makeValidName(Lr);
    Lc=matlab.lang.makeValidName(Lc);
    
end

% Extract labels for rows and columns
if istable(N)
    Lc=N.Properties.VariableNames;
    Lr=N.Properties.RowNames;
    Ntable=N;
    N=table2array(N);
else
    if isempty(Lr)
        Lr=cellstr(num2str((1:I)'));
    else
        % Check that the length of Lr is equal to I
        if length(Lr)~=I
            error('FSDA:balloonplot:WrongInputOpt','Wrong length of row labels');
        end
    end
    
    if isempty(Lc)
        Lc=cellstr(num2str((1:J)'));
    else
        % Check that the length of Lc is equal to J
        if length(Lc)~=J
            error('FSDA:balloonplot:WrongInputOpt','Wrong length of column labels');
        end
    end
    Ntable=array2table(N,'RowNames',Lr,'VariableNames',Lc);
end


[I,J]=size(N);
ycoo=meshgrid((I:-1:1),J:-1:1)';
ycoo=ycoo(:);
xcoo=repmat(1:J,I,1);
xcoo=xcoo(:);
Nvector=N(:);
if isempty(ax)
    h=bubblechart(xcoo,ycoo,Nvector,Nvector);
else
    h=bubblechart(ax,xcoo,ycoo,Nvector,Nvector);
end
axes1=gca;
jall=1:J;
% If there are more than 70 columns just show a sistematic sample of
% (approximately) 70 of them
if J<=70
    set(axes1,'XTick',jall,'XTickLabel',Lc,'TickLabelInterpreter','none');
else
    step=ceil(J/70);
    sel=1:step:J;
    set(axes1,'XTick',jall(sel),'XTickLabel',Lc(sel),'TickLabelInterpreter','none');
end
set(axes1,'YTick',1:I,'YTickLabel',flip(Lr),'TickLabelInterpreter','none');
bubblesize([3 20])
colorbar(axes1)
grid('on')

end
%FScategory:VIS-Mult