function h  = funnelchart(x, varargin)
%funnelchart displays a funnel chart
%
%<a href="matlab: docsearchFS('funnelchart')">Link to the help function</a>
%
% Funnel charts provide a graphic representation of data values across
% multiple stages in a process. The chart displays progressively decreasing
% values in proportions amounting to 100 percent in total. For example, you
% could use a funnel chart to show the number of sales prospects at each
% stage in a sales pipeline. Typically, the values decrease gradually,
% allowing the bars to resemble a funnel. This type of chart can also be
% useful in identifying potential problem areas in an organization’s sales
% processes. A funnel chart is similar to a stacked percent bar chart. For
% more details see https://en.wikipedia.org/wiki/Funnel_chart
%
%  Required input arguments:
%
%           x : input data. Vector or matrix.
%                If x is a vector, a single funnel chart is displayed.
%                If x is a matrix, the function displays one funnel chart
%                for each column of x.
%          Data Types - double
%
%  Optional input arguments:
%
%
% Labels  :  Box labels. character array | string array | cell array |
%            numeric vector. Box labels, specified as the comma-separated
%            pair consisting of 'Labels' and a character array, string
%            array, cell array, or numeric vector containing the box label
%            names. Specify one label per x value.
%                 Example - 'Labels',false
%                 Data Types - char | string | cell | single | double
%
%
% Color  :   Color of the boxes. Character | RGB triplet vector. The color
%            specified by the user. 
%                 Example - 'Color',[0.12 0.6 0.15]
%                 Data Types - char | array
%
% Title  :   Title. Character array. The title of the Funnel Chart
%            specified by the user. 
%                 Example - 'Title','The Funnel Chart'
%                 Data Types - char
%
%   h     :   Target axes. If you do not specify the axes the plot function
%             uses the current axes.
%                 Example - h,ax, where ax=subplot(3,2,5);
%                 Data Types - axes object | graphics object
%
%  Output:
%
%         h:   Graphics handle to the plot. Graphics handle. Graphics
%               handle which is produced on the screen.
%
%
%
% See also: carbikeplot
%
% References:
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('funnelchart')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% funnelchart with all default options.
    x=[500 425 200 150 100 90];
    funnelchart(x)
%}

%{
    %% funnelchart with Labels option.
    x=[500 425 200 150 100 90];
    labels={'Prospects', 'Qualified prospects', 'Needs analysis', 'Price quotes', ...
        'Negotiations', 'Closed sales'};
    funnelchart(x,'Labels',labels)
%}

%{
    %% funnelchart when x is a matrix.
    x=100*abs(randn(10,4));
    x=sort(x,1,'descend');
    labels={'A', 'B', 'C', 'D', 'E', ...
        'F', 'G' 'H', 'I' 'J'};
    funnelchart(x,'Labels',labels)
%}

%{
    %% funnelchart with a non-default color.
    x=100*abs(randn(10,4));
    x=sort(x,1,'descend');
    labels={'A', 'B', 'C', 'D', 'E', ...
        'F', 'G' 'H', 'I' 'J'};
    funnelchart(x,'Labels',labels,'Color',FSColors.greysh.RGB)
%}

%{
    %% funnelchart with a title.
    x=100*abs(randn(10,4));
    x=sort(x,1,'descend');
    labels={'A', 'B', 'C', 'D', 'E', ...
        'F', 'G' 'H', 'I' 'J'};
    funnelchart(x,'Labels',labels,'Title','A Funnel Chart')
%}

%{
    % Example of the use of option h.
    % Create a subplot
    ax=subplot(3,2,5);
    % Create the funnelchart inside the axes specified by ax instead of in
    % the current axes (gca)
    funnelchart(rand(5,1),'h',ax)
%}

%% Beginning of code

if ~isnumeric(x)
    error('FSDA:funnelchart:WrongInput','First input argument must be numeric.');
end

x(x<0)=0;
if isvector(x)
    n=length(x);
    p=1;
    x=x(:); % Make sure that x is a column vector
else
    [n,p]=size(x);
end

% some settings
space    = 0.15;   % space between rectangles
fontsize = 14;    % font of numbers and axes

% default options
Labels   = cellstr(num2str((1:n)'));
BoxColor = 'c';
Title    = [];
h='';
options  = struct('Labels',{Labels},'Color',BoxColor,'Title',Title,'h',h);


% user options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:funnelchart:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk = isfield(options,UserOptions);
    WrongOptions = UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:funnelchart:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin>1
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end  
    Labels    = options.Labels;  
    BoxColor  = options.Color;
    Title     = options.Title;
    h=options.h;
end

Labels = flipud(Labels(:));

% Create figure
if isempty(h)
figure;
h=gca;
emptyh=true;
else
    emptyh=false;
end

set(h,'YTick','','XTick','');

switch p
    case 1
        nr=1;
        nc=1;
    case 2
        nr=2;
        nc=1;
    case {3,4}
        nr=2;
        nc=2;
    case {5,6}
        nr=3;
        nc=2;
    case {7,8,9}
        nr=3;
        nc=3;
    case {10,11,12}
        nr=3;
        nc=4;
    otherwise
        nr=4;
        nc=4;
end

% line below is to have more space for the subplots
setappdata(gcf, 'SubplotDefaultAxesLocation',[0.11,0.05,0.795,0.85]);
kk = 1.05;
for j=1:p
    if emptyh==true
    h=subplot(nr,nc,j);
    end
    % Create axes
    % axes1 = axes('Parent',figure1);
    Xj=x(:,j);
    axis([-max(Xj)*kk  max(Xj)*kk 1.5 n+0.5]);
    grid off;
    set(h,'YTick',1.5:(n+0.5),'XTick','');
    ylim([0.5 n+1.3]);
    for i=1:n
        if Xj(i)>0
            rectangle(h,'position',[-Xj(i) n+1-i 2*Xj(i)  1-space],...
                'facecolor',BoxColor,'Curvature',[0.02 0.02],...
                'AlignVertexCenters','on');
            text(h,0, n+1.5-i ,num2str(Xj(i),'%.1f'),'HorizontalAlignment','center',...
                'FontSize',fontsize,'VerticalAlignment','middle');
        end
    end
    set(h,'YTickLabel',Labels);
    box('on');
end

sgtitle(Title,'FontSize',fontsize+6);

h=gcf;
allAxesInFigure = findall(h,'type','axes');
set(allAxesInFigure,'FontSize',fontsize+2);
%ActivePositionProperty or  PositionConstraint
%set(allAxesInFigure,'ActivePositionProperty','outerposition','FontSize',fontsize+2);

end
%FScategory:VIS-Mult
