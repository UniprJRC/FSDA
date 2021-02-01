function h = waterfallchart(Y,varargin)
%waterfallchart creates a waterfall chart
%
%<a href="matlab: docsearchFS('waterfallchart')">Link to the help function</a>
%
% A waterfall chart is a form of data visualization that helps in
% understanding the cumulative effect of sequentially introduced positive
% or negative values. These intermediate values can either be time based or
% category based. The waterfall chart is also known as a flying bricks
% chart or Mario chart due to the apparent suspension of columns (bricks)
% in mid-air. Often in finance, it will be referred to as a bridge.
% Waterfall chart. For more details see
% https://en.wikipedia.org/wiki/Waterfall_chart
%
%  Required input arguments:
%
%            Y: Input data. Vector of table with one column.
%               Data matrix containing n observations
%                 Data Types -  double
%
%  Optional input arguments:
%
%   SetAsTotal  : elements which which have to be set as total.
%                  Logical or numeric index vector.
%                  Logical or numeric index vector indicating which
%                  elements of Y have to be set as total.
%               Example - 'SetAsTotal',[2,3] | logical([0 1 1 0 0 0])
%               Data Types - double | logical
%
%   BarWidth    : width of the bars. Scalar. A number in the interval 0 1
%                 which specifies the width of the bars. The deafult value is 0.6.
%                 For example if bardwidth is 0.5, than the first bar has
%                 coordinated 0.75 1.25 the second 1.75 2.25 ...
%               Example - 'barwidth',0.7
%               Data Types - double
%
%    titl     : plot title. character or cell. string containing plot title
%                    If Y is a table, the variableName of the table is used
%                    as title
%               Example - 'titl',{'demographic movements'}
%               Data Types - cell
%
%    Labels     : x labels. cell. cell of length n containing the labels
%                   for the n elements to add to the xTickLabels of the
%                   plot.
%                   If Y is a table, the rowsnames
%                   of the table are used as labels. The default is to use
%                   the labels 1, 2, ..., n
%               Example - 'labels',{'aa' 'bb' 'cc' ddf'}
%               Data Types - cell
%
% DisplayValueOnTopOfPatches : display values on top of patches. Boolean.
%               If this option is set to true, than the Y values are
%               displayed on top of the patches. The default is false, that
%               is values on top of the patches are not shown.
%               Example - 'DisplayValueOnTopOfPatches',true
%               Data Types - boolean
%
% ShowConnectorLines : Show connector lines. Boolean.
%               If this option is set to true (default), connector lines are shown in
%               the plots to connect the patches.
%               Example - 'ShowConnectorLines',false
%               Data Types - boolean
%
%  Output:
%
%      h :    handles to the patches. Vector of graphic handles.
%               3-by-1 vector of graphic handles.
%               h(1) is the handle to the patches which are SetAsTotal.
%               h(2) is the handle to the patches which have positive values.
%               h(3) is the handle to the patches which have positive values.
%
%
% See also funnelchart
%
% References:
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('waterfallchart')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % waterfall chart with all default options.
    X=[245631 -2412 243219 -114899 -18731 -6244 103345 -26745 ...
    -11279 -36000 29321 -4400 24921];
    waterfallchart(X);
%}

%{
    %% waterfall chart with options SetAsTotal and Labels.
    X=[245631 -2412 243219 -114899 -18731 -6244 103345 -26745 ...
    -11279 -36000 29321 -4400 24921];
    % Define the elements of X which contain the 
    SetAsTotal=[1 3 7 11 13];
    Labels={'Gross revenues' 'Ordinary expenses' ... 
        'Net revenues' 'Initial inventories' 'Merchandising', 'Other selling costs', ...
        'Gross Income' 'Wages' 'Marketing advertising' 'Insurance bank charges' ...
        'Operating income' 'Taxes' 'Net income'};
    waterfallchart(X,'SetAsTotal',SetAsTotal,'Labels',Labels);
%}

%{
    % waterfall chart with options ShowConnectorLines and BarWidth. 
    X=[245631 -2412 243219 -114899 -18731 -6244 103345 -26745 ...
    -11279 -36000 29321 -4400 24921];
    % Define the elements of X which contain the 
    SetAsTotal=[1 3 7 11 13];
    Labels={'Gross revenues' 'Ordinary expenses' ... 
        'Net revenues' 'Initial inventories' 'Merchandising', 'Other selling costs', ...
        'Gross Income' 'Wages' 'Marketing advertising' 'Insurance bank charges' ...
        'Operating income' 'Taxes' 'Net income'};
    BarWidth=0.8;
    waterfallchart(X,'SetAsTotal',SetAsTotal,'Labels',Labels,'ShowConnectorLines',false,'BarWidth',BarWidth);
%}

%{
    %% waterfall chart with input as a table.
    figure
    X=[515 133 -65 583 159 -70 672 189 -100 761]';
    rownam={'Population 2015' 'Births2015' 'Deaths2015' ...
    'Population 2016' 'Births2016' 'Deaths2016' ...
    'Population 2017' 'Births2017' 'Deaths2017' ...
    'Population 2018'};
    sel=[1 4 7 10];
    Xtable=array2table(X,'RowNames',rownam','VariableNames',{'Detailed Analysis of the population growth'});
    waterfallchart(Xtable,'SetAsTotal',sel);
%}

%{
    %%  waterfall called with h ouptut.
    X=[515 133 -65 583 159 -70 672 189 -100 761]';
    rownam={'Population 2015' 'Births2015' 'Deaths2015' ...
    'Population 2016' 'Births2016' 'Deaths2016' ...
    'Population 2017' 'Births2017' 'Deaths2017' ...
    'Population 2018'};
    sel=[1 4 7 10];
    Xtable=array2table(X,'RowNames',rownam','VariableNames',{'Detailed Analysis of the population growth'});
    h=waterfallchart(Xtable,'SetAsTotal',sel);
    % Change the colors of the patches  which are set as total
    h(1).FaceColor='k';
    % Change the colors of the patches  which have positive values
    h(2).FaceColor='c';
    % Change the colors of the patches  which have negative values
    h(3).FaceColor='b';
%}

%{
    % waterfall chart with option DisplayValueOnTopOfPatches.
    X=[245631 -2412 243219 -114899 -18731 -6244 103345 -26745 ...
    -11279 -36000 29321 -4400 24921];
    waterfallchart(X,'DisplayValueOnTopOfPatches',true);
%}


%{
    % waterfall chart with option SetAsTotal passed as logical.
    X=[245631 -2412 243219 -114899 -18731 -6244 103345 -26745 ...
    -11279 -36000 29321 -4400 24921];
    SetAsTotal=[true false true false false false true false false false true false true];
    Labels={'Gross revenues' 'Ordinary expenses' ... 
        'Net revenues' 'Initial inventories' 'Merchandising', 'Other selling costs', ...
        'Gross Income' 'Wages' 'Marketing advertising' 'Insurance bank charges' ...
        'Operating income' 'Taxes' 'Net income'};
    waterfallchart(X,'SetAsTotal',SetAsTotal,'Labels',Labels);
%}

%% Beginning of code

ShowConnectorLines = true;
SetAsTotal = '';
Labels     = '';
titl       = '';
BarWidth   = 0.5;
DisplayValueOnTopOfPatches = false;

options=struct('ShowConnectorLines',ShowConnectorLines,...
    'SetAsTotal',SetAsTotal,'BarWidth',BarWidth,'Labels',Labels,...
    'titl',titl,'DisplayValueOnTopOfPatches',DisplayValueOnTopOfPatches);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:waterfallchart:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present in
    % structure options Remark: the nocheck option has already been dealt
    % by routine chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:waterfallchart:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the choices of the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    ShowConnectorLines  = options.ShowConnectorLines;
    SetAsTotal          = options.SetAsTotal;
    BarWidth            = options.BarWidth;
    Labels              = options.Labels;
    titl                = options.titl;
    DisplayValueOnTopOfPatches = options.DisplayValueOnTopOfPatches;
end

% transform SetAsTotal from logical to numbers from 1 to n which identify
% the elements which have to be set as total
if islogical(SetAsTotal)
    seq        = 1:length(SetAsTotal);
    SetAsTotal = seq(SetAsTotal);
end

if istable(Y)
    Labels = Y.Properties.RowNames;
    titl   = Y.Properties.VariableNames;
    Y      = Y{:,1};
else
    if ~isvector(Y)
        error('FSDA:waterfallchart:WrgInput','Y must be a vector or a table')
    end
end

% make sure that Y is a column vector
Y = Y(:)';
Ywithborders = [0 Y 0];
n = length(Ywithborders);

% colum j of x refers to the x coordinates of the j-th patch vertices
% the four elements of x respectively refer to
% x(1,j) = lower left x coordinate
% x(2,j) = upper left x coordinate
% x(3,j) = upper right x coordinate
% x(4,j) = lower right x coordinated
x = meshgrid(1:n-1, 1:4);

horizontalWidthBars = [-BarWidth/2; -BarWidth/2; BarWidth/2; BarWidth/2];
x = bsxfun(@plus, x, horizontalWidthBars);

ycoordlower = 1:n-1;
yoordupper  = 2:n;
% Build 4-by-length(Y)+1 matrix y which will contain the vertical
% coordinates of the patches
y = Ywithborders([ycoordlower;yoordupper;yoordupper;ycoordlower]);

% Make ssure that the elements of SetAsTotal are integers in 1:n-2
if max(SetAsTotal)>n-2
    error('FSDA:waterfallchart:WrgInput',['max(SetasTotal)must be ' num2str(n-2)])
elseif min(SetAsTotal)<1
    error('FSDA:waterfallchart:WrgInput','min(SetasTotal)must be 1')
else
end

% Set the upper coordinated of the columns just after SetAsTotal
y(2:3,SetAsTotal+1) = y([1 4],SetAsTotal+1) + y(2:3,SetAsTotal+1);

% For subtotals the lower y coordinate are 0
y([1 4],SetAsTotal) = 0;

notsel = setdiff(1:size(x,2),[SetAsTotal SetAsTotal+1]);
% Set the coordinates of the columns defined by notsel
for j=1:length(notsel)
    jj=notsel(j);
    if jj>1
        y([1 4 ],jj)=y(2:3,jj-1);
        y(2:3,jj)=y([1 4],jj)+y(2:3,jj);
    end
end

y(2,SetAsTotal) = Ywithborders(SetAsTotal+1);
y(3,SetAsTotal) = Ywithborders(SetAsTotal+1);

x=x(:,1:end-1);
y=y(:,1:end-1);

% Gray color for the bars which are set as totals Color (211/255) is light
% gray dark gray is 168/255 while silver is 192/255
%lightgray = (211/255)*ones(1,3); 
lightgrey = FSColors.lightgrey.RGB; 
h         = gobjects(3,1);

if ~isempty(SetAsTotal)
    h(1) = patch(x(:,SetAsTotal),y(:,SetAsTotal),lightgrey);
    if DisplayValueOnTopOfPatches == true
        xtips2  = h(1).XData(1,:);
        ytips2  = h(1).YData(2,:);
        labels2 = string(ytips2);
        text(xtips2,ytips2,labels2,'HorizontalAlignment','left',...
            'VerticalAlignment','bottom')
    end
end
notSetAsTotal = setdiff(1:size(x,2),SetAsTotal);
pos = Y(notSetAsTotal)>0;
if max(pos)>0
    % green color for positive values
    h(2)=patch(x(:,notSetAsTotal(pos)),y(:,notSetAsTotal(pos)),'g');
    if DisplayValueOnTopOfPatches == true
        xtips2 = h(2).XData(1,:);
        ytips2 = h(2).YData(2,:);
        labels2 = string(ytips2);
        text(xtips2,ytips2,labels2,'HorizontalAlignment','left',...
            'VerticalAlignment','bottom')
    end
end

if max(~pos)>0
    % red color for negative values
    h(3)=patch(x(:,notSetAsTotal(~pos)),y(:,notSetAsTotal(~pos)),'r');
    if DisplayValueOnTopOfPatches == true
        xtips2  = h(3).XData(1,:);
        ytips2  = h(3).YData(1,:);
        labels2 = string(ytips2);
        text(xtips2,ytips2,labels2,'HorizontalAlignment','left',...
            'VerticalAlignment','bottom')
    end
end

if ShowConnectorLines==true
    xcoo=[x(3,1:end-1); x(3,1:end-1)+BarWidth];
    line(xcoo,y(2:3,1:end-1),'color',lightgrey);
end
set(gca,'XTick',1:(n-2))

if ~isempty(Labels)
    nl = max(strlength(Labels(:)));
    if nl > 10
        rot = 45;
    else
        rot = 90;
    end
    set(gca,'XTickLabel',Labels,'XTickLabelRotation',rot)
end

title(titl,'FontSize',16);

allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'ActivePositionProperty','outerposition','FontSize',14);

end
%FScategory:VIS-Mult