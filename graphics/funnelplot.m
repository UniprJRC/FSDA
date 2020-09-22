function h  = funnelplot(x, varargin)
%carbikeplot produces the carbike plot to find best relevant clustering solutions
%
%<a href="matlab: docsearchFS('funnelplot')">Link to the help function</a>
%
% Funnel charts show values across multiple stages in a process. For
% example, you could use a funnel chart to show the number of sales
% prospects at each stage in a sales pipeline. Typically, the values
% decrease gradually, allowing the bars to resemble a funnel.
%
%  Required input arguments:
%
%           x : input data. Vector or matrix.
%                If y is a vector, funnelplot plots one funnel plot. If y
%                is a matrix, funnelplot plots one funnel plot for each
%                column of y.
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
%  Output:
%
%         h:   graphics handle to the plot. Graphics handle. Graphics
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
%<a href="matlab: docsearchFS('funnelplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% funnelplot with all default options.
    x=[500 425 200 150 100 90];
    funnelplot(x)
%}

%{
    %% funnelplot with Labels option.
    x=[500 425 200 150 100 90];
    labels={'Prospects', 'Qualified prospects', 'Needs analysis', 'Price quotes', ...
        'Negotiations', 'Closed sales'};
    funnelplot(x,'Labels',labels)
%}



%% Beginning of code

if ~isnumeric(x)
    error('FSDA:carbikeplot:WrongInput','First input argument must be numeric.');
end

n=length(x);
Labels=cellstr(num2str((1:n)'));

    options=struct('Labels',{Labels});
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:funnelplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:funnelplot:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
 if nargin>1   
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    Labels=options.Labels;
    
end


Labels=flipud(Labels(:));

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);

kk=1.05;
 axis([-max(x)*kk  max(x)*kk 1.5 n+0.5]);
grid off;
set(axes1,'YTick',1.5:(n+0.5),'XTick','');
ylim([0.5 n+1.3])
for i=1:n
    rectangle('position',[-x(i) n+1-i 2*x(i)  1],'facecolor','c','Curvature',[0.02 0.02]);
    text(0, n+1.5-i ,num2str(x(i)),'HorizontalAlignment','center','FontSize',15,'VerticalAlignment','middle');
end

set(axes1,'YTickLabel',Labels);

box('on');
h=gcf;
end
%FScategory:VIS-Clu
