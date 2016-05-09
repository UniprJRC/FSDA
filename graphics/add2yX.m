function add2yX(varargin)
%add2yX adds objects to the yXplot.
%
%<a href="matlab: docsearchFS('add2yx')">Link to the help page for this function</a>
%
%
% Required input arguments:
%
%
% Optional input arguments:
%
%           bivarfit:   Add a line fit. Char.
%                       '0' fit a line to each group; 
%                       '1' fit 1 line to all data, regardless the groups; 
%                       '2' fit a line on all data and a line to relevant
%                       data; 
%                       ''  the default, nothing is added.
%                       Example - 'bivarfit','1'
%                       Data Types - char
%
%           multivarfit: Add a multivariate fit. Char.
%                       '1' one multivariate fit on all units
%                       '2' one multivariate fit on all units and one on relevant data
%                       ''  the default, nothing is added
%                       Example - 'multivarfit','1'
%                       Data Types - char
%
%           labeladd:   Add labels. Char.
%                       '1' add labels to relevant units
%                       ''  the default, nothing is added
%                       Example - 'labeladd','1'
%                       Data Types - char
%
%           intercept:  Indicator for constant term. Scalar. 
%                       intercept = 1 (default) assumes the intercept for the
%                       bivarfit and multivarfit. When intercept = 0, 
%                       the intercept is not used for the fits.
%                       Example - 'intercept',1 
%                       Data Types - double
%
% Output:
%
%
% More About:
%
% Note that the function extracts the data from the graphical
% objects in the plot. At the current stage the objects that can be 
% added to yXplot using add2yX are:
% - a bivariate fit on each panel of yXplot (see olsline.m); 
% - a multivariate fit that is shown on each panel of yXplt; 
% - the labels of relevant observations, e.g. outliers or brushed groups.
%
%
% See also olsline
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('add2yX')">Link to the help function</a>
% Last modified 06-Feb-2015


%% Beginning of code

% user options
options= struct('intercept',1,'bivarfit','','multivarfit','','labeladd','');

% FontSizelabeladd= height of text labels
FontSizelabeladd=12;

%get optional user options
if nargin>0
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:add2yX:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

intercept   = options.intercept;
bivarfit    = options.bivarfit;
multivarfit = options.multivarfit;
labeladd    = options.labeladd;

%% get the values of the observations from the yXplot

% assume the current figure is yXplot and get its handle
fig = gcf;

% the children of the yXplot
hChildren  = get(fig,'Children');

% BigAx: the handle to big (invisible) axes framing the entire plot matrix.
% remark: fig = ancestor(BigAx,'figure');
% BigAx = hChildren(nhChildren);

% AX: the handle(s) of the axes of the panels in the yXplot

% In other words, we have to extract from hchildren just the handles whose
% Type is axes but not bigaxes (therefore Tag must be empty
% In Matlab <2014b these tags where all but the first and the last (the
% first was the legend and the last was bigaxes). From Matlab 2014b the
% order has changed so it is necessary to use instruction findobj in order
% to extract just the subhandles axes (excluding PlotMatrixBigAx)
% AX  = hChildren;
% nhChildren = numel(hChildren);
% AX([1 nhChildren])=[];

AX= findobj(hChildren,'Type','axes','-and','Tag','');

nAX = numel(AX);

% For bivariate data no need of multivarfit
if nAX ==1
    if isempty(bivarfit) && ~isempty(multivarfit)
        bivarfit = multivarfit;
    else multivarfit = '';
    end
end

% hPlotMatrixAxC: the children of AX are in a cell of nAX elements (one
% for each panel of the yXplot). Each element of the cell is a vector
% containing the handles of the groups within each panel.
hPlotMatrixAxC = get(AX,'Children');
if iscell(hPlotMatrixAxC)
    ngroups = size(hPlotMatrixAxC{1,1},1);
else ngroups = size(hPlotMatrixAxC,1);
end

% H: the handles of the groups within each panel are rearranged in a
%    three-dimensional matrix H having the structure of that produced by
%    gplotmatrix. Then, matrix H is used to extract the values of the
%    yXplot data points.
% Remark 1: the points in each group are MATLAB line objects.
% Remark 2: the first group of H, i.e. H(1,i,1)  for each panel i, is the
% group of the unselected units (Xigood ygood).
% Remark 3: the last group of H, i.e. H(1,i,end) for each panel i, is the
% group of the last selected units (Xilast ylast).
% Remark 4: in yXplot the call to gplotmatrix is so that the unselected
% units are always identified by the "group" parameter with value 1.
% Remark 5: The order of the observations as passed to gplotmatrix is not
% the same of the Xi extracted here. As a consequence, the observations
% identified by "group=1" may not be equivalent to (Xigood ygood).

verMatlab=verLessThan('matlab','8.4.0');
if verMatlab
    H = NaN(1,nAX,ngroups);
else
    H=gobjects(1,nAX,ngroups);
end

if nAX > 1
    for i=1:nAX
        H(1,i,:) = fliplr(hPlotMatrixAxC{i,:}');
        if i==1
            if ngroups == 1
                y = get(H(1,i,:),'YData')';
            else
                y = cell2mat(get(H(1,i,:),'YData')')';
                ygood  = get(H(1,i,1),'YData')';
                ylast  = get(H(1,i,end),'YData')';
                Xi     = zeros(numel(y),nAX);
                Xigood = zeros(numel(ygood),nAX);
                Xilast = zeros(numel(ylast),nAX);
            end
        end
        if ngroups == 1
            Xi(:,i)     = get(H(1,i,:),'XData')';
        else
            Xi(:,i)     = cell2mat(get(H(1,i,:),'XData')')';
            Xigood(:,i) = get(H(1,i,1),'XData')';
            Xilast(:,i) = get(H(1,i,end),'XData')';
        end
    end
elseif nAX==1
    H(1,1,:) = fliplr(hPlotMatrixAxC');
    if ngroups == 1
        y        = get(H(1,1,:),'YData')';
        Xi       = get(H(1,1,:),'XData')';
    else
        y        = cell2mat(get(H(1,1,:),'YData')')';
        ygood    = get(H(1,1,1),'YData')';
        ylast    = get(H(1,1,end),'YData')';
        Xi       = cell2mat(get(H(1,1,:),'XData')')';
        Xigood   = get(H(1,1,1),'XData')';
        Xilast   = get(H(1,1,end),'XData')';
    end
    %     y  = get(hPlotMatrixAxC(numel(hPlotMatrixAxC)),'YData')';
    %     Xi = get(hPlotMatrixAxC(numel(hPlotMatrixAxC)),'XData')';
else
    disp('Error: nAX cannot be negative or zero.');
end



% Get the labels of the last selected group of units.
nbrush = get(H(:,1,end), 'UserData');

% Get the legenda of the last selected group of units.
%DisplayNameLast = get(H(1,1,end), 'DisplayName');
DisplayNameFirst = get(H(1,1,1), 'DisplayName');

if intercept == 1
    X = cat(2,ones(numel(y),1),Xi);
    if ngroups > 1
        Xgood = cat(2,ones(numel(ygood),1),Xigood);
    end
else
    X = Xi;
    if ngroups > 1
        Xgood = Xigood;
    end
end


[~,p]=size(X);

intcolumn = find(max(X,[],1)-min(X,[],1) == 0);

%% Add the objects

% We need to add objects to the scatterplots of y|X
set(fig,'NextPlot','add');
set(AX,'NextPlot','add');

% Now, for each scatter of y|X do:
for i = 1:length(AX)
    % Make the axes of the panel with handle AX(i) the current axes.
    set(fig,'CurrentAxes',AX(i));
    % Remark: axes(AX(i)) would also do the job, but would
    % restack the axes above all other axes in the figure.
    
    % Fit least square line(s) to the scatterplot AX(i).
    switch bivarfit
        case ''
            %do nothing: no line is fit.
        case '0'
            h=olsline(0,intercept); % fit a line to each group
            if length(h)==1
                set(h,'DisplayName',['bivarfit on ' DisplayNameFirst] ); %unbrushed units of y|Xi
            else
                for brugrp = 1:size(h,1)-1
                    set(h(brugrp),'DisplayName',['bivarfit on ' get(H(1,1,size(h,1)-brugrp+1), 'DisplayName')]); % brushed units num2str(size(h,1)-brugrp)
                end
                set(h(size(h,1)),'DisplayName',['bivarfit on ' DisplayNameFirst]);%unbrushed units of y|Xi
            end
        case '1'
            h=olsline(-1,intercept);    % fit 1 line to all data, regardless the groups
            set(h,'DisplayName','bivarfit on all units');
        case '2'
            h=olsline(-1,intercept);    % fit a line to all data
            set(h,'DisplayName','bivarfit on all units');
            if ngroups > 1
                h=olsline(size(H,3),intercept); % fit a line to the unselected data, i.e. the last group among the handles H
                set(h,'DisplayName',['bivarfit on ' DisplayNameFirst]);%'fit on unbrushed units of y|Xi'
            end
        otherwise
            if strncmp('i',bivarfit,1)
                token = strtok(bivarfit, 'i');
                selgroup=str2double(token);
                if ~isnan(selgroup) && selgroup > 0.5 && selgroup <= size(H,3)
                    h=olsline(size(H,3)-round(selgroup)+1 , intercept); % fit one group only: the one with index round(bivarfit)
                    set(h,'DisplayName','bivarfit on a group of y|Xi');
                else
                    %do nothing: no line is fit.
                end
            else
                error('FSDA:add2yX:WrongBivarifit','Valid values for option ''bivarfit'' are: '''', ''0'', ''1'', ''2'', ''i1'', ''i2'', ... , ''ig'', ... being ''g'' the index of a selected group.')
                %do nothing: no line is fit
            end
    end
    
    % Add the labels for the last selected group.
    if strcmp('1',labeladd) && ngroups > 1
        xlimits = get(AX(i),'Xlim'); ylimits = get(AX(i),'Ylim');
        dx = (xlimits(2)-xlimits(1))*0.01*length(AX); dy = (ylimits(2)-ylimits(1))*0.01*length(AX)/2; % displacement
        %       text(Xi(nbrush,i)+dx,y(nbrush)+dy,numtext(nbrush),'HorizontalAlignment', 'Left');
        text(Xilast(:,i) + dx,ylast + dy,cellstr(num2str(nbrush,'%d')),'HorizontalAlignment', 'Left','FontSize',FontSizelabeladd);
    end
    
    % Add to each plot AX(i) the line(s) based on the hyperplane fit to y|X
    switch multivarfit
        case {'1' , '2'}
            coef=regress(y,X);
            % indcoef = vector which contains the indexes of columns of X
            % except that which is about to be plotted
            if intercept==1;
                indcoef = setdiff(1:p,i+1);
            else
                indcoef = setdiff(1:p,i);
            end
            % The mean value of the columns of matrix X using all the units
            meaot=mean(X(:,indcoef))*coef(indcoef);
            
            coef(intcolumn)=[]; % the other coefficients
            xlimits = get(AX(i),'Xlim');
            hline1 = line(xlimits , meaot + coef(i).*xlimits);
            datacolor = [.3 .3 .3]; %Dark grey;  [1 .62 .40] == %Copper
            set(hline1,'Color',datacolor,'LineWidth',2,...
                'DisplayName','multivarfit on all units');
            if strcmp(multivarfit, '2') && ngroups > 1
                coef = regress(ygood,Xgood);
                meaot=mean(Xgood(:,indcoef))*coef(indcoef); %Mean of all other variables considering just good units
                coef(intcolumn)=[];
                hline2 = line(xlimits , meaot  + coef(i).*xlimits);
                
                %'fit on unbrushed units of y|X'
                set(hline2,'Color',get(H(1,i,1),'Color'),'LineWidth',2,...
                    'DisplayName',['multivarfit on ' DisplayNameFirst]);
            end
        otherwise
            %do nothing
    end
end

%% Update the legends with the new objects and make them clickable
hLines = findobj(AX(end), 'type', 'line');
eLegend=get(hLines, 'DisplayName');
% Create clickable multilegend if there is more than one group in the spm
if iscell(eLegend)
    clickableMultiLegend(hLines, eLegend{:});
end
end
%FScategory:VIS-Reg