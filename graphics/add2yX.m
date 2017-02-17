function add2yX(H,AX,BigAx,varargin)
%add2yX adds objects to the yXplot.
%
%
%
% Required input arguments:
%
%            H   : handles to the lines on the graphs. 2D or 3D array. H is an array of handles to the
%                   lines on the graphs. The array's third dimension
%                   corresponds to groups in the grouping variable.
%                   For more details see the output arguments of gplotmatrix.
%           AX   : handles to the axes of the individual plots. Matrix.  AX is a
%                   matrix of handles to the axes of the individual plots.
%                   For more details see the output arguments of gplotmatrix.
%          BigAx :  handle to big (invisible) axes framing the
%                   entire plot matrix. Scalar.
%                   For more details see the output arguments of gplotmatrix.
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
% More About:
%
%   add2yX essentially used within FSDA function yXplot.m. However its
%   logic can be also demonstrated with MATLAB function gplotmatrix.m, as in
%   the examples below.
%
%   As default add2yX makes legends in the existing yX plot matrix
%   clickable or creates a clickable multilegend if the legend does not
%   exist.
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
% References:
%
%   Tufte E.R. (1983). The visual display of quantitative information. Graphics Press, Cheshire
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('add2spm')">Link to the help function</a>
% Last modified 14-06-2016

% Examples:

%
%{
    %% add2yX with all default options.
    n=100;
    p=3;
    X=randn(n,p);
    y=100+randn(n,1);
    sel=51:100;
    y(sel)=y(sel)+2;
    group=ones(n,1);
    group(sel)=2;
    [H,AX,BigAx] = gplotmatrix(X,y,group);
    % The legengs are not clickable
    add2yX(H,AX,BigAx)
    % Now the legends become clickable
%}

%{
    %% Compare gplotmatrix + addyX with yXplot for IRIS data.
    load fisheriris;
    % Create scatter plot matrix with specific legends
    % plot Sepal length (y) againt the other variables
    y=meas(:,1);
    X=meas(:,2:4);
    [H,AX,BigAx]=gplotmatrix(X,y,species,[],[],[],'on');
    % The legends are not clickable
    add2yX(H,AX,BigAx)
    % Now the legends become clickable
    % It is easier to call directly function yXplot
    [H,AX,BigAx]=yXplot(meas(:,1),meas(:,2:4),species)
%}

%{
    %% Example of use of option bivarfit.
    load fisheriris;
    % Create scatter plot matrix with specific legends
    % plot Sepal length (y) againt the other variables
    y=meas(:,1);
    X=meas(:,2:4);
    [H,AX,BigAx]=yXplot(meas(:,1),meas(:,2:4),species)
    % add a regression line to each group
    add2yX(H,AX,BigAx,'bivarfit','0')
%}
    
%% Beginning of code

% FontSizelabeladd= height of text labels
FontSizelabeladd=12;

intercept=1;
bivarfit='';
multivarfit='';
labeladd='';
userleg='';

%get optional user options
if nargin>3
    % user options
    options= struct('intercept',1,'bivarfit','','multivarfit','','labeladd','','userleg','');

    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:add2yX:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    intercept   = options.intercept;
bivarfit    = options.bivarfit;
multivarfit = options.multivarfit;
labeladd    = options.labeladd;
userleg  = options.userleg;

end


% To take account a change in property names of the legend object in 2016b
if verLessThan('matlab','9.1')
    legstring='LegendPeerHandle';
else
    legstring='LayoutPeers';
end



%% get the values of the observations from the yXplot

% get its handle
fig = ancestor(BigAx,'figure');

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

% AX= findobj(hChildren,'Type','axes','-and','Tag','');

% force to build default legends if there are no legends in the scatterplot
if ~isappdata(AX(1,end),legstring)
    userleg = '1';
end

% These are the legends already in the plot
legplot = get(getappdata(AX(1,end),legstring),'String');

% if 'userleg' is empty, use the legend already in the plot.
if isempty(userleg)
    legnew = legplot;
end

% if 'userleg' is a cell of string, use such strings as user-defined legends
if ~isempty(userleg) && iscell(userleg)
    legnew = userleg;
end


nAX = numel(AX);

% For bivariate data no need of multivarfit
if nAX ==1
    if isempty(bivarfit) && ~isempty(multivarfit)
        bivarfit = multivarfit;
    else
        multivarfit = '';
    end
end

% hPlotMatrixAxC: the children of AX are in a cell of nAX elements (one
% for each panel of the yXplot). Each element of the cell is a vector
% containing the handles of the groups within each panel.
hPlotMatrixAxC = get(AX,'Children');
if iscell(hPlotMatrixAxC)
    ngroups = size(hPlotMatrixAxC{1,1},1);
else
    ngroups = size(hPlotMatrixAxC,1);
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



if isempty(userleg) || (~isempty(userleg) && iscell(userleg))
    nleg = numel(legnew);
    
    for i = 1 : nleg
        HH = H(:,:,i) ;
        set(HH(:),'DisplayName',legnew{i});
    end
    set(gcf,'Name','Scatter plot matrix with groups highlighted');
end


% if 'userleg' is '1', set context sensitive group-specific legends.
% The context is determined by the occcurence of specific words in the Tag
% of the current figure. The currently addressed strings/contexts are
% 'outlier' (for outliers/normal units), 'brush' (for Brushed units 1,
% Brushed units 2, etc.) and 'group' (for 'Group 1, Group 2, etc.).
if ~isempty(userleg) && ischar(userleg) && strcmp(userleg,'1')
    
    % add multilegend
    v = size(AX,2);
    leg = get(getappdata(AX(1,end),legstring),'String');
    nleg = numel(leg);
    
    if ndims(H) == 3
        % The third dimension of H is to distinguish the groups. In the next
        % 'if' statement we use two equivalent ways to deal with H, considering
        % that the diagonal of the scatter matrix is dedicated to the
        % histograms.
        if nleg == 2 && ~isempty(strfind(lower(get(gcf, 'Tag')),'outlier')) %#ok<*STREMP>
            set(H(H(:,:,2)~=0),'DisplayName','Normal units');
            linind      = sub2ind([v v],1:v,1:v);
            outofdiag   = setdiff(1:v^2,linind);
            lin2ind     = outofdiag+v^2;
            set(H(lin2ind),'DisplayName','Outliers');
        else
            % Assign to this figure a name
            set(gcf,'Name','yXplot with groups highlighted');
            % Reset the handles of the main diagonal (histograms) to zero.
            
            
            % Now reshape the handles array to make it more manageable: while H
            % is a 3-dimensional array with the third dimension associated to
            % the groups, newH is 2-dimensional with columns associated to the
            % lines of the scatterplot and lines associated to the groups.
            newH = reshape(H,numel(H)/nleg,nleg);
            if ~isempty(strfind(lower(get(gcf, 'Tag')),'brush'))
                %            if strcmp(get(gcf, 'Tag'),'pl_yX') || ~isempty(strfind(lower(get(gcf, 'Tag')),'brush'))
                % set the legend of the unbrushed units
                set(newH(newH(:,1)~=0),'DisplayName','Unbrushed units');
                % set the legend of the brushed units
                for i = 2 : nleg
                    set(newH(newH(:,1)~=0,i),'DisplayName',['Brushed units ' num2str(i-1)]);
                end
            elseif ~isempty(strfind(lower(get(gcf, 'Tag')),'group'))
                for i = 1 : nleg
                    set(newH(newH(:,1)~=0,i),'DisplayName',['Group ' num2str(i)]);
                end
            else
                % here the tag is empty: in this case take the legends
                % provided by the user
                for i = 1 : nleg
                    %leguser = get(getappdata(AX(1,end),'LegendPeerHandle'),'String');
                    set(newH(newH(:,1)~=0,i),'DisplayName',legplot{i});
                end
            end
        end
    else
        % In this case there are no groups in the data
        set(setdiff(H(:),diag(H)),'DisplayName','Units')
    end
    
    % Get the final legends
    legnew = get(getappdata(AX(1,end),legstring),'String');
end

% Now update the legends in the plot and make them clickable.
hLines  = findobj(AX(1,end), 'type', 'line');
%  [~,sorindlegnew]=sort({hLines.DisplayName});


if ~isempty(legnew)
    % The loop which follows has been done in order to guarantee a one to
    % one correspondence between the order of the names of the groups in
    % legnew and the names contained in {hLines.DisplayName}
    % For example it may happen that
    % legnew = {'Unbrushed units'    'Brushed units 1'} and
    % {hLines.DisplayName} is  {'Brushed units 1'    'Unbrushed units'}
    hNames={hLines.DisplayName};
    sorindlegnew=zeros(nleg,1);
    seqnames=1:nleg;
    for j=1:nleg
        sorindlegnew(j)=seqnames(strcmp(hNames(j),legnew));
    end
    % make the legend entries clickable
    clickableMultiLegend(hLines(sorindlegnew), legnew{:});
    
end


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
            if intercept==1
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

%
% legnew = get(getappdata(AX(1),legstring),'String');
% % Now update the legends in the plot and make them clickable.
% hLines  = findobj(AX(1), 'type', 'line');
% if ~isempty(legnew)
%     clickableMultiLegend(sort(double(hLines)), legnew{:});
% end

%% Update the legends with the new objects and make them clickable
% hLines = findobj(AX(end), 'type', 'line');
% hLines = findobj(AX(1), 'type', 'line');
%
% eLegend=get(hLines, 'DisplayName');
% % Create clickable multilegend if there is more than one group in the spm
% if iscell(eLegend)
%     hleg = findobj(fig,'Tag','legend');
%     set(hleg,'Visible','off');
%     clickableMultiLegend(hLines, eLegend{:},'Location','northwest');
% elseif ischar(eLegend)
%     hleg = findobj(fig,'Tag','legend');
%      set(hleg,'Visible','off');
%     clickableMultiLegend(hLines, eLegend,'Location','northwest');
% else
%
% end
end
