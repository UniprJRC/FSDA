function [varargout] = clickableMultiLegend(varargin)
%clickableMultiLegend hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend.
%
%<a href="matlab: docsearchFS('clickableMultiLegend')">Link to the help page for this function</a>
%
% It is typically applied to gplotmatrix figures. By clicking on a text
% label in the legend, the graphics (line or patch) objects associated to
% that label in all subplots are turned on and off (hide/show).
%
% The extention to multiple plots is realised by looking for graphics
% objects with the same DisplayName property of the one associated to the
% legend label. Therefore, the function should work also through plots in
% different figures.
%
% clickableMultiLegend accepts the same parameters of the legend function
% and can be used in the same way.
%
% Required input arguments:
%
% Optional input arguments:
%
% Output:
%
% Optional Output:
%
%     HLEG : handle to legend. Graphics handle. This is the handle to
%            legend on the current axes or empty if none exists.
%
% See also: legend, yXplot
%
% References:
%
% Deoras A. (2008),
% http://www.mathworks.com/matlabcentral/fileexchange/21799-clickablelegend-interactive-highlighting-of-data-in-figures/content/clickableLegend.m
% [clickableMultiLegend extends the clickableLegend by Ameya Deoras to
% figures with one legend for several subplots].
%
% Copyright 2008-2024.
% clickableMultiLegend has been adapted to this toolbox by FSDA team
%
%<a href="matlab: docsearchFS('clickableMultiLegend')">Link to the help page for this function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples

%{
   % ClickableMultilegend applied to a single scatter with several groups.
    % Simulate 3 groups
    n1=30; n2=20; n3=50;
    y = [rand(n1,1); rand(n2,1)+1; rand(n3,1)+2];
    group= [2*ones(n1,1); ones(n2,1); zeros(n3,1)];
    X = rand(n1+n2+n3,1);
    gscatter(X,y,group);

    % Make the legends clickable
    clickableMultiLegend();
%}

%{
     % clickableMultiLegend applied to two plots with lines with same legend names.
     % Clicking on a legend line, will show/hide the selected line in both plots. 
     z = peaks(100);
     plot(z(:,26:5:50),'LineWidth',3)
     grid on;

     % define the names of the legend labels for both plots
     names = {'Line1','Line2','Line3','Line4','Line5'};
     clickableMultiLegend(names , 'Location', 'NorthWest','FontSize',14);

     figure;
     plot(z(:,32:5:56),'LineWidth',3)
     grid on;
     hlegend=clickableMultiLegend(names, 'Location', 'NorthWest','FontSize',14);
%}

%{
    % clickableMultiLegend applied without arguments to multiple subplots.
    
    % Simulate X and y with 3 groups
    X = rand(100,4);
    y = [rand(10,1); rand(20,1)+1; rand(70,1)+2];
    group= [2*ones(10,1); ones(20,1); zeros(70,1)];

    % Generate the scatter matrix
    gplotmatrix(X,y,group);

    % Withut arguments, the effect is on the default legend
    clickableMultiLegend();

%}

%{
    % clickableMultiLegend applied with legend style arguments to multiple subplots.
    
    % Simulate X and y with 3 groups
    X = rand(100,4);
    y = [rand(10,1); rand(20,1)+1; rand(70,1)+2];
    group= [2*ones(10,1); ones(20,1); zeros(70,1)];

    % Generate the scatter matrix
    gplotmatrix(X,y,group);

    % Now we just want to change the font size of the legend
    clickableMultiLegend('FontSize',14);
%}

%{
    % clickableMultiLegend applied with legend names to multiple subplots.

    % Simulate X and y with 3 groups
    X = rand(100,4);
    y = [rand(10,1); rand(20,1)+1; rand(70,1)+2];
    group= [2*ones(10,1); ones(20,1); zeros(70,1)];

    % Generate the scatter matrix
    gplotmatrix(X,y,group);

    % Update the legend and make them clickable
    clickableMultiLegend({'group 1' , 'group 2' , 'group 3'},'FontSize',14);

%}

%{
    % clickableMultiLegend applied with the handles of the line objects in multiple subplots.
    
    % This example shows that clickableMultiLegend can also receive the 
    % handles returned by a multiple panel figure such as gplotmatrix.

    % Simulate X and y with 3 groups
    X = rand(100,4);
    y = [rand(10,1); rand(20,1)+1; rand(70,1)+2];
    group= [2*ones(10,1); ones(20,1); zeros(70,1)];

    % Generate the scatter matrix
    [H,AX,bigax] = gplotmatrix(X,y,group);

    % Set the DisplayName property (i.e. the texts of the legend) in all panels.
    % Note that in the gplotmatrix only one legend is visible.
    set(H(:,:,1),'DisplayName','group 10');
    set(H(:,:,2),'DisplayName','group 20');
    set(H(:,:,3),'DisplayName','group 30');

    % Get the handles of the legend to update
    hLines  = findobj(AX(1,end), 'type', 'line');

    % Update the legend and make them clickable
    clickableMultiLegend(hLines,'FontSize',14);

%}

%{
    % clickableMultiLegend again with the line objects handles, with a change in the legend order.

    % Simulate X and y with 3 groups
    X = rand(100,4);
    y = [rand(10,1); rand(20,1)+1; rand(70,1)+2];
    group= [2*ones(10,1); ones(20,1); zeros(70,1)];

    % Generate the scatter matrix
    [H,AX,bigax] = gplotmatrix(X,y,group);

    % Get the new legend texts directly from the plot
    % (takes into account a change in the legend property names in 2016b)
    if verLessThan('matlab','9.1')
        legstring='LegendPeerHandle';
    else
        legstring='LayoutPeers';
    end
    leg = get(getappdata(AX(1,end),legstring),'String');

    % Get the handles of the legend to update
    hLines  = findobj(AX(1,end), 'type', 'line');

    % Change the order of the legend strings, for exaple by sorting
    hLines = sort(hLines);

    % Update the legend and make them clickable
    clickableMultiLegend(hLines, leg{:},'FontSize',14);

%}

%{
    % clickableMultiLegend again with the line objects handles, with changes in the legend text.

    % Here we make a gplotmatrix plot clickable, then we change the labels
    % after a convenient reshape of its handles array: while H is a 3-dimensional
    % array with the third dimension associated to the groups, newH is 2-dimensional 
    % with lines associated to the subplots of the scatterplot and columns
    % associated to the groups. This simplifies the redefinition of the 
    % DisplayName property.

    % Simulate X and y with 3 groups
    X = rand(100,4);
    y = [rand(10,1); rand(20,1)+1; rand(70,1)+2];
    group= [2*ones(10,1); ones(20,1); zeros(70,1)];

    % Generate the scatter matrix
    [H,AX,bigax] = gplotmatrix(X,y,group);

    % make the legend clickable
    clickableMultiLegend('FontSize',14);

    % reshape the line handles array
    nleg = numel(hLines);
    newH = reshape(H,numel(H)/nleg,nleg);

    % redefine the legend texts
    for i = 1 : nleg
        set(newH(:,i),'DisplayName',['Redefined group n. ' num2str(i)]);
    end

    % If the legend texts were clickable before the re-definition, they
    % will remain clickable.
%}


%% Beginning of code
 
drawnow;

xlim manual;

% fix the DisplayName property in case of multiple panels before applying the callback function
if nargin > 0 

    % ensure to work on the axes handle to which a legend belongs
    hLeg = findobj(gcf,'tag','legend');
    if ~isempty(hLeg)
        axes(hLeg.Axes);
    end

    % 
    lgd = legend(varargin{:});
    fixDisplayName(gcf);

else

    leg = fixDisplayName(gcf);
    lgd = legend(leg);

end

% apply the callback for making the legend clickable
set(lgd, 'ItemHitFcn', @(src, event)togglevisibility(event.Peer));
varargout={lgd};

drawnow;

axis manual;

    function togglevisibility(plotHandle)

        % whichFigs contains the handle(s) of the target figures. When it
        % is set to groot, clickableMultiLegend switches on/off all figures
        % containing objects with the same DisplayName. If it is set to,
        % for example, parentFigs (the ancestor), then clickableMultiLegend
        % switches on/off only the active figure. 
        %parentFig = ancestor(plotHandle,'Figure');
        whichFigs  = groot;

        % Toggle the visibility of the plot handle (groot vs parentFig)
        allPlots = findall(whichFigs, 'DisplayName',plotHandle.DisplayName);
        for i = 1:length(allPlots)
            if (allPlots(i).Visible)
                allPlots(i).Visible = 'off';
            else
                allPlots(i).Visible = 'on';
            end
        end

        h1 = findall(groot, '-not','Type','Line','-not','Type','Axes','-and','Tag',plotHandle.DisplayName);
        if ~isempty(h1)
            ccur = get(h1(1),'FaceColor');
            %ccur = cell2mat(ccur(1,:));
            if ~isequal(ccur , [1 1 1])
                set(h1, 'UserData',get(h1,'FaceColor'));
                set(h1, 'FaceColor','w', 'EdgeColor','k');
            else
                cori = get(h1(1),'UserData');
                cori = cell2mat(cori(1,:));
                set(h1, 'FaceColor',cori, 'EdgeColor','k');
            end
        end

    end

    function out = fixDisplayName(hfig)

        if ~isempty(findobj(hfig, 'Tag', 'PlotMatrixBigAx'))

            % Get the handle to the legend and the legend entries
            hLegend        = findobj(hfig, 'Type', 'Legend');
            %legendEntries  = fliplr(hLegend.String);
            legendEntries  = hLegend.String;
            nlegendEntries = numel(legendEntries);

            % Get the handles of the panels
            axesObjects = findobj(hfig, 'Type', 'axes');
            
            % Go over the axes
            for i = 1:numel(axesObjects)
                % Not all the axes refer to a panel
                if isa(axesObjects(i).Children, 'matlab.graphics.chart.primitive.Line')
                    % disp(['Panel ', num2str(i), ' contains a line plot'])
                    % Assign the same DisplayName property to each group in each panel
                    for g = 1:nlegendEntries
                        axesObjects(i).Children(g).DisplayName = legendEntries{g};
                    end
                else
                    % disp(['Panel ', num2str(i), ' does not contain a line plot'])
                end
            end            

            % Get the input for the subsequent legend/clickableMultiLegend
            out  = findobj(axesObjects(1), 'type', 'line');

        else

            out = [];

        end

    end

end

%FScategory:UTIGEN
