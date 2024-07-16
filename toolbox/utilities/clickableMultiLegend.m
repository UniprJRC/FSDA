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
    gscatter(X,y,group)
    % Make the legends clickable
    clickableMultiLegend();
%}

%{
     % clickableMultiLegend applied to two plots with lines with same legend names.
     % Clicking on a legend line, will show/hide the selected line in both plots. 
     z = peaks(100);
     plot(z(:,26:5:50))
     grid on;
     % define the names of the legend labels for both plots
     names = {'Line1','Line2','Line3','Line4','Line5'};
     clickableMultiLegend(names , 'Location', 'NorthWest');
     figure;
     plot(z(:,32:5:56))
     grid on;
     hlegend=clickableMultiLegend(names, 'Location', 'NorthWest');
%}

%{
    % clickableMultiLegend applied to multiple subplots.
    % For example let us start with a gplotmatrix.
    % Simulate X
    X = rand(100,4);

    % Simulate y with 3 groups
    y = [rand(10,1); rand(20,1)+1; rand(70,1)+2];

    group= [2*ones(10,1); ones(20,1); zeros(70,1)];
    % Generate the scatter matrix
    [H,AX,bigax] = gplotmatrix(X,y,group);

    % Set the DisplayName property (i.e. the texts of the legend) in all panels.
    % Note that in the gplotmatrix only one legend is visible.
    set(H(:,:,1),'DisplayName','group 1');
    set(H(:,:,2),'DisplayName','group 2');
    set(H(:,:,3),'DisplayName','group 3');

    % Get the handles of the legend to update
    hLines  = findobj(AX(1,end), 'type', 'line');

    % Update the legend and make them clickable
    clickableMultiLegend(hLines);

    % % % Get the new legend texts directly from the plot
    % % % To take account a change in property names of the legend object in 2016b
    % % if verLessThan('matlab','9.1')
    % %     legstring='LegendPeerHandle';
    % % else
    % %     legstring='LayoutPeers';
    % % end
    % % legnew = get(getappdata(AX(1,end),legstring),'String');
    % % 
    % % % Get the handles of the legend to update
    % % hLines  = findobj(AX(1,end), 'type', 'line');
    % % 
    % % % Update the legend and make them clickable
    % % clickableMultiLegend(sort(hLines), legnew{:});

    % Now, it is possible to click with the mouse on the different entries
    % to hide/show a particular group of units. For example, clicking on the
    % entry "group 2" in the legend we hide group 2.

    % Function gplotmatrix generates the legend texts automatically, based on
    % the values in the vector defined by option 'group'. In the example above we
    % have re-defined manually the legend texts set by option 'group' (which are
    % '1', '2' and '3') as "group 1", "group 2" and "group 3". More conveniently,
    % especially when the number of groups is not known in advance, one may
    % re-define the legend texts in a more general way as follows:

    % it is convenient to reshape the gplotmatrix handles array to make it
    % more manageable: while H is a 3-dimensional array with the third
    % dimension associated to the groups, newH is 2-dimensional with lines
    % associated to the subplots of the scatterplot and columns associated
    % to the groups.
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

% This preample documents an issue that concerns the legend function, which
% can be very slow because of a 'drawnow' call. Please report to the FSDA
% team any issue that might be related to this problem.
%
% We report below few guidelines on the problem mainly taken from
% http://undocumentedmatlab.com/articles/plot-performance of Yan Altman.
%
% Force the legend to be static.
%ax=gca;
% % set(ax,'LegendColorbarListeners',[]); This does not work anymore
%setappdata(ax,'LegendColorbarManualSpace',1);
%setappdata(ax,'LegendColorbarReclaimSpace',1);

% DrawMode is to avoid checking which objects need to be displayed on top
% of the others; Matlab will redraw objects following the order in which
% they were created.
% Clipping is to not waste time to check whether data beyond xlim/ylim need
% to be excluded from display.
% NextPlot is to avoid many automatic checking and property reset.
%set(ax,'DrawMode','fast','Clipping','off','NextPlot','replacechildren');

% Additional intervention: disabling legend for specific plot lines using:
% hasbehavior(hPlotLineToDisable,'legend',false);

%% Make the legend clickable

drawnow;

xlim manual;

lgd = legend(varargin{:});
set(lgd, 'ItemHitFcn', @(src, event)togglevisibility(event.Peer));
varargout={lgd};

drawnow;

axis manual;

    function togglevisibility(plotHandle)

        % whichFigs contains the handle(s) of the target figures. When it
        % is set to groot, clickableMultiLegend switches on/off all figures
        % containing objects with the same DisplayName. If it is set to,
        % for example, parentFigs (the ancestor), then clickableMultiLegend 
        % switches on/off only the active figure. Try the first example.
        whichFigs  = groot;
        %parentFig = ancestor(plotHandle,'Figure');

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

end

%FScategory:UTIGEN
