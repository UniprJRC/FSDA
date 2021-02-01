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
% Copyright 2008-2021.
% clickableMultiLegend has been adapted to this toolbox by FSDA team
%
%<a href="matlab: docsearchFS('clickableMultiLegend')">Link to the help page for this function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples

%{
     % clickableMultiLegend applied to a single plot.
     z = peaks(100);
     plot(z(:,26:5:50))
     grid on;
     clickableMultiLegend({'Line1','Line2','Line3','Line4','Line5'}, 'Location', 'NorthWest');
     axis manual;
     figure;
     z = peaks(100);
     plot(z(:,26:5:50))
     grid on;
     hlegend=clickableMultiLegend({'Line1','Line2','Line3','Line4','Line5'}, 'Location', 'NorthWest');
     axis manual;
     % legend(hlegend,'off');
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

    % Get the new legend texts directly from the plot
    % To take account a change in property names of the legend object in 2016b
    if verLessThan('matlab','9.1')
        legstring='LegendPeerHandle';
    else
        legstring='LayoutPeers';
    end
    legnew = get(getappdata(AX(1,end),legstring),'String');

    % Get the handles of the legend to update
    hLines  = findobj(AX(1,end), 'type', 'line');

    % Update the legend and make them clickable
    clickableMultiLegend(sort(hLines), legnew{:});

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

xlim manual;

% Create legend as if it was called directly
[varargout{1:nargout(@legend)}] = legend(varargin{:});

[~, objhan, plothan] = varargout{1:4};
varargout = varargout(1:nargout);

% Set the callbacks
for i = 1:length(plothan)
    if ~isempty(objhan)
        set(objhan(i), 'HitTest', 'on', 'ButtonDownFcn',...
            @(varargin)togglevisibility(objhan(i),plothan(i)),...
            'UserData', true);
    end
end
end

%% The callback function that shows/hides legends

function togglevisibility(hObject, obj)
if verLessThan('matlab','8.4.0')
    histobj='patch';
else
    histobj='Bar';
end

% hObject is the handle of the text of the legend
if get(hObject, 'UserData') % It is on, turn it off
    set(obj,'HitTest','off','Visible','off','handlevisibility','off');
    pause(0.001); %artificially introduced from 2016a to leave time to refresh HG2 object.
    set(hObject, 'Color', (get(hObject, 'Color') + 1)/1.5, 'UserData', false);
    
    similar_obj_h = findobj('DisplayName',get(obj,'DisplayName'));
    similar_obj_h(logical(similar_obj_h==obj)) = [];
    %similar_obj_h(find(similar_obj_h==obj)) = []; %slower than line before
    set(similar_obj_h,'HitTest','off','Visible','off','handlevisibility','on');
    
    % This is to make the patches of a group histogram white
    h = findobj('Type',histobj,'Tag',get(obj,'DisplayName'));
    
    % THIS SHOULD ALSO WORK
    %  h1 = findobj('-not','Type','Line','-not','Type','Axes','-and','Tag',get(obj,'DisplayName'));
    
    if ~isempty(h)
        set(h, 'UserData',get(h,'FaceColor'));
        set(h, 'FaceColor','w', 'EdgeColor','k');
    end
    
    % this is to hide the labels possibly associated to a group of units
    h_plo_labeladd = findobj('-regexp','Tag','plo_labeladd');
    if ~isempty(h_plo_labeladd)
        color_labeladd = get(h_plo_labeladd, 'Color');
        color_labeladd_1 = color_labeladd{1,:};
        if color_labeladd_1 == get(obj, 'Color')
            set(h_plo_labeladd,'HitTest','off','Visible','off','handlevisibility','on');
        end
    end
    
else
    
    set(hObject, 'Color', get(hObject, 'Color')*1.5 - 1, 'UserData', true);
    set(obj, 'HitTest','on','visible','on','handlevisibility','on');
    
    similar_obj_h = findobj('DisplayName',get(obj,'DisplayName'));
    similar_obj_h(logical(similar_obj_h==obj)) = [];
    %similar_obj_h(find(similar_obj_h==obj)) = []; %slower than line before
    set(similar_obj_h,'HitTest','on','Visible','on','handlevisibility','on');
    
    % This is to re-establish the color of the white patches of a group histogram
    h = findobj('Type',histobj,'Tag',get(obj,'DisplayName'));
    
    % THIS SHOULD ALSO WORK
    % h = findobj('-not','Type','Line','-not','Type','Axes','-and','Tag',get(obj,'DisplayName'));
    
    if ~isempty(h)
        cori = get(h(1),'UserData'); cori = cori{1};
        set(h, 'FaceColor',cori, 'EdgeColor','k');
    end
    
    % this is to show the labels possibly associated to a group of units
    h_plo_labeladd = findobj('-regexp','Tag','plo_labeladd');
    if ~isempty(h_plo_labeladd)
        color_labeladd = get(h_plo_labeladd, 'Color');
        color_labeladd_1 = color_labeladd{1,:};
        if color_labeladd_1 == get(obj, 'Color')
            set(h_plo_labeladd,'HitTest','on','Visible','on','handlevisibility','on');
        end
    end
    
end
end
%FScategory:UTIGEN