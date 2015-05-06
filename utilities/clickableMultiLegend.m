function varargout = clickableMultiLegend(varargin)
%clickableMultiLegend extends the clickableLegend by Ameya Deoras to figures with one legend for several subplots
%
%<a href="matlab: docsearchFS('clickableMultiLegend')">Link to the help page for this function</a>
%
% An example could be a gplotmatrix figure. In this case, by clicking on a
% text label in the legend, we want to turn on and off (hide or show) the
% graphics object (line or patch) associated to that label in all subplots.
%
% The extention to multiple plots is realised by looking for graphics
% objects with the same DisplayName property of the one associated to the
% legend label. Therefore, the function should work also through plots in
% different figures.
%
% See also
% clickableLegend by Ameya Deoras:
% http://www.mathworks.com/matlabcentral/fx_files/21799/1/clickableLegend.m
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('clickableMultiLegend')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples

%{
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


%% Create legend as if it was called directly
[varargout{1:nargout(@legend)}] = legend(varargin{:});

[~, objhan, plothan] = varargout{1:4};
varargout = varargout(1:nargout);

% Set the callbacks
for i = 1:length(plothan)
    set(objhan(i), 'HitTest', 'on', 'ButtonDownFcn',...
        @(varargin)togglevisibility(objhan(i),plothan(i)),...
        'UserData', true);
end


% The callback function that shows/hides legends

function togglevisibility(hObject, obj)
if verLessThan('matlab','8.4.0')
    histobj='patch';
else
    histobj='Bar';
end

% hObject is the handle of the text of the legend
if get(hObject, 'UserData') % It is on, turn it off
    set(hObject, 'Color', (get(hObject, 'Color') + 1)/1.5, 'UserData', false);
    set(obj,'HitTest','off','Visible','off','handlevisibility','off');
    
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
