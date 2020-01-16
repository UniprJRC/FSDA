function varargout = dsxy2figxy(varargin)
%dsxy2figxy transforms points or positions from 'axes data units' to 'normalized figure units'
%
%
% Required input arguments:
%
%   [figx figy] = dsxy2figxy([x1 y1],[x2 y2])
%   figpos = dsxy2figxy([x1 y1 width height])
%
%       x1,x2,y1,y2:  scalars (in plot coordinates)
%       x1 y1 width height: scalars (in plot coordinates)
%
% Output:
%
% dsxy2figxy transforms [axx axy] or [xypos] from axes hAx (data) coords
% into coords wrt GCF for placing annotation objects that use figure coords
% into data space. The annotation objects that use just normalized coordinates are:
%    arrow, doublearrow, textarrow
%    ellipses (in this case coordinates must be transformed to [x, y,
%    width, height] where x and y, width and height are numbers between 0
%    and 1)
%
% Consider that line, text, and rectangle anno objects already are placed
% on a plot using axes coordinates and must be located within an axes.
%
%In other words, for example
%    rectangle('Position',[x,y,w,h]) 
%   draws a rectangle from the
% point x,y  having a width of w and a height of h.
% Note that in this case the coordinates are specified in axes data units.
% On the other hand
%   annotation('line',x,y) 
% creates a line annotation object that
% extends from the point defined by x(1),y(1) to the point
% defined by x(2),y(2), specified in normalized figure units
%
% See also 
%
%
% Copyright 2008-2019.
% REMARK: function  dsxy2figxy has been written by MATLAB programmers but for some
% strange reason has not included yet into the standard documentation
%
%
%
%<a href="matlab: docsearchFS('dsxy2figxy')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

%{
    %% Plot of sin between -pi and pi 
    % Add an arrow from (0,0) to say (pi,1):
    x = -pi:pi/10:pi;
    plot(x,sin(x));
    % set(gcf,'Units','normalized');
    % set(gcf,'Units','points');
    [figx figy] = dsxy2figxy([0 pi],[0 1]);
    annotation('textarrow',figx,figy)
%}

%% Obtain arguments (only limited argument checking is performed).

% Determine if axes handle is specified
if length(varargin{1})== 1 && ishandle(varargin{1}) && ...
  strcmp(get(varargin{1},'type'),'axes')	
	hAx = varargin{1};
	varargin = varargin(2:end);
else
	hAx = gca;
end

% Parse either a position vector or two 2-D point tuples
if length(varargin)==1	% Must be a 4-element POS vector
	pos = varargin{1};
else
	[x,y] = deal(varargin{:});  % Two tuples (start & end points)
end

%% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Need normaized units to do the xform
axpos = get(hAx,'Position');
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));

%% Transform data from figure space to data space
if exist('x','var')     % Transform a and return pair of points
	varargout{1} = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
	varargout{2} = (y-axlim(3))*axpos(4)/axheight + axpos(2);
else                    % Transform and return a position rectangle
	pos(1) = (pos(1)-axlim(1))/axwidth*axpos(3) + axpos(1);
	pos(2) = (pos(2)-axlim(3))/axheight*axpos(4) + axpos(2);
	pos(3) = pos(3)*axpos(3)/axwidth;
	pos(4) = pos(4)*axpos(4)/axheight;
	varargout{1} = pos;
end

%% Restore axes units
set(hAx,'Units',axun)
end

