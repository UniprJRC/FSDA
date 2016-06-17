function h = olsline(varargin)
%olsline adds least-squares fit line(s) to scatter plot(s).
%   It is a modified version of the standard function lsline, especially
%   conceived to work with gplotmatrix.
%   OLSLINE can be called without arguments or with a single argument,
%   which can be 0 or get an integer between 1 and the number of
%   data groups in each scatter plot.
%   The fit on each scatterplot in the figure is done as follows:
%   - if called without argument, OLSLINE fits a single ordinary least
%   squares line to all data, independently from the data groups in the
%   scatter plot.
%   - if the argument is 0, OLSLINE fits a line to each data group
%   in the scatter plot;
%   - if the argument is a positive index i >= 1, OLSLINE fits a single line
%   to the i-group in the scatter plot.
%
%   H = LSLINE returns the handle to the line object(s) in H.
%
% See also mdrplot.m, resplot.m
%
% References:
%
%   Riani M., Cerioli A., Perrotta D., Torti F. (2008). Fitting Mixtures of
%   Regression Lines with the Forward Search. Mining Massive Data Sets for
%   Security F. Fogelman-Soulié et al. EDS. (pp. 271-286). IOS Press,
%   Amsterdam (The Netherlands).
%
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('olsline')">Link to the help function</a>
% Last modified 31-05-2016

% Examples:

%{
    % A bivariate dataset with two groups of units
    load('fishery.txt','fishery.txt');
    y=fishery(:,2);
    X=fishery(:,1);
    group=ones(677,1); group(181:216)=2;

    % plot the data
    [H,AX,BigAx]=gplotmatrix(X,y,group);

    set(gcf,'CurrentAxes',AX)

    % fit a line to all data
    h=olsline;
    set(h,'DisplayName','fit on all units');

    % fit a line to one group
    h=olsline(size(H,3));
    set(h,'DisplayName','fit on group 1');

    % fit a line to another group
    h=olsline(size(H,2));
    set(h,'DisplayName','fit on group 2');

    % this just updates all legends
    hLines = findobj(gca, 'type', 'line');
    legend(hLines);
%}


%% Beginning of code

defaultColor = [.3 .3 .3] ; %Grey  [.75 .75 .75] == Light Gray
%      [.85 .64 .85] == Violet

if nargin == 0
    dataToFit = -1;
    with_intercept = 1;
elseif nargin == 1
    dataToFit = varargin{1};
    with_intercept = 1;
elseif nargin == 2
    dataToFit = varargin{1};
    with_intercept = varargin{2};
    if (with_intercept ~= 0) && (with_intercept ~= 1)
        with_intercept = 1;
        warning('FSDA:olsline:WrongIntercept','The intercept option can be either 0 or 1: it was set to 1 (default)');
    end
else
    error('FSDA:olsline:TooManyArgs','Too many arguments. First argument: \n 0 to add a line for each group; \n an index i >= 1 to add a single ols line to group i.');
end

% Find any line objects that are descendents of the axes.
AxCh = get(gca,'Children');
lh = findobj(AxCh,'Type','line');

% Ignore certain continuous lines.
if ~isempty(lh)
    style = get(lh,'LineStyle');
    if ~iscell(style)
        style = cellstr(style);
    end
    ignore = strcmp('-',style) | strcmp('--',style) | strcmp('-.',style);
    lh(ignore) = [];
end

% Find hggroups that are immediate children of the axes, such as plots made
% using SCATTER.
hgh = findobj(AxCh,'flat','Type','hggroup');
% Ignore hggroups that don't expose both XData and YData.
if ~isempty(hgh)
    ignore = ~isprop(hgh,'XData') | ~isprop(hgh,'YData');
    hgh(ignore) = [];
end

hh = [lh;hgh];
numlines = length(hh);
if numlines > 0
    hlslines = zeros(numlines,1);
    
    if dataToFit == -1
        % Fit all points in the scatter plot and plot the line in light Gray
        xdatone = []; ydatone = [];
        for k = 1:length(hh)
            % Extract data from the points we want to fit.
            xdat = get(hh(k),'XData'); xdat = xdat(:);
            ydat = get(hh(k),'YData'); ydat = ydat(:);
            xdatone = cat(1,xdatone,xdat);
            ydatone = cat(1,ydatone,ydat);
        end
        ok = ~(isnan(xdatone) | isnan(ydatone));
        if with_intercept
            beta = [ones(size(xdatone(ok,1))) xdatone(ok,:)]\ydatone(ok,:);
            intercept  = beta(1); slope = beta(2);
        else
            beta = xdatone(ok,:)\ydatone(ok,:);
            intercept  = 0; slope = beta;
        end
        xlimits = get(get(hh(1),'Parent'),'Xlim');
        xdat2 = xlimits;
        ydat2 = intercept + slope.*xdat2;
        hlslines(1) = line(xdat2,ydat2);
        datacolor = defaultColor;
        set(hlslines(1),'Color',datacolor,'LineWidth',2,'LineStyle','--');
        
    elseif dataToFit == 0
        % Fit the points of a given group and plot a line for each group.
        for k = 1:length(hh)
            % Extract data from the points we want to fit.
            xdat = get(hh(k),'XData'); xdat = xdat(:);
            ydat = get(hh(k),'YData'); ydat = ydat(:);
            if length(xdat)>1 % no line to fit if there is only 1 point.
                ok = ~(isnan(xdat) | isnan(ydat));
                if isprop(hh(k),'Color')
                    datacolor = get(hh(k),'Color');
                else
                    datacolor = defaultColor;
                end
                if with_intercept
                    beta = [ones(size(xdat(ok,1))) xdat(ok,:)]\ydat(ok,:);
                    intercept  = beta(1); slope = beta(2);
                else
                    beta = xdat(ok,:)\ydat(ok,:);
                    intercept  = 0; slope = beta;
                end
                
                %xlimits = get(gca,'Xlim');
                xlimits = get(get(hh(k),'Parent'),'Xlim');
                xdat2 = xlimits;
                ydat2 = intercept + slope.*xdat2;
                hlslines(k) = line(xdat2,ydat2);
                set(hlslines(k),'Color',datacolor,'LineWidth',2,'LineStyle','--');
            end
        end
    else
        % Fit the points of the group of index dataToFit and plot the line
        % for this group.
        xdat = get(hh(dataToFit),'XData'); xdat = xdat(:);
        ydat = get(hh(dataToFit),'YData'); ydat = ydat(:);
        if length(xdat)>1 % no line to fit if there is only 1 point.
            ok = ~(isnan(xdat) | isnan(ydat));
            if isprop(hh(dataToFit),'Color')
                datacolor = get(hh(dataToFit),'Color');
            else
                datacolor = defaultColor;
            end
            if with_intercept
                beta = [ones(size(xdat(ok,1))) xdat(ok,:)]\ydat(ok,:);
                intercept  = beta(1); slope = beta(2);
            else
                beta = xdat(ok,:)\ydat(ok,:);
                intercept  =0; slope = beta;
            end
            %xlimits = get(gca,'Xlim');
            xlimits = get(get(hh(dataToFit),'Parent'),'Xlim');
            xdat2 = xlimits;
            ydat2 = intercept + slope.*xdat2;
            hlslines(1) = line(xdat2,ydat2);
            set(hlslines(1),'Color',datacolor,'LineWidth',2,'LineStyle','--');
        end
    end
    hlslines(hlslines==0) = [];
    set(hlslines,'Tag','lsline');
else
    warning('FSDA:olsline:NoLinesFound','No allowed line types or scatterplots found. Nothing done.');
    hlslines = [];
end

if nargout == 1
    h = hlslines;
end
end
