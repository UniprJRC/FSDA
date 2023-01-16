function centers=augStarplot(X,obslabs,varlabs, varargin)
%augStarplot creates the Augmented star plot.
%
% This plot is useful to indicate which options are important in a particular
% data analysis. The rays in individual plots
% are of equal length for those features used in an analysis.
%
%
%  Required input arguments:
%
%       X   :  Input data. Matrix.
%               n x p data matrix; n observations and p variables. Rows of
%               X represent observations, and columns represent variables.
%
%   obslabs : labels for each star. String array, character array or cell
%             array of strings. Vector of length n containing the labels of
%             the stars.
%
%   varlabs : labels for the spokes of each star. String array, character
%             array or cell array of strings. Vector of length p containing
%             the labels of each slike of each star.
%
%
% Optional input arguments:
%
%  addPolygons : polygons around the outside. Boolean.
%               if addPolygons is true (default) , polygons around the
%               outside of the radial lines are added to the augmented star
%               plot. If addPolygons is false just the radial lines from
%               the center are shows.
%           Example - 'addPolygons',false
%           Data Types - logical
%
%      BestSols :  Best solutions. table.
%                 A table containing the details of the admissible
%                 solutions which have been found. 
%                 First column of BestSols contains R2       
%                 Second column of BestSols contains nused/n    
%                 Third column of BestSols contains  pvalDW   
%                 Fourth column of BestSols contains  pvalJB
%           Data Types - table


%% Beginning of code
% Create a figure to host the Augmented star plot
h=figure;
set(h,'Name', 'Augmented star plot', 'NumberTitle', 'off');

if nargin < 3
    error(message('FSDA:augStarplot:TooFewInputs'));
end
BestSols='';
addPolygons=true;

if nargin>3
    options=struct('BestSols',BestSols,'addPolygons',addPolygons);
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:augStarplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    BestSols=options.BestSols;
    addPolygons=options.addPolygons;
end
if isempty(BestSols)
    showBars=false;
else
    showBars=true;
end

n = size(X,1);

limX=[0 1];  % [min (max-min)] spoke length

obslabs = obslabs(:);
varlabs = varlabs(:);

if n == 0
    newplot
    % h = [];
    return
end

% Transform data

minX = min(X(:));
maxX = max(X(:));
if maxX > minX
    Xstd = limX(1) + limX(2)*(X - minX) ./ (maxX - minX);
else
    Xstd(:) = limX(1) + limX(2)*(X) ./ (maxX );
end


% Figure out where everything will go on the plot, and how big.  Unless
% this is a scrolling window, we'll plot in the current axes.
axesh = newplot;
% By default, put all the glyphs on one page.  ngridx-by-ngridy is big
% enough to fit everything.
ngridy = floor(sqrt(n));
ngridx = ceil(n./ngridy);
ngrid = ngridx*ngridy;
ctrx = repmat(1:ngridx,1,ngridy);
ctry = ngridy - reshape(repmat(1:ngridy,ngridx,1),1,ngrid) + 1;
radius = 0.4;
% xlim = [1-1.5*radius ngridx+1.5*radius];
% ylim = [1-1.5*radius ngridy+1.5*radius];


% Find the rows of X that will be displayed on this page.  On the last
% page, there may be fewer glyphs to plot than centers.
if ngrid <= n
    nglyphs = ngrid;
else % final page not full
    nglyphs = mod(n,ngrid);
    ctrx = ctrx(1:nglyphs);
    ctry = ctry(1:nglyphs);
end
pagerows = 1:nglyphs;

% Plot the grid of glyphs for the current page.
centers=plotStars(Xstd(pagerows,:),ctrx,ctry,radius,axesh,varlabs,addPolygons);

% Clipping 'on' sets the text to the axes boundaries
cy=ctry-1.1*radius;
h_axes = get(gcf,'CurrentAxes');    %get axes handle.

if showBars==true
    axesoffsets = get(h_axes,'Position');%get axes position on the figure.
    y_axislimits = get(gca, 'ylim');     %get axes extremeties.
    x_axislimits = get(gca, 'xlim');     %get axes extremeties.

    %get axes length

    for i=1:length(ctrx)
        y_axislength_lin = abs(y_axislimits(2) - y_axislimits(1));
        x_axislength_lin = abs(x_axislimits(2) - x_axislimits(1));
        y1 = axesoffsets(2)+axesoffsets(4)*abs((cy(i)-y_axislimits(1))/y_axislength_lin);
        x1 = axesoffsets(1)+axesoffsets(3)*abs((ctrx(i)-x_axislimits(1))/x_axislength_lin);

        sp = uipanel('Parent',gcf,'Title','','FontSize',12,...
            'Position',[x1-0.1 y1-0.10 .20 .1],'units','normalized');

        leftAxis = axes(sp, 'Position', [0 0 1 1],'units','normalized');
        % The order is  R2       nused/n    pvalDW     pvalJB
        bar(leftAxis,BestSols{i,[1 4 2 3]})
        ylim([0 1])
        texth = text(ctrx(i),cy(i)+0.1,['Sol' num2str(i)],...
            'Clipping','on', 'VerticalAlignment','top', ...
            'Parent',axesh);
    end
    % pause(0.1)
else
    % add text to each star
    texth = text(ctrx,cy,obslabs(pagerows),...
        'Clipping','on', 'HorizontalAlignment','Center', ...
        'Parent',axesh);
    axis(h_axes,'equal')
end
set(texth,'Tag','obs label');

end


% -----------------------------------------
function centers=plotStars(X,ctrx,ctry,radius,axesh,varlabs,addPolygons)
%PLOTSTARS Plot a grid of stars.

[n,p] = size(X);
theta = 2*pi*(0:(p-1))' ./ p;

ctrx = ctrx(ones(1,p),:);
ctry = ctry(ones(1,p),:);

% Set up the perimeter of each star as a column.
tipx = radius * repmat(cos(theta),1,n).*X' + ctrx;
tipy = radius * repmat(sin(theta),1,n).*X' + ctry;

% Set up spokes for each star as column with NaNs separating each spoke.
spokesx = cat(1, reshape(ctrx,[1 p n]), reshape(tipx,[1 p n]), NaN(1,p,n));
spokesx = reshape(spokesx, [3*p n]);
spokesy = cat(1, reshape(ctry,[1 p n]), reshape(tipy,[1 p n]), NaN(1,p,n));
spokesy = reshape(spokesy, [3*p n]);

mycolors=[1 0 0;  % red
    0 0 1;   % blue
    0 0 0;   % black
    1 0 1;   % magenta
    0.4660 0.6740 0.1880]; % dark green


if addPolygons==true
    % The following line plots perimeter
    % plot(axesh, tipx([1:p 1],:), tipy([1:p 1],:),'-','Color',colors(1,:),'Visible','off', plotArgs{:});
    plot(axesh, tipx([1:p 1],:), tipy([1:p 1],:),'-', ...
        'Color',mean(mycolors(:,:)),'Visible','on','LineWidth',2');
end

hold('on')
for i=1:p
    plot(axesh, ...
        spokesx(i*3-2:i*3-1,:), spokesy(i*3-2:i*3-1,:), '-', ...
        'Color',mycolors(i,:), 'LineWidth',2);
end
% Add center
plot(spokesx(1,:),spokesy(1,:),'o','MarkerFaceColor','b')

% Store centers
centers=[spokesx(1,:);spokesy(1,:)]';

for j=1:n
    boo=abs(tipy(:,j)-spokesy(1,j))>1e-05 | abs(tipx(:,j)-spokesx(1,j));
    text(tipx(boo,j),tipy(boo,j),varlabs(boo));
end
aa=gca;
aa.YLim(1)=0.3;
aa.XTickLabel='';
aa.XMinorTick='off';
aa.YMinorTick='off';

aa.YTickLabel='';

end % plotStars

% -----------------------------------------


