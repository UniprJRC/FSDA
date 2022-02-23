function augStarplot(X,obslabs,varlabs)
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

%% Beginning of code

if nargin < 1
    error(message('FSDA:augStarplot:TooFewInputs'));
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
plotStars(Xstd(pagerows,:),ctrx,ctry,radius,axesh,varlabs);

% Clipping 'on' sets the text to the axes boundaries
texth = text(ctrx,ctry-1.1*radius,obslabs(pagerows),...
    'Clipping','on', 'HorizontalAlignment','Center', ...
    'Parent',axesh);
set(texth,'Tag','obs label');

end 


% -----------------------------------------
function plotStars(X,ctrx,ctry,radius,axesh,varlabs)
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

% The following line plots perimeter
% plot(axesh, tipx([1:p 1],:), tipy([1:p 1],:),'-','Color',colors(1,:),'Visible','off', plotArgs{:});

hold('on')
mycolors=[1 0 0;  % red
          0 0 1;   % blue
          0 0 0;   % black
          1 0 1;   % magenta 
          0.4660 0.6740 0.1880]; % dark green 
for i=1:p
    plot(axesh, ...
        spokesx(i*3-2:i*3-1,:), spokesy(i*3-2:i*3-1,:), '-', ...
        'Color',mycolors(i,:), 'LineWidth',2);
end
% Add center
plot(spokesx(1,:),spokesy(1,:),'o','MarkerFaceColor','b')

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


