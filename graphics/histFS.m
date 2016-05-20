function [ng, hb] = histFS(y,nbins,gy,gylab,ax,barcolors)
%histFS plots a histogram with the elements in each bin grouped according to a vector of labels. 
%
%<a href="matlab: docsearchFS('histFS')">Link to the help function</a>
%
%  Required input arguments
%
%     y        : vector of n elements.
%     nbins    : the number of bins.
%     gy       : vector of n numeric identifiers for the group of each 
%                element in y. If there are k groups in y, unique(gy) = k.
%
%  Optional input arguments
%
%     gylab     : legend labels identifying the groups of each element in y.
%                 size(gylab) = size(gy) must hold.
%                 If not specified, gylab is set to '' and legends are not
%                 displayed. 
%                 If gylab = {}, default legend labels are generated.
%                 If gylab is a cell of strings of lenght n, e.g. gylab =
%                 {'G1' 'G2' 'G3'}, such strings are used for the legends.
%     ax        : the axis handle where to plot the grouped histogram (e.g.
%                 a traditional histogram plot to be superimposed). Default
%                 gca.
%     barcolors : Personalised RGB matrix of the colors used for the groups.
%
%  Output
%
%     ng        : a matrix with a column for each group and a row for 
%                 each bin. In ng, the column i contains the the number
%                 of elements of group i in each bin.
%     hb        : a vector hb with the handles to the barseries objects.
%
%<a href="matlab: docsearchFS('histFS')">Link to the help function</a>
%
% Example:
%{
      y = randn(500,1);
      % three groups
      groups = randi(4,500,1);
      % number of bins
      nbins = 10;
      [ng, hb] = histFS(y,nbins,groups);
%}
%{
      % Now with default clickable legends
      [ng, hb] = histFS(y,nbins,groups,{});
%}
%{
      % Now with personalised clickable legends
      myleg = {'my group 1' 'my group 2' 'my group 3' 'my group 4'};
      [ng, hb] = histFS(y,nbins,groups,myleg);
%}
%{
      % Now apply to the grouped histogram the legends of a different plot, 
      
      hs = gscatter(1:numel(y),y,groups);     % e.g. a scatterplot

      hfs = gcf;                              % get the handle of the scatterplot
      has = get(hfs,'CurrentAxes');           % it is the same as has = gca
      hlegends  = get(has,'Children');        % get the legend entries
      getleg = get(hlegends,'DisplayName');   % get the names of the legend entries
      getcol = get(hlegends,'Color');         % get the color of the legend entries
      getcolm = cell2mat(getcol);             % arrange the RGB vectors into a matrix

      figure;                                 
      [ng, hb] = histFS(y,nbins,groups,getleg,gca,getcolm); 

%}
%   See also hist, bar.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
% Last modified 06-Feb-2015

%% Beginning of code
if nargin < 6
    barcolors = 'brcmykgbrcmykgbrcmykg';
end
if nargin < 5
    ax = gca;
end
if nargin < 4
    gylab = '';
end

% Set groups and number of groups.
% REMARK: the  vector resulting from unique is sorted in ascending order:
% in histFS this is required to ensure to have always the same order of
% colors
groups = unique(gy); 
ngroups = numel(groups);

if nargin >= 4 
    if ischar(gylab) && isempty(gylab)
        doleg = 0;
    end
    if iscell(gylab) 
        if isempty(gylab)
            gylab = cell(1,ngroups); 
            for i = 1 : ngroups
                gylab(1,i) = {['Group n. ' num2str(i)]};
            end
            doleg = 1;
        elseif numel(gylab) == ngroups
            doleg = 1;
        else
            doleg = 0;
            disp('Warning: number of labels and groups do not match');
            disp('Legends are not displayed');
        end
    end
else
    doleg = 0;
end

% size of vector
y = y(:); 
n = size(y,1);

[~,x]   = hist(y,nbins);       % bins locations
xc      = repmat(x,n,1);       % repeated bin locations
yr      = repmat(y,1,nbins);   % repeated data for each bin
dist    = abs(yr-xc);          % distance between data and bin centers

% use dist to assign data to bins
[~,bins] = min(dist,[],2);

% build matrix of labels of groups in each bin
ng  = zeros(nbins,ngroups);
for i = 1 : nbins
    bing = gy(bins==i); % groups in bin i
    for j=1:ngroups
        ng(i,j)=numel(find(bing==groups(j)));
    end
end

% prepare the colormap for the histogram
if isnumeric(barcolors)
    Crgb = barcolors;
else
    C = textscan(barcolors, '%1c');
    C = C{:};
    Crgb = zeros(size(C,1),3);
    for i=1:size(C,1)
        switch C(i)
            case 'y'
                rgb = FSColors.yellow.RGB;
            case 'm'
                rgb = FSColors.magenta.RGB;
            case 'c'
                rgb = FSColors.cyan.RGB;
            case 'r'
                rgb = FSColors.red.RGB;
            case 'g'
                rgb = FSColors.green.RGB;
            case 'b'
                rgb = FSColors.blue.RGB;
            case 'w'
                rgb = FSColors.white.RGB;
            case 'k'
                rgb = FSColors.black.RGB;
        end
        Crgb(i,:)=rgb;
    end
end
colormap(ax,Crgb(1:ngroups,:));

% Draw a bar for each element in bm at locations specified in x.
% Notice that x is a vector defining the x-axis intervals for the vertical
% bars. Function bar groups the elements of each row in bm at corresponding
% locations in x.

hb = bar(ax,x,ng,'stacked','FaceColor','flat','BarWidth',1);
% REMARK: as an alternative to the above Colormap instruction, one can
% create an empty bar and color it by setting with a loop the FaceColor
% property. This is shown in the next 4 lines.
% hb = bar(ax,x,ng,'stacked','FaceColor','none','BarWidth',1);
% for i=1:numel(hb)
%     set(hb(i),'FaceColor',Crgb(i,:));
% end

% add clickable legends, so that to show/hide grouped bins
if doleg
    caxis(caxis);
    clickableMultiLegend(hb, gylab{:});
    axis(axis);
    for i=1:numel(hb)
         set(findobj(hb(i),'Type','patch'),'DisplayName',gylab{i});
    end
end

end
