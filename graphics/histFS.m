function [ng, hb] = histFS(y,nbins,gy,gylab,ax,barcolors,W)
%histFS plots a histogram with the elements in each bin grouped according to a vector of labels.
%
%<a href="matlab: docsearchFS('histFS')">Link to the help function</a>
%
%  Required input arguments:
%
%     y        : vector of n elements to bin. Vector. Vector which contains
%                the elements which have to be binned.
%     nbins    : the number of bins. Scalar. The elements of of input vector
%                y are binned into nbins equally spaced containers.
%     gy       : identifier vector. Vector. Vector of n numeric identifiers
%                for the group of each element in y. If there are k groups
%                in y, unique(gy) = k.
%
%  Optional input arguments:
%
%     gylab     : legend labels. String | cell of strings.
%                 Legend labels identifying the groups of each element in y.
%                 length(gylab) = length(unique(gy)) must hold.
%                 If not specified, gylab is set to '' and legends are not
%                 displayed.
%                 If gylab = {}, default legend labels are generated.
%                 If gylab is a cell of strings of length 3, e.g. gylab =
%                 {'G1' 'G2' 'G3'}, such strings are used for the legends.
%               Example - {'G1' 'G2'}
%               Data Types - cell array of strings or char
%
%     ax        : plots into ax instead of gca. Axis handle. The axis handle
%                 where to plot the grouped histogram (e.g. a traditional
%                 histogram plot to be superimposed). Default is gca.
%               Example - gca
%               Data Types - graphics handle
%
%     barcolors : colors of the bars. char or matrix.
%                Vector containing the strings of the colors to use 
%                (e.g. 'rgy') or RGB matrix of the colors used for the
%                groups (e.g. [1 0 0; 0 0 1]). If the number of colors
%                supplied is smaller than the number of groups the program
%                displays an error.
%               Example - 'rgy'
%               Data Types - character or numeric matrix
%
%     W         : Weights. Vector. Vector which contains optional weights
%                 associated to the elements of y.
%               Example - [10 20 3000]
%               Data Types - numeric vector
%  Output:
%
%       ng      : number of elements in each container for each group.
%                 Matrix. A matrix with a column for each group and a row for
%                 each bin. In ng, column i contains the the number
%                 of elements of group i in each bin.
%                 For example ng(2,3) contains the number of elements
%                 belonging to the third group included into the second
%                 container.
%        hb     : Bar array handles. Vector. A vector containing the
%                 handles to the barseries objects.
%
% See also hist
%
% References:
%
%   Tufte E.R. (1983), "The visual display of quantitative information",
%   Graphics Press, Cheshire.

%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('histFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
      %% An example with 4 groups.
      close all;
      y = randn(500,1);
      % four groups
      groups = randi(4,500,1);
      % number of bins
      nbins = 10;
      [ng, hb] = histFS(y,nbins,groups);
      title('Histogram with default color rotation (brcmykg)','interpreter','latex','FontSize',18);

%}

%{
      %% The same histogram is now plotted with different legends.
      figure;
      [ng, hb] = histFS(y,nbins,groups,{'BMW','FIAT','VOLVO','FERRARI'});
      title('The same, with personalized legends','interpreter','latex','FontSize',18);
%}

%{
      %% The same histogram is now plotted with different colors.
      figure;
      [ng, hb] = histFS(y,nbins,groups,{'BMW','FIAT','VOLVO','FERRARI'},gca,'kgbr');
      title('FERRARI must be red! Color sequence changed to kgbr','interpreter','latex','FontSize',18);
%}

%{
      % Apply to the grouped histogram the legends of a different plot.
      % Create a scatterplot
      close all;
      hs = gscatter(1:numel(y),y,groups);
      hfs = gcf;                              % get the handle of the scatterplot
      has = get(hfs,'CurrentAxes');           % it is the same as has = gca
      hlegends  = get(has,'Children');        % get the legend entries
      getleg = get(hlegends,'DisplayName');   % get the names of the legend entries
      getcol = get(hlegends,'Color');         % get the color of the legend entries
      getcolm = cell2mat(getcol);             % arrange the RGB vectors into a matrix

      figure;
      [ng, hb] = histFS(y,nbins,groups,getleg,gca,getcolm);
      title('Color sequence is taken from the scatterplot','interpreter','latex','FontSize',18);
%}

%{
      % Example with personalised clickable legends.
      myleg = {'my group 1' 'my group 2' 'my group 3' 'my group 4' };
      [ng, hb] = histFS(y,nbins,groups,myleg);
      title('Example with personalised clickable legends','interpreter','latex','FontSize',18);
%}

%{
        %% Example of weighted histogram.
        close all;clear all;
        X = crosstab2datamatrix([10 20 30]); X=X(:,2);
        X1 = 444*ones(60,1); X1(40:60)= 999; X1 = shuffling(X1);
        X = X+5;
        if verLessThan('matlab','8.5')
            groups = num2str(X1);
        else
            groups = categorical(X1);
        end
        W = ones(size(X,1),1); W(12) = 100;

        histFS(X,10,groups,[],[],[],W);
        title({'Weighted histograms:' , 'please check with tabulateFS(W) and tabulateFS(X)'},'interpreter','latex','FontSize',18);
        figure;
        histFS(X,10,groups);
        title('Standard histogram','interpreter','latex','FontSize',18);
        cascade;
        disp('tabulateFS(X)');
        tabulateFS(X);
        disp('tabulateFS(W)');
        tabulateFS(W);
%}

%% Beginning of code
if nargin < 7 || not(isvector(W) && ~isscalar(W))
    W = [];
end
if nargin < 6 || isempty(barcolors)
    barcolors = 'brcmykgbrcmykgbrcmykg';
end
if nargin <5 || isempty(ax)
    ax = gca;
end
if nargin < 4 || isempty(gylab)
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

ng  = zeros(nbins,ngroups);
Wempty = isempty(W);

for i = 1 : nbins
    bing = gy(bins==i); % groups in bin i
    if ~Wempty
        bingW=W(bins==i);
    end
    
    for j=1:ngroups
        
        if Wempty
            % build matrix of labels of groups in each bin
            ng(i,j)=sum(bing==groups(j));
        else
            % build matrix of weights of groups in each bin
            ng(i,j)=sum(bingW(bing==groups(j)));
        end
    end
end

% prepare the colormap for the histogram
if isnumeric(barcolors)
    Crgb = barcolors;
    assert(size(Crgb,1) >=ngroups,'Number of supplied colors smaller than number of groups')
else
    C = textscan(barcolors, '%1c');
    C = C{:};
    assert(size(C,1) >= ngroups ,'Number of supplied colors smaller than number of groups')
    
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

% the use of colormap with the bar object seems to not work in recent
% matlab releases. Read the comments of the next code segments for details.
% colormap(ax,Crgb(1:ngroups,:));


% Draw a bar for each element in bm at locations specified in x.
% Notice that x is a vector defining the x-axis intervals for the vertical
% bars. Function bar groups the elements of each row in bm at corresponding
% locations in x. 
hb = bar(ax,x,ng,'stacked','FaceColor','flat','BarWidth',1);

% In the next 4 lines we color the bar by setting with a loop the FaceColor
% property. An alternative which was working with previous MATLAB versions
% was to set the colormap with the instruction commented above.
%hb = bar(ax,x,ng,'stacked','FaceColor','none','BarWidth',1);
for i=1:numel(hb)
    set(hb(i),'FaceColor',Crgb(i,:));
    %set(hb(i),'FaceColor',C(i));
end

% add clickable legends, so that to show/hide grouped bins
if doleg
    caxis(caxis);
    clickableMultiLegend(hb, gylab{:});
    axis(axis);
    for i=1:numel(hb)
        set(findobj(hb(i),'Type','patch'),'DisplayName',gylab{i},'EdgeColor','k');
    end
end

end
%FScategory:VIS-Reg