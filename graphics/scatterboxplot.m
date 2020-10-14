function h = scatterboxplot(x,y,varargin)
%
%<a href="matlab: docsearchFS('scatterboxplot')">Link to the help page for this function</a>
%
%  scatterboxplot displays a 2D scatter plot with marginal boxplots. It
%  receives data in two vectors X and Y, and puts a univariate boxplot on
%  the horizontal and vertical axes of the plot. x and y must have the same
%  length.
%
%  Required input arguments:
%
%   x :          Input data. Vector or 1-column matrix. The structure
%                contains the univariate data to display in the x axis.
%
%   y :          Input data. Vector or 1-column matrix. The structure
%                contains the univariate data to display in the y axis.
%
%  Optional input arguments:
%
%   Group :     A grouping variable. Vector. Vector of group identifiers.
%                 With this option scatterboxplot creates a 2D GSCATTER
%                 plot instead of a SCATTER plot, and the marginal
%                 boxplots are replaced by grouped boxplots.
%                 Data Types - scalar
%                 Example - Group,[1,1,1,2,2,2,2,2,2]
%
%
%  Output:
%
%   h :         Vector of handles. Vector. Contains the three axes handles 
%               for the scatterplot, the plot along the horizontal axis, 
%               and the plot along the vertical axis, respectively..
%               Data Types - single | double.
% 
%
%  See also scatterhist, scatterhistogram
%
%
% References:
%
% Copyright 2008-2020.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('scatterboxplot')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2019-11-19 19:15:24 #$: Date of the last commit
%
%
% Examples:
%
%{
    %% A 2D scatter plot with marginal boxplots. 
    n=100;
    scatterboxplot(randn(n,1),randn(n,1))
%}
%
%{
    %% A 2D scatter plot with marginal boxplots, for grouped data.
    n=100; group=ones(n,1);
    group(1:50)=2; scatterboxplot(randn(n,1),randn(n,1),'Group',group)
%}

%% Start calling function scatterhist
if nargin<3
    h=scatterhist(x,y);
    group='';
else
    h=scatterhist(x,y,varargin{:,:});
    
    UserOptions=varargin(1:2:length(varargin));
    checkGroup = strcmp(UserOptions,'Group');
    if sum(checkGroup)
        group = varargin(2*find(checkGroup));
    else
        group='';
    end
end
if iscell(group)
    group=group{:};
end

hold on;
clr = get(h(1),'colororder');
boxplot(h(2),x,group,'orientation','horizontal',...
    'label','','color',clr);
set(h(2:3),'XTickLabel','');
boxplot(h(3),y,group,'orientation','horizontal',...
    'label', '','color',clr);
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto');    % Sync axes
hold off;

end
%FScategory:VIS-Mult
