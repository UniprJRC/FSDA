function [haxis,hlabel]=suplabel(text,whichaxis,possuperaxes)
%suplabel places text as a title, xlabel, or ylabel on a group of subplots.
%
%
%<a href="matlab: docsearchFS('suplabel')">Link to the help function</a>
%
% Required input arguments:
%
%   text       : any string. Character or cell array of string.
%               If text is a Character it containg the string which has to put as title,
%               xlabel or ylabel on a group of subplots. If text is a cell
%               array of strings each element of the cell is put on a
%               different line. In other words, in this case suplabel will
%               display multiline titles contained in a cell array
%
%
% Optional input arguments:
%
%   whichaxis   : where to put the string. Character.
%                 String equal to any of 'x', 'y', 'yy', or 't', specifying
%                 whether the text is to be the xlabel, ylabel, right
%                 side y-label, or title respectively. If whichaxis is not
%                 specified it is set to 'x'
%                 Example - 'y'
%                 Data Types - character
%
% possuperaxes  : super axes position. double. 
%                 vector of length 4 which specifies the
%                 position of the "super" axes surrounding the subplots.
%                 The default values of possuperaxes is [.08 .08 .84 .84];
%                 The meaning of the four elements of possuperaxes are
%                 [left bottom width height] where left and bottom define
%                 the distance from the lower-left corner of the container
%                 to the lower-left corner of the rectangle. width and
%                 height are the dimensions of the rectangle. The Units
%                 property specifies the units for all measurements.
%                 Example - [.08 .10 .84 .84]
%                 Data Types - numeric vector of length 4
%
%
% Output:
%
%      haxis: handle to the axis. Graphics handle.
%             Graphics handle to the axis.
%      hlabel: handle to the label. Graphics handle.
%             Graphics handle to the label.
%
% See also: subplot, spmplot, yXplot
%
% References:
%
% Acknowledgements: 
%  This code is based on the code suplabel by Ben Barrowes
% <barrowes@alum.mit.edu>
% https://www.mathworks.com/matlabcentral/fileexchange/7772-suplabel?s_tid=srchtitle
%
% 
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('suplabel')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples

  
%{
  % Add a top, let, right and bottom titles to a gplotmatrix.  
  % Creare a scatter plot matrix with three variables and add gloabal X
  % label on the X axis, global Y label on the left and right, and global
  % title.
  n=200;
  p=3;
  state1=123498;
  randn('state', state1);
  X=randn(n,p);
  gplotmatrix(X);
  % Add a common label on the x axis
  [ax1,h1]=suplabel('super X label');
  % Add a common label on the y axis
  [ax2,h2]=suplabel('super Y label','y');
  % Add a common label on the y axis
  [ax3,h2]=suplabel('super Y label (right)','yy');
  % Add a common label on top of the plot
  [ax4,h3]=suplabel('super Title'  ,'t');
  % set the fontsize of the string on top of the plot
  set(h3,'FontSize',30);
%}

%{
    % Two panels with a common y label.
    figure
    subplot(2,1,1);
    plot((1:10).^2)
    subplot(2,1,2);
    plot((1:10).^3)
    suplabel('Population growth','y')
%}

%{
    % Example with 6 panels.
    % The three panels of the left have a common xlabel, ylabel and
    % right ylabel
    figure
    subplot(3,2,1);
    plot((1:10).^2)
    subplot(3,2,3);
    plot((1:10).^2)
    subplot(3,2,5);
    plot((1:10).^2)
    possuperaxes=[0.1 0.1 0.35 0.8];
    suplabel('Population growth','y',possuperaxes)
    suplabel('Right label','yy',possuperaxes)
    suplabel('Years','x',possuperaxes)
    
    % The three panels of the right have a common xlabel, ylabel and
    % right ylabel
    subplot(3,2,2);
    plot((1:10).^2)
    subplot(3,2,4);
    plot((1:10).^2)
    subplot(3,2,6);
    plot((1:10).^2)
    possuperaxes=[0.55 0.1 0.3 0.8];
    suplabel('Population growth 2','y',possuperaxes)
    suplabel('Right label','yy',possuperaxes)
    suplabel('Months','x',possuperaxes)
%}


%{
    % Example of suplabel with output arguments.
    % A label is added to the spmplot on the x axis
    load fisheriris;
    plo=struct;
    plo.nameY={'SL','SW','PL','PW'};
    spmplot(meas,species,plo,'hist');
    % insert text 'any string' as x label
    % and return both the handle to the axis (inside haxis) and the handle to the label (inside hlabel).
    [haxis,hlabel]=suplabel('Title added to the x axis')  
%}

%{
    % Call to subplot using as text a cell array of strings.
    figure
    subplot(3,2,1);
    plot((1:10).^2)
    subplot(3,2,3);
    plot((1:10).^2)
    subplot(3,2,5);
    plot((1:10).^2)
    possuperaxes=[0.1 0.1 0.35 0.8];
    text={'Population growth','3 countries'} ;
    suplabel(text,'y',possuperaxes)
%}

%% Beginning of code
if nargin<1
    error('FSDA:suplabel:Missingtext','text to add is missing');
end

% Set default position where to put the text
if nargin < 2 
    whichaxis = 'x';  
end

% Check that text is a character or a cell array of strings and text is a
% character
if ((iscell(text) && ~ischar([text{:}])) || (~iscell(text) &&~ischar(text) )) || ~ischar(whichaxis)
    error('FSDA:suplabel:WrongInput','text must be a string and whichLabel must be either a text or a cell array of strings')
end

whichaxis=lower(whichaxis);

% Find all graphics objects of type axes exclusing those which have tag
% suplabel
currax=findobj(gcf,'type','axes','-not','tag','suplabel');

if nargin < 3
    possuperaxes=[.08 .08 .84 .84];
    ah=findall(gcf,'type','axes');
    if ~isempty(ah)
        leftMin=inf;  bottomMin=inf;  leftMax=0;  bottomMax=0;
        axBuf=.04;
        set(ah,'units','normalized')
        ah=findall(gcf,'type','axes');
        for ii=1:length(ah)
            if strcmp(get(ah(ii),'Visible'),'on')
                thisPos=get(ah(ii),'Position');
                leftMin=min(leftMin,thisPos(1));
                bottomMin=min(bottomMin,thisPos(2));
                leftMax=max(leftMax,thisPos(1)+thisPos(3));
                bottomMax=max(bottomMax,thisPos(2)+thisPos(4));
            end
        end
        possuperaxes=[leftMin-axBuf,bottomMin-axBuf,leftMax-leftMin+axBuf*2,bottomMax-bottomMin+axBuf*2];
    end
end

% creates Cartesian axes in the current figure in position specified by
% possuperaxes and make them invisible. This is the trick to add xlabel,
% ylabels, title... to a series of multiple plots
haxis=axes('Units','Normal','Position',possuperaxes,'Visible','off','tag','suplabel');

if strcmp('t',whichaxis)   % Global title
    set(get(haxis,'Title'),'Visible','on')
    title(text);
elseif strcmp('x',whichaxis) % Global x label
    set(get(haxis,'XLabel'),'Visible','on')
    xlabel(text);
elseif strcmp('y',whichaxis) % Global y label
    set(get(haxis,'YLabel'),'Visible','on')
    ylabel(text);
elseif strcmp('yy',whichaxis) % Global y label (on the right(
    set(get(haxis,'YLabel'),'Visible','on')
    ylabel(text);
    set(haxis,'YAxisLocation','right')
else
    error('FSDA:suplabel:WrongInput','String label must be any of ''x'', ''y'', ''yy'', or ''t''')
end

% The loop seems to be unnecessary and in recent MATLAB (version greater
% than 2014a) can simply be replaced by. 
% axes(currax);
% However, for compatibility with older version we leave it
for k=1:length(currax)
    axes(currax(k)) %#ok<LAXES>
end % restore all other axes

if (nargout < 2)
    return
end

if strcmp('t',whichaxis)
    hlabel=get(haxis,'Title');
    set(hlabel,'VerticalAlignment','middle')
elseif strcmp('x',whichaxis)
    hlabel=get(haxis,'XLabel');
elseif strcmp('y',whichaxis) || strcmp('yy',whichaxis)
    hlabel=get(haxis,'YLabel');
end

end
%FScategory:UTIGEN