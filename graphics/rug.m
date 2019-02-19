function varargout=rug(varargin)
% Adds "rug" tick marks to an existing plot
%
% function handles=rug(ax,tickScaleFactor)
%
%
% Purpose
% Creates a so-called "rug" plot similar to those produced by the
% free stats package, R. Rug plots use tick marks along the axes as a
% compact way of illustrating the marginal distributions of a variable 
% in a scatter plot. These tick marks are reminiscent of the tassels on a 
% rug. 
%
% Inputs
% - ax [optional] - handle defining the axes on which rug will work. If ax
%    is missing, rug is applied to the current axis. ax can also be a 
%    vector of axis handles. 
% - tickScaleFactor [optional] - a scalar defining the length of the tick 
%    marks as proportion of each axis' length. tickScaleFactor must have a 
%    value between 0 and 1. Defaults to 0.015. 
% Note, that input arguments can be supplied in any order. 
%
%
% Outputs
% - handles [optional] - plot handles for each set of rug objects.
%
% Examples
% - Rug plot with two populations denoted by different colours
%{
clf
  plot(randn(1,100)-3,randn(1,100)+3,'ko','markerfacecolor',[.5,.5,.5])
  hold on
  plot(randn(1,100),randn(1,100),'ro','markerfacecolor',[1,.5,.5])
 rug
%}
%
% - Add rugs to defined subplots
%{
     clf
     for ii=1:4
     subplot(2,2,ii), plot(randn(1,50),randn(1,50),'.k')
     end
     c=get(gcf,'children');
     rug([c(1),c(4)])
     rug(0.1,c(2))
%}

% Known Issues:
% Overlays rug tick marks onto the same axes as the original
% data. Changing axis dimensions after calling rug will therefore
% cause the tick marks to become disassociated from the axes. 
%
%
% Rob Campbell - May 2010


% Check input arguments and use default values if needed
error(nargchk(0,2,nargin))

ax=gca;
defaultTickScale=0.015;
tickScaleFactor=defaultTickScale;

if nargin>0
    for ii=1:length(varargin)
        if length(varargin{ii})==1 && ~ishandle(varargin{ii}) && ~isempty(varargin{ii})
            tickScaleFactor=varargin{ii};
        else
            ax=varargin{ii};
        end    
    end    
end

if tickScaleFactor>1 || tickScaleFactor<0
    fprintf('tickScaleFactor should be between 0 and 1. Defaulting to %0.3f\n',...
        defaultTickScale)
    tickScaleFactor=defaultTickScale;
end


%Check that plot axes are of a sort which we can work with (e.g. we
%don't want to process a figure legend). 
for ii=length(ax):-1:1

    if isempty(strmatch(get(gca,'Type'),'axes'))
        ax(ii)=[];
    end
    if isfield(get(ax(ii)),'Location')
        ax(ii)=[];
    end

end

if isempty(ax)
    error('No suitable axes found')
end




%Add rug tick marks to each pair of axes
n=1;
for ii=1:length(ax)
    set(gca,'TickDir','out') %otherwise the rug plot won't make sense

    axes(ax(ii)) %Focus on this pair of axes
    %To ensure we don't change the hold status by running this function
    holdStatus=ishold;
    if ~holdStatus, hold on,  end


    %Add rug tick marks to all axis elements of type "line"
    chil=get(ax(ii),'children'); 
    for jj=1:length(chil)
        H=get(chil(jj));
        if strmatch(H.Type,'line')
            handles(n)=addRug(ax(ii),H,tickScaleFactor);
            n=n+1;
        end
    end


    if ~holdStatus, hold off, end

end



%Return handles if user asks for them
if nargout==1
    varargout{1}=handles;
end








%----------------------------------------------------------------
function h=addRug(ax,H,tickScaleFactor)
%Explicitly setting the axis limits seems to stop them being    
%automatically re-set.
X=xlim; xlim(X)
Y=ylim; ylim(Y)

h.ax=ax;

tickProps={'-','color',H.Color};

%Tick marks along x axis
tickSize=diff(Y)*tickScaleFactor;
xd=repmat(H.XData,3,1);
xd(1,:)=nan;
yd=repmat([Y(1),Y(1),Y(1)+tickSize],1,length(H.XData));
h.xRug=plot(xd(:),yd(:),tickProps{:});

%Tick marks along y axis
tickSize=diff(X)*tickScaleFactor;
yd=repmat(H.YData,3,1);
yd(1,:)=nan;
xd=repmat([X(1),X(1),X(1)+tickSize],1,length(H.YData));
h.yRug=plot(xd(:),yd(:),tickProps{:});
