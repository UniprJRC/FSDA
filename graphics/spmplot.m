function [H,AX,BigAx] = spmplot(Y,group,plo,dispopt)
%spmplot: scatterplot matrix with dynamic scatters and boxplots or histograms on the main diagonal
%
%<a href="matlab: docsearch('spmplot')">Link to the help function</a>
%
%  Required input arguments:
%     Y : Data matrix containing n observations on v variables.
%         Rows of Y represent observations, and columns represent variables.
%  group: vector with n elements, grouping variable that determines the
%         marker and color assigned to each point. It can be a categorical
%         variable, vector, string matrix, or cell array of strings.
%
%  Optional input arguments:
%
%    plo: spmplot(Y,group,plo) enables to specify the names which are
%         displayed in the margins of the scatter-plot matrix and the
%         labels of the legend. plo can be an empty value [], a scalar,
%         or a structure.
%         If plo is set to the empty vector [], then nameY and labeladd are
%         both set to the empty string '' (default), and no labels nor
%         names are added to the plot.
%         If plo = 1 the names Y1,..., Yv are added to the margins of the
%         the scatter plot matrix else nothing is added.
%         If plo is a structure it may contain the following fields:
%         - labeladd: if this option is '1', the elements belonging to the
%                max(group) in the spm are labelled with their unit row index.
%                The default value is labeladd = '', i.e. no label is added.
%         - nameY: cell array of strings containing the labels of the
%                variables. As default value, the labels which are added
%                are Y1, ..., Yv.
%         - clr: a string of color specifications. By default, the colors
%                are 'brkmgcy'.
%         - sym: a string or a cell of marker specifications. For example,
%                if sym = 'o+x', the first group will be plotted with a
%                circle, the second with a plus, and the third with a 'x'.
%                This is obtained with the assignment plo.sym = 'o+x'
%                or equivalently with plo.sym = {'o' '+' 'x'}.
%                By default the sequence of marker types is:
%                '+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'
%         - siz: scalar, a marker size to use for all plots. By default the
%                marker size depends on the number of plots and the size of
%                the figure window. Default is siz = '' (empty value).
%       - doleg: a string to control whether legends are created or not.
%                Set doleg to 'on' (default) or 'off'.
%
%   dispopt: spmplot(Y,group,[],dispopt) enables to add on the main diagonal
%       of the scatter plot matrix stacked histograms, or univariate boxplots
%       for each of the groups.
%       dispopt is a string which lets you control how to fill the diagonals
%       in a plot of Y vs Y. Set dispopt to 'hist' (default) to plot histograms,
%       or 'box' to plot boxplots.
%       Remark: to set dispopt without changing the defaults for plo use, e.g.,
%          spmplot(Y,group,[],'box');
%          REMARK: The style which is used for univariate boxplots is
%          'traditional' if the number of groups is <=5, else it is compact
%
%  Output:
%
%   spmplot has the same output of gplotmatrix in the statistics toolbox:
%   [H,AX,BigAx] = spmplot(...) returns an array of handles H to the
%   plotted points; a matrix AX of handles to the individual subaxes; and a
%   handle BIGAX to big (invisible) axes framing the subaxes.  The third
%   dimension of H corresponds to groups in G. AX contains one extra row of
%   handles to invisible axes in which the histograms are plotted. BigAx is
%   left as the CurrentAxes so that a subsequent TITLE, XLABEL, or YLABEL
%   will be centered with respect to the matrix of axes.
%
%
% Copyright 2008-2014.
% Written by FSDA team
%
%<a href="matlab: docsearch('spmplot')">Link to the help function</a>
% Last modified 08-Dec-2013

% Examples:

%
%{

    % Generate contaminated data
    state=100;
    randn('state', state);
    n=200;
    Y=randn(n,3);
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;

    % spmplot is called automatically by all outlier detection methods, e.g. FSM
    [out]=FSM(Ycont,'plots',1);

%}

%{

    % Now test the direct use of FSM. Set two groups, e.g. those obtained
    % from FSM
    % Generate contaminated data
    state=100;
    randn('state', state);
    n=200;
    Y=randn(n,3);
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
   
    close all;
    [out]=FSM(Ycont,'plots',1);

    group = zeros(n,1);
    group(out.outliers)=1;
    plo=struct; plo.labeladd='1'; % option plo is used to label the outliers

    % By default, the legend identifies the groups with 'Group 1', Group 2', etc.
    spmplot(Ycont,group,plo,'box');

%}

%{

    % With two groups, and if the tag figure includes the word 'outlier',
    % the legend will identify one group for outliers and the other for
    % normal units
    figure('tag','This is an outlier scatterplot');
    spmplot(Ycont,group,plo);

    % If the tag figure includes the word 'brush', the legend will identify
    % one group for 'Unbrushed units' and the others for 'Brushed units 1',
    % 'Brushed units 2', etc.
    figure('Tag','Scatterplot with brushed units');
    spmplot(Ycont,group,plo);

%}

%{
    % Iris data: scatter plot matrix with univariate boxplots on the main
    % diagonal
    load fisheriris;
    plo=struct;
    plo.nameY={'SL','SW','PL','PW'}
    spmplot(meas,species,plo,'box');

%}


%{
    % Iris data: scatter plot matrix with univariate boxplots on the main
    % diagonal and personalized options for symbols, colors and symbol
    % size and no legend
    load fisheriris;
    plo=struct;
    plo.nameY={'SL','SW','PL','PW'}; % Name of the variables
    plo.clr='kbr'; % Colors of the groups
    plo.sym={'+' '+' 'v'}; % Symbols of the groups (inside a cell)
    % Symbols can also be specified as characters
    % plo.sym='++v'; % Symbols of the groups
    plo.siz=3.4; % Symbol size
    plo.doleg='off'; % Remove the legend
    spmplot(meas,species,plo,'box');
%}

%{
% An example with 5 groups
    n1=100;
    n2=80;
    n3=50;
    n4=80;
    n5=70;
    v=5;
    Y1=randn(n1,v)+5;
    Y2=randn(n2,v)+3;
    Y3=rand(n3,v)-2;
    Y4=rand(n4,v)+2;
    Y5=rand(n5,v);

    group=ones(n1+n2+n3+n4+n5,1);
    group(n1+1:n1+n2)=2;
    group(n1+n2+1:n1+n2+n3)=3;
    group(n1+n2+n3+1:n1+n2+n3+n4)=4;
    group(n1+n2+n3+n4+1:n1+n2+n3+n4+n5)=5;


    Y=[Y1;Y2;Y3;Y4;Y5];
    spmplot(Y,group,[],'box');
%}



%% Beginning of code
ngroups = numel(unique(group,'stable'));

[n,v]=size(Y);
if nargin<3;
    plo=1;
end

if nargin<4;
    dispopt='hist';
end

% Specify default values for colors, symbols, size of symbols and precence
% of legend
clrdef='brkmgcybrkmgcybrkmgcybrkmgcybrkmgcybrkmgcy';
symdef={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.';'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
sizdef='';
dolegdef='on';

if isstruct(plo)
    fplo=fieldnames(plo);
    
    d=find(strcmp('nameY',fplo));
    if d>0
        nameY=plo.nameY;
    else
        nameY=cellstr(num2str((1:v)','Y%d'));
    end
    d=find(strcmp('labeladd',fplo));
    if d>0
        labeladd=plo.labeladd;
    else
        labeladd='';
    end
    d=find(strcmp('clr',fplo));
    if d>0
        clr=plo.clr;
        
        if length(clr) ~= ngroups
            warning('Number of colors which have been supplied is not equal to the number of groups')
            disp(['Number of groups =' num2str(ngroups)])
            disp(['Number of colors =' num2str(length(clr))])
            if length(clr)< ngroups
                disp('Supplied colors will be duplicated')
                if isrow(clr)
                    clr=repmat(clr,1,ngroups);
                else
                    clr=repmat(clr,ngroups,1);
                end
            else
                disp(['Just the first ' num2str(ngroups) ' colors will be used'])
            end
        end
        
    else
        clr=clrdef;
    end
    d=find(strcmp('sym',fplo));
    if isempty(d)
        d=0;
    end
    if d>0 % && (ngroups == numel(plo.sym))
        sym=plo.sym;
        if length(sym) ~= ngroups
            warning('Number of symbols which have been supplied is not equal to the number of groups')
            disp(['Number of groups =' num2str(ngroups)])
            disp(['Number of symbols =' num2str(length(sym))])
            if length(sym)< ngroups
                disp('Supplied symbols will be duplicated')
                if isrow(sym)
                    sym=repmat(sym,1,ngroups);
                else
                    sym=repmat(sym,ngroups,1);
                end
            else
                disp(['Just the first ' num2str(ngroups) ' symbols will be used'])
            end
        end
    else
        sym=symdef;
    end
    d=find(strcmp('siz',fplo));
    if d>0
        siz=plo.siz;
    else
        siz=sizdef;
    end
    d=find(strcmp('doleg',fplo));
    if d>0
        doleg=plo.doleg;
    else
        doleg='on';
    end
    
    
else
    if ischar(plo)
        error('FSDA: Third argument must be a structure, or a scalar or an empty value []')
    end
    
    if plo==1
        nameY = cellstr(num2str((1:v)','Y%d'));
    else
        nameY = '';
    end
    labeladd='';
    
    clr=clrdef;
    sym=symdef;
    siz=sizdef;
    doleg=dolegdef;
end

if iscell(group)
    groupv = zeros(numel(group),1);
    guni = unique(group,'stable');
    for ii=1:numel(guni)
        groupv(strcmp(group,guni(ii))) = ii;
    end
else
    groupv = group;
end

unigroup = 1:ngroups;




%% The scatterplot matrix with histograms or boxplots (on the main diagonal) generalised to groups

% sym can be either a cell array or a character
if iscell(sym)
    charsym=char(sym{unigroup});
else
    charsym=sym(unigroup);
end
[H,AX,BigAx] = gplotmatrix(Y,[],group,clr(unigroup),charsym,siz,doleg,'hist',nameY,nameY);

% The third dimension of H distinguishes the groups. If there are no groups
% then ndims(H) = 2.
if ndims(H) == 3
    for i=1:size(AX,2)
        hold('on');
        
        if strcmp(dispopt,'hist')==1
            % Add the histograms generalised to groups
            
            ax = AX(i,i);
            Xlim=get(ax,'Xlim');
            
            ax = AX(end,i);
            Ylim=get(ax,'Ylim');
            
            %set(gcf,'CurrentAxes',ax);
            
            % the strings used to label the tick marks
            XTickLabel = get(ax,'XTickLabel');
            YTickLabel = get(ax,'YTickLabel');
            XLabel = get(get(ax,'XLabel'),'String');
            YLabel = get(get(ax,'YLabel'),'String');
            
            [~, ~] = histFS(Y(:,i),10,groupv,'',ax,clr(unigroup)); %'br'
            % Prevent from changing the limits when the figure is resized:
            % 1.Freeze the current limits
            set(ax,'XLimMode','manual','YLimMode','manual');
            set(ax,'xlim',Xlim)
            set(ax,'ylim',Ylim)
            % 2.Freeze the current tick values
            set(ax,'XTickMode','manual','YTickMode','manual');
            %Now restore the labels of the gplotmatrix
            set(ax,'XTickLabel',XTickLabel,'YTickLabel',YTickLabel);
            set(get(ax,'XLabel'),'String',XLabel);
            set(get(ax,'YLabel'),'String',YLabel);
            
        else
            % Add the boxplots generalised to groups
            
            ax = AX(end,i);
            axPosition = get(ax,'position');
            
            if length(unigroup) <= 5
                plotstyle = 'traditional';
            else
                plotstyle = 'compact';
            end
            
            hbp = boxplot(ax,Y(:,i),groupv,'plotstyle',plotstyle,'colors',clr(unigroup),'labelverbosity','minor','symbol','+');
            
            % Use boxplot without grouping variable (less efficient):
            % Y_n=[];
            % for groups  = 1:length(unigroup)
            %     Y_tmp   = Y(groupv == groups-1 , i);
            %     Y_tmp   = [Y_tmp ; NaN(length(Y)-length(Y_tmp),1)];
            %     Y_n     = [Y_n Y_tmp];
            % end
            % boxplot(ax, Y_n, ...);
            
            % Remove the x tick labels from the graph containing boxplots
            set(ax,'XTickLabel',{' '});
            
            % Remove the y tick labels from the graph containing boxplots
            set(ax,'YTickLabel',{' '});
            
            % Put the graph containing boxplots in the correct position
            set(ax,'position',axPosition);
            
            % Adjust the vertical scale (using the y scale of a scatter in
            % the same row of the scatter plot matrix)
            if i < size(Y,2);
                ax1 = AX(i,end);
                Ylim=get(ax1,'Ylim');
            else
                ax1 = AX(i,1);
                Ylim=get(ax1,'Ylim');
            end
            set(ax,'Ylim',Ylim)
            
            % The tag is set for later use in add2spm by
            % clickableMultiLegend
            for gg=1:numel(unigroup)
                set(hbp(:,gg),'Tag',['boxplot' num2str(gg)]);
            end
            
        end
    end
    
    % Put in the figure UserData field the list of units in the last group,
    % i.e. (depending on context) the outliers or the units last brushed.
    seq = (1:n)';
    if isnumeric(group)
        set(H(:,:,end), 'UserData' , seq(group==max(group)));
    else
        set(H(:,:,end), 'UserData' , seq(groupv==max(groupv)));
    end
end

if strcmp(doleg,'on')
    % Add to the spm the clickable multilegend and eventually the text labels
    % of the selections
    if isnumeric(group)
        % use context sensitive legends
        add2spm(H,AX,BigAx,'labeladd',labeladd,'userleg','1');
    else
        % use legends in guni
        add2spm(H,AX,BigAx,'labeladd',labeladd,'userleg',guni);
    end
end

end