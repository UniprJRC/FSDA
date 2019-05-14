function cdsplot(outms,varargin)
%cdsplot produces the candlestick plot for robust model selection in linear regression
%
%<a href="matlab: docsearchFS('cdsplot')">Link to the help function</a>
%
% Required input arguments:
%
%  outms :  plot data. Structure. Structure containing the following fields
%        outms.stor = k x 9 matrix containing statistics which are used to create
%               the candles. 
%               1st col: max Cp values; 
%               2nd col: min Cp values; 
%               3rd col: averages Cp values; 
%               4nd col: median Cp values; 
%               5th col: x coordinates (or size of
%               submodel); 
%               6th col: number of explanatory variables of the submodel; 
%               7th col: y coordinate of final Cp; 
%               8th col: units entering the final step of the search; 
%               9th col: maximum Cp value during the search (This is
%               information is used to print the labels on top of each
%               model). 
%        outms.outl = r x 4 matrix containing information about 'influential
%               units' or empty matrix. 
%               Influential units in this context are defined as the units
%               which enter the subset in the final part of the search and
%               bring the value of the Cp below the minimum or
%               above the maximum value of the central part of the search. 
%               1st col: x coordinates; 
%               2nd col: y coordinates; 
%               3rd col: step of entry into subset; 
%               4nd col: unit number. 
%               If matrix outl contains more columns they are ignored. 
%         outms.siz = vector of length 2 containing information about n (number of
%               units of the sample and bigP, number of explanatory variables,
%               including the constant, in the full model). This
%               information is necessary to compute the envelopes.
%      Data Types - struct
%
% Optional input arguments:
%
% CandleWidth: width. Scalar. Scalar defining the width of the boxes which represent central part of the
%              candles. The default width is 0.05.
%                   Example - 'CandleWidth',0.1
%                   Data Types - double
%      color : Color. Vector. Three elements color vector, [R G B], or a string specifying the
%              color name. MATLAB supplies a default color if none is specified
%              or if it is empty. The default color differs depending on the
%              background color of the figure window. See COLORSPEC in the
%              MATLAB Reference Guide for color names.
%                   Example - 'color',[0.1 0.2 0.5]
%                   Data Types - double
%  LineWidth : Line Width. Scalar. Line Width (in points) for the vertical lines outside the boxes of the
%              candles. The default LineWidth is 0.5 points.
%                   Example - 'LineWidth',0.3
%                   Data Types - double
%    ylimy    : y axis scale. Vector. Vector with two elements controlling minimum and maximum
%              on the y axis. Default value is [-2 50] (automatic scale).
%                   Example - 'ylimy',[0 100]
%                   Data Types - double
%    xlimx    : x axis scale. Vector. Vector with two elements controlling minimum and maximum
%              on the x axis. Default value is '' (automatic scale).
%                   Example - 'xlimx',[0 100]
%                   Data Types - double
%   label   : Labels of the selected models. Cell array of strings. Cell array of strings of length k (number of rows of matrix stat)
%              containing the labels of the selected models. Default value is ''
%              that is no label is plotted on the screen.
%                   Example - 'label',{'a' 'b' 'c'}
%                   Data Types - char
%    quant   : Quantiles. Vector. Vector containing quantiles for the horizontal lines
%              associated witht the confidence bands of Cp. 
%              The default is to plot 2.5% and
%              97.5% envelopes. In other words the default is
%              quant=[0.025;0.975].
%                   Example - 'quant',[0.01 0.025 0.975 0.99]
%                   Data Types - double
%   lablast  : Label for the last unit entered. Scalar. Scalar which specifies whether to add the label of the unit
%              which enters the final step of the search close to its
%              symbol. If lablast=1 label is added else (default) no label
%              is added.
%                   Example - 'lablast',0
%                   Data Types - double
%   laboutl  : Label for the influential units. Scalar. Scalar which specifies whether to add the labels of the
%              'influential units'
%              if laboutl=1 the unit number is added close to its symbol.
%              if laboutl=2 the unit number together with the entry step is
%              added close to its symbol else (default) no label is added.
%                   Example - 'laboutl',1
%                   Data Types - double
%   labbold  : Models to highliht. Cell array of strings. Cell array of strings which specifies the models which have
%              to be highlighted (the linewidth of the vertical lines
%              outside the boxes of the models specified in labbold is
%              considerably increased). As default labbold=''.
%                   Example - 'labbold',{'a' 'b'}
%                   Data Types - char
%    labenv  :  Quantiles labels. Scalar. If labelv=1 labels of the quantiles used to generate the
%               horizontal lines associated with the envelopes are added,
%               else if labelv=0 (default) no label is added.  
%                   Example - 'labenv',1
%                   Data Types - double
%    barend  : Adding horizontal lines. Scalar. Scalar which specifies whether to add small horizontal lines
%              at the end of the vertical lines representing the whiskers.
%              If barend=1 horizontal lines are added else (default) no
%              additional line is drawn. 
%                   Example - 'barend',1
%                   Data Types - double
%    cpbrush : Brushing. Empty value or matrix. 
%              If cpbrush is an empty value (default), no brushing is done.
%              The activation of this option (cpbrush is a scalar) enables
%              the user  to select a set of candles in the candlestick plot
%              and to monitor the corresponding forward searches in a new
%              plot. 
%              If cpbrush is not an empty value the user has
%              to supply the matrix which in the first column contains the
%              fwd search index and in the other k columns the values of Cp
%              associated with the k models displayed in the candlestick
%              plot.
%                   Example - 'cpbrush',''
%                   Data Types - double
%
% Output:
%
% See also:
%
% References:
%
%   Riani M. and Atkinson A.C. (2010), Robust Model Selection with Flexible Trimming,
%   "Computational Statistics and Data Analysis", Vol. 54, p. 3300-3312.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('cdsplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % cdsplot with all default options.
    % Load Ozone data (reduced data)
    X=load('ozone.txt');
    % Tranform the response using logs
    X(:,end)=log(X(:,end));
    % Add a time trend
    X=[(-40:39)' X];
    % Define y
    y=X(:,end);
    % Define X
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
    % Robust model selection using Cp
    [Cpms]=FSRms(y,X,'labels',labels);
    % Candlestick plot
    cdsplot(Cpms);
%}


%{
    % Interactive_example
    % cdsplot with optional arguments.        
    % Load Ozone data (full data)
    X=load('ozone_330_obs.txt');
    y=log(X(:,9));
    Time1=[(1:165)';(165:-1:1)'];
    X=[Time1 X(:,1:8)];
    labels={'Time','1','2','3','4','5','6','7','8'};
    outms=FSRms(y,X,'labels',labels,'smallpint',5:6);
    cdsplot(outms,'cpbrush',1,'laboutl',1);
%}

%% Beginning of code

if nargin < 1
    error('FSDA:cdsplot:missingInputs', ...
        'Missing input structure containing required arguments to compute cdsplot.')
end

stat=outms.stor;

if size(stat,2)<8
    error('FSDA:cdsplot:missingInputs', ...
        'Missing matrix containing: max, min, averages, or medians of Cp values.')
end

% Retrieve matrix which contains information about the influential units
outl=outms.outl;

% Retrieve matrix which contains information about size of dataset
siz=outms.siz;

hi=stat(:,1);
lo=stat(:,2);
av=stat(:,3);
me=stat(:,4);
xcoord=stat(:,5);

% Extract the information about the number of units and the number of
% explanatory variables in the full model
n=siz(1);
p=siz(2);


CandleWidthdef=0.05;
options=struct('CandleWidth',CandleWidthdef,...
    'color','','barend',0,'cpbrush','',...
    'LineWidth',0.5,'lablast',0,'laboutl',0,'labbold','',...
    'xlimx','','ylimy','','label',1,'quant',[0.025 0.975],'labenv',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:tclust:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 1
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

% Determine if current plot is held or not
if ishold
    hldflag = 1;
else
    hldflag = 0;
end

width=options.CandleWidth;
color=options.color;
LineWidth=options.LineWidth;
xlimx=options.xlimx;
ylimy=options.ylimy;
% Labels of the models which will be put vertically on top of each candle
label=options.label;
% Confidence bands of the F distribution
quant=options.quant;
% If quant has been supplied as a column vector it is transposed
if size(quant,1)>1
    quant=quant';
end

%Information about adding label associated with the last unit which enters the
%search
lablast=options.lablast;
%Information about adding labels associated with the 'outliers'
laboutl=options.laboutl;

%Information about adding labels of the quantiles (horizontal bands)
labenv=options.labenv;

%Information about the models which have to be highlighted in the plot
labbold=options.labbold;

%Information about the horizontal lines at the end of the whiskers
barend=options.barend;

%Information about the possibility of brushing the candles
cpbrush=options.cpbrush;

index=xcoord';

if isempty(color)
    avs = get(gca, 'colororder');
    color = avs(1, :);
end


back = get(gca, 'color');



m = length(hi(:));

% Need to pad all inputs with NaN's to leave spaces between day data
tmp = nan;
nanpad = tmp(1, ones(1, m));
hilo = [hi'; lo'; nanpad];


indhilo = index(ones(3, 1), :);
plot(indhilo(:), hilo(:), 'color', color,'LineWidth',LineWidth)


avpad = [av(:)';nanpad];
avpad = avpad(:)';
mepad = [me(:)'; nanpad];
mepad = mepad(:)';

% Create boundaries for filled regions
xbottom = index - width;
xbotpad = [xbottom(:)'; nanpad];
xbotpad = xbotpad(:)';
xtop = index + width;
xtoppad = [xtop(:)'; nanpad];
xtoppad = xtoppad(:)';
ybottom = min(avpad, mepad);
ytop = max(avpad, mepad);

% Plot lines between high and low value of Cp for each model
hold on

% If the median is less than the average, box is empty
i = find(mepad(:) <= avpad(:));
boxes(i) = patch([xbotpad(i); xbotpad(i); xtoppad(i); xtoppad(i)],...
    [ytop(i); ybottom(i); ybottom(i); ytop(i)],...
    back, 'edgecolor', color);

% If the median price is greater than the average, box is filled
i = find(mepad(:) > avpad(:));
boxes(i) = patch([xbotpad(i); xbotpad(i); xtoppad(i); xtoppad(i)],...
    [ytop(i); ybottom(i); ybottom(i); ytop(i)],...
    color, 'edgecolor', color); %#ok


% Add the outliers
plot(outl(:,1),outl(:,2),'*r');


% Plot final value of Cp as a blue circle
% Fifth column of stat contains x coordinates
% Seventth column of stat contains final value of Cp
plot(stat(:,5),stat(:,7),'ob');


% Add vertical lines in order to separate the different values of p using black color
xl1=min(stat(:,6))-0.5;
xl2=max(stat(:,6))+0.5;
seql=xl1:xl2;
ons=ones(1,length(seql));
ycoord=vertcat(-100*ons,100*ons);
xcoord=repmat(seql,2,1);
% plot the vertical lines
line(xcoord,ycoord,'Color','k','tag','bar');

% Define limits for y axis
if isempty(ylimy)
    ylim([-2 60]);
else
    ylim(ylimy);
end

% Define limits for x axis
if isempty(xlimx)
    xlim([xl1 xl2]);
else
    xlim(xlimx);
end


if label
    label=outms.LAB;
    % Add labels associated with the different selected models
    % 9th col of matrix stor contains overall maximum of Cp values
    text(stat(:,5),stat(:,9)+0.5,label,'Rotation',90,'FontSize',14);
end

smallp=unique(stat(:,6));
Env=zeros(length(smallp)*3-1,length(quant)+1);
Env(:)=NaN;
% ienv = index which is linked to the rows of matrix Env
ienv=1;

for i=1:length(smallp)
    ast=finv(quant,p-smallp(i),n-p);
    ast=ast*(p-smallp(i))+2*smallp(i)-p;
    
    Env(ienv,:)=horzcat(smallp(i)-0.5,ast);
    Env(ienv+1,:)=horzcat(smallp(i)+0.5,ast);
    ienv=ienv+3;
end


% Superimpose bands (horizontal lines for each smallp) based on the F distribution
plot(Env(:,1),Env(:,2:end),'b','tag','bar');

% Write quantiles which are used to compute the envelopes
% The last element of envesel contains the index of matrix Env which
% specifies the y coordinates of the labels
if labenv
    envsel=find(Env(:,1)<=xl2);
    text(repmat(xl2+0.05,1,length(quant)),Env(envsel(end),2:end),num2str(quant'));
end

% Label for the last unit which entered the search
if lablast
    % Add the label of the units which entered the last step of the search
    text(stat(:,5)+0.03,stat(:,7),num2str(stat(:,8)));
end

% Label for the influential units
if laboutl==1
    % Add labels of the influential units
    text(outl(:,1)+0.03,outl(:,2),num2str(outl(:,4)));
elseif laboutl==2
    % Add labels of the influential units with the information about the
    % entry step
    aa=cellstr(num2str(outl(:,[4 3]),'%d,'));
    text(outl(:,1)+0.03,outl(:,2),aa);
else
    % No label is added
end

% Check if some models have to be highlighted
if ~isempty(labbold)
    [~, ~, ib] = intersect(labbold, label);
    % The strings contained in vector labbold are in rows ib of
    % cell array label
    % Extract from matrices indhilo and hilo the required candles
    % (vertical lines)
    indhilos=indhilo(:,ib);
    hilos=hilo(:,ib);
    plot(indhilos(:), hilos(:), 'color', 'k','LineWidth',LineWidth+5)
    
end

% if barend=1 add small horizontal lines at the end of the whiskers
if(barend)
    barendx=[(stat(:,5)-0.5*width) (stat(:,5)+0.5*width)];
    barendlowy=repmat(stat(:,2),1,2);
    barendupy=repmat(stat(:,1),1,2);
    line(barendx',barendlowy','color','b','tag','bar');
    line(barendx',barendupy','color','b','tag','bar');
else
    % no additional line is plotted
end
% If original figure was not held, turn hold off
if ~hldflag
    hold off
end

% Interactive brushing part
if ~isempty(cpbrush)
    storFWD=[(n-size(outms.MAL,1)+1:n)' outms.MAL];
    
    disp('Interactive part: select particular trajectories with the mouse');
    [~,xselect]=selectdataFS('selectionmode','Brush','brushsize',width,'brushshape','vrect','Ignore',findobj(gcf,'tag','bar'),'FlagColor','k');
    
    % From the x selected coordinates find the selected searches
    [~, ~, ib] = intersect(unique(cell2mat(xselect)), stat(:,5));
    
    % The selected fwd searches can be found in column ib of matrix storFWD
    figure;
    plot(storFWD(:,1),storFWD(:,ib+1));
    
    % Add labels for the brushed searches
    if ~isempty(label)
        % Add the proper forward envelopes based on the F distribution
        hold('on');
        text(repmat(n,length(ib),1),stat(ib,7),label(ib));
        % If the user has selected searches associated with only one small p then
        % the exact confidence bands are drawn
    end
    
    if max(stat(ib,6))==min(stat(ib,6))
        smallp=stat(ib(1),6);
        title(['Number of selected search =', num2str(length(ib)),'  small p=',num2str(smallp)]);
    else
        % draw confidence bands based on the maximum value of the small p
        % which are selected
        title(['Number of selected search =', num2str(length(ib)),'  Range of small p=',num2str(min(stat(ib,6))) '-' num2str(max(stat(ib,6)))]);
        smallp=max(stat(ib,6));
    end
    
    % Compute and superimpose envelopes based on the F distribution
    EnvF=[storFWD(:,1) zeros(size(storFWD,1),length(quant)+1)];
    
    for i=1:size(storFWD,1)
        ast=finv([quant 0.5],p-smallp,storFWD(i,1)-p);
        ast=ast*(p-smallp)+2*smallp-p;
        EnvF(i,2:end)=ast;
    end
    % Plot the envelopes
    line(EnvF(:,1),EnvF(:,2:end),'LineStyle','--','Color','r');
    
    % Write the labels for the quantiles which have been used to produce
    % the envelopes
    text(repmat(n,length(quant)+1,1),EnvF(end,2:end)',num2str([quant 0.5]'));
    
    
    xlabel('Subset size m')
    
end

%FScategory:VIS-Reg