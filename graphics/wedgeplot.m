function hf = wedgeplot(RES,varargin)
%wedgeplot generates the double wedge plot of a time series
%
%<a href="matlab: docsearchFS('wedgeplot')">Link to the help function</a>
%
%  Required input arguments:
%
%        RES  :  absolute scaled residuals. Matrix or structure.
%                Matrix of size T-by-(T-lshift) containing scaled residuals
%                (in absolute value) for all the T units of the original time
%                series monitored in steps lshift+1, lshift+2, ..., T-lshift,
%                where lshift+1 is the first tentative level shift position,
%                lshift +2 is the second level shift position, and so on. This
%                matrix can be created by funtion LTSts (Least Trimmed Squares
%                in time series). If RES is a structure, it must contain field:
%                RES.RES = matrix containing scaled residuals.
%                Data Types - double.
%
%       Optional input arguments:
%
%     transpose: option determining the posiiton of the index number or tentative
%                level shift. Boolean. If transpose is true (default) the x-axis
%                contains the tentative level shift position and the y-axis the
%                index number else if it is false the axes are interchanged.
%                When transpose is true, it is possible with option extradata to
%                add on a separate panel a subplot of the original time series
%                (and possibly the series of fitted values). See extradata
%                option for details.
%                Example - 'transpose',false
%                Data Types - Boolean
%
%   extradata :  extra data to plot in a separate panel in association to the
%                wedge plot. Matrix. Matrix of size T-by-1 or T-by-p containing
%                the data which have to be plotted in the separate panel.
%                Generally extradata is a matrix of size T-by-2 containing the
%                original time series and the corresponding fitted values in
%                order to link the irregularities shown by the wedgeplot with
%                the original time series.
%                - If extradata is empty (default) the double wedge plot will be
%                  shown in a single panel.
%                - If extradata is not empty a two panel plot will be created:
%                  one will contain the double wedge plot and extradata will be
%                  plot in the other panel. This options makes sense only if
%                  transpose is true, that is if the x axis of the double wedge
%                  plot contains the index number.
%                When option transpose is left by the user unspecified, the
%                default position of the extradata subplot is at the bottom.
%                Otherwise, the position of the two panels depends on the order
%                with which the user specifies the two options: if extradata is
%                specified first, the corresponding subplot will be at the top,
%                otherwse it will fall at the bottom.
%                Example - 'extradata', [y yhat]
%                Data Types - double
%
%      cmapname: color map. Character. Character which indicates the type of
%                colormmap to use in the wedge plot. The accepted values are
%                'hot', 'autumn', 'spring', 'pink', 'summer', 'winter', 'gray'.
%                The default is 'hot'.
%                Example - 'cmapname','summer'
%                Data Types - Character
%
%        labls : label of the axis which contains the level shift position.
%                Character. Character containing the label to put on the axis
%                which contains the level shift position. This axis could be
%                either the horizontal or vertical depending on the option
%                transpose. The default label is 'Tentative level shift
%                position'.
%                Example - 'labls','Position of level shift'
%                Data Types - Character
%
%        labin : label of the axis which contains the index number.
%                Character. Character containing the label to put on the axis
%                which contains the index number of the units of the time
%                series. This axis could be either the horizontal or vertical
%                depending on the option transpose. The default label is 'Index
%                number'.
%                Example - 'labin','unit number'
%                Types - Character
%
%     titl     : Title. String. A label for the title 
%                (default: 'Double wedge plot').
%                Example - 'titl','Plot with two wedges'
%                Data Types - char
%
%     FontSize:  Font size of the labels. Scalar. Scalar which controls the font
%                size of the labels of the axes and of the labels inside the
%                plot. Default value is 12.
%                Example - 'FontSize',12
%                Data Types - double
%
%  SizeAxesNum:  Size of the numbers of the axis. Scalar. Scalar which controls
%                the size of the numbers of the axes. Default value is 12.
%                Example - 'SizeAxesNum',10
%                Data Types - double
%
% Output:
%
%       hf  : handle to the figure. Graphics handle. Handle to the figure
%               which has just been created.
%
% See also: plot
%
% References:
%
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [RPRH]
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('wedgeplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Double wedge plot with simulated data with linear trend and level shift.
    % No seasonal component.
    pwd
    n=45;
    a=1;
    b=0.8;
    sig=1;
    seq=(1:n)';
    y=a+b*seq+sig*randn(n,1);
    y(round(n/2):end)=y(round(n/2):end)+10;
    % model with a quadratic trend, non seasonal and level shift
    model=struct;
    model.trend=2;
    model.seasonal=0;
    % Potential level shift position is investigated in positions:
    % t=10, t=11, ..., t=T-10.
    model.lshift=10:n-1;
    out=LTSts(y,'model',model);
    wedgeplot(out,'transpose',true,'extradata',[y out.yhat]);  
%}

%{
    % Example of double wedge plot in series with level shift.
    % Analysis of contaminated airline data.
    % Load the airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
        118  126  150  180  196  188  233  277  301  318  342  391    % Feb
        132  141  178  193  236  235  267  317  356  362  406  419    % Mar
        129  135  163  181  235  227  269  313  348  348  396  461    % Apr
        121  125  172  183  229  234  270  318  355  363  420  472    % May
        135  149  178  218  243  264  315  374  422  435  472  535    % Jun
        148  170  199  230  264  302  364  413  465  491  548  622    % Jul
        148  170  199  242  272  293  347  405  467  505  559  606    % Aug
        136  158  184  209  237  259  312  355  404  404  463  508    % Sep
        119  133  162  191  211  229  274  306  347  359  407  461    % Oct
        104  114  146  172  180  203  237  271  305  310  362  390    % Nov
        118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl
    y=(y(:));
    % Add a level shift contamintion plus some outliers.
    y(68:end)=y(68:end)+1300;
    y(67)=y(67)-600;
    y(45)=y(45)-800;
    y(68:69)=y(68:69)+800;
    % Create structure specifying model
    model=struct;
    model.trend=2;              % quadratic trend
    model.s=12;                 % monthly time series
    model.seasonal=204;         % number of harmonics
    model.lshift=40:120;        % position where to start monitoring level shift
    model.X='';
    % Create structure lts specifying lts options
    lts=struct;
    lts.bestr=20; % number of best solutions to bring to full convergence
    % h = dimension of the h subset (75 per cent of the data, bdp=0.25)
    h=round(0.75*length(y));
    [out, varargout]=LTSts(y,'model',model,'nsamp',500,...
        'lts',lts,'h',h,'plots',0,'msg',1);
    % Create the double wedge plot.
    wedgeplot(out);
%}

%{
    %% Example of double wedge plot in series with level shift with option transpose.
    % Analysis of contaminated airline data.
    % Load the airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960.
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
        118  126  150  180  196  188  233  277  301  318  342  391    % Feb
        132  141  178  193  236  235  267  317  356  362  406  419    % Mar
        129  135  163  181  235  227  269  313  348  348  396  461    % Apr
        121  125  172  183  229  234  270  318  355  363  420  472    % May
        135  149  178  218  243  264  315  374  422  435  472  535    % Jun
        148  170  199  230  264  302  364  413  465  491  548  622    % Jul
        148  170  199  242  272  293  347  405  467  505  559  606    % Aug
        136  158  184  209  237  259  312  355  404  404  463  508    % Sep
        119  133  162  191  211  229  274  306  347  359  407  461    % Oct
        104  114  146  172  180  203  237  271  305  310  362  390    % Nov
        118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl
    y=(y(:));
    % Add a level shift contamintion plus some outliers.
    y(50:55)=y(50:55)-300;
    y(68:end)=y(68:end)-700;
    y(70:75)=y(70:75)+300;
    y(90:90)=y(90:90)+300;
    % Create structure specifying model
    model=struct;
    model.trend=2;              % quadratic trend
    model.s=12;                 % monthly time series
    model.seasonal=204;         % number of harmonics
    model.lshift=40:120;        % position where to start monitoring level shift
    model.X='';
    % Create structure lts specifying lts options
    lts=struct;
    lts.bestr=20; % number of best solutions to bring to full convergence
    % h = dimension of the h subset (75 per cent of the data, bdp=0.25)
    [out, varargout]=LTSts(y,'model',model,'nsamp',500,...
        'lts',lts,'plots',0,'msg',1);
    % Create the double wedge plot.
    % Remember to remove the last column of the matrix of the residuals
    % obtained for each level shift position if you want to avoid the
    % top orange band (just execute RES(:,64)=[] before line 258).
    wedgeplot(out,'transpose',true,'extradata',[y out.yhat]);
%}

%{
    % Same double wedge plot as before, but with the time series at the top
    % subplot. This is obtained simply by specifying extradata before transpose.
    % Load the airline data.
    %   1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960.
    y = [112  115  145  171  196  204  242  284  315  340  360  417    % Jan
        118  126  150  180  196  188  233  277  301  318  342  391    % Feb
        132  141  178  193  236  235  267  317  356  362  406  419    % Mar
        129  135  163  181  235  227  269  313  348  348  396  461    % Apr
        121  125  172  183  229  234  270  318  355  363  420  472    % May
        135  149  178  218  243  264  315  374  422  435  472  535    % Jun
        148  170  199  230  264  302  364  413  465  491  548  622    % Jul
        148  170  199  242  272  293  347  405  467  505  559  606    % Aug
        136  158  184  209  237  259  312  355  404  404  463  508    % Sep
        119  133  162  191  211  229  274  306  347  359  407  461    % Oct
        104  114  146  172  180  203  237  271  305  310  362  390    % Nov
        118  140  166  194  201  229  278  306  336  337  405  432 ]; % Dec
    % Source:
    % http://datamarket.com/data/list/?q=provider:tsdl
    y=(y(:));
    % Add a level shift contamintion plus some outliers.
    y(50:55)=y(50:55)-300;
    y(68:end)=y(68:end)-700;
    y(70:75)=y(70:75)+300;
    y(90:90)=y(90:90)+300;
    % Create structure specifying model
    model=struct;
    model.trend=2;              % quadratic trend
    model.s=12;                 % monthly time series
    model.seasonal=204;         % number of harmonics
    model.lshift=40:120;        % position where to start monitoring level shift
    model.X='';
    % Create structure lts specifying lts options
    lts=struct;
    lts.bestr=20; % number of best solutions to bring to full convergence
    % h = dimension of the h subset (75 per cent of the data, bdp=0.25)
    [out, varargout]=LTSts(y,'model',model,'nsamp',500,...
        'lts',lts,'plots',0,'msg',1);
    % Create the double wedge plot.
    % Remember to remove the last column of the matrix of the residuals
    % obtained for each level shift position if you want to avoid the
    % top orange band (just execute RES(:,64)=[] before line 258).
    wedgeplot(out,'extradata',[y out.yhat],'transpose',true);
%}

%% Beginning of code

options=struct('extradata',[],'cmapname','hot',...
    'labls','Level shift position','labin','Index number',...
    'titl','Double wedge plot',...
    'FontSize',14,'SizeAxesNum',14,'transpose',true);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:wedgeplot:WrongInputOpt',...
            'Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

[extradata, cmapname, labls, labin, titl, FontSize, SizeAxesNum, transpose] = deal(...
    options.extradata, options.cmapname, options.labls, options.labin, ...
    options.titl, options.FontSize, options.SizeAxesNum, options.transpose);

%% Initialise key structures

% Check if input is a structure
if isstruct(RES)
    if isfield(RES,'posLS')
        posLS = RES.posLS;
    else
        posLS = [];
    end
    if isfield(RES,'RES')
        residuals = RES.residuals;
    else
        residuals = [];
    end
    if isfield(RES,'RES')
        RES=abs(RES.RES);
    else
        error('FSDA:wedgeplot:WrongInput','Input structure must contain a field named ''RES''')
    end
else
    % Take absolute values of RES
    RES = abs(RES);
    residuals = [];
    posLS = [];
end
% ADcont2 data; to avoid orange band, execute      RES(:,64)=[]
[T, l]  = size(RES);

% LSH = vector of integers associated with tentative level shift positions
lshift = (T-l)/2;
LSH    = (lshift+1):(T-lshift);

%% colormap
scmap = T*l;
switch cmapname
    case 'hot'
        cmap = hot(scmap);
    case 'autumn'
        cmap = autumn(scmap);
    case 'spring'
        cmap = spring(scmap);
    case 'pink'
        cmap = pink(scmap);
    case 'summer'
        cmap = summer(scmap);
    case 'winter'
        cmap = winter(scmap);
    case 'gray'
        cmap = gray(scmap);
    otherwise
        error('FSDA:wedgeplot:WrongInputOpt','The colormap name supplied is invalid.');
end

%Small values should go with light colors
cmap = flipud(cmap);

% add a section to the colormap to obtain a grey background (gray colormaps are
% unaffected, of course)
if ~strcmp(cmapname,  'gray')
    gray_levels = 50;
    gcmap = flipud(gray(gray_levels));
    gray_levels = 10;
    gcmap = gcmap(1:gray_levels,:);
    cmap = [gcmap ; cmap];
end

%% initialise figure with a colorbar

hf   = figure;
colormap(cmap);

%% colors that will be used for the wedge plot

% The colors of the wedge plot will be proportional to the scaled residuals.
% Very large values (> 50) are set to 50.
Cres=RES;
Cres(Cres>50)=50;

% a linear sequence of colormap weights from the minimum to the maximum
% residuals values (0 to 50, in practice).
scres_lin = (min(Cres(:)):max(Cres(:)))+1; scres_lin(1) = 0;

%Partition the residuals in three parts: 'large', 'medium' and 'small'
thtmp=50;
rlarge   = find(Cres>thtmp);            %#ok<NASGU>  may be used in future releases, e.g. Cres(rlarge)=thtmp;
rmedium  = find(Cres>2.58 & Cres<=10);  %#ok<NASGU>  may be used in future releases
rsmall   = Cres<=2.58;

%The colors are rescaled so that to map the elements of the colormap
Cres = Cres ./ max(max(Cres));
Cres = Cres * size(cmap,1);
Cres = round(Cres);
%Colors that point to 0 must be moved to 1
Cres(Cres==0) = 1;

% Control backgroud
if ~strcmp(cmapname,  'gray')
    %associate the small residuals to a gray color
    %for no gray,     set Cres(rsmall) to 2
    %for faint gray,  set Cres(rsmall) to floor(gray_levels/2)
    %for strong gray, set Cres(rsmall) to gray_levels
    Cres(rsmall)  = 1;
end
% if 0
%     %this has to be discussed
%     Cres(rmedium) = round(Cres(rmedium) *3);
% end

%% generates the double wedge plot

% In MATLAB releases before R2012b properties of surface and colorbar
% objects were different
vlt15 = verLessThan('matlab', '7.15');

Cres = [Cres; nan(1,length(LSH))];
Cres = [Cres nan(T+1,1)];

if transpose == false
    
    % the surface of the wedgeplot
    surface(zeros(size(Cres)),Cres,...
        'EdgeColor','none','Xdata',[LSH nan]','CDataMapping','direct');
    % axes labels
    xlabel(labls,'Fontsize',FontSize,'interpreter','none');
    ylabel(labin,'Fontsize',FontSize,'interpreter','none');
    
    % Colorbar and properties of the surface axes. Note the -1 and +1 in
    % the Ylim settings, i.e. in min(LSH-1) and max(LSH+1); this is needed
    % because the line of the Box would be covered by the surface of the
    % wedgeplot, at least in the bottom part of the plot.
    % The if statement addresses properties in different MATLAB releases.
    if  ~vlt15
        set(gca,'Box','on','Boxstyle','full','LineWidth',1,...
            'Xlim',[min(LSH-1), max(LSH+1)],'Fontsize',SizeAxesNum);
        
        colorbar('Ticks' , prctile(1:size(cmap,1),[1 20 40 60 80 100]),...
            'TickLabels' , round(prctile(scres_lin , [1 20 40 60 80 100])*100)/100,...
            'Fontsize',FontSize);
    else
        set(gca,'Box','on','LineWidth',1,...
            'Xlim',[min(LSH-1), max(LSH+1)],'Fontsize',SizeAxesNum);
        
        colorbar('YTick' , prctile(1:size(cmap,1),[1 20 40 60 80 100]),...
            'YTickLabels' , round(prctile(scres_lin , [1 20 40 60 80 100])*100)/100,...
            'Fontsize',FontSize);
    end
    
    title(titl,'interpreter','none','FontSize',FontSize+2);
    
else
    
    % figure formed by two panels: the wedge plot and the time series
    if ~isempty(extradata)
        
        trapos = find(strcmpi('transpose',UserOptions));
        extpos = find(strcmpi('extradata',UserOptions));
        if extpos > trapos
            dps = 2; wps = 1;
        else
            dps = 1; wps = 2;
        end
        % subplot hosting the wedgeplot
        A(wps) = subplot(2,1,wps);
        
    else
        xlabel(labin,'Fontsize',FontSize,'interpreter','none');
    end
    
    % the surface of the wedgeplot
    surface(zeros(size(Cres))',Cres',...
        'EdgeColor','none','Ydata',[LSH nan]','CDataMapping','direct');
    % axes labels
    % xlabel(labin,'Fontsize',FontSize);
    ylabel(labls,'FontSize',FontSize,'interpreter','none');
    
    % Colorbar and properties of the surface axes. Note the -1 and +1 in
    % the Ylim settings, i.e. in min(LSH-1) and max(LSH+1); this is needed
    % because the line of the Box would be covered by the surface of the
    % wedgeplot, at least in the bottom part of the plot.
    % The if statement addresses properties in different MATLAB releases.
    if ~vlt15
        set(gca,'Box','on','BoxStyle','full','LineWidth',1,...
            'Ylim',[min(LSH-1), max(LSH+1)],'Fontsize',SizeAxesNum);
        
        colorbar('eastoutside','Ticks' , prctile(1:size(cmap,1),[1 20 40 60 80 100]),...
            'TickLabels' , round(prctile(scres_lin , [1 20 40 60 80 100])*100)/100,...
            'Fontsize',FontSize-2);
    else
        set(gca,'Box','on','LineWidth',1,...
            'Ylim',[min(LSH-1), max(LSH+1)],'Fontsize',SizeAxesNum);
        
        colorbar('eastoutside','YTick' , prctile(1:size(cmap,1),[1 20 40 60 80 100]),...
            'YTickLabels' , round(prctile(scres_lin , [1 20 40 60 80 100])*100)/100,...
            'Fontsize',FontSize-2);
    end
    
    % the subplots have to be rescaled for leaving space to the colorbar
    if ~isempty(extradata)
        
        mine = min(extradata(:));
        maxe = max(extradata(:));
        delta = (maxe-mine)*0.1;
        yaxlim = [mine - delta ; maxe + delta];
        
        A(dps) = subplot(2,1,dps);
        
        hold('on');
        
        % mark outliers with their severity
        if ~isempty(residuals)
            seq = 1:T;
            quant = sqrt(chi2inv(1-0.01,1));
            boolean=abs(residuals)>quant;
            resboo=residuals(boolean);
            modres=resboo;
            th=8;
            modres(abs(resboo)>th)=th;
            %Rescale residuals in the interval [0 3]
            sizeout=3*(abs(modres)-quant)/(th-quant);
            outliers=seq(boolean);
            for i=1:length(sizeout)
                plot(seq(outliers(i)),extradata(outliers(i),1),'x','LineWidth',sizeout(i),'Color','r', 'MarkerFaceColor','k');
            end
        end
        
        % plot the vertical line of the level shift position and the
        % associated label on the X axis
        if ~isempty(posLS)
            line(posLS*ones(2,1) , yaxlim , 'LineStyle' , ':' , 'LineWidth' , 1.5 , 'Color' , 'k');
            text(posLS , yaxlim(1) , num2str(posLS) , 'HorizontalAlignment' , 'Center' , 'VerticalAlignment' ,  'Top');
        end
        
        % plot the time series
        dd  = size(extradata,2);
        clr = 'bkrgmcy';
        syb = {'-','--','-.',':','-','--','-.'};
        for d=1:dd
            plot(extradata(:,d),'Color',clr(d),'LineStyle',syb{d},'LineWidth',1);
        end
        
        xlabel(A(2),labin,'FontSize',FontSize,'interpreter','none');
        if ~vlt15
            set(gca,'FontSize',SizeAxesNum,'Ylim' , yaxlim,'Box','on','BoxStyle','full');
        else
            set(gca,'FontSize',SizeAxesNum,'Ylim' , yaxlim,'Box','on');
        end
        for i=1:2
            pos=get(A(i), 'Position');
            axes(A(i)) ; %#ok<LAXES>
            set(A(i), 'Position', [pos(1) pos(2) .6626 pos(4)]);
        end
        title(A(1),titl,'interpreter','none','FontSize',FontSize+2);
        box('on');
    else
        title(titl,'interpreter','none','FontSize',FontSize+2);
    end
    
end

% get the current axes position and increases the bottom margin to ensure
% that the x label is not cut off
position = get(gca,'OuterPosition');
[left, bottom, width, height] = deal(position(1),position(2),position(3),position(4));
set(gca,'OuterPosition',[left 0.7*bottom width height]);

end
%FScategory:VIS-Reg
