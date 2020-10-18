function tclustICplotGPCM(IC,varargin)
%tclustICplotGPCM plots information criterion as a function of  $c_{det}$, $c_{shw}$,  $c_{shb}$ and $k$
%
%<a href="matlab: docsearchFS('tclustICplotGPCM')">Link to the help function</a>
%
%   tclustICplotGPCM takes as input the output of function tclustICgpcm
%   (that is a series of matrices which contain the values of the
%   information criteria BIC/ICL/CLA for different values of $k$,
%   $c_{det}$, $c_{shw}$ $c_{shb}$ and type or rotation The plot enables
%   interaction in the sense that, if option databrush has been activated,
%   it is possible to click on a point in the plot and to see the
%   associated classification in the scatter plot matrix.
%
%  Required input arguments:
%
%           IC : Information criterion to use. Structure.
%                It contains the following fields.
%                IC.MIXMIX  (or IC.MIXCLA or IC.CLACLA) = 3D array of size
%                   length(kk)-times-length(IC.ccdet)-times-length(IC.shw)
%                   containinig the value of the penalized classification
%                   likelihood.
%                IC.IDXMIX (or IC.IDXCLA) = cell of size
%                   length(kk)-times-length(pa.ccdet)-times-length(pa.ccshw).
%                   Each element of the cell is a vector of
%                   length n containinig the assignment of each unit using
%                   the classification model.
%                IC.kk = vector containing the values of k (number of
%                   components) which have been considered.
%                IC.ccdet = scalar or vector containing the values of
%                   $c_{det}$ (values of the restriction factor for the
%                   determinants) which have been considered.
%                IC.ccshw = scalar or vector containing the values of
%                   $c_{shw}$ (values of the restriction factor for the
%                   withing groups shape elements) which have been considered.
%                IC.alpha = scalar value of
%                   trimming which has been considered.
%                IC.BICbest = scalar containing optimal value of
%                   BIC (smallest value of BIC).
%                IC.nameY=  cell of length(size(Y,2)) containing the names
%                   of the variables of original matrix Y
%                 Data Types - struct
%
%  Optional input arguments:
%
%
%
%       tag     :   Personalized tag. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag 'pl_ICgpcm'.
%                   Note that if the program finds a plot which has a tag
%                   equal to the one specified by the user, then the output
%                   of the new plot overwrites the existing one in the same
%                   window else a new window is created.
%                   Example - 'tag','myplot'
%                   Data Types - char
%
%   datatooltip :   interactive clicking. Empty value (default) or
%                   structure. The default is datatooltip=''.
%                   If datatooltip = 1, the user can select with the
%                   mouse a solution in order to
%                   have the following information:
%                   1) value of $k$ which has been selected
%                   2) value of $c_{det}$ which has been selected
%                   3) value of $c_{shw}$ which has been selected
%                   4) value of $c_{shb}$ which has been selected
%                   5) values of the information criterion
%                   6) frequency distribution of the associated
%                   classification
%                   If datatooltip is a structure it may contain the
%                   following the fields
%                   datatooltip.DisplayStyle = Determines how the data
%                       cursor displays. datatip | window.
%                       - datatip displays data cursor
%                       information in a small yellow text box attached to a
%                       black square marker at a data point you interactively
%                       select.
%                       - window displays data cursor information for the
%                       data point you interactively select in a floating
%                       window within the figure.
%                   datatooltip.SnapToDataVertex=  Specifies whether the
%                       data cursor snaps to the nearest data value or is
%                       located at the actual pointer position.  on | off.
%                       - on data cursor snaps to the nearest data value
%                       - off data cursor is located at the actual pointer
%                       position.
%                   (see the MATLAB function datacursormode or the examples
%                   below). Default values are datatooltip.DisplayStyle =
%                   'Window' and datatooltip.SnapToDataVertex = 'on'.
%                   Example - 'datatooltip',''
%                   Data Types - scalar double or struct
%
%    databrush  :   interactive mouse brushing. Empty value, scalar or structure.
%                   If databrush is an empty value (default), no brushing
%                   is done. The activation of this option
%                   (databrush is a scalar or a structure) enables the user
%                   to select a set of values of IC in the current plot and
%                   to see the corresponding classification highlighted in
%                   the scatter plot matrix (spm). If spm does not exist it
%                   is automatically created. Please, note that the window
%                   style of the other figures is set equal to that which
%                   contains the IC plot. In other words, if the IC plot is
%                   docked all the other figures will be docked too.
%                   DATABRUSH IS A SCALAR. If databrush is a scalar the
%                   default selection tool is a rectangular brush and it is
%                   possible to brush only once (that is persist='').
%                   DATABRUSH IS A STRUCTURE. If databrush is a structure,
%                   it is possible to use all optional arguments of
%                   function selectdataFS.m and the following optional
%                   arguments: - databrush.persist = repeated brushing
%                   enabled. Persist is an empty value or a scalar
%                     containing the strings 'on' or 'off'.
%                     The default value of persist is '', that is brushing
%                     is allowed only once.
%                     If persist is 'on' or 'off' brushing can be done as
%                     many time as the user requires.
%                     If persist='on' then the unit(s) currently brushed
%                     are added to those previously brushed. it is
%                     possible, every time a new brushing is done, to use a
%                     different color for the brushed solutions.
%                     If persist='off' every time a new brush is performed
%                     units previously brushed are removed.
%                   - databrush.Label = add labels (i.e. x=value of k and
%                      y=values of IC) of brushed solutions in the ICplot.
%                     Character. [] (default) | '1'.
%                   - dispopt = string which controls how to fill the
%                      diagonals in the scatterplot matrix of the brushed
%                      solutions. Set dispopt to 'hist' (default) to plot
%                      histograms, or 'box' to plot boxplots.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%
%         nameY  : variable labels. Cell array. Cell array of strings
%                   containing the labels of the
%                   variables. As default value, the labels which are added
%                   are Y1, ..., Yv.
%                   Example - 'nameY',{'myY1', 'myY2'}
%                   Data Types - cell
%
%  Output:
%
%
%
% See also: tclustICplot, tclustICgpcm, tclust, tclustICsolGPCM
%
% References:
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustICplotGPCM')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Plot BIC, ICL and CLA for for Geyser data with all default options.
    Y=load('geyser2.txt');
    % Make sure (whenever possible) that units 15, 30 and 69 are inside
    % groups which have labels respectively equal to 1, 2 and 3.
    UnitsSameGroup=[15 30 69];
    nsamp=50;
    out=tclustICgpcm(Y,'cleanpool',false,'plots',0,'alpha',0.1,...
        'UnitsSameGroup',UnitsSameGroup,'nsamp',nsamp);
    tclustICplotGPCM(out);
%}

%{
    %   Example of the use of option datatooltip (all default options).
    %   Gives the user the possibility of clicking on the different points
    %   and have information about
    %   1) value of k which has been selected
    %   2) value of restriction factor which has been selected
    %   3) values of the information criterion
    %   4) frequency distribution of the associated classification
    Y=load('geyser2.txt');
    UnitsSameGroup=[15 30 69];
    nsamp=50;
    out=tclustICgpcm(Y,'cleanpool',false,'plots',0,'alpha',0.1,...
        'UnitsSameGroup',UnitsSameGroup,'nsamp',nsamp);
    tclustICplotGPCM(out,'datatooltip',1);
%}

%{
    % Example of the use of option datatooltip (personalized options).
    % Gives the user the possibility of clicking on the different points
    % and have information about the selected, the step of entry
    % into the subset and the associated label.
    datatooltip = struct;
    % In this example the style of the datatooltip is 'datatip'. Click on a
    % point when the ICplot is displayed.
    datatooltip.DisplayStyle = 'datatip';
    Y=load('geyser2.txt');
    UnitsSameGroup=[15 30 69];
    nsamp=50;
    out=tclustICgpcm(Y,'cleanpool',false,'plots',0,'alpha',0.1,...
        'UnitsSameGroup',UnitsSameGroup,'nsamp',nsamp);
    tclustICplotGPCM(out,'datatooltip',datatooltip);
%}


%{
    % Interactive_example
    % databrushing from any of the plot which was produced.
    % Use all default options for databrush (brush just once)
    close all
    Y=load('geyser2.txt');
    UnitsSameGroup=[15 30 69];
    out=tclustICgpcm(Y,'cleanpool',false,'plots',0,'alpha',0.1,...
        'UnitsSameGroup',UnitsSameGroup,'nsamp',nsamp);
    tclustICplotGPCM(out,'databrush',1);
%}

%{
    % Interactive_example
    % Repeated databrushing.
    % enable repeated brushing and show boxplots of groups inside diag of spm
    databrush=struct;
    % Set the shape of the brush
    databrush.selectionmode='Rect';
    % Enable repeated brushing
    databrush.persist='on';
    % Include x and y coordinates of brushed solutions inside ICplot
    databrush.Label='on';
    % Remove x and y coordinated just after btushing
    databrush.RemoveLabels='on';
    % show boxplots of the groups instead of histograms on the main
    % diagonal of the spm
    databrush.dispopt='box';
    Y=load('geyser2.txt');
    UnitsSameGroup=[15 30 69];
    out=tclustICgpcm(Y,'cleanpool',false,'plots',0,'alpha',0.1,...
        'UnitsSameGroup',UnitsSameGroup,'nsamp',nsamp);
    tclustICplotGPCM(out,'databrush',databrush);
%}



%% Beginning of code

% Specify information criterion for which you want to see best solutions
% no brushsing as default
databrush='';
datatooltip=1;
tag='';

if nargin>1
    options=struct('datatooltip',datatooltip, ...
        'tag',tag,'databrush', databrush,'nameY','');
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:tclustICplotGPCM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:tclustICplotGPCM:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    
    datatooltip=options.datatooltip;
    databrush=options.databrush;
    tag=options.tag;
end

cdet=IC.ccdet;
cshw=IC.ccshw;
BICbest=IC.BICbest;

if isfield(IC,'MIXMIX')
    selIC=IC.MIXMIX;
    nameselIC='MIXMIX';
    IDX=IC.IDXMIX;
elseif isfield(IC,'MIXCLA')
    selIC=IC.MIXCLA;
    nameselIC='MIXCLA';
    IDX=IC.IDXMIX;
elseif isfield(IC,'CLACLA')
    selIC=IC.CLACLA;
    nameselIC='CLACLA';
    IDX=IC.IDXCLA;
else
    error('FSDA:tclustICplotGPCM:NonExistInputOpt','Missing input');
end
[lkk,lcdet,lcshw]=size(selIC);

% Extract the values of k (number of groups)
kk=IC.kk;
% Extract the values of c (number of restriction factors)


[~,posmin]=min(selIC,[],'all','linear');

% bestk, bestcdet and best cshw are respectively best values for number of
% groups (bestk), restriction among determinants (bestcdet) and restriction
% among the shape elements inside each group (bestcshw)
[bestk,bestcdet,bestcshw]=ind2sub([lkk lcdet lcshw],posmin);

% Global variables to declare is datatooltip is not empty
if ~isempty(datatooltip)
    hTarget=[];
    hTargetlwd=[];
    hTargetcol=[];
    hTargetms=[];
end

%% Prepate the figure to display the ICplot
% Create a figure to host the plot or clear the existing one
% if isempty(tag)
%     h=figure;
% else
%     h=findobj('-depth',1,'tag',tag);
%     if (~isempty(h))
%         clf(h);
%         figure(h);
%         axes;
%     else
%         h=figure;
%     end
% end
% set(h,'Name', 'ICgpcm plot', 'NumberTitle', 'off');

% figure('Name','BIC')
% set line width of the trajectories of BIC
LineWidth=1;
% Define marker type
styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
lcc=max([lcdet,lcshw]);
styp=repmat(styp,ceil(200/lcc),1);
% Define line type
slintyp={'-';'--';':';'-.'};
slintyp=repmat(slintyp,ceil(200/lcc),1);
% Define legend entries
xkk=0:(1/(length(kk)-1)):1;

%% Beginning  of refinement plot for cdet
h=figure;
set(h,'Name', 'ICgpcmcdet plot', 'NumberTitle', 'off');
selIC2D=selIC(:,:,bestcshw);
plot1=plot(kk',selIC2D,'LineWidth',LineWidth);
title([nameselIC ' | $c_{shw}$=' num2str(cshw(bestcshw))],'Interpreter','latex')
% Add labels for the best value of cdet for each k
cmin=zeros(lcdet,1);
for j=1:lkk
    [~,posj]=min(selIC2D(j,:));
    cmin(j)=cdet(posj);
    text(xkk(j),0.98,['c_{det}=' num2str(cmin(j))],'Units','Normalized')
end

% Set line type and markers
set(plot1,{'LineStyle'},slintyp(1:lcdet));
set(plot1,{'Marker'},styp(1:lcdet))
xlabel('Number of groups')
set(gca,'xtick',kk)

a=cell(lcdet,1);
a(:)={'c_{det}='};
if isrow(cdet)
    legstrcdet=strcat(a, cellstr(num2str(cdet')));
else
    legstrcdet=strcat(a, cellstr(num2str(cdet')));
end
legend(legstrcdet,'location','best');

plot1cdet=gcf;
% Datatooltip mode (call to function ICplotLbl)
if ~isempty(datatooltip)
    PrepareDatatooltip(selIC,IDX)
end

%% Beginning  of refinement plot for cshw
h=figure;
set(h,'Name', 'ICgpcmcshw plot', 'NumberTitle', 'off');
selIC2D=squeeze(selIC(:,bestcdet,:));
plot1=plot(kk',selIC2D,'LineWidth',LineWidth);
title([nameselIC ' | $c_{det}=$' num2str(cdet(bestcdet))],'Interpreter','latex')
% Add labels for the best value of cshw for each k
cmin=zeros(lcshw,1);
for j=1:lkk
    [~,posj]=min(selIC2D(j,:));
    cmin(j)=cshw(posj);
    text(xkk(j),0.98,['c_{shw}=' num2str(cmin(j))],'Units','Normalized')
end

% Set line type and markers
set(plot1,{'LineStyle'},slintyp(1:lcshw));
set(plot1,{'Marker'},styp(1:lcshw))
xlabel('Number of groups')
set(gca,'xtick',kk)

a=cell(lcshw,1);
a(:)={'c_{shw}='};
if isrow(cshw)
    legstrcshw=strcat(a, cellstr(num2str(cshw')));
else
    legstrcshw=strcat(a, cellstr(num2str(cshw')));
end
legend(legstrcshw,'location','best');
% set(plot1,'Tag','BIC')
disp('The labels of in the top part of Figure named BIC denote the values of $c_{det}$ ($c_{shw}$) for which IC is minimum')

plot1cshw=gcf;
% Datatooltip mode (call to function ICplotLbl)
if ~isempty(datatooltip)
    PrepareDatatooltip(selIC,IDX)
end


%% Beginninig of heatmap
if isempty(databrush)
    figure
    selIC2Dbestk=squeeze(selIC(bestk,:,:));
    heatmap(cshw,cdet,selIC2Dbestk)
    xlabel('c_{shw}')
    ylabel('c_{det}')
    title(['Heatmap for k=' num2str(kk(bestk)) '. Best c_{shw}=' ...
        num2str(cshw(bestcshw))  ', best c_{det}=' num2str(cdet(bestcdet))...
        ' min(' nameselIC ')=' num2str(BICbest,'%1.1f') ])
end

%% Beginnning of refinemenet plot for cshb
candshb=IC.BICshb;
IDXshb=IC.IDXshb;
if size(candshb,1)>1
    h=figure;
    set(h,'Name', 'ICgpcmcshb plot', 'NumberTitle', 'off');
    hold('on')
    for i=1:size(candshb,1)
        plot(i,candshb(i,2),'o','LineWidth',5)
    end
    xlabel('c_{shb}')
    ylabel('BIC to select best c_{shb}')
    title(['Best c_{shb}=' num2str(IC.cshbbest)])
    
    set(gca,'xtick',1:size(candshb,1),'xticklabel',candshb(:,1))
    % Datatooltip mode (call to function ICplotLbl)
    if ~isempty(datatooltip)
        PrepareDatatooltip(candshb(:,2),IDXshb)
    end
else
    % if cshw is <=2 there is just one point and the plot is not shown
end
plot1cshb=gcf;

%% Beginnning of refinement plot for type or rotation
typerot={'I';'E';'V'};
h=figure;
set(h,'Name', 'ICgpcmROT plot', 'NumberTitle', 'off');

BICrot=IC.BICrot;
IDXrot=IC.IDXrot;

[~,indminrot]=min(BICrot);
hold('on')
for i=1:3
    plot(i,BICrot(i),'o','LineWidth',5)
end

xlabel('Type of rotation')
ylabel('BIC to select best type of rotation')
title(['Best rot =' typerot{indminrot}])
xlim([0.5 3.5])
set(gca,'XTick',1:3)
set(gca,'XTickLabel',char(typerot));

% Datatooltip mode (call to function ICplotLbl)
if ~isempty(datatooltip)
    PrepareDatatooltip(BICrot,IDXrot)
end
plot1crot=gcf;


% Make the main BIC plot the current figure
plBIC=findobj('type','figure','Name','BIC');
if ~isempty(plBIC)
    figure(plBIC(1))
end


%% Interactivity
if isstruct(databrush)
    fdatabrush=fieldnames(databrush);
    
    % persist option
    dpers=find(strcmp('persist',fdatabrush));
    if dpers>0
        %         persist=options.databrush(d+1);
        %         options.databrush(d:d+1)=[];
        persist=databrush.persist;
        databrush=rmfield(databrush,'persist');
        fdatabrush=fieldnames(databrush);
        
        ColorOrd=[1 0 0;0 1 1; 1 0 1; 1 1 0; 0 0 0; 0 1 0; 0 0 1];
        ColorOrd=repmat(ColorOrd,4,1);
    else
        persist='';
        ColorOrd=[1 0 0];
    end
    
    % FlagColor option Initialize colors for the brushing option: default
    % colors are blue (unbrushed unit) and red (brushed units)
    
else
    persist='';
    ColorOrd=[1 0 0];
end

%% Brush mode (call to function selectdataFS)
if ~isempty(databrush) || isstruct(databrush)
    
    if isempty(tag)
        set(gcf,'Tag','pl_IC');
    else
        set(gcf,'Tag',tag);
    end
    if isfield(IC,'Y')
        Y=IC.Y;
    else
        error('FSDA:tclustICplotGPCM:InvalidArg','Input struct must contain field ''Y''')
    end
    
    if isstruct(databrush)
        
        % Choose what to put on the main diagonal of the spm
        d=find(strcmp('dispopt',fieldnames(databrush)),1);
        if  isempty(d)
            dispopt='hist';
        else
            dispopt=databrush.dispopt;
            databrush=rmfield(databrush,'dispopt');
            fdatabrush=fieldnames(databrush);
        end
        
        
        % transform the input structure databrush into a cell array
        cv=[fdatabrush struct2cell(databrush)]';
        sele=cv(:)';
        % Add the FlagSize of the brushed points if it has not been
        % previously specified by the user
        d=find(strcmp('FlagSize',fdatabrush));
        if d==0
            sele=[sele 'FlagSize' '3'];
        end
    else
        sele={'selectionmode' 'Rect' };
        dispopt='hist';
    end
    
    % Set the labels of the axes in the spmplot which shows the brushed solutions.
    d=find(strcmp('nameY',fieldnames(IC)),1);
    if  isempty(d)
        v=size(Y,2);
        p1=1:v;
        nameY=cellstr(num2str(p1','y%d'));
    else
        nameY=IC.nameY;
    end
    
    plo=struct;
    plo.nameY=nameY;
    
    
    % some local variables
    but=0; brushcum=[]; ij=1;
    
    sele=[sele 'FlagColor' ColorOrd(ij,:) 'FlagMarker' char(styp(ij+1))];
    
    % loop brushing
    while but<=1
        
        % Remark: function selectdataFS cannot be used on the current
        % figure if the "selection mode" or the "zoom tool" are on. Setting
        % the plotedit mode initially to on and then to off, has the effect
        % to deselect plotedit mode.
        plotedit on
        plotedit off
        
        if strcmp(persist,'on')
            % add to cell sele option FlagColor (color of selection) and
            % FlagMarker (symbol to be used for selection)
            
            if ij>1
                chkexist=find(strcmp('FlagColor',sele)==1);
                sele{chkexist+1}=ColorOrd(ij,:);
                sele{chkexist+3}=char(styp(ij+1));
            end
        end
        
        
        % CALL to selectdataFS
        cascade
        prompt = ['\bf Enter type of plot you wish to brush:' newline ...
            '\bf 1= Analysis of restriction factor for determinants $c_{det}$' newline ...
            '\bf 2= Analysis of within groups restriction factor  $c_{shw}$' newline ...
            '\bf 3= Analysis of between groups restriction factor $c_{shb}$' newline ...
            '\bf 4= Analysis of type of rotation $c_{rot}$' newline ...
            '\bf else stop brushing'];
        
        dlgtitle = 'Choose the plot you wish to brush';
        dims = [1 80] ;
        definput = {'1'};
        opts.Interpreter='latex';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        answerd=str2double(cell2mat(answer));
        
        if  answerd ==1
            % The user has decided to investigage cdet
            Xcla=selIC(:,:,bestcshw);
            IDXselected=IDX(:,:,bestcshw);
            plot1=plot1cdet;
            cselected=IC.ccdet;
            disp('Select a region to brush in the IC plot for cdet');
            lab='cdet=';
            
        elseif answerd ==2
            % The user has decided to investigage cshw
            Xcla=squeeze(selIC(:,bestcdet,:));
            IDXselected=squeeze(IDX(:,bestcdet,:));
            plot1=plot1cshw;
            cselected=IC.ccdet;
            lab='cshw=';
            disp('Select a region to brush in the IC plot for cshw');
            
        elseif answerd ==3
            % the user has selected the figure which contains the
            % refinement for cshb
            Xcla= IC.BICshb(:,2)';
            IDXselected=IDXshb;
            plot1=plot1cshb;
            lab='cshb=';
            disp('Select a region to brush in the IC plot for cshb');
            cselected=IC.BICshb(:,1);
            
        elseif answerd ==4
            % the user has selected the figure which contains the
            % refinement for type of rotation
            Xcla=IC.BICrot';
            IDXselected=IDXrot;
            plot1=plot1crot;
            lab='type of rotation=';
            disp('Select a region to brush in the IC plot for crot');
            cselected=typerot;
        else
            return
        end
        
        
        % Remark: function selectdataFS cannot be used on the current
        % figure if the "selection mode" or the "zoom tool" are on. Setting
        % the plotedit mode initially to on and then to off, has the effect
        % to deselect plotedit mode.
        figure(plot1);
         plotedit on
        plotedit off
       
        [~,~,pl] = selectdataFS(sele{:});
        
        % exit if the ICplot figure was closed before selection
        if isnumeric(pl) && ~isempty(pl) && (pl == -999)
            return
        end
        
        if ~isempty(cell2mat(pl))
            
            nbrush=1;
            bic=cell2mat(pl);
            disp([num2str(length(bic)) ' solution(s) have been selected'])
            
        else
            disp('Wrong selection: Try again');
            disp('Select a region to brush in the IC plot');
            figure(plot1);
            nbrush='';
        end
        
        %% For each brushing operation, do the following:
        if ~isempty(nbrush)
            % brushcum = - the list of selected observations in all
            % iterations if
            %   persist=on
            % - the list of selected observations in the current iteration
            %   if persist=off
            if strcmp(persist,'on')
                brushcum=unique([brushcum; nbrush]);
            else
                brushcum=nbrush;
            end
            
            for r=1:length(bic)
                [Iwithk,Jwithc] = ind2sub(size(Xcla),find(Xcla==bic(r)));
                % It is necessary to put inside the tag the word group in
                % order to let spmplot understand that we are not dearling
                % with brushed units but simply with groups.
                
                figure('Tag','pl_spm');
                set(gcf,'WindowStyle',get(plot1,'WindowStyle'));
                % tabulate(IDX{Iwithk,Jwithc})
                if answerd<=2
                    spmplot(Y,cellstr(num2str(IDXselected{Iwithk,Jwithc})),plo,dispopt);
                    title([nameselIC ', k='  num2str(kk(Iwithk)) ' '  lab num2str(cselected(Jwithc) )])
                elseif answerd==3
                    spmplot(Y,IDXselected(:,Jwithc),plo,dispopt);
                    title([nameselIC ', k='  num2str(IC.kbest) ' '  lab num2str(cselected(Jwithc) )])
                elseif answerd==4
                    spmplot(Y,IDXselected(:,Jwithc),plo,dispopt);
                    title([nameselIC ', k='  num2str(IC.kbest) ' '  lab cselected{Jwithc} ])
                else
                end
                % spmplot(Y,IDX{Iwithk,Jwithc},plo)
            end
            
            %% - check condition to exit from the brush mode
            % If the option persistent is not equal off or on than get out
            % of the loop
            if strcmp('on',persist) || strcmp('off',persist)
                if strcmp('on',persist)
                    ij=ij+1;
                    % but=1 makes that previous highlightments in other
                    % figures are not deleted
                    but=1;
                end
                
                % Before waitforbuttonpress:
                % - the IC plot is highlighted again
                figure(plot1);
                % - and a function to be executed on figure close is set
                set(gcf,'CloseRequestFcn',@closereqFS);
                
                % Lay down the plots before continuing
                position(plot1);
                disp('Highlight the IC plot then: click on it to continue brushing or press a keyboard key to stop');
                ss=waitforbuttonpressFS;
                disp('------------------------');
                
                % After waitforbuttonpress:
                % - the standard MATLAB function to be executed on figure
                %   close is recovered
                set(gcf,'CloseRequestFcn','closereq');
                
                Open_spm = findobj(0, 'type', 'figure','tag','pl_spm');
                
                Open_mal = findobj(0, 'type', 'figure','tag','pl_IC');
                if isempty(Open_mal)  % User closed the main brushing window
                    if ~isempty(Open_spm)
                        delete(Open_spm);
                        % just in case Open_spm has length greater than 1
                        Open_spm=Open_spm(1); %#ok<NASGU>
                    end % spmplot  (yXplot) is deleted
                    delete(get(0,'CurrentFigure')); % deletes Figure if still one left open
                end
                
                % - and the 'but' variable is set if keyboard key was
                % pressed
                if ss==1
                    but=2;
                end
            else
                but=2;
            end % close loop associated with persist 'on' or persist 'off'
            position(plot1)
        end % for each brushing operation do ...
    end % close loop associated with but (loop brushing)
else
    % Apply cascade to existing plots in case databrush is not invoked
    position(0);
    
end

    function PrepareDatatooltip(IC,IDX)
        try
            chkgpu=gpuDevice; %#ok<NASGU>
            % datacursormode on;
            hdt = datacursormode;
            set(hdt,'Enable','on');
            % If options.datatooltip is not a struct then use our default options
            if ~isstruct(datatooltip)
                set(hdt,'DisplayStyle','window','SnapToDataVertex','on');
            else
                % options.datatooltip contains a structure where the user can set
                % the properties of the data cursor
                set(hdt,datatooltip);
            end
            
            LineColor=[1 0 0];
            % Declare a custom datatooltip update function to display additional
            % information about the selected unit
            
            %if contains(hdt.Figure.Name,'ICgpcm plot')
            set(hdt,'UpdateFcn',{@ICplotLbl,IC,IDX,LineColor});
            %elseif contains(hdt.Figure.Name,'ICgpcmROT plot')
            %    set(hdt,'UpdateFcn',{@ICplotLblROT,IC,IDX,LineColor});
            %else
            %end
        catch
            disp('No graphical device, interactive datatooltip not enabled')
        end
    end

    function output_txt = ICplotLbl(~,event_obj,selIC,IDX,~)
        % ICplotLbl provides information about
        %
        % Required input arguments:
        %
        %       obj =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object
        %               (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the
        %               sense that these arguments are automatically passed
        %               to the function when it executes.
        %       out =   a structure containing the following fields
        %               IC = a matrix containing the IC
        %
        % Output:
        %
        %   output_txt=  Datatip text (string or string cell array) which
        %                informs about the value of k and c which have been
        %                selected and the frequency distribution of the
        %                associated classification
        %
        % REMARK: this function is called by function tclustICplotGPCM
        %
        %
        % Written by FSDA team
        
        % find the plot the user has clicked on
        titl=get(gca,'title');
        titl=titl.String;
        
        if contains(titl, 'shw')
            % If the user has selected left panel (that is given bestcshw)
            ICsel=selIC(:,:,bestcshw);
            ICIDXsel=IDX(:,:,bestcshw);
            cloop=2;
        elseif contains(titl, 'det')
            % the user has selected right panel (that is given bestcdet)
            ICsel=squeeze(selIC(:,bestcdet,:));
            ICIDXsel=squeeze(IDX(:,bestcdet,:));
            cloop=1;
        elseif contains(titl, 'shb')
            % the user has selected the figure which contains the
            % refinement for cshb
            ICsel=selIC';
            ICIDXsel=IDX;
            cloop=3;
            
        elseif contains(titl, 'Best rot')
            % the user has selected the figure which contains the
            % refinement for type of rotation
            ICsel=selIC';
            ICIDXsel=IDX;
            cloop=4;
        else
            warning('FSDA:tclustICplotGPCM:InvalidArg','Supplied plot is not supported.')
            error('FSDA:tclustICplotGPCM:WrongIC','title of plot must be ''MIXMIX'' , ''MIXCLA'', ''CLACLA''')
        end
        
        pos = get(event_obj,'Position');
        
        % x and y, plot coordinates of the mouse
        x1 = pos(1); y1 = pos(2);
        
        % Find index to retrieve value of k (number of groups) and value of c (restriction factor)
        % Consider that find return the
        % linear indexing of matrix xydata
        idx = find(ICsel == y1,1);
        
        
        % Linear indexing is transformed into normal indexing using
        % function ind2sub row and column contain the column and row
        % indexed of the observation which has been selected with the mouse
        [row,col] = ind2sub(size(ICsel),idx);
        
        
        if ~isempty(hTarget)
            % set old line width and old color for old selection
            set(hTarget,'LineWidth',hTargetlwd,'Color',hTargetcol,'MarkerSize',hTargetms);
        else
        end
        
        % Store line width and color of selected trajectory
        % Notice that changing event_obj.Target in subsequent lines seems
        % to affect also hTarget
        hTarget=event_obj.Target;
        hTargetlwd=get(hTarget,'LineWidth');
        hTargetcol=get(hTarget,'Color');
        hTargetms=get(hTarget,'MarkerSize');
        
        % Increase Line width and keep line color
        set(hTarget,'LineWidth',hTargetlwd+1.5,'Color',hTargetcol);
        set(hTarget,'MarkerSize',hTarget.MarkerSize+5)
        
        if isempty(row)
            output_txt{1}=['no IC x,y' num2str(x1) '' num2str(y1)] ;
        else
            output_txt=cell(4,1);
            
            % output_txt is what it is shown on the screen
            output_txt(1,1) = {[ titl ' equal to: ',num2str(y1,4)]};
            
            if cloop== 1
                % Add information about k and cshw (user has selected right panel)
                output_txt{2,1} = ['$k= ' num2str(x1) ', c_{shw}=' num2str(IC.ccshw(col)) '$'];
                % Add information about the corresponding frequency
                % distribution  of associated classification
                clas=tabulate(ICIDXsel{x1,col});
                
            elseif cloop==2
                % Add information about k and cdet (user has selected left panel)
                output_txt{2,1} = ['$k= ' num2str(x1) ', c_{det}=' num2str(IC.ccdet(col)) '$'];
                % Add information about the corresponding frequency
                % distribution  of associated classification
                clas=tabulate(ICIDXsel{x1,col});
            elseif cloop==3
                output_txt{2,1} = ['$k= ' num2str(IC.kbest) ', c_{shb}=' num2str(candshb(col)) '$'];
                clas=tabulate(ICIDXsel(:,col));
                
            elseif cloop==4
                if col==1
                    typerotation='I';
                elseif col==2
                    typerotation='E';
                else
                    typerotation='V';
                end
                output_txt{2,1} = ['$k= ' num2str(IC.kbest) ', rot=' typerotation '$'];
                % Add information about the corresponding frequency
                % distribution  of associated classification
                clas=tabulate(ICIDXsel(:,col));
            end
            
            output_txt{3,1} = 'Classification';
            
            output_txt{4,1} = num2str(clas);
            
            set(0,'ShowHiddenHandles','on');    % Show hidden handles
            % hText = findobj('Type','text','Tag','DataTipMarker');
            % set(hText(1),'Interpreter','latex');
        end
    end

end

%FScategory:VIS-Clu


