function tclustICplot(IC,varargin)
%tclustICplot plots information criterion as a function of c and k
%
%<a href="matlab: docsearchFS('tclustICplot')">Link to the help function</a>
%
%   tclustICplot takes as input the output of function tclustIC (that is a
%   series of matrices which contain the values of the information criteria
%   BIC/ICL/CLA for different values of k and c) and plots them as function
%   of c or of k. The plot enables interaction in the sense that, if option
%   databrush has been activated, it is possible to click on a point in the
%   plot and to see the associated classification in the scatter plot
%   matrix.
%
%  Required input arguments:
%
%           IC : Information criterion to use. Structure.
%                It contains the following fields.
%                IC.CLACLA = matrix of size length(kk)-times
%                   length(cc) containinig the values of the penalized
%                   classification likelihood (CLA).
%                   This field is linked with
%                   out.IDXCLA.
%                IC.IDXCLA = cell of size length(kk)-times
%                   length(cc). Each element of the cell is a vector of
%                   length n containinig the assignment of each unit using
%                   the classification model.
%                Remark: fields CLACLA and IDXCLA are linked together.
%                   CLACLA and IDXCLA are compulsory just if optional input
%                   argument 'whichIC' is 'CLACLA' or 'ALL'
%                IC.MIXMIX = matrix of size length(kk)-times
%                   length(cc) containinig the value of the penalized
%                   mixture likelihood (BIC). This field is linked with
%                   out.IDXMIX.
%                IC.MIXCLA = matrix of size length(kk)-times
%                   length(cc) containinig the value of the ICL. This field
%                   is linked with out.IDXMIX.
%                IC.IDXMIX = cell of size length(kk)-times
%                   length(cc). Each element of the cell is a vector of
%                   length n containinig the assignment of each unit using
%                   the mixture model.
%                Remark 1: fields MIXMIX and IDXMIX are linked together.
%                   MIXMIX and IDXMIX are compulsory just if optional input
%                   argument 'whichIC' is 'CLACLA' or 'ALL'.
%                Remark 2: fields MIXCLA and IDXMIX are linked together.
%                   MIXCLA and IDXMIX are compulsory just if optional input
%                   argument 'whichIC' is 'MIXCLA' or 'ALL'.
%                IC.kk = vector containing the values of k (number of
%                   components) which have been considered.
%                IC.cc = vector containing the values of c (values of the
%                   restriction factor) which have been considered.
%                IC.Y =  original n-times-v data matrix on which the IC
%                   (Information criterion) has
%                    been computed
%                IC.nameY=  cell of length(size(Y,2)) containing the names
%                   of the variables of original matrix Y
%
%                 Data Types - struct
%
%  Optional input arguments:
%
%
%   whichIC  : character which specifies the information criterion to use
%               in the plot. Character.
%               Possible values for whichIC are:
%               'CLACLA' = in this case best solutions are referred to
%                   the classification likelihood.
%               'MIXMIX'  = in this case in this case best solutions are
%                   referred to the mixture likelihood (BIC).
%               'MIXCLA'  = in this case in this case best solutions are
%                   referred to ICL.
%               'ALL'  = in this case best solutions both three solutions using
%                     classification and mixture likelihood are produced.
%                   In output structure out all the three matrices
%                   out.MIXMIXbs, out.CLACLAbs and out.MIXCLAbs are given.
%               The default value of 'whichIC' is 'ALL'
%                 Example - 'whichIC','ALL'
%                 Data Types - character
%       tag     :   Personalized tag. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag 'pl_IC'.
%                   Note that if the program finds a plot which has a tag
%                   equal to the one specified by the user, then the output
%                   of the new plot overwrites the existing one in the same
%                   window else a new window is created.
%                   Example - 'tag','myplot'
%                   Data Types - char
%   datatooltip :   interactive clicking. Empty value (default) or
%                   structure. The default is datatooltip=''.
%                   If datatooltip = 1, the user can select with the
%                   mouse a solution in order to
%                   have the following information:
%                   1) value of k which has been selected
%                   2) value of c which has been selected
%                   3) values of the information criterion
%                   4) frequency distribution of the associated
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
%    databrush  :   interactive mouse brushing. empty value, scalar or structure.
%                   If databrush is an empty value (default), no brushing
%                   is done.
%                   The activation of this option (databrush is a scalar or
%                   a structure) enables the user  to select a set of
%                   values of IC in the current plot and to see the
%                   corresponding classification highlighted in the scatter
%                   plot matrix (spm).
%                   If spm does not exist it is automatically created.
%                   Please, note that the window style of the other figures is set
%                   equal to that which contains the IC
%                   plot. In other words, if the IC plot
%                   is docked all the other figures will be docked too.
%                   DATABRUSH IS A SCALAR.
%                   If databrush is a scalar the default selection tool is
%                   a rectangular brush and it is possible to brush only
%                   once (that is persist='').
%                   DATABRUSH IS A STRUCTURE.
%                   If databrush is a structure, it is possible to use all
%                   optional arguments of function selectdataFS.m and the
%                   following fields
%                   - databrush.persist = repeated brushng enabled. Persist is an empty value or a scalar
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
%  Output:
%
%
%
% See also: tclustIC, tclust
%
% References:
%
% A. Cerioli, L.A. Garcia-Escudero, A. Mayo-Iscar and M. Riani (2016). Finding
% the Number of Groups in Model-Based Clustering via Constrained
% Likelihoods, submitted.
% L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of
% Classification 2:193-218
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustICplot')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    %% Plot BIC, ICL and CLA for for Geyser data with all default options.
    Y=load('geyser2.txt');
    out=tclustIC(Y,'cleanpool',false,'plots',0,'alpha',0.1);
    tclustICplot(out)

%}

%{
    %   Example of the use of option datatooltip (all default options).
    %   Gives the user the possibility of clicking on the different points
    %   and have information about
    %   1) value of k which has been selected
    %   2) value of c which has been selected
    %   3) values of the information criterion
    %   4) frequency distribution of the associated classification
    tclustICplot(out,'datatooltip',1);
%}
%
%{
    % Example of the use of option datatooltip (personalized options).
    % Gives the user the possibility of clicking on the different points
    % and have information about the selected, the step of entry
    % into the subset and the associated label.
    datatooltip = struct;
    % In this example the style of the datatooltip is 'datatip'. Click on a
    % point when the ICplot is displayed.
    %
    datatooltip.DisplayStyle = 'datatip';
    tclustICplot(out,'datatooltip',datatooltip);
%}

%{
    % Simultaneous datatooltip with all 3 plots (MIXMIX, MIXCLA and CLACLA).
    tclustICplot(out,'whichIC','ALL')
%}

%{
    % Interactive_example (databrushing from the ICplot)
    % Use all default options for databrush (brush just once)
    tclustICplot(out,'databrush',1)
%}

%{
    % Interactive_example (databrushing from the ICplot)
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
    tclustICplot(out,'databrush',databrush)
%}



%% Beginning of code

% Specify information criterion for which you want to see best solutions
whichIC='MIXMIX';
% no brushsing as default
databrush='';
datatooltip=1;
tag='';

if nargin>1
    options=struct('whichIC',whichIC,...
        'datatooltip',datatooltip, 'tag',tag,'databrush', databrush,'nameY','');
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:tclustBICsol:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:tclustBICsol:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
    
    
    datatooltip=options.datatooltip;
    databrush=options.databrush;
    tag=options.tag;
    whichIC=options.whichIC;
end

% Extract the values of k (number of groups)
kk=IC.kk;
% Extract the values of c (number of components)
cc=IC.cc;
lcc=length(cc);

% set line width of the trajectories of BIC
LineWidth=1;
% Define marker type
styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
styp=repmat(styp,ceil(lcc/length(styp)),1);
% Define line type
slintyp={'-';'--';':';'-.'};
slintyp=repmat(slintyp,ceil(lcc/length(slintyp)),1);
% Define legend entries
a=cell(lcc,1);
a(:)={'c='};
if isrow(cc)
    legstr=strcat(a, cellstr(num2str(cc')));
else
    legstr=strcat(a, cellstr(num2str(cc')));
end
xkk=0:(1/(length(kk)-1)):1;


if strcmp(whichIC,'ALL')
    typeIC=3;
elseif strcmp(whichIC,'MIXMIX')
    typeIC=2;
elseif strcmp(whichIC,'MIXCLA')
    typeIC=1;
elseif strcmp(whichIC,'CLACLA')
    typeIC=0;
else
    warning('Supplied string for whichIC is not supported.')
    error('FSDA:tclustICplot:WrongIC','Specified information criterion is not supported: possible values are ''MIXMIX'' , ''MIXCLA'',  ''CLACLA'', ''ALL''')
end

% Global variables to declare is datatooltip is not empty
if ~isempty(datatooltip)
    hTarget=[];
    hTargetlwd=[];
    hTargetcol=[];
end

%% Prepate the figure to display the ICplot
% Create a figure to host the plot or clear the existing one
if typeIC<3
    if isempty(tag)
        h=figure;
    else
        h=findobj('-depth',1,'tag',tag);
        if (~isempty(h))
            clf(h);
            figure(h);
            axes;
        else
            h=figure;
        end
    end
    set(h,'Name', 'IC plot', 'NumberTitle', 'off');
    hold('all');
end

% CLACLA
if typeIC==0 || typeIC==3
    if typeIC==3
        figure
    end
    plot1CLACLA=plot(kk',IC.CLACLA,'LineWidth',LineWidth);
    title('CLACLA')
    % Add labels for the bet value of c for each k
    cmin=zeros(length(cc),1);
    for j=1:length(kk)
        [~,posj]=min(IC.CLACLA(j,:));
        cmin(j)=cc(posj);
        text(xkk(j),0.98,['c=' num2str(cmin(j))],'Units','Normalized')
    end
    
    % Set line type and markers
    set(plot1CLACLA,{'LineStyle'},slintyp(1:lcc));
    set(plot1CLACLA,{'Marker'},styp(1:lcc))
    xlabel('Number of groups')
    set(gca,'xtick',kk)
    legend(legstr,'location','best')
    plot1=gcf;
    % Datatooltip mode (call to function ICplotLbl)
    if ~isempty(datatooltip)
        PrepareDatatooltip(IC)
    end
end


% MIXMIX
if typeIC==2 || typeIC==3
    if typeIC==3
        figure
    end
    plot1MIXMIX=plot(kk',IC.MIXMIX,'LineWidth',LineWidth);
    title('MIXMIX')
    % Add labels for the bet value of c for each k
    cmin=zeros(length(cc),1);
    for j=1:length(kk)
        [~,posj]=min(IC.MIXMIX(j,:));
        cmin(j)=cc(posj);
        text(xkk(j),0.98,['c=' num2str(cmin(j))],'Units','Normalized')
    end
    % Set line type and markers
    set(plot1MIXMIX,{'LineStyle'},slintyp(1:lcc));
    set(plot1MIXMIX,{'Marker'},styp(1:lcc))
    xlabel('Number of groups')
    set(gca,'xtick',kk)
    legend(legstr,'location','best')
    plot1=gcf;
    % Datatooltip mode (call to function ICplotLbl)
    if ~isempty(datatooltip)
        PrepareDatatooltip(IC)
    end
end

%MIXCLA
if typeIC==1 || typeIC==3
    if typeIC==3
        figure
    end
    plot1MIXCLA=plot(kk',IC.MIXCLA,'LineWidth',LineWidth);
    title('MIXCLA')
    
    % Add labels for the best value of c for each k
    cmin=zeros(length(cc),1);
    for j=1:length(kk)
        [~,posj]=min(IC.MIXCLA(j,:));
        cmin(j)=cc(posj);
        text(xkk(j),0.98,['c=' num2str(cmin(j))],'Units','Normalized')
    end
    % Set line type and markers
    set(plot1MIXCLA,{'LineStyle'},slintyp(1:lcc));
    set(plot1MIXCLA,{'Marker'},styp(1:lcc))
    xlabel('Number of groups')
    set(gca,'xtick',kk)
    legend(legstr,'location','best')
    plot1=gcf;
    % Datatooltip mode (call to function ICplotLbl)
    if ~isempty(datatooltip)
        PrepareDatatooltip(IC)
    end
    
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
    if  typeIC==3
        sd=input(['Please choose the plot you wish to brush:\n Please type MIXMIX, MIXCLA or CLACLA\n' ...
            'Remark: to avoid this message please choose directly MIXMIX, MIXCLA or CLACLA in the optional input arguments\n'],'s');
        if strcmp(sd,MIXMIX)
            Xcla=IC.MIXMIX;
            IDX=IC.IDXMIX;
        elseif strcmp(sd,MIXCLA)
            Xcla=IC.MIXCLA;
            IDX=IC.IDXMIX;
        elseif strcmp(sd,CLACLA)
            Xcla=IC.CLACLA;
            IDX=IC.IDXCLA;
        else
            error('wrong option')
        end
    elseif typeIC==2
        Xcla=IC.MIXMIX;
        IDX=IC.IDXMIX;
    elseif typeIC==1
        Xcla=IC.MIXCLA;
        IDX=IC.IDXMIX;
    elseif typeIC==0
        
        Xcla=IC.CLACLA;
        IDX=IC.IDXCLA;
    else
    end
    
    if isempty(tag)
        set(gcf,'Tag','pl_IC');
    else
        set(gcf,'Tag',tag);
    end
    
    Y=IC.Y;
    
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
    
    % sele={sele{:} 'Tag' options.tag}; OLD inefficient code
    % MMMMMM
    % sele=[sele 'Tag' {options.tag}];
    
    
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
        
        figure(plot1);
        
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
        disp('Select a region to brush in the IC plot');
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
                spmplot(Y,cellstr(num2str(IDX{Iwithk,Jwithc})),plo,dispopt);
                % spmplot(Y,IDX{Iwithk,Jwithc},plo)
                title([whichIC '=' num2str(bic(r),'%10.2f') ', k='  num2str(kk(Iwithk))    ' c=' num2str(cc(Jwithc) )])
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
                    if ~isempty(Open_spm), delete(Open_spm); end % spmplot is deleted
                    delete(get(0,'CurrentFigure')); % deletes Figure if still one left open
                end
                
                % - and the 'but' variable is set if keyboard key was
                % pressed
                if ss==1;
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



    function PrepareDatatooltip(IC)
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
        set(hdt,'UpdateFcn',{@ICplotLbl,IC,LineColor});
        
    end

    function output_txt = ICplotLbl(~,event_obj,IC,~)
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
        % REMARK: this function is called by function tclustICplot
        %
        %
        % Written by FSDA team
        
        % find the plot the user has clicked on
        titl=get(gca,'title');
        titl=titl.String;
        
        if strcmp(titl,'MIXMIX')
            ICsel=IC.MIXMIX;
            ICIDXsel=IC.IDXMIX;
        elseif strcmp(titl,'MIXCLA')
            ICsel=IC.MIXCLA;
            ICIDXsel=IC.IDXMIX;
        elseif strcmp(titl,'CLACLA')
            ICsel=IC.CLACLA;
            ICIDXsel=IC.IDXCLA;
        else
            warning('Supplied plot is not supported.')
            error('FSDA:tclustICplot:WrongIC','title of plot must be ''MIXMIX'' , ''MIXCLA'', ''CLACLA''')
        end
        
        
        if ~isempty(hTarget)
            % set old line width and old color for old selection
            set(hTarget,'LineWidth',hTargetlwd,'Color',hTargetcol);
        else
        end
        
        % Store line width and color of selected trajectory
        % Notice that changing event_obj.Target in subsequent lines seems
        % to affect also hTarget
        hTarget=event_obj.Target;
        hTargetlwd=get(hTarget,'LineWidth');
        hTargetcol=get(hTarget,'Color');
        
        % Increase Line width and keep line color
        set(hTarget,'LineWidth',hTargetlwd+1.5,'Color',hTargetcol);
        
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
        
        if isempty(row)
            output_txt{1}=['no IC x,y' num2str(x1) '' num2str(y1)] ;
        else
            output_txt=cell(4,1);
            
            % output_txt is what it is shown on the screen
            output_txt(1,1) = {[ titl ' equal to: ',num2str(y1,4)]};
            
            % Add information about k and c
            % investigation
            output_txt{2,1} = ['k= ' num2str(x1) ', c=' num2str(IC.cc(col))];
            
            
            output_txt{3,1} = 'Classification';
            
            % Add information about the corresponding frequency
            % distribution  of associated classification
            clas=tabulate(ICIDXsel{x1,col});
            output_txt{4,1} = num2str(clas);
            
            set(0,'ShowHiddenHandles','on');    % Show hidden handles
            hText = findobj('Type','text','Tag','DataTipMarker');
            set(hText,'Interpreter','latex');
            
            
        end
    end

end




